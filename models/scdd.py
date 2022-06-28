import tensorflow as tf
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from tqdm import tqdm
import random
from scipy.stats import spearmanr
l = 0
tf.compat.v1.disable_v2_behavior()


class SC_Denoising:
    def __init__(self, data, A, Omega, Target, batch_size=1024):
        global l
        self.l = l
        self.data = np.array(data, dtype=np.float32)
        self.batch_size = batch_size
        self.Omega = Omega
        self.N = np.sum(Omega > 0)
        self.Target = Target
        self.A = np.array(A, dtype=np.float32)
        self.preprocessing(self.data)
        self.build()
        self.sess = tf.compat.v1.Session()
        self.sess.run(tf.compat.v1.initialize_all_variables())

    def preprocessing(self, data):
        self.zeroarg = np.argwhere(np.all(data == 0, axis=0))[:, 0]
        self.data = np.delete(data, self.zeroarg, axis=1)
        self.factor = np.linalg.norm(self.data, ord=2, axis=1, keepdims=True)
        self.data = self.data / self.factor
        self.M = np.max(self.factor)
        self.factor = self.factor / self.M
        self.Target = self.Target / np.linalg.norm(self.Target, ord=2, axis=1, keepdims=True)
        if self.data.shape[0] <= 20000:
            print("the cells are less than 20000, using PCA with mode 'auto'...")
            pca = PCA(n_components=0.99, whiten=True)
        else:
            print("the cells are more than 20000, using PCA with mode 'arpack'...")
            if self.data.shape[1] <= 5000:
                pca = PCA(n_components=int(self.data.shape[1] / 2), whiten=True, svd_solver='arpack')
            else:
                pca = PCA(n_components=2500, whiten=True, svd_solver='arpack')
        self.input = pca.fit_transform(self.data)

        # calculate the total batches
        if self.input.shape[0] < self.batch_size:
            self.total_batch = 1
        else:
            self.total_batch = np.int((self.input.shape[0] - 1) / self.batch_size + 1)
        print(self.input.shape)
        print("total batch:{0}".format(self.total_batch))

    def bn_layer(self, x, scope, is_training, epsilon=0.001, decay=0.99, reuse=None):
        """
        Performs a batch normalization layer
        Args:
            x: input tensor
            scope: scope name
            is_training: python boolean value
            epsilon: the variance epsilon - a small float number to avoid dividing by 0
            decay: the moving average decay
        Returns:
            The ops of a batch normalization layer
        """
        with tf.compat.v1.variable_scope(scope, reuse=reuse):
            shape = x.get_shape().as_list()
            # gamma: a trainable scale factor
            gamma = tf.compat.v1.get_variable("SCDD_pca_" + scope + "_gamma", shape[-1],
                                              initializer=tf.constant_initializer(1.0), trainable=True)
            # beta: a trainable shift value
            beta = tf.compat.v1.get_variable("SCDD_pca_" + scope + "_beta", shape[-1],
                                             initializer=tf.constant_initializer(0.0), trainable=True)
            avg, var = tf.nn.moments(x, np.arange(len(shape) - 1), keepdims=True)
            avg = tf.reshape(avg, [avg.shape.as_list()[-1]])
            var = tf.reshape(var, [var.shape.as_list()[-1]])
            output = tf.nn.batch_normalization(x, avg, var, offset=beta, scale=gamma, variance_epsilon=epsilon)
        return output

    def bn_layer_top(self, x, is_training, epsilon=1e-7):
        """
        Returns a batch normalization layer that automatically switch between train and test phases based on the
        tensor is_training
        Args:
        x: input tensor
        scope: scope name
        is_training: boolean tensor or variable
        epsilon: epsilon parameter - see batch_norm_layer
        decay: epsilon parameter - see batch_norm_layer
        Returns:
        The correct batch normalization layer based on the value of is_training
        """
        # assert isinstance(is_training, (ops.Tensor, variables.Variable)) and is_training.dtype == tf.bool
        scope = 'scope' + str(self.l)
        self.l += 1
        return tf.cond(
            is_training,
            lambda: self.bn_layer(x=x, scope=scope, epsilon=epsilon, is_training=True, reuse=None),
            lambda: self.bn_layer(x=x, scope=scope, epsilon=epsilon, is_training=False, reuse=True),
        )

    def build(self):
        global l
        self.x = tf.compat.v1.placeholder(tf.float32, shape=(None, None))

        # if total_batch == 1, load A, Target, Omega as constant to accelerate
        if self.total_batch == 1:
            self.a = self.A
            self.t = self.Target
            self.o = self.Omega
        else:
            self.a = tf.compat.v1.placeholder(tf.float32, shape=(None, None))
            self.t = tf.compat.v1.placeholder(tf.float32, shape=(None, None))
            self.o = tf.compat.v1.placeholder(tf.float32, shape=(None, None))
        self.is_training = tf.compat.v1.placeholder(tf.bool)
        self.dropout_rate = tf.compat.v1.placeholder(tf.float32)
        W0 = tf.Variable(
            tf.random.uniform([self.input.shape[1], 128], minval=-(1.5 / (self.input.shape[1] + 128)) ** 0.5,
                              maxval=(1.5 / (self.input.shape[1] + 128)) ** 0.5))
        W1 = tf.Variable(tf.random.uniform([128, 64], minval=-(6 / 192) ** 0.5, maxval=(6 / 192) ** 0.5))
        W2 = tf.Variable(tf.random.uniform([64, 16], minval=-(6 / 80) ** 0.5, maxval=(6 / 80) ** 0.5))
        self.layer1 = tf.nn.relu(
            self.bn_layer_top(tf.matmul(tf.matmul(self.a, self.x), W0), is_training=self.is_training))
        self.layer2 = tf.nn.relu(
            self.bn_layer_top(tf.matmul(tf.matmul(self.a, self.layer1), W1), is_training=self.is_training))
        self.z = tf.matmul(tf.matmul(self.a, self.layer2), W2)
        self.embedding = self.z
        self.embedding = self.bn_layer_top(self.embedding, is_training=self.is_training)
        W3 = tf.Variable(tf.random.uniform([16, 64], minval=-(6 / 80) ** 0.5, maxval=(6 / 80) ** 0.5))
        W4 = tf.Variable(tf.random.uniform([64, 128], minval=-(6 / 192) ** 0.5, maxval=(6 / 192) ** 0.5))
        W5 = tf.Variable(tf.random.uniform([128, 1024], minval=-(6 / 1152) ** 0.5, maxval=(6 / 1152) ** 0.5))
        W6 = tf.Variable(
            tf.random.uniform([1024, self.Target.shape[1]], minval=-(1.5 / (self.Target.shape[1] + 1024)) ** 0.5,
                              maxval=(1.5 / (self.Target.shape[1] + 1024)) ** 0.5))
        self.layer3 = tf.nn.relu(self.bn_layer_top(tf.matmul(self.embedding, W3), is_training=self.is_training))
        self.layer4 = tf.nn.relu(self.bn_layer_top(tf.matmul(self.layer3, W4), is_training=self.is_training))
        self.layer5 = tf.nn.relu(self.bn_layer_top(tf.matmul(self.layer4, W5), is_training=self.is_training))
        self.layer6 = tf.nn.dropout(self.layer5, rate=self.dropout_rate)
        self.res = tf.sigmoid(tf.matmul(self.layer6, W6))
        self.contractive_loss = tf.reduce_sum(tf.square(tf.gradients(self.z, self.x, stop_gradients=[self.x])))
        self.real_loss = tf.reduce_mean(tf.reduce_sum(((self.t - self.res) * self.o) ** 2, axis=1))
        self.loss = tf.reduce_mean(
            tf.reduce_sum(((self.t - self.res)) ** 2, axis=1)) + self.contractive_loss + self.real_loss
        self.opt = tf.compat.v1.train.RMSPropOptimizer(1e-3).minimize(self.loss)
        l += self.l

    def shuffle_set(self):
        train_row = list(range(self.input.shape[0]))
        random.shuffle(train_row)
        return train_row

    def get_batch(self, train_row, now_batch):
        if now_batch < self.total_batch - 1:
            batch = train_row[now_batch * self.batch_size:(now_batch + 1) * self.batch_size]
        else:
            batch = train_row[now_batch * self.batch_size:]
        cell_batch = self.input[batch]
        target_batch = self.Target[batch]
        omega_batch = self.Omega[batch]
        A_batch = self.A[batch, :]
        A_batch = spearmanr(A_batch.transpose())[0]
        A_batch = np.where(np.diag([1] * len(batch)) == 1, np.diag([np.min(A_batch)] * len(batch)), A_batch)
        A_batch = (A_batch - np.min(A_batch)) / (np.max(A_batch) - np.min(A_batch))
        A_batch = np.where(np.diag([1] * len(batch)) == 1, np.diag([1] * len(batch)), A_batch)
        D = np.diag(np.where(np.sum(A_batch, axis=0) == 0, [1] * len(A_batch), np.sum(A_batch, axis=0)) ** -0.5)
        A_batch = D * np.mat(A_batch) * D
        A_batch = np.array(A_batch)
        return cell_batch, A_batch, omega_batch, target_batch

    def train(self, n):
        print("Using model:SCDD...")
        min_loss = np.inf
        for _ in tqdm(range(n)):
            Loss = []
            if self.total_batch == 1:
                a, b, c = self.sess.run((self.opt, self.loss, self.contractive_loss),
                                        feed_dict={self.x: self.input, self.is_training: True, self.dropout_rate: 0.5})
                Loss.append(b)
                # print(b-c,c)
            else:
                train_row = self.shuffle_set()
                for now_batch in range(self.total_batch):
                    cell_batch, A_batch, omega_batch, target_batch = self.get_batch(train_row, now_batch)
                    a, b, c = self.sess.run((self.opt, self.loss, self.contractive_loss),
                                            feed_dict={self.x: cell_batch, self.a: A_batch, self.o: omega_batch,
                                                       self.t: target_batch, self.is_training: True,
                                                       self.dropout_rate: 0.5})
                    Loss.append(b)
            print("Epoch {0} Loss: {1}".format(_, np.mean(Loss)))

    def impute(self):
        if self.total_batch == 1:
            result = self.sess.run((self.res),
                                   feed_dict={self.x: self.input, self.is_training: False, self.dropout_rate: 0})
        else:
            reslist = []
            train_row = self.shuffle_set()
            s = pd.DataFrame(train_row)
            s[1] = s.index
            s.index = s[0]
            s = s.sort_index()
            rev_row = list(s[1])
            for now_batch in range(self.total_batch):
                cell_batch, A_batch, omega_batch, target_batch = self.get_batch(train_row, now_batch)
                res = self.sess.run((self.res), feed_dict={self.x: cell_batch, self.a: A_batch, self.o: omega_batch,
                                                           self.t: target_batch, self.is_training: False,
                                                           self.dropout_rate: 0})
                reslist.append(res)
            result = np.vstack(reslist)
            result = result[rev_row]
            print(result.shape)
        result = np.array(result) * self.factor * self.M
        return np.array(result)
