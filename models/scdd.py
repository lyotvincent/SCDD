import tensorflow as tf
tf.compat.v1.disable_v2_behavior()
import numpy as np
from sklearn.decomposition import PCA
from tqdm import tqdm

l = 0
class SC_Denoising:
    def __init__(self,data, A, Omega, Target):
        global l
        self.l=l
        self.data=np.array(data,dtype=np.float32)
        self.Omega=Omega
        self.N = np.sum(Omega > 0)
        self.Target = Target
        self.A=np.array(A,dtype=np.float32)
        self.preprocessing(self.data)
        self.build()
        self.sess = tf.compat.v1.Session()
        self.sess.run(tf.compat.v1.initialize_all_variables())
    def preprocessing(self,data):
        self.zeroarg = np.argwhere(np.all(data == 0, axis=0))[:, 0]
        self.data = np.delete(data, self.zeroarg, axis=1)
        self.factor = np.linalg.norm(self.data, ord=2, axis=1, keepdims=True)
        self.data = self.data / self.factor
        self.M = np.max(self.factor)
        self.factor = self.factor / self.M
        self.Target = self.Target / np.linalg.norm(self.Target, ord=2, axis=1, keepdims=True)
        pca = PCA(n_components=0.99,whiten=True)
        self.input = pca.fit_transform(self.data)
        print(self.input.shape)
    def bn_layer(self,x, scope, is_training, epsilon=0.001, decay=0.99, reuse=None):
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
            gamma = tf.compat.v1.get_variable("SCDD_pca_"+scope+"_gamma", shape[-1], initializer=tf.constant_initializer(1.0), trainable=True)
            # beta: a trainable shift value
            beta = tf.compat.v1.get_variable("SCDD_pca_"+scope+"_beta", shape[-1], initializer=tf.constant_initializer(0.0), trainable=True)
            avg, var = tf.nn.moments(x, np.arange(len(shape)-1), keepdims=True)
            avg=tf.reshape(avg, [avg.shape.as_list()[-1]])
            var=tf.reshape(var, [var.shape.as_list()[-1]])
            output = tf.nn.batch_normalization(x, avg, var, offset=beta, scale=gamma, variance_epsilon=epsilon)
        return output
    def bn_layer_top(self,x, is_training, epsilon=1e-7):
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
        #assert isinstance(is_training, (ops.Tensor, variables.Variable)) and is_training.dtype == tf.bool
        scope='scope'+str(self.l)
        self.l+=1
        return tf.cond(
            is_training,
            lambda: self.bn_layer(x=x, scope=scope, epsilon=epsilon, is_training=True, reuse=None),
            lambda: self.bn_layer(x=x, scope=scope, epsilon=epsilon, is_training=False, reuse=True),
            )
    def build(self):
        global l
        self.x = tf.compat.v1.placeholder(tf.float32, shape=(None, None))
        self.is_training = tf.compat.v1.placeholder(tf.bool)
        self.dropout_rate = tf.compat.v1.placeholder(tf.float32)
        W0 = tf.Variable(tf.random.uniform([self.input.shape[1], 128], minval=-(1.5/(self.input.shape[1]+128))**0.5,maxval=(1.5/(self.input.shape[1]+128))**0.5))
        W1 = tf.Variable(tf.random.uniform([128, 64], minval=-(6/192)**0.5,maxval=(6/192)**0.5))
        W2 = tf.Variable(tf.random.uniform([64, 16], minval=-(6/80)**0.5,maxval=(6/80)**0.5))
        self.layer1=tf.nn.relu(self.bn_layer_top(tf.matmul(tf.matmul(self.A,self.x),W0),is_training=self.is_training))
        self.layer2=tf.nn.relu(self.bn_layer_top(tf.matmul(tf.matmul(self.A,self.layer1),W1),is_training=self.is_training))
        self.z = tf.matmul(tf.matmul(self.A,self.layer2),W2)
        self.embedding = self.z
        self.embedding = self.bn_layer_top(self.embedding,is_training=self.is_training)
        W3 = tf.Variable(tf.random.uniform([16, 64], minval=-(6/80)**0.5,maxval=(6/80)**0.5))
        W4 = tf.Variable(tf.random.uniform([64, 128], minval=-(6/192)**0.5,maxval=(6/192)**0.5))
        W5 = tf.Variable(tf.random.uniform([128, 1024], minval=-(6/1152)**0.5,maxval=(6/1152)**0.5))
        W6 = tf.Variable(tf.random.uniform([1024,self.Target.shape[1]], minval=-(1.5/(self.Target.shape[1]+1024))**0.5,maxval=(1.5/(self.Target.shape[1]+1024))**0.5))
        self.layer3=tf.nn.relu(self.bn_layer_top(tf.matmul(self.embedding,W3),is_training=self.is_training))
        self.layer4=tf.nn.relu(self.bn_layer_top(tf.matmul(self.layer3,W4),is_training=self.is_training))
        self.layer5=tf.nn.relu(self.bn_layer_top(tf.matmul(self.layer4,W5),is_training=self.is_training))
        self.layer6 = tf.nn.dropout(self.layer5, rate=self.dropout_rate)
        self.res=tf.sigmoid(tf.matmul(self.layer6,W6))
        self.contractive_loss = tf.reduce_sum(tf.square(tf.gradients(self.z, self.x, stop_gradients = [self.x])))
        self.real_loss = tf.reduce_mean(tf.reduce_sum(((self.Target-self.res) * self.Omega)**2, axis=1))
        self.loss = tf.reduce_mean(tf.reduce_sum(((self.Target-self.res))**2, axis=1))+ self.contractive_loss + self.real_loss
        self.opt = tf.compat.v1.train.RMSPropOptimizer(1e-3).minimize(self.loss)
        l +=self.l
    def train(self,n, label=None):
        print("Using model:SCDD...")
        for _ in tqdm(range(n)):
            a,b,c=self.sess.run((self.opt,self.loss,self.contractive_loss),feed_dict={self.x:self.input,self.is_training:True, self.dropout_rate:0.5})
            # print(b-c,c)
    def impute(self):
        res=self.sess.run((self.res),feed_dict={self.x:self.input,self.is_training:False, self.dropout_rate:0})
        res = np.array(res) * self.factor * self.M
        return np.array(res)