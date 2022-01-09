import numpy as np
from scipy.stats import spearmanr


def getA(trainData, method="TFIDF", M:np.array=None):
    data = trainData>0
    print("using method:{0} to get A...".format(method))
    if method == "TFIDF":
        data = data[:, (np.sum(data, axis=0) > 0.1 * len(data)) & (np.sum(data, axis=0) < 0.9 * len(data))]
        data = np.log(1.01 + data)
        IDF = np.log(len(data) / (np.sum(data > 0, axis=0) + 1))
        TF = data / np.sum(data, axis=1, keepdims=True)
        A = spearmanr((TF * IDF).transpose())[0]
        A = np.where(np.diag([1] * len(data)) == 1, np.diag([np.min(A)] * len(data)), A)
        A = (A - np.min(A)) / (np.max(A) - np.min(A))
        A = np.where(np.diag([1] * len(data)) == 1, np.diag([1] * len(data)), A)
        D = np.diag(np.where(np.sum(A, axis=0) == 0, [1] * len(A), np.sum(A, axis=0)) ** -0.5)
        A = D * np.mat(A) * D
    A = np.array(A)
    return A