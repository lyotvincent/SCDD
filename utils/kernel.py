import numpy as np
from scipy.stats import spearmanr


def getA(trainData, method="TFIDF", filter=True, structure="Array"):
    data = trainData > 0
    print("using method:{0} to get A...".format(method))
    if method == "TFIDF":
        if filter:
            data = data[:, (np.sum(data, axis=0) > 0.1 * len(data)) & (np.sum(data, axis=0) < 0.9 * len(data))]
        else:
            data = trainData
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

def getA_batch(C, index):
    C = C[:, index]
    batch_size = len(index)
    A = []
    for i in range(batch_size):
        tp = C[:, i:i+1]
        tp2 = (C == tp)
        A.append(np.sum(tp2, axis=0))
    D = np.diag(np.where(np.sum(A, axis=0) == 0, [1] * len(A), np.sum(A, axis=0)) ** -0.5)
    A = D * np.mat(A) * D
    A = np.array(A)
    return A