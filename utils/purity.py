import numpy as np
from sklearn.metrics import accuracy_score
def purity_score(y_true, y_pred):
    y_true = np.array(y_true).T
    y_pred = np.array(y_pred).T
    # matrix which will hold the majority-voted labels
    y_voted_labels = np.zeros(y_true.shape)
    # Ordering labels
    ## Labels might be missing e.g with set like 0,2 where 1 is missing
    ## First find the unique labels, then map the labels to an ordered set
    ## 0,2 should become 0,1
    labels = np.unique(y_true)
    ordered_labels = np.arange(labels.shape[0])
    for k in range(labels.shape[0]):
        y_true[y_true == labels[k]] = ordered_labels[k]
    # Update unique labels
    labels = np.unique(y_true)
    # We set the number of bins to be n_classes+2 so that
    # we count the actual occurence of classes between two consecutive bins
    # the bigger being excluded [bin_i, bin_i+1[
    bins = np.concatenate((labels, [np.max(labels) + 1]), axis=0)

    for cluster in np.unique(y_pred):
        hist, _ = np.histogram(y_true[y_pred == cluster], bins=bins)
        # Find the most present label in the cluster
        winner = np.argmax(hist)
        y_voted_labels[y_pred == cluster] = winner

    return accuracy_score(y_true, y_voted_labels)


# Calculate purity from `predicts`, which is the cluster results of each methods
# and the `perfs` represents other performances such as ARI, NMI, F1-score, etc..
def cal_purity(perfs, predicts, mcounts):
    label = predicts.iloc[:, 0]
    pur = []
    for i in range(1, mcounts+2):
        pred = predicts.iloc[:, i]
        pur.append(purity_score(label, pred))
    perfs['Purity'] = pur
    return perfs