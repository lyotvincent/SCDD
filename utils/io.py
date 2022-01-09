import numpy as np
import pandas as pd

cells = np.array([])
genes = np.array([])


def LoadData(expName, dataPath, labelPath=None, needTrans=True, labelHeader=None):
    """
    :rtype: tuple(np.array, np.array)
    :param expName: experiment name
    :param dataPath: the origin data path
    :param labelPath: the label path
    :param needTrans: if true, transpose the origin data.
    :param labelHeader: need to read the first line or not
    :return: training data and training label
    """
    df = pd.read_csv(dataPath, sep='\t', index_col=0)
    print("Experiment:{0}".format(expName))
    print("Data path:{0}".format(dataPath))
    print("Label path:{0}".format(labelPath))
    # set global values to save cells and genes, which will be used in SaveData
    global cells
    global genes
    cells = df.columns
    genes = df.index
    if needTrans:
        train_data = np.array(df).transpose()
    else:
        train_data = np.array(df)
    print("Data size:{0}".format(train_data.shape))
    label = list()
    if labelPath:
        train_label = np.array(pd.read_csv(labelPath, header=labelHeader))
        train_label = [elem[0] for elem in train_label]
        tp, l = dict(), 0
        cell_type = list()
        for elem in train_label:
            if elem in tp.keys():
                label.append(tp[elem])
            else:
                cell_type.append(elem)
                tp[elem] = l
                label.append(l)
                l += 1
    print("Load Data OK.")
    return train_data, label


def SaveData(expName, modelName, result, needTrans=True, batch = None):
    """
    This function will call `impute' member in model and save the result,
        which will use global cells and genes.
    Please make sure that genes are in index and cells are in rows.
    We will save the result in `results' dictionary
    :param needTrans:  if true, transpose the origin data.
    :param expName:  the experiment name
    :param modelName:  the model name
    :param result:  impute result
    """
    Batch = str("")
    if batch:
        Batch = str(batch)
    res = result
    if needTrans:
        res = pd.DataFrame(res.transpose())
    res.index = genes
    res.columns = cells
    path = "results/"+modelName+"/"+expName+"_"+modelName+Batch+"_impute.tsv"
    res.to_csv(path, sep='\t')
    print("Write to {0} successfully!".format(path))


def SaveTargets(M, Omega, Target, dropout_rate, null_genes):
    np.savez("temp/clusters_M", M)
    np.savez("temp/Omega", Omega)
    np.savez("temp/dropout_rate", dropout_rate)
    np.savez("temp/null_genes", null_genes)
    np.savez("temp/Target", Target)

def LoadTargets():
    M = np.load("temp/clusters_M.npz")
    M = M['arr_0']
    print("Load related Matrix `M' from temp...")
    Omega = np.load("temp/Omega.npz")
    Omega = Omega['arr_0']
    print("Load Omega from temp...")
    Target = np.load("temp/Target.npz")
    Target = Target['arr_0']
    print("Load Target from temp...")
    dropout_rate = np.load("temp/dropout_rate.npz")
    dropout_rate = dropout_rate['arr_0']
    print("Load dropout_rate from temp...")
    null_genes = np.load("temp/null_genes.npz")
    null_genes = null_genes['arr_0']
    print("Load null_genes from temp...")
    return M, Omega, Target, dropout_rate, null_genes