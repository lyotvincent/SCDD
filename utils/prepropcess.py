import numpy as np
import pandas as pd
from scipy.stats import gamma, norm
from scipy.optimize import root
from scipy.special import digamma
from tqdm import tqdm
import rpy2.robjects as robjects
import os
import sys

if sys.platform == 'win32':
    print("using Windows system...")
    os.environ["R_HOME"] = r"D:\R\R-4.1.0"
    os.environ["PATH"] = r"D:\R\R-4.1.0\bin\x64" + ";" + os.environ["PATH"]
elif sys.platform == 'linux':
    print("using Linux system...")


def RobjKmeans(train_data: np.array) -> np.array:
    """
    use Kmeans in R object, by reading training data in files, and returns
    the related Matrix in M
    :param train_data: a Matrix with genes as columns and cells as rows
    :param store: if true, save M as .npz file
    """
    np.savetxt('temp/train_data.txt', train_data, fmt="%f," * len(train_data[0]))
    print("Converted data to R-object and start R-engine...")
    M, C = robjects.r('''
    library('SC3')
    a<-read.table('temp/train_data.txt', header=F,sep=',')
    tp<-unlist(a)
    tp<-array(tp,dim=c(%d,%d))
    tp<-t(tp)
    tp<-tp[which(rowSums(tp)>0),]
    tplog<-log2(tp+1)
    library(SingleCellExperiment)
    sce<-SingleCellExperiment(assays=list(counts=tp,logcounts=tplog))
    ks<-sc3_estimate_k(sce)@metadata$sc3$k_estimation
    rowData(sce)$feature_symbol=1:dim(tp)[1]
    t<-sc3_prepare(sce, svm_max=10000)
    t<-sc3_calc_dists(t)
    t<-sc3_calc_transfs(t)
    t<-sc3_kmeans(t,ks=ks:ks)
    t<-sc3_calc_consens(t)
    info <- metadata(t)$sc3$consensus[[as.character(ks), exact = FALSE]]
    con <- info$consensus
    clu <- info$silhouette[,1]
    return (list(con, clu))
    ''' % (train_data.shape[0], train_data.shape[1]))
    M = np.array(M)
    return M

# 计算q函数
def calculate_weight(xdata, params):
    rate, alpha, scale = params[0], params[1], 1 / params[2]
    mu, std = params[3], params[4]
    pz1 = params[0] * gamma.pdf(xdata, a=alpha, scale=scale)    # scale = 1 / beta
    pz2 = (1 - params[0]) * norm.pdf(xdata, mu, std)
    pz = pz1 / (pz1 + pz2)
    pz[pz1 == 0] = 0
    wt = [pz, 1-pz]
    return wt    # 表示gamma分布和normal分布占的权重


# 更新gamma函数的参数
def update_gamma_pars(xdata, wt):
    xdata = np.array(xdata)
    wt = np.array(wt)
    tp_s = np.sum(wt)
    tp_t = np.sum(wt * xdata)
    tp_u = sum(wt * np.log(xdata))
    tp_v = -tp_u / tp_s - np.log(tp_s / tp_t)
    alpha = 0
    if tp_v <= 0:
        alpha = 20
    else:
        alpha0 = (3 - tp_v + np.sqrt((tp_v - 3) ** 2 + 24 * tp_v)) / 12 / tp_v
        if alpha0 >= 20:
            alpha = 20
        else:
            alpha = root(lambda alpha: np.log(alpha) - digamma(alpha) - tp_v, np.array([0.9, 1.1])*alpha0)
            alpha = alpha.x[0]
    ## need to solve log(x) - digamma(x) = tp_v
    ## We use this approximation to compute the initial value
    beta = tp_s / tp_t * alpha
    return alpha, beta


# 用新的参数更新模型
def dmix(xdata, params):
    rate, alpha, beta, mu, std = params
    return rate * gamma.pdf(xdata, a=alpha, scale=(1 / beta)) + (1 - rate) * norm.pdf(xdata, mu, std)


# 返回混合模型的参数
def get_mix_internal(xdata, point):
    # "rate", "alpha", "beta", "mu", "sigma"
    params = [0, 0, 0, 0, 0]
    params[0] = np.sum(xdata == point) / len(xdata)
    if params[0] > 0.95:    # 当某种基因有95%的缺失值时，我们认为该基因无效
        params = np.repeat(np.nan, 5)
        return params
    if params[0] == 0:
        params[0] = 0.01    # 初始化时将rate设置为0.01
    params[1], params[2] = 1.5, 1   # α，β分别为0.5, 1
    xdata_rm = xdata[xdata > point]
    params[3], params[4] = np.mean(xdata_rm), np.std(xdata_rm)  # μ， σ分别为均值和标准差
    if params[4] == 0:     # 处理方差为0的情况
        params[4] = 0.01
    eps = 10    # 损失值
    iter = 0    # 迭代次数
    loglik_old = 0

    while eps > 0.5:    # 不断执行EM算法知道损失减小至0.5以下
        wt = calculate_weight(xdata, params)
        tp_sum = np.sum(wt, axis=1)     # 分别将每个细胞的gamma权重和norm权重相加
        rate = tp_sum[0] / len(wt[0])
        mu = np.sum(wt[1] * xdata) / np.sum(wt[1])
        std = np.sqrt(np.sum(wt[1] * ((xdata - mu ) ** 2)/sum(wt[1])))
        alpha, beta = update_gamma_pars(xdata, wt[0])
        params = [rate, alpha, beta, mu, std]
        new = dmix(xdata, params)
        loglik = np.sum(np.log10(new))
        eps = (loglik - loglik_old) ** 2
        loglik_old = loglik
        iter = iter + 1
        if iter > 100:
            break
    return params


# 返回混合分布的参数，以及无效的基因序号
def get_mix_parameters(count, point = np.log10(1.01)):
    # count 表示每一类细胞的矩阵，行代表细胞的变化，列代表基因的变化
    count = np.array(count)
    # 筛选出表达量极低的基因
    genes_expr = abs(np.mean(count, axis=0) - point)
    null_genes = np.argwhere(genes_expr < 1e-2)
    # paralist 表示对于每一个基因，对应的5个参数
    paralist = np.zeros([count.shape[1], 5])
    # 循环该类的细胞矩阵的每一个基因
    for i in tqdm(range(count.shape[1])):
        if i in null_genes:     # 无效基因的参数都用nan来填充
            paralist[i] = np.repeat(np.nan, 5)
        else:
            xdata = count[:,i]    # 将该列数据挑选出来
            params = get_mix_internal(xdata, point)     # 针对该列数据，计算混合模型参数
            paralist[i] = params
    return paralist, null_genes


# 获取xdata的dropout率
def dhat(xdata, params):
    rate, alpha, beta, mu, std = params
    gam = rate * gamma.pdf(xdata, a=alpha, scale=(1 / beta))
    nor = (1 - rate) * norm.pdf(xdata, mu, std)
    dropout = gam / (gam + nor)
    return dropout


# 获取该聚类每一个细胞的每一个基因是dropout的概率
def get_dropout_rate(count, point = np.log(1.01)):
    paralist, null_genes = get_mix_parameters(count, point)     # 获取混合分布模型参数
    df = pd.DataFrame(paralist, columns=['alpha', 'beta', 'mu', 'sigma', 'gamma'])
    print(df)
    dropout_rate = np.zeros((count.shape[1], count.shape[0]))    # 初始化矩阵，行表示基因的变化，列表示细胞的变化
    print("Calculating dropout rate...")
    for i in tqdm(range(count.shape[1])):
        if paralist[i][0] == np.nan:    # 当该列基因被认为无效，dropout概率也为nan
            dropout_rate[i] = np.repeat(dropout_rate.shape[1], np.nan)
        else:
            xdata = count[:, i]     # 获取该基因对应的细胞
            dropout_rate[i] = dhat(xdata, paralist[i])  # 获取dropout率
    null_genes = null_genes.flatten()
    return dropout_rate.T, null_genes


# 获取神经网络框架的输入项和监督项
def get_supervise(data,dropout,null_genes, M):
    nullarg = null_genes
    # 删除无效的基因
    data = np.delete(data, nullarg, axis=1)
    data1 = data
    Omega = np.where(data > 0.5, 1, 0)  # 输入项
    dropout = np.delete(dropout, nullarg, axis=1)
    # 获得关于data的M矩阵，在M矩阵当中，值越大说明越相似
    dist = np.array(M)
    # 对于每一个dropout点位，先用与其距离最近的10个点位的均值进行填充
    neibor = np.argsort(-dist, axis=1)
    neibor_data = map(lambda i:np.mean(data[neibor[i, 1:13]], axis=0), range(len(neibor)))
    neibor_data = list(neibor_data)
    for i in tqdm(range(data.shape[0])):
        for j in range(data.shape[1]):
            if dropout[i][j] > 0.2:     # 当dropout_rate 大于 0.5，我们就认为是dropout点位
                # 这里采用了原始值和10个邻近点的中位数值的加权和，作为简单填充
                data[i][j] = dropout[i][j] * neibor_data[i][j] + (1 - dropout[i][j]) * data[i][j]
                Omega[i][j] = 0.5
    return Omega, data1


# 补全神经网络训练后的矩阵
def makeup_results(result, data, null_genes, dropout, conservative=False):
    # 先根据null_genes, 找到有效基因在data中的列号
    colarg = range(data.shape[1])
    nullarg = null_genes.flatten()
    validarg = list(set(colarg) - set(nullarg))
    # 再将result的每一列填到data的对应列中
    for i in range(result.shape[0]):
        for j in range(result.shape[1]):
            # 这里对于每一列的回填时，需要查看dropout率，只有当dropout率过高时才回填
            if conservative:
                if dropout[i, validarg[j]] > 0.2:
                    data[i, validarg[j]] = result[i][j]#  * dropout[i, validarg[j]] + (1 - dropout[i, validarg[j]]) * data[i, validarg[j]]
            elif dropout[i, validarg[j]] > 0:
                data[i, validarg[j]] = result[i][j]# * dropout[i, validarg[j]] + (1 - dropout[i, validarg[j]]) * data[i, validarg[j]]
    # 将填充后的data作为result返回
    return data

# 补全神经网络训练后的矩阵
def makeup_results_all(result, data, null_genes, dropout):
    # 先根据null_genes, 找到有效基因在data中的列号
    colarg = range(data.shape[1])
    nullarg = null_genes.flatten()
    validarg = list(set(colarg) - set(nullarg))
    cdata = np.array(data, copy = True)
    # 再将result的每一列填到data的对应列中
    for i in range(result.shape[0]):
        for j in range(result.shape[1]):
            # 这里对于每一列的回填时，需要查看dropout率，只有当dropout率过高时才回填
            if dropout[i, validarg[j]] > 0.2:
                cdata[i, validarg[j]] = result[i][j]# * dropout[i, validarg[j]] + (1 - dropout[i, validarg[j]]) * data[i, validarg[j]]
            data[i, validarg[j]] = result[i][j]# * dropout[i, validarg[j]] + (1 - dropout[i, validarg[j]]) * data[i, validarg[j]]
    # 将填充后的data作为result返回
    return data, cdata