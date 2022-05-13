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
    t<-sc3_prepare(sce, gene_filter = F, svm_max=10000)
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

# calculate the q()
def calculate_weight(xdata, params):
    rate, alpha, scale = params[0], params[1], 1 / params[2]
    mu, std = params[3], params[4]
    pz1 = params[0] * gamma.pdf(xdata, a=alpha, scale=scale)    # scale = 1 / beta
    pz2 = (1 - params[0]) * norm.pdf(xdata, mu, std)
    pz = pz1 / (pz1 + pz2)
    pz[pz1 == 0] = 0
    wt = [pz, 1-pz]
    return wt    # Represents the weight of gamma distribution and normal distribution


# Update parameters of gamma function
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


# Update the model with new parameters
def dmix(xdata, params):
    rate, alpha, beta, mu, std = params
    return rate * gamma.pdf(xdata, a=alpha, scale=(1 / beta)) + (1 - rate) * norm.pdf(xdata, mu, std)


# Returns the parameters of the mix model
def get_mix_internal(xdata, point):
    # "rate", "alpha", "beta", "mu", "sigma"
    params = [0, 0, 0, 0, 0]
    params[0] = np.sum(xdata == point) / len(xdata)
    if params[0] > 0.95:    # When a gene has 95% dropouts, we consider it invalid
        params = np.repeat(np.nan, 5)
        return params
    if params[0] == 0:
        params[0] = 0.01    # init rate = 0.01
    params[1], params[2] = 1.5, 1   # α，β = 0.5, 1
    xdata_rm = xdata[xdata > point]
    params[3], params[4] = np.mean(xdata_rm), np.std(xdata_rm)  # μ， σ
    if params[4] == 0:     # the variance is 0
        params[4] = 0.01
    eps = 10    # loss
    iter = 0    # iterations
    loglik_old = 0

    # The EM algorithm is continuously executed until the loss is reduced to less than 0.5
    while eps > 0.5:    
        wt = calculate_weight(xdata, params)
        tp_sum = np.sum(wt, axis=1)   
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


# Return the parameters of mixed distribution and invalid gene sequence number
def get_mix_parameters(count, point = np.log10(1.01)):
    # count represents each type of cell, with rows representing cell changes and columns representing gene changes
    count = np.array(count)
    # The genes with very low expression were screened out
    genes_expr = abs(np.mean(count, axis=0) - point)
    null_genes = np.argwhere(genes_expr < 1e-2)
    # paralist represents 5 parameters of each gene
    paralist = np.zeros([count.shape[1], 5])
    for i in tqdm(range(count.shape[1])):
        if i in null_genes:     # use nan to fill in invalid genes
            paralist[i] = np.repeat(np.nan, 5)
        else:
            xdata = count[:,i]    
            params = get_mix_internal(xdata, point)     # calculate the params of each columns
            paralist[i] = params
    return paralist, null_genes


# get the dropout rate of xdata
def dhat(xdata, params):
    rate, alpha, beta, mu, std = params
    gam = rate * gamma.pdf(xdata, a=alpha, scale=(1 / beta))
    nor = (1 - rate) * norm.pdf(xdata, mu, std)
    dropout = gam / (gam + nor)
    return dropout


#  Obtain the dropout rate of each gene of each cell 
def get_dropout_rate(count, point = np.log(1.01)):
    paralist, null_genes = get_mix_parameters(count, point)     
    df = pd.DataFrame(paralist, columns=['alpha', 'beta', 'mu', 'sigma', 'gamma'])
    print(df)
    dropout_rate = np.zeros((count.shape[1], count.shape[0]))    # init the matrix
    print("Calculating dropout rate...")
    for i in tqdm(range(count.shape[1])):
        if paralist[i][0] == np.nan:    # dropout rate should be `nan` if the gene is invalid
            dropout_rate[i] = np.repeat(dropout_rate.shape[1], np.nan)
        else:
            xdata = count[:, i]     
            dropout_rate[i] = dhat(xdata, paralist[i])  
    null_genes = null_genes.flatten()
    return dropout_rate.T, null_genes


# Obtain the input items and supervision items for denoising
def get_supervise(data,dropout,null_genes, M):
    nullarg = null_genes
    data = np.delete(data, nullarg, axis=1)
    data1 = data
    Omega = np.where(data > 0.5, 1, 0)  # 输入项
    dropout = np.delete(dropout, nullarg, axis=1)
    # In the M matrix, the larger the value, the more similar the description is
    dist = np.array(M)
    neibor = np.argsort(-dist, axis=1)
    neibor_data = map(lambda i:np.mean(data[neibor[i, 1:13]], axis=0), range(len(neibor)))
    neibor_data = list(neibor_data)
    for i in tqdm(range(data.shape[0])):
        for j in range(data.shape[1]):
            if dropout[i][j] > 0.2: 
                data[i][j] = dropout[i][j] * neibor_data[i][j] + (1 - dropout[i][j]) * data[i][j]
                Omega[i][j] = 0.5
    return Omega, data1


# Complete the matrix after neural network training
def makeup_results(result, data, null_genes, dropout, conservative=False):
    colarg = range(data.shape[1])
    nullarg = null_genes.flatten()
    validarg = list(set(colarg) - set(nullarg))
    for i in range(result.shape[0]):
        for j in range(result.shape[1]):
            if conservative:
                if dropout[i, validarg[j]] > 0.2:
                    data[i, validarg[j]] = result[i][j]
            elif dropout[i, validarg[j]] > 0:
                data[i, validarg[j]] = result[i][j]
    return data

# Complete the matrix after neural network training for two results
def makeup_results_all(result, data, null_genes, dropout):
    colarg = range(data.shape[1])
    nullarg = null_genes.flatten()
    validarg = list(set(colarg) - set(nullarg))
    cdata = np.array(data, copy = True)
    for i in range(result.shape[0]):
        for j in range(result.shape[1]):
            if dropout[i, validarg[j]] > 0.2:
                cdata[i, validarg[j]] = result[i][j]
            data[i, validarg[j]] = result[i][j]
    return data, cdata