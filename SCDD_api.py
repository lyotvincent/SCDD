import numpy as np
import scanpy as sc
from utils.io import *
from utils.kernel import *
from models.scdd import *
from utils.prepropcess import get_dropout_rate, get_supervise, makeup_results, makeup_results_all, RobjKmeans, \
    AnnDataSNN

modelName = "SCDD"


class SCDD:
    def __init__(self, name=None, raw=None, label=None,
                 Tran=True, id=0, dropout=True, structure='AnnData', max_epoch=500,
                 method="TFIDF", neighbor_method="SC3", filter=True, format="tsv", conservative=True,
                 neighbors=20, threshold=0.2, batch_size=5000, slot=None):
        """
        :param name: the name of experiments or dataset to identify the imputed results. You can use the name
                    in our default experiments;
        :param raw: the raw data path, except format: `.tsv` and `.h5ad`, make sure that your delimiter in
                    .tsv file is `\t` and your raw data is .X in .h5ad file.
        :param label: the path to the annotation of datas, default none.
        :param Tran: if your data layout is: columns are cells and rows are genes, set this param `True`, else `False`.
        :param id: the id number to identify different results from the same methods, default None.
        :param dropout: whether need dropout in Denoising, default `True`;
        :param structure: the main structure of the pipline, default `AnnData`;
        :param max_epoch: the max epoch of Denoising, default 500.
        :param method: the method to generate graph embedding item, default `TFIDF`.
        :param neighbor_method: the method to generate concensus matrix and neighbor relations, default `SC3`. if the
                    number of cells too large, it will automatically turn to `SNN`.
        :param filter: whether need to filter genes when generating concensus mantrix, default `True`.
        :param format: the format of input raw data file, support `.tsv` and `.h5ad`
        :param conservative: whether needs conservative results, default `True`.
        :param neighbors: the neighbors of each cell, default 20.
        :param threshold: the threshold of whether the point is a drop-out point, the drop-out rate higher than this
                        value will be treated as a drop-out point, default 0.2.
        :param batch_size: the batch_size of each input in the nerual network when denoising, default 5000, which means
                        if cell number are less than 5000, the total batch will be 1.
        :param slot: for format `h5ad`, which slot contains raw data, default adata.X.
        """
        self.name = name
        self.raw = raw
        self.label = label
        self.Tran = Tran
        self.log_data = None
        self.data = None
        self.id = id
        self.dropout = dropout
        self.structure = structure
        self.filter = filter
        self.format = format
        self.method = method
        self.max_epoch = max_epoch
        self.neighbor_method = neighbor_method
        self.conservative = conservative
        self.neighbors = neighbors
        self.threshold = threshold
        self.batch_size = batch_size
        self.point = 1.01
        self.slot=slot
        if self.name == "Cellcycle":
            self.raw = "data/Cellcycle.raw.txt"
            self.label = "data/Cellcycle.label.txt"
            self.Tran = True
        elif self.name == "Li":
            self.raw = "data/Li.raw.tsv"
            self.label = "data/Li.label.txt"
            self.Tran = True
        elif self.name == "Bone":
            self.raw = "data/Bone.raw.tsv"
            self.label = "data/Bone.label.txt"
            self.Tran = True
        elif self.name == "Timecourse":
            self.raw = "data/Timecourse.raw.tsv"
            self.label = "data/Timecourse.label.txt"
            self.Tran = True
        elif self.name == "liver":
            self.raw = "data/liver.raw.tsv"
            self.label = "data/liver.label.txt"
            self.Tran = True
        elif self.name == "fish":
            self.raw = "data/fish.raw.tsv"
            self.label = ""
            self.Tran = True
        elif self.name == "CIDR":
            self.raw = "data/CIDR.raw.tsv"
            self.label = "data/CIDR.label.txt"
            self.Tran = True
        elif self.name == "sc_10x":
            self.raw = "data/sc_10x.raw.tsv"
            self.label = "data/sc_10x.label.txt"
            self.Tran = True
        elif self.name == "sc_celseq2":
            self.raw = "data/sc_celseq2.raw.tsv"
            self.label = "data/sc_celseq2.label.txt"
            self.Tran = True
        elif self.name == "sc_dropseq":
            self.raw = "data/sc_dropseq.raw.tsv"
            self.label = "data/sc_dropseq.label.txt"
            self.Tran = True

    def run(self, store=False):
        """
        Do Diffusion and Denoising, and save the results in ./results
        :param store: if `True`, it will load the temp data from ./temp. This options is for the running SCDD twice
        and more.
        """
        fm = self.raw.split('.')[-1]
        if self.format is None:
            self.format = fm    # automatically identify the format
            print("Identified raw data format: {0}".format(self.format))
        self.data, self.Label = LoadData(self.name, self.raw,
                                         labelPath=self.label,
                                         format=self.format,
                                         structure=self.structure,
                                         needTrans=self.Tran,
                                         slot=self.slot)
        if self.structure == "AnnData":
            if isinstance(self.data.X, np.ndarray) is False:
                self.log_data = self.data.X.toarray()
                self.log_data = np.log(self.log_data + self.point)
            else:
                self.log_data = np.log(self.data.X + self.point)
        else:
            self.log_data = np.log(self.data + self.point)
        print("Using neighbors:{0}.".format(self.neighbors))
        print("Using threshold:{0}.".format(self.threshold))
        print("Using max epoch:{0}.".format(self.max_epoch))
        print("Using batch size:{0}.".format(self.batch_size))

        # only support SNN if the number of cells are more than 10000
        if self.data.shape[0] >= 10000:
            self.neighbor_method = "SNN"
        if store:
            if self.neighbor_method == "SNN":
                M, Omega, Target, dropout_rate, null_genes = LoadTargets(True)
            else:
                M, Omega, Target, dropout_rate, null_genes = LoadTargets(False)
        else:
            if self.neighbor_method == "SC3":
                M = RobjKmeans(self.data)
                self.sparse = False
            elif self.neighbor_method == "SNN":
                M = AnnDataSNN(self.data)
                self.sparse = True
            dropout_rate, null_genes = get_dropout_rate(self.log_data)
            Omega, Target = get_supervise(self.log_data, dropout_rate, null_genes, M, self.neighbors, self.threshold,
                                          self.sparse)
            SaveTargets(M, Omega, Target, dropout_rate, null_genes, sparse=self.sparse)
        # whether need to sperate A to batch
        if self.batch_size >= M.shape[0]:
            if self.structure == "AnnData":
                A = getA(self.data.X.toarray(), method=self.method, filter=self.filter)
            else:
                A = getA(self.data, method=self.method, filter=self.filter)
        else:
            A = self.data.X.todense() > 0
        self.md = SC_Denoising(self.log_data, A, Omega, Target, batch_size=self.batch_size)
        for i in range(4):
            self.md.train(self.max_epoch)
            self.result = self.md.impute()
            self.result, self.cresult = makeup_results_all(self.result, self.log_data, null_genes, dropout_rate,
                                                           self.threshold)
            self.cresult = np.exp(self.cresult) - 1.01 + 0.5
            self.cresult = self.cresult.astype(np.int)
            # SaveData(self.name, modelName, self.result, format=self.format, id=self.id)
            SaveData(self.name, modelName, self.cresult, format=self.format, id=self.id + i)
        self.m_i = 6

    def continous_running(self, n):
        if self.neighbor_method == "SNN":
            M, Omega, Target, dropout_rate, null_genes = LoadTargets(True)
        else:
            M, Omega, Target, dropout_rate, null_genes = LoadTargets(False)
        self.md.train(n)
        self.result = self.md.impute()
        self.result, self.cresult = makeup_results_all(self.result, self.log_data, null_genes, dropout_rate,
                                                       self.threshold)
        self.cresult = np.exp(self.cresult) - 1.01 + 0.5
        self.cresult = self.cresult.astype(np.int)
        # SaveData(self.name, modelName, self.result, format=self.format, id=self.id)
        SaveData(self.name, modelName, self.cresult, format=self.format, id=self.id + self.m_i)
        self.m_i += 1

    def run_Diffusion(self, store=False):
        self.data, self.Label = LoadData(self.name, self.raw,
                                         labelPath=self.label,
                                         format=self.format,
                                         needTrans=self.Tran)
        self.log_data = np.log(self.data + 1.01)
        if store:
            M, Omega, Target, dropout_rate, null_genes = LoadTargets()
        else:
            M = RobjKmeans(self.data)
            dropout_rate, null_genes = get_dropout_rate(self.log_data)
            Omega, Target = get_supervise(self.log_data, dropout_rate, null_genes, M)
            SaveTargets(M, Omega, Target, dropout_rate, null_genes)
        self.result = Target
        self.result = makeup_results(self.result, self.log_data, null_genes, dropout_rate, self.threshold)
        self.result = np.exp(self.result) - 1
        self.result = np.round(self.result)
        self.result[self.result < 0] = 0
        SaveData(self.name, "Diffusion", self.result, format=self.format)
