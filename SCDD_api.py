import numpy as np
from utils.io import *
from utils.kernel import *
from models.scdd import *
from utils.prepropcess import get_dropout_rate, get_supervise, makeup_results, makeup_results_all, RobjKmeans
modelName = "SCDD"

class SCDD:
    def __init__(self, name=None, raw=None, label=None,
                 Tran=True, batch=None, dropout=True,
                 method="TFIDF", filter=True, format="tsv", conservative=False):
        self.name = name
        self.raw = raw
        self.label = label
        self.Tran = Tran
        self.log_data = None
        self.data = None
        self.batch = batch
        self.dropout = dropout
        self.filter = filter
        self.format = format
        self.method = method
        self.conservative = conservative
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
        elif self.name == "guo":
            self.raw = "data/guo.raw.tsv"
            self.label = "data/guo.label.txt"
            self.Tran = True

    def run(self, store = False):
        self.data, self.Label = LoadData(self.name, self.raw,
                                         labelPath=self.label,
                                         format=self.format,
                                         needTrans=self.Tran)
        self.log_data = np.log(self.data + 1.01)
        A = getA(self.data, method=self.method, filter=self.filter)
        if store:
            M, Omega, Target, dropout_rate, null_genes = LoadTargets()
        else:
            M = RobjKmeans(self.data)
            dropout_rate, null_genes = get_dropout_rate(self.log_data)
            Omega, Target = get_supervise(self.log_data , dropout_rate, null_genes, M)
            SaveTargets(M, Omega, Target, dropout_rate, null_genes)
        md = SC_Denoising(self.log_data, A, Omega, Target)
        md.train(200)
        self.result = md.impute()
        self.result, self.cresult = makeup_results_all(self.result, self.log_data, null_genes, dropout_rate)
        self.result = np.exp(self.result) - 1.01 + 0.5
        self.result = self.result.astype(np.int)
        self.cresult = np.exp(self.cresult) - 1.01 + 0.5
        self.cresult = self.cresult.astype(np.int)
        SaveData(self.name, modelName, self.result, format=self.format, batch=self.batch)
        SaveData(self.name, modelName, self.cresult, format=self.format, batch=self.batch+1)
    
    def run_Diffusion(self, store = False):
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
            Omega, Target = get_supervise(self.log_data , dropout_rate, null_genes, M)
            SaveTargets(M, Omega, Target, dropout_rate, null_genes)
        self.result = Target
        self.result = makeup_results(self.result, self.log_data, null_genes, dropout_rate)
        self.result = np.exp(self.result) - 1
        self.result = np.round(self.result)
        self.result[self.result < 0] = 0
        SaveData(self.name, "Diffusion", self.result, format=self.format)
