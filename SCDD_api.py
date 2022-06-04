import numpy as np
from utils.io import *
from utils.kernel import *
from models.scdd import *
from utils.prepropcess import get_dropout_rate, get_supervise, makeup_results, makeup_results_all, RobjKmeans
modelName = "SCDD"

class SCDD:
    def __init__(self, name=None, raw=None, label=None,
                 Tran=True, id=None, dropout=True,
                 method="TFIDF", filter=True, format="tsv", conservative=False,
                 neighbors=20, threshold=0.2, batch_size=1024):
        self.name = name
        self.raw = raw
        self.label = label
        self.Tran = Tran
        self.log_data = None
        self.data = None
        self.id = id
        self.dropout = dropout
        self.filter = filter
        self.format = format
        self.method = method
        self.conservative = conservative
        self.neighbors = neighbors
        self.threshold = threshold
        self.batch_size = batch_size
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

    def run(self, store = False):
        self.data, self.Label = LoadData(self.name, self.raw,
                                         labelPath=self.label,
                                         format=self.format,
                                         needTrans=self.Tran)
        self.log_data = np.log(self.data + 1.01)
        A = getA(self.data, method=self.method, filter=self.filter)
        print("Using neighbors:{0}.".format(self.neighbors))
        print("Using threshold:{0}.".format(self.threshold))
        if store:
            M, Omega, Target, dropout_rate, null_genes = LoadTargets()
        else:
            M = RobjKmeans(self.data)
            dropout_rate, null_genes = get_dropout_rate(self.log_data)
            Omega, Target = get_supervise(self.log_data , dropout_rate, null_genes, M, self.neighbors, self.threshold)
            SaveTargets(M, Omega, Target, dropout_rate, null_genes)
        md = SC_Denoising(self.log_data, A, Omega, Target, batch_size=self.batch_size)
        md.train(2000)
        self.result = md.impute()
        self.result, self.cresult = makeup_results_all(self.result, self.log_data, null_genes, dropout_rate, self.threshold)
        self.result = np.exp(self.result) - 1.01 + 0.5
        self.result = self.result.astype(np.int)
        self.cresult = np.exp(self.cresult) - 1.01 + 0.5
        self.cresult = self.cresult.astype(np.int)
        SaveData(self.name, modelName, self.result, format=self.format, id=self.id)
        SaveData(self.name, modelName, self.cresult, format=self.format, id=self.id+1)
    
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
        self.result = makeup_results(self.result, self.log_data, null_genes, dropout_rate, self.threshold)
        self.result = np.exp(self.result) - 1
        self.result = np.round(self.result)
        self.result[self.result < 0] = 0
        SaveData(self.name, "Diffusion", self.result, format=self.format)
