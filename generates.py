# This scripts generates all comparisons including MAGIC, 
# SAVER, DCA, SCDD and SCDD(Diffusion) for default experiments in paper
# All of the results will be saved in `results` fold
from utils.comps import *
from utils.purity import cal_purity
from SCDD_api import *
import os
import pandas as pd

def RunComps(expName, dataPath, labelPath=None, format="tsv", id=0):
    dd = SCDD(expName, raw=dataPath, id=0, format=format)
    dd.run_Diffusion()
    dd.run(store=True)
    imputeByscGNN(expName, dataPath, labelPath, format, id)
    # # # imputeByscImpute(expName, dataPath, labelPath, format)
    # # # imputeByscTSSR(expName, dataPath, labelPath, format)
    imputeByDeepImpute(expName, dataPath, labelPath, format, id)
    imputeByscIGANs(expName, dataPath, labelPath, format, id)
    imputeByDrImpute(expName, dataPath, labelPath, format, id)
    imputeByVIPER(expName, dataPath, labelPath, format, id)
    imputeByMAGIC(expName, dataPath, labelPath, format)
    imputeByDCA(expName, dataPath, labelPath, format)
    imputeByALRA(expName, dataPath, labelPath, format, id)
    imputeBySAVER(expName, dataPath, labelPath, format)


def Generate_all_imputations():
    # Cellcycle
    expName = "Cellcycle"
    dataPath = "data/Cellcycle.raw.txt"
    labelPath = "data/Cellcycle.label.txt"
    RunComps(expName, dataPath, labelPath)

    # Li
    expName = "Li"
    dataPath = "data/Li.raw.tsv"
    labelPath = "data/Li.label.txt"
    RunComps(expName, dataPath, labelPath)

    # Timecourse
    expName = "Timecourse"
    dataPath = "data/Timecourse.raw.tsv"
    labelPath = "data/Timecourse.label.txt"
    RunComps(expName, dataPath, labelPath)

    # Bone
    expName = "Bone"
    dataPath = "data/Bone.raw.tsv"
    labelPath = "data/Bone.label.txt"
    RunComps(expName, dataPath, labelPath)


def Generate_Bone(has_results=False):
    if has_results == False:
        expName = "Bone"
        dataPath = "data/Bone.raw.tsv"
        labelPath = "data/Bone.label.txt"
        RunComps(expName, dataPath, labelPath)
    os.system("Rscript paper/Bone/Bone.R")
    predicts = pd.read_csv("temp/Bone_predicts1.tsv", sep='\t')
    perfs = pd.read_csv("temp/Bone_perfs1.tsv", sep='\t')
    mcounts = len(perfs) - 1
    perfs = cal_purity(perfs, predicts, mcounts)
    perfs.to_csv("paper/Bone/Bone_perfs1.tsv", sep='\t')
    print("Generate Bone OK.")


def Generate_Cellcycle(has_results=False):
    if has_results == False:
        # Cellcycle
        expName = "Cellcycle"
        dataPath = "data/Cellcycle.raw.txt"
        labelPath = "data/Cellcycle.label.txt"
        RunComps(expName, dataPath, labelPath)
    os.system("Rscript paper/Cellcycle/Cellcycle.R")
    print("Generate Cellcycle OK.")

def Generate_Li(has_results=False):
    if has_results == False:
        # Li
        expName = "Li"
        dataPath = "data/Li.raw.tsv"
        labelPath = "data/Li.label.txt"
        RunComps(expName, dataPath, labelPath)
    os.system("Rscript paper/Li/Li.R")
    predicts = pd.read_csv("temp/Li_predicts1.tsv", sep='\t')
    perfs = pd.read_csv("temp/Li_perfs1.tsv", sep='\t')
    mcounts = len(perfs) - 1
    perfs = cal_purity(perfs, predicts, mcounts)
    perfs.to_csv("paper/Li/Li_perfs1.tsv", sep='\t')
    print("Generate Li OK.")

def Generate_Timecourse(has_results=False):
    if has_results == False:
        # Timecourse
        expName = "Timecourse"
        dataPath = "data/Timecourse.raw.tsv"
        labelPath = "data/Timecourse.label.txt"
        RunComps(expName, dataPath, labelPath)
    os.system("~/anaconda3/bin/python paper/Timecourse/Timecourse.py")

def Generate_GUO(has_results=False):
    if has_results == False:
        # Timecourse
        expName = "guo"
        dataPath = "data/guo.raw.tsv"
        labelPath = "data/guo.label.txt"
        RunComps(expName, dataPath, labelPath)
    os.system("~/anaconda3/bin/python paper/guo/guo.py")

def Generate_liver(has_results=False):
    if has_results == False:
        expName = "liver"
        dataPath = "data/liver.raw.tsv"
        labelPath = "data/liver.label.txt"
        format = "tsv"
        RunComps(expName, dataPath, labelPath=labelPath, format=format)
    os.system("~/anaconda3/bin/python paper/liver/liver.py")

def Generate_fish(has_results=False):
    if has_results == False:
        expName = "fish"
        dataPath = "data/fish.raw.tsv"
        format = "tsv"
        RunComps(expName, dataPath, labelPath=None, format=format)
        # SCDD:150mins
        # scIGANs: 2.5days
        # scGNN: 192mins
        # deepimpute: 10min

def Generate_CIDR(has_results=False):
    if has_results == False:
        expName = "CIDR"
        dataPath = "data/CIDR.raw.tsv"
        labelPath = "data/CIDR.label.txt"
        format = "tsv"
        RunComps(expName, dataPath, labelPath, format=format)
    os.system("Rscript paper/CIDR/CIDR.R")
    predicts = pd.read_csv("temp/CIDR_predicts1.tsv", sep='\t')
    perfs = pd.read_csv("temp/CIDR_perfs1.tsv", sep='\t')
    mcounts = len(perfs) - 1
    perfs = cal_purity(perfs, predicts, mcounts)
    perfs.to_csv("paper/CIDR/CIDR_perfs1.tsv", sep='\t')
    print("Generate CIDR OK.")

def Generate_sc_10x(has_results=False):
    if has_results == False:
        expName = "sc_10x"
        dataPath = "data/sc_10x.raw.tsv"
        labelPath = "data/sc_10x.label.txt"
        format = "tsv"
        RunComps(expName, dataPath, labelPath, format=format)
    os.system("Rscript paper/sc_10x/sc_10x.R")
    predicts = pd.read_csv("temp/sc_10x_predicts1.tsv", sep='\t')
    perfs = pd.read_csv("temp/sc_10x_perfs1.tsv", sep='\t')
    mcounts = len(perfs) - 1
    perfs = cal_purity(perfs, predicts, mcounts)
    perfs.to_csv("paper/sc_10x/sc_10x_perfs1.tsv", sep='\t')
    print("Generate Li OK.")

def Generate_sc_celseq2(has_results=False):
    if has_results == False:
        expName = "sc_celseq2"
        dataPath = "data/sc_celseq2.raw.tsv"
        labelPath = "data/sc_celseq2.label.txt"
        format = "tsv"
        RunComps(expName, dataPath, labelPath, format=format)
    os.system("Rscript paper/sc_celseq2/sc_celseq2.R")
    predicts = pd.read_csv("temp/sc_celseq2_predicts1.tsv", sep='\t')
    perfs = pd.read_csv("temp/sc_celseq2_perfs1.tsv", sep='\t')
    mcounts = len(perfs) - 1
    perfs = cal_purity(perfs, predicts, mcounts)
    perfs.to_csv("paper/sc_celseq2/sc_celseq2_perfs1.tsv", sep='\t')
    print("Generate Li OK.")

def Generate_sc_dropseq(has_results=False):
    if has_results == False:
        expName = "sc_dropseq"
        dataPath = "data/sc_dropseq.raw.tsv"
        labelPath = "data/sc_dropseq.label.txt"
        format = "tsv"
        RunComps(expName, dataPath, labelPath, format=format)
    os.system("Rscript paper/sc_dropseq/sc_dropseq.R")
    predicts = pd.read_csv("temp/sc_dropseq_predicts1.tsv", sep='\t')
    perfs = pd.read_csv("temp/sc_dropseq_perfs1.tsv", sep='\t')
    mcounts = len(perfs) - 1
    perfs = cal_purity(perfs, predicts, mcounts)
    perfs.to_csv("paper/sc_dropseq/sc_dropseq_perfs1.tsv", sep='\t')
    print("Generate Li OK.")

def Generate_Bladder(has_results=False):
    if has_results == False:
        expName = "sc_dropseq"
        dataPath = "data/TS_Bladder.h5ad"
        labelPath = ""
        format = "h5ad"
        RunComps(expName, dataPath, labelPath, format=format)
    print("Generate Bladder OK.")

# Generate_CIDR(has_results=False)
# Generate_sc_10x(has_results=False)
Generate_Timecourse(has_results=True)
# Generate_sc_dropseq(has_results=True)
# imputeByDCA("immune", "data/TS_immune.h5ad", format='h5ad')