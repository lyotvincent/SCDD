# This scripts generates all comparisons including MAGIC, 
# SAVER, DCA, SCDD and SCDD(Diffusion) for default experiments in paper
# All of the results will be saved in `results` fold
from utils.comps import *
from utils.purity import cal_purity
from SCDD_api import *
import os
import pandas as pd

def RunComps(expName, dataPath, labelPath=None, format="tsv"):
    dd = SCDD(expName, raw=dataPath, batch=0, format=format)
    dd.run_Diffusion()
    dd.run(store=True)
    # imputeByscGNN(expName, dataPath, labelPath, format)
    # # imputeByscImpute(expName, dataPath, labelPath, format)
    # # imputeByscTSSR(expName, dataPath, labelPath, format)
    imputeByDeepImpute(expName, dataPath, labelPath, format)
    # imputeByscIGANs(expName, dataPath, labelPath, format)
    imputeByDrImpute(expName, dataPath, labelPath, format)
    imputeByVIPER(expName, dataPath, labelPath, format)
    imputeByMAGIC(expName, dataPath, labelPath, format)
    imputeBySAVER(expName, dataPath, labelPath, format)
    imputeByDCA(expName, dataPath, labelPath)
    imputeByALRA(expName, dataPath, labelPath, format)


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

def Generate_CIDR(has_results=False):
    if has_results == False:
        expName = "CIDR"
        dataPath = "data/CIDR.raw.tsv"
        labelPath = "data/CIDR.label.txt"
        format = "tsv"
        RunComps(expName, dataPath, labelPath, format=format)

def Generate_sc_10x(has_results=False):
    if has_results == False:
        expName = "sc_10x"
        dataPath = "data/sc_10x.raw.tsv"
        labelPath = "data/sc_10x.label.txt"
        format = "tsv"
        RunComps(expName, dataPath, labelPath, format=format)

def Generate_sc_celseq2(has_results=False):
    if has_results == False:
        expName = "sc_celseq2"
        dataPath = "data/sc_celseq2.raw.tsv"
        labelPath = "data/sc_celseq2.label.txt"
        format = "tsv"
        RunComps(expName, dataPath, labelPath, format=format)

def Generate_sc_dropseq(has_results=False):
    if has_results == False:
        expName = "sc_dropseq"
        dataPath = "data/sc_dropseq.raw.tsv"
        labelPath = "data/sc_dropseq.label.txt"
        format = "tsv"
        RunComps(expName, dataPath, labelPath, format=format)

# Generate_CIDR(has_results=False)
Generate_sc_10x(has_results=False)
Generate_sc_celseq2(has_results=False)
Generate_sc_dropseq(has_results=False)
