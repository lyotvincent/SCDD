# This scripts generates all comparisons including MAGIC, 
# SAVER, DCA, SCDD and SCDD(Diffusion) for default experiments in paper
# All of the results will be saved in `results` fold
from utils.comps import *
from utils.purity import cal_purity
from SCDD_api import *
import os
import pandas as pd

def RunComps(expName, dataPath, labelPath=None):
    imputeByMAGIC(expName, dataPath, labelPath)
    imputeBySAVER(expName, dataPath, labelPath)
    imputeByDCA(expName, dataPath, labelPath)
    dd = SCDD(expName, batch=0)
    dd.run_Diffusion()
    dd.run(store=True)


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
    predicts = pd.read_csv("temp/Bone_predicts.tsv", sep='\t')
    perfs = pd.read_csv("temp/Bone_perfs.tsv", sep='\t')
    mcounts = len(perfs) - 1
    perfs = cal_purity(perfs, predicts, mcounts)
    perfs.to_csv("paper/Bone/Bone_perfs.tsv", sep='\t')
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
    predicts = pd.read_csv("temp/Li_predicts.tsv", sep='\t')
    perfs = pd.read_csv("temp/Li_perfs.tsv", sep='\t')
    mcounts = len(perfs) - 1
    perfs = cal_purity(perfs, predicts, mcounts)
    perfs.to_csv("paper/Li/Li_perfs.tsv", sep='\t')
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
        dataPath = "data/guo.raw.txt"
        labelPath = "data/guo.label.txt"
        RunComps(expName, dataPath, labelPath)
    os.system("~/anaconda3/bin/python paper/guo/guo.py")

Generate_GUO(True)
