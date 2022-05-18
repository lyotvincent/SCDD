from SCDD_api import *

if __name__ == '__main__':
    DataSplit("data/TS_Bladder.h5ad", "data/Bladder", chunksize=5000)