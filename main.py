from SCDD_api import *

if __name__ == '__main__':
    e = SCDD(raw="data/Bladder/chunk5000_part0.h5ad", format="h5ad")
    e.run(store=False)