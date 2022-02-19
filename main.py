from SCDD_api import *
from generates import *
if __name__ == '__main__':
    e = SCDD("Timecourse", batch=224, dropout=True, method="TFIDF")
    e.run(False)
    Generate_Cellcycle()

