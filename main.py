from SCDD_api import *

if __name__ == '__main__':
    e = SCDD("guo", batch=4, dropout=True, method="TFIDF")
    e.run(True)