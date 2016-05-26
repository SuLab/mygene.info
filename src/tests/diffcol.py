import sys
import utils.diff as diff
from databuild.backend import GeneDocMongoDBBackend
from pymongo import MongoClient


def main(s1, s2, p=False):
    conn1, dbstr1, colstr1 = s1
    conn2, dbstr2, colstr2 = s2
    col1 = MongoClient(conn1)[dbstr1][colstr1]
    col2 = MongoClient(conn2)[dbstr2][colstr2]

    b1 = GeneDocMongoDBBackend(col1)
    b2 = GeneDocMongoDBBackend(col2)

    changes = diff.diff_collections(b1, b2, use_parallel=p, step=1000)

    return changes


if __name__ == "__main__":
    conn1 = sys.argv[1]
    dbstr1, colstr1 = sys.argv[2].split("/")
    conn2 = sys.argv[3]
    dbstr2, colstr2 = sys.argv[4].split("/")
    p = "-p" in sys.argv

    res = main([conn1, dbstr1, colstr1], [conn2, dbstr2, colstr2], p)
    import pickle
    pickle.dump(res, open("zeres%s" % colstr1, "wb"))
    print("done see pickle")

