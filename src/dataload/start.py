import sys, os
import os.path
import dataload
import time
import random
import importlib

src_path = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
sys.path.append(src_path)


def main(source):
    '''
    Example:
        python -m dataload/start ensembl
        python -m dataload/start entrez
        python -m dataload/start pharmgkb

    '''
    # discover sources
    try:
        src_mod = importlib.import_module("dataload.sources.%s.upload" % source)
        mods = []
        for pyfile in os.listdir(src_mod.__path__[0]):
            if pyfile.startswith("_"):
                continue
            if not pyfile.endswith(".py"):
                continue
            mods.append("%s.upload.%s" % (source, pyfile.replace(".py","")))

        # alphanum sorted, in case sub-sources have dependencies
        dataload.__sources__ = sorted(mods)
        dataload.register_sources()
        dataload.load_all()

    except ImportError as e:
        raise ValueError('Unknown source "%s" (%s)' % (source,e))


def main_test(src):
    t0 = time.time()
    i = 0
    limit = 0 + random.randint(0, 50)
    while True:
        random.random() * random.random()
        j = int(round(time.time() - t0, 0))
        if j > i:
            if j > limit:
                break
            else:
                print((src, j))
                i = j

if __name__ == '__main__':
    main(sys.argv[1])
