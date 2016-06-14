import time
from optparse import OptionParser

from databuild.builder import DataBuilder, timesofar
# from utils.es import es_clean_indices
from utils.common import ask
from .tunnel import open_tunnel, es_local_tunnel_port
from config import ES_HOST


def validate(build_config=None):
    from pprint import pprint
    from utils.diff import diff_collections
    from databuild.backend import GeneDocMongoDBBackend, GeneDocESBackend
    from utils.mongo import get_src_build, get_target_db
    from utils.es import ESIndexer

    src_build = get_src_build()
    _cfg = src_build.find_one({'_id': build_config})
    last_build = _cfg['build'][-1]
    print("Last build record:")
    pprint(last_build)
    target_name = last_build['target']

    mongo_target = get_target_db()
    b1 = GeneDocMongoDBBackend(mongo_target[target_name])
    b2 = GeneDocESBackend(ESIndexer(es_index_name=target_name,
                                    es_host='127.0.0.1:' + str(es_local_tunnel_port)))
    changes = diff_collections(b1, b2, use_parallel=True, step=10000)
    return changes


def test(build_config):
    with open_tunnel():
        validate(build_config)


def main():
    parser = OptionParser()
    parser.add_option("-c", "--conf", dest="config",
                      action="store", default=None,
                      help="ES indexing building config name")
    parser.add_option("-b", "--noconfirm", dest="noconfirm",
                      action="store_true", default=False,
                      help="do not ask for confirmation")
    parser.add_option("-e", "--es-index", dest="es_index_name",
                      action="store", default=None,
                      help="provide an alternative ES index name")
    # parser.add_option("", "--no-cleanup", dest="nocleanup",
    #                   action="store_true", default=False,
    #                   help="do not clean up old ES indices")
    (options, args) = parser.parse_args()

    with open_tunnel() as tunnel:
        if tunnel.ok:
            es_host = '127.0.0.1:' + str(es_local_tunnel_port)
        else:
            es_host = ES_HOST
        # if not options.nocleanup:
        #     es_clean_indices(noconfirm=options.noconfirm)
        t00 = time.time()
        bdr = DataBuilder(backend='es')
        if options.config:
            config_li = [options.config]
        else:
            config_li = ['mygene', 'mygene_allspecies']

        if not options.noconfirm:
            print('\n'.join(["Ready to build these ES indices on %s (tunnel=%s):" % (es_host, tunnel.ok)] +
                            ['\t' + conf for conf in config_li]))
            if ask('Continue?') != 'Y':
                print("Aborted")
                return

        for _conf in config_li:
            t0 = time.time()
            print('>"{}">>>>>>'.format(_conf))
            bdr.build_index2(_conf,
                             es_index_name=options.es_index_name,
                             es_host=es_host,noconfirm=options.noconfirm)
            print('<<<<<<"{}"...done. {}'.format(_conf, timesofar(t0)))

        print('=' * 20)
        print("Finished.", timesofar(t00))


if __name__ == '__main__':
    main()
