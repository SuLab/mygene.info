'''
A thin python layer for accessing GeneDoc ElasticSearch host.

Currently available URLs:

    /query?q=cdk2      gene query service
    /gene/<geneid>     gene annotation service

'''
import sys

#import os.path
import os, logging
import subprocess
#import json

import tornado.options
import tornado.web
from tornado.options import options


#TODO: move to a "common api" module (not specific to an API version)
#*******#
from biothings.www.api.handlers import StatusHandler, BiothingHandler
from utils.es import ESQuery
class MyGeneStatusHandler(StatusHandler):
    ''' This class is for the /status endpoint. '''
    esq = ESQuery()

class MainHandler(BiothingHandler):
    def get(self):
        if INCLUDE_DOCS:
            self.render(os.path.join(DOCS_STATIC_PATH, 'index.html'))

DEMO_STATIC_FILE = '../docs/demo/index.html'
class DemoHandler(BiothingHandler):
    def get(self):
        with open(DEMO_STATIC_FILE,'r') as demo_file:
            self.write(demo_file.read())
#*******#

def _get_rev():
    '''return current git rev number.
        e.g.
           9d9811e21a878bed4c8d7de5ab20bf0baa5aa87f
    '''
    pipe = subprocess.Popen(["git", "rev-parse", "HEAD"],
                            stdout=subprocess.PIPE)
    output = pipe.stdout.read().strip().decode()
    return ':'.join(reversed(output.replace('+', '').split(' ')))
__revision__ = _get_rev()
os.environ["MYGENE_REVISION"] = __revision__


from biothings.www.index_base import main, settings, get_app
from biothings.www.helper import add_apps
from biothings.settings import BiothingSettings
btsettings = BiothingSettings()

# build API routes
from www.api.handlers import APP_LIST as api_v2_app_list
#from api.handlers_async import APP_LIST as api_v2_async_app_list
from demo.handlers import APP_LIST as demo_app_list
#from auth.handlers import APP_LIST as auth_app_list
from www.api.handlers import MyGeneMetaDataHandler
from www.api.handlers import MyGeneFieldsHandler


from config import INCLUDE_DOCS
__USE_WSGI__ = False
DOCS_STATIC_PATH = '../docs/_build/html'
if INCLUDE_DOCS and not os.path.exists(DOCS_STATIC_PATH):
    raise IOError('Run "make html" to generate sphinx docs first.')


APP_LIST = [
    (r"/", MainHandler),
    (r"/status", MyGeneStatusHandler),
    (r"/metadata", MyGeneMetaDataHandler),
    #TODO: what is v2a ?
    (r"/v2a/metadata", MyGeneMetaDataHandler),
    (r"/metadata/fields", MyGeneFieldsHandler),
    (r"/demo/?$", DemoHandler),
]
APP_LIST += add_apps('', api_v2_app_list)
APP_LIST += add_apps('v2', api_v2_app_list)
#APP_LIST += add_apps('demo', demo_app_list)
#APP_LIST += add_apps('v2a', api_v2_async_app_list)
#APP_LIST += add_apps('auth', auth_app_list)

if options.debug:
    APP_LIST += [
        #/widget/* static path
        (r"/widget/(.*)", tornado.web.StaticFileHandler, {'path': os.path.join(btsettings.static_path, 'widget')}),
        #this should be the last one
        (r"/?(.*)", tornado.web.StaticFileHandler, {'path': DOCS_STATIC_PATH}),
    ]


# TODO: in config ? (prod = wsgi app)
if __USE_WSGI__:
    import tornado.wsgi
    wsgi_app = tornado.wsgi.WSGIApplication(APP_LIST)


if __name__ == '__main__':
    main(APP_LIST)



