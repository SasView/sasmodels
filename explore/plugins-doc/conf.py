import os
import sys

# Location of the base conf.py
#_base_conf_dir = '../../../sasview/docs/sphinx-docs/source-temp'
_base_conf_dir = '../../doc'  # from sasmodels/explore/plugins-doc
sys.path.insert(0, os.path.abspath(_base_conf_dir))
from conf import *

# overriding theme so I can see it is different
#html_theme = 'haiku'

# Add the intersphinx extensions
extensions = [
    *extensions,
    'sphinx.ext.intersphinx',
    'sasmodels.sphinx.remoteimage',
    'sasmodels.sphinx.intersphinx_nav',
    # 'sphinx_external_toc',
]

html_static_path = [ f'{_base_conf_dir}/_static']

# Define the external document location
_remote_root = 'https://www.sasview.org/docs'
#_remote_root = 'file:///Users/pkienzle/Source/sasview/build/doc/html'
intersphinx_mapping = {
    # requests library fails for file:// schema, so strip it from the URI
    'sasview': (_remote_root, f'{_remote_root.replace("file://","")}/objects.inv'),
}
remote_image_url = f'{_remote_root}/_images'


# Place the plugin at the following point in the SasView documentation:
#
#     sasview > user-guide > models > plugins
#
# Need to specify the siblings of the parents to resolve next/prev links:
#
#     sasview: ... user-guide ...
#     user: ... models menu-bar ...
#     models: ... structure-factor plugins
#
_user_docs = "user/user"
_models_root = "user/qtgui/Perspectives/Fitting/models"
_models_page = f"{_models_root}/index"
_structure_factors = f"{_models_root}/structure-factor"
remote_relations = {
    'index': [_user_docs],
    _user_docs: [_models_page, 'user/menu_bar'],
    _models_page: [_structure_factors, 'model/plugins'],
}
