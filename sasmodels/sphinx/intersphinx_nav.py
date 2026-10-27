"""
======================
Intersphinx navigation
======================

.. role:: code-py(code)
  :language: Python

If you are making a document which virtually resides in an existing remote document
then you can set up an intersphinx navigation plugin to tie the root of your document
to a point in the remote document.

Your document will operate as normal, except that the html navbar will contain the
all the ancestors up to the root of the remote document.  Your document will be
invisible to the remote document. You can use previous and next to navigate around
your tree, but as soon as you leave for the remote you will not be able to get back.

Configuration
-------------

To use Intersphinx navigation, add ``'sasmodels.sphinx.intersphinx_nav'`` to your
*extensions* config value, and use these config values to activate
linking:

**remote_relations**

.. container:: confval

    **Type**: :code-py:`dict[str, list[str]]`

    **Default**: :code-py:`{{}}`

    A toctree dictionary giving all parents keyed by all parents
    of the insertion point, plus the parents of the page that comes next after the
    final page of local toctree.

    For each local toctree we need to include its parent documents all the way to
    the root (these are displayed on the navigation bar). As well we need enough
    information to determine the next and previous document, which is the "uncle"
    pages directly preceding and following the parent.

    Be sure to include the target master_doc somewhere in the toctree.

    **Example**

    ::

        # The plugins document defined by plugins/index.rst is located here
        #
        #     sasview > users guide > models > plugins
        #
        # The document before the first plugin is the user/models/ellipsoids page.
        # The document afterward the last plugin is user/menubar. The ... below are
        # remote documents that we don't care about.
        #
        #     sasview: ... user ...
        #     user: ... user/models user/menubar ...
        #     user/models: ... user/models/ellipsoids plugins
        #
        # If we wanted to put plugins at the start of user/models rather than the
        # end, then we wouldn't need to mention user/menubar. If in addition,
        # navigation were set up to skip toc pages, then we would need the preceding
        # sibling to user model and additional entries tracing the path to the final
        # child in the prior tree.
        #
        remote_toctree_includes = {
            'index': ['user/index'],
            'user/index': ['user/models/index', 'user/menu_bar'],
            'user/models/index': ['user/models/structure_factors', 'plugins/index'],
        }

When setting up intersphinx using a local build of the external sphinx document,
you need to use ``file://`` as the external target, but you must strip ``file://``
from the inventory location. For example:

.. code:: python

    _remote_root = file:///path/to/external
    intersphinx_mapping = {
        'external': (_remote_root, f'{_remote_root.replace("file://","")}/objects.inv'),
    }

Notes
-----

This extension monkeypatches ``sphinx.util.osutil.relative_uri`` (returning the
target directly if it is a fully qualified URI) and
``shpinx.builders.html.StandaloneHTMLBuilder.get_target_uri`` (resolving to the
remote link if the target docname is in the intersphinx inventory). It adds
``intersphinx_pages`` to the the builder attributes on the ``env-updated`` event.
"""

# If this were in its own package I would use the sphinx_toolbox.confval extension
# from pip sphinx-toolbox, but it pulls in a lot of stuff. Another option is to
# vendor in the extension from:
#
#   https://github.com/sphinx-toolbox/sphinx-toolbox/blob/master/sphinx_toolbox/confval.py
#
#
# .. confval:: remote_relations
#     :type: :code-py:`dict[str, list[str]]`
#    :default: :code-py:`{{}}`
#

import contextlib
from urllib.parse import quote

import docutils.nodes
from sphinx.application import Sphinx
from sphinx.locale import _
from sphinx.util import logging
from sphinx.util.inventory import _InventoryItem

# TODO: replace hardcoded name with __name__ when this moves to its own module
logger = logging.getLogger(__name__)

def inject_remote_relations(app, env):
    """
    Inject remote parents into the toctree structure before writing.
    """
    remote_relations = env.config.remote_relations
    app.env.toctree_includes.update(remote_relations)
    # print(f"inject parents: {env.tocs=}\n{env.toc_num_entries=}\n{env.toctree_includes=}\n{env.current_document=}")

    # print(f"^^^^ inject remote {env.titles=}")

    # Intersphinx remote inventory of remote links (mostly std:doc but some std:label)
    # We may be missing some...
    intersphinx_pages = _get_intersphinx_pages(app)

    # Set the titles for each remote document
    # Note: as of the env-updated event, the env.titles structure appears to be empty
    titles = {}
    for parent, children in remote_relations.items():
        for docname in (parent, *children):
            if docname not in titles:
                if docname in intersphinx_pages:
                    titles[docname] = _build_title(intersphinx_pages[docname].display_name)
                elif docname not in app.env.titles:
                    logger.warning(f"missing {docname} in config.remote_relations for {__name__}")
    app.env.titles.update(titles)

    # Save the intersphinx pages in the build attributes for later link hijacking
    app.env.intersphinx_pages = intersphinx_pages

def _get_intersphinx_pages(app: Sphinx) -> dict[str, _InventoryItem]:
    """
    Retrieve the inventory from the intersphinx extension.
    """
    # TODO: maybe need to include py:modules links, etc.
    # Get inventory of remote pages from intersphinx
    inventory = getattr(app.env, 'intersphinx_inventory', {})
    pages = inventory.get("std:doc", {}).copy()
    # Include genindex from the std:label group.
    labels = inventory.get("std:label", {})
    if 'genindex' in labels:
        pages['genindex'] = labels['genindex']
    return pages


def _build_title(title: str) -> docutils.nodes.title:
    """
    Compose the title node needed to label to relations.
    """
    node = docutils.nodes.title()
    node += docutils.nodes.Text(title)
    return node


_patched = False
def _monkeypatch_link_hijack():
    from sphinx import builders
    from sphinx.builders.html import StandaloneHTMLBuilder
    from sphinx.util import osutil

    global _patched
    if _patched:
        return
    _patched = True

    get_target_uri_orig = StandaloneHTMLBuilder.get_target_uri
    def get_target_uri(self, docname: str, typ: str | None = None) -> str:
        # Intersphinx link hijacking
        remote = getattr(self.app.env, 'intersphinx_pages', {})
        if docname in remote:
            # print(f"   {docname} is remote")
            return remote[docname].uri
        return get_target_uri_orig(self, docname, typ)
        #return quote(docname) + self.link_suffix
    StandaloneHTMLBuilder.get_target_uri = get_target_uri

    # Keep existing relative_uri so that we can fallback to it for local checks
    # relative_uri was already loaded into the html builders, so replace it there
    relative_uri_orig = osutil.relative_uri
    def relative_uri(base: str, to: str) -> str:
        if '://' in to:
            # print(f"   ... target {to} is absolute")
            return to
        # print(f"   ... target {to} is relative")
        return relative_uri_orig(base, to)
    osutil.relative_uri = relative_uri
    builders.relative_uri = relative_uri
    builders.html.relative_uri = relative_uri


# *** Unused code ***
# If we don't monkeypatch link hijacking in the html builder then we need to
# duplicate a lot of code that referencess relative_uri and get_target_uri,
# either directly or indirectly. In this case, the monkeypatch is less brittle.
def add_remote_nav_context(app: Sphinx, pagename, templatename, context, doctree):
    from sphinx.util.osutil import relative_uri
    # TODO: do nothing if intersphinx inventory is not available

    # Mimic original code as much as possible by using the same variable names.
    # - The self variable points to the builder instance.
    # - StandaloneHTMLBuilder.handle_page uses pagename
    # - StandaloneHTMLBuilder.get_doc_context uses docname
    self = app.builder
    docname = pagename

    # Note: Subclassing builder would almost be good enough, with an additional
    # intersphinx_pages attribute added in the inject_remote_relations so that
    # get_target_uri() could resolve remote links. It won't work unless we also
    # monkeypatch relative_uri() to return 'to' if it is a fully qualified name.

    # Cribbed from StandaloneHTMLBuilder.get_target_uri
    def get_target_uri(docname: str, typ: str | None = None) -> str:
        """
        Replacement StandaloneHTMLBuilder.get_target_uri, with link hijacking for
        remote pages.
        """
        # Intersphinx link hijacking
        if docname in self.intersphinx_pages:
            return self.intersphinx_pages[docname].uri
        # Original code
        return quote(docname) + self.link_suffix

    # Cribbed from Builder.get_relative_uri
    def get_relative_uri(from_: str, to: str, typ: str | None = None) -> str:
        """
        Replacement Builder.get_relative_uri, with absolute links for remote targets.
        """
        from_= get_target_uri(from_)
        to = get_target_uri(to, typ)
        return to if '://' in to else relative_uri(from_, to)

        # Alternative specific to intershinx link hijacking
        # if to in intersphinx_pages:
        #     return pages[to].uri
        # # Original, but with replacement get_target_uri
        # return relative_uri(
        #     get_target_uri(from_),
        #     get_target_uri(to, typ),
        # )

    # Cribbed from StandaloneHTMLBuilder.handle_page
    default_baseuri = get_target_uri(pagename)
    def pathto(
        otheruri: str,
        resource: bool = False,
        baseuri: str = default_baseuri,
    ) -> str:
        """
        pathto() macro for jinja2 templates.

        Replaces StandaloneHTMLBuilder.handle_page.<local>.pathto() with
        intersphinx link hijacking.
        """
        # Intersphinx link hijacking
        if otheruri in self.intersphinx_pages:
            return self.intersphinx_pages[otheruri].uri
        # Original code, with replacement get_target_uri
        if resource and '://' in otheruri:
            # allow non-local resources given by scheme
            return otheruri
        elif not resource:
            otheruri = get_target_uri(otheruri)
        uri = relative_uri(baseuri, otheruri) or '#'
        if uri == '#' and not self.allow_sharp_as_current_path:
            uri = baseuri
        return uri
    context['pathto'] = pathto

    # Cribbed from StandaloneHTMLBuilder.get_doc_context
    # Uses replacement get_relative_uri() rather than self.get_relative_uri()
    # to define prev, next, parents and rellinks.
    prev = next = None
    parents = []
    rellinks = self.globalcontext['rellinks'][:]
    related = self.relations.get(docname)
    titles = self.env.titles
    if related and related[2]:
        try:
            next = {
                'link': get_relative_uri(docname, related[2]),
                'title': self.render_partial(titles[related[2]])['title'],
            }
            rellinks.append((related[2], next['title'], 'N', _('next')))
        except KeyError:
            next = None
    if related and related[1]:
        try:
            prev = {
                'link': get_relative_uri(docname, related[1]),
                'title': self.render_partial(titles[related[1]])['title'],
            }
            rellinks.append((related[1], prev['title'], 'P', _('previous')))
        except KeyError:
            # the relation is (somehow) not in the TOC tree, handle
            # that gracefully
            prev = None
    while related and related[0]:
        with contextlib.suppress(KeyError):
            parents.append({
                'link': get_relative_uri(docname, related[0]),
                'title': self.render_partial(titles[related[0]])['title'],
            })

        related = self.relations.get(related[0])
    if parents:
        # remove link to the master file; we have a generic
        # "back to index" link already
        parents.pop()
    parents.reverse()

    context['prev'] = prev
    context['next'] = next
    context['parents'] = parents
    context['rellinks'] = rellinks

    #print("   == rellinks=%(rellinks)s\n   &= parents=%(parents)s"%context)


def setup(app):
    f"""
    Registers the {__name__} extension with sphinx.

    Requires *sphinx.ext.intersphinx*.

    Transforms on *env-updated* and *html-page-context*.
    """
    app.add_config_value(
        name='remote_relations',
        default={},
        rebuild='html', # all html files change if you modify the remote
        types=[dict[str: list[str]]],
        description='Dictionary of parent pages and their siblings up to document root',
        )

    # We require intersphinx, so make sure it is included
    app.setup_extension('sphinx.ext.intersphinx')
    app.connect('env-updated', inject_remote_relations)
    _monkeypatch_link_hijack()
    # Don't need to tweak the context if we've applied the monkeypatch
    #app.connect('html-page-context', add_remote_nav_context)
