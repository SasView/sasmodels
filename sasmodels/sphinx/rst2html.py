r"""
Convert a restructured text document to html.

Inline math markup can uses the *math* directive, or it can use latex
style *\$expression\$*.  Math can be rendered using simple html and
unicode, or with mathjax.
"""

# CRUFT: locale.getlocale() fails on some versions of OS X
# See https://bugs.python.org/issue18378
import locale
import re
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

# TODO: remove _parse_localename cruft
# CRUFT: Old versions of locale did not support UTF-8 as a locale name.
if hasattr(locale, '_parse_localename'):
    try:
        locale._parse_localename('UTF-8')
    except ValueError:
        _old_parse_localename = locale._parse_localename
        def _parse_localename(localename):
            code = locale.normalize(localename)
            if code == 'UTF-8':
                return None, code
            else:
                return _old_parse_localename(localename)
        locale._parse_localename = _parse_localename

from docutils.core import publish_parts
from docutils.nodes import SkipNode, literal
from docutils.parsers.rst import Directive
from docutils.writers.html4css1 import HTMLTranslator, Writer

from .dollarmath import replace_dollar

# TODO: Get the files using importlib.resources
# TODO: The prolog is sasview/sasmodels specific... it doesn't belong in rst2html
# TODO: Remove all extra copies of rst_prolog
# from importlib import resources
THEME_PATH = Path(__file__).expanduser().resolve().parent
#STYLESHEET = THEME_PATH / "classic.css"
TEMPLATE = THEME_PATH / "template.txt"
RST_PROLOG_PATH = THEME_PATH / "prolog.rst"

#MATHJAX_PATH = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_CHTML"
MATHJAX_PATH = "https://cdn.jsdelivr.net/npm/mathjax@4/tex-mml-chtml.js" # recommended version as of 2026-11


if RST_PROLOG_PATH.exists():
    RST_PROLOG = RST_PROLOG_PATH.read_text()
else:
    # CRUFT: Fallback in case the resources aren't available in the install.
    RST_PROLOG = r"""
.. |Ang| unicode:: U+212B
.. |Ang^-1| replace:: |Ang|\ :sup:`-1`
.. |Ang^2| replace:: |Ang|\ :sup:`2`
.. |Ang^-2| replace:: |Ang|\ :sup:`-2`
.. |1e-6Ang^-2| replace:: 10\ :sup:`-6`\ |Ang|\ :sup:`-2`
.. |Ang^3| replace:: |Ang|\ :sup:`3`
.. |Ang^-3| replace:: |Ang|\ :sup:`-3`
.. |Ang^-4| replace:: |Ang|\ :sup:`-4`
.. |nm^-1| replace:: nm\ :sup:`-1`
.. |cm^-1| replace:: cm\ :sup:`-1`
.. |cm^2| replace:: cm\ :sup:`2`
.. |cm^-2| replace:: cm\ :sup:`-2`
.. |cm^3| replace:: cm\ :sup:`3`
.. |1e15cm^3| replace:: 10\ :sup:`15`\ cm\ :sup:`3`
.. |cm^-3| replace:: cm\ :sup:`-3`
.. |sr^-1| replace:: sr\ :sup:`-1`

.. |cdot| unicode:: U+00B7
.. |deg| unicode:: U+00B0
.. |g/cm^3| replace:: g\ |cdot|\ cm\ :sup:`-3`
.. |mg/m^2| replace:: mg\ |cdot|\ m\ :sup:`-2`
.. |fm^2| replace:: fm\ :sup:`2`
.. |Ang*cm^-1| replace:: |Ang|\ |cdot|\ cm\ :sup:`-1`
"""

# TODO: make a better sphinx stubs

def noop_directive(
    required_arguments=0,
    optional_arguments=0,
    final_argument_whitespace=False,
    has_content=False,
):
    class _NoOpDirective(Directive):
        """A minimal stub for Sphinx ``py:`` directives.

        The real directives expect a single required argument (the module,
        class, function name, …) and no content.  By declaring ``required_arguments
        = 1`` the parser will treat the argument correctly and will not try to
        interpret the following line as content, avoiding the "no content
        permitted" error.
        """
        def run(self):
            # Returning an empty list means the directive produces no output.
            return []
    _NoOpDirective.required_arguments = required_arguments          # e.g. ``sasmodels`` in ``.. py:currentmodule:: sasmodels``
    _NoOpDirective.optional_arguments = optional_arguments
    _NoOpDirective.final_argument_whitespace = final_argument_whitespace
    _NoOpDirective.has_content = has_content
    return _NoOpDirective

def noop_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    return [literal(text, text)], []

@contextmanager
def sphinx_stubs():
    """
    Temporarily register stubs for sphinx directives and roles. The originals
    will be restored when exiting the context.
    """
    from docutils.parsers.rst import directives, roles

    # Adding :orphan: as a role or a directive still left "orphan:" at the top of the rendered html.
    # Instead, lines containing :orphan: are automatically stripped.
    # See sasview/src/sas/qtgui/Perspectives/CorFunc/media/fdr-pdfs.rst

    _directives = directives._directives.copy()
    _role_registry = roles._role_registry.copy()
    _roles = roles._roles.copy()

    for name in "eq ref numref func class meth attr mod download".split():
        roles.register_canonical_role(name, noop_role)
    directives.register_directive('toctree', noop_directive(has_content=True))
    directives.register_directive('currentmodule', noop_directive(required_arguments=1))
    directives.register_directive('py:currentmodule', noop_directive(required_arguments=1))
    for name in "py:module py:class py:meth py:func".split():
        directives.register_directive(name, noop_directive())
    yield
    directives._directives = _directives
    roles._role_registry = _role_registry
    roles._roles = _roles

class WriterShim(Writer):
    """
    HTML writer that allows extra substitution arguments for the template.
    """
    def __init__(self, template_args=None):
        super().__init__()
        self.template_args = template_args
    def interpolation_dict(self):
        subs = super().interpolation_dict()
        if self.template_args:
            subs.update(self.template_args)
        #print("available subs", "\n".join(f"{k}: {v[:40]}" for k,v in subs.items()))
        return subs


class HTMLTranslatorShim(HTMLTranslator):
    # Suppress HTML error messages
    def visit_system_message(self, node):
        raise SkipNode

@dataclass
class URI:
    """
    Sphinx layout templates use uri.title and uri.link in various locations.

    Initialize with the base path of the reference, including filename without
    extension. Before using in a layout substitution set *uri.link = pathto(uri.base)*.
    This creates the full URI, including an implicit *.html* extension.

    We are not resolving the link in the constructor since we don't know the local
    and remote URIs when we create the link.
    """
    base: str
    title: str


def pseudo_sphinx(rst, path=None, title=None, context=(), rst_prolog=None, doc_root=None):
    """
    renders the document, linking to the sasview doc tree.
    """
    path = Path(path)
    if title is None:
        title = path.stem

    # TODO: fix prev, next, parents
    # TODO: use loops to render navigation parents
    # TODO: Pull DOC_ROOT/objects.inv to resolve links in the rst.
    # TODO: Remove the assumption that we are within plugins (needs parents, next and previous)
    # TODO: Allow the html to go into a cache directory.
    # TODO: Copy images to the correct directory relative to cache.
    # TODO: Write the rst file to cache/_sources/filename.rst
    # TODO: Use sphinx themes instead of hardcoding basic layout.html navigation and sidebar
    # TODO: Guess the sasview version from DOC_ROOT, which probably encodes the version number.
    # TODO: When modifying an existing file, allow linking to remote images
    # TODO: won't handle sphinx extensions
    # TODO: need full context; get it directly from sphinx?
    # TODO: If no cache, cache in auto-reaped tempdir so we don't have to delete.
    # TODO: Scan for figure/image tags to copy to cache

    # We should have a plugins directory cache with all the plugin help prebuilt,
    # and the associated img files copied. The rellink list should point to the next
    # and previous in the plugins directory.
    rst_root_path = None

    this_uri = title # base uri for document we are rendering
    local_uri = [this_uri] # list of uri base names that are resolved locally

    # ==== context used by the sphinx templates ===

    # Sphinx theme variables: see StandaloneHTMLBuilder.global_context
    #    https://github.com/sphinx-doc/sphinx/blob/master/sphinx/builders/html/__init__.py
    def pathto(uri, resource=False):
        """
        Resolve a uri base into a full URI.

        If uri base refers to a *local_uri*, expand into a reference relative to
        *path*, otherwise it is a reference relative to *doc_root*.

        if *resource* is True, then return the uri as is, possibly substituting
        *_sources/* with the parent of the reStructuredText sources tree. Sphinx
        is opinionated about the location of the rst files when *html_copy_source*
        is True in the sphinx configuration.

        The variables *rst_root_path*, *local_uri*, *path* and *doc_root* are
        local to pseudo_sphinx.
        """
        nonlocal rst_root_path, local_uri, path, doc_root
        if resource:
            # If we need to redirect rst files from _sources, do it here
            if rst_root_path and uri.startswith('_sources/'):
                uri = uri.replace('_sources', str(rst_root_path))
            return uri
        # if uri is a sister plugin model in the cache, then return that path; if it is
        # coming from the inventory, return relative to docroot.
        res = f"{str(path.parent)}/{uri}.html" if uri in local_uri else f"{doc_root}/{uri}.html"
        return res

    # Expand uri.base into uri.link using pathto()
    for _uri in context: # resolve links
        _uri.link = pathto(_uri.base)

    # Navigation state
    prev, next, _root_doc, *parents = context

    def hasdoc(uri):
        # Used by the template engine to make decisions about what to include.
        # We could implement it using requests, trying to load the files from DOC_ROOT.
        # If we start using the sphinx themes, though, I think it will be okay to drop some parts
        # of the page when we are rendering the local document.
        return True

    def accesskey(key):
        return f'accesskey="{key}"' if key else ''

    # Note: "next" is the name in the context, so overriding the builtin next locally (bad form!)
    prev, next = context[0], context[1]
    #   [page name, link title, accesskey, link text],
    rellinks: list[tuple[str, str, str, str]] = []
    if hasdoc('genindex'):
        rellinks.append(('genindex', 'General Index', 'I', 'index'))
    if hasdoc('py-modindex'):
        rellinks.append(('py-modindex', 'Python Module Index', '', 'modules'))
    if next is not None:
        rellinks.append((next.base, next.title, 'N', 'next'))
    if prev is not None:
        rellinks.append((prev.base, prev.title, 'P', 'previous'))
    root_doc, shorttitle = _root_doc.base, _root_doc.title
    link = pathto(this_uri)
    reldelim1 = " »" # the default if not defined
    reldelim2 = " |" # the default if not defined
    # Ick! The default layout assumes sources are in `_sources/{sourcename}`, but we can
    # override this in pathto().
    sourcename = path.name
    copyright = "2026, The SasView Project" # from conf.py
    last_updated = None # If defined, expands to "Last updated on {last_updated}." in the footer.
    show_sphinx = False # Not created using sphinx and we don't know the sphinx version.
    sphinx_version = "9.1.0" # from sphinx
    sidebars = ["localtoc.html", "relations.html", "sourcelink.html", "searchbox.html"] # from theme.toml
    stylesheets = ["classic.css"] # from theme.toml; I didn't check if the template sees this symbol.

    css_files = [
        (THEME_PATH / css).relative_to(path.parent, walk_up=True)
        for css in stylesheets
    ]

    # likely quite a bit more we could include, but our current sphinx theme isn't using them.

    # ==== end of context ====

    # Copy of the classic layout navigation, sidebar and footers, with macros expanded.
    # https://github.com/sphinx-doc/sphinx/blob/master/sphinx/themes/basic/layout.html
    # This is good enough to show the rst/latex rendering but it may not match the
    # rest of the docs.
    _rendered_rellinks = "\n".join(
        f"""<li class="right"{' style="margin-right: 10px"' if k == 0 else ''}>
          <a href="{pathto(rellink[0])}" title="{pathto(rellink[1])}" {accesskey(rellink[2])}>{rellink[3]}</a>{reldelim2 if k > 0 else ''}</li>"""
        for k, rellink in enumerate(rellinks) if rellink
    )
    _rendered_parents = "\n".join(
        f"""<li class="nav-item nav-item-{1}"><a href="{parent.link}">{parent.title}</a>{reldelim1}</li>"""
        for parent in parents
    )
    navigation = f"""\
<div class="related" role="navigation" aria-label="Related">
      <h3>Navigation</h3>
      <ul>
        {_rendered_rellinks}
        <li class="nav-item nav-item-0"><a href="{pathto(root_doc)}">{shorttitle}</a>{reldelim1}</li>
          {_rendered_parents}
        <li class="nav-item nav-item-this"><a href="{link}">{title}</a></li>
      </ul>
    </div>
"""
    sidebar = f"""\
<div class="sphinxsidebar" role="navigation" aria-label="Main">
<div class="sphinxsidebarwrapper">
  <div>
    <h4>Previous topic</h4>
    <p class="topless"><a href="{prev.link}" title="previous chapter">{prev.title}</a></p>
  </div>
  <div>
    <h4>Next topic</h4>
    <p class="topless"><a href="{next.link}" title="next chapter">{next.title}</a></p>
  </div>
  <div role="note" aria-label="source link">
    <h3>This Page</h3>
    <ul class="this-page-menu">
      <li><a href="{pathto('_sources/' + sourcename, True)}" rel="nofollow">Show Source</a></li>
    </ul>
   </div>
<search id="searchbox" style="display: block;" role="search">
  <h3 id="searchlabel">Quick search</h3>
    <div class="searchformwrapper">
    <form class="search" action="{pathto('search')}" method="get">
      <input type="text" name="q" aria-labelledby="searchlabel" autocomplete="off" autocorrect="off" autocapitalize="off" spellcheck="false">
      <input type="submit" value="Go">
    </form>
    </div>
</search>
<script>document.getElementById('searchbox').style.display = "block"</script>
        </div>
      </div>
"""
    # Original footer. It's problematic because it needs copyright year and sphinx version.
    footer = f"""\
<div class="footer" role="contentinfo">
    © Copyright {copyright}.
    Created using <a href="https://www.sphinx-doc.org/">Sphinx</a> {sphinx_version}.
    </div>"
"""
    footer = "" # suppressing because the page might not be from SasView and the generator was not Sphinx.

    # Extra blocks available to the docutils template engine, used by TEMPLATE
    template_args = dict(navigation=navigation, sidebar=sidebar, sphinx_footer=footer)

    return rst2html(
        rst=rst, rst_prolog=RST_PROLOG_PATH, template=TEMPLATE, template_args=template_args, css_list=css_files,
    )

def rst2html(
        rst, part="whole", math_output="mathjax",
        rst_prolog=None, template=None, css_list=None,
        template_args=None,
        ):
    r"""
    Convert restructured text into simple html.

    Valid *math_output* formats for formulas include:
    - HTML
    - MathML
    - MathJax

    See `<http://docutils.sourceforge.net/docs/user/config.html#math-output>`_
    for details. The MathJax library url is defined in rst2html.MATHJAX_PATH.

    The following *part* choices are available:
    - whole: the entire html document
    - html_body: document division with title and contents and footer
    - body: contents only

    There are other parts, but they don't make sense alone:

        subtitle, version, encoding, html_prolog, header, meta,
        html_title, title, stylesheet, html_subtitle, html_body,
        body, head, body_suffix, fragment, docinfo, html_head,
        head_prefix, body_prefix, footer, body_pre_docinfo, whole

    *rst_prolog* is the path to the reStructureText prolog file
    """
    if rst_prolog and Path(rst_prolog).exists():
        prolog = Path(rst_prolog).read_text()
        rst = f"{prolog}\n{rst}"

    # Ick! mathjax doesn't work properly with math-output, and the
    # others don't work properly with math_output!
    if math_output == "mathjax":
        # TODO: this is copied from docs/conf.py; there should be only one
        settings = {"math_output": math_output + " " + MATHJAX_PATH}
    else:
        settings = {"math-output": math_output}

    if css_list:
        settings["embed_stylesheet"] = False
        # The stylesheet_path setting makes paths relative to current directory.
        # Clear it out so that we can instead use the stylesheet paths given by the caller.
        settings["stylesheet_path"] = None
        settings["stylesheet"] = css_list

    # math2html and mathml do not support \frac12
    # mathml, html do not support \tfrac
    if math_output in ("mathml", "html"):
        rst = replace_compact_fraction(rst)
        rst = rst.replace(r'\tfrac', r'\frac')

    if template:
        settings["template"] = template

    # TODO: docutils doesn't support :orphan: metadata
    # TODO: docutils math doesn't support :nowrap: or :label:
    # The :nowrap: option is used when the equation has its own \begin{align*}...\end{align*}
    # The docutils helper pick_math_environment() looks for \\ in the text, and if it
    # sees it, it wraps the block in "align" rather than "equation".
    # We will strip label, nowrap, \begin{align*} and \end{align*}.
    # This is too simple: if the author has :nowrap: in sphinx for a multiline equation, but
    # doesn't include their own \begin...\end block then it will display fine in the preview
    # but fail when rendering with sphinx. Doing this correctly is too much work.
    pattern = r"^\s*:(?:nowrap|no[-_]wrap|label|orphan):.*\n?"
    rst = re.sub(pattern, "", rst, flags=re.MULTILINE)
    rst = re.sub(r"\\(?:begin|end){align\*?}", "", rst, flags=re.MULTILINE)
    #print(f"=== rst ===\n", rst)

    rst = replace_dollar(rst)
    writer = WriterShim(template_args=template_args)
    #writer.translator_class = HTMLTranslatorShim
    #with suppress_html_errors(),
    with sphinx_stubs():
        parts = publish_parts(
            source=rst,
            writer=writer,
            #writer_name='html',
            settings_overrides=settings,
            )
    return parts[part]

@contextmanager
def suppress_html_errors():
    r"""
    Context manager for keeping error reports out of the generated HTML.

    Within the context, system message nodes in the docutils parse tree
    will be ignored.  After the context, the usual behaviour will be restored.
    """
    visit_system_message = HTMLTranslator.visit_system_message
    HTMLTranslator.visit_system_message = _skip_node
    yield None
    HTMLTranslator.visit_system_message = visit_system_message

def _skip_node(self, node):
    raise SkipNode


_compact_fraction = re.compile(r"(\\[cdt]?frac)([0-9])([0-9])")
def replace_compact_fraction(content):
    r"""
    Convert \frac12 to \frac{1}{2} for broken latex parsers
    """
    return _compact_fraction.sub(r"\1{\2}{\3}", content)


def load_rst_as_html(filename):
    """Load rst from file and convert to html"""
    # TODO: Take a configuration from elsewhere rather than assuming sasview docs
    from ..generate import DOC_ROOT  # Ick! Circular import of sasmodels specific stuff

    sasview_version = "" # Can pull this from DOC_TREE
    sasview_doc = URI("index", f"SasView {sasview_version} Documentation")
    user_doc = URI("user/user", "SasView User Documentation")
    # prev next root parents
    context = (user_doc, user_doc, sasview_doc, user_doc)

    path = Path(filename).resolve()
    with open(path) as fid:
        rst = fid.read()
    return pseudo_sphinx(rst, path=path, context=context, doc_root=DOC_ROOT)
    #return rst2html(rst=f"{RST_PROLOG}\n{rst}", css_list=[stylesheet])

def wxview(html, url="", size=(850, 540)):
    # type: (str, str, tuple[int, int]) -> "wx.Frame"
    """View HTML in a wx dialog"""
    import wx
    from wx.html2 import WebView

    frame = wx.Frame(None, -1, size=size)
    view = WebView.New(frame)
    view.SetPage(html, url)
    frame.Show()
    return frame

def view_html_wxapp(html, url=""):
    # type: (str, str) -> None
    """HTML viewer app in wx"""
    import wx  # type: ignore

    app = wx.App()
    frame = wxview(html, url)  # pylint: disable=unused-variable
    app.MainLoop()

def view_url_wxapp(url):
    # type: (str) -> None
    """URL viewer app in wx"""
    import wx  # type: ignore
    from wx.html2 import WebView

    app = wx.App()
    frame = wx.Frame(None, -1, size=(850, 540))
    view = WebView.New(frame)
    view.LoadURL(url)
    frame.Show()
    app.MainLoop()

def can_use_wx() -> bool:
    """Return True if wx web viewer is available."""
    try:
        import wx
        return True
    except ImportError:
        return False

def qtview(html, url=""):
    # type: (str, str) -> "QWebView"
    """View HTML in a Qt dialog"""
    from PySide6.QtCore import QUrl
    from PySide6.QtWebEngineWidgets import QWebEngineView

    helpView = QWebEngineView()
    helpView.setHtml(html, QUrl(url))
    helpView.show()
    return helpView

def view_html_qtapp(html, url=""):
    # type: (str, str) -> None
    """HTML viewer app in Qt"""
    import sys

    from PySide6.QtWidgets import QApplication

    app = QApplication([])
    frame = qtview(html, url)  # pylint: disable=unused-variable
    sys.exit(app.exec_())

def view_url_qtapp(url):
    # type: (str) -> None
    """URL viewer app in Qt"""
    import sys

    from PySide6.QtCore import QUrl
    from PySide6.QtWebEngineWidgets import QWebEngineView
    from PySide6.QtWidgets import QApplication

    app = QApplication([])
    frame = QWebEngineView()
    frame.load(QUrl(url))
    frame.show()
    sys.exit(app.exec_())

def can_use_qt() -> bool:
    """Return True if Qt web viewer is available."""
    try:
        from PySide6.QtWebEngineWidgets import QWebEngineView
        return True
    except ImportError:
        return False

def view_html_browser(html, url, overwrite=False):
    # Show the docs in the default browser
    import time
    import webbrowser

    # Write the html to the target file, found by removing file:// from the url
    html_file = Path(url[7:])
    if html_file.exists() and not overwrite:
        delete_after = False
        print(f"Warning: showing existing {str(html_file)}")
    else:
        delete_after = True
        html_file.write_text(html)
    webbrowser.open(url)

    # Give the file time to open in the browser, then delete and exit the program.
    # This may fail on Windows if it holds the file open while viewing
    #delete_after = False
    if delete_after:
        time.sleep(1.0)
        html_file.unlink()

def view_url_browser(url, autoraise=False):
    import webbrowser

    webbrowser.open(url, autoraise=autoraise)


# Set default html viewer
view_html = view_html_browser

def view_help(filename, viewer="browser"):
    # type: (str, bool) -> None
    """
    View rst or html file.
    If *qt* use q viewer, otherwise use wx.
    """
    if viewer == "qt" and not can_use_qt():
        viewer = "browser"
        print("QtWebKit is not available. Use browser to view help")
    if viewer == "wx" and not can_use_wx():
        viewer = "browser"
        print("wb is not available. Use browser to view help")

    url = Path(filename).expanduser().resolve().as_uri()  # file://{absolute path}
    if filename.endswith('.rst'):
        # TODO: this fails without stubs for sphinx specific roles and directives
        html = load_rst_as_html(filename)
        if viewer == "browser":
            view_html_browser(html, url + ".html", overwrite=True)
        elif viewer == "qt":
            view_html_qtapp(html, url + ".html")
        else:
            view_html_wxapp(html, url + ".html")
    else:
        if viewer == "browser":
            view_url_browser(url)
        elif viewer == "qt":
            view_url_qtapp(url)
        else:
            view_url_wxapp(url)

def main():
    # type: () -> None
    """Command line interface to rst or html viewer."""
    import sys

    view_help(sys.argv[1])

if __name__ == "__main__":
    main()
