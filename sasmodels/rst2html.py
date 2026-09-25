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
from docutils.parsers.rst import Directive, directives, roles
from docutils.writers.html4css1 import HTMLTranslator

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

def sphinx_stubs():
    for name in "ref numref mod func class meth".split():
        roles.register_canonical_role(name, noop_role)
    directives.register_directive('toctree', noop_directive(has_content=True))
    directives.register_directive('currentmodule', noop_directive(required_arguments=1))
    directives.register_directive('py:currentmodule', noop_directive(required_arguments=1))
    for name in "py:module py:class py:func".split():
        directives.register_directive(name, noop_directive())
sphinx_stubs()

#MATHJAX_PATH = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_CHTML"
MATHJAX_PATH = "https://cdn.jsdelivr.net/npm/mathjax@4/tex-mml-chtml.js" # recommended version as of 2026-11

def rst2html(rst, part="whole", math_output="mathjax", rst_prolog=None, css_list=None):
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
    if rst_prolog:
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
        # The sytlesheet_path setting makes paths relative to current directory.
        # Clear it out so that we can instead use the stylesheet paths given by the caller.
        settings["stylesheet_path"] = None
        settings["stylesheet"] = css_list

    # math2html and mathml do not support \frac12
    # mathml, html do not support \tfrac
    if math_output in ("mathml", "html"):
        rst = replace_compact_fraction(rst)
        rst = rst.replace(r'\tfrac', r'\frac')

    # TODO: docutils math doesn't support :nowrap: or :label:
    pattern = r'^\s\s*:(?:nowrap|no[-_]wrap|label):.*\n?'
    rst = re.sub(pattern, '', rst, flags=re.MULTILINE)

    rst = replace_dollar(rst)
    with suppress_html_errors():
        parts = publish_parts(
            source=rst, writer_name='html',
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


_dollar = re.compile(r"(?:^|(?<=\s|[-(]))[$]([^\n]*?)(?<![\\])[$](?:$|(?=\s|[-.,;:?\\)]))")
_notdollar = re.compile(r"\\[$]")
def replace_dollar(content):
    r"""
    Convert dollar signs to inline math markup in rst.
    """
    content = _dollar.sub(r":math:`\1`", content)
    content = _notdollar.sub("$", content)
    return content


def test_dollar():
    """
    Test substitution of dollar signs with equivalent RST math markup
    """
    assert replace_dollar("no dollar") == "no dollar"
    assert replace_dollar("$only$") == ":math:`only`"
    assert replace_dollar("$first$ is good") == ":math:`first` is good"
    assert replace_dollar("so is $last$") == "so is :math:`last`"
    assert replace_dollar("and $mid$ too") == "and :math:`mid` too"
    assert replace_dollar("$first$, $mid$, $last$") == ":math:`first`, :math:`mid`, :math:`last`"
    assert replace_dollar("dollar\\$ escape") == "dollar$ escape"
    assert replace_dollar("dollar \\$escape\\$ too") == "dollar $escape$ too"
    assert replace_dollar("spaces $in the$ math") == "spaces :math:`in the` math"
    assert replace_dollar("emb\\ $ed$\\ ed") == "emb\\ :math:`ed`\\ ed"
    assert replace_dollar("$first$a") == "$first$a"
    assert replace_dollar("a$last$") == "a$last$"
    assert replace_dollar("$37") == "$37"
    assert replace_dollar("($37)") == "($37)"
    assert replace_dollar("$37 - $43") == "$37 - $43"
    assert replace_dollar("($37, $38)") == "($37, $38)"
    assert replace_dollar("a $mid$dle a") == "a $mid$dle a"
    assert replace_dollar("a ($in parens$) a") == "a (:math:`in parens`) a"
    assert replace_dollar("a (again $in parens$) a") == "a (again :math:`in parens`) a"

def load_rst_as_html(filename):
    """Load rst from file and convert to html"""
    from .generate import RST_PROLOG, STYLESHEET  # Ick! Circular import of sasmodels specific stuff

    # Make stylesheet path relative to the html file
    filename = Path(filename).expanduser().absolute()
    stylesheet = STYLESHEET.relative_to(filename.parent, walk_up=True)

    with open(filename) as fid:
        rst = fid.read()
    return rst2html(rst=f"{RST_PROLOG}\n{rst}", css_list=[stylesheet])

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
        time.sleep(0.5)
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

    url = Path(filename).expanduser().absolute().as_uri()  # file://{absolute path}
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
