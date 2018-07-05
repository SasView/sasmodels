r"""
Convert a restructured text document to html.

Inline math markup can uses the *math* directive, or it can use latex
style *\$expression\$*.  Math can be rendered using simple html and
unicode, or with mathjax.
"""

import re
from contextlib import contextmanager

# CRUFT: locale.getlocale() fails on some versions of OS X
# See https://bugs.python.org/issue18378
import locale
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
from docutils.writers.html4css1 import HTMLTranslator
from docutils.nodes import SkipNode

def rst2html(rst, part="whole", math_output="mathjax"):
    r"""
    Convert restructured text into simple html.

    Valid *math_output* formats for formulas include:
    - html
    - mathml
    - mathjax
    See `<http://docutils.sourceforge.net/docs/user/config.html#math-output>`_
    for details.

    The following *part* choices are available:
    - whole: the entire html document
    - html_body: document division with title and contents and footer
    - body: contents only

    There are other parts, but they don't make sense alone:

        subtitle, version, encoding, html_prolog, header, meta,
        html_title, title, stylesheet, html_subtitle, html_body,
        body, head, body_suffix, fragment, docinfo, html_head,
        head_prefix, body_prefix, footer, body_pre_docinfo, whole
    """
    # Ick! mathjax doesn't work properly with math-output, and the
    # others don't work properly with math_output!
    if math_output == "mathjax":
        # TODO: this is copied from docs/conf.py; there should be only one
        mathjax_path = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_CHTML"
        settings = {"math_output": math_output + " " + mathjax_path}
    else:
        settings = {"math-output": math_output}

    # TODO: support stylesheets
    #html_root = "/full/path/to/_static/"
    #sheets = [html_root+s for s in ["basic.css","classic.css"]]
    #settings["embed_styesheet"] = True
    #settings["stylesheet_path"] = sheets

    # math2html and mathml do not support \frac12
    rst = replace_compact_fraction(rst)

    # mathml, html do not support \tfrac
    if math_output in ("mathml", "html"):
        rst = rst.replace(r'\tfrac', r'\frac')

    rst = replace_dollar(rst)
    with suppress_html_errors():
        parts = publish_parts(source=rst, writer_name='html',
                              settings_overrides=settings)
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
    assert replace_dollar(u"no dollar") == u"no dollar"
    assert replace_dollar(u"$only$") == u":math:`only`"
    assert replace_dollar(u"$first$ is good") == u":math:`first` is good"
    assert replace_dollar(u"so is $last$") == u"so is :math:`last`"
    assert replace_dollar(u"and $mid$ too") == u"and :math:`mid` too"
    assert replace_dollar(u"$first$, $mid$, $last$") == u":math:`first`, :math:`mid`, :math:`last`"
    assert replace_dollar(u"dollar\\$ escape") == u"dollar$ escape"
    assert replace_dollar(u"dollar \\$escape\\$ too") == u"dollar $escape$ too"
    assert replace_dollar(u"spaces $in the$ math") == u"spaces :math:`in the` math"
    assert replace_dollar(u"emb\\ $ed$\\ ed") == u"emb\\ :math:`ed`\\ ed"
    assert replace_dollar(u"$first$a") == u"$first$a"
    assert replace_dollar(u"a$last$") == u"a$last$"
    assert replace_dollar(u"$37") == u"$37"
    assert replace_dollar(u"($37)") == u"($37)"
    assert replace_dollar(u"$37 - $43") == u"$37 - $43"
    assert replace_dollar(u"($37, $38)") == u"($37, $38)"
    assert replace_dollar(u"a $mid$dle a") == u"a $mid$dle a"
    assert replace_dollar(u"a ($in parens$) a") == u"a (:math:`in parens`) a"
    assert replace_dollar(u"a (again $in parens$) a") == u"a (again :math:`in parens`) a"

def load_rst_as_html(filename):
    from os.path import expanduser
    with open(expanduser(filename)) as fid:
        rst = fid.read()
    html = rst2html(rst)
    return html

def wxview(html, url="", size=(850, 540)):
    import wx
    from wx.html2 import WebView
    frame = wx.Frame(None, -1, size=size)
    view = WebView.New(frame)
    view.SetPage(html, url)
    frame.Show()
    return frame

def view_html_wxapp(html, url=""):
    import wx  # type: ignore
    app = wx.App()
    frame = wxview(html, url)
    app.MainLoop()

def view_url_wxapp(url):
    import wx  # type: ignore
    from wx.html2 import WebView
    app = wx.App()
    frame = wx.Frame(None, -1, size=(850, 540))
    view = WebView.New(frame)
    view.LoadURL(url)
    frame.Show()
    app.MainLoop()

def qtview(html, url=""):
    try:
        from PyQt5.QtWebKitWidgets import QWebView
        from PyQt5.QtCore import QUrl
    except ImportError:
        from PyQt4.QtWebKit import QWebView
        from PyQt4.QtCore import QUrl
    helpView = QWebView()
    helpView.setHtml(html, QUrl(url))
    helpView.show()
    return helpView

def view_html_qtapp(html, url=""):
    import sys
    try:
        from PyQt5.QtWidgets import QApplication
    except ImportError:
        from PyQt4.QtGui import QApplication
    app = QApplication([])
    frame = qtview(html, url)
    sys.exit(app.exec_())

def view_url_qtapp(url):
    import sys
    try:
        from PyQt5.QtWidgets import QApplication
    except ImportError:
        from PyQt4.QtGui import QApplication
    app = QApplication([])
    try:
        from PyQt5.QtWebKitWidgets import QWebView
        from PyQt5.QtCore import QUrl
    except ImportError:
        from PyQt4.QtWebKit import QWebView
        from PyQt4.QtCore import QUrl
    frame = QWebView()
    frame.load(QUrl(url))
    frame.show()
    sys.exit(app.exec_())

# Set default html viewer
view_html = view_html_qtapp

def can_use_qt():
    """
    Return True if QWebView exists.

    Checks first in PyQt5 then in PyQt4
    """
    try:
        from PyQt5.QtWebKitWidgets import QWebView
        return True
    except ImportError:
        try:
            from PyQt4.QtWebKit import QWebView
            return True
        except ImportError:
            return False

def view_help(filename, qt=False):
    import os

    if qt:
        qt = can_use_qt()

    url = "file:///"+os.path.abspath(filename).replace("\\", "/")
    if filename.endswith('.rst'):
        html = load_rst_as_html(filename)
        if qt:
            view_html_qtapp(html, url)
        else:
            view_html_wxapp(html, url)
    else:
        if qt:
            view_url_qtapp(url)
        else:
            view_url_wxapp(url)

def main():
    import sys
    view_help(sys.argv[1], qt=False)

if __name__ == "__main__":
    main()
