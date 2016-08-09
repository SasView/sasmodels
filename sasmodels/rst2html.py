r"""
Convert a restructured text document to html.

Inline math markup can uses the *math* directive, or it can use latex
style *\$expression\$*.  Math can be rendered using simple html and
unicode, or with mathjax.
"""

import re
from contextlib import contextmanager

from docutils.core import publish_parts
from docutils.writers.html4css1 import HTMLTranslator
from docutils.nodes import SkipNode


def rst2html(rst, part="whole", math_output="html"):
    r"""
    Convert restructured text into simple html.

    Valid *math_output* formats for formulas include:
    - html
    - mathml
    - mathjax
    See `http://docutils.sourceforge.net/docs/user/config.html#math-output`_
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
        settings = {"math_output": math_output}
    else:
        settings = {"math-output": math_output}

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


_dollar = re.compile(r"(?:^|(?<=\s|[(]))[$]([^\n]*?)(?<![\\])[$](?:$|(?=\s|[.,;)\\]))")
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
    assert replace_dollar(ur"dollar\$ escape") == u"dollar$ escape"
    assert replace_dollar(ur"dollar \$escape\$ too") == u"dollar $escape$ too"
    assert replace_dollar(u"spaces $in the$ math") == u"spaces :math:`in the` math"
    assert replace_dollar(ur"emb\ $ed$\ ed") == ur"emb\ :math:`ed`\ ed"
    assert replace_dollar(u"$first$a") == u"$first$a"
    assert replace_dollar(u"a$last$") == u"a$last$"
    assert replace_dollar(u"$37") == u"$37"
    assert replace_dollar(u"($37)") == u"($37)"
    assert replace_dollar(u"$37 - $43") == u"$37 - $43"
    assert replace_dollar(u"($37, $38)") == u"($37, $38)"
    assert replace_dollar(u"a $mid$dle a") == u"a $mid$dle a"
    assert replace_dollar(u"a ($in parens$) a") == u"a (:math:`in parens`) a"
    assert replace_dollar(u"a (again $in parens$) a") == u"a (again :math:`in parens`) a"

if __name__ == "__main__":
    test_dollar()
