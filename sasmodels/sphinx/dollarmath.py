# This program is public domain
# Author: Paul Kienzle
r"""
Allow \$math\$ and \\\$\$math\\\$\$ markup in text and docstrings, ignoring \\\$.

For inline math, \$...\$ is translated into a math role \:math\:\`...\`.
The math markup should be separated from the surrounding text by spaces,
To embed markup within a word, place backslash-space before and after.
For convenience, \$ can be used beside some punctuation. To put \$ in the
text it needs to be escaped with a backslash.

Display math "... \\\$\$ ... \\\$\$ ..." is translated into a math block:

.. code:: restructuredtext

    ...

    .. math::

        ...

    ...


**Examples**

* \$a\$-\$b\$ becomes $a$-$b$
* (\$a\$, \$b\$) becomes ($a$, $b$); to remove the space (\$a\$,\\ \$b\$) becomes ($a$,\ $b$)
* \$f(x)\$. becomes $f(x)$. Similarly for , : ; and ?
* K\\ \$\alpha\$ becomes K\ $\alpha$
* \\\\\$3.56 becomes \$3.56; normally you won't need the escape since $3.56 doesn't have a closing dollar sign.
* \\\$\$ x &= 1 \\\\ y &= 2 \\\$\$ becomes $$ x &=1 \\ y &= 2 $$
* \\\\\$\\\\\$...\\\\\$\\\\\$ becomes \\\$\$...\\\$\$ which gets expanded by MathJax. To
  show \\\$\$ in your text you need to use \\\\\\\\\\\\\\$\\\\\$.

If you are using *sphinx.ext.autodoc* then it needs to appear in the extension
list before *dollarmath*.

The *dollarmath* extension is incompatible with *myst_parser*.

The *dollarmath* extension transforms the reStructureText source before it is parsed, so the
substitution happens even when it occurs in a comment block.
"""

# TODO: Support myst-parser.
# TODO: requires sphinx.ext.autodoc before dollarmath on the extensions list

import re
import textwrap

# Match $...$
_inline_math = re.compile(
r"""
    (?:           # Non-captured look-behind before $
      ^           # Allow $ at start of line
      | (?<=\s|[-(]) # Allow leading space, dash, comma or open parenthesis
    )
    [$]           # Opening $  (one of ^$, \s$, -$ or \($, but not \\$)
    ([^\n]*?)     # Capture everything on the line up to the next $, non-greedy
    (?<![\\])     # Allow embedded \$ using negative look-behind
    [$]           # closing $
    (?:           # Non-captured look-ahead after $
      $           # Allow $ at the end of line
      | (?=\s|[-.,;:?\\)]) # Allow trailing space or punctuation -.,;:?\)
    )
""", re.VERBOSE)

# Match $$...$$
_display_math = re.compile(
r"""
    (?:
      ^           # Allow $$ at the start of line
      | \s+       # Eat leading space
    )
    [$][$]        # Opening $$
    (.*?)         # Capture everything up to the next $$, non-greedy
    (?<![\\])     # Allow embedded \$ using negative look-behind
    [$][$]        # closing $$
    (?:
      $           # Allow $$ at the end of line
      | \s+       # Eat trailing space
    )
""", re.VERBOSE|re.DOTALL)

# Match \$
_escaped_dollar = re.compile(r"\\[$]") # Match \$ so it can be replaced by $ after prior transform

def replace_dollar(content):
    r"""
    Converts latex math markup using dollar signs into sphinx math markup.

    For inline math, \$...\$ becomes \:math\:\`...\`. For display math, \\\$\$ ... \\\$\$
    becomes a math block:

    .. code:: restructuredtext

        .. math::

            ...

    """
    # original = content
    # print("text:", repr(content))
    content = _display_math.sub(_display_math_sub, content)
    content = _inline_math.sub(r":math:`\1`", content)
    content = _escaped_dollar.sub("$", content)
    # print("==>", repr(content))
    # # For debugging within sphinx, write directly to stdout
    # if '$' in content:
    #     import sys
    #     sys.stdout.write("\n========> not converted\n")
    #     sys.stdout.write(content)
    #     sys.stdout.write("\n")
    # elif '$' in original:
    #     import sys
    #     sys.stdout.write("\n========> converted\n")
    #     sys.stdout.write(content)
    #     sys.stdout.write("\n")
    return content

def _display_math_sub(match):
    return rf"""

.. math::

{textwrap.indent(match.group(1), "    ")}

"""

def _rewrite_rst(app, docname, source):
    source[0] = replace_dollar(source[0])

def _rewrite_autodoc(app, what, name, obj, options, lines):
    lines[:] = replace_dollar("\n".join(lines)).split("\n")

def setup(app):
    """
    Register sphinx event listeners for source-read (for rst source) and
    autodoc-process-docstring (for python docstrings). These transform the
    input text, replacing the dollar signs with the equivalent restructuredtext
    math commands.
    """
    from sphinx.errors import SphinxError

    if "myst_parser" in app.config.extensions:
        raise SphinxError(f"The {__name__} extension is incompatible with myst_parser")

    app.connect('source-read', _rewrite_rst)
    if 'autodoc-process-docstring' in app.events.events:
        app.connect('autodoc-process-docstring', _rewrite_autodoc)


def test_dollar():

    # display math tests
    assert replace_dollar("$$only$$") == "\n\n.. math::\n\n    only\n\n"
    assert replace_dollar("$$first$$ is good") == "\n\n.. math::\n\n    first\n\nis good"
    assert replace_dollar("so is $$last$$") == "so is\n\n.. math::\n\n    last\n\n"
    assert replace_dollar("and $$mid$$ too") == "and\n\n.. math::\n\n    mid\n\ntoo"
    assert replace_dollar("   $$indented$$") == "\n\n.. math::\n\n    indented\n\n"
    assert replace_dollar("   $$indented$$ ") == "\n\n.. math::\n\n    indented\n\n"
    assert replace_dollar("new line\n\n   $$indented$$") == "new line\n\n.. math::\n\n    indented\n\n"
    assert replace_dollar("$$only$$   \n\nmore") == "\n\n.. math::\n\n    only\n\nmore"
    assert replace_dollar("$$multiline\nequation$$""") == "\n\n.. math::\n\n    multiline\n    equation\n\n"
    assert replace_dollar("\n$$\n   math\n   more math\n$$\n") == "\n\n.. math::\n\n\n       math\n       more math\n\n\n\n"

    # inline math tests
    assert replace_dollar("no dollar") == "no dollar"
    assert replace_dollar("$only$") == ":math:`only`"
    assert replace_dollar("$first$ is good") == ":math:`first` is good"
    assert replace_dollar("so is $last$") == "so is :math:`last`"
    assert replace_dollar("and $mid$ too") == "and :math:`mid` too"
    assert replace_dollar("$first$, $mid$, $last$") == ":math:`first`, :math:`mid`, :math:`last`"
    assert replace_dollar(r"dollar\$ escape") == "dollar$ escape"
    assert replace_dollar(r"dollar \$escape\$ too") == "dollar $escape$ too"
    assert replace_dollar("spaces $in the$ math") == "spaces :math:`in the` math"
    assert replace_dollar(r"emb\ $ed$\ ed") == r"emb\ :math:`ed`\ ed"
    assert replace_dollar("$first$a") == "$first$a"
    assert replace_dollar("a$last$") == "a$last$"
    assert replace_dollar("$37") == "$37"
    assert replace_dollar("($37)") == "($37)"
    assert replace_dollar("$37 - $43") == "$37 - $43"
    assert replace_dollar("($37, $38)") == "($37, $38)"
    assert replace_dollar("a $mid$dle a") == "a $mid$dle a"
    assert replace_dollar("a ($in parens$) a") == "a (:math:`in parens`) a"
    assert replace_dollar("a (again $in parens$) a") == "a (again :math:`in parens`) a"


if __name__ == "__main__":
    test_dollar()
