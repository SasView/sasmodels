"""
Read a relaxed JSON file.

Relaxed JSON allows comments (introduced by // and going to the end
of the line), optional quotes on key names in dictionaries and optional
trailing commas at the ends of lists and dictionaries.  It also strips
the leading characters up to the first '{'.  Multiline strings can be
formatted using "\n\" at the end of each line.

If the file contains e.g., "var _ = {...}", then it can be edited with
a JavaScript aware editor such as VJET for eclipse, and it will be easier
to locate errors and adjust formatting.
"""
import re
import json
from contextlib import contextmanager

try:
    from collections import OrderedDict
except:
    from ordered_dict import OrderedDict


_LEADING_TEXT = re.compile(r'^.*?[{]', re.DOTALL)
_LINE_CONTINUATION = re.compile(r'\\\s*\n')
_TRAILING_COMMENT = re.compile(r'(?P<comment>\s*//.*?)\n')
_MULTILINE_COMMENT = re.compile(r'(?P<comment>/\*.*?\*/)', re.DOTALL)
_UNQUOTED_FIELDNAME = re.compile(r'(?P<prefix>[,{]\s*)(?P<key>[^\s,{}:"]+)(?P<tail>\s*:)')
_TRAILING_COMMA = re.compile(r',(?P<tail>\s*[]}])')

def relaxed_load(path, **kw):
    return relaxed_loads(open(path).read(), **kw)

def relaxed_loads(text, **kw):
    """
    Parse and return a relaxed JSON string.
    """
    ordered = kw.pop('ordered', False)
    if ordered:  kw['object_pairs_hook'] = OrderedDict
    # TODO: need a little state machine that performs the translation so that
    # TODO: line and column numbers are preserved, and so that we can have
    # TODO: http:// in a string (instead of it being treated like a comment). 
    #print "== raw text\n", text
    text = _LINE_CONTINUATION.sub('', text)
    #print "== joined lines\n", text
    text = _TRAILING_COMMENT.sub(r'\n', text)
    text = _MULTILINE_COMMENT.sub(r'', text)
    #print "== stripped comments\n", text
    text = _LEADING_TEXT.sub('{', text)
    #print "== trimmed text\n", text
    text = _UNQUOTED_FIELDNAME.sub(r'\g<prefix>"\g<key>"\g<tail>', text)
    #print "== quoted field names\n", text
    text = _TRAILING_COMMA.sub(r'\g<tail>', text)
    #print "== processed text\n", text
    try:
        obj = json.loads(text, object_hook=decode_dict_as_str, **kw)
    except ValueError, e:
        msg = [str(e)]
        M = re.findall('line ([0-9]*) column ([0-9]*)', msg[0])
        if M:
            line,col = int(M[0][0]), int(M[0][1])
            lines = text.split("\n")
            if line>=2: msg.append(lines[line-2])
            if line>=1: msg.append(lines[line-1])
            msg.append(" "*(col-1) + "^")
            if line<len(lines): msg.append(lines[line])
        msg = "\n".join(msg)
        raise e.__class__(msg)
    return obj

@contextmanager
def float_format(formatstr='.15g'):
    """
    Allow the float format to be changed for a json encoding action.

    This is a context manager, and should be used for example as::

        >>> with float_format('.2g'):
        >>>    print json.dumps(sqrt(2))
        1.41 
    """
    formatter = json.encoder.FLOAT_REPR
    json.encoder.FLOAT_REPR = lambda o: format(o, formatstr)
    yield
    json.encoder.FLOAT_REPR = formatter

def numpy_encoder(o):
    """
    JSON encoder for numpy data.

    To automatically convert numpy data to lists when writing a datastream
    use json.dumps(object, default=numpy_json).
    """
    try:
        return o.tolist()
    except AttributeError:
        raise TypeError


def _decode_list(lst):
    newlist = []
    for i in lst:
        if isinstance(i, unicode):
            i = i.encode('utf-8')
        elif isinstance(i, list):
            i = _decode_list(i)
        newlist.append(i)
    return newlist

def decode_dict_as_str(dct):
    newdict = {}
    for k, v in dct.iteritems():
        if isinstance(k, unicode):
            k = k.encode('utf-8')
        if isinstance(v, unicode):
            v = v.encode('utf-8')
        elif isinstance(v, list):
            v = _decode_list(v)
        newdict[k] = v
    return newdict


def test():
    """
    Verify that the translation from pseudo-JSON to JSON works.
    """
    good = """\
// This is a source definition with no errors
var entry = {
field : { // A comment about the field
  "field" : "te\\
x\\  
t",
  other$field : 56,
  },
/*
multiline comment
*/
secondfield : {
  content: ["string", "string"], /* a second comment */
  content: [{name:"good", value:3, URL:"http:\\/\\/my.url.com"},]
  },
}
"""
    broken = """\
// This is a source definition with a missing comma
{
field : { // A comment about the field
  field : "te\\   
x\\  
t"
  other$field : 56,
  },
/*
multiline comment
*/
secondfield : {
  content: ["string", "string"], /* a second comment */
  content: [{name:"good", value:3},]
  },
}
"""
    result = relaxed_loads(good)
    assert result['field']['field'] == "text"
    assert result['field']['other$field'] == 56
    assert result['secondfield']['content'][0]['name'] == 'good'
    try: relaxed_loads(broken)
    except ValueError, _: pass
    else: raise Exception("No exception raised in broken")

if __name__ == "__main__":
    test()
