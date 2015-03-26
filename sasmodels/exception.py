"""
Utility to add annotations to python exceptions.
"""

def annotate_exception(exc, msg):
    """
    Add an annotation to the current exception, which can then be forwarded
    to the caller using a bare "raise" statement to reraise the annotated
    exception.
    Example::
        >>> D = {}
        >>> try:
        ...    print D['hello']
        ... except Exception,exc:
        ...    annotate_exception(exc, "while accessing 'D'")
        ...    raise
        Traceback (most recent call last):
            ...
        KeyError: "hello while accessing 'D'"
    """
    args = exc.args
    if not args:
        exc.args = (msg,)
    else:
        try:
            arg0 = " ".join((args[0],msg))
            exc.args = tuple([arg0] + list(args[1:]))
        except:
            exc.args = (" ".join((str(exc),msg)),)

