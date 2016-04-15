"""
Utility to add annotations to python exceptions.
"""
import sys

# Platform cruft: WindowsError is only defined on Windows.
try:
    WindowsError
except NameError:
    class WindowsError(Exception):
        """
        Fake WindowsException when not on Windows.
        """
        pass

def annotate_exception(msg, exc=None):
    """
    Add an annotation to the current exception, which can then be forwarded
    to the caller using a bare "raise" statement to raise the annotated
    exception.  If the exception *exc* is provided, then that exception is the
    one that is annotated, otherwise *sys.exc_info* is used.

    Example::

        >>> D = {}
        >>> try:
        ...    print(D['hello'])
        ... except:
        ...    annotate_exception("while accessing 'D'")
        ...    raise
        Traceback (most recent call last):
            ...
        KeyError: "hello while accessing 'D'"
    """
    if not exc:
        exc = sys.exc_info()[1]

    # Can't extend WindowsError exceptions; instead raise a new exception.
    # TODO: try to incorporate current stack trace in the raised exception
    if isinstance(exc, WindowsError):
        raise OSError(str(exc) + " " + msg)

    args = exc.args
    if not args:
        exc.args = (msg,)
    else:
        try:
            arg0 = " ".join((args[0], msg))
            exc.args = tuple([arg0] + list(args[1:]))
        except:
            exc.args = (" ".join((str(exc), msg)),)
