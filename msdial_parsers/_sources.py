"""Private helpers for parser input sources."""

from contextlib import contextmanager


@contextmanager
def open_binary_source(source):
    """Open a path in binary mode or borrow a caller-owned binary stream."""
    if hasattr(source, "read"):
        yield source
    else:
        with open(source, "rb") as stream:
            yield stream
