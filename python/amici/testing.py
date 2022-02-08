"""Test support functions"""

import sys
from tempfile import TemporaryDirectory


class TemporaryDirectoryWinSafe(TemporaryDirectory):
    """TemporaryDirectory that will not raise if cleanup fails.

    If any extension was loaded from the temporary directory, cleanup would
    otherwise fail on Windows with a ``PermissionError``. This class ignores
    such failures.
    """
    def cleanup(self):
        try:
            super().cleanup()
        except PermissionError as e:
            if sys.platform not in {'win32', 'cygwin'}:
                raise e
        except NotADirectoryError:
            # Ignore exception on Windows for pyd files:
            #  NotADirectoryError: [WinError 267] The directory name is
            #  invalid: '....pyd'
            pass
