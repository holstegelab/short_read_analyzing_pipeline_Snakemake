import os
import sys
import subprocess
import importlib

def ensure_fastcheck(logger=None):
    """
    Try to import the local 'fastcheck' C-extension from the scripts folder.
    If not found, build it in-place with setup.py and retry.
    Returns (has_fastcheck: bool, fastcheck_module or None)
    """
    try:
        import fastcheck  # type: ignore
        return True, fastcheck
    except ModuleNotFoundError:
        here = os.path.dirname(os.path.abspath(__file__))
        if logger:
            logger(f"fastcheck module not found; attempting local build in {here}")
        try:
            subprocess.check_call([sys.executable, os.path.join(here, 'setup.py'), 'build_ext', '--inplace'], cwd=here)
            importlib.invalidate_caches()
            if here not in sys.path:
                sys.path.insert(0, here)
            import fastcheck  # type: ignore
            if logger:
                logger('fastcheck built in-place and loaded')
            return True, fastcheck
        except Exception as e:
            if logger:
                logger(f'fastcheck build failed; falling back to Python implementation: {e}')
            return False, None
    except Exception as e:
        if logger:
            logger(f'fastcheck import failed; falling back to Python implementation: {e}')
        return False, None
