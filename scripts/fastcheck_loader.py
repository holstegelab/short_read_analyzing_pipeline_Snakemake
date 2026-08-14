import os
import sys
import subprocess
import importlib
import sysconfig


def _build_and_import(here, logger=None):
    subprocess.check_call(
        [
            sys.executable,
            os.path.join(here, 'setup.py'),
            'build_ext',
            '--inplace',
            '--force',
        ],
        cwd=here,
    )
    importlib.invalidate_caches()
    sys.modules.pop('fastcheck', None)
    if here not in sys.path:
        sys.path.insert(0, here)
    import fastcheck  # type: ignore
    if logger:
        logger('fastcheck built in-place and loaded')
    return fastcheck


def _local_extension_is_stale(here):
    extension_suffix = sysconfig.get_config_var('EXT_SUFFIX')
    if not extension_suffix:
        return False
    source = os.path.join(here, 'fastcheck.c')
    extension = os.path.join(here, f'fastcheck{extension_suffix}')
    return os.path.exists(extension) and os.path.getmtime(source) > os.path.getmtime(extension)


def ensure_fastcheck(logger=None):
    """
    Try to import the local 'fastcheck' C-extension from the scripts folder.
    If not found, build it in-place with setup.py and retry.
    Returns (has_fastcheck: bool, fastcheck_module or None)
    """
    here = os.path.dirname(os.path.abspath(__file__))
    if _local_extension_is_stale(here):
        if logger:
            logger('fastcheck C source is newer than the local extension; rebuilding')
        try:
            return True, _build_and_import(here, logger=logger)
        except Exception as e:
            if logger:
                logger(f'fastcheck rebuild failed; falling back to Python implementation: {e}')
            return False, None

    try:
        import fastcheck  # type: ignore
        return True, fastcheck
    except ModuleNotFoundError:
        if logger:
            logger(f"fastcheck module not found; attempting local build in {here}")
        try:
            return True, _build_and_import(here, logger=logger)
        except Exception as e:
            if logger:
                logger(f'fastcheck build failed; falling back to Python implementation: {e}')
            return False, None
    except Exception as e:
        if logger:
            logger(f'fastcheck import failed; falling back to Python implementation: {e}')
        return False, None
