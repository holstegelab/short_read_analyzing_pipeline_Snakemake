import os
import sys
import subprocess
import importlib
from typing import Tuple, Optional


def ensure_fastcheck_hts(logger=None) -> Tuple[bool, Optional[object]]:
    """
    Try to import local 'fastcheck_hts' C-extension from the scripts folder.
    If not found, build it in-place with setup_hts.py and retry.
    Returns (has_module: bool, module or None)
    """
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)

    try:
        import fastcheck_hts  # type: ignore
        return True, fastcheck_hts
    except ModuleNotFoundError:
        if logger:
            logger(f"fastcheck_hts not found; attempting local build in {here}")
        try:
            # Propagate conda hints to setup so pkg-config can find htslib
            env = os.environ.copy()
            conda_prefix = env.get('CONDA_PREFIX') or env.get('VIRTUAL_ENV')
            if conda_prefix:
                pkgp = os.path.join(conda_prefix, 'lib', 'pkgconfig')
                env['PKG_CONFIG_PATH'] = (pkgp + (':' + env['PKG_CONFIG_PATH'] if 'PKG_CONFIG_PATH' in env and env['PKG_CONFIG_PATH'] else ''))
                inc = os.path.join(conda_prefix, 'include')
                lib = os.path.join(conda_prefix, 'lib')
                env.setdefault('HTSLIB_INCLUDE', inc)
                env.setdefault('HTSLIB_LIBDIR', lib)
            subprocess.check_call([sys.executable, os.path.join(here, 'setup_hts.py'), 'build_ext', '--inplace'], cwd=here, env=env)
            importlib.invalidate_caches()
            import fastcheck_hts  # type: ignore
            if logger:
                logger('fastcheck_hts built in-place and loaded')
            return True, fastcheck_hts
        except Exception as e:
            if logger:
                logger(f'fastcheck_hts build failed; not using HTS path: {e}')
            return False, None
    except Exception as e:
        if logger:
            logger(f'fastcheck_hts import failed; not using HTS path: {e}')
        return False, None
