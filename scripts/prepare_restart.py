#!/usr/bin/env python3
"""Recover local archived QC and freeze the completed samples for a restart."""
import argparse
import json
import os
from pathlib import Path
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--workdir', required=True)
parser.add_argument('--replace-newer', action='store_true',
                    help='Also restore archived stats overwritten after sample completion')
parser.add_argument('--backup-dir', type=Path)
parser.add_argument('--rebuild-samples', default='')
parser.add_argument('--caller', choices=['Deepvariant', 'HaplotypeCaller', 'BOTH'], default='Deepvariant')
parser.add_argument('--chrm', choices=['Yes', 'No'], default='Yes')
args = parser.parse_args()
os.chdir(args.workdir)
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import SAMPLEINFO, level0_regions, level1_regions
from restart_state import prepare

path, state = prepare(Path.cwd(), SAMPLEINFO,
                      {'restart_rebuild_samples':args.rebuild_samples,
                       'caller':args.caller, 'chrM':args.chrm},
                      (level0_regions, level1_regions),
                      replace_newer=args.replace_newer, backup_root=args.backup_dir)
print(json.dumps({'manifest':str(path), 'reused_samples':len(state['reused_samples']),
                  'restored_files':sum(map(len,state['restored_files'].values()))}), flush=True)
