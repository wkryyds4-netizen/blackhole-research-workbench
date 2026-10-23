"""Pinned, read-only inventory; no repository code is executed here."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import inspect_projects

inspect_projects.ROOT=Path(__file__).resolve().parent/'generalization-snapshots'
REPOS=['kyleyhw/black_hole','samrht/blackhole-sim','gyoto/Gyoto',
       'sxs-collaboration/sxs','pesummary/pesummary','ColmTalbot/gwpopulation',
       'bilby-dev/bilby','einsteinpy/einsteinpy','duetosymmetry/qnm']
if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(inspect_projects.fetch,REPOS))
