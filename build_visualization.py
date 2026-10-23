from pathlib import Path
import json
import argparse

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,default=ROOT/'results',help='Output directory (default: repository results directory)')
OUT=parser.parse_args().output
def compact(x):
    if isinstance(x,float): return float(f'{x:.9g}')
    if isinstance(x,list): return [compact(v) for v in x]
    if isinstance(x,dict): return {k:compact(v) for k,v in x.items()}
    return x
data=json.loads((ROOT/'results/model.json').read_text(encoding='utf-8'))
template=(ROOT/'visualization.html').read_text(encoding='utf-8')
output=template.replace('__MODEL_DATA__',json.dumps(compact(data),ensure_ascii=False,separators=(',',':')))
assert '__MODEL_DATA__' not in output
assert len(output.encode('utf-8'))<1_000_000
OUT.mkdir(exist_ok=True,parents=True)
dest=OUT/'black-hole-robustness.html'
dest.write_text(output,encoding='utf-8')
print(dest, len(output.encode('utf-8')))
