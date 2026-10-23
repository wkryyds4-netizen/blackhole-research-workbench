"""Audit existing population inputs without modifying selection or posteriors."""
from pathlib import Path
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'research/population-audit'

def main():
    events=json.loads((ROOT/'data/sample_audit.json').read_text(encoding='utf-8'))
    result=json.loads((ROOT/'results/model.json').read_text(encoding='utf-8'))
    precision=json.loads((ROOT/'results/generalization-summary.json').read_text(encoding='utf-8'))
    selected=[e for e in events if e['included']]
    baseline=result['scenarios'][result['default']]
    with np.load(ROOT/'data/injections.npz',allow_pickle=False) as z:
        fields={k:{'shape':list(z[k].shape),'dtype':str(z[k].dtype)} for k in z.files}
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as z:
        event_fields=sorted(k.removeprefix('0_') for k in z.files if k.startswith('0_'))
    report={
        'schema':'blackhole-population-input-audit/1','date':'2026-10-23',
        'status':'exploratory_only; provenance and selection audit incomplete',
        'event_count':len(events),'selected_count':len(selected),
        'selected_missing_far':sum(e.get('far') is None for e in selected),
        'selected_catalogs':{c:sum(e['catalog']==c for e in selected) for c in sorted({e['catalog'] for e in selected})},
        'identity_method':'global assignment on posterior summary similarity plus chronology; not original event identifiers',
        'identity_match_max_cost':max(e['match_cost'] for e in selected),
        'posterior_array_fields':event_fields,'injection_array_fields':fields,
        'baseline':{'mixture_definition':'phenomenological aligned-tilt fraction; not formation-channel fraction','corrected_quantiles_05_50_95':baseline['corrected'],'uncorrected_quantiles_05_50_95':baseline['uncorrected'],'posterior_mass_passing_exploratory_mc_thresholds':precision['baseline_posterior_mass_passing']},
        'blockers':[
            {'id':'event_identity','finding':'Event labels inferred from numerical summaries','required':'Recover identifiers from original PE files or verified preprocessing manifest'},
            {'id':'selection','finding':'Current code accepts missing FAR and uses posterior median secondary mass >3; injection search thresholds/run labels absent','required':'Specify event detection/classification policy and reconstruct matching injection selection by observing run; model any data-dependent classification cut'},
            {'id':'prior_measure','finding':'A numeric prior column exists but original PE analysis/prior/Jacobian metadata are not preserved','required':'Verify prior density in the same source-mass, q, redshift and spin coordinates used in population weights'},
            {'id':'mc_precision','finding':'Only 86.7 percent of baseline posterior mass passes current exploratory MC screen','required':'Increase effective posterior/injection sample coverage and verify numerical stability over posterior support'},
            {'id':'physical_interpretation','finding':'Two tilt components are not uniquely tied to formation mechanisms','required':'Test identifiable population features before adding forward formation-model predictions'}
        ],
        'sources':{'tutorial':'https://colmtalbot.github.io/gwpopulation/examples/GWTC3.html','population_paper':'https://arxiv.org/abs/2111.03634'},
        'input_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['data/posteriors.npz','data/injections.npz','data/sample_audit.json','results/model.json']},
        'changes_to_analysis':'none; existing inference preserved'
    }
    OUT.mkdir(exist_ok=True)
    (OUT/'input-audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['event_count','selected_count','selected_missing_far','selected_catalogs','baseline']},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
