"""Validated conditional-channel summaries and shared publication figures."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[1]
PAPER=ROOT/'research/paper/channels-2026-10-27'
FAMILIES=['uniform_magnitude','beta2_5_magnitude'];BANKS=['pure_bbh','official_mixture'];POLICIES=['matched_filter','with_cwb'];EPS=[0,.001,.01]

def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    fits=json.loads((BASE/'branching-results.json').read_text());mc=json.loads((BASE/'observed-mcmc-validation.json').read_text())
    products=json.loads((BASE/'channel-product-validation.json').read_text());components=json.loads((BASE/'detected-component-results.json').read_text())
    assert len(fits['cases'])==24 and mc['status']=='pass' and products['status']=='pass' and components['status']=='computed'
    for report in [mc,products,components]:
        for name,digest in report['input_sha256'].items():assert sha(BASE/name)==digest,name
    channels=fits['channels']
    def key(f,b,p,e):return f'{f}:{b}:{p}:epsilon={e:g}'
    def case(f,b,p,e):return fits['cases'][key(f,b,p,e)]
    def med(f,b,p,e):return np.array(case(f,b,p,e)['intrinsic_quantiles'])[:,1]
    contrasts={name:[] for name in ['regularization','hidden_spin','bank','search_policy']}
    for f in FAMILIES:
        for b in BANKS:
            for p in POLICIES:
                values=np.array([med(f,b,p,e) for e in EPS]);span=np.ptp(values,axis=0)
                for c in range(5):contrasts['regularization'].append(dict(channel=channels[c],value_pp=float(100*span[c]),family=f,bank=b,policy=p))
    for b in BANKS:
        for p in POLICIES:
            for e in EPS:
                delta=100*(med(FAMILIES[1],b,p,e)-med(FAMILIES[0],b,p,e))
                for c in range(5):contrasts['hidden_spin'].append(dict(channel=channels[c],value_pp=float(abs(delta[c])),signed_value_pp=float(delta[c]),bank=b,policy=p,epsilon=e))
    for f in FAMILIES:
        for p in POLICIES:
            for e in EPS:
                delta=100*(med(f,BANKS[1],p,e)-med(f,BANKS[0],p,e))
                for c in range(5):contrasts['bank'].append(dict(channel=channels[c],value_pp=float(abs(delta[c])),signed_value_pp=float(delta[c]),family=f,policy=p,epsilon=e))
        for b in BANKS:
            for e in EPS:
                delta=100*(med(f,b,POLICIES[1],e)-med(f,b,POLICIES[0],e))
                for c in range(5):contrasts['search_policy'].append(dict(channel=channels[c],value_pp=float(abs(delta[c])),signed_value_pp=float(delta[c]),family=f,bank=b,epsilon=e))
    maxima={name:max(rows,key=lambda r:r['value_pp']) for name,rows in contrasts.items()}
    all_rows=list(fits['cases'].values());median=np.array([r['intrinsic_quantiles'] for r in all_rows])[:,:,1]
    synthesis=dict(schema='blackhole-conditional-channel-synthesis/1',status='conditional_results_with_explicit_remaining_gaps',
        events=69,posterior_rows=11180578,cases=24,channels=channels,context=fits['context'],
        scientific_statement='Within the fixed-context restricted-domain emulator experiment, intrinsic branch medians depend more strongly on regularization and declared hidden-spin completion than on the tested bank/search replacements. This is not a measurement of cosmic formation fractions or general robustness.',
        contrast_maxima=maxima,contrasts=contrasts,ce_median_range=[float(median[:,0].min()),float(median[:,0].max())],
        maximum_qmc_replicate_quantile_span_pp=float(100*max(r['replicate_quantile_span'] for r in all_rows)),
        maximum_independent_mcmc_quantile_difference_pp=float(100*max(r['maximum_qmc_quantile_difference'] for r in mc['cases'].values())),
        maximum_mcmc_rhat=max(r['maximum_rhat'] for r in mc['cases'].values()),
        minimum_mcmc_bulk_ess=min(r['minimum_bulk_ess'] for r in mc['cases'].values()),
        minimum_mcmc_tail_ess=min(r['minimum_tail_ess'] for r in mc['cases'].values()),
        maximum_normalizer_half_quantile_shift_pp=float(100*max(r['maximum_normalizer_half_quantile_shift'] for r in all_rows)),
        minimum_event_concentration=min(r['concentration']['minimum_event_weight_concentration'] for r in all_rows),
        minimum_selection_concentration=min(r['concentration']['selection_weight_concentration'] for r in all_rows),
        maximum_decimal_log_product_error=max(r['max_log_product_error'] for r in products['cases']),
        limits=fits['limits']+['Low weight concentration is retained: fixed-bank channel posteriors are not calibrated campaign or PE-chain estimates.','Regularized detected labels include a shared artificial density term; the decomposition must not be reported as observed astrophysical membership.'],
        input_sha256={p.name:sha(p) for p in [Path(__file__),BASE/'branching-results.json',BASE/'observed-mcmc-validation.json',BASE/'channel-product-validation.json',BASE/'detected-component-results.json',BASE/'channel-log-normalization.json',BASE/'likelihood-design.json']})
    (BASE/'branching-synthesis.json').write_text(json.dumps(synthesis,indent=2)+'\n',encoding='utf-8')
    PAPER.mkdir(parents=True,exist_ok=True);(PAPER/'figures').mkdir(exist_ok=True)
    facts=dict(MaxReg=maxima['regularization']['value_pp'],MaxHidden=maxima['hidden_spin']['value_pp'],MaxBank=maxima['bank']['value_pp'],MaxPolicy=maxima['search_policy']['value_pp'],MaxIndependent=synthesis['maximum_independent_mcmc_quantile_difference_pp'],MaxReplicate=synthesis['maximum_qmc_replicate_quantile_span_pp'],MaxNormHalf=synthesis['maximum_normalizer_half_quantile_shift_pp'],MinEvent=synthesis['minimum_event_concentration'],MinSelection=synthesis['minimum_selection_concentration'],CEMin=100*median[:,0].min(),CEMax=100*median[:,0].max())
    (PAPER/'facts.tex').write_text('\n'.join('\\newcommand{\\'+name+'}{'+(f'{value:.3f}' if name=='MaxNormHalf' else f'{value:.2f}')+'}' for name,value in facts.items())+'\n',encoding='utf-8')
    for lang in ['en','zh']:
        header=r'Base $Q$ & $\epsilon$ & $\beta_{\rm CE}$ [5,50,95]\% & $\eta_{\rm CE}^{\rm label}$\% & $w_U$\% \\'
        if lang=='zh':header=header.replace('Base $Q$','基准 $Q$')
        rows=[]
        for f,label in zip(FAMILIES,['U','B']):
            for e in EPS:
                r=case(f,'official_mixture','matched_filter',e);q=100*np.array(r['intrinsic_quantiles'])[0]
                eta=100*np.array(r['detected_quantiles'])[0,1]
                u=100*np.array(components['cases'][key(f,'official_mixture','matched_filter',e)]['uniform_term_quantiles'])[1]
                rows.append(f'{label} & {e:g} & [{q[0]:.2f}, {q[1]:.2f}, {q[2]:.2f}] & {eta:.2f} & {u:.2f} '+r'\\')
        table='\\begin{tabular}{llrrr}\n\\toprule\n'+header+'\n\\midrule\n'+'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}\n'
        (PAPER/f'table-main-{lang}.tex').write_text(table,encoding='utf-8')
    plt.rcParams.update({'font.size':8.5,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    colors=['#333f53','#2765a6','#ba5825'];markers=['o','D']
    fig,(left,right)=plt.subplots(1,2,figsize=(6.5,5.0),gridspec_kw={'width_ratios':[1,1.25]})
    labels=[]
    for fi,f in enumerate(FAMILIES):
        for ei,e in enumerate(EPS):
            j=fi*3+ei;q=100*np.array(case(f,'official_mixture','matched_filter',e)['intrinsic_quantiles'])
            left.errorbar(q[0,1],5-j,xerr=np.array([[q[0,1]-q[0,0]],[q[0,2]-q[0,1]]]),fmt=markers[fi],color=colors[ei],capsize=3)
            labels.append(f'$Q_{"U" if fi==0 else "B"},\\ \\epsilon={e:g}$')
            for ci in range(1,5):
                y=4-ci+(2.5-j)*.10
                right.errorbar(q[ci,1],y,xerr=np.array([[q[ci,1]-q[ci,0]],[q[ci,2]-q[ci,1]]]),fmt=markers[fi],color=colors[ei],capsize=2,ms=4)
    left.set(yticks=np.arange(6),yticklabels=labels[::-1],xlim=(50,96),ylim=(-.6,5.6),xlabel=r'Intrinsic CE weight $\beta_{\rm CE}$ / %',title='A  Common envelope')
    right.set(yticks=[0,1,2,3],yticklabels=['SMT','NSC','GC','CHE'],xlim=(0,25),ylim=(-.5,3.5),xlabel=r'Intrinsic weight $\beta_c$ / %',title='B  Other channels')
    for ax in [left,right]:ax.grid(axis='x',alpha=.18)
    fig.suptitle('Conditional branch weights: birth spin 0.1, CE efficiency 2',fontsize=10.5)
    fig.text(.5,.025,'Points and bars: marginal medians and model-conditional 90% intervals.\nOfficial mixture; matched-filter searches; 69 events; not cosmic rates.\nMedians need not sum to1. Q_U/Q_B denote conditional base laws.',ha='center',fontsize=7.5)
    fig.subplots_adjust(left=.18,right=.98,top=.86,bottom=.23,wspace=.35)
    for ext in ['svg','png','pdf']:fig.savefig((BASE if ext!='pdf' else PAPER/'figures')/f'channel-posteriors.{ext}',dpi=160)
    plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(6.5,4.1),gridspec_kw={'width_ratios':[1.05,1,1]})
    axes[0].barh(np.arange(4),[maxima[k]['value_pp'] for k in ['search_policy','bank','hidden_spin','regularization']],color=['#8ca0b5','#71869d','#2765a6','#ba5825'])
    axes[0].set(yticks=np.arange(4),yticklabels=['Search policy','O3 bank','Hidden-spin law','Regularizer'],xlabel='Largest median shift / pp',title='A  Specified grid\nsensitivity')
    for fi,f in enumerate(FAMILIES):
        x=100*np.array(EPS);eff=[];u=[];bounds=[]
        for e in EPS:
            r=case(f,'official_mixture','matched_filter',e);eff.append(100*r['efficiency'][0])
            uq=100*np.array(components['cases'][key(f,'official_mixture','matched_filter',e)]['uniform_term_quantiles']);u.append(uq[1]);bounds.append([uq[1]-uq[0],uq[2]-uq[1]])
        color=['#2765a6','#ba5825'][fi];label=['Q_U','Q_B'][fi]
        axes[1].plot(x,eff,markers[fi]+'-',color=color,label=label)
        axes[2].errorbar(x+(fi-.5)*.006,u,yerr=np.array(bounds).T,fmt=markers[fi]+'-',color=color,capsize=3,label=label)
    for ax in axes[1:]:
        ax.set(xticks=[0,.1,1],xticklabels=['0','0.1','1'],xlabel=r'Regularizer $\epsilon$ / %')
        plt.setp(ax.get_xticklabels(),rotation=45,ha='right')
        ax.legend(frameon=False,fontsize=7.5);ax.grid(alpha=.15)
    axes[1].set(ylabel='Computed CE efficiency / %',title='B  Low channel\nefficiency')
    axes[2].set(ylabel='Uniform term / %',title='C  Shared term in\ndetected mixture')
    fig.text(.5,.025,'A: maxima over specified comparisons, not uncertainty bars.\nC: conditional 90% intervals for a mathematical density term.\nB-C: fixed domain, official mixture, matched-filter searches; no campaign calibration.',ha='center',fontsize=7)
    fig.subplots_adjust(left=.18,right=.98,top=.82,bottom=.26,wspace=.65)
    for ext in ['svg','png','pdf']:fig.savefig((BASE if ext!='pdf' else PAPER/'figures')/f'channel-sensitivity.{ext}',dpi=160)
    plt.close(fig)
    manifest=dict(schema='blackhole-channel-paper-evidence/1',synthesis_sha256=sha(BASE/'branching-synthesis.json'),facts=facts,
                  source_sha256=synthesis['input_sha256'],generated_sha256={p.relative_to(PAPER).as_posix():sha(p) for p in [PAPER/'facts.tex',PAPER/'table-main-en.tex',PAPER/'table-main-zh.tex',*sorted((PAPER/'figures').glob('*.pdf'))]})
    (PAPER/'evidence-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(maxima=maxima,ce_median_range=synthesis['ce_median_range'],independent_error_pp=synthesis['maximum_independent_mcmc_quantile_difference_pp']),indent=2))


if __name__=='__main__':main()
