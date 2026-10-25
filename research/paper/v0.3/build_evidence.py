"""Build bilingual numerical tables and publication-size vector figures.

Read validated research products only. No population fit or original data is
changed; v0.2 manuscript and build products are deliberately untouched.
"""
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import expit

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
RESEARCH=ROOT/'research'
INV=RESEARCH/'population-invariant'
CAT=RESEARCH/'population-catalog'


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def read(path):return json.loads(path.read_text(encoding='utf-8'))


def table(headers,rows,columns):
    return '\n'.join([r'\begin{tabular}{'+columns+'}',r'\toprule',' & '.join(headers)+r' \\',r'\midrule',
        *[' & '.join(row)+r' \\' for row in rows],r'\bottomrule',r'\end{tabular}',''])


def main():
    files={
        'native':CAT/'native-catalog-results.json','native_validation':CAT/'native-catalog-validation.json',
        'control':INV/'control-results.json','control_validation':INV/'control-validation.json',
        'o3':INV/'paired-o3-results.json','o3_validation':INV/'paired-o3-validation.json',
        'paired':INV/'paired-pe-results.json','paired_validation':INV/'paired-pe-validation.json',
        'support':INV/'support-results.json','support_validation':INV/'support-validation.json',
        'tail_plan':INV/'tail-validation-plan.json','milestone':INV/'finite-bank-milestone.json',
        'alias_validation':CAT/'alias-heldout-validation.json','historical_defaults':CAT/'historical-defaults-validation.json',
        'transport':INV/'reference-invariance-summary.json'}
    data={k:read(path) for k,path in files.items()}
    for key in ['native','control','o3','paired','support']:
        assert data[key+'_validation']['result_sha256']==sha(files[key])
    assert data['native']['completed_full_source_events']==69
    assert data['paired_validation']['original_step_derivative_gate_passed'] is False
    m=data['milestone']
    for name,digest in m['file_sha256'].items():assert sha(INV/name)==digest,name
    for name,digest in m['frozen_hashes'].items():assert sha(ROOT/name)==digest,name
    before={str(path):sha(path) for path in files.values()}
    rows={r['name']:r for r in data['paired']['contrasts']}
    facts={
        'DataEvents':'69','DataRows':'11,180,578','DataBytes':'23,554,981,402',
        'DataReplicates':'512','DataEventDraw':'6,913','DataOThreeRows':'284,672','DataOThreeDraw':'128,454',
        'DataNativeLowGap':f"{data['native']['width_gaps_pp'][0]:.5f}",
        'DataNativeBaseGap':f"{data['native']['width_gaps_pp'][1]:.5f}",
        'DataMaxSourceShift':f"{max(abs(x) for x in data['native']['median_shift_pp']):.5f}",
        'DataJointLowNonpositive':f"{100*rows['matched:native_width:low_spin']['joint_nonpositive_fraction']:.2f}",
        'DataOThreeLowNonpositive':f"{100*next(r['nonpositive_fraction'] for r in data['o3']['contrasts'] if r['name']=='matched:native_width:low_spin'):.2f}",
        'DataInteraction':f"{rows['matched:native_width:low_spin']['interaction_rms_pp']:.5f}",
        'DataTailMass':f"{100*data['support']['known_o3_omitted_domain'][0]['population_redshift_tail_mass']:.4f}",
        'DataTailBound':f"{data['support']['known_o3_omitted_domain'][0]['o3_time_weighted_selection_contribution_upper']:.5f}",
        'DataCompletionMax':f"{m['maximum_matched_infinite_completion_shift_pp']:.2f}",
        'DataTailPlan':f"{data['tail_plan']['total_direct_draws_if_all_zero']:,}",
        'DataEventZeroRows':f"{data['support']['completion_event']['rows']:,}"}
    # General counters use commas deliberately; no locale-dependent formatting.
    fact_source='\n'.join('\\newcommand{\\'+key+'}{'+value+'}' for key,value in facts.items())+'\n'
    table_outputs={}
    for lang in ['en','zh']:
        project=BASE/f'latex-{lang}';assert (project/'main.tex').exists()
        (project/'facts.tex').write_text(fact_source,encoding='utf-8')
        header=['Magnitude prior',r'$\lambda$',r'$g$ median [5th, 95th] (\%)'] if lang=='en' else ['幅度先验',r'$\lambda$',r'$g$ 中位数[5\%,95\%分位数]（\%）']
        q=np.array(data['control']['scenarios']['matched_filter']['quantiles'])*100
        values=[]
        for k in range(6):
            a,b=[(1.5,8.5),(2,5)][k//3]
            values.append([rf'$\mathrm{{Beta}}({a},{b})$',str([0,4,16][k%3]),f'{q[k,1]:.2f} [{q[k,0]:.2f}, {q[k,2]:.2f}]'])
        table_outputs[lang,'g-table.tex']=table(header,values,'lrr')
        names=['matched:native_width:low_spin','matched:native_width:baseline_spin',*[f'matched:effective_amplitude:lambda{l}' for l in [0,4,16]]]
        labels=([r'$\Delta r$, low spin',r'$\Delta r$, baseline spin',r'$\Delta g$, $\lambda=0$',r'$\Delta g$, $\lambda=4$',r'$\Delta g$, $\lambda=16$'] if lang=='en' else
                [r'$\Delta r$，低自旋',r'$\Delta r$，基线幅度',r'$\Delta g$，$\lambda=0$',r'$\Delta g$，$\lambda=4$',r'$\Delta g$，$\lambda=16$'])
        values=[]
        for name,label in zip(names,labels):
            row=rows[name];vals=[label,f"{row['full_bank_value_pp']:.3f}"]
            for field in ['event_only_5_50_95_pp','selection_only_5_50_95_pp','joint_5_50_95_pp']:
                a,_,b=row[field];vals.append(f'[{a:.3f}, {b:.3f}]')
            values.append(vals)
        header=(['Contrast','Full banks','Event only','O3 only','Joint'] if lang=='en' else ['对比','全银行','仅事件','仅O3','联合'])
        table_outputs[lang,'finite-table.tex']=table(header,values,'lrrrr')
        values=[]
        for lam in [0,4,16]:
            a,b=[next(x for x in data['support']['ordering_certificate_thresholds'] if x['policy']==p and x['lambda_shape']==lam) for p in [0,1]]
            values.append([str(lam),f"{a['critical_log_correction_oscillation']:.5f}",f"{b['critical_log_correction_oscillation']:.5f}",f"{100*a['uniform_relative_selection_completion_bound_at_threshold']:.3f}"])
        header=([r'$\lambda$',r'$B_*$ matched',r'$B_*$ with cWB',r'Equivalent bound (\%)'] if lang=='en' else [r'$\lambda$',r'$B_*$ 匹配滤波',r'$B_*$ 含cWB',r'等效上界（\%）'])
        table_outputs[lang,'support-table.tex']=table(header,values,'rrrr')
        values=[]
        for label in ['0.0','0.01','1.0','infinity_limit']:
            row=next(r for r in data['support']['event_completions'] if r['eta_label']==label)
            shifts=row['median_shift_pp'][0];gap=row['amplitude_median_difference_pp'][0]
            eta=r'$\infty$' if label=='infinity_limit' else label
            values.append([eta,f'[{min(shifts):.3f}, {max(shifts):.3f}]',f'[{min(gap):.3f}, {max(gap):.3f}]'])
        header=([r'$\eta$','Median shifts across six models','Three amplitude gaps'] if lang=='en' else [r'$\eta$','六模型中位数位移范围','三项幅度差范围'])
        table_outputs[lang,'completion-table.tex']=table(header,values,'rrr')
    for (lang,name),content in table_outputs.items():(BASE/f'latex-{lang}'/name).write_text(content,encoding='utf-8')
    figures=BASE/'figures';figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'legend.fontsize':7,
        'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):
        fig.savefig(figures/(name+'.pdf'));fig.savefig(figures/(name+'.svg'));fig.savefig(figures/(name+'.png'),dpi=190);plt.close(fig)
    # Native mixed-reference diagnostic, including the five-source baseline.
    fig,axes=plt.subplots(1,2,figsize=(6.3,3.4),layout='constrained');y=np.arange(6)
    labels=[f'{a}, s={s}' for a in ['low','base'] for s in [.3,.6,1.2]]
    for quant,offset,color,label in [(data['native']['baseline_quantiles'],-.12,'.6','Five-source baseline'),(data['native']['quantiles'],.12,'#277e9e','All 69 full sources')]:
        q=np.array(quant)*100
        axes[0].errorbar(q[:,1],y+offset,xerr=[q[:,1]-q[:,0],q[:,2]-q[:,1]],fmt='o',ms=3,capsize=2,color=color,label=label)
    axes[0].set(yticks=y,yticklabels=labels,xlabel='Native-reference r (%)',title='Conditional 90% credible intervals')
    axes[0].set_ylim(6.8,-.5);axes[0].legend(frameon=False,loc='lower left')
    axes[1].barh(y,data['native']['median_shift_pp'],color=['#277e9e']*3+['#cc6b36']*3)
    axes[1].set(yticks=y,yticklabels=[str(s) for s in [.3,.6,1.2]*2],xlabel='Median shift (pp)',ylabel='Tilt width',title='Change after source restoration')
    axes[1].axvline(0,c='.4',lw=.7);axes[1].invert_yaxis()
    for ax in axes:ax.grid(axis='x',alpha=.15)
    save(fig,'native-restoration')
    with np.load(INV/'control-checkpoints.npz',allow_pickle=False) as cp:g=cp['g'];density=cp['densities'][0]
    fig,axes=plt.subplots(1,2,figsize=(6.3,2.8),sharey=True,layout='constrained')
    colors=['#277e9e','#57916b','#cc6b36']
    for a,ax in enumerate(axes):
        for k,lam in enumerate([0,4,16]):ax.plot(g,density[a*3+k],color=colors[k],label=rf'$\lambda={lam}$')
        ax.set(xlabel=r'$g=P(\chi_{\rm eff}<0)$',xlim=(0,.8),title=['Beta(1.5, 8.5)','Beta(2, 5)'][a]);ax.grid(alpha=.15)
    axes[0].set_ylabel('Conditional posterior density');axes[1].legend(frameon=False)
    save(fig,'effective-spin')
    fig,axes=plt.subplots(1,2,figsize=(6.3,3.2),layout='constrained')
    cases=[('event_only_5_50_95_pp',-.16,'#57916b','Event'),('selection_only_5_50_95_pp',0,'#277e9e','O3'),('joint_5_50_95_pp',.16,'#cc6b36','Joint')]
    for ax,names,labs,title in [(axes[0],names[:2],['Low spin','Baseline spin'],'Native width contrasts'),(axes[1],names[2:],[r'$\lambda=0$',r'$\lambda=4$',r'$\lambda=16$'],'Effective-spin amplitude contrasts')]:
        yy=np.arange(len(names))
        for field,offset,color,label in cases:
            q=np.array([rows[n][field] for n in names]);ax.errorbar(q[:,1],yy+offset,xerr=[q[:,1]-q[:,0],q[:,2]-q[:,1]],fmt='o',ms=3,capsize=2,color=color,label=label)
        ax.scatter([rows[n]['full_bank_value_pp'] for n in names],yy,c='k',s=15,marker='D',label='Full')
        ax.set(yticks=yy,yticklabels=labs,xlabel='Median difference (pp)',title=title);ax.invert_yaxis();ax.axvline(0,c='.4',lw=.6);ax.grid(axis='x',alpha=.15)
    axes[0].set_ylim(1.6,-.5);axes[0].legend(frameon=False,ncol=2,loc='lower left')
    save(fig,'finite-banks')
    with np.load(INV/'support-checkpoints.npz',allow_pickle=False) as cp:g=cp['g'];cdf=cp['cdf'][0]
    fig,ax=plt.subplots(figsize=(6.3,3.1),layout='constrained');B=np.linspace(0,1.5,301)
    for k,(lam,color) in enumerate(zip([0,4,16],colors)):
        lower=100*np.array([np.interp(expit(-b),cdf[k+3],g)-np.interp(expit(b),cdf[k],g) for b in B])
        upper=100*np.array([np.interp(expit(b),cdf[k+3],g)-np.interp(expit(-b),cdf[k],g) for b in B])
        ax.fill_between(B,lower,upper,color=color,alpha=.08);ax.plot(B,lower,color=color,label=rf'$\lambda={lam}$: lower bound')
    ax.axhline(0,c='.4',lw=.7);ax.set(xlabel='Assumed log-correction oscillation B',ylabel='Amplitude-model median difference (pp)',title='Conditional outer envelopes; B is not measured')
    ax.set_xlim(0,1.5);ax.set_xticks([0,.25,.5,.75,1,1.25,1.5])
    ax.legend(frameon=False,ncol=3,loc='upper left');ax.grid(alpha=.15);save(fig,'support-envelope')
    # Preserve the failed release-agreement gate alongside its historical replay.
    fig,ax=plt.subplots(figsize=(6.3,2.8),layout='constrained');x=np.arange(2);width=.24
    a=data['alias_validation']['summaries'];b=data['historical_defaults']['summaries']
    ax.bar(x-width,[r['max_target_cosine_difference'] for r in a],width,color='#cc6b36',label='Independent vs release')
    ax.bar(x,[r['max_final_tolerance_difference'] for r in a],width,color='#277e9e',label='Tolerance-pair change')
    ax.bar(x+width,[r['maximum_released_cosine_residual'] for r in b],width,color='#57916b',label='Historical defaults vs release')
    ax.axhline(1e-6,c='.3',ls='--',label='Original failed gate');ax.set(xticks=x,xticklabels=['16 Hz / Phenom','20 Hz / Mixed'],ylabel='Maximum cosine residual',ylim=(0,1.8e-6))
    ax.legend(frameon=False,ncol=2);ax.grid(axis='y',alpha=.15);save(fig,'reference-residual')
    for lang in ['en','zh']:
        out=BASE/f'latex-{lang}'/'figures';out.mkdir(exist_ok=True)
        for path in figures.glob('*.pdf'):shutil.copyfile(path,out/path.name)
    assert before=={str(path):sha(path) for path in files.values()},'Research products were modified'
    manifest=dict(schema='blackhole-bilingual-v03-evidence/1',date='2026-10-25',facts=facts,
        source_sha256={path.relative_to(ROOT).as_posix():sha(path) for path in files.values()},
        generated_sha256={p.relative_to(BASE).as_posix():sha(p) for p in sorted(BASE.glob('latex-*/*.tex')) if p.name.endswith('-table.tex') or p.name=='facts.tex'},
        figure_sha256={p.name:sha(p) for p in figures.glob('*.pdf')},producer_sha256=sha(Path(__file__)),
        scope='Bilingual fact/table identity and data-derived figures. Not an independent scientific validation or a peer review.')
    (BASE/'evidence-manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(dict(facts=facts,figures=len(manifest['figure_sha256'])),indent=2))


if __name__=='__main__':main()
