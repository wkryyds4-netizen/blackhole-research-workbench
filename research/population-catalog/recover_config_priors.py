"""Recover declared bounds and the source-supported LALInference default family."""
from pathlib import Path
import hashlib
import json
import math

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
SOURCES=ROOT/'research/population-selection-audit/sources'
PREFIX='gitlab-1438-af4bc21a9f92-'


def main():
    manifest={v['file']:v for v in json.loads((SOURCES/'manifest.json').read_text())}
    texts={};sources=[]
    for suffix in ['lalinference--lib--LALInferencePrior.c','lalinference--lib--LALInferenceInitCBC.c']:
        path=SOURCES/(PREFIX+suffix);raw=path.read_bytes();assert hashlib.sha256(raw).hexdigest()==manifest[path.name]['sha256']
        texts[suffix]=raw.decode();sources.append(dict(file=path.name,sha256=manifest[path.name]['sha256'],url=manifest[path.name]['url']))
    prior=texts['lalinference--lib--LALInferencePrior.c'];initial=texts['lalinference--lib--LALInferenceInitCBC.c']
    assert 'logPrior+=3.0* log_dist;' in prior and 'logPrior+=2.0*log(dist);' in prior
    assert 'logPrior+=log(m1*m1/mc);' in prior and 'logPrior+=log(m1*m1);' in prior
    assert 'Usual case has uniform in a' in prior and '--spin-volumetricprior' in prior
    assert 'REAL8 Dmin=1.0;' in initial and 'REAL8 a1min=0.0,a1max=1.0;' in initial and 'REAL8 qMax=1.0;' in initial
    configs=json.loads((BASE/'missing-prior-configurations.json').read_text());events=[]
    alternatives=['distance-prior-uniform','distance-prior-comoving-volume','spin-volumetricprior','alignedspin-zprior','disable-spin','noSpin','singleSpin','spinAligned','aligned-spin','fix-distance']
    for row in configs:
        selected=row['selected_prior_configuration'];paths=set(row['configuration_paths'])
        found=[flag for flag in alternatives if 'engine/'+flag in paths];assert not found,(row['event'],found)
        def value(key):
            v=selected['engine/'+key];assert isinstance(v,list) and len(v)==1;return float(v[0])
        bounds=dict(chirp_mass=dict(minimum=value('chirpmass-min'),maximum=value('chirpmass-max')),
            mass_1=dict(minimum=value('comp-min'),maximum=value('comp-max')),
            mass_2=dict(minimum=value('comp-min'),maximum=value('comp-max')),
            mass_ratio=dict(minimum=value('q-min'),maximum=1.),luminosity_distance=dict(minimum=1.,maximum=value('distance-max')),
            a_1=dict(minimum=0.,maximum=value('a_spin1-max')),a_2=dict(minimum=0.,maximum=value('a_spin2-max')),
            tilt_1=dict(minimum=0.,maximum=math.pi),tilt_2=dict(minimum=0.,maximum=math.pi))
        assert row['metadata']['approximant']==['IMRPhenomXPHM']
        configured=selected['engine/approx'];assert configured==['IMRPhenomXPHMpseudoFourPN']
        events.append(dict(index=row['index'],event=row['event'],bounds=bounds,configured_approximant=configured,
            standardized_metadata_approximant=row['metadata']['approximant'],
            waveform_name_scope='Preserve the configured pseudoFourPN suffix; no assertion that all waveform settings match other events.',
            explicit_bounds=['chirp_mass_min_max','component_mass_min_max','mass_ratio_min','spin_magnitude_max','luminosity_distance_max'],
            documented_defaults=['mass_ratio_max=1','spin_magnitude_min=0','full_isotropic_tilts','distance_min=1 Mpc','flat_detector_component_masses','distance_squared_density','uniform_spin_magnitudes'],
            incompatible_prior_flags_present=found,source_supported_default_family=True,
            complete_historical_runtime_verified=False,
            scope='Official released engine configuration plus fixed historical source defaults; not Bilby analytic strings or complete historical command-line provenance. Retained-row checks are performed separately.'))
    assert len(events)==6
    report=dict(schema='blackhole-lalinference-config-prior-recovery/1',source_commit='af4bc21a9f9246ef5403052ec3b66ed2b52d90ed',verified_sources=sources,events=events,
        coordinate_density_up_to_event_constant='D_L(z)^2 * dD_L/dz * m1_source * (1+z)^2 / 4',
        mass_jacobian='Historical prior is m1_detector^2/Mchirp in (Mchirp,q), or m1_detector^2 in (log Mchirp,q); transforming to (m1_source,q,z) gives m1_source*(1+z)^2.',
        limitations='Configuration/default-family recovery does not constrain likelihood outside support, recover unrecorded overrides, establish every historic build, or fix common spin references.')
    (BASE/'configuration-prior-recovery.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(dict(events=[dict(index=e['index'],event=e['event'],bounds=e['bounds']) for e in events],density=report['coordinate_density_up_to_event_constant']),indent=2))


if __name__=='__main__':main()
