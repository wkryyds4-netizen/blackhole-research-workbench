"""Download GWOSC GW150914 v3 strain and prepare an auditable display product."""
from pathlib import Path
import gzip, hashlib, io, json
import numpy as np
import requests
from scipy import signal

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
FS = 4096
GPS = 1126259462.4

def whiten_filter(x, psd):
    x = signal.detrend(x) * signal.windows.tukey(len(x), .15)
    w = np.fft.irfft(np.fft.rfft(x) / np.sqrt(psd * FS / 2), n=len(x))
    return signal.sosfiltfilt(signal.butter(4, [35, 350], fs=FS, btype='bandpass', output='sos'), w)

def main():
    DATA.mkdir(exist_ok=True)
    endpoint = 'https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/'
    meta_path = DATA/'event_metadata.json'
    if not meta_path.exists():
        r=requests.get(endpoint,timeout=40);r.raise_for_status();meta_path.write_text(r.text,encoding='utf-8')
    event=next(iter(json.loads(meta_path.read_text(encoding='utf-8'))['events'].values()))
    result={'event':'GW150914','gps':GPS,'sample_rate':FS,'processing':'32 s strain; off-source Welch PSD (first/last 10 s, 2 s segments); whitening; zero-phase 35–350 Hz Butterworth order 4','detectors':{}}
    manifest=[]
    for detector in ['H1','L1']:
        item=next(x for x in event['strain'] if x['detector']==detector and x['duration']==32 and x['sampling_rate']==FS and x['format']=='txt')
        path=DATA/(detector+'.txt.gz')
        if not path.exists():
            r=requests.get(item['url'],timeout=60);r.raise_for_status();path.write_bytes(r.content)
        raw=path.read_bytes()
        x=np.loadtxt(io.BytesIO(gzip.decompress(raw)))
        if x.ndim!=1 or len(x)!=FS*32 or not np.all(np.isfinite(x)):raise ValueError('Unexpected strain shape or missing samples')
        f,p1=signal.welch(x[:FS*10],fs=FS,nperseg=FS*2,average='median')
        _,p2=signal.welch(x[-FS*10:],fs=FS,nperseg=FS*2,average='median')
        psd=np.maximum(np.interp(np.fft.rfftfreq(len(x),1/FS),f,(p1+p2)/2),1e-49)
        y=whiten_filter(x,psd)
        t=np.arange(len(x))/FS+item['GPSstart']-GPS
        np.savez_compressed(DATA/(detector+'.npz'),time=t,raw=x,white=y,psd=psd)
        take=(t>=-.5)&(t<=.25)
        ff,tt,power=signal.spectrogram(y,fs=FS,nperseg=128,noverlap=112)
        tt+=t[0];sel_t=(tt>=-.5)&(tt<=.25);sel_f=(ff>=32)&(ff<=512)
        result['detectors'][detector]={'t':np.round(t[take],7).tolist(),'strain':np.round(y[take],6).tolist(),'spectrogram':{'t':np.round(tt[sel_t],6).tolist(),'f':ff[sel_f].tolist(),'power':np.round(10*np.log10(np.maximum(power[np.ix_(sel_f,sel_t)],1e-12)),3).tolist()}}
        manifest.append({'detector':detector,'url':item['url'],'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'gps_start':item['GPSstart'],'sample_rate':FS,'duration':32})
        print(detector,len(x),'samples',flush=True)
    (DATA/'event.json').write_text(json.dumps(result,separators=(',',':')),encoding='utf-8')
    (DATA/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')

if __name__=='__main__':main()
