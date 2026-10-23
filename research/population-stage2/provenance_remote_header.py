"""Read HDF5 attributes via verified HTTP byte ranges, bounded to 2 MiB/file.

Never fetch full large datasets. Requests returning 200 instead of 206 abort
without consuming the body. The range audit records bytes actually transferred.
"""
from pathlib import Path
import io, json, hashlib
import requests, h5py, numpy as np
OUT=Path(__file__).resolve().parent/'provenance-sources'
META=json.loads((OUT/'zenodo-5636816.json').read_text())

class RangeReader(io.RawIOBase):
    def __init__(self,url,size):
        self.url=url; self.size=size; self.pos=0; self.cache={}; self.log=[]
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self,offset,whence=0):
        self.pos=offset if whence==0 else (self.pos+offset if whence==1 else self.size+offset)
        return self.pos
    def readinto(self,b):
        result=self.read(len(b)); b[:len(result)]=result; return len(result)
    def read(self,n=-1):
        if n<0 or n>1024*1024: raise ValueError('Unbounded read refused')
        end=min(self.size,self.pos+n); result=bytearray()
        while self.pos<end:
            start=self.pos//65536*65536; stop=min(self.size,start+65536)-1
            if start not in self.cache:
                if sum(len(v) for v in self.cache.values())+stop-start+1>2*1024*1024:
                    raise ValueError('2 MiB per-file metadata budget exhausted')
                with requests.get(self.url,headers={'Range':f'bytes={start}-{stop}'},stream=True,timeout=40) as r:
                    if r.status_code!=206:raise ValueError(f'Range unsupported HTTP {r.status_code}; body not consumed')
                    assert r.headers.get('Content-Range')==f'bytes {start}-{stop}/{self.size}',r.headers
                    data=r.content
                    assert len(data)==stop-start+1
                self.cache[start]=data
                self.log.append(dict(start=start,end=stop,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
            amount=min(end-self.pos,stop+1-self.pos)
            result.extend(self.cache[start][self.pos-start:self.pos-start+amount]); self.pos+=amount
        return bytes(result)

def serial(v):
    if isinstance(v,bytes):return v.decode(errors='replace')
    if isinstance(v,np.ndarray):return [serial(x) for x in v.tolist()]
    if isinstance(v,np.generic):return serial(v.item())
    return v

results=[]
for file in META['files']:
    if not file['key'].endswith('.hdf5'):continue
    reader=RangeReader(file['links']['self'],file['size'])
    rec=dict(file=file['key'],url=reader.url,size=reader.size,official_checksum=file['checksum'])
    try:
        with h5py.File(reader,'r') as h:
            rec['root_attrs']={k:serial(v) for k,v in h.attrs.items()}
            group=h['injections'];rec['injection_attrs']={k:serial(v) for k,v in group.attrs.items()}
            rec['dataset_metadata']={k:dict(shape=group[k].shape,dtype=str(group[k].dtype)) for k in group}
    except Exception as e:
        rec['error']=f'{type(e).__name__}: {e}'
    rec['ranges']=reader.log; rec['bytes_transferred']=sum(r['bytes'] for r in reader.log)
    results.append(rec)
    print(rec['file'],rec.get('root_attrs',{}),rec.get('error',''),rec['bytes_transferred'],flush=True)
    (OUT/'remote-headers.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
