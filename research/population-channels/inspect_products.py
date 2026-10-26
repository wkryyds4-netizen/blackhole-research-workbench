"""Inspect verified public channel products without unpickling model objects."""
from pathlib import Path,PurePosixPath
import hashlib
import io
import json
import pickletools
import tarfile
import zipfile
import h5py
import numpy as np

BASE=Path(__file__).resolve().parent
FOLDER=BASE/'sources'


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    inventory=json.loads((BASE/'archive-inventory-14967688.json').read_text())
    extracted=FOLDER/'inspected-products';extracted.mkdir(exist_ok=True)
    report=dict(schema='blackhole-physical-channel-product-inspection/1',archives=[],json_metadata={},numeric_mappings={},hdf={},torch_archives={})
    for product in inventory:
        record=product['source'];path=FOLDER/record['file'];assert sha(path)==record['sha256']
        report['archives'].append(record)
        with tarfile.open(path,'r:gz') as archive:
            for member in archive:
                p=PurePosixPath(member.name)
                assert not p.is_absolute() and '..' not in p.parts and not member.issym() and not member.islnk()
                if not member.isfile() or not p.name.endswith(('.hdf5','.json','.npy','.pt')):continue
                assert member.size<80000000
                data=archive.extractfile(member).read()
                assert len(data)==member.size
                target=extracted/p.name;target.write_bytes(data)
                if p.suffix=='.json':report['json_metadata'][p.name]=json.loads(data)
                elif p.suffix=='.npy':
                    array=np.load(io.BytesIO(data),allow_pickle=False)
                    assert array.dtype.kind in 'iuf'
                    report['numeric_mappings'][p.name]=array.tolist()
                elif p.suffix=='.hdf5':
                    fields=[]
                    with h5py.File(target,'r') as h:
                        def inspect(name,obj):
                            if isinstance(obj,h5py.Dataset):
                                fields.append(dict(path=name,shape=list(obj.shape),dtype=str(obj.dtype)))
                        h.visititems(inspect)
                    report['hdf'][p.name]=dict(sha256=sha(target),fields=fields)
                elif p.suffix=='.pt':
                    with zipfile.ZipFile(io.BytesIO(data)) as zipped:
                        pickles=[n for n in zipped.namelist() if n.endswith('/data.pkl')]
                        assert len(pickles)==1
                        raw=zipped.read(pickles[0]);assert len(raw)<1000000
                        # Bytecode inspection only. pickletools does not resolve
                        # or execute referenced globals or reconstruct objects.
                        globals_seen=[str(arg) for op,arg,_ in pickletools.genops(raw) if op.name=='GLOBAL']
                        strings=[str(arg) for op,arg,_ in pickletools.genops(raw) if op.name in ['BINUNICODE','SHORT_BINUNICODE']]
                        report['torch_archives'][p.name]=dict(sha256=sha(target),globals=sorted(set(globals_seen)),
                                                             state_key_examples=[s for s in strings if 'weight' in s or 'mask' in s][:18],zip_entries=len(zipped.namelist()),executed=False)
    report['inspection_tool_sha256']=sha(Path(__file__))
    report['limits']='Numeric arrays and metadata inspected; model pickle bytecode examined without execution. No channel fractions or validated model evaluations yet.'
    (BASE/'product-inspection.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['hdf','archives']},indent=2))
    for name,entry in report['hdf'].items():
        print(name,len(entry['fields']),'datasets')
        print(json.dumps(entry['fields'][:60],indent=2))


if __name__=='__main__':main()
