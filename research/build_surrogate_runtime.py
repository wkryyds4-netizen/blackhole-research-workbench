"""Compile the release with explicit LLP64 integer ABI fix; no formula edits."""
from pathlib import Path
import sys,sysconfig,subprocess,tarfile,hashlib,json,os,difflib
import numpy as np
root=Path(__file__).resolve().parent/'adapter-build'
runtime=root/'zig-runtime'
archive=root.parent/'adapter-sources/gwsurrogate-1.2.0.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest()=='316a640576d32aaa45aaf9a6d876b240af356c89131e6af8257039b3818baae8'
runtime.mkdir(exist_ok=True)
with tarfile.open(archive) as tf:
    tf.extractall(runtime,filter='data')
src=runtime/'gwsurrogate-1.2.0'
# HDF5 basis orders are int64. Windows LLP64 long is 32-bit; Linux LP64
# long is 64-bit. Preserve the stored integer width explicitly, without
# changing the fitted polynomial or the floating-point arithmetic.
cpath=src/'gwsurrogate/precessing_utils/src/precessing_utils.c'
old=cpath.read_text(); new=old.replace('long *bf_order_data','npy_int64 *bf_order_data').replace('long   *bf_order_data','npy_int64 *bf_order_data').replace('const long *orders','const npy_int64 *orders').replace('(long *)   PyArray_DATA(bf_orders)','(npy_int64 *) PyArray_DATA(bf_orders)').replace('(long *) PyArray_DATA(bf_orders)','(npy_int64 *) PyArray_DATA(bf_orders)')
cpath.write_text(new)
(root/'gwsurrogate-zig-int64.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/precessing_utils.c',tofile='b/precessing_utils.c')))
zig=root/'zig-x86_64-windows-0.14.1/zig.exe'
base=Path(sys.base_prefix)
common=['-target','x86_64-windows-gnu','-shared','-O2','-I'+str(base/'include'),'-I'+np.get_include()]
suffix=sysconfig.get_config_var('EXT_SUFFIX')
(src/'spline_exports.def').write_text('EXPORTS\nspline_interp\nspline_interp_multi\nspline_interp_multi_complex\n')
targets=[('cc','gwsurrogate/precessing_utils/src/precessing_utils.c','gwsurrogate/precessing_utils/_utils', ['-std=c11','-DCMPLX(x,y)=__builtin_complex((double)(x),(double)(y))','-Igwsurrogate/precessing_utils/include']),('c++','gwsurrogate/spline_interp_Cwrapper/_spline_interp.cpp','gwsurrogate/spline_interp_Cwrapper/_spline_interp',['-std=c++14','-Wl,--export-all-symbols'])]
records=[]
env=dict(os.environ,ZIG_GLOBAL_CACHE_DIR=str(root/'zig-cache'),ZIG_LOCAL_CACHE_DIR=str(root/'zig-local-cache'))
for lang,source,out,flags in targets:
    cmd=[str(zig),lang,*common,*flags,source,*(['spline_exports.def'] if lang=='c++' else []),str(base/'libs/python312.lib'),'-o',out+suffix]
    result=subprocess.run(cmd,cwd=src,env=env,capture_output=True,text=True)
    (root/(lang.replace('+','p')+'-build.log')).write_text(result.stdout+result.stderr,encoding='utf-8')
    print(result.stdout+result.stderr,flush=True)
    result.check_returncode()
    p=src/(out+suffix)
    records.append({'source':source,'source_sha256':hashlib.sha256((src/source).read_bytes()).hexdigest(),'binary_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'command':cmd})
report={'compiler':'zig 0.14.1 x86_64-windows-gnu','source_archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'formula_edits':False,'portability_edits':['bfOrders long pointers replaced with npy_int64 for Windows LLP64','CMPLX macro supplied via compiler builtin complex constructor','three ctypes spline functions explicitly exported via .def file'],'build':records}
(root/'zig-runtime-build.json').write_text(json.dumps(report,indent=2))
results=root.parent/'surrogate-results';results.mkdir(exist_ok=True)
(results/'runtime-build.json').write_text(json.dumps(report,indent=2))
(results/'windows-int64.patch').write_text((root/'gwsurrogate-zig-int64.patch').read_text())
print('Built ABI-adapted runtime',src)
