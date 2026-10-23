"""Explicit build-only portability patch for the verified 1.2.0 source archive."""
from pathlib import Path
import difflib, json, hashlib
root=Path(__file__).resolve().parent
base=root/'adapter-build/gwsurrogate-1.2.0'
patches=[]
for relative in ['setup.py','gwsurrogate/spline_interp_Cwrapper/_spline_interp.cpp']:
    p=base/relative; old=p.read_text(encoding='utf-8'); new=old
    if relative=='setup.py':
        new=new.replace('extra_compile_args=["-std=c++14", "-O3"] + _arch_flags,','extra_compile_args=["/std:c++14", "/O2", "/utf-8"],\n    export_symbols=["spline_interp", "spline_interp_multi", "spline_interp_multi_complex"],')
        new=new.replace('extra_compile_args=["-std=c99", "-fPIC", "-O3", \'-Wcpp\'] + _arch_flags,','extra_compile_args=["/std:c11", "/O2", "/utf-8"],')
    else:
        new=new.replace('#define restrict __restrict__','#ifdef _MSC_VER\n#define restrict __restrict\n#else\n#define restrict __restrict__\n#endif')
    assert new != old, 'Apply only once to clean release source'
    patches.extend(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+relative,tofile='b/'+relative))
    p.write_text(new,encoding='utf-8')
(root/'adapter-sources/gwsurrogate-windows.patch').write_text(''.join(patches),encoding='utf-8')
print('MSVC flags, restrict spelling, ctypes symbol exports; numerical formulas unchanged')
