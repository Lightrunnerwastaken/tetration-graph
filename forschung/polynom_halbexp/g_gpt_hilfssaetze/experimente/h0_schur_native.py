"""Public FLINT C API bridge: opaque storage, no Python-object layout access.

Numerical midpoint arithmetic only. Matrix windows do not own their data.
The DLL comes from the installed python-flint wheel. See FLINT arb_mat.h.
"""
import ctypes as C
from pathlib import Path
import flint

P=C.c_void_p
S=C.c_longlong # FLINT's 64-bit slong, including its Windows build
U=C.c_ulonglong

class Native:
    def __init__(self):
        libdir=Path(flint.__file__).parent.parent/'python_flint.libs'
        self.lib=C.CDLL(str(next(libdir.glob('libflint*.dll'))))
        signatures={
            'arb_mat_init':[P,S,S], 'arb_mat_clear':[P],
            'arb_mat_window_init':[P,P,S,S,S,S], 'arb_mat_entry_ptr':[P,S,S],
            'arb_mat_set':[P,P], 'arb_mat_get_mid':[P,P],
            'arb_mat_transpose':[P,P], 'arb_mat_neg':[P,P],
            'arb_mat_approx_mul':[P,P,P,S], 'arb_mat_add':[P,P,P,S],
            'arb_mat_scalar_div_arb':[P,P,P,S],
            'arb_mat_scalar_mul_arb':[P,P,P,S],
            'arb_mat_scalar_addmul_arb':[P,P,P,S],
            'arb_mat_frobenius_norm':[P,P,S],
            'arb_set_str':[P,C.c_char_p,S], 'arb_get_str':[P,S,U],
            'arb_init':[P], 'arb_clear':[P], 'arb_set':[P,P],
            'arb_get_mid_arb':[P,P], 'arb_one':[P], 'arb_sub':[P,P,P,S],
            'arb_div':[P,P,P,S], 'arb_contains_zero':[P], 'flint_free':[P],
        }
        self.functions={}
        for name,args in signatures.items():
            fn=getattr(self.lib,name);fn.argtypes=args
            fn.restype=P if name in ('arb_mat_entry_ptr','arb_get_str') else (C.c_int if name in ('arb_set_str','arb_contains_zero') else None)
            self.functions[name]=fn
    def __getattr__(self,name):return self.functions[name]
    def matrix(self,r,c):
        buf=C.create_string_buffer(128);self.arb_mat_init(buf,r,c);return buf
    def scalar(self):
        buf=C.create_string_buffer(128);self.arb_init(buf);return buf
    def window(self,mat,r1,c1,r2,c2):
        buf=C.create_string_buffer(128);self.arb_mat_window_init(buf,mat,r1,c1,r2,c2);return buf
    def entry(self,mat,i,j):return self.arb_mat_entry_ptr(mat,i,j)
    def string(self,x,digits=60):
        p=self.arb_get_str(x,digits,2)
        try:return C.string_at(p).decode('ascii')
        finally:self.flint_free(p)
    def set_string(self,x,s,bits):
        if self.arb_set_str(x,str(s).encode('ascii'),bits):raise ValueError(s)
        self.arb_get_mid_arb(x,x)

class Border:
    def __init__(self,J,bits):
        self.api=F=Native();self.N=J.nrows()-1
        self.J=F.matrix(self.N+1,self.N+1)
        self.inv=F.matrix(self.N+1,self.N+1)
        self.tmp=F.matrix(self.N+1,self.N+1)
        self.z=F.matrix(self.N+1,1);self.w=F.matrix(1,self.N+1)
        self.Q=F.matrix(self.N+1,1);self.dot=F.matrix(1,1)
        self.h=F.scalar();self.one=F.scalar();self.qnorm=F.scalar();self.pivot=F.scalar()
        F.arb_one(self.one)
        digits=int(bits*.30103)+12
        for i in range(self.N+1):
            for j in range(self.N+1):
                F.set_string(F.entry(self.J,i,j),J[i,j].mid().str(digits,radius=False),bits)
        F.arb_div(F.entry(self.inv,0,0),self.one,F.entry(self.J,0,0),bits)
        F.arb_mat_get_mid(self.inv,self.inv)
    def step(self,n,bits):
        F=self.api
        if n<1 or n>self.N:raise ValueError(n)
        I=F.window(self.inv,0,0,n,n)
        b=F.window(self.J,0,n,n,n+1);r=F.window(self.J,n,0,n+1,n)
        z=F.window(self.z,0,0,n,1);w=F.window(self.w,0,0,1,n)
        F.arb_mat_approx_mul(z,I,b,bits)
        F.arb_mat_approx_mul(w,r,I,bits)
        F.arb_mat_approx_mul(self.dot,r,z,bits)
        F.arb_sub(self.h,F.entry(self.J,n,n),F.entry(self.dot,0,0),bits)
        F.arb_get_mid_arb(self.h,self.h)
        if F.arb_contains_zero(self.h):raise ArithmeticError(f'zero Schur pivot n={n}')
        F.arb_div(self.pivot,self.one,self.h,bits)
        F.arb_get_mid_arb(self.pivot,self.pivot)
        Q=F.window(self.Q,0,0,n+1,1)
        F.arb_mat_neg(F.window(self.Q,0,0,n,1),z)
        F.arb_one(F.entry(self.Q,n,0))
        F.arb_mat_frobenius_norm(self.qnorm,Q,bits)
        F.arb_get_mid_arb(self.qnorm,self.qnorm)
        col=F.window(self.inv,0,n,n,n+1);row=F.window(self.inv,n,0,n+1,n)
        F.arb_mat_scalar_mul_arb(col,z,self.pivot,bits)
        F.arb_mat_get_mid(col,col)
        T=F.window(self.tmp,0,0,n,n)
        F.arb_mat_approx_mul(T,col,w,bits)
        F.arb_mat_add(I,I,T,bits)
        F.arb_mat_neg(col,col)
        F.arb_mat_scalar_mul_arb(row,w,self.pivot,bits);F.arb_mat_neg(row,row)
        F.arb_set(F.entry(self.inv,n,n),self.pivot)
        F.arb_mat_get_mid(F.window(self.inv,0,0,n+1,n+1),F.window(self.inv,0,0,n+1,n+1))
        return self
    def close(self):
        F=self.api
        for name in ('J','inv','tmp','z','w','Q','dot'):F.arb_mat_clear(getattr(self,name))
        for name in ('h','one','qnorm','pivot'):F.arb_clear(getattr(self,name))
