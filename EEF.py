import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spint
from scipy import integrate
import math
import const1

F1 = lambda x: x/np.log(x)

F2 = lambda x: math.sqrt(x)/np.log(x)


def E(r,t):
    if r == 0:
        return 0.
    else:
        N=274
        y=np.linspace(1.,const1.E_m,N+1)
        sum=0.
        h=(const1.E_m-1.)/N
        for i in range(N):
            up1=r1(y[i],t)
            up2=r1(y[i+1],t)
            secont_int1=spint.quad(fun,0.,up1,args=(r,y[i]))[0]
            secont_int2=spint.quad(fun,0.,up2,args=(r,y[i+1]))[0]
            sum+=((secont_int1/y[i]**2)+(secont_int2/y[i+1]**2))*0.5*h
        return const1.E_0*r*sum


def fun(x,r,y):
    r_0=math.sqrt((1./y)-(1./const1.E_m))
    a=x-(r-r_0)**2
    b=(r+r_0)**2-x
    c=1.-0.5*(r**2+r_0**2-x)/r**2
    if (a <= 0.) or (b <= 0.):
        return 0.
    else:
        return c/(x*math.sqrt(a*b))


def r1(y2,t):
    q=(const1.I**2/(const1.po_max*math.pi*const1.Z_2*const1.N_e*const1.charge**4))*math.sqrt(1.-y2/const1.E_m)
    if y2 == 1.:
        return 0.
    else:
        v=E_1(t,y2)
        return (q*(spint.quad(F1,v,y2)[0]))**2

def E_1(t,y3):
    d=const1.I**2/(math.pi*const1.Z_2*const1.N_e*const1.charge**4*const1.V_I*t)
    a=1.+1e-11
    b=y3
    z=0.5*(b+a)
    eps=1e-10
    while(math.fabs(b-a) > eps):
        w_m=1.-d*(spint.quad(F2,z,y3)[0])
        w_l=1.-d*(spint.quad(F2,a,y3)[0])
        w_r=1.-d*(spint.quad(F2,b,y3)[0])
        if w_m*w_l < 0.:
            b=z
        elif w_m*w_r < 0.:
            a=z
        #else:
        #   break
        
        z=0.5*(b+a)
    
    return z

T=1e-16
R=5.
M=100
r=np.linspace(0.,R,M)
E_E=np.zeros(M)

for i in range(M):
    E_E[i]=E(r[i],T)
    print(E_E[i], '   ', i)
E_E=E_E*(3e-6)
r=r*4.5
plt.plot(r,E_E)
plt.show()

