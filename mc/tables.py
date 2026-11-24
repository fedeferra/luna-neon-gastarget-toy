import numpy as np, xraylib
NA=6.02214e23
a=[8.7842e+2,-1.9625e+3,1.2949e+3,-2.0028e+2,1.2575e+1,-2.8333e-1]
b=[-1.0342e+1,1.7692e+1,-8.2381,1.3063,-9.0815e-2,2.3586e-3]
c=[-4.5263e+2,1.1161e+3,-8.6749e+2,2.1773e+2,-2.0467e+1,6.5372e-1]
me=0.510999
def pair_ub_raw(Z,E):
    X=np.log(E/me);P=lambda k:sum(k[i]*X**i for i in range(6));return (Z+1)*Z*(P(a)+P(b)*Z+P(c)/Z)
def pair_b(Z,E):   # barn per atom, E MeV
    if E<=2*me: return 0.0
    if E<1.5: return pair_ub_raw(Z,1.5)*1e-6*((E-2*me)/(1.5-2*me))**2
    return pair_ub_raw(Z,E)*1e-6
def kn_b(E):
    k=E/me;r0=2.8179403e-13
    s=2*np.pi*r0**2*((1+k)/k**2*(2*(1+k)/(1+2*k)-np.log(1+2*k)/k)+np.log(1+2*k)/(2*k)-(1+3*k)/(1+2*k)**2);return s*1e24
Z={'Na':11,'I':53,'Bi':83,'Ge':32,'O':8,'Fe':26,'Cr':24,'Ni':28}
A={'Na':22.98977,'I':126.90447,'Bi':208.9804,'Ge':72.63,'O':15.999,'Fe':55.845,'Cr':51.9961,'Ni':58.6934}
EG=np.logspace(np.log10(0.010),np.log10(20),1700)   # MeV
def el(sym):
    z=Z[sym];ph=np.zeros_like(EG);co=np.zeros_like(EG);pr=np.zeros_like(EG)
    Elim=0.78   # MeV, safe limit of xraylib
    # photoelectric: xraylib up to Elim, then power law with local slope (1..3), slope 1 above 2 MeV
    p1=xraylib.CS_Photo(z,Elim*1000);p0=xraylib.CS_Photo(z,500.0);n=np.clip(np.log(p0/p1)/np.log(Elim/0.5),1,3)
    for i,E in enumerate(EG):
        if E<=Elim:
            ph[i]=xraylib.CS_Photo(z,E*1000);co[i]=xraylib.CS_Compt(z,E*1000)
        else:
            ph[i]=p1*(E/Elim)**(-n) if E<2 else p1*(2/Elim)**(-n)*(E/2)**(-1)
            co[i]=z*kn_b(E)*1e-24*NA/A[sym]
        pr[i]=pair_b(z,E)*1e-24*NA/A[sym]
    # Compton matching at Elim
    kn_l=z*kn_b(Elim)*1e-24*NA/A[sym];xc=xraylib.CS_Compt(z,Elim*1000)
    return ph,co,pr,xc/kn_l
mats={'nai':{'rho':3.67,'w':{'Na':22.98977/(22.98977+126.90447),'I':126.90447/(22.98977+126.90447)},'Z':50},
      'bgo':{'rho':7.13,'w':None,'Z':70},'ge':{'rho':5.323,'w':{'Ge':1.0},'Z':32},'steel':{'rho':8.0,'w':{'Fe':0.715,'Cr':0.19,'Ni':0.095},'Z':26}}
mb=4*A['Bi']+3*A['Ge']+12*A['O'];mats['bgo']['w']={'Bi':4*A['Bi']/mb,'Ge':3*A['Ge']/mb,'O':12*A['O']/mb}
E=[];out=[]
elc={s:el(s) for s in Z}
for s in Z: print(s,'Compton xraylib/KN at 0.78MeV',round(elc[s][3],4))
with open('tables.bin','wb') as f:
    EG.astype('<f8').tofile(f)
    for m in ['nai','bgo','ge','steel']:
        d=mats[m];ph=sum(d['w'][s]*elc[s][0] for s in d['w'])*d['rho'];co=sum(d['w'][s]*elc[s][1] for s in d['w'])*d['rho'];pr=sum(d['w'][s]*elc[s][2] for s in d['w'])*d['rho']
        for arr in (ph,co,pr): arr.astype('<f8').tofile(f)
        np.array([d['Z']],dtype='<f8').tofile(f)
        i=np.searchsorted(EG,1.0);print(m,'mu(1MeV)/rho',(ph[i]+co[i]+pr[i])/d['rho'],'mu(10MeV)/rho',(ph[-1]+co[-1]+pr[-1])/d['rho'],'mu(0.1MeV)/rho',(ph[np.searchsorted(EG,0.1)]+co[np.searchsorted(EG,0.1)])/d['rho'])
