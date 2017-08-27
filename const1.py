import math

Z_2=26
Z_1=92
#N_e=Z_2*8.47e+22
N_e=8.47e+22 
m_e=9.1e-28                     #грамм
M=238*1.66e-24                  #грамм
I=8*1.6e-12                     #эрг
po_max=4.5e-8                   #см
charge=4.8e-10                  #ед.заряда СГС
E_1=(238*10e+6)*1.6e-12         #эрг
V_I=math.sqrt(2*I/m_e)          #см/с
E_m=(4*m_e*E_1/M)/I             #эрг
E_0=2*charge*N_e*Z_2*po_max     #эрг
M_2=56*1.66e-24
k_boltz=1.38e-16

#print(E_0)
#print(E_m)
#print(math.sqrt((1./1.)-(1./E_m)))
#print(math.sqrt((1./E_m)-(1./E_m)))



