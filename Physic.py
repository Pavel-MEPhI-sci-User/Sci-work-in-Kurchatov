# -*- coding: utf-8 -*-
"""
Here I realise all phisic thermodinamic quantities
for my quest


Physic
    |
    |---С_e
    |---C_i
    |---K_e
    |---K_i
    |
    |

NM --перевод всего в СГС! NB
"""

class Physic:    
    """
    В этом классе я введу все физические свойства интересующего материала
    и свойства электронов
    
    пока тут Fe85B15
    """
    
    #def __init__(self, electric_prop = 'Fe85B15', lattice = 'Fe85B15'):
        #if electric_prop == 'Fe85B15':
    def __init__(self):
        self.dzeta = 8.97e-5 #J/cm**3 *K
        self.T_melt = 1538
                  
        
    def C_e(self, T_e):
        # J g-1 K-1
        #самый простой вариант возьму мажоранту
        #но что то тут не то...
        # k = 1.38e-16
        #Ne = 8.48e+22 * 26
        return 1.5 * (1.38e-16) * (26*8.48e22)
    
    def C_i(self, T_i):
        # J g-1 K-1
        if 10 < T_i < 100 :
            return  0.17 - 0.002*T_i +7.3e-5*T_i**2 - 3.3e-7*T_i**3
        elif 100 < T_i < 300 :
            return -0.20 + 6e-3*T_i -2e-5*T_i**2 +2.5e-8*T_i**3
        elif 300 < T_i < self.T_melt :
            return 0.79 + 5.4e-6*T_i
        elif T_i > self.T_melt :
            return 0.8
    
    def K_e(self, T_e):
        # W cm-1 K-1
        pass
    
    def K_i(self, T_i):
        # W cm-1 K-1
        T = T_i
        if 1 < T < 20 :
            return -0.45 + 0.97*T - 22e-3*T**2
        elif 20 < T < 100 :
            return 16.5 - 0.35*T + 2e-3*T**2
        elif 100 < T < self.T_melt :
            return 1.24 - 17e-4*T + 8.8e-7*T**2 - 1.3e-10*T**3
        elif T > self.T_melt :
            return 0.33
        
        


































