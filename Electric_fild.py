# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 21:13:48 2017

@author: Alpharius

        Ion_Class
      /
     /
    /
Electric_Field
|
|---root_eps_til_1
|---r1_til
|---f_right
|---E

Data_creator_for_FILD
|
|---create_df
|---read_csv
|---create_plot


"""
from Ionclass import *
import pandas as pd
from matplotlib import animation
import numpy as np
from matplotlib import pyplot as plt


class Electric_Field(Ion_Class):
    
    def __init__(self, t, q_til, E, phy_prop_i, phy_prop_target ):
        super().__init__( E, phy_prop_i, phy_prop_target)
        self.t = t
        self.q_til = q_til
        
        #СГС напряжение  http://www.akin.ru/spravka/s_rel2.htm
        self.Eo = 2*self.e*self.N*self.Z2*self.q_max()  * 3e-6
        
        #my functions for numerical calculation
        self.F1 = lambda x: x/np.log(x)
        self.F2 = lambda x: sqrt(x)/np.log(x)
    
    def root_eps_til_1(self, eps_til_0):
        v_i = sqrt(2*self.I*self.ev_to_erg/self.m_e)
        term1 = (self.I * self.ev_to_erg)**2
        term2 = pi*self.Z2*self.e**4*self.t*v_i*self.N
               
        a = 1. + 1e-11
        b = eps_til_0
        m = (a+b)/2
        EPS = 1e-9
        while (abs(b - a)>EPS):
            v_a =1 - term1/term2 * scipy.integrate.quad(self.F2,a, eps_til_0)[0]
            v_m =1 - term1/term2 * scipy.integrate.quad(self.F2,m, eps_til_0)[0]
            v_b =1 - term1/term2 * scipy.integrate.quad(self.F2,b, eps_til_0)[0]
            #print('a = ',a, v_a)
            #print('m = ',m, v_m)
            #print('b = ',b, v_b)
            if v_m * v_a < 0:
                b = m
            else:
                a = m
                
            m = (a + b)/2
        #print(eps_til_0 - m)
        return m
    
    def r1_til(self, eps_til_0):
        term1 = (self.I * self.ev_to_erg)**2 / self.q_max()#сразу обезразмериваю
        term2 = pi*self.Z2*self.N*self.e**4
        term3 = sqrt(1 - eps_til_0* self.I/self.eps_m())
        c = term1*term3/term2        
        return c * scipy.integrate.quad(self.F1,
                                    self.root_eps_til_1(eps_til_0),
                                    eps_til_0)[0]
    
    def f_right(self, x, eps_til_0):
        r_til_0 = sqrt(1./eps_til_0 - self.I/self.eps_m())
        term1 = 1 - (self.q_til**2 + r_til_0**2 - x)/(2*self.q_til**2)
        term2 = x - (self.q_til - r_til_0)**2
        term3 = (self.q_til + r_til_0)**2 - x
        if (term2 <= 0 ) or (term3 <=0) or (x == 0) :
            return 0
        else:
            return term1/(x * sqrt(term2) * sqrt(term3))
        
    def E(self):
        a = time.time()
        if self.q_til == 0:
            return 0
        else:
            N = 200
            mas_eps_0 = np.linspace(11, self.eps_m()/self.I, N+1)
            h = (self.eps_m()/self.I - 1.1)/N
            integral = 0.
            for i in range(N):
                up_1 = self.r1_til(mas_eps_0[i])**2
                up_2 = self.r1_til(mas_eps_0[i+1])**2
                int1 = scipy.integrate.quad(self.f_right,
                                           0., up_1,
                                           args =(mas_eps_0[i],) )[0]
                int2 = scipy.integrate.quad(self.f_right,
                                           0., up_2,
                                           args =(mas_eps_0[i+1],) )[0]
                integral += (int1/(mas_eps_0[i]**2) +\
                             int2/(mas_eps_0[i+1]**2))*0.5*h
            b = time.time()
            print('time = ',b-a)
            return self.Eo*self.q_til*integral






