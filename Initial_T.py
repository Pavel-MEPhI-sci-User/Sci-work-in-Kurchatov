# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 20:07:59 2017

@author: Alpharius

here I made class which make initial condition for my SOLUTOR class

        Data_creator_for_FILD
        /
       /
Initial_T
|
|---make initial
|
|
|
"""

from Data_creator_for_FILD import *
import pandas as pd

def Initial_T():
    
    @staticmethod
    def T_init():
        pass
    
    def __init__(self):
        DcfF = Data_creator_for_FILD()
        self.mas_field = DcfF[ 9.52631578947e-17]   #подбираем такое время которое совпадает со временем релаксации
        print(self.mas_field)
        

a = Initial_T()