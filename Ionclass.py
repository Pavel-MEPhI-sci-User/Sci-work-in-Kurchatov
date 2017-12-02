# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.

Ion_Class
    |
    |---eps_m
    |---Z1_star
    |---q_max
    |---dEdz
    |---interpolated_dEdz
    |---coordinat_as_function_from_energi
    |---interpolated_coordinat_as_function_from_energi
    |---velocity_how_function_from_E
    |---coordinate_from_time
"""

from IPython.display import Image
from math import *
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.integrate as integrate
import time
from scipy.interpolate import interp1d 
#%matplotlib inline


from scipy.special import expi as Ei
from scipy.optimize import newton, brenth


def vect(f):
    return np.vectorize(f)

U_ion = (98, 238)
Fe = (26, 56, 7.8)
Fe_ion = (26, 56)
E_init_Fe = Fe_ion[1]*10e6
E_init_U = 2380e6


class Ion_Class:
    """
    класс характеризующий состояние иона в нутри металла-вводим сюда характеристики иона 
    в данный момент на выходе получаем все необходимые его свойства в этот 
    же момент, +некоторое описание среды которое влияет на тормозные потери
    иона.
    
    В итоге мы можем посчитать энергопотерии
    
    такие как:
        eps_m - максимально возможная энергия дельта электрона
        Z_star1 - эффективный заряд движущегося иона
        q_max - область ионизации    
    
    """
    
    def __init__(self,E, phy_prop_i, phy_prop_target ):
        """
        Вводные параметры:
        E - Энергия движущегося иона [эв]
        phy_prop_i - кортедж который хранит в себе физические свойства тяжелого иона
                   (Z1,A) -заряд иона , число его нуклонов, плотность материала
                   #self.Z1  --порядковый номер иона  #безразмерно   
                   #self.A   --число нуклонов иона    #безразмерно
                   
                   
        phy_prop_target - кортедж который хранит в себе физические свойства среды
                   #self.Z2 --порядковый номер среды
                   #self.A_o   --число нуклонов вмешенях    #безразмерно    
                   #self.rho --плотность материала иона  # грамм на см**3
                   
        U = (98, 238)
        Fe = (26, 56, 7.8)
        Fe_ion = (26, 56)
        
        """     
        
        self.Z1, self.A = phy_prop_i                                      
        self.Z2, self.A_o, self.rho = phy_prop_target
        
        self.m_e = 9.1e-28                      #гр
        self.m_ion = self.A * 1.6e-24           #грамм
        self.E1 = E                             #в эв
        self.N = self.rho/(1.66e-24*self.A_o)   #1|cm3 
        self.e = 4.8e-10                        #эд СГСЭ  [гр1/2 * см 3/2 с-1]
        self.I = 8                              #эв энергия иониз Fe
        
        self.ev_to_erg = 1.6e-12                #коэффициент перевода
    
    def eps_m(self):
        """ 
        максимально возможная энергия дельта электрона в [эв]----тут я уверен
        возвращает:
                    [эв]
        """
        #возвращает эв
        return 4  * self.m_e *self.E1/self.m_ion
    
    def Z1_star(self):
        """
        эффективный заряд движущегося иона [безразмерная величина]         
        примечание:
            скорость есть по сути энергия но здесь я сознательно не усложняю функцию
            Работает так как надо нареканий нет
        """
        #v в метрах на секунду
        v = sqrt(2000*self.E1*1.6e-19/(self.m_ion))    #метры в секунду
        с = 3e8     #метров в секунду
        return self.Z1*(1 - exp((-137*v)/(с * self.Z1**(2./3))))
    
    def q_max(self):
        """
        функция возвращающая размер области ионизации   ----тут я уверен
        принимает:
            [E] энергию иона в [эв] (внутри я перевожу в эрг)
        возвращает:
            длину в [см]
        замечание:
            произвел проверку на данных из статьи --полностью удовлетворен работой модуля        
        """       
        term1 = 2*self.Z1_star()*(self.e**2)/(self.eps_m() * self.ev_to_erg)
        term2 = sqrt(self.eps_m()/self.I)
        return term1 * term2
    
    
    
    def dEdz(self):
        """
        Возвращает кэв на ангстрем . но это не окончяательно. в общем это надо еще проверить --еще не проверял!
        по табличным результатам все подходит
        """
        term1 = 4*pi*self.e**4 * self.N * 16 * self.Z1_star()**2  #
        term2 = self.eps_m() #эв
        term3 = log(self.eps_m()/self.I)
        return abs(term1*term3/term2   *3.8e12) *1000 #для того чтобы в эв а не кэв считать ибо косячно
    
    #Расширяю свой класс --функциями которые позволят установить энергию для каждой глубины
    
    def interpolated_dEdz(self,val):
        #эв на ангстрем
        x = np.linspace(4.e-4,1,1000)
        mas = []
        for i in x:
            mas.append(Ion_Class(i*E_init_Fe,Fe_ion,Fe).dEdz())
        return interp1d(x,mas)(val)
    
    def coordinat_as_function_from_energi(self, E_finish):
        f = lambda x: 1./self.interpolated_dEdz(x)
        return scipy.integrate.quad(f, E_finish,1.)[0] * self.E1
    
    def interpolated_coordinate_from_energi(self):    #поменять название на правильное координата от энергии!!!
        """
        Область определеня функции:
                    val is the part of non dimention energi from 0.001 to 1
        
        возвращает:
            функцию которая интерплированна
        """
        a = Ion_Class(E_init_Fe,Fe_ion,Fe)
        part = np.linspace(0.001,1.,10)    #примерная точка остановки 58.57 мкм
        mas_cord = []
        for i in part:
            mas_cord.append(a.coordinat_as_function_from_energi(i)*1e-4) #для микрометров

        #plt.plot(part,mas_cord) #выводит графиек функции с которой наинтерволировали
        return interp1d(part,mas_cord)
    
    def velocity_how_function_from_E(self,E):
        """
        E-безразмерное от 4e-4 to 1
        """
        v = sqrt(2*E*self.E1*self.ev_to_erg/self.m_ion)  #см сек-1
        return v*1e-4 #мкм сек-1
        #visualisation of this function
        #mas =[]
        #for i in np.linspace(0.001,1,100):
        #    mas.append(a.velocity_how_function_from_E(i))
        #plt.plot(np.linspace(0.001,1,100), mas)
    
    def coordinate_from_time(self):
        """
        time begin of process
        return time in seconds
        """
        NUM = 100
        mas_E = np.linspace(1.,0.001, NUM)
        mas_x = []
        mas_v = []
        mas_t = [0,]
        t = 0     
        fun = self.interpolated_coordinate_from_energi()
        for i in mas_E:
            a = time.time()
            mas_x.append(fun(i))
            mas_v.append(self.velocity_how_function_from_E(i))
            b = time.time()
            print(b - a)
        print(mas_x,'\n',mas_v)   
        for i in range(NUM - 1):
            dx = mas_x[i+1] - mas_x[i]
            avg_v = (mas_v[i] + mas_v[i+1])/2. #чтобы был нужный знак
            dt = dx/avg_v
            print('dx = ',dx)
            print('v = ',avg_v)
            print('dt = ',dt)
            
            t += dt
            mas_t.append(t)
        return (mas_t, mas_x)
        