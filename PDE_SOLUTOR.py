# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 14:06:03 2017

@author: Alpharius

SOLUTOR
    |
    |---explisit_differential_scheme
    |---A
    |---create_df_for_temp
    |---save
"""
import numpy as np
import pandas as pd
from math import *
from Electric_fild import *
from Ionclass import *
    

class SOLUROR(Phisic_priperties):
    """
    решает систему уравнений теплопроводности
    и выводит данные для 3д картинки.
    """
    def __init__(self, grid, init_T, bound_condition, phy_property   ):
        """
        создание сетки с нулевыми значениями
        формирование начального и 2 краевых условий
        введение в оборот физических параметров для изучаемого образца
        
        вводим:
            grid - кортедж из 2 массивов отвечающих координатам и времени
                координаты будут безразмерными (размерность введу по ходу)
                grid[0]-координаты
                grid[1]-время
                grid[2]--тип данных
                
            в иоге у нас прямоугольная матрица в которой строки это координаты
            а столбцы времена
            
            init_T - массив начальных пемператур
        
        комментарий:
            все расчеты проходят в СГС и ни как иначе
            
            координаты и времена будут безразмерными но с учетом ножителей
            scale_t, scale_coord;
            
            в сетке столбцы это разные моменты времени а строки это разные 
            координаты. 
                1.первый столбец это начальные условия
                2.первая строка это граничные условия на оси трека
                3.последняя строка это граничное условие в немозмущенном удалении от оси трека
        """
        #сетка
        self.GRID = np.zeros((grid[0].shape[0],grid[1].shape[0]), dtype = grid[2])
        self.GRID[:,0] = init_T[:]   # их размеры равны
        self.GRID[0,:] = bound_condition[0][:]
        self.GRID[grid[0].shape[0],:] = bound_condition[1][:]
        #термодинамические функции
        self.C_e = phy_property[0]
        self.C_i = phy_property[1]
        self.K_e = phy_property[2]
        self.K_i = phy_property[3]
        self.g = 0.01                #электрон фононное взаимодействие мне не известно!!!
        #физические масштабы решаемой задачи
        self.skale_t = phy_property[4]              
        self.scale_coord = phy_property[5]
        
        self.CURRANT = 0.1
        
    def explisit_differential_scheme(self,time = 1, dx = 1):
        """
        time -  безразмерная постоянная для измерения долей масштаба 
        (вариируем для устойчивости)
        dx - то же самое но для координаты
        
        те это расстояния между соседними узлами сетки в безразмерных параметрах
        вся размерность учтена в self.skale...
        """
        deltha_t = self.skale_t * time          
        deltha_x = self.scale_coord * dx
        for j in range(0,self.GRID.shape[1]-1):  #time
            for i in range(1,self.GRID.shape[0]-1): #space
                #половинчатые температуры
                T_n_minus_half = 0.5*(self.GRID[i-1,j] + self.GRID[i,j])  
                T_n_plus_half = 0.5*(self.GRID[i+1,j] + self.GRID[i,j])
                
                #температура Иона
                self.GRID[i,j+1]['i'] = self.GRID[i,j]['i'] + \
                    deltha_t/(self.C_i(self.GRID[i,j]['i'])*deltha_x**2) *\
                    (
                    self.K_i(T_n_minus_half)*(self.GRID[i,j]['i'] - self.GRID[i-1,j]['i'])\
                    -\
                    self.K_i(T_n_plus_half_half)*(self.GRID[i+1,j]['i'] - self.GRID[i.j]['i'])\
                    ) \
                    +\
                    (self.g*deltha_t)/(self.C_i(self.GRID[i,j]['i'])*deltha_x) *\
                    (self.GRID[i,j]['e'] - self.GRID[i,j]['i'])
                #температура Электрона
                self.GRID[i,j+1]['e'] = self.GRID[i,j]['e'] + \
                    deltha_t/(self.C_e(self.GRID[i,j]['e'])*deltha_x**2) *\
                    (
                    self.K_e(T_n_minus_half)*(self.GRID[i,j]['e'] - self.GRID[i-1,j]['e'])\
                    -\
                    self.K_e(T_n_plus_half_half)*(self.GRID[i+1,j]['e'] - self.GRID[i.j]['e'])\
                    ) \
                    -\
                    (self.g*deltha_t)/(self.C_e(self.GRID[i,j]['e'])*deltha_x) *\
                    (self.GRID[i,j]['e'] - self.GRID[i,j]['i']) \
                    +\
                    deltha_t/(self.C_e(self.GRID[i,j])*deltha_x)# *\
                    #self.A(deltha_x*i, deltha_t*j)
                    
                    
    def A(self,r,t):
        """
        Эффективный источник энергии в электронную подсистему
        """
        r_0 = 2e-7                          #характерная область существования поля
        t_0 = 1e-15                         #характерное время термализации электронов
        C = 1./(4*sqrt(2)*pi*r_0**2*t_0)   #нормировочная кончстанта
        if t>2e-15:
            return 0
        else:
            
            
    def create_df_for_temp(self):
        """
        Сохранение нашей расчитанной сетки (сохранена вся история эволюции 
        температурных профилей) в экселевском файле для дальнейшего удобного
        воспроизведения.
        """
        self.df = pd.DataFrame(self.GRID)
        self.df.to_csv('temperatur evalution')
                
                        
                
                
            
    def save(self):
        """
        наша прямоугольная матрица должна быть сохранена в csv файл с использованием
        pandas
        """
        pass
        
        
class Experimental_cls():
    self.phy_property = (1,2,3,4)
    self.bound_condition = ([],[])
    self.grid
        
        
        
        
        