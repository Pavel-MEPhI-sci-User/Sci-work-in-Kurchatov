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

#a = Electric_Field(1e-18,1, E_init_Fe,Fe_ion, Fe)        
#print(a.E())
"""
mas_q_til = np.linspace(0.,4, 6)
mas_t = np.linspace(1., 250, 40)
mas_f_from_t = []
mas_f = []
file = open('table.txt', 'w')
for j in mas_q_til:
    mas_f = []
    for i in mas_t:
        a = Electric_Field(1e-18 * i,j, E_init_Fe,Fe_ion, Fe)  
        f = a.E()
        print(i,j,f)
        mas_f.append(f)         #массив полей в конкретный момент времени по разным координатам
    mas_f_from_t.append(mas_f)  #массив полей в различные времена и в разных координатах
                                #по строкам идут разные времена по столбцам разные координаты 
file.write(mas_f_from_t)
plt.plot(mas_t, mas_f)
"""

class Data_creator_for_FILD():
    """
    этот класс предназначен для того чтобы:
        a.  осуществлять расчет электирического поля в различных координатах
            и временах
        b.  сохранять результаты в csv файл чтобы потом не производить долгих расчетов 
            а просто загружать необходимое
        c.  создавать анимацию которая как мне кажется будет полезна в моем дипломе
        
    здесь применяется фреймворк пандас
    """
    def __init__(self,mas_q_til, mas_time):
        """
        конструктор принимающий в себя массивы отвечающие за расчет поля
        в определенных координатах в определенное время-первый аргумент
        координаты второй время
        """
        self.mas_q = mas_q_til
        self.mas_t = mas_time
        self.mas_index =['time = ' + str(1e-18 * i) for i in self.mas_t]  #для заполнения таблицы
        self.mas_col = [str(i) for i in self.mas_q]   #для заполения колонок
    def create_df(self):
        """
        функция создающая эту таблицу --в ней мы вызываем экземпляр другого
        класса и при помощи метода E() расчитываем поле. которое заносим в 
        таблицу. 
        По завершению функци у нас появляется таблица --это поможет избежать
        длительных расчетов в дальнейшем.
        """
        mas_f_from_t = []    
        for i in self.mas_t:
            mas_f = []
            for j in self.mas_q:
                a = Electric_Field(1e-18 * i,j, E_init_Fe,Fe_ion, Fe) 
                f = a.E()
                print(i,j,f)
                mas_f.append(f)
            mas_f_from_t.append(mas_f)
        self.df = pd.DataFrame(mas_f_from_t, index = self.mas_index, 
                               columns = self.mas_col)
        self.df.to_csv('El_field.csv')
        
    def read_csv(self):
        """
        здесь мы просто считываем эту таблицу и вводим е в свой оборот
        """
        self.df = pd.read_csv('El_field.csv')
        
    def create_plot(self,time, use_save_file = True):
        """
        создание одномерной картинки (в момент времени time) и анимации-
        показывающей эволюцию электрического поля со временем
        принимает:
            time - целочисленный параметр - номер ячейки во временном
            по поводу второго аргумента это просто элемент гибкости-
            для возможного расширения функционала
        """
        if use_save_file:
            self.read_csv()
        fig = plt.figure()
        ax = plt.axes(xlim=(0, 5), ylim=(0, 50))
        line, = ax.plot([], [], lw=2)
        
        def init():
            line.set_data([], [])
            return line,
        def animate(i):
            x = self.mas_q
            y = self.df.iloc[i].values[1:]
            line.set_data(x, y)
            return line,
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=20, interval=100, blit=True)
        anim.save('basic_animation.html', fps=30, extra_args=['-vcodec', 'libx264'])

        plt.show()
if __name__ == "__main__":
    a = Data_creator_for_FILD(np.linspace(0,3,5), np.linspace(1,200,20))
    a.create_plot(2)


