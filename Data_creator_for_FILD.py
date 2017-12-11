# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 20:03:27 2017

@author: Alpharius

Electric_fild
    |
    |---create_df
    |---read_csv
    |---create_plot
    |

"""
from Electric_fild import *


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
        По завершению функци у нас появляется таблица (поле в ращличных
        координатах для разных времен . По  строкам время а столбцы расстояние)
        --это поможет избежать
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