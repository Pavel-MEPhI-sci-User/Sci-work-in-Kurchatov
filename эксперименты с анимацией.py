import numpy as np

dt = np.dtype([('electron','f8'),('proton','f8')])
mas = np.array([(1,2),(3,4)], dtype = dt)
print(mas['proton'])

x, y = np.meshgrid([1,2,3],[10,20,30])
print(x,'\n',y)


class A():
    def f(self):
        print ('hello gandon')
        
    def g(self, arg):
        print(arg)
        self.f()
        return 0
    
    
a = A()
a.g(4)
a.f()