import numpy as np
from math import sqrt
from Method.Integrator import Integrator #导入数值积分
class Insert(object):
    '''
    @brief:   计算线性插值和二次插值误差类
    @param:   debug 调试标志                                                
    '''
    debug = False
    '''
    @brief:   类初始化函数
    @param:   fx 被插值函数地址
    @param    x  隔根区间  x = [a,b]
    @param    n  等分区间数
    '''
    def __init__(self,fx,x,n):
        self.F = fx
        self.x = np.linspace(x[0],x[1],n+1)
        self.y = np.zeros(len(self.x))
        for i in range(len(self.x)):
            self.y[i] = self.F(self.x[i])
        if self.debug:
            print(self.x)
            print(self.y)
    '''
    @brief:   线性插值函数
    @param  xi 插值区间xi =[xk,xk+1]
    @param  yi 插值函数值yi = [yk,yk+1]
    @param  x  插值函数自变量
    @return 返回插值函数在插值区间xi上在x点处的值
    '''
    def P1n(self,xi,yi,x):
        return yi[0]*(x-xi[1])/(xi[0]-xi[1]) + yi[1]*(x-xi[0])/(xi[1]-xi[0])
    '''
    @brief:   二次插值函数
    @param  xn 插值区间xi =[xk,xk+1]
    @param  yn 插值函数值yi = [yk,yk+1]
    @param  x  插值函数自变量
    @return 返回插值函数在插值区间xi上在x点处的值
    '''
    def P2n(self,xn,yn,x):
        xi = [xn[0],xn[0]/2+xn[1]/2,xn[1]]
        yi = [yn[0],self.F(xn[0]/2+xn[1]/2),yn[1]]
        return yi[0]*(x-xi[1])*(x-xi[2])/(xi[0]-xi[1])/(xi[0]-xi[2]) + \
             yi[1]*(x-xi[0])*(x-xi[2])/(xi[1]-xi[0])/(xi[1]-xi[2]) + \
             yi[2]*(x-xi[0])*(x-xi[1])/(xi[2]-xi[0])/(xi[2]-xi[1])
    '''
    @brief:   计算线性插值误差 ||f(x)-P1n||2 = sqrt(Integrator{[f(x)-P1n]*[f(x)-P1n])}积分区间为【a,b】)
    @return   返回误差函数的2范数
    '''
    def LinearInsert(self):
        '''
        @brief:线性插值误差函数
        '''
        def err(x):
            return self.F(x) - self.P1n(xi,yi,x)
        sum = 0
        for i in range(len(self.x)-1):#分区间计算误差函数的内积,便于后面计算误差函数的2范数
            xi = [self.x[i],self.x[i+1]]
            yi = [self.y[i],self.y[i+1]]
            if self.debug:
                print(xi)
                print(yi)
                print(sum)
            sum += self.fanshu(err,xi)
        return sqrt(sum)#返回误差函数的2范数

    '''
    @brief:   计算线性插值误差 ||f(x)-P2n||2 = sqrt(Integrator{[f(x)-P2n]*[f(x)-P2n])}积分区间为【a,b】)
    @return   返回误差函数的2范数
    '''  
    def QuadraticInsert(self):
        '''
        @brief:二次插值误差函数
        '''
        def err(x):
            return self.F(x) - self.P2n(xi,yi,x)
        sum = 0
        for i in range(len(self.x)-1):#分区间计算误差函数的内积,便于后面计算误差函数的2范数
            xi = [self.x[i],self.x[i+1]]
            yi = [self.y[i],self.y[i+1]]
            if self.debug:
                print(xi)
                print(yi)
                print(sum)
            sum += self.fanshu(err,xi)
        return sqrt(sum)#返回误差函数的2范数

    '''
    @brief:  计算函数在区间x=[a,b]上的内积,函数命名写错了，改太麻烦了，，，
    @param： fx  被计算函数的地址
    @param： x   被计算的区间x=[a,b]
    @return  返回函数内积值
    '''
    def fanshu(self,fx,x):
        '''
        @brief:  函数相乘
        '''
        def f(x):
            return fx(x)*fx(x)
        myInt = Integrator(f,x)
        return myInt.cal()