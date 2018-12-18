import numpy as np

class RungeKutta4(object):
    """Runge-Kutta 4阶微分方程求解法"""
    '''
    @brief:   初始化Runge-Kutta 算法
    @param:   F  一阶微分方程描述函数的地址                                                
    @param:   y0 微分方程求解初始值
    @param:   x 微分方程求解区间 x = [a,b]
    @param:   n  微分方程求解次数  
    '''
    def __init__(self,fx,y0,x,n):
        self.F = fx
        self.n =n
        self.h  = (x[1]-x[0])/n #计算求解步长
        self.x  = np.linspace(x[0],x[1],n+1)#构造横坐标序列
        self.y  = np.zeros(len(self.x))#申请存错解的内存
        self.y[0] = y0#赋值解初值
    '''
    @brief:  实现Runge-Kutta 算法,为算法类内部函数，不建议外部直接调用
    '''
    def cal(self):
        for i in range(1,self.n+1,1): 
            K1=self.F(self.x[i-1],self.y[i-1])
            K2=self.F(self.x[i-1]+self.h/2,self.y[i-1]+(self.h*K1)/2)
            K3=self.F(self.x[i-1]+self.h/2,self.y[i-1]+(self.h*K2)/2)
            K4=self.F(self.x[i-1]+self.h,self.y[i-1]+self.h*K3)
            self.y[i]=self.y[i-1]+self.h*(K1+2*K2+2*K3+K4)/6
    '''
    @brief:  Runge-Kutta 算法，结果输出函数，调用此函数输出微分方程解
    @return  返回值为微分方程解的横坐标和纵坐标
    '''
    def get_result(self):
        self.cal()
        return self.x,self.y