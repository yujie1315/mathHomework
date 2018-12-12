import numpy as np
class Polynomial(object):
    """
    @brief 通用多项式拟合程序
           构造的多项式为f(x) = a0*P(x,0) +a1*P(x,1)+a2*P(x,2)...
           其中P(x,n)为正交多项式
    @param debug 程序调试标志
    """
    debug = False
    """
    @brief  多项式拟合程序初始化
    @param  xi  输入散点序列的横坐标
    @param  yi  输入散点序列的纵坐标
    @param  n   拟合多项式的最高次数
    """
    def __init__(self,xi,yi,n):
        self.xi = xi
        self.yi = yi 
        self.n  = n+1
        self.a  = np.zeros(n+1);
        self.nomial = []
        if len(xi) != len(yi):#判断输入横纵坐标长度是否相同
            raise ValueError("输入参数不正确")
        self.coefficientCal()
    """
    @brief 输出显示拟合的正交多项式
           输出为字符串形式，程序不能调用
    """
    def showMultinomial(self):
        print("正交多项式为：")
        for k in range(len(self.a)):
           if k==0:
               self.nomial.append("1")
           elif k==1:
               self.nomial.append("(x - %f)"%(self.A(1)))
           else:
               self.nomial.append("(x-%f)*"%self.A(k)+self.nomial[k-1]+"-%f*"%self.B(k-1)+self.nomial[k-2])
        for _ in self.nomial:
            print(_)
    """
    @brief   利用递归的方式构造正交多项式
    @param   x 正交多项式的自变量，传入的是具体的数值
    @param   k 正交多项式的次数
    @return  返回正交多项式的计算结果
    """
    def P(self,x,k):
        if k==0:
            return 1
        elif k==1:
            return (x - self.A(1))
        else:
            return ((x - self.A(k))*self.P(x,k-1) - self.B(k-1)*self.P(x,k-2))
    """
    @brief  返回构造多项式系数alpha()的值 
            alpha() = (xi*P(k-1),P(k-1))/(P(k-1),P(k-1))
    @param  k  正交多项式的次数
    @return 返回值构造的正交多项式的系数alpha
    """
    def A(self,k):
        num = 0
        den = 0
        for i in range(len(self.xi)):
            num += self.xi[i]*self.P(self.xi[i],k-1)**2
            den += self.P(self.xi[i],k-1)**2
        if self.debug:
            print("A(%d)为："%k)
            print(num/den)
        return num/den
    """
    @brief  返回构造多项式系数beta()的值 
            beta() = (xi*P(k),P(k))/(P(k-1),P(k-1))
    @param  k  正在构造的正交多项式的次数 减1
    @return 返回值构造的正交多项式的系数beta
    """
    def B(self,k):
        num = 0
        den = 0
        for i in range(len(self.xi)):
            num += self.P(self.xi[i],k)**2
            den += self.P(self.xi[i],k-1)**2
        if self.debug:
            print("B(%d)为："%k)
            print(num/den)
        return num/den
    """
    @brief   输出拟合多项式的系数a0,a1,a2...
             f(x) = a0*P(x,0) +a1*P(x,1)+a2*P(x,2)...
             其中P(x,n)为正交多项式
    """
    def coefficientCal(self):
        for _ in range(len(self.a)):
            num = 0
            den = 0
            for i in range(len(self.xi)):
                num += self.yi[i]*self.P(self.xi[i],_)
                den += self.P(self.xi[i],_)*self.P(self.xi[i],_)
            self.a[_] = num/den
        if self.debug:
            print("多项式系数为："   )
            print(self.a)
    """
    @brief  拟合的多项式函数
    @return 返回拟合多项式的地址
    """
    def fitFunction(self):
        """
        @brief  拟合的多项式函数计算，嵌套内部函数
        @param  x 需要计算的拟合多项式的横坐标
        @return 返回拟合多项式的计算纵坐标的结果
        """
        def fitfun(x):
            sum = 0
            for _ in range(len(self.a)):
                sum += self.a[_]*self.P(x,_)
            return sum
        return fitfun