class Integrator(object):
    """利用复化梯形公式计算数值积分"""
    '''                                              
    @param:   __MAXSTEP 数值积分最大迭代次数，最大细分区间为2**__MAXSTEP等份
    @param:   __ret     数值积分迭代误差
    @param:   debug     程序书写过程中的调试标志  
    '''
    __MAXSTEP = 20
    __ret = 0.000001
    debug = False
    '''
    @brief:   初始化梯形公式数值积分算法
    @param:   fx   被积函数地址,调用被积函数格式为fx(x)
    @param    x    函数被积分区间，x = [a,b]

    '''
    def __init__(self,fx,x):
        self.F = fx
        self.a = x[0]
        self.b = x[1]

    '''
    @brief:   梯形公式数值积分算法
    @return   T2n 返回被积函数数值积分结果

    '''
    def cal(self):
        hn = (self.b -self.a)/2 #初始化细分区间长度
        Tn = 0
        T2n = 0
        Tn = hn/2*(self.F(self.a)+self.F(self.b)+2*self.F(self.a+hn))# 计算二等分区间数值积分值
        step = 1
        while step < self.__MAXSTEP:
            h2n = hn/2 #细分区间长度
            sum = 0   
            for i in range(2**step):
                sum = sum + self.F(self.a+h2n+hn*i)
            T2n = Tn/2 + hn/2*sum  #梯形数值积分迭代公式 T2n = Tn/2 + hn/2*sum(f(x_k+1/2))
            if abs(T2n-Tn)<self.__ret:
                break;
            if self.debug:
                print("第%d次迭代"%step)
                print("Tn = %10f"%Tn)
                print("T2n = %10f"%T2n)
            hn = h2n
            Tn = T2n
            step +=1
        return T2n