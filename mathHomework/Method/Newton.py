from sympy import diff #导入方程求导方法

class Newton(object):

    '''
    @brief:   使用牛顿法求解非线性方程,牛顿迭代法的参数
    @param:   __MAXSTEP  牛顿迭代法最大迭代次数                                                
    @param:   __ret      牛顿法迭代求解精度	  
    @param    debug      程序调试标志
    '''
    __MAXSTEP = 1000
    __ret = 0.000001
    debug = False

    '''
	@brief:   初始化牛顿法牛顿迭代法的参数
    @param:   func  需要使用牛顿迭代法求解的方程
    @param:   x0    牛顿法迭代的初始值
    @param:   x     被求解方程的自变量，符号变量
	@param:	  xLim  牛顿迭代法的隔根区间,长度为2的列表表示，如 xLim = [a,b]
    @param:   self 表示这个类所拥有的函数，python 语法决定的
    '''
    def __init__(self, func, x0, x,xLim):
        self.func = func
        self.dfunc = diff(func,x)#求解方程的导数
        print(self.dfunc)
        self.x0 = x0
        self.x = x
        self.xmin = xLim[0]
        self.xmax = xLim[1]
    '''
	@brief:   使用牛顿法求解非线性方程
	@return： 返回求解方程的根和方程求解精度
	@notice： 求解过程中输出迭代次数和迭代求解的值
    '''
    def cal(self):
        stepCount = 0
        while stepCount < self.__MAXSTEP:#迭代循环
            x1 = self.x0 - self.func.subs(self.x,self.x0)/self.dfunc.subs(self.x,self.x0)#牛顿法迭代格式x1 = x0-f(x)/f'(x)
            stepCount += 1
            if self.debug:
                print("迭代次数为：%d,本次迭代解为%f,前一次迭代解为%f"%(stepCount,x1,self.x0))
            if abs(x1-self.x0)<self.__ret:
                return [x1,self.__ret]#返回求解方程的根和方程求解精度
            if x1<self.xmin:
                x1 = self.xmin
            if x1 > self.xmax:
                x1 = self.xmax
            self.x0 = x1