import numpy as np

class NewtonInsert(object):
    __debug = False
    """
	@brief:	   初始化牛顿法插值多项式
	@param:   xi   所有插值节点的横坐标集合                                                 
    @param:   fi   所有插值节点的纵坐标集合
	@param:   a    所有插值多项式的系数
	@param:   Px   存储多项式插值多项式函数的地址
	"""
    def __init__(self,xi,fi):
        self.xi = xi
        self.fi = fi
        self.a = np.zeros(len(xi))
        self.a[0] = fi[0]
        if self.__debug:
            print("节点数为：%d"%len(xi))

    """
    @brief:   计算n阶差商 f[x0, x1, x2 ... xn]，采用递归的方法，未考虑牛顿插值法的优点，计算较慢
    @param:   xi   所有插值节点的横坐标集合                                                 
    @param:   fi   所有插值节点的纵坐标集合                                                      
    @return:  返回xi的i阶差商(i为xi长度减1)                                                     
    @notice:  a. 必须确保xi与fi长度相等                                                       
              b. 使用递归的方法求解n阶差商.                                         
              c. 递归减递归(每层递归包含两个递归函数), 每层递归次数呈二次幂增长                                                                                
    """
    def get_order_diff_quot(self, xi = [], fi = []):
        if self.__debug:
            print("求%d阶差商，，，"%(len(xi)-1))
            print(xi)
            print(fi)

        if len(xi) > 2 and len(fi) > 2:
            return (self.get_order_diff_quot(xi[:len(xi)-1], fi[:len(fi)-1]) - self.get_order_diff_quot(xi[1:len(xi)], fi[1:len(fi)])) / float(xi[0] - xi[-1])
        return (fi[0] - fi[1]) / float(xi[0] - xi[1])

    """
    @brief:  获得Wi(x)函数;
             Wi的含义举例 W1 = (x - x0); W2 = (x - x0)(x - x1); W3 = (x - x0)(x - x1)(x - x2)
    @param:  i  i阶(i次多项式)
    @param:  xi  所有插值节点的横坐标集合
    @return: 返回Wi(x)函数
    """
    def get_Wi(self, i = 0):
        def Wi(x):
            result = 1.0
            for each in range(i):
                result *= (x - self.xi[each])
            return result
        return Wi
     
    """
    @brief:   获得牛顿插值多项式函数  y=P(x)
    @param:   xi   所有插值节点的横坐标集合
	@param:   fi   所有插值节点的纵坐标集合
	@return： 返回牛顿法插值多项式函数的地址
    """
    def get_Newton_inter(self):
        def Newton_inter(x):
            result = self.fi[0]
            for i in range(2,len(self.xi)):
                if self.__debug:
                    print("求插值多项式第%d次项"%i)
                self.a[i] = self.get_order_diff_quot(self.xi[:i], self.fi[:i])
                if self.__debug:
                    print("打印系数**********************")
                    print(self.a[i])
                #get_Wi()返回的是函数的地址，get_Wi(i-1, xi)(x)表示调用返回的函数
                result += (self.a[i] * self.get_Wi(i-1)(x))
            return result
        return Newton_inter