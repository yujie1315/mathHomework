import numpy as np

class ClassGauss(object):
    '''
    高斯列主元法
    @brief:   初始化高斯列主元法
    @param:   a   线性方程组的系数矩阵                                                 
    @param:   b   方程组右端矩阵  
	@param:   n   方程组的维数	                                                                          
    '''
    def __init__(self, a, b):	#初始化高斯列主元法
        super(ClassGauss, self).__init__()
        self.a = np.array(a)#赋值系数矩阵
        self.b = np.array(b)#方程组右端矩阵
        self.n = len(self.b)#计算方程维数
        self.x = np.zeros(self.n)#构造方程解容器
    '''
    @brief:   寻找列最大值和列最大值所在行号
    @param:   j   需要寻找最大值的系数矩阵列号
	@return： 返回寻找列的最大值和最大值所在行号
    '''
    def max(self,j):
        max = abs(self.a[j][j])
        maxLine = j
        for i in range(j,self.n-1):
            if max < abs(self.a[i][j]):
                max = abs(self.a[i][j])
                maxLine = i
        return max,maxLine
    '''
    @brief:   根据主元选择交换方程
    @param:   j         需要交换的行号
	@param：  maxLine   列主元所在行号
    '''
    def swapLine(self,maxLine,j):
        self.a[[j,maxLine],:] =self.a[[maxLine,j],:]#进行系数矩阵行交换
        self.b[j],self.b[maxLine] = self.b[maxLine],self.b[j]#交换b
        print("SWAP***************")#输出交换后的结果
        print(self.a)
        print(self.b)

    '''
    @brief:   采用高斯列主元法求解线性方程组的解
    @notice:  直接打印线性方程组的解
    '''
    def Guass(self):#进行高斯列主元法求解线性方程
        maxLine = 0
        maxValue = self.a[0][0]
        for j in range(0,self.n):
            maxValue,maxLine = self.max(j)#寻找主元
            if maxValue == 0:
                raise ValueError('no unique solution')#方程奇异，无解退出
            if maxLine != j:#根据主元选择进行方程交换
                self.swapLine(maxLine,j)
            '''-----------进行消元------------'''
            for i in range(j+1,self.n):
                l = self.a[i][j] / self.a[j][j]
                self.b[i] -= l * self.b[j]
                for p in range(0,self.n):
                    self.a[i][p] -= l*self.a[j][p]
        '''----------打印消元后的系数矩阵---------'''
        print("Guass************")
        print(self.a)
        '''-----------回带求解线性方程组-----------'''
        self.x[self.n-1] = self.b[self.n-1]/self.a[self.n-1][self.n-1]
        for i in range(self.n-2,-1,-1):
            sum = 0.0
            for j in range(i,self.n):
                sum += self.a[i][j]*self.x[j]
            self.x[i] = (self.b[i]-sum)/float(self.a[i][i])
        print("线性方程组的解为：")
        print(self.x)
            