#导入自己编写的算法，验证程序是否正确运行
from Method.Newton import Newton as nt #导入牛顿迭代法
from Method.classGauss import ClassGauss #导入高斯列主元法
from Method.NewtonInsert import NewtonInsert #导入牛顿插值多项式解法
from Method.Integrator import Integrator #导入数值积分方法
from Method.Insert import Insert #导入线性插值和二次插值
from Method.Polynomial import Polynomial #导入通用多项式拟合算法
from Method.RungeKutta4 import RungeKutta4 #导入RungeKutta4算法
#导入库函数，对计算结果进行验证
from sympy import * #导入符号运算库，对计算结果进行验证
import numpy as np	#矩阵运算支持库
import math  #导入常用的数学函数，如exp(x)
import matplotlib.pyplot as plt #导入绘图函数

if __name__ == "__main__":

    """
    print("------牛顿迭代法调用的例子-------")
    x = symbols("x") #定义符号变量x
    fx = x**4 - 1.4*x**3 - 0.48*x**2 + 1.048*x-0.512 #定义需要求解的方程
    mynewton = nt(fx,1,x,[1,10]) #初始化牛顿迭代法
    result = mynewton.cal() #调用牛顿迭代法求解方程
    print("方程的根为：%f\n求解精度为:%f\n"%(result[0],result[1]))
    print("使用sympy.solve()求解方程对比求解是否正确")
    print(solve(fx,x))
    

    print("----高斯列主元法求解线性方程组的例子----")
    a= [ [ 1, 0, 0, 0, 1],
         [-1, 1, 0, 0, 1],
         [-1,-1, 1, 0, 1],
         [-1,-1,-1, 1, 1],
         [-1,-1,-1,-1, 1]]#系数矩阵
    b = [0.1,10,1,5,4]#
    
    g = ClassGauss(a,b)#初始化高斯列主元法
    g.Guass()#调用方法求解兴县方程组的解
   
    
    print("----牛顿法插值多项式求解的例子----")
    sr_x = np.linspace(0,3,6)
    sr_fx = [sin(i)+1 for i in sr_x] #构造sin(x)+1函数的插值点
    newtonInsert = NewtonInsert(sr_x,sr_fx)#初始化牛顿插值法
    fx = newtonInsert.get_Newton_inter()
    tmp_x =  np.linspace(-1,4,40)         #构造插值函数的绘图点
    tmp_y = [fx(i) for i in tmp_x] # 根据插值函数获得测试用例的纵坐标
    print("-----插值多项式的系数----")
    print(newtonInsert.a)#打印插值多项式的系数

	#牛顿法插值多项式画图
    plt.figure("牛顿法插值多项式绘图")
    ax1 = plt.subplot(111)
    plt.sca(ax1)
    plt.plot(sr_x, sr_fx, linestyle = '', marker='o', color='b')
    plt.plot(tmp_x, tmp_y, linestyle = '--', color='r')
    plt.show()
    
    
    print("-------分段线性插值和分段二次插值的例子------")
    def f(x):#被插值函数表达式
        return sin(x)
    for i in range(4):
        print("等分区间数为%d时"%2**(i+1))
        myIns = Insert(f,[0,1],2**(i+1))#构造插值误差计算实例
        print("线性插值误差||f(x)-P1n||2为：")
        print(myIns.LinearInsert())#调用线性插值误差计算，计算插值误差函数的2范数
        print("线性插值误差||f(x)-P2n||2为：")
        print(myIns.QuadraticInsert())#调用二次插值误差计算，计算插值误差函数的2范数
    
    print("----------通用多项式函数拟合-----------")
    s = 0.000000001
    x = [1-2*s,1-s,1,1+s,1+2*s]
    y = [1,2,3,4,5]
    mypoly = Polynomial(x,y,3)
    print("构造的正交多项式为：")
    print(mypoly.showMultinomial())
    print("正交多项式系数为：")
    print(mypoly.a)
    fit = mypoly.fitFunction()#获取拟合函数的地址，调用格式 y = fit(x)
    #计算实际误差
    err = np.zeros(len(y))
    for i in range(len(y)):
        err[i] = y[i] - fit(x[i]) 
    print("实际误差为：")
    print(err)
    #绘制被拟合的散点图和拟合后的曲线图
    plt.figure("通用多项式拟合绘图")
    ax1 = plt.subplot(111)
    plt.sca(ax1)
    plt.scatter(x, y)
    tmp_x = np.linspace(x[0],x[-1],20)
    tmp_y = np.zeros(20)
    for i in range(len(tmp_x)):
        tmp_y[i] = fit(tmp_x[i])
    plt.plot(tmp_x, tmp_y, linestyle = '--', color='r')
    plt.show()
    """
    
    print("-------Runge-Kutta4阶算法例子------")
    #-----@brief:  微分方程描述函数
    #-----@return  返回值为微分方程
    def f(x,y):
        return -3*y

    myRungeKutta = RungeKutta4(f,1,[0,2],10) #调用Runge-kutta算法
    rk_x,rk_y = myRungeKutta.get_result()   #获取微分方程求解结果
    fxreal = np.zeros(len(rk_x))
    error = np.zeros(len(rk_x)) #构造误差存储器
    #计算微分仿真准确解
    for _ in range(len(rk_x)):
        fxreal[_] = exp(-3*rk_x[_])
    error[:] = fxreal[:] - rk_y[:]#求解误差相量
    #计算误差的2范数
    Err = 0.0
    for _ in error:
        Err = abs(_)**2
    print("***误差****")
    print(sqrt(Err))#输出误差的2范数
    #绘制Runge-Kutta4算法求解的微分方程解和微分方程准确解以及误差变化图像
    plt.figure("Runge-Kutta4阶绘图")
    ax2 = plt.subplot(121)
    plt.sca(ax2)
    plt.plot(rk_x,rk_y, linestyle = '-', marker='o', color='b')
    plt.plot(rk_x,fxreal, linestyle = '-', marker='', color='r')
    ax3 = plt.subplot(122)
    plt.sca(ax3)
    plt.plot(rk_x,abs(error), linestyle = '-', marker='', color='r')
    plt.title("error")
    plt.show()
    
