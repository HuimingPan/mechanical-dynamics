
import matplotlib.pyplot as plt
from math import e,sqrt,atan2,pi
import numpy as np
#math包里的sin函数只接受一个标量，在程序中，使用np.sin代替math.sin

def VTB1(m,c,k,x0,v0,tf):
    '''
        VTB1用来计算单自由度自由振动系统的响应；
        VTB1绘出单自由度有阻尼自由振动系统的相应图
        m：mass；
        c:阻尼；
        k：刚度；
        x0:初始位移；
        v0:初始速度；
        tf：仿真时间
        程序中：
        zeta为阻尼比；wn为固有频率；A为振动幅度；phi为初相位；
    '''
    wn=sqrt(k/m)
    xi=c/2/m/wn
    wd=wn*sqrt(1-xi**2)
    print("固有频率为{} rad/s\n阻尼比为{}\n有阻尼的固有频率为{:.5}".format(wn,xi,wd))
    t=np.arange(0,tf,tf/1000)
    if xi<1:#弱阻尼状态
        A=sqrt(((v0+xi*wn*x0)**2+(x0+wd)**2)/wd**2)
        phi=atan2(x0*wd,v0+xi*wn*x0)
        x=A*e**(-xi*wn*t)*np.sin(wd*t+phi)
        print("A={:.3}".format(A))
        print("phi={:.3}".format(phi))
    elif xi==1:#临界阻尼状态
        c1=x0
        c2=v0+wn*x0
        print("c1={:.3}".format(c1))
        print("c2={:.3}".format(c2))
        x=e**(-wn*t)*(c1+c2*t)
    else :#过度阻尼状态
        c1=(v0+(-xi+sqrt(xi*2-1))*wn*x0)/(2*wn*sqrt(xi*2-1))
        c2=(v0+(xi+sqrt(xi*2-1))*wn*x0)/(2*wn*sqrt(xi*2-1))
        print("c1={:.3}".format(c1))
        print("c2={:.3}".format(c2))
        x=e**(-wn*t)*(c1+c2*t)
    fig = plt.figure()                      # 创建一个没有 axes 的 figure
    fig.suptitle('vibration curve')   # 添加标题以便我们辨别
    plt.xlim(0, tf)
    plt.xlabel("time/s")
    plt.ylabel("displacement/mm")
    plt.grid(ls='--')
    plt.plot(t,x) 
    plt.show()

def VTB2(m,c,k,x0,v0,tf,w,f0):
    wn=sqrt(k/m)
    xi=c/2/m/wn
    wd=wn*sqrt(1-xi**2)
    lam=w/wn
    A=sqrt(((v0+xi*wn*x0)**2+(x0+wd)**2)/wd**2)
    phi=atan2(x0*wd,v0+xi*wn*x0)
    B=f0/k/sqrt((1-lam**2)**2+(2*lam*xi)**2)
    Phi=atan2(2*xi*lam,1-lam**2)
    t=np.arange(0,tf,tf/1000)
    x=A*e**(-xi*wn*t)*np.sin(sqrt(1-xi**2)*wn*t+phi-pi/2)+B*np.sin(w*t-Phi-pi/2)
    fig = plt.figure()                      # 创建一个没有 axes 的 figure
    fig.suptitle('vibration curve')  # 添加标题以便我们辨别
    plt.xlim(0, tf)
    plt.xlabel("time/s")
    plt.ylabel("displacement/mm")
    plt.grid(ls='--')
    plt.plot(t,x) 
    plt.show()

def matrix_iteration1():
    '''
        This function is programmed for the example on page 43.
    '''
    n=3
    alpha=0

    M=np.array([(2,0,0),(0,1.5,0),(0,0,1)])
    K=np.array([(5,-2,0),(-2,3,-1),(0,-1,1)])
    D=np.linalg.inv(K+alpha*M)@M
    A=np.array([[1,1,1]]).T
    delta=10**(-12)
    
    
    for i in range(1,4):
        B=D@A
        p20=0
        p2=1/B[n-1]
        A=p2*B
        while abs(p2-p20)/p2>delta:
            B=D@A
            p20=p2
            p2=1/B[n-1]
            A=B*p2
        
        wn=sqrt(p20)
        f=wn/2/pi
        
        print("系统的第{}阶固有频率为{:.5}，主振型为\n{}\n".format(i,f,A))
        M1=A.T@M@A
        D=D-A@A.T@M/(p2*M1)
        
            
matrix_iteration1()
#VTB2(18.2,1.49,43.8,1,1,100,15,44.5)
    

