
import matplotlib.pyplot as plt
from math import e,sqrt,atan2,pi,sin,cos
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
    '''
    This function is programmed for example 1-3 on page 12

    The main function for this function is:
    #VTB2(18.2,1.49,43.8,1,1,100,15,44.5)
    '''
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

def matrix_iteration1(M,K,alpha=0,delta=10**(-12)):
    '''
        This function is programmed for the example on page 43.
        M is the mass matrix
        K is the stiffness matrix
        alpha is a non-negative number to prevent the stiffness matrix K from being a non-positive definite matrix
        delta controls the accuracy.

        The main function for this function is:
        #M=np.array([(2,0,0),(0,1.5,0),(0,0,1)])     #mass matrix
        #K=np.array([(5,-2,0),(-2,3,-1),(0,-1,1)])   #stiffness matrix      
        #     
        # n=13
        # X=np.identity(n)
        # Y=np.array([[1,1,1,1,1,1,1,1,1,1,1,1,1]]).T
        # M=np.array([[35.52,19.92,19.92,19.92,19.92,19.92,19.92,19.92,19.92,6.35,312.87,10.77,157.05]])
        # M=Y@M*X

        # matrix_iteration1(M,K)
    '''
    n=M.shape[0]                                    #The DOF of the system
    
    D=np.linalg.inv(K+alpha*M)@M    # D is dynamic matrix
    A=np.array([[1,1,1]]).T    #This is A1
    delta=10**(-12)
    
    for i in range(1,n+1): 
        B=D@A
        p20=0
        p2=1/B[n-1]
        A=p2*B
        while abs(p2-p20)/p2>delta:#begin to iterate
            B=D@A
            p20=p2
            p2=1/B[n-1]
            A=p2*B
        
        wn=sqrt(p2)
        f=wn/2/pi
        print("系统的第 {} 阶固有频率为 {:.5} Hz，主振型为\n{}\n".format(i,f,A))
        
        M1=A.T@M@A
        D=D-A@A.T@M/(p2*M1)#clear data of the last natural frequency from dynamic matrix

def transfer_matrix(): v
    '''
    This function is programmed for example 2-5 on page 51
    The main function for this function is:
    
    transfer_matrix() 
    '''
    J1=1
    J2=1
    J3=2
    K2=1100000
    K3=1200000
    K4=100000
    
    i=1
    theta_T_R4_show=[]
    for WN in np.arange(0,2000,0.01):
        U1=np.array([[1,0],[-WN**2*J1,1]])
        U2=np.array([[1,1/K2],[-WN**2*J2,1-WN**2*J2/K2]])
        U3=np.array([[1,1/K3],[-WN**2*J3,1-WN**2*J3/K3]])
        U4=np.array([[1,1/K4],[0,1]])#These are transfer matrixs.
        
        T1L=0
        theta1R=1
        T1R=-WN**2*J1
        theta4R_real=0#These are boundary conditions

        theta_T_L1=np.array([[theta1R,T1L]]).T
        theta_T_R1=U1@theta_T_L1 #Theta_T_1 is {\theta,T}.T_R1
        theta_T_R1[1][0]=T1R
        theta_T_R2=U2@theta_T_R1
        theta_T_R3=U3@theta_T_R2
        theta_T_R4=U4@theta_T_R3 #These sentences are to transfer.

        if abs(theta_T_R4[0]-theta4R_real)<0.0002:
            theta=np.array([[theta_T_R1[0],theta_T_R2[0],theta_T_R3[0],theta_T_R4[0]]]).T
            print("系统的第 {} 阶固有频率为 {} ,其振型为\n {}\n".format(i,WN,theta))
            i=i+1
        theta_T_R4_show.append(theta_T_R4[0][0])
    
    theta_T_R4_show=np.array(theta_T_R4_show)
    WN=np.arange(0,2000,0.01)
    fig = plt.figure()                     
    fig.suptitle('vibration curve')   # Add a title
    plt.xlim(0, 2000)
    plt.ylim(-200, 100)
    plt.xlabel("time/s") #Add a title for X axis
    plt.ylabel("displacement/mm")#Add a title for Y axis
    plt.grid(ls='--')
    plt.plot(WN,theta_T_R4_show) 
    plt.show()  

def VTB5(tf,delta):
    '''
    This function is programmed for example 3-2 on page 63.

    tf is the time to simulate;
    delta controls the accuracy.

    The main function for this function is:
    #VTB5(100,0.1)
    '''
    
    M=np.array([[2,0,0],[0,2,0],[0,0,2]]) #mass matrix
    C=np.array([[3,-1.5,0],[-1.5,3,-1.5],[0,-1.5,3]]) #damping matirx
    K=np.array([[100,-50,0],[-50,100,-50],[0,-50,100]]) #stiffness matrix
    def Force(t):
        F=np.array([[2.0*sin(3.754*t),-2.0*cos(2.2*t),sin(2.8*t)]]).T
        return F
    
    V0=np.array([[1,1,1]]).T  #V0 is the initial speed
    X0=np.array([[1,1,1]]).T  #X0 is the initial displacment
    A0=np.linalg.inv(M)@(Force(0)-C@V0-K@X0)#A0 is the initial accelaration
    MD=np.linalg.inv(M+1/2*delta*C+1/6*delta**2*K)
    X_record=X0
    t_record=np.arange(0,tf,delta)
    T=t_record.tolist()

    #####################

    #This is direct solution 1

    for t in T:
        F=Force(t+delta)
        A=MD@(F-C@(V0+delta/2*A0)-K@(X0+delta*V0+delta**2/3*A0))
        V=V0+delta/2*(A0+A)
        X=X0+delta*V0+delta**2/3*A0+delta**2/6*A
        X0=X
        V0=V
        A0=A
        X_record=np.concatenate((X_record,X0),axis=1)

    ######################

    #This is direct solution2

    # for t in T:
    #     F=Force(t)
    #     X=MD@(M@(X0+delta*V0+1/3*delta**2*A0)+C@(1/2*delta*X0+1/3*delta**2*V0+1/12*delta**3*A0)+1/6*delta**2*F)
    #     A=6/delta**2*(X-(X0+delta*V0+1/3*delta**2*A0))
    #     V=V0+delta/2*(A0+A)
    #     X0=X
    #     V0=V
    #     A0=A
    #     X_record=np.concatenate((X_record,X0),axis=1)

    ########################

    #This is New-Mark beta method

    # beta=1/6
    # for t in T:
    #     F=Force(t+delta)
    #     A=MD@(F-C@(V0+delta/2*A0)-K@(X0+delta*V0+(1/2-beta)*delta**2*A0))
    #     V=V0+delta/2*(A0+A)
    #     X=X0+delta*V0+1/2*delta**2*A0+beta*delta**3*(A-A0)/delta
    #     X0=X
    #     V0=V
    #     A0=A
    #     X_record=np.concatenate((X_record,X0),axis=1)

    ############################

    #This is Wilson theta method for direct-solution2

    # theta=1.4
    # MD=np.linalg.inv(K+(3*C)/(theta*delta)+6*M/(theta*delta)**2)
    # for t in T:
    #     F=Force(t)
    #     X_theta=MD@(M@(2*A0+6*V0/(theta*delta)+6*X0/(theta*delta)**2)+C@(1/2*delta*theta*A0+2*V0+3*X0/(theta*delta))+F)
    #     A_theta=6/(theta*delta)**2*(X_theta-X0)-6*V0/(theta*delta)-2*A0
        
    #     A=(1-1/theta)*A0+1/theta*A_theta
    #     V=V0+delta/2*(A0+A)
    #     X=X0+delta*V0+delta**2/6*(2*A0+A)
    #     X0=X
    #     V0=V
    #     A0=A
    #     X_record=np.concatenate((X_record,X0),axis=1)

    ################################

    t_record=np.append(t_record,tf)
    fig = plt.figure()
    fig.suptitle('vibration curve')   # Add a title
    plt.xlim(0, tf+1)
    plt.ylim(-1, 1)
    plt.xlabel("time/s") #Add a title for X axis
    plt.ylabel("displacement/mm")#Add a title for Y axis
    plt.grid(ls='--')
    k=2
    plt.plot(t_record,X_record[k])
    plt.show()

def VTB8(tf,delta):
    '''
    This function is programmed for R-K method on page70
    '''
    M=np.array([[2,0,0],[0,2,0],[0,0,2]]) #mass matrix
    C=np.array([[3,-1.5,0],[-1.5,3,-1.5],[0,-1.5,3]]) #damping matirx
    K=np.array([[100,-50,0],[-50,100,-50],[0,-50,100]]) #stiffness matrix
    MT=M.T

    def Force(t):
        F=np.array([[2.0*sin(3.754*t),-2.0*cos(2.2*t),sin(2.8*t)]]).T
        return F
    
    Z0=V0=np.array([[1,1,1]]).T  #V0 is the initial speed
    X0=np.array([[1,1,1]]).T  #X0 is the initial displacment
    
    X_record=X0
    t_record=np.arange(0,tf,delta)
    T=t_record.tolist()
    
    for t in T:
        K1=MT@(Force(t)-K@X0-C@V0)
        L1=Z0
        
        K2=MT@(Force(t+delta/2)-K@(X0+delta/2*L1)-C@(Z0+delta/2*K1))
        L2=Z0+delta/2*K1

        K3=MT@(Force(t+delta/2)-K@(X0+delta/2*L2)-C@(Z0+delta/2*K2))
        L3=Z0+delta/2*K2

        K4=MT@(Force(t+delta)-K@(X0+delta*L3)-C@(Z0+delta*K3))
        L4=Z0+delta*K3

        Z=Z0+delta/6*(K1+2*K2+2*K3+K4)
        X=X0+delta/6*(L1+2*L2+2*L3+L4)

        Z0=Z
        X0=X

        X_record=np.concatenate((X_record,X0),axis=1)
    
    t_record=np.append(t_record,tf)

    fig = plt.figure()
    fig.suptitle('vibration curve')   # Add a title
    plt.xlim(0, tf+1)
    plt.ylim(-1, 1)
    plt.xlabel("time/s") #Add a title for X axis
    plt.ylabel("displacement/mm")#Add a title for Y axis
    plt.grid(ls='--')
    k=1
    plt.plot(t_record,X_record[k])
    plt.show()

def problem_2_5(M,K,alpha=0,delta=10**(-12)):
    '''
        This function is programmed for the problem2-5 on page 54.
        M is the mass matrix
        K is the stiffness matrix
        alpha is a non-negative number to prevent the stiffness matrix K from being a non-positive definite matrix
        delta controls the accuracy.

        The main function for this function is:
        # m=1
        # k=1

        # M=np.array([(m,0,0),(0,m,0),(0,0,m)])     #mass matrix
        # K=np.array([(k,-k,0),(-k,2*k,-k),(0,-k,k)])   #stiffness matrix      

        # problem_2_5(M,K,alpha=1)
    '''
    n=M.shape[0]                    #The DOF of the system    
    D=np.linalg.inv(K+alpha*M)@M    # D is dynamic matrix
    A=np.array([[1,1,1]]).T         #This is A1
    
    for i in range(1,n+1): 
        
        B=D@A
        p20=0 
        p2=1/B[n-1]
        A=p2*B
        while abs(p2-p20)/p2>delta:#begin to iterate
            B=D@A
            p20=p2
            p2=1/B[0]
            A=p2*B

        wn=sqrt(p2-alpha)
        f=wn/2/pi
        print("系统的第 {} 阶固有频率为 {:.5} Hz，主振型为\n{}\n".format(i,f,A))
        
        M1=A.T@M@A
        D=D-A@A.T@M/(p2*M1)#clear data of the last natural frequency from dynamic matrix

def problem_2_14():
    '''
    This function is programmed for problem 2-14 on page.57
    '''
    K0=1
    K1=K0
    K2=K0
    J=1
    J1=J
    J2=J/2

    i=1
    theta_T_R0_show=[]

    WNf=4
    delta=0.001
    for WN in np.arange(0,WNf,delta):
        U1=np.array([[1-WN**2*J1/K1,1/K1],[-WN**2*J1,1]])
        U2=np.array([[1-WN**2*J2/K2,1/K2],[-WN**2*J2,1]])#These are transfer matrixs.

        theta0R=0
        theta2R=1
        T2R=0#These are boundary conditions

        

        theta_T_R2=np.array([[theta2R,T2R]]).T
        theta_T_R1=U2@theta_T_R2 #Theta_T_1 is {\theta,T}.T_R1
        theta_T_R0=U1@theta_T_R1#These sentences are to transfer.

        if abs(theta_T_R0[0][0]-theta0R)<0.002:
            theta=np.array([[theta_T_R1[0,0],theta_T_R2[0,0]]]).T
            print("系统的第 {} 阶固有频率为 {} ,其振型为\n {}\n".format(i,WN,theta))
            i=i+1
        theta_T_R0_show.append(theta_T_R0[0][0])
    
    theta_T_R0_show=np.array(theta_T_R0_show)
    WN=np.arange(0,WNf,delta)
    fig = plt.figure()                     
    fig.suptitle('vibration curve')   # Add a title
    plt.xlim(0, 4)
    plt.ylim(-1, 1)
    plt.xlabel("time/s") #Add a title for X axis
    plt.ylabel("displacement/mm")#Add a title for Y axis
    plt.grid(ls='--')
    plt.plot(WN,theta_T_R0_show) 
    plt.show()

def problem_3_2(tf,delta):
    '''
    This function is programmed for problem 3-2 on page 73.
    It applies NewMark-beta method.
    Its main function is:
    # problem_3_2(4,0.002)
    '''
    
    M=np.array([[2,0],[0,1]]) #mass matrix
    C=np.array([[250,-250],[-250,250]]) #damping matirx
    K=np.array([[1000,-500],[-500,500]]) #stiffness matrix
    def Force(t):
        F=np.array([[2.0*sin(70*t),0]]).T
        return F
    
    V0=np.array([[0,0]]).T  #V0 is the initial speed
    X0=np.array([[0.2,0]]).T  #X0 is the initial displacment
    A0=np.linalg.inv(M)@(Force(0)-C@V0-K@X0)#A0 is the initial accelaration
    MD=np.linalg.inv(M+1/2*delta*C+1/6*delta**2*K)
    X_record=X0
    t_record=np.arange(0,tf,delta)
    T=t_record.tolist()

    beta=1/6
    for t in T:
        F=Force(t+delta)
        A=MD@(F-C@(V0+delta/2*A0)-K@(X0+delta*V0+(1/2-beta)*delta**2*A0))
        V=V0+delta/2*(A0+A)
        X=X0+delta*V0+1/2*delta**2*A0+beta*delta**3*(A-A0)/delta
        X0=X
        V0=V
        A0=A
        X_record=np.concatenate((X_record,X0),axis=1)

    t_record=np.append(t_record,tf)
    fig = plt.figure()
    fig.suptitle('vibration curve of m1')   # Add a title
    plt.xlim(0, tf+1)
    plt.ylim(-0.5, 0.5)
    plt.xlabel("time/s") #Add a title for X axis
    plt.ylabel("displacement/mm")#Add a title for Y axis
    plt.grid(ls='--')
    k=0
    plt.plot(t_record,X_record[k])
    plt.show()

def problem_3_3(tf,delta):
    '''
    This function is programmed for problem 3-2 on page 73.
    It applies Wilson-theta method.
    Its main function is:
    # problem_3_3(40,0.002)
    '''
    
    M=np.array([[3,0,0],[0,1,0],[0,0,2]]) #mass matrix
    C=np.array([[0,0,0],[0,0,0],[0,0,0]]) #damping matirx  C=np.array([[0,0,0],[0,0,0],[0,0,0]]) 
    K=np.array([[140,-60,-20],[-60,200,-60],[-20,-60,140]]) #stiffness matrix
    def Force(t):
        F=np.array([[1,0,0]]).T
        return F
    
    V0=np.array([[1,1,1]]).T  #V0 is the initial speed
    X0=np.array([[1,1,1]]).T  #X0 is the initial displacment
    A0=np.linalg.inv(M)@(Force(0)-C@V0-K@X0)#A0 is the initial accelaration
    MD=np.linalg.inv(M+1/2*delta*C+1/6*delta**2*K)
    X_record=X0
    t_record=np.arange(0,tf,delta)
    T=t_record.tolist()

    theta=1.4
    MD=np.linalg.inv(K+(3*C)/(theta*delta)+6*M/(theta*delta)**2)
    for t in T:
        F=Force(t)
        X_theta=MD@(M@(2*A0+6*V0/(theta*delta)+6*X0/(theta*delta)**2)+C@(1/2*delta*theta*A0+2*V0+3*X0/(theta*delta))+F)
        A_theta=6/(theta*delta)**2*(X_theta-X0)-6*V0/(theta*delta)-2*A0
        
        A=(1-1/theta)*A0+1/theta*A_theta
        V=V0+delta/2*(A0+A)
        X=X0+delta*V0+delta**2/6*(2*A0+A)
        X0=X
        V0=V
        A0=A
        X_record=np.concatenate((X_record,X0),axis=1)

    t_record=np.append(t_record,tf)
    fig = plt.figure()
    fig.suptitle('vibration curve of m3')   # Add a title
    plt.xlim(0, tf+1)
    plt.ylim(-5, 5)
    plt.xlabel("time/s") #Add a title for X axis
    plt.ylabel("displacement/mm")#Add a title for Y axis
    plt.grid(ls='--')
    k=2
    plt.plot(t_record,X_record[k])
    plt.show()

def example_6_1():
    '''
    This function is programmed for example 6-1 on page 120 which applies thransfer matrix method.
    The main function for this function is:
    
    example_6_1()
    '''
    J0=0.155
    J1=0.0079
    J2=0.0065
    J3=0.024
    J4=0.0036
    J5=0.0088
    J6=0.0084
    J7=0.017
    K1=1/130*10**6
    K2=1/190*10**6
    K3=1/50*10**6
    K4=1/39*10**6
    K5=1/280*10**6
    K6=1/300*10**6
    K7=1/39*10**6

    i=1
    theta_T_R7_show=[]

    WN_range=np.arange(0,2000,0.1)

    for WN in WN_range:
        U0=np.array([[1,0],[-WN**2*J0,1]])
        U1=np.array([[1,1/K1],[-WN**2*J1,1-WN**2*J1/K1]])
        U2=np.array([[1,1/K2],[-WN**2*J2,1-WN**2*J2/K2]])
        U3=np.array([[1,1/K3],[-WN**2*J3,1-WN**2*J3/K3]])
        U4=np.array([[1,1/K4],[-WN**2*J4,1-WN**2*J4/K4]])#These are transfer matrixs.
        U5=np.array([[1,1/K5],[-WN**2*J5,1-WN**2*J5/K5]])
        U6=np.array([[1,1/K6],[-WN**2*J6,1-WN**2*J6/K6]])
        U7=np.array([[1,1/K7],[-WN**2*J7,1-WN**2*J4/K7]])
        
        T0L=0
        theta0R=1
        theta7R_real=1#These are boundary conditions

        theta_T_L0=np.array([[theta0R,T0L]]).T
        theta_T_R0=U0@theta_T_L0 #Theta_T_1 is {\theta,T}.T_R1
        theta_T_R1=U1@theta_T_R0
        theta_T_R2=U2@theta_T_R1
        theta_T_R3=U3@theta_T_R2
        theta_T_R4=U4@theta_T_R3 #These sentences are to transfer.
        theta_T_R5=U5@theta_T_R4
        theta_T_R6=U6@theta_T_R5
        theta_T_R7=U7@theta_T_R6

        if abs(theta_T_R7[0]-theta7R_real)<0.000001:
            theta=np.array([[theta_T_R1[0],theta_T_R2[0],theta_T_R3[0],theta_T_R4[0],theta_T_R5[0],theta_T_R6[0],theta_T_R7[0]]]).T
            print("系统的第 {} 阶固有频率为 {} ,其振型为\n {}\n".format(i,WN,theta))
            i=i+1
        theta_T_R7_show.append(theta_T_R7[0][0])
    
    theta_T_R7_show=np.array(theta_T_R7_show)
    fig = plt.figure()                     
    fig.suptitle('vibration curve')   # Add a title
    plt.xlim(0, 2000)
    plt.ylim(-200, 1500)
    plt.xlabel("Wn/(1/s)") #Add a title for X axis
    plt.ylabel("displacement/mm")#Add a title for Y axis
    plt.grid(ls='--')
    plt.plot(WN_range,theta_T_R7_show) 
    plt.show()

def example_6_1_2(M,K,alpha=0,delta=10**(-12)):
    '''
    This function is programmed for example6-1on page 120 whilch applies matrix iteration method.

    # M=np.array([[0.155,0.0079,0.0065,0.024,0.0036,0.0088,0.0084,0.017]]).T
    # M=M@np.ones(8).reshape(1,8)
    # K1=1/130*10**6
    # K2=1/190*10**6
    # K3=1/50*10**6
    # K4=1/39*10**6
    # K5=1/280*10**6
    # K6=1/300*10**6
    # K7=1/39*10**6
    # K=np.array([[K1,-K1,0,0,0,0,0,0],[-K1,K1+K2,-K2,0,0,0,0,0],[0,-K2,K2+K3,-K3,0,0,0,0],[0,0,-K3,K3+K4,-K4,0,0,0],[0,0,0,-K4,K4+K5,-K5,0,0],[0,0,0,0,-K5,K5+K6,-K6,0],[0,0,0,0,0,-K6,K6+K7,-K7],[0,0,0,0,0,0,-K7,K7]])
    # example_6_1_2(M,K,alpha=10)
    '''
    n=M.shape[0]                                    #The DOF of the system
    
    D=np.linalg.inv(K+alpha*M)@M    # D is dynamic matrix
    A=np.ones(n).reshape(1,n).T    #This is A1
    
    for i in range(1,n+1): 
        B=D@A
        p20=0
        p2=1/B[n-1]
        A=p2*B
        while abs(p2-p20)/p2>delta:#begin to iterate
            B=D@A
            p20=p2
            p2=1/B[n-1]
            A=p2*B
        
        wn=sqrt(p2)
        f=wn/2/pi
        print("系统的第 {} 阶固有频率为 {:.5} Hz，主振型为\n{}\n".format(i,f,A))
        
        M1=A.T@M@A
        D=D-A@A.T@M/(p2*M1)#clear data of the last natural frequency from dynamic matrix

