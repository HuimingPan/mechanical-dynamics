clc,clear all.close all
%输入质量矩阵
M=[2 0 0;0 1.5 0;0 0 1];
%输入刚度矩阵
K=[5 -2 0;-2 3 -1;0 -1 1];
%计算特征值与特征向量
D=inv(K+0.1*M)*M;
A=ones(3,1);
for i=1:3
    pp0=0;
    i
    B=D*A;
    pp=1/B(3);
    A=B/B(3);
    while abs((pp-pp0)/pp)>1e-10
        pp0=pp;
        B=D*A;
        pp=1/B(3);
        A=B*pp;
    end
       f=sqrt(pp-0.1)/2/pi;
       A=vpa(A,5);
       f=vpa(f,5);
       disp(A)
       disp(f)
    D=D-A*A'*M/(A'*M*A*pp);
end