clc
close all
clear all
global A B  C F  r Q R K
%******************************************************
%matrix coefficient
A=[-1.01887 0.90506   -0.00215   0    0;...
    0.82225  -1.07741 -0.17555   0    0;...
    0        0        -20.2      0    0;...
    10       0         0         -10  0;...
    -16.26   -0.9788   0.4852    0   0];
B=[0 0 20.2 0 0]';
G=[0 0 0 0 1]';
C=[0 0 0 57.2958  0;...
   0 57.2958 0 0 0;...
   -16.26 -0.9788 0.4852 0 0;...
   0 0 0 0 1];
F=[0 0 0 1 0];
H1=[16.26 0.9788 -0.04852 0 0];
D=0;
%******************************************************
%******************************************************
%estimate the size of matrix
[rowA columnA]=size(A);
[rowB columnB]=size(B);
[rowC columnC]=size(C);
[rowD columnD]=size(D);
[rowF columnF]=size(F);
[rowG columnG]=size(G);
[rowr columnr]=size(r);
%******************************************************
%******************************************************
Statesnum=columnA;
Inputsnum=columnB;
Outputsnum=rowC;
%******************************************************
%******************************************************
a1=min(rowC,columnC);
if rank(C)<a1
    sprintf('The rank of C is not full row rank',2,3)
    return;
end
%******************************************************

%******************************************************
%Q,R,X
Q=[264 16 1 0 0;16 60 0 0 0;1 0 0 0 0;0 0 0 0 0;0 0 0 0 100];
rou=0.1;
R=rou*eye(Inputsnum);
X=eye(Statesnum);
%******************************************************

%******************************************************
L0=0.001*[1 1 1 1 1];%gave an initial value of L closed to 0
P0=are(A,B*inv(R)*B',Q+L0'*inv(R)*L0);%get the initial P0
K0=inv(R)*(B'*P0+L0)*C'*inv(C*C')%get the initial K0
L0=R*K0*C-B'*P0;%iteration L0
%***********************************
% Check if K0 stable, if not return
damp(A-B*K0*C);
    Eig1=eig(A-B*K0*C);%calclate the eig of A-B*K0*C
    REig1=real(Eig1);%get the real part
    for i=1:Statesnum
        if REig1(i)>=0
           return;
        end
    end
%***********************************
    
for i=1:1000
    J0=abs(trace(P0*X));%%calculate the performance value
    P1=are(A,B*inv(R)*B',Q+L0'*inv(R)*L0);%calculate P for the new L0
    K1=inv(R)*(B'*P1+L0)*C'*inv(C*C');%calculate K for the new P
    L1=R*K1*C-B'*P1;%calculate the L for the new P,K
    J1=abs(trace(P1*X));%calculate the performance value
    %
    Eig2=eig(A-B*K1*C);%calclate the eig of A-B*K0*C
    REig2=real(Eig2);%get the real part
    if abs(K1-K0)<=0.0001
          P=P1;
          K=K1;
          L=L1;
          break
       else
         P0=P1;
         K0=K1;
         L0=L1;
    end
       
%     if REig2<0
%         if (J1<=J0)%if J1<=J0 then iterate the value of P,K,L
%            P0=P1;
%            K0=K1;
%            L0=L1;
%         else
%             break
%         end
%     end
    
end
P
L
K
x0=[1 0 0 0 0];%the initial value
T=10; 
tspan=[0 T];
[t,Y]=ode45(@ABCD,tspan,x0);
YY=C*Y';
%figure
figure(1)
plot(t,Y(:,1),'k','LineWidth',2)
xlabel('time response')
ylabel('Angle of attack')
grid on

figure(2)
plot(t,Y(:,2),'k','LineWidth',2);
xlabel('time response')
ylabel('Pitch Rate')
grid on

figure(3)
plot(t,Y(:,3),'k','LineWidth',2);
xlabel('time response')
ylabel('elevator actuator')
grid on

figure(4)
plot(t,Y(:,4),'k','LineWidth',2);
xlabel('time response')
ylabel('filtered measurement of angle of attack')
grid on

figure(5)
plot(t,Y(:,4),'k','LineWidth',2);
xlabel('time response')
ylabel('integral controller')
grid on

figure(6)
plot(t,YY(3,:),'k','LineWidth',2);
xlabel('time response')
ylabel('tracking error')
grid on

figure(7)
poles=eig(A-B*K*C)
plot(poles,'.')
xlabel('real axis')
ylabel('imaginery axis')
title('poles map')
grid on






