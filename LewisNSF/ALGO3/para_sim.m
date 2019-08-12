%Program to call algorithm for calculating output feedback gains

%Algorithm by J. Gadewadikar,F. L. Lewis, and M Abu Khalaf.

%Writtern By J. Gadewadikar.

clear all
close all
clc
%----------------------------------------------------------------
%8.1.1


A=[-.3220   +0.0640  +0.0364    -0.9917   +0.0003   +0.0008   0;
    0           0       1       +0.0037       0           0   0;
    -30.6492    0    -3.6784    0.6646      -0.7333  0.1315   0;
    8.5396      0    -0.0254    -0.4764     -0.0319 -0.0620   0;
    0           0       0           0       -20.2       0     0;
    0           0       0           0           0       -20.2 0;
    0           0       0         57.2958       0       0     -1];  

B=[0    0;
   0    0;
   0    0;
   0    0;
   20.2 0;
   0    20.2;
   0    0];

C=[0        0      0       57.2958     0       0   -1;
   0        0       57.2958 0           0       0     0;
   57.2958  0       0       0           0       0      0;
   0        57.2958 0       0           0       0       0];   
    

D=[ 0   0   0   0   1   0   0;
    0   0   0   0   0   1   0;
    0   0   0   0   0   0   1]';

gamma=1.09989;
    
R=1*[1 0;
    0 1]; % control weighing matrix

display('Closed loop poles')
[P,K,F,L,i] = closed_form_kv2(A,B,C,R,D,gamma);
lewis1_poles=eig(A-B*F*C) % closed loop poles
i
d=[1 1 1]';
ti=0;
tf=10;
tspan=[ti tf];
x0=[1 0 0 0 0 0 0]'; % initial state vector
[t,x]= ode23('f16_1',tspan,x0,[],F,d,A,B,C,D);
figure
plot(t(:,1),x(:,1),t(:,1),x(:,4),'LineWidth',2)
legend('Beta','r');
xlabel('Time')
ylabel('Time response')
title('Dutch-roll states beta and r')
grid on
figure
plot(t(:,1),x(:,2),t(:,1),x(:,3),'LineWidth',2)
legend('Phi','p');
xlabel('Time')
ylabel('Time response')
title('Roll mode states phi and p')
grid on


