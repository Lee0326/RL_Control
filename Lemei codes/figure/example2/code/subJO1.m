function xdot=subJO1(t,x)
global P R A B C D Q rowA;
global u;
global k
global T
global L
global K
 
A1=[0  1   0 0;...
   0 -10  0 0 ;...
  -1  0   0 0;...
   0  200 0 -200];
L=[10 0 0 ;-1 1 0 ; 1 10 -1;1 10 1];
B=[0; 1; 0; 0];
  %calculating the control signal
%  u=-K*C*x(rowA+2:rowA*2+1);
 u=-K*C*x(1:rowA);
% u=-K*C*x(1:rowA)+0.02*sin(1.0*tm)^5+0.02*cos(20*tm)^3+0.02*sin(0.1*tm)*cos(5*tm);
%  u=-inv(R)*B'*P*x(1:rowA)+0.02*sin(1.0*tm)^5;%+2*cos(20*tm)^3+2*sin(0.1*tm)*cos(5*tm);
  %updating the derivative of the state=[x(1:2) V]
  Whut=[ 0.0525 -0.0087 -0.0510 0.1177 -0.0206; 0.0634 -0.0151 -0.0263 0.0673 0.0130; 0.0319 0.0724 0.1477 0.0853 0.1802; 0.0108 0.0889 -0.0068 0.1030 0.0967];
%   xdot=[A1*x(1:rowA)+Whut*tanh([x(1:rowA);u])+B*u+L*(C*x(rowA+2:rowA+1+rowA)-C*x(1:rowA))
%       x(rowA+2:rowA*2+1)'*Q*x(rowA+2:rowA*2+1)+u'*R*u
%       A*x(rowA+2:rowA+rowA+1)+B*u];
   xdot=[A*x(1:rowA)+B*u
        x(1:rowA)'*Q*x(1:rowA)+u'*R*u
      A1*x(rowA+2:rowA+rowA+1)+B*u+L*(C*x(1:rowA)-C*x(rowA+2:rowA+rowA+1))];