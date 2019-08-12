function xdot=subF16(t,x)
global P R A B C D Q rowA;
global u;
global k
global T
global L
global K
 tm=t+(k-1)*T;
A1=[-1        1   0   0    0;...
    1  -1 -0.2   0    0;...
    0        0        -20     0    0;...
    10       0         0         -10  0;...
    -16   -1   0.5    0   0];%A3

L=[10 0 0 0 ;-1 1 0 0; 1 10 -1 0;0 1 10 1;0 1 10 1];%this Gain is good with A2
Ac=A1-L*C;
  %calculating the control signal
%  u=-K*C*x(rowA+2:rowA*2+1);
 u=-K*C*x(1:rowA);
% u=-K*C*x(1:rowA)+0.02*sin(1.0*tm)^5+0.02*cos(20*tm)^3+0.02*sin(0.1*tm)*cos(5*tm);
%  u=-inv(R)*B'*P*x(1:rowA)+0.02*sin(1.0*tm)^5;%+2*cos(20*tm)^3+2*sin(0.1*tm)*cos(5*tm);
  %updating the derivative of the state=[x(1:2) V]
%   Whut=[ 0.0015    0.0141   -0.0001   -0.0007    0.0007  
%            -0.0028    0.0236    0.0152     0.0037   -0.0160
%            0.0018    0.0240    0.0063   -0.0025   -0.0111
%         -0.0048   0.0164   -0.0023   -0.0045   -0.0111
%        -0.0104    0.0240    0.0132    0.0084  -0.0002]';
Whut=[0.0014    0.0141   -0.0001   -0.0004    0.0002   
           0.0029    0.0198   -0.0002   -0.0119   -0.0074 
           0.0087    0.0209   -0.0133   -0.0112   -0.0077  
           -0.0043   0.0319    0.0009   -0.0078   -0.0084  
           -0.0022    0.0201   -0.0018   -0.0190   -0.0087];
%   xdot=[A1*x(1:rowA)+Whut*tanh([x(1:rowA)])+B*u+L*(C*x(rowA+2:rowA+1+rowA)-C*x(1:rowA))
%   x(rowA+2:rowA*2+1)'*Q*x(rowA+2:rowA*2+1)+u'*R*u
%       A*x(rowA+2:rowA+rowA+1)+B*u];
   xdot=[A*x(1:rowA)+B*u
        x(1:rowA)'*Q*x(1:rowA)+u'*R*u
        A1*x(rowA+2:rowA+rowA+1)+Whut*tanh(x(rowA+2:rowA+rowA+1))+B*u+L*(C*x(1:rowA)-C*x(rowA+2:rowA+rowA+1))];