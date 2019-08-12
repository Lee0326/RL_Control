function xdot=subreintwojob(t,x)
global P R A B C D Q rowA;
global u;
global k
global T
global L
global K
  tm=t+(k-1)*T;
% u=0.02*sin(0.1*t)+0.03*sin(0.2*t)+30*sin(20*t);
  %calculating the control signal
%  u=-K*C*x(1:rowA);
 u=-K*C*x(rowA+2:rowA*2+1);
% u=-K*C*x(1:2);
% u=3*sin(2*tm);
% u=-K*C*x(1:rowA)+0.02*sin(1.0*tm)^5+0.02*cos(20*tm)^3+0.02*sin(0.1*tm)*cos(5*tm);
%  u=-inv(R)*B'*P*x(1:rowA)+0.02*sin(1.0*tm)^5;%+2*cos(20*tm)^3+2*sin(0.1*tm)*cos(5*tm);
  %updating the derivative of the state=[x(1:2) V]
 
% xdot=[A*x(1:rowA)+B*u  
%     x(1:rowA)'*Q*x(1:rowA)+u'*R*u
%     [-1.0051 1;-1.3851 0]*x(rowA+2:rowA+1+rowA)+B*u+[7.9950 ;18.6149]*C*(x(1:rowA)-x(rowA+2:rowA+rowA+1))
%      x(rowA+3:rowA+rowA+2)'*Q*x(rowA+3:rowA+rowA+2)+u'*R*u];
% xdot=[[-1.0051 1;-1.3851 0]*x(1:rowA)+B*u-[7.9950;18.6149]*C*(x(1:rowA)-x(rowA+2:rowA*2+1))%Khut=10,20
%        x(1:rowA)'*Q*x(1:rowA)+u'*R*u
%      A*x(rowA+2:rowA*2+1)+B*u  
%      x(rowA+2:rowA+rowA+1)'*Q*x(rowA+2:rowA+rowA+1)+u'*R*u];
%  xdot=[[ -1.0201  1.0000; -1.3847   0]*x(1:rowA)+B*u-[198.9799; 198.6153]*C*(x(1:rowA)-x(rowA+2:rowA*2+1))%Khut=200,200,C=[1 2]
xdot=[[-1.0652    1.0000;   -1.2171         0]*x(1:rowA)+B*u-[ 198.9348;198.7829]*C*(x(1:rowA)-x(rowA+2:rowA*2+1))%Khut=200,200
% xdot=[[-1.0652    1.0000;   -1.2171         0]*x(1:rowA)+B*u-[ 198.9348;198.7829]*C*(x(1:rowA)-x(rowA+2:rowA*2+1))%Khut=200,200,C=[1 0]
%        x(1:rowA)'*Q*x(1:rowA)+u'*R*u
 x(rowA+2:rowA*2+1)'*Q*x(rowA+2:rowA*2+1)+u'*R*u
     A*x(rowA+2:rowA*2+1)+B*u  
     x(rowA+2:rowA*2+1)'*Q*x(rowA+2:rowA*2+1)+u'*R*u];