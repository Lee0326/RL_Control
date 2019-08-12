
%××××××××××××××××××××××××××××××××××××××××××××××
function xdot=cal(t,x)

global A B G C F D H1 r Q R K 

%calculating the control signal

u=-K*(C*x);

  %updating the derivative of the state=[x(1:2) V]
% A
% x
% B
% u
% G
% r
xdot=[A*x+B*u];
%××××××××××××××××××××××××××××××××××××××××××××××××?