% Simulation of f-16
function xdot=f16(t,x,options,F,d,A,B,C,D);
% F16 dynamics
Ac=A-B*F*C;
xdot=Ac*x+D*d;
