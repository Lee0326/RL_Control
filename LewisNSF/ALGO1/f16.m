% Simulation of f-16
function xdot=f16(t,x,options,K,d,A,B,C,D);
% F16 dynamics
Ac=A-B*K*C;
xdot=Ac*x+D*d;
