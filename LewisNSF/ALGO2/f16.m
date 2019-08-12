% Simulation of f-16
function xdot=f16(t,x,options,K,d,A,B,C,D);
% F16 dynamics
r=1;
G=[0 0 0 0 1]';

Ac=A-B*K*C;

F=[0 0 1 0]';

xdot=Ac*x+D*d+(G-B*K*F)*r;
