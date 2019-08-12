function xdot=slovelinear(t,x)
global k n m l C

%The matrix of the system
A=[0  1   0 0;...
   0 -11  0 0 ;...
  -1  0   0 0;...
   0  1/0.005 0 -1/0.005];
A1=[0  1   0 0;...
   0 -10  0 0 ;...
  -1  0   0 0;...
   0  100 0 -100];%A2
% A1=[0  1   0  0;...
%     0  0   1  0 ;...
%     0  0   0  1;...
%    -12 -5 -10 -7];%A3
%  A1=[0  1   0  0;...
%     0  0   1  0 ;...
%     0  0   0  1;...
%    0 0 0 0];
B=[0; 1; 0; 0];
C=[0 0 1 0;...
   -1 0 0 0;...
    0 0 0 1];
L=[10 0 0 ;-1 1 0 ; 1 10 -1;1 10 1];%this Gain is good with A2
%L=[10 0 -1 ;-1 1 0 ; 1 10 -1;1 10 1];%this Gain is good with A2
%  L=[280 0 0;0 0 10;-30 0 150 ;-100 -10 200];%this gain is with A3
C1=eye(3,3);

Ac=A1-L*C;

%The learning rates 
  yita1=180;
  yita2=180;
  
%The design papameters
  xita1=150;
  xita2=150;
  
%The input
  u=10*sin(200*t)+cos(0.5*t)+3*sin(5*t)+sin(40*t)+5*sin(30*t)+40*sin(4*t)+15*sin(3*t)+100*sin(20*t)+1.3*cos(5.5*t)+31*sin(115*t)+0.9*sin(4*t)+15*sin(36*t)+45*sin(9.4*t)+10*sin(30*t);
%  u=-50;
%restore Whut
 Whut=zeros(n,k);
 for i=1:k
     Whut(1:n,i)=x((i+1)*n+1:(i+2)*n);
 end

Whutdot=-yita1*((C*x(1:n)-C*x(n+1:2*n))'*C*inv(Ac))'*tanh([x(n+1:2*n);u])'-xita1*norm(C*x(1:n)-C*x(n+1:2*n))*Whut;
    
 
  
  Whutdot0=zeros(n*k,1);
for i=1:k
    Whutdot0((i-1)*n+1:n*i)=Whutdot(1:n,1);
end

xdot=[A*x(1:n)+B*u
       A1*x(n+1:2*n)+Whut*tanh([x(n+1:2*n);u])+B*u+L*(C*x(1:n)-C*x(n+1:2*n))
      Whutdot0];

    
    
    
    
    
    
    
    
    
    