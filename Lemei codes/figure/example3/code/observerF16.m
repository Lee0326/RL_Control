function xdot=observerF16(t,x)
global A B C k n m l xxx number
number=number+1
%The matrix of the observer system
% A1=[0  1   0   0   0;...
%     0  0   1   0   0;...
%     0  0   0   1   0;...
%     0  0   0   0   1;...
%   -12 -7  -16 -8  -10];%A2
A1=[-1        1   0   0    0;...
    1  -1 -0.2   0    0;...
    0        0        -20     0    0;...
    10       0         0         -10  0;...
    -16   -1   0.5    0   0];%A3
% A1=[-1.01887 0.90506   -0.00215   0    0;...
%     0.82225  -1.07741 -0.17555   0    0;...
%     0        0        -20.2      0    0;...
%     10       0         0         -10  0;...
%     -16.26   -0.9788   0.4852    0   0];
L=[10 0 0 0 ;-1 1 0 0; 1 10 -1 0;0 1 10 1;0 1 10 1];%this Gain is good with A2
Ac=A1-L*C;
%  damp(Ac)
%The learning rates 
  yita1=1.8;
  yita2=1.8;
  
%The design papameters
  xita1=1.5;
  xita2=1.5;
  
%The input
% u=-100;
u=-[0   0    0  10]*C*x(1:n)+30*sin(20*t)+cos(0.5*t)+3*sin(5*t)^5+sin(40*t)+10*sin(200*t)*cos(5*t)+cos(0.5*t)+3*sin(5*t)^5+sin(40*t)+5*sin(30*t)+40*sin(4*t)+15*sin(3*t)+100*sin(20*t)+1.3*cos(5.5*t)+31*sin(115*t)+0.9*sin(4*t)+15*sin(36*t)+45*sin(9.4*t)+10*sin(30*t)+1*sin(1*t)+2*sin(2*t)+3*sin(3*t)^3+4*sin(4*t)+5*sin(5*t)+6*sin(6*t)*7*sin(7*t)^110+8*sin(8*t)*9*cos(9*t)+10*sin(10*t)+11*sin(11*t)+12*sin(12*t)+13*sin(13*t)+14*sin(14*t)+15*sin(15*t)+16*sin(16*t)+17*sin(17*t)+18*cos(18*t)+19*sin(19*t)+20*sin(20*t)+21*sin(21*t)+22*sin(22*t)+23*sin(23*t)++24*sin(24*t)+25*sin(25*t)++26*sin(26*t)+27*sin(27*t)++28*sin(28*t)+29*sin(29*t)+30*sin(30*t)+31*sin(31*t)+32*sin(32*t)+33*sin(330*t);
%   u=10*sin(200*t)*cos(5*t)+cos(0.5*t)+3*sin(5*t)^5+sin(40*t)+5*sin(30*t)+40*sin(4*t)+15*sin(3*t)+100*sin(20*t)+1.3*cos(5.5*t)+31*sin(115*t)+0.9*sin(4*t)+15*sin(36*t)+45*sin(9.4*t)+10*sin(30*t)+1*sin(1*t)+2*sin(2*t)+3*sin(3*t)^3+4*sin(4*t)+5*sin(5*t)+6*sin(6*t)+7*sin(7*t)+8*sin(8*t)+9*cos(9*t)+10*sin(10*t)+11*sin(11*t)+12*sin(12*t)+13*sin(13*t)+14*sin(14*t)+15*sin(15*t)+16*sin(16*t)+17*sin(17*t)+18*cos(18*t)+19*sin(19*t)+20*sin(20*t)+21*sin(21*t)+22*sin(22*t)+23*sin(23*t)++24*sin(24*t)+25*sin(25*t)++26*sin(26*t)+27*sin(27*t)++28*sin(28*t)+29*sin(29*t)+30*sin(30*t)+31*sin(31*t)+32*sin(32*t)+33*sin(330*t);
% u=1*sin(1*t)+2*sin(2*t)+3*sin(3*t)^3+4*sin(4*t)+5*sin(5*t)+6*sin(6*t)+7*sin(7*t)+8*sin(8*t)+9*cos(9*t)+10*sin(10*t)+11*sin(11*t)+12*sin(12*t)+13*sin(13*t)+14*sin(14*t)+15*sin(15*t)+16*sin(16*t)+17*sin(17*t)+18*cos(18*t)+19*sin(19*t)+20*sin(20*t)+21*sin(21*t)+22*sin(22*t)+23*sin(23*t)++24*sin(24*t)+25*sin(25*t)++26*sin(26*t)+27*sin(27*t)++28*sin(28*t)+29*sin(29*t)+30*sin(30*t)+31*sin(31*t)+32*sin(32*t)+33*sin(330*t)+34*sin(34*t)+35*sin(350*t)+360*sin(360*t)+37*sin(37*t); 
%restore Whut
  Whut=zeros(n,k);
  for i=1:k
      Whut(1:n,i)=x((i+1)*n+1:(i+2)*n);
  end
%restore Vhut
%Vhut=zeros(k,n+m);
for i=1:k
    Vhut(1:k,i)=x((k+2)*n+(i-1)*k+1:(k+2)*n+i*k);
end
  
  
%   Whutdot=-yita1*((C*x(1:n)-C*x(n+1:2*n))'*C*inv(Ac))'*tanh(Vhut*[x(n+1:2*n);u])'-xita1*norm(C*x(1:n)-C*x(n+1:2*n))*Whut;
%     Whutdot=-yita1*((C*x(1:n)-C*x(n+1:2*n))'*C*inv(Ac))'*tanh([x(n+1:2*n);u])'-xita1*norm(C*x(1:n)-C*x(n+1:2*n))*Whut;
  Whutdot=-yita1*((C*x(1:n)-C*x(n+1:2*n))'*C*inv(Ac))'*tanh(x(n+1:2*n))'-xita1*norm(C*x(1:n)-C*x(n+1:2*n))*Whut;
%   sigma=tanh(Vhut*[x(n+1:2*n);u]);
  sigma=tanh(Vhut*[x(n+1:2*n)]);
  for i=1:k
      sigmagamma(i)=sigma(i)^2;
  end
  gamma=diag(sigmagamma);

%   Vhutdot=-yita2*((C*x(1:n)-C*x(n+1:2*n))'*C*inv(Ac)*Whut*(eye(k)-gamma))'*sign([x(n+1:2*n);u])'-xita2*norm(C*x(1:n)-C*x(n+1:2*n))*Vhut;
  Vhutdot=-yita2*((C*x(1:n)-C*x(n+1:2*n))'*C*inv(Ac)*Whut*(eye(k)-gamma))'*sign([x(n+1:2*n)])'-xita2*norm(C*x(1:n)-C*x(n+1:2*n))*Vhut;
  Whutdot0=zeros(n*k,1);
  for i=1:k
      Whutdot0((i-1)*n+1:n*i)=Whutdot(1:n,1);
  end

  Vhutdot0=zeros(k*k,1);
  for i=1:k
      %Vhutdot0((i-1)*k+1:k*i)=Vhutdot(1:k,1);
      Vhutdot0((i-1)*k+1:k*i)=0;
  end

%    Whutknown=[ 0.0015    0.0141   -0.0001   -0.0007    0.0007  
%            -0.0028    0.0236    0.0152     0.0037   -0.0160
%            0.0018    0.0240    0.0063   -0.0025   -0.0111
%         -0.0048   0.0164   -0.0023   -0.0045   -0.0111
%        -0.0104    0.0240    0.0132    0.0084  -0.0002]';
% Whutknown=[0.0015    0.0144   -0.0001   -0.0006    0.0005 
%            0.0003    0.0014    0.0063   -0.0038   -0.0057 
%            0.0033    0.0083   -0.0086   -0.0110   -0.0122  
%            0.0062    0.0059   -0.0035   -0.0040   -0.0022   
%            0.0020    0.0079   -0.0096   -0.0040  -0.0076];
Whutknown=[0.0014    0.0141   -0.0001   -0.0004    0.0002   
           0.0029    0.0198   -0.0002   -0.0119   -0.0074 
           0.0087    0.0209   -0.0133   -0.0112   -0.0077  
           -0.0043   0.0319    0.0009   -0.0078   -0.0084  
           -0.0022    0.0201   -0.0018   -0.0190   -0.0087];
xdot=[A*x(1:n)+B*u
%       A1*x(n+1:2*n)+Whut*tanh(Vhut*[x(n+1:2*n);u])+B*u+L*(C*x(1:n)-C*x(n+1:2*n))
% A1*x(n+1:2*n)+Whut*tanh([x(n+1:2*n);u])+B*u+L*(C*x(1:n)-C*x(n+1:2*n))
A1*x(n+1:2*n)+Whut*tanh([x(n+1:2*n)])+B*u+L*(C*x(1:n)-C*x(n+1:2*n))
      Whutdot0
      Vhutdot0
      A1*x(n+1:2*n)+Whutknown*tanh([x(n+1:2*n)])+B*u+L*(C*x(1:n)-C*x(n+1:2*n))];

    
    
    
    
    
    
    
    
    
    