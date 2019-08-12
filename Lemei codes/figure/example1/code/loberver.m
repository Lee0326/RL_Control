function xdot=loberver(t,x) 
global A B C

%   u=0.02*sin(0.1*t)+0.03*sin(0.2*t)+30*sin(20*t);

Ahut=[-x(11) 1 ;-x(12) 0];
Bhut=[x(9);x(10)];
Khut=[9-x(11); 20-x(12)];

y=C*x(1:2);
yhut=C*x(3:4);
% u=-1.406*y+30*sin(20*t);
 u=30*sin(20*t);
% u=-1.406*y+0.03*sin(0.2*t);
A1=[-9 -20;1 0];
B1=[1;0];
B2=[-1;0];

Z=y+[5,6]*x(7:8);
Zhut=x(9:10)'*x(5:6)+x(11:12)'*x(7:8);
yta=Z-Zhut;

gmma1=0.1;
gmma2=0.2;

xdot=[ A*x(1:2)+B*u
       Ahut*x(3:4)+B*u+Khut*(y-yhut)
         A1*x(5:6)+B1*u
      A1*x(7:8)+B2*y
     gmma1*yta*x(5:6)
      gmma2*yta*x(7:8)
      [ -1.0201 1.0000;-1.3847  0]*x(13:14)+B*u+[78.9799;198.6153]*C*(x(1:2)-x(13:14))];
