clear all
close all
clc
%Model from text file

%x= [u           v   p           q       phi         theta   as          bs w       r   rbf]';

A= [  -0.1778         0         0         0         0   -9.7807   -9.7807;
         0   -0.3104         0         0    9.7807         0         0;
   -0.3326   -0.5353         0         0         0         0   75.7640;
    0.1903   -0.2940         0         0         0         0  172.6200;
         0         0    1.0000         0         0         0         0;
         0         0         0    1.0000         0         0         0;
         0         0         0   -1.0000         0         0   -8.1222;
         0         0   -1.0000         0         0         0   -0.0921;
         0         0         0         0         0         0   17.1680;
         0         0   -0.2834         0         0         0         0;
         0         0         0         0         0         0         0];
     
     A(:,8:11)=[0         0         0         0;
                9.7807         0         0         0;
                343.8600         0         0         0;
                -59.9580         0         0         0;
                0         0         0         0;
                0         0         0         0;
                4.6535         0         0         0;
                -8.1222         0         0         0;
                7.1018   -0.6821   -0.1070         0;
                0   -0.1446   -5.5561  -36.6740;
                0         0    2.7492  -11.1120];
eigA=eig(A);
%U=[Del_lati    Del_longi   Del_col     Del_ped]';

   B1=[  0         0         0         0;
         0         0         0         0;
         0         0         0         0;
         0         0         0         0;
         0         0         0         0;
         0         0         0         0;
    0.0632    3.3390         0         0;
    3.1739    0.2216         0         0;
         0         0   19.9250         0;
         0         0    2.0816  -74.3640;
         0         0         0         0];

B=[B1(:,1) B1(:,2) B1(:,4)];

%x= [u           v   p           q       phi         theta   as          bs w       r   rbf]';
%y=[phi theta r p q]'

C=zeros(5,11);

C(1,5)=1;C(2,6)=1;C(3,10)=1;C(4,3)=1,
C(5,4)=1;

Dsys=zeros(5,3);
G=ss(A,B,C,Dsys);
 
%  num11=.6*[1 1];
%  den11=[1 0];
%  B1=tf(num11,den11);
%  num22=.5*[1 .5];
%  den22=[1 0];
%  num33=.5*[1 .6];
%  den33=[1 0];

% Filter values as on 07/25/06
%num11= [1 .1];
%den11= [1 2 0];

num11= 3*[1 0.5];
den11= [1 5 0];

B1=tf(num11,den11);
 
 
num22= num11;%[1 .1];
den22=den11;%[1 1.25 0];

num33=num11;%[1 0.1];
den33=den11;%[1 1.25 0]; 


 W1 = [ tf(num11,den11)     0              0;
        0              tf(num22,den22)    0;
        0               0                  tf(num33,den33)];
 
W1Gt= series(W1,G);
%y=[phi theta r p q]'

W2=diag([1 1 1 1 1]);
tauc=3.2;
Aw=[-1/tauc 0; 0 -1/tauc];
Bw=.5*eye(2);Cw=eye(2);Dw=zeros(2,2);
WindG=ss(Aw,Bw,Cw,Dw);

W1G=series(W1Gt,W2);
figure
sigma(W1G,'-',G,'*', WindG,':');
xlim([.01 50])
hold on
legend('Loop-Shaped Plant','Open Loop Plant','Wind Model')
grid on


%=sigma(W1G,'r',G,'y--',WindG,'gx')
%[SV,W]=SIGMA(W1G,{.001,10})
Q=zeros(size(W1G.a));

%x= [u           v              p           q       phi         theta   as          bs            w       r               rbf]';
    Q(1,1)=2500;Q(2,2)=2500;Q(3,3)=100;   Q(4,4)=100;  Q(5,5)=1e6;  Q(6,6)=1e6; Q(7,7)=1; Q(8,8)=1; Q(9,9)=2500;Q(10,10)=100; Q(11,11)=100;
     %                       Q(3,5)=100;    Q(4,6)=1;   Q(5,3)=100; Q(6,4)=1;
                            
 %Q(12,12)=1000;Q(13,13)=10000;Q(14,14)=10000;

R=zeros(3,3);

R=[1690000         0        0;
      0        1690000      0;  
       0         0           7800];

Q=.0001*Q;

gamma=.62;

R=.0001*R;

D=[(W1G.a(:,1)) (W1G.a(:,2))] 



K=Hinfopfb(W1G.a,W1G.b,W1G.c,D,gamma,Q,R);

W1Gc=ss(W1G.a-W1G.b*K*W1G.c,W1G.b,W1G.c,W1G.d);
figure
sigma(G,W1G,'-x',W1Gc,'-+',WindG,':');
%setoptions(j,'FreqUnits','Hz');
legend('W1G','G','W1Gc','WindG')
grid on

Phicorr=1;
Thetacorr=1;

%Phicorr=1/1.18906;
%Thetacorr=1/1.1815;

qu=-20*0.3048;
qv=-20*0.3048;
sim withshapingprgmdistsim

%x= [u           v   p           q       phi         theta   as          bs w       r   rbf]';
%y=[phi theta r p q]'

% Plot Results for Phi Command Variables
figure
subplot(221)
plot(y(:,1),y(:,2),y(:,1),ones(length(y(:,1)),1),'--r','LineWidth',2);ylabel(' \phi in radians')
%title('Lat - Dir state histories for \phi_{command}')
subplot(223)
plot(y(:,1),y(:,5),'LineWidth',2);ylabel(' P in radian/s');xlabel('time in seconds')
subplot(222)
plot(x1(:,1),x1(:,3),'LineWidth',2);ylabel(' V in m/s');
subplot(224)
plot(x1(:,1),x1(:,10),'LineWidth',2);ylabel(' R in rad/s');xlabel('time in seconds')

figure
subplot(221)
plot(y(:,1),y(:,3),'LineWidth',2);ylabel(' \theta in radians')
%title('Long. state histories for \phi_{command}')
subplot(223)
plot(y(:,1),y(:,6),'LineWidth',2);ylabel(' Q in radian/s');xlabel('time in seconds')
subplot(222)
plot(x1(:,1),x1(:,2),'LineWidth',2);ylabel(' U in m/s')
subplot(224)
plot(x1(:,1),x1(:,10),'LineWidth',2);ylabel(' W in m/s');xlabel('time in seconds')

figure
subplot(311)
plot(InertialPos(:,1),InertialPos(:,2),'LineWidth',2);ylabel('X in m');grid;
title('Inertial Positions for the \phi_{command}')
subplot(312)
plot(InertialPos(:,1),InertialPos(:,3),'LineWidth',2);ylabel('Y in m');grid;
subplot(313)
plot(InertialPos(:,1),InertialPos(:,4),'LineWidth',2);ylabel('Z in m');grid;
xlabel('time in seconds')

% Plot Controls & Misc variables
load control_lat.mat
ulat = ulat';
figure
subplot(311),plot(ulat(:,1),ulat(:,2),'LineWidth',2); ylabel('Cyc. pitch - lat.');title('Control histories for \phi_{command}')
subplot(312),plot(ulat(:,1),ulat(:,3),'LineWidth',2); ylabel('Cyc. pitch - long.')
subplot(313),plot(ulat(:,1),ulat(:,4),'LineWidth',2); ylabel('Pedal');xlabel('time in seconds')

figure
subplot(311),plot(x1(:,1),x1(:,8),'LineWidth',2); ylabel('as')
subplot(312),plot(x1(:,1),x1(:,9),'LineWidth',2); ylabel('bs')
subplot(313),plot(x1(:,1),x1(:,12),'LineWidth',2); ylabel('rbf');xlabel('time in seconds');


% Plot Pitch angle command
figure
subplot(221)
plot(y2(:,1),y2(:,3),y2(:,1),ones(length(y2(:,1)),1),'--r','LineWidth',2);ylabel(' \theta in radians')
%title('Long. state histories for \theta_{command}')
subplot(223)
plot(y2(:,1),y2(:,6),'LineWidth',2);ylabel(' Q in radian/s');xlabel('time in seconds')
subplot(222)
plot(x2(:,1),x2(:,2),'LineWidth',2);ylabel(' U in m/s')
subplot(224)
plot(x2(:,1),x2(:,10),'LineWidth',2);ylabel(' W in m/s');xlabel('time in seconds')

figure
subplot(221)
plot(y2(:,1),y2(:,2),'LineWidth',2);ylabel(' \phi in radians')
%title('Lat - Dir state histories for \theta_{command}')
subplot(223)
plot(y2(:,1),y2(:,5),'LineWidth',2);ylabel(' P in radian/s');xlabel('time in seconds')
subplot(222)
plot(x2(:,1),x2(:,3),'LineWidth',2);ylabel(' V in m/s');
subplot(224)
plot(x2(:,1),x2(:,10),'LineWidth',2);ylabel(' R in rad/s');xlabel('time in seconds')

figure
subplot(311)
plot(InertialPos1(:,1),InertialPos1(:,2),'LineWidth',2);ylabel('X in m');grid;
title('Inertial Positions for the \theta_{command}')
subplot(312)
plot(InertialPos1(:,1),InertialPos1(:,3),'LineWidth',2);ylabel('Y in m');grid;
subplot(313)
plot(InertialPos1(:,1),InertialPos1(:,4),'LineWidth',2);ylabel('Z in m');grid;
xlabel('time in seconds')

% Plot Controls & Misc variables
load control_long.mat
ulong = ulong';
figure
subplot(311),plot(ulong(:,1),ulong(:,2),'LineWidth',2); ylabel('Cyc. pitch - lat.');title('Control histories for \theta_{command}')
subplot(312),plot(ulong(:,1),ulong(:,3),'LineWidth',2); ylabel('Cyc. pitch - long.')
subplot(313),plot(ulong(:,1),ulong(:,4),'LineWidth',2); ylabel('Pedal');xlabel('time in seconds')

figure
subplot(311),plot(x2(:,1),x2(:,8),'LineWidth',2); ylabel('as')
subplot(312),plot(x2(:,1),x2(:,9),'LineWidth',2); ylabel('bs')
subplot(313),plot(x2(:,1),x2(:,12),'LineWidth',2); ylabel('rbf');xlabel('time in seconds');

figure
plot(d(:,1),d(:,2),d(:,1),d(:,3),'--','LineWidth',2); xlabel('time in seconds'); 
title('Disturbance velocities in the Body fixed X & Y axes (m/s)'); legend('du','dv');