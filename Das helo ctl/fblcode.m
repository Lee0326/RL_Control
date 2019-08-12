% Simulation of backstepping method using ODE23
close all;clear all;
clc;
%
global m M1 J g Ix Iy Iz l
global C1 C2 C3 K1 K2 K3 lamda1 lamda2  
global AA BB CC DD P1 Q1
global x0 y0 z0 xtf ytf ztf tf
%%
dt=0.1;
tf=20;
tt=0:dt:tf;
ini1=[0 0 0]; % Initial condition for x, y, and z as x0, y0, and z0
ini2=zeros(1,21); % Initial condition for rest of the states
ini=[ini1 ini2];
%% Quadrotor Parameters
m=1; % Quadrotor mass
l=1; % Distance between center of the quadrotor and rotor
M1 = m*eye(3); % Mass matrix 
Ix = 5; Iy = 5; Iz = 5; % Innertia components
% J=[Ix 0 0;0 Iy 0;0 0 Iz];
J=2;
g=9.81; % Gravity
mg=[0;0;-m*g];
%% Constant Gains
lamda1 = 10*eye(4);
lamda2 = 10*eye(2);
R0 = 80*eye(4);
S0 = 1*eye(2);
K1 = lamda1+R0;
K2 = lamda1*R0;
K3 = 80*eye(4);
C1 = lamda2+S0;
C2 = lamda2*S0;
C3 = 10*eye(2);
%% Lyapunov gains
P1 = eye(4,4);
Q1 = eye(2,2);
%% Differentiator Parameters
AA = -100*eye(3,3);
BB = eye(3,3);
CC = -100*eye(3,3);
DD = 0.01*eye(3,3); 
%% Trajectory Generation data
% Initial and desired final position
x0=ini1(1,1); y0=ini1(1,2); z0=ini1(1,3);
xtf=20; ytf=15; ztf=10;
%% Integration with ODE45
[t,S] = ode45(@quadsysfblc,tt,ini);
U = []; Tau = [];
for ct = 1:length(t),
    [dummy, u, tau, etad, zetad] = quadsysfblc(t(ct),S(ct,:)');
    U1(ct, 1) = u;
    Tau(ct, :) = tau';
    Zetad(ct,:) = zetad';
    Etad(ct,:)  = etad';
end 
%% Plot simulation result
close all;
figure; 
subplot(2,1,1),plot(t,Zetad(:,1),'b-',t,S(:,1),'r:',t,Zetad(:,2),'c--',t,S(:,2),'k-.',t,Zetad(:,3),'g.',t,S(:,3),'m+','LineWidth',2);grid;zoom;
xlabel('Time in Sec')
ylabel('Positions in meter')
legend('x_d','x','y_d','y','z_d','z')
subplot(2,1,2),plot(t,Zetad(:,1)-S(:,1),'b-',t,Zetad(:,2)-S(:,2),'r:',t,Zetad(:,3)-S(:,3),'k-.','LineWidth',2);grid;zoom;
xlabel('Time in Sec')
ylabel('Position erros in meter')
legend('x_d - x','y_d - y','z_d - z')
%
figure
subplot(2,1,1),plot(t,Etad(:,1),'b-',t,S(:,4),'r:',t,Etad(:,2),'c-.',t,S(:,5),'k--','LineWidth',2);grid;zoom;
xlabel('Time in Sec')
ylabel('Angular Positions in rad')
legend('\phi_d','\phi','\theta_d','\theta')
subplot(2,1,2),plot(t,Etad(:,1)- S(:,4),'b-',t, Etad(:,2)- S(:,5),'r:','LineWidth',2);grid;zoom;
xlabel('Time in Sec')
ylabel('Angular Position errors in rad')
legend('\phi_d - \phi','\theta_d - \theta')
figure
subplot(211),plot(t,U1,'LineWidth',2);grid;zoom;xlabel('Time in Sec');ylabel('Control input U in Newton')
subplot(212),plot(t,Tau(:,1),'-',t,Tau(:,2),':',t,Tau(:,3),'-.','LineWidth',2);grid;zoom;xlabel('Time in Sec');ylabel('Control Torques In N-m');legend('tao_{phi}','tao_{theta}','tao_{si}')
