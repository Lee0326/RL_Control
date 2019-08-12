%closed loop form algorithm function for output feedback TIMC paper
%Version 1, 5/13/05 
%Version 2, 5/14/05
%Version 3, 6/25/05
clear all
close all
clc
%--------------------------------------------------------
%A, B, C : State space description 

A=[-1.01887     0.090506        -0.00215    0           0;
    0.82225     -1.07741        -0.17555    0           0; 
    0           0               -20.2       0           0;
    10          0               0           -10         0;
    -16.26      -0.9788         .04852      0           0];

B=[0            0               20.2        0           0]';

C=[0            0               0           57.2958     0;
    0           57.2958         0           0           0;
    -16.26      -0.9788         0.04852     0           0;
    0           0               0           0           1];
%--------------------------------------------------------
%G, F, H : State space formulation for Feedback design Stewens and Lewis 

G=[0            0               0           0           1]';


F=[0            0               1           0]';

H=[16.26        0.9788          -0.04852    0           0];
%--------------------------------------------------------
%Stewens and Lewis Section 5.4
%Example 5.4.1
%Q, R: Performance index matrices 
%K Optimal output feedback Gain matrix for specified Q and R
Q=[264 16 1 0 0;
    16 60 0 0 0;
    1 0 0 0 0;
    0 0 0 0 0;
    0 0 0 0 100];
R=.01;

K=[-1.629       -1.316          18.56       77.6];

%for simulations without disturbance
D=[0            0               0           0           0;
   0            0               0           0           0 ]';

d=[0 0]';
ti=0;tf=10;
tspan=[ti tf];
x0=[0 0 0 0 0]'; % initial state vector
[t,x]= ode23('f16',tspan,x0,[],K,d,A,B,C,D);
zoptwodist=H*x(:,:)';
plot(t,H*x(:,:)');
z1=H*x(:,:)';
title('without disturbance')
grid on
figure

% Optimal output feedback Gain matrix  taken from Stevens and Lewis
K=[-1.629       -1.316          18.56       77.6];

%Introducing disturbance in simulations
D=[0            0               1           0           0;
   0            0               0           1           0 ]';
d=10*[1 1]';

[t,x]= ode23('f16',tspan,x0,[],K,d,A,B,C,D);
zoptdist=H*x(:,:)';
hold on
plot(t,zoptdist,'k');
grid on
%***************************************************************
%H-Infinity static output feedback closed form with initial 
%stabilizing gains
%J. Gadewadikar, F. Lewis, Murad Abu-Khlaf.
%Copyright F. L. Lewis.
%***************************************************************

% Optimal output feedback Gain matrix  taken from Stevens and Lewis
% These gains will be used as initial stabilizing gains
K=[-1.629       -1.316          18.56       77.6];
%Ko :Initial stabilizing gain
%Define Ko such that Re(eig(A-B.Ko))<0
%K=Ko
% Defining disturbance attenuation
gamma =2;
Ac=A-B*K*C;
B_are=-inv(gamma^2)*D*D';
L=zeros(size(B')); %L_o=0;
i=1;            % Starting Index
Id=eye(size(C'*inv(C*C')*C));
Cplus=C'*inv(C*C');% Pseudo inverse

 % First run
     %------------------------------------------
     P=are(Ac,B_are,Q+C'*K'*R*K*C); % calculate P
     Pold=P
     Pnew=zeros(size(P));
     Lold=L; % Save initial L
     K=inv(R)*(B'*P+L)*Cplus;
     Lnew=R*K*C-B'*P;% Save new L
     L=Lnew
     Ac=A-B*K*C;
     i=i+1; 
     %------------------------------------------
        % Iteration
        %--------------------------------------------
            while norm((Pnew-Pold),2)>.0001
            Pold=P;    
            P=are(Ac,B_are,Q+C'*K'*R*K*C);
            Pnew=P;
            Lold=L;% Save old L before calculating new L
            K=inv(R)*(B'*P+L)*Cplus;
            Lnew=R*K*C-B'*P;% Save new L
            L=Lnew; 
            Ac=A-B*K*C;
            i=i+1; 
            end
        %------------------------------------------
        % Iteration finished
      %--------------------------------------------
i

   K=inv(R)*(B'*P+L)*C'*inv(C*C') % calculate K
   poles=eig(A-B*K*C) % closed loop poles

[t,x]= ode23('f16',tspan,x0,[],K,d,A,B,C,D);
hold on
zhinf=H*x(:,:)';
plot(t,zhinf,'r');
xlabel('Performance Output');
title('H-infinity "red";Optimal OPFB with dis "Black"')
grid on




