%
close all;
clear;
clc;
%Example 5.4.1 Aircraft Control and Simulation [Stevens and Lewis]
%Version 1, 11/23/04
clear all
close all
%----------------------------------------------------------------
% xdot=Ax+Bu+Gr, y=Cx+Fr, z=Hx
A=[-1.01887     0.090506        -0.00215    0   0;
    0.82225     -1.07741        -0.17555    0   0; 
    0           0               -20.2       0   0;
    10          0               0           -10 0;
    -16.26      -0.9788         .04852      0   0];

B=[0 0 20.2 0 0]';

C=[0 0 0 57.2958 0;
    0 57.2958 0 0 0;
    -16.26 -0.9788 0.04852 0 0;
    0 0 0 0 1];

D=[0 0 1 0 0]';
% help are
% 
%  ARE  Algebraic Riccati Equation solution.
%  
%      X = ARE(A, B, C) returns the stablizing solution (if it
%      exists) to the continuous-time Riccati equation:
%  A'*X + X*A - X*B*X + C = 0
%  
%      assuming B is symmetric and nonnegative definite and C is
%      symmetric.

gamma=.2;
R=[.1];
L=0*[1 1 1 1 1];

Q=[264     16      1       0       0; 
    16     60      0       0       0; 
    1       0       0       0       0; 
    0       0       0       0       0; 
    0       0       0       0       100];
B_are=(-(1/gamma^2)*D*D'+B*inv(R)*B');
Q_are=Q+L'*inv(R)*L;
i=1;            % Starting Index
Id=eye(size(C'*inv(C*C')*C));
Cplus=C'*inv(C*C');% Pseudo inverse
 % First run
     %------------------------------------------
     P=are(A,B_are,Q+L'*inv(R)*L); % calculate P
     Lold=L; % Save initial L
     L=L*Cplus*C-B'*P*(Id-Cplus*C);
     Lnew=L; % Save new L
     i=i+1; 
     %------------------------------------------
        % Iteration
        %--------------------------------------------
            while norm((Lnew-Lold),2)>.0001
            P=are(A,B_are,Q+L'*inv(R)*L); % calculate P
            Lold=L;% Save old L before calculating new L
            i
            lcplus=L*Cplus % Contribution because of Ln*Cplus
            L=L*Cplus*C-B'*P*(Id-Cplus*C);
            Lnew=L;% update new L for comparision  
            i=i+1; 
            end
i

   K=inv(R)*(B'*P+L)*C'*inv(C*C') % calculate K
poles=eig(A-B*K*C)
d=1; % define disturbance input


res=P*A+A'*P+Q+(1/gamma^2)*P*D*D'*P-P*B*inv(R)*B'*P+L'*inv(R)*L;
norm(res)
 
%  poles =
% 
%  -66.2764          
%   -1.6218 + 1.0922i
%   -1.6218 - 1.0922i
%   -8.3558          
%  -10.0000      
% Time span
ti=0;
tf=10;
tspan=[ti tf];
x0=[.1 1 0 0 0]'; % initial state vector
[t,x]= ode23('f16',tspan,x0,[],K,d,A,B,C,D);
tinf=t;
xinf=x;
figure
plot(t(:,1),x(:,2),'k','LineWidth',2)
xlabel('time response')
ylabel('Pitch Rate')
grid on
figure
plot(t(:,1),x(:,1),'k','LineWidth',2)
xlabel('time response')
ylabel('Angle of Attack')
grid on


% %K=[-1.629 -1.316 18.56 77.6]
% K=[0.006 -0.152 1.17 .996]
% [t,x]= ode23('f16',tspan,x0,[],K,d,A,B,C,D);
% th2=t;
% xh2=x;
% figure
% plot(t(:,1),x(:,2),'k','LineWidth',2)
% xlabel('time response')
% ylabel('Pitch Rate')
% title('H2 with disturbance')
% grid on
% figure
% plot(t(:,1),x(:,1),'k','LineWidth',2)
% xlabel('time response')
% ylabel('Angle of Attack')
% title('H2 with disturbance')
% grid on
%      
% figure
% plot(tinf(:,1),xinf(:,2),th2(:,1),xh2(:,2),'k','LineWidth',1)
% legend('Hinf','H2');
% title('pitch rate')
% grid on
% figure
% plot(tinf(:,1),xinf(:,1),th2(:,1),xh2(:,1))
% legend('Hinf','H2');
% title('Angle of Attack')
% grid on



