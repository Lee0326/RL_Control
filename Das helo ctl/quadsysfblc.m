function [sdot, u, tau, etad, zetad] = quadsysfblc(t,s)
%
global m J g Ix Iy Iz l
global C1 C2 C3 K1 K2 K3 lamda1 lamda2  
global AA BB CC DD P1 Q1
global x0 y0 z0 xtf ytf ztf tf

%% Quadrotor System 
zeta=s(1:3,1);
eta=s(4:6,1);
zetadot=s(7:9,1);
etadot=s(10:12,1);
%
x=zeta(1,1);
y=zeta(2,1);
z=zeta(3,1);
%
xdot=zetadot(1,1);
ydot=zetadot(2,1);
zdot=zetadot(3,1);
%
phidot=etadot(1,1);
thetadot=etadot(2,1);
sidot=etadot(3,1);
% 
phi=eta(1);
theta=eta(2);
si=eta(3);
%%Diff parameters
sd1=s(13:15,1);
sd2=s(16:18,1);
%%integration parameter
r1i=s(19:22,1);
r2i=s(23:24,1);
%% Trajectory Generation (see the function trajectoryplaning)
[zetad, zetaddot1, zetaddot2] = trajectoryplaning(t,tf,x0,y0,z0,xtf,ytf,ztf);
%
xd = zetad(1,1);
yd = zetad(2,1);
zd = zetad(3,1);
%
xddot=zetaddot1(1,1);
yddot=zetaddot1(2,1);
zddot=zetaddot1(3,1);
%% Thrust Calculation required for calculating 'm/u' while m=1 
vn1=zetaddot2(3,1)+K1(1,1)*(zddot-zdot)+K2(1,1)*(zd-z)+K3(1,1)*r1i(1,1);
u_t = (vn1+g)/(cos(theta)*cos(phi));
um=u_t/m; % um denotes u/m and m=1
U_M=[um,0;0,um];
M_U=inv(U_M);
% 
e2=[xd-x;yd-y];
e2dot=[xddot-xdot;yddot-ydot];
r2=e2dot+lamda2*e2;
thetad1 = (m/u_t)*(-zetaddot2(1,1)-10*C1(1,1)*e2dot(1,1)-10*C2(1,1)*e2(1,1));
phid1 = (m/u_t)*(zetaddot2(2,1)+C1(2,2)*e2dot(2,1)+C2(2,2)*e2(2,1));
etad00 = [thetad1;phid1]+C3*r2i;
thetad=etad00(1,1);
phid=etad00(2,1);
sid=0;
% desired attitude
etad=[phid;thetad;sid];
%% differentiating etad to obtain etaddot1
diffdot1 = AA*sd1 + BB*etad;
etaddot1 = CC*sd1 + DD*etad;
%
phiddot=etaddot1(1,1);
thetaddot=etaddot1(2,1);
siddot=etaddot1(3,1);
%% differentiating etaddot1 to obtain etaddot2
diffdot2 = AA*sd2 + BB*etaddot1;
etaddot2 = CC*sd2 + DD*etaddot1;
%% Calculating the Mh and Eh matrix required to compute FBLC input
% U = inv(Eh) * (vnew - Mh);
Mh = [-g; thetadot*sidot*(Iy-Iz)/Ix; phidot*sidot*(Iz-Ix)/Iy; phidot*thetadot*(Ix-Iy)/Iz];
Eh = [cos(theta)*cos(phi)/m 0 0 0;0 l/Ix 0 0;0 0 l/Iy 0; 0 0 0 l/Iz];
% e1 = [zd - z;phid - phi;thetad - theta;sid - si];
e1dot=[zddot - zdot;phiddot - phidot;thetaddot - thetadot;siddot - sidot];
e1=[zd - z;phid - phi;thetad - theta;sid - si];
r1=e1dot+lamda1*e1;
% Computing A_tilda 
Atil=[sin(theta)-thetad;phid-cos(theta)*sin(phi)];
% Comuting the robustifying term. Note that the robustifying component are
% only available for phi and theta channels.
%
vrp = um*r2(2,1)*Q1(2,2)*Atil(2,1)/(r1(2,1)*P1(2,2)+0.001);
vrt = um*r2(1,1)*Q1(1,1)*Atil(1,1)/(r1(3,1)*P1(3,3)+0.001);
Vr = [0;vrp;vrt;0];
% Vr = r1*r2'*U_M*Atil/(norm(r1)^2+0.001);
% Computing new input for FBLC or the linear PID like controller for
% feedback linearized system
vnew=[zetaddot2(3,1);etaddot2]+K1*e1dot+K2*e1+K3*r1i+Vr;
% Computing the Feedback linearizing input
U=inv(Eh)*(vnew-Mh);
%% Distributing the inputs into thrust 'u' and angular torque components
% 'tauphi', 'tautheta' and 'tausi'
u=U(1,1);
tauphi=U(2,1);
tautheta=U(3,1);
tausi=U(4,1);
tau=[tauphi;tautheta;tausi];
% Computing the angular velocities of each rotor
w4=abs(sqrt(0.25*(u-tausi+2*tauphi)));
w2=abs(sqrt(w4^2 - tauphi));
w3=abs(sqrt(0.25*(u+tausi+2*tautheta)));
w1=abs(sqrt(w3^2 - tautheta));
% Added disturbance to the quadrotor as 'OMEGA'
omega=w2+w4-w1-w3;
%% Plant Dynamics
% Transational dynamics (xdot2, ydot2, zdot2)
zetadot2=[-(u/m)*sin(theta);(u/m)*cos(theta)*sin(phi);(u/m)*cos(theta)*cos(phi)-g];
% Rotational Dynamics (phidot2, thetadot2 and sidot2)
etadot2=[thetadot*sidot*(Iy-Iz)/Ix - J*thetadot*omega/Ix; phidot.*sidot*(Iz-Ix)/Iy + J*phidot*omega/Iy; phidot.*thetadot*(Ix-Iy)/Iz]+...
    [l/Ix 0 0;0 l/Iy 0;0 0 l/Iz]*[tauphi;tautheta;tausi];
% Collecting all states
sdot=[zetadot;etadot;zetadot2;etadot2;diffdot1;diffdot2;r1;r2];
%
%% Displaying time and torque in command window
disp([t tau']);






    