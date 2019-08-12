function[K]=Hinfopfb(A,B,C,D,gamma,Q,R) 
%-------------------------------------------------------------------------
%This algorithm solves two coupled necessary and sufficient 
%design equations of H-Infinity output-feedabck problem posed in
%Gadewadikar, Lewis and Abu-Khalaf, "Necessary and Sufficient Conditions for
%H-Infinity Static Output-Feedback Control", Journal of Guidance, Dynamics
%and Control
%-------------------------------------------------------------------
% Design Equation 1 
%PA+A'P+Q+(1/gamma^2)PDD'P-PBinv(R)B'P+L'inv(R)L =0    (1)
%Design Equation 2
%KC=invR(B'P+L)                                         (2)
%-------------------------------------------------------------------

%-------------------------------------------------------------------
% System Definitions
%System or Plant Matrix A
%Control Input matrix B
%Mesurement matrix C

%Disturbance matrix D Note: D is not a direct feed matrix

%gamma= disturbance Attenuation; % Define L2 gain gamma

%-------------------------------------------------------------------
%ARE  routine is used to solved equation 1
%  ARE  Algebraic Riccati Equation solution.
%  
%      X = ARE(A, B, C) returns the stablizing solution (if it
%      exists) to the continuous-time Riccati equation:
%  A'*X + X*A - X*B*X + C = 0
%  
%      assuming B is symmetric and nonnegative definite and C is
%      symmetric.
%-------------------------------------------------------------------
%   Author   : J. Gadewadikar 11-21-2004
%   Revised    J. Gadewadikar  04-26-06
%   jyotir@arri.uta.edu

L= zeros(size(B')); %Start with zero initial gain L with dimensions m by n

B_are=(-(1/gamma^2)*D*D'+B*inv(R)*B');

i=1;            % Starting Index

Cplus=C'*inv(C*C');% Define Pseudo inverse Cplus
CplusC=Cplus*C;
Id=eye(size(Cplus*C));
Prj=Id-CplusC;

%-------------------------------------------------------------------
%start with zero index n=0
%Start with L(0)=0;
%1. Solve ARE
%P(n)A+A'P+Q+(1/gamma^2)P(n)DD'P(n)-P(n)Binv(R)B'P(n)+L(n)'inv(R)L(n) =0
%2. Update L
%L(n+1)=B'P(n)(I-CplusC)+L(n)CplusC
%3. Check convergence 
%If not converged increment index and go to 1.
%If converged 
%calculate and display output-feedback gain matrix and closed loop poles
%-------------------------------------------------------------------

 % First run
     %-------------------------------------------------------------
     P=are(A,B_are,Q+L'*inv(R)*L); % calculate P
     Lold=L; % Save initial L
     L=-B'*P*Prj+L*CplusC;
     Lnew=L; % Save new L
     i=i+1; 
     %------------------------------------------------------------
        % Iteration checks convergence on L using Frobenius norm
        %---------------------------------------------------------
            while norm((Lnew-Lold),'fro')>.01 
            P=are(A,B_are,Q+L'*inv(R)*L); % calculate P
            Lold=L;% Save old L before calculating new L
            L=-B'*P*Prj+L*CplusC;
            Lnew=L;% update new L for comparision  
            i=i+1;
            end
    %---------------------------------------------------------
 
%-----------------------------------------------------
   K=inv(R)*(B'*P)*C'*inv(C*C') % calculate K
%--------------------------------------------------------------------------
