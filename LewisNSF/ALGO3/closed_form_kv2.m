%closed loop form algorithm function
%Version 1, 5/13/05 
%Version 2, 5/14/05
%Version 3,10/16/05
%Writtern By J. Gadewadikar.
%MATLAB routine A, B, C : State space description 
%Q=C'*C, R: Performance index matrices

function [P,K,F,L,i] = closed_form_kv2(A,B,C,R,D,gamma)
    [K P]=lqr(A,B,C'*C,R);
    B_are=(-(1/gamma^2)*D*D');
    Cplus=C'*inv(C*C'); % pseudo inverse
    
    [U,S,V] = svd(C);
    C_dim=size(C);
    r=C_dim(1,1); %no of outputs
    n=C_dim(1,2); %states
    V1=V(1:n,1:r);%(n)x(r) first r columns
    
    V2=V(1:n,r+1:n);%(n)x(n-r) last n-r columns
    Pold=P;
    Pnew=zeros(size(P)); % To run iteration first time
    % Update Ac
    Ac=A-B*K;
    i=1;
    
%ITERATION    
    while norm((Pnew-Pold),'fro')>.00001 % Check convergence will run first time
        Pold=P;
        %Step 1 : Solve lyapunov equation for P
        P=are(Ac,B_are,C'*C+K'*R*K);
        Pnew=P;
        % Step 3: Update L
        L=R*K-B'*P;
        %Step 2 : Update K 
        K=inv(R)*(B'*P+L);
        %Step 4: Update Ac
        Ac=A-B*K-B*K*V2*V2';
        i=i+1
        eig(P);
    end 
F=K*Cplus;


