%function odestart
%xytable 
clear all;close all;clc;
global P R A B C D Q rowA K F r G;
global u;
global k
global T
%system matrices
A=[0  1   0 0;...
   0 -11  0 0 ;...
  -1  0   0 0;...
   0  1/0.005 0 -1/0.005];

B=[0 1 0 0]';
C=[0 0 1 0;...
   -1 0 0 0;...
    0 0 0 1];

D=0;
%r=[0 0 0 1 0];


figure(1);hold on;

[rowA columnA]=size(A);
[rowB columnB]=size(B);

statesnum=columnA;
inputsnum=columnA;
unknownnum=rowA*(rowA+1)/2;


%***********************************
%Need to reset
%initial conditions
x0=[20,-10,10,-20,0,20,-10,10,-20];%size(x0)=statesnum+1
P=zeros(rowA,rowA);
% P=eye(rowA,rowA);
% P=[55.5177    4.9043  -11.5240   -0.0002
%     4.9043    0.5219   -1.0000    0.0024
%   -11.5240   -1.0000    4.9044    0.0000
%    -0.0002    0.0024    0.0000    0.0025];
% P=[55    5  -10   0
%     5    0.5   -1    0
%   -10   -1    5    0.0000
%    0    00    0.0000    0.0025];
% L=0.1*eye(rowA,rowA);
L=0.001*[1 1 1 1];
NN=20;
N=30;
e=0.001;
Fsamples=N*unknownnum; %length of the simulation in samples 
T=0.02; % sample time
%*************************************888


Q=[10 0 0 0;0 10 0 0;0 0 105 0;0 0 0 10];
R=1;

uu=[]; % saving the control signal for plot

%parameters of the critic rearanged as it would be returned by the least
%squares
%WW=[P(1,1); 2*P(1,2); 2*P(1,3); P(2,2); 2*P(2,3); P(3,3)];

WW=zeros(unknownnum,1);
k=1;
    for i=1:rowA 
        for j=i:rowA
            if i==j
                WW(k)=P(i,j);
            else
                WW(k)=2*P(i,j);
            end
             k=k+1;
        end
    end 
   
WWP=[WW; 0];



%parameters for the batch least squares
L1=0.001*[1 1 1 1];
x1c=[];
x11c=[];
Xpi=[];
Xpi1=[];
K=[-1 -10 0];
% K=inv(R)*(B'*P+L1)*C'*inv(C*C')
E=[eig(A-B*K*C)] % saves the poles of the closed loop system over the duration of the simulation
damp(A-B*K*C)



upd=[]; % stores information relative to updates of the critic parameters
ch=0;
j=0;
mm=1;
for k=1:Fsamples,
    j=j+1;
    l=1;
    for m=1:rowA
        for n=m:rowA
            X0(j,l)=x0(m)*x0(n);
            l=l+1;
        end
    end

       l=1;
    for m=1:rowA
        for n=m:rowA
            X01(j,l)=x0(m+rowA+1)*x0(n+rowA+1);
            l=l+1;
        end
    end  

    before_cost=x0(1:rowA)*P*x0(1:rowA)';
    tspan=[0 T];
    [t,x]= ode45('subJO1',tspan,x0);
    
   x1=zeros(1,rowA);
   x1=x(length(x),1:rowA);
   mm=length(x);
   for mmm=1:mm
        x1c((k-1)*mm+mmm,:)=x(mmm,1:rowA);
   end

   Y(j,:)=x(length(x),rowA+1);
    
   x11=zeros(1,rowA);
   x11=x(length(x),rowA+2:rowA*2+1);
   for mmm=1:mm
        x11c((k-1)*mm+mmm,:)=x(mmm,rowA+2:rowA*2+1);
    end

    after_cost=x(length(x),rowA+1)+x1*P*x1';
    
    l=1;
    for m=1:rowA
        for n=m:rowA
            X1(j,l)=x1(m)*x1(n);
            l=l+1;
        end
    end
    
    l=1;
    for m=1:rowA
        for n=m:rowA
            X11(j,l)=x11(m)*x11(n);
            l=l+1;
        end
    end
    
    Xpi(j,:)=X0(j,:)-X1(j,:);
    Xpi1(j,:)=X01(j,:)-X11(j,:);
      
          
   
    plot(t+T*(k-1),x(:,6),'b','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,7),'r','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,8),'g','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,9),'black','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,1),':b','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,2),':r','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,3),':g','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,4),':black','LineWidth',2); hold on 
    
    for i=1:rowA
        x0(i)=x(length(t),i);
    end
   for i=rowA+2:rowA+rowA+1
        x0(i)=x(length(t),i);
   end
   
% plot(t(length(t))+T*(k-1),x1(1:rowA),'o');hold on
    
% plot(t(length(t))+T*(k-1),x11(1:rowA),'o');

    uu=[uu u];
    
    if (abs(after_cost-before_cost)>e)&&(ch==0)&&(mod(j,unknownnum+NN)~=1),
        j=0;
        ch=ch+1;
    else
        if abs(after_cost-before_cost)>e,
        ch=ch+1;
        end
    end

    % the batch least squares is made on 6 values
    if mod(j,unknownnum+NN)==0,
       if (abs(after_cost-before_cost)>e)&&(ch==unknownnum+NN),
          weights=Xpi\Y ;%calculation of the weights
          upd=[upd 1];
       else 
        %there is no reason to update
        upd=[upd 0];       
       end
         WWP=[WWP [weights; k*T]];
         WW=[WW weights];
         %calculating the matrix P
         %P=[weights(1) weights(2)/2 weights(3)/2 ; weights(2)/2 weights(4) weights(5)/2; weights(3)/2 weights(5)/2 weights(6)]
            
         l=1;
       for m=1:rowA 
           for n=m:rowA
               if m==n
                   P(m,n)=weights(rowA*(m-1)-(m-1)*m/2+n);
               else
                   P(m,n)=weights(rowA*(m-1)-(m-1)*m/2+n)/2;
            end
            
            P(n,m)=P(m,n);
            l=l+1;
        end
    end 
 
        j=0;
        ch=0;
     
        K=inv(R)*(B'*P+L1)*C'*inv(C*C');
        damp(A-B*K*C);
        E=[E eig(A-B*K*C)];
         mm=mm+1;
%          L=R*K*C-B'*P;
%         Q=C'*Q*C+L'*inv(R)*L;
    end
end
%%

% figure(1); title('System states'); xlabel('Time(s)');legend('x(1)','x(2)','x(3)','x(4)');
figure(1); title('System states'); xlabel('Time(s)');
str1='$$\hat{x}$$';
str2='$$\hat{\dot{x}}$$'; 
str3='$$\hat{\epsilon}$$';
str4='$$\hat{\dot{x}}_F$$'; 
str5='$$x$$';
str6='$$\dot{x}$$'; 
str7='$$\epsilon$$';
str8='$$\dot{x}_F$$'; 
legend({str1,str2,str3,str4,str5,str6,str7,str8},'interpreter','latex','fontsize',20)
K=inv(R)*(B'*P+L1)*C'*inv(C*C')
P
%%
K
K1=[ -10.2470  -16.8947    0.0234]
Ac1=A-B*K1*C;
Ac=A-B*K*C;
damp(Ac)
x01=x0(1:rowA);
t=(0:0.05:10);
u=zeros(size(t));
[y,x]=lsim(Ac,B,C,D,u,t,x01);
[y1,x1]=lsim(Ac1,B,C,D,u,t,x01);

figure(2);

plot(t,y(:,1),'b','LineWidth',2);hold on
plot(t,y(:,2),'r','LineWidth',2);hold on
plot(t,y(:,3),'g','LineWidth',2);hold on
plot(t,y1(:,1),':b','LineWidth',2);hold on
plot(t,y1(:,2),':r','LineWidth',2);hold on
plot(t,y1(:,3),':g','LineWidth',2);hold on
title('Outputs y'); xlabel('Time(s)');
str1='$$\epsilon (Algotirhm 3)$$';
str2='$$e (Algotirhm 3)$$'; 
str3='$$ \dot{x}_F (Algotirhm 3)$$'; 
str4='$$ \epsilon (Algotirhm 2)$$';
str5='$$e (Algotirhm 2)$$'; 
str6='$$\dot{x}_F (Algotirhm 2)$$';
legend({str1,str2,str3,str4,str5,str6},'interpreter','latex','fontsize',20)
hold off

%printing for comparison
% sol=care(A,B,eye(rowA),1);
% 
% WW=[WW [sol(1,1); 2*sol(1,2); 2*sol(1,3); sol(2,2); 2*sol(2,3); sol(3,3)]];

figure(3); plot([0:T:T*(length(uu)-1)],uu);
title('Control signal'); xlabel('Time(s)');

%plotting the poles of the closed loop system


figure(5);
hold on;
plot(WWP(unknownnum+1,:)',WWP(1:unknownnum,:)','.-');
title('P matrix parameters'); xlabel('Time(s)');legend('P_1_1','P_1_2','P_1_3','P_1_4','P_2_2','P_2_3','P_2_4','P_3_3','P_3_4','P_4_4');

% figure;
% WWP=[WWP [[sol(1,1); 2*sol(1,2); 2*sol(1,3); sol(2,2); 2*sol(2,3); sol(3,3)]; T*(Fsamples-1)]];
% xy=size(WWP);
% plot(WWP(7,1:(xy(2)-1)),WWP(2:2:6,1:(xy(2)-1)),'.');
% hold on;
% plot(WWP(7,xy(2)),WWP(2:2:6,xy(2)),'*');
% title('P matrix parameters'); xlabel('Time(s)'); axis([0 3 -0.2 2.5]);

figure(6);
plot(upd,'*'); title('P parameters updates'); xlabel('Iteration number');axis([0 10 -0.1 1.1]);
%*************************************************************
%track problem
x0=[0 10 0 0]; 
T=20;
tspan=[0 T];
[t,XXXX]= ode45('cal',tspan,x0);
[M3,N3]=size(XXXX);

YYY=zeros(3,M3);
for i=1:M3
    %Y(:,i)=C*XX(i,:)'+F*r;
    YYY(:,i)=C*XXXX(i,:)'+F*r;%canclute the output state
end
%UU=zeros()

for i=1:M3
    UUUU(i)=-K*(C*XXXX(i,:)'-F*r);%canclute the input command
end

%figure
figure(7)
plot(YYY(1,:),'r','LineWidth',2);
xlabel('time response')
ylabel('error intergral')
grid on

%hold on
figure(8)
plot(YYY(2,:),'r','LineWidth',2);
xlabel('time response')
ylabel('error')
grid on

%hold on
figure(9)
plot(YYY(3,:),'r','LineWidth',2);
xlabel('time response')
ylabel('filter velocity')
grid on
%hold on

figure(10)
plot(XXXX(:,1),'r','LineWidth',2);
xlabel('time response')
ylabel('position')
grid on

figure(11)
plot(XXXX(:,2),'r','LineWidth',2);
xlabel('time response')
ylabel('velocity')
grid on

figure(12)
plot(UUUU,'r','LineWidth',2)
xlabel('time response')
ylabel('input')
grid on



%figure(2)
%plot(Y(4,:),'r');
%-------------------------------------------------------------------------------------

