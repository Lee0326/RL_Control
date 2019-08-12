%F-16
clear all;close all;clc;
global P R A B C D Q rowA K F r G nnn mmm kkk lll;
global u;
global k
global T
%******************************************************
%matrix coefficient
A=[-1.01887 0.90506   -0.00215   0    0;...
    0.82225  -1.07741 -0.17555   0    0;...
    0        0        -20.2         0    0;...
    10       0         0         -10  0;...
    -16.26   -0.9788   0.4852    0   0];
B=[0 0 20.2 0 0]';
G=[0 0 0 0 1]';
C=[0 0 0 57.2958  0;...
   0 57.2958 0 0 0;...
   -16.26 -0.9788 0.4852 0 0;...
   0 0 0 0 1];
F=[0 0 0 1 0];
H1=[16.26 0.9788 -0.04852 0 0];
D=0;
%******************************************************
%******************************************************
%estimate the size of matrix
[rowA columnA]=size(A);
[rowB columnB]=size(B);
[rowC columnC]=size(C);
[rowD columnD]=size(D);
[rowF columnF]=size(F);
[rowG columnG]=size(G);
[rowr columnr]=size(r);
%******************************************************
%******************************************************
Statesnum=columnA;
Inputsnum=columnB;
Outputsnum=rowC;
unknownnum=rowA*(rowA+1)/2;
%******************************************************

%***************************************************
a1=min(rowC,columnC);
if rank(C)<a1
    sprintf('The rank of C is not full row rank',2,3)
    return;
end
%******************************************************

%******************************************************
%Q,R,X
Q=[264 16 1 0 0;16 60 0 0 0;1 0 0 0 0;0 0 0 0 0;0 0 0 0 100];
rou=0.1;
R=rou*eye(Inputsnum);
X=eye(Statesnum);
%******************************************************

%***********************************
%Need to reset
%initial conditions
x0=[1 2 1 -1 1,0,1 2 1 -1 1];%size(x0)=statesnum+1
% P=zeros(rowA,rowA);
P=eye(rowA,rowA);
% L1=0.001*[1 1 1 1 1];
%  L1=0.00*[1 1 1 1 1];
L1=[0.0043    0.0010    0.1126    0.0010    0.0010];
%*************************************888
%*************************************888
NN=50;
N=30;
e=0.0000001;
Fsamples=N*unknownnum; %length of the simulation in samples 
T=0.02; % sample time
%*************************************888
K=[0   0    0   100]
%***********************************
% Check if K0 stable, if not return
damp(A-B*K*C);
    Eig1=eig(A-B*K*C);%calclate the eig of A-B*K0*C
    REig1=real(Eig1);%get the real part
    for i=1:Statesnum
        if REig1(i)>=0
           return;
        end
    end
%***********************************

%88888888888888888888888888888888888888888888
uu=[]; % saving the control signal for plot
x1c=[];
x11c=[];
Xpi=[];
Xpi1=[];
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
%8888888888888888888888888888888888888888888

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
    [t,x]= ode45('subF16',tspan,x0);
%     x(end,:)
%    pause
   x1=zeros(1,rowA);
   x1=x(length(x),1:rowA);
%    x1(1)=x(length(x),1);
%    x1(2)=x(length(x),8);
%    x1(3)=x(length(x),3);
%    x1(4)=x(length(x),10);
%    x1(5)=x(length(x),11);
   x11=zeros(1,rowA);
   x11=x(length(x),rowA+2:rowA+1+rowA);
   
   Y(j,:)=x(length(x),rowA+1);
    
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
          
  
    plot(t+T*(k-1),x(:,7),'b','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,8),'r','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,9),'g','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,10),'y','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,11),'black','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,1),':b','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,2),':r','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,3),':g','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,4),':y','LineWidth',2); hold on 
    plot(t+T*(k-1),x(:,5),':black','LineWidth',2); hold on 
    
    for i=1:rowA
        x0(i)=x(end,i);
    end

   for i=rowA+2:rowA+rowA+1
        x0(i)=x(end,i);
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
%            Xpi-Xpi1
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
     P
        K=inv(R)*(B'*P+L1)*C'*inv(C*C')
        damp(A-B*K*C)
%         pause
%         E=[E eig(A-B*K*C)];
         mm=mm+1;
%          L=R*K*C-B'*P;
%          Q=C'*Q*C+L'*inv(R)*L;
    end
end

%%

% figure(1); title('System states'); xlabel('Time(s)');legend('x(1)','x(2)','x(3)','x(4)');
figure(1); title('System states'); xlabel('Time(s)');
str1='$$\hat{\alpha}$$';
str2='$$\hat{\dot{q}}$$'; 
str3='$$\hat{\delta}_e$$';
str4='$$\hat{\alpha}_F$$'; 
str5='$$\hat{\epsilon}$$'; 
str6='$$\alpha$$';
str7='$$\dot{q}$$'; 
str8='$$\delta_e$$';
str9='$$\alpha_F$$'; 
str10='$$\epsilon$$'; 
legend({str1,str2,str3,str4,str5,str6,str7,str8,str9,str10},'interpreter','latex','fontsize',20)
K=inv(R)*(B'*P+L1)*C'*inv(C*C')
P
%%
K
K1=[ 0.0002 -0.1441 4.3917 31.6328]
Ac1=A-B*K1*C;
Ac=A-B*K*C;
damp(Ac)
% x01=x0(1:rowA);
 x01=[1 2 1 -1 1];
t=(0:0.05:7);
u=zeros(size(t));
[y,x]=lsim(Ac,B,C,D,u,t,x01);
[y1,x1]=lsim(Ac1,B,C,D,u,t,x01);

figure(2);
y
plot(t,y(:,1),'b','LineWidth',2);hold on
plot(t,y(:,2),'r','LineWidth',2);hold on
plot(t,y(:,3),'g','LineWidth',2);hold on
plot(t,y(:,4),'black','LineWidth',2);hold on
plot(t,y1(:,1),':b','LineWidth',2);hold on
plot(t,y1(:,2),':r','LineWidth',2);hold on
plot(t,y1(:,3),':g','LineWidth',2);hold on
plot(t,y1(:,4),':black','LineWidth',2);hold on
title('Outputs y'); xlabel('Time(s)');
str1='$$\alpha_F (Algotirhm 3)$$';
str2='$$q (Algotirhm 3)$$'; 
str3='$$e (Algotirhm 3)$$'; 
str4='$$ \epsilon (Algotirhm 3)$$';
str5='$$\alpha_F (Algotirhm 2)$$'; 
str6='$$q (Algotirhm 2)$$';
str7='$$e (Algotirhm 2)$$';
str8='$$\epsilon (Algotirhm 2)$$';
legend({str1,str2,str3,str4,str5,str6,str7,str8},'interpreter','latex','fontsize',20)
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
title('P matrix parameters'); xlabel('Time(s)');legend('P_1_1','P_1_2','P_1_3','P_1_4','P_1_5','P_2_2','P_2_3','P_2_4','P_2_5','P_3_3','P_3_4','P_3_5','P_4_4','P_4_5','P_5_5');

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

