 clc
close all
 clear all
global A B G C F D Q R r K
%******************************************************
%matrix coefficient
A=[0  1   0 0;...
   0 -11  0 0 ;...
  -1  0   0 0;...
   0  1/0.005 0 -1/0.005];
B=[0; 1; 0; 0];
C=[0 0 1 0;...
   -1 0 0 0;...
    0 0 0 1];
D=0;
F=[0 1 0]';
G=[0 0 1 0]';
r=[1];
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
%******************************************************
%******************************************************
a1=min(rowC,columnC);
if rank(C)<a1
    sprintf('The rank of C is not full row rank',2,3)
    return;
end
%******************************************************

%******************************************************
%Q R X
rou=1;
% the size of Q is (rowA,rowA)
Q=[10 0 0 0;...
   0 10 0 0;...
   0 0 105 0;...
   0 0 0 10];
R=rou*eye(Inputsnum);
X=eye(Statesnum);
%******************************************************

%******************************************************
L0=0.00*[1 1 1 1];%gave an initial value of L closed to 0
P0=are(A,B*inv(R)*B',Q+L0'*inv(R)*L0)%get the initial P0
K0=inv(R)*(B'*P0+L0)*C'*inv(C*C')%get the initial K0
L0=R*K0*C-B'*P0;%iteration L0
%***********************************
% Check if K0 stable, if not return
% damp(A-B*K0*C);
    Eig1=eig(A-B*K0*C);%calclate the eig of A-B*K0*C
    REig1=real(Eig1);%get the real part
    for i=1:Statesnum
        if REig1(i)>=0
           return;
        end
    end
%***********************************

for i=1:1000
    J0=abs(trace(P0*X));%%calculate the performance value
    P1=are(A,B*inv(R)*B',Q+L0'*inv(R)*L0);%calculate P for the new L0
    K1=inv(R)*(B'*P1+L0)*C'*inv(C*C');%calculate K for the new P
    L1=R*K1*C-B'*P1;%calculate the L for the new P,K

    J1=abs(trace(P1*X));%calculate the performance value
    %
    Eig2=eig(A-B*K1*C);%calclate the eig of A-B*K0*C
    REig2=real(Eig2);%get the real part
    if abs(K1-K0)<=0.0001
          P=P1;
          K=K1;
          L=L1;
          break
       else
         P0=P1;
         K0=K1;
         L0=L1;
       end
%     if REig2<0
%         if (J1<=J0)%if J1<=J0 then iterate the value of P,K,L
%            P0=P1;
%            K0=K1;
%            L0=L1;
%         else
%             break
%         end
%     end
%     
end
L
P
K
%  damp(A-B*K*C);
 Eig1=eig(A-B*K*C);%calclate the eig of A-B*K0*C
 REig1=real(Eig1);
%  figure(9)
%  plot(Eig1)
%figure
%*************************************************************
%track problem
x0=[0 10 0 0]; 
T=20;
tspan=[0 T];
[t,XXX]= ode45('cal',tspan,x0);
[M3,N3]=size(XXX);
YY=zeros(3,M3);
for i=1:M3
    %Y(:,i)=C*XX(i,:)'+F*r;
    YY(:,i)=C*XXX(i,:)'+F*r;%canclute the output state
end
%UU=zeros()

for i=1:M3
    UUU(i)=-K*(C*XXX(i,:)'-F*r);%canclute the input command
end

% %figure
% figure(1)
% %  plot(YYY(1,:),'b','LineWidth',2);hold on
%  plot(YY(1,:),'r','LineWidth',2);hold on
% % plot(Y(1,:),'g','LineWidth',2);
% xlabel('time response')
% ylabel('\epsilon')
% legend('jyovrabie','OPFB','H-infinite SOFC')
% title('error intergral \epsilon ')
% grid on
% 
% %hold on
% figure(2)
% % plot(YYY(2,:),'b','LineWidth',2);hold on
% plot(YY(2,:),'r','LineWidth',2);hold on
% % plot(Y(2,:),'g','LineWidth',2);
% 
% xlabel('time response')
% ylabel('e')
% 
% legend('jyovrabie','OPFB','H-infinite SOFC')
% title('error e')
% grid on
% 
% %hold on
% figure(3)
% % plot(YYY(3,:),'b','LineWidth',2);hold on
% plot(YY(3,:),'r','LineWidth',2);hold on
% % plot(Y(3,:),'g','LineWidth',2);
% xlabel('time response')
% ylabel('d(x_F)/dt')
% 
% legend('jyovrabie','OPFB','H-infinite SOFC')
% title('filter velocity d(x_F)/dt')
% grid on
% %hold on
% 
% figure(4)
% % plot(XXXX(:,1),'b','LineWidth',2);hold on
% plot(XXX(:,1),'r','LineWidth',2);hold on
% % plot(XX(:,1),'g','LineWidth',2);
% xlabel('time response')
% ylabel('x')
% legend('jyovrabie','OPFB','H-infinite SOFC')
% title('Position x')
% grid on
% 
% figure(5)
% % plot(XXXX(:,2),'b','LineWidth',2);hold on
% plot(XXX(:,2),'r','LineWidth',2);hold on
% % plot(XX(:,2),'g','LineWidth',2);
% xlabel('time response')
% ylabel('Postion velocity dx/dt')
% legend('jyovrabie','OPFB','H-infinite SOFC')
% title('Position velocity dx/dt')
% grid on
% 
% 
% figure(6)
% % plot(UUUU,'b','LineWidth',2); hold on
% plot(UUU,'r','LineWidth',2); hold on
% % plot(UU,'g','LineWidth',2);
% xlabel('time response')
% ylabel('u')
% legend('jyovrabie','OPFB','H-infinite SOFC')
% title('input u')
% grid on


%××××××××××××××××××××××××××××××××××××××××××××××××××××










