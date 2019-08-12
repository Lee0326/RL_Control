clc
clear all
close all
global A L B C Ac yhut y k n m l xxx

%******************************************************
%matrix coefficient
A=[-1.01887 0.90506   -0.00215   0    0;...
    0.82225  -1.07741 -0.17555   0    0;...
    0        0        -20.2      0    0;...
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
%******************************************************
%******************************************************
Statesnum=columnA;
Inputsnum=columnB;
Outputsnum=rowC;
unknownnum=rowA*(rowA+1)/2;
%******************************************************
%******************************************************
a1=min(rowC,columnC);
if rank(C)<a1
    sprintf('The rank of C is not full row rank',2,3)
    return;
end
%******************************************************

%initial condition
%The time 
n=Statesnum%the number of states 
m=Inputsnum%The number of the input
l=Outputsnum%the number of the output
k=n%the number of the hidden layer nurons


%The initial weighs of W and V
Whut0=0.02*rand(n,k);
Vhut0=0.00*rand(k,k);


%The initial value
x0=[2,-1,1,-2,1]';
xhut0=[2,-1,1,-2,1]';
xhutknown=[2,-1,1,-2,1]';

Whut00=zeros(n*k,1);
for i=1:k
    Whut00(n*(i-1)+1:n*i)=Whut0(1:n,i);
end
Vhut00=zeros(k*k,1);
for i=1:k
    Vhut00(k*(i-1)+1:k*i)=Vhut0(1:k,i);
end


xxhutwv0=[x0;xhut0;Whut00;Vhut00;xhutknown];
%size(xxhutwv0)


%Slove the equation
T=20;%The time 
tspan=[0,T]; 
[t,xxhutwv]=ode45(@observerF16,tspan,xxhutwv0);
y=C*xxhutwv(:,1:n)';
x11=xxhutwv(:,1:n)'-xxhutwv(:,2*n+n*k+k*k+1:2*n+n*k+k*k+n)';
yhut=C*xxhutwv(:,n+1:n+n)';
yhutknown=C*xxhutwv(:,2*n+n*k+k*k+1:2*n+n*k+k*k+n)';
error=C*x11;
%plot
figure(1)
for i=1:n
    %figure(i)
    plot(t,xxhutwv(:,1),'r'); hold on
    plot(t,xxhutwv(:,n+1),'b'); hold on
%     plot(t,xxhutwv(:,2*n+n*k+k*k+i),'g'); hold on
    xlabel('time response')
    ylabel('states')
    legend('x','xhut','xhutknown')
    title('states compare')
end


figure(n+1)
for i=2*n+1:2*n+n*k    
    plot(t,xxhutwv(:,i),'r','LineWidth',i*0.01); hold on
end
    xlabel('time response')
    ylabel('Value W')
    title('weights W')
 
figure(n+2)
for i=2*n+n*k+1:2*n+n*k+k*k    
    plot(t,xxhutwv(:,i),'r','LineWidth',i*0.01); hold on
end
    xlabel('time response')
    ylabel('Value V ')
    title('weights V')
    
 figure(n+3)
    plot(t,y,'r'); hold on
    plot(t,yhut,'b'); hold on
%     plot(t,yhutknown,'g'); hold on
    xlabel('time response')
    ylabel('outputs')
    legend('y(1)','y(2)','y(3)','yhut(1)','yhut(2)','yhut(3)','yhutknown(1)','yhutknown(2)','yhutknown(3)')
    title('ouputs compare')

    
 figure(n+4)
%     plot(t,error,'r'); hold on
     plot(t,x11,'b'); hold on
    xlabel('time response')
    ylabel('error')
%     legend('y(1)','y(2)','y(3)')
    title('error')
    






