clc
clear all
close all
global A L C Ac yhut y k n m l
%initial condition
%The time 
n=4%the number of states 
m=1%The number of the input
l=3%the number of the output
k=5%the number of the hidden layer nurons


%The initial weighs of W and V
Whut0=0.2*rand(n,k);


%The initial value
x0=[1 -1 1 -1]';
xhut0=[0 0 0 -0.8]';

Whut00=zeros(n*k,1);
for i=1:k
    Whut00(n*(i-1)+1:n*i)=Whut0(1:n,i);
end


xxhutwv0=[x0;xhut0;Whut00];
%size(xxhutwv0)


%Slove the equation
T=200;%The time 
tspan=[0,T]; 
[t,xxhutwv]=ode45(@slovelinear,tspan,xxhutwv0);

y=C*xxhutwv(:,1:n)';

yhut=C*xxhutwv(:,n+1:n+n)';

%plot
figure(1)
for i=1:3
   
    %figure(i)
    plot(t,xxhutwv(:,i),'r'); hold on
    plot(t,xxhutwv(:,n+i),'b'); hold on
    xlabel('time response')
    ylabel('states')
%     legend('x','xhut')
    title('states compare')
end
legend('x(1)','x(2)','x(3)','x(4)','x(5)','x(6)','x(7)','x(8)')


figure(n+1)
for i=2*n+1:2*n+n*k    
    plot(t,xxhutwv(:,i),'r','LineWidth',i*0.01); hold on
end
    xlabel('time response')
    ylabel('Value W')
    title('weights W')

    
 figure(n+3)
    plot(t,y,'r'); hold on
    plot(t,yhut,'b'); hold on
    xlabel('time response')
    ylabel('outputs')
    legend('y','yhut')
    title('ouputs compare')

    

    






