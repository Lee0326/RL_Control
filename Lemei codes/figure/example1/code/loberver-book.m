%linear observer
clc
clear all
close all

global A B C

A=[-4 1;-4 0];
B=[1 3]';
C=[1 2];

xy0=[2 -1 0 0 1 1 1 1 1 1 1 1 2 -1]';

T=20;%The time 
tspan=[0,T]; 
[t,xy]=ode45(@loberver,tspan,xy0);

y=C*xy(:,1:2)';
yhut=C*xy(:,3:4)';
yhutknown=C*xy(:,13:14)';

plot
figure(4)
for i=1:2
    %figure(i)
    plot(t,xy(:,6),'r'); hold on
    plot(t,xy(:,7),'b'); hold on
    plot(t,xy(:,8),'g');hold on
    plot(t,xy(:,9),'g');hold on
    xlabel('time response')
    ylabel('A')
    legend('A11','A12','A21','A22')
    title('A compare')
end

figure(1)
for i=1:1
    %figure(i)
    plot(t,xy(:,i),'r'); hold on
    plot(t,xy(:,2+i),'b'); hold on
    plot(t,xy(:,12+i),'g');hold on
    xlabel('time response')
    ylabel('states')
    legend('x','xhut','xhutknown')
    title('states compare')
end
figure(1+2)

    %figure(i)
    plot(t,xy(:,2),'r'); hold on
    plot(t,xy(:,2+2),'b'); hold on
    plot(t,xy(:,12+2),'g');hold on
    xlabel('time response')
    ylabel('states')
    legend('x','xhut','xhutknown')


    
 figure(2)
    plot(t,y,'r'); hold on
    plot(t,yhut,'b'); hold on
    plot(t,yhutknown,'g');hold on
    xlabel('time response')
    ylabel('outputs')
    legend('y','yhut','yhutknown')
    title('ouputs compare')