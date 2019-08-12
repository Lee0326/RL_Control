clear all;clc; close all
S=[ 0.5100    0.0057   -0.0038   -0.0065         0         0   -0.0064   -0.0072
    0.0057    0.5032   -0.0067   -0.0112         0         0   -0.0074   -0.0082
         0         0    0.5015    0.0025   -0.0042   -0.0047         0         0
         0         0    0.0025    0.5042   -0.0046   -0.0051         0         0
   -0.0109   -0.0062         0         0    0.5028    0.0031         0         0
   -0.0061   -0.0035         0         0    0.0031    0.5034         0         0
         0         0         0         0         0         0    0.5047    0.0052
         0         0         0         0         0         0    0.0052    0.5058];
%% System
a=[0.995 0.09983;-0.09983 0.995];
b1=[0.2047;0.08984];b2=[0.2147;0.2895];b3=[0.2097;0.1897];b4=[0.2;0.1];
g1=0;g2=0;g3=0;g4=1;
R11=1;R22=1;R33=1;R44=1;
e12=0.8;   e14=0.7;    e13=0;    e21=0;     e23=0.6;   e31=0.8;    e32=0;               e41=0;
R12=1;     R14=1;      R13=0;    R21=0;     R23=1;     R31=1;      R32=0;               R41=0;
Q11=eye(2);Q22=eye(2);Q33=eye(2);Q44=eye(2);
EE=[0   e12 e13 e14
    e21 0   e23 0
    e31 e32 0   0
    e41 0   0   0 ];
di=EE*[1;1;1;1];
d1=di(1);d2=di(2);d3=di(3);d4=di(4);
Di=blkdiag(d1,d2,d3,d4);
LG=Di-EE+blkdiag(g1,g2,g3,g4);
wc1(:,:,1)=zeros(2,6);wa1(:,:,1)=rand(1,6);
wc2(:,:,1)=zeros(2,4);wa2(:,:,1)=rand(1,4);
wc3(:,:,1)=zeros(2,4);wa3(:,:,1)=rand(1,4);
wc4(:,:,1)=zeros(2,2);wa4(:,:,1)=rand(1,2);
ua=0.1;uc=0.1;
A=blkdiag(a,a,a,a);
B=blkdiag(b1,b2,b3,b4);
R=blkdiag(R11,R22,R33,R44);
g=blkdiag((g1+d1),(g2+d2),(g3+d3),(g4+d4));
E=[(g1+d1)*eye(2) -1*e12*eye(2) -1*e13*eye(2) -1*e14*eye(2);-1*e21*eye(2) (g2+d2)*eye(2) -1*e23*eye(2) zeros(2,2);-1*e31*eye(2) -1*e32*eye(2) (g3+d3)*eye(2) zeros(2,2);-1*e41*eye(2) zeros(2,2) zeros(2,2) (g4+d4)*eye(2)];
%% global
    k=1;l=1;

    u1(l)=rand;u2(l)=rand;u3(l)=rand;u4(l)=rand;

    x(:,1)=rand(8,1);x0(:,1)=rand(2,1);
    x0bar=kron([1;1;1;1],eye(2))*x0(:,1);Z(:,1)=-(kron(LG,eye(2)))*(x-x0bar); 
    x1(:,k)=Z(1:2,1);x2(:,k)=Z(3:4,1);x3(:,k)=Z(5:6,1);x4(:,k)=Z(7:8,1);
while (l <= 4000)
    Z1(:,k)=[x1(:,k);x2(:,k);x4(:,k)];
    Z2(:,k)=[x2(:,k);x3(:,k)];
    Z3(:,k)=[x1(:,k);x3(:,k)];
    Z4(:,k)=[x4(:,k)];
    u1(l)=wa1(:,:,l)*Z1(:,k);u2(l)=wa2(:,:,l)*Z2(:,k);u3(l)=wa3(:,:,l)*Z3(:,k);u4(l)=wa4(:,:,l)*Z4(:,k);
    x1(:,k+1)=(a*x1(:,k))-((g1+d1)*b1*u1(l))+(e12*b2*u2(l))+(e13*b3*u3(l))+(e14*b4*u4(l));
    x2(:,k+1)=(a*x2(:,k))-((g2+d2)*b2*u2(l))+(e21*b1*u1(l))+(e23*b3*u3(l));
    x3(:,k+1)=(a*x3(:,k))-((g3+d3)*b3*u3(l))+(e31*b1*u1(l))+(e32*b2*u2(l));
    x4(:,k+1)=(a*x4(:,k))-((g4+d4)*b4*u4(l))+(e41*b1*u1(l));
    
    ZZ1(:,k)=[2*x1(:,k+1);x2(:,k+1);x4(:,k+1)];
    ZZ2(:,k)=[2*x2(:,k+1);x3(:,k+1)];
    ZZ3(:,k)=[x1(:,k+1);2*x3(:,k+1)];
    ZZ4(:,k)=[2*x4(:,k+1)];
    
    SS(:,:,k)=[wc1(:,1:4,k) zeros(2,2) wc1(:,5:6,k);zeros(2,2) wc2(:,:,k) zeros(2,2);wc3(:,1:2,k) zeros(2,2) wc3(:,3:4,k) zeros(2,2);zeros(2,6) wc4(:,:,k)];
    
    x0bar(:,l)=kron([1;1;1;1],eye(2))*x0(:,l);
    x0(:,l+1)=a*x0(:,l);
    Z(:,l)=[x1(:,k);x2(:,k);x3(:,k);x4(:,k)]; % varepsilon updat
    XZ(:,l)=x0bar(:,l)+inv(-(kron(LG,eye(2))))*Z(:,l);
    
    lym1(:,:,l)=wc1(:,:,l)*ZZ1(:,k);lym2(:,:,l)=wc2(:,:,l)*ZZ2(:,k);lym3(:,:,l)=wc3(:,:,l)*ZZ3(:,k);lym4(:,:,l)=wc4(:,:,l)*ZZ4(:,k);
    C1=(g1+d1)*(R11^-1)*b1'*lym1(:,:,l);
    C2=(g2+d2)*(R22^-1)*b2'*lym2(:,:,l);
    C3=(g3+d3)*(R33^-1)*b3'*lym3(:,:,l);
    C4=(g4+d4)*(R44^-1)*b4'*lym4(:,:,l);
    
    cu1=wa1(:,:,l)*ZZ1(:,k);
    cu2=wa2(:,:,l)*ZZ2(:,k);
    cu3=wa3(:,:,l)*ZZ3(:,k);
    cu4=wa4(:,:,l)*ZZ4(:,k);

    V1=(0.5*x1(:,k)'*Q11*x1(:,k))+(0.5*(cu1'*R11*cu1+cu2'*R12*cu2+cu4'*R14*cu4))+(x1(:,k+1)'*wc1(:,:,l)*ZZ1(:,k));
    V2=(0.5*x2(:,k)'*Q22*x2(:,k))+(0.5*(cu2'*R22*cu2+cu3'*R23*cu3))+(x2(:,k+1)'*wc2(:,:,l)*ZZ2(:,k));
    V3=(0.5*x3(:,k)'*Q33*x3(:,k))+(0.5*(cu3'*R33*cu3+cu1'*R31*cu1))+(x3(:,k+1)'*wc3(:,:,l)*ZZ3(:,k));
    V4=(0.5*x4(:,k)'*Q44*x4(:,k))+(0.5*(cu4'*R44*cu4))+(x4(:,k+1)'*wc4(:,:,l)*ZZ4(:,k));

    wa1(:,:,l+1)=wa1(:,:,l)-((ua*Z1(:,k)')*((wa1(:,:,l)*Z1(:,k))-(C1)));
    wa2(:,:,l+1)=wa2(:,:,l)-((ua*Z2(:,k)')*((wa2(:,:,l)*Z2(:,k))-(C2)));
    wa3(:,:,l+1)=wa3(:,:,l)-((ua*Z3(:,k)')*((wa3(:,:,l)*Z3(:,k))-(C3)));
    wa4(:,:,l+1)=wa4(:,:,l)-((ua*Z4(:,k)')*((wa4(:,:,l)*Z4(:,k))-(C4)));

    wc1(:,:,l+1)=wc1(:,:,l)-(uc*(-V1+((wc1(:,:,l)*Z1(:,k))'*x1(:,k)))*x1(:,k)*Z1(:,k)');
    wc2(:,:,l+1)=wc2(:,:,l)-(uc*(-V2+((wc2(:,:,l)*Z2(:,k))'*x2(:,k)))*x2(:,k)*Z2(:,k)');
    wc3(:,:,l+1)=wc3(:,:,l)-(uc*(-V3+((wc3(:,:,l)*Z3(:,k))'*x3(:,k)))*x3(:,k)*Z3(:,k)');
    wc4(:,:,l+1)=wc4(:,:,l)-(uc*(-V4+((wc4(:,:,l)*Z4(:,k))'*x4(:,k)))*x4(:,k)*Z4(:,k)');
        
    if (l >= 1500) && (((norm(wc1(:,:,l)-wc1(:,:,l-1))+norm(wc2(:,:,l)-wc2(:,:,l-1))+norm(wc3(:,:,l)-wc3(:,:,l-1))+norm(wc4(:,:,l)-wc4(:,:,l-1)))/4)<= 0.000001)
    break;
    end
    l=l+1;k=k+1;
end

for i=1:500
    m1(i)=wc1(1,1,i);
    m2(i)=wc1(1,2,i);
    m3(i)=wc1(1,3,i);
    m4(i)=wc1(1,4,i);
    m5(i)=wc1(1,5,i);
    m6(i)=wc1(1,6,i);
    m9(i)=wc1(2,1,i);
    m10(i)=wc1(2,2,i);
    m11(i)=wc1(2,3,i);
    m12(i)=wc1(2,4,i);
    m13(i)=wc1(2,5,i);
    m14(i)=wc1(2,6,i);
end
figure
plot(1:i,m1,1:i,m2,1:i,m3,1:i,m4,1:i,m5,1:i,m6,1:i,m9,1:i,m10,1:i,m11,1:i,m12,1:i,m13,1:i,m14)

xlabel('Iteration Steps')
Ylabel ('Critic Weights')
for i=1:500
    mm1(i)=wa1(1,1,i);
    mm2(i)=wa1(1,2,i);
    mm3(i)=wa1(1,3,i);
    mm4(i)=wa1(1,4,i);
    mm5(i)=wa1(1,5,i);
    mm6(i)=wa1(1,6,i);
end
figure
plot(1:i,mm1,1:i,mm2,1:i,mm3,1:i,mm4,1:i,mm5,1:i,mm6)
%title('Actor Update for agent 1' )
xlabel('Iteration Steps')
Ylabel ('Actor Weights')
figure
plot(1:1500,Z(:,1:1500)')
xlabel('Iteration Steps')
ylabel('Tracking Error Dynamics')
figure
plot(XZ(1,:),XZ(2,:),XZ(3,:),XZ(4,:),XZ(5,:),XZ(6,:),XZ(7,:),XZ(8,:))
xlabel('Phase 1')
ylabel('Phase 2')
figure
plot(1:1500,XZ(:,1:1500)')
xlabel('Iteratin Steps')
ylabel('Dynamics')
K=inv(eye(8)+(E*B*g*inv(R)*B'*S))*A; 
Zm(:,1)=[x1(:,1);x2(:,1);x3(:,1);x4(:,1)];
for l=1:l
  x0bar(:,l)=kron([1;1;1;1],eye(2))*x0(:,l);
  Xm(:,l)=x0bar(:,l)+inv(-(kron(LG,eye(2))))*Zm(:,l);
  x0(:,l+1)=a*x0(:,l);
  Zm(:,l+1)=K*Zm(:,l); 
end
for j=1:1500
vric(j)=Zm(:,j)'*S*Zm(:,j);
vcrit(j)=Z(:,j)'*SS(:,:,j)*Z(:,j);
end
figure
plot(1:j,vcrit,1:j,vric,'-.')
xlabel('Time index')
ylabel('The Value Function')
%title('Value of Initial conditions by Riccati and Neural Netwroks')
legend('NN','Riccati')
figure
plot(1:1500,XZ(1,1:1500)',1:1500,XZ(2,1:1500)',1:1500,Xm(1,1:1500)','-.',1:1500,Xm(2,1:1500)','-.')
xlabel('Iteratin Steps')
ylabel('Dynamics of agent 1')
