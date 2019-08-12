  % File dynamics.m
  %   
  % Simulation for paper: Deadzone Compensation in Motion Control
  % Systems Using Neural Networks.
  % Authors: R. R. Selmic and F. L. Lewis

  % Author: Rastko Selmic, Nov. 97

  % This file contains 2-link Robot Arm & NEURAL NETWORK DYNAMICS


    function xp = dynamics(t,x)

    global   N1 N2 N3 N1i N2i N3i Nai n Kv;
    global   V v0 Vi v0i;
    
    Lambda = 7;

  
  % COMPUTATION OF NONLINEAR CONTROL LAW

  % Auxiliary Control Input

%    Vc(1,1) = sin(t);           % Desired trajectory for q1;
%    Vc(2,1) = cos(t);           % Desired trajectory for q2; 
%    Vc(3,1) = cos(t);           % Desired q1dot;
%    Vc(4,1) = -sin(t);          % Desired q2dot;
%    Vc(5,1) = -sin(t);          % Desired q1doubledot;
%    Vc(6,1) = -cos(t);          % Desired q2doubledot;


    Vc(1,1) = 1;                % Desired trajectory for q1;
    Vc(2,1) = 1;                % Desired trajectory for q2; 
    Vc(3,1) = 0;                % Desired q1dot;
    Vc(4,1) = 0;                % Desired q2dot;
    Vc(5,1) = 0;                % Desired q1doubledot;
    Vc(6,1) = 0;                % Desired q2doubledot;


    %  Pulse train is the input

    % if (t>=0)&(t<5) 
    % Vc(1,1) = 1;                % Desired trajectory for q1;
    % Vc(2,1) = 1;                % Desired trajectory for q2; 

    % elseif (t>=5)&(t<10) 
    % Vc(1,1) = 0;                % Desired trajectory for q1;
    % Vc(2,1) = 0;                % Desired trajectory for q2; 

    % elseif (t>=10)&(t<15)
    % Vc(1,1) = 1;                % Desired trajectory for q1;
    % Vc(2,1) = 1;                % Desired trajectory for q2;   
    % else
    % Vc(1,1) = 0;                % Desired trajectory for q1;
    % Vc(2,1) = 0;                % Desired trajectory for q2; 
    % end;

    % Vc(3,1) = 0;                % Desired q1dot;
    % Vc(4,1) = 0;                % Desired q2dot;
    % Vc(5,1) = 0;                % Desired q1doubledot;
    % Vc(6,1) = 0;                % Desired q2doubledot;
 

  % CALCULATION OF THE CONTROL LAW
    for i = 1:n,
       e(i,1) = Vc(i) - x(i);
    end
    for i = 1:N3,
       r(i,1) = Lambda*e(i) + e(i+N3);
    end
    norm_r = norm(r);               % Finding the two-norm of r;


  % SIGNAL w
    w = Kv*r;
    
  % NN II  
  % NN II Input Vector (NN used as compensator of the deadzone)
    Y(1,1) = 1;
    Y(2,1) = w(1);
    Y(3,1) = w(2);

  % Computation of sigmoid function
    aux = [v0i; Vi]'*Y;
    sigmai = 1./(1 + exp(-aux));
    sigmai = [1; sigmai]; 
    p = zeros (1, N2i);

  % Generation of the weight matrix W1Ti
    for i = 1:N3i,
        for j = 1:(N2i+1),
            W1Ti(i,j) = x(n + N3*(N2+1) + (i-1)*(N2i+1) + j);
        end
    end
    W1i = W1Ti';
    
  % Generation of the weight matrix W2Ti
    for i = 1:N3i,
        for j = 1:(Nai+1),
            W2Ti(i,j) = x(n + N3*(N2+1) + N3i*(N2i+1) + (i-1)*(Nai+1) + j);
        end
    end
    W2i = W2Ti';

  % Generation of vector function fi:
    if w(1)>0 
       fi1 = 1;
       fi2 = 1-exp(-w(1));
    else
       fi1 = 0;
       fi2 = 0; 
    end

    if w(2)>0 
       fi3 = 1;
       fi4 = 1-exp(-w(2));
    else
       fi3 = 0;
       fi4 = 0;
    end
    fi = [1; fi1; fi2; fi3; fi4];
    
    Wi = [W1Ti, W2Ti]';
    sigmai = [sigmai; fi];
    
  % NN II output
  
    NNII = Wi'*sigmai;

  % Deadzone compensation  
  
    u = w + NNII;
    
%    u = w;
  
   
  % NN I Input Vector (NN used as estimator of the deadzone)
    X(1,1) = 1;
    X(2,1) = u(1);
    X(3,1) = u(2);

  % Computation of sigmoid function and derivativ of sigmoid
    aux = [v0; V]'*X;
    sigma = 1./(1 + exp(-aux));
    sigmap = diag(sigma.*(1 - sigma));
    sigma = [1; sigma]; 
    p = zeros (1, N2);
    sigmap = [p; sigmap];

  % Generation of the weight matrix WT
    for i = 1:N3,
      for j = 1:(N2+1),
        WT(i,j) = x(n + (i-1)*(N2+1) + j);
      end
    end
    W = WT';

  % NN I output
    NNI = WT*sigma;
    

  % DEADZONE IN THE MECHANICAL SYSTEM
  % deadzone parameters:
    d_plus = 40; d_minus = -30;
    m_plus = 1; m_minus = 1;
    
    if u(1) >= d_plus
        tau(1) = m_plus*sqrt(u(1)-d_plus);
    elseif u(1) <= d_minus;
        tau(1) = m_minus*(-sqrt(u(1)-d_minus));
    else
        tau(1) = 0;
    end
     
    if u(2) >= d_plus
        tau(2) = m_plus*sqrt(u(2)-d_plus);
    elseif u(2) <= d_minus;
        tau(2) = m_minus*(-sqrt(u(2)-d_minus));
    else
        tau(2) = 0;
    end
       
    
  % 2-LINK ARM DYNAMICS

    m1=0.8; m2=2.3; le1=1; le2=1; g=9.8; 

    %  if t >= 5
    %      m2 = 4;
    %  end

    m11=(m1+m2)*le1*le1 + m2*le2*le2 + 2.*m2*le1*le2*cos(x(2));
    m22=m2*le2*le2;
    m12=m22 + m2*le1*le2*cos(x(2));

    v1=-m2*le1*le2*(2.*x(3)*x(4) + x(4)*x(4))*sin(x(2));
    v2=m2*le1*le2*x(3)*x(3)*sin(x(2));

    g1=(m1+m2)*g*le1*cos(x(1)) + m2*g*le2*cos(x(1) + x(2));
    g2=m2*g*le2*cos(x(1) + x(2));
    
    det=m11*m22 - m12*m12;
    
    a1 = tau(1) - v1 -g1;
    a2 = tau(2) - v2 -g2;

    
  % STATE EQUATIONS
  % System equations (n)

    xp(1) = x(3);
    xp(2) = x(4);
    xp(3) = ( m22*a1 - m12*a2 )/det;
    xp(4) = ( m11*a2 - m12*a1 )/det;


  % NN tuning rules
  % Design Parameters
    S = 25; T = 60000000;
    k1 = 0.0000001; k2 = 0.00000005;
    
       
    WT_prime = (-S*sigmap*V'*Wi'*sigmai*r' - k1*S*norm_r*W)';    % Dif. eq. for W;
    for i = 1:N3,
        for j = 1:(N2+1),
            xp(n + (i-1)*(N2+1) + j) = WT_prime(i,j);
        end;
    end;

    WTi_prime = (T*sigmai*r'*WT*sigmap*V' - k1*T*norm_r*Wi...
                 - k2*T*norm_r*(trace(Wi'*Wi))*Wi)';             % Dif. eq. for Wi;
    for i = 1:N3i,
        for j = 1:(N2i+1),
            xp(n + N3*(N2+1) + (i-1)*(N2i+1) + j) = WTi_prime(i,j);
        end;
    end;
    for i = 1:N3i,
        for j = 1:(Nai+1),
            xp(n + N3*(N2+1) + N3i*(N2i+1) + (i-1)*(Nai+1) + j) = WTi_prime(i,j);
        end;
    end;
    
    