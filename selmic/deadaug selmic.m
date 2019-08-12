    %  File DEADAUG.M
    %   
    %  Simulation of Deadzone Compensation Using Augmeted NN.

    %  Simulation for paper: Deadzone Compensation in Motion Control
    %  Systems Using Neural Networks.
    %  Authors: R. Selmic and F. L. Lewis

    %  Author: Rastko Selmic, Nov. 97

    %  Definition of global variables, coeffictions in simulation, etc.     

    
       clear all;
    
    %  Definition of Global Variables

       global N1 N2 N3 N1i N2i N3i Nai n Kv;
       global V v0 Vi v0i;
       

    %  PARAMETERS OF THE NN I (Neural Network used as estimator)
       N1 = 2;         %  Input layer
       N2 = 25;        %  Hidden layer
       N3 = 2;         %  Output layer
       
    %  Generate the random values for V
    %  random numbers are uniformly distributed between v_min and v_max.
       V_min = -0.001; V_max = 0.001;
       v0_min = -35; v0_max = 35;
              
       V = rand(N1, N2);                % first layer are random numbers. 
       v0 = rand(1, N2);                % thresholds are random numbers. 
       V = (V.*(V_max - V_min))+V_min;
       v0 = (v0.*(v0_max - v0_min))+v0_min;
       
       
    %  PARAMETERS OF THE NN II (Neural Network used as compensator)
       N1i = 2;        %  Input layer
       N2i = 20;       %  Hidden layer
       N3i = 2;        %  Output layer
       Nai = 4;        %  Jump layer 

       Vi = rand(N1i, N2i);             % first layer are random numbers. 
       v0i = rand(1, N2i);              % thresholds are random numbers. 
       Vi = (Vi.*(V_max - V_min))+V_min;
       v0i =(v0i.*(v0_max - v0_min))+v0_min;
       
       
       n  = 4;         %  # of equations of the system
       
       %  x1 = q1;
    %  x2 = q2;
    %  x3 = q1dot;
    %  x4 = q2dot; 

    %  Initialize states
       x0(1:(n + N3*(N2+1) + N3i*(N2i+1) + N3i*(Nai+1)),1)...        
          =zeros((n + N3*(N2+1) + N3i*(N2i+1) + N3i*(Nai+1)),1);  %system states and NN I & NN II
        
    %  Generate random initial values for the matrix WT
    %  random numbers are uniformly distributed between w_min and w_max.
       w_min = -2; w_max = 2;
       for i = 1:N3,
           for j = 1:(N2+1),
               x0(n + (i-1)*(N2+1) + j) = rand(1)*(w_max - w_min)+w_min;
           end
       end
       
%       save data.mat V v0 Vi v0i x0;    % what ever is randomly selected, save it.

       load data.mat;

    %  Velocity Following Control Gain

    %  Kv is given as gobal varable;
    
       Kv = 40;
    
       
    %  Definition of parameters for simulation

       T0 = 0;         %  Initial time
       Tf = 10;        %  Final time
       
    %  Integrate using the "DYNAMICS.M" file      
	
       Tolerance = 0.005;
                   
    %  This part of the code speeds up the integration of the 
    %  differential equation for long simulations. ODE23 is now
    %  called with short vectors x and t instead with very long variables
    
       T = 0.5;
       while T0 < Tf
           [t,x] = ode23('dynamics', T0, T0+T, x0, Tolerance);
           Xall = [Xall; x];
           Tall = [Tall; t];
           save temp.mat Xall, Tall;
           T0 = T0 + T;
           x0 = x(length(x), :);
           clear x, t;
       end;
       
       x = Xall;
       t = Tall;


%       figure;
%       plot (t, x(:,1)-sin(t));
%       title ('Position error for the 1. joint');
       
%       figure;
%       plot (x(:,1)-sin(t), x(:,3)-cos(t));
%       title ('Error space for the 1. joint');
       
%       figure;
%       plot (t, x(:,2)-cos(t));
%       title ('Position error for the 2. joint');

       load pd.mat     % load values for PD controller for comparison
                     
       figure;
       plot (t, x(:,1)-1);
       hold;
       plot (pdt, pd(:,1)-1, 'r');
       title ('Position error for the 1. joint');
       
%       figure;
%       plot (t, x(:,3));
%       hold;
%       plot (pdt, pd(:,3), 'r');
%       title ('Velocity error for the 1. joint');
       
       figure;
       plot (x(:,1)-1, x(:,3));
       hold;
       plot (pd(:,1)-1, pd(:,3), 'r');
       title ('Error space for the 1. joint');
       
       figure;
       plot (t, x(:,2)-1);
       hold;
       plot (pdt, pd(:,2)-1, 'r');
       title ('Position error for the 2. joint');
       
 %      figure;
 %      plot (t, x(:,4));
 %      hold;
 %      plot (pdt, pd(:,4), 'r');
 %      title ('Velocity error for the 2. joint');
       
       figure;
       plot (x(:,2)-1, x(:,4));
       hold;
       plot (pd(:,2)-1, pd(:,4), 'r');
       title ('Error space for the 2. joint');
       
       save temp.mat
       
      
       TRALALA!!!!!!
