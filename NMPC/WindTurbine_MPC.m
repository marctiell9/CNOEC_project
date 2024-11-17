
%% WindTurbine_MPC.m
%   Authors:
%
% - Lupatini, Alessandro
% - Malusardi, Diego
% - Puchades Ibanez, Mar
% - Rositani Marco
clc
clear 
close all
%% Model parameters

% All numerical data can be found in 'data_init'
data_init
load("data.mat") 

M = [   m_t+N_blades*m_b   N_blades*m_b*r_b       0           0;
        N_blades*m_b*r_b   N_blades*m_b*(r_b^2)   0           0;
        0                  0                      J_r         0;
        0                  0                      0           J_g];

C = [   damp_t      0               0           0;
        0           damp_b*(r_b^2)  0           0;
        0           0               damp_s      -damp_s;
        0           0               -damp_s     damp_s];

L = [   eye(3)  [0  0   -1]'];

K_tilde = [ stiff_t     0               0;
            0           stiff_b*(r_b^2) 0;
            0           0               stiff_s;
            0           0               -stiff_s];

Q = [   N_blades       0       0;
        N_blades*r_b   0       0;
        0              1       0;
        0              0       -1];

A_tilde = [ zeros(3)            L              zeros(3,1);
            -(M^(-1))*K_tilde   -(M^(-1))*C    zeros(4,1);
            zeros(1,7)                        -1/1];

B_tilde = [ zeros(3,4);
            (M^(-1))*Q zeros(4,1);
            zeros(1,3) 1/1];

% Array of coefficients needed for the MatLab wind turbine model
fun_data = [r_b;rotor_radius;rho;Bg;synch_speed];

%% FHOCP parameters - single shooting
z0          =   [y_t_0;             % Initial state
                 zeta_0;
                 theta_s_0;
                 y_t_dot_0;
                 zeta_dot_0;
                 Omega_r_0;
                 Omega_g_0;
                 beta_0];  
Ts          =   0.005;               % Sampling time
Tend        =       1;               % Time horizon
N           =   Tend/Ts;             % Prediction steps
K           =   100;                 % Downsampling factor of the control variable
coll        =   0;                   % 1: collocation enabled, else: collocation not enabled
Q=eye(2*N);
Q(1:N,1:N)=1*eye(N);                % omega_g weight
Q(1*N+1:2*N,1*N+1:2*N)=1*eye(N);    % yt_dot  weight 
R=1e-4*eye(N);                      % beta_dot weight 
%% REFERENCE DEFINITIONS
% Defining the reference profile for the controlled variables
Omegarated=rated_gen_speed/gearbox_ratio;
Yref(1:N,1)=Omegarated;

%% WIND PROFILE GENERATION
% Generation of wind profile used in the simulation
T_wind       =   650; %length of the wind speed profile
pole_lpf=0.005; %pole of the low pass filter for the wind estimation, higher values leads to a better estimation
out=sim('wind_simulink_constr_1',"StopTime","T_wind");
V_real=out.V_real;
V_estimated=out.V_estimated;

%% Constraints 
% Inequality constraints imposed on the optimization variable beta
C       =       [-eye(N/K);
                eye(N/K)];

d       =       [-45*ones(N/K,1); 
                 0*ones(N/K,1)];




% Number of nonlinear inequality constraints
q       =     6*N; % inequality constraints are on omega_g, y_t, and beta_dot
%% Solution -  BFGS
% Initialize solver options
myoptions               =   myoptimset;
myoptions.Hessmethod  	=	'BFGS';
myoptions.gradmethod  	=	'CD'; 
myoptions.graddx        =	2^-17;
myoptions.tolgrad    	=	1e-9;
myoptions.ls_beta       =	0.8; 
myoptions.ls_c          =	0.05;
myoptions.ls_nitermax   =	20;
myoptions.tolfun    	=	1e-6;  
myoptions.nitermax      =	20;

%% Simulate with MPC
Tmpc=100;
Tmpc=ceil(Tmpc/(Ts*K))*Ts*K;
Nsim=Tmpc/(Ts*K);
zt                 =  z0;  
beta_d0            =  20*ones(N/K,1);
z_MPC=z0;

% We can introduce a flag to run the ideal and the real case
estim=0; % 1: real case, else: ideal case

tic
Uprec=20;
for ind=1:Nsim
    if estim == 1
        [Ustar,fxstar,k,exitflag,xsequence] = myfmincon(@(beta_d)WindTurbine_cost_constr(beta_d,zt,Yref,V_estimated(K*(ind-1)+1:N+K*(ind-1),1),Ts,Q,R,A_tilde,B_tilde,fun_data,K,Uprec,coll),beta_d0,[],[],C,d,0,q,myoptions);
    else
        [Ustar,fxstar,k,exitflag,xsequence] = myfmincon(@(beta_d)WindTurbine_cost_constr(beta_d,zt,Yref,V_real(K*(ind-1)+1:N+K*(ind-1),1),Ts,Q,R,A_tilde,B_tilde,fun_data,K,Uprec,coll),beta_d0,[],[],C,d,0,q,myoptions);
    end
    if coll == 1
        for i=1:N-1
            per = rem(i, K)/K;
            if floor(i/K)==0
                pre = Uprec;
            else
                pre = Ustar(floor(i/K),1);
            end
            post = Ustar(ceil(i/K),1);
            Ustar2(i,1) = (post-pre)*(3*per^2-2*per^3)+pre;
        end
    else
        for i=1:N-1
            Ustar2(i,1)=Ustar(floor((i-1)/K)+1,1);
        end
    end
    for i=1:K
        [zdot]               =   WindTurbine(z_MPC(:,K*(ind-1)+i),Ustar2(i,1),V_real(K*(ind-1)+i,1),A_tilde,B_tilde,fun_data);
        z_MPC(:,K*(ind-1)+i+1)         =   z_MPC(:,K*(ind-1)+i)+Ts*zdot;
        beta_d_opt(K*(ind-1)+i,1)=Ustar2(i,1);
    end
    zt=z_MPC(:,K*ind);     % initial state update
    beta_d0=[Ustar(2:N/K,1);Ustar(end,1)]; %FHOCP warm starting, initial guess update
    Uprec=Ustar(1,1);
    fprintf("Time: "+ind*K*Ts+"s\n");
end
toc

%% FIGURES

figure()
plot(0:Ts:(K*Nsim-1)*Ts,beta_d_opt)
grid
title('beta_desired')
hold on
plot(0:Ts:K*Nsim*Ts,z_MPC(8,:),'LineWidth',3)
legend('beta_d','beta'),xlabel('Time [s]'),ylabel('Angle degrees [°]')

figure()
plot(0:Ts:K*Nsim*Ts,V_estimated(1:K*Nsim+1,1))
grid
title('wind speed')
hold on
plot(0:Ts:K*Nsim*Ts,V_real(1:K*Nsim+1,1))
legend('estimated','real'),xlabel('Time [s]'),ylabel('Velocity [m/s]')

figure(),plot(0:Ts:K*Nsim*Ts,z_MPC(7,:)),grid,title('omega_g'),ylabel('Angular velocity [rad/s]'),xlabel('Time [s]')
figure(),plot(0:Ts:K*Nsim*Ts,Bg*(z_MPC(7,:)-synch_speed/97*ones(7,K*Nsim+1))),grid,title('Tg'),xlabel('Time [s]'),ylabel('Torque [Nm]')
figure(), plot(0:Ts:K*Nsim*Ts,z_MPC(1,:),'LineWidth',3), grid, title('yt'),xlabel('Time [s]'),ylabel('Displacement [m]')
figure() ,plot(0:Ts:K*Nsim*Ts,z_MPC(4,:)) ,grid ,title('yt_dot'),xlabel('Time [s]'),ylabel('Velocity [m/s]')
figure(),plot(0:Ts:(K*Nsim-2)*Ts,z_MPC(8,2:K*Nsim)/Ts-z_MPC(8,1:K*Nsim-1)/Ts),grid,title('beta_derivative'),xlabel('Time [s]'),ylabel('Pitch angle rate[°/s]')










