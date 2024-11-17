function data_init
% This function return all the data needed by the script to work
% NREL 5-MW Baseline Wind Turbine

    m_tower= 347460;
    m_t = m_tower + 240000;         % Mass of the tower and nacelle         [kg]
    m_b = 110000;                   % Mass of each blade                    [kg]
    gearbox_ratio=97;
    J_r = 115926+38759226;          % Inertia of the hub + rotor            [kg*m^2]
    J_g = 534.116*gearbox_ratio^2;  % Inertia of the generator              [Kg*m^2], we consider even the gear box ratio
    omega_n_t = 0.3240;             % Natural frequency of the tower        [rad/s]
    omega_n_b = 0.66;               % Natural frequency of each blade       [rad/s]
    omega_n_s = 0.62;               % Natural frequency of the transmission [rad/s]
    damp_ratio_t = 1e-2;            % Damping ratio of the tower            [rad/(s^2)]
    damp_ratio_b = 0.477465e-2;     % Damping ratio of the blade            [rad/(s^2)]
    damp_ratio_s = 0.5e-2;          % Damping ratio of the transmission     [rad/(s^2)]
    rated_rotor_speed=12.1*2*pi/60; % [rad/s]
    N_blades = 3;                   % Number of blades              []
    r_b = 61.5/2;                   % Point of application of zeta  [m]
    rotor_radius=63;         
    rho = 1.293;                    % Air density                   [kg/(m^3)]
    
    rated_gen_speed=123;            % [rad/s]
    rated_gen_T=43093.55;           % [Nm]
    synch_speed=1000*2*pi/60;       % [rad/s] 
    Bg=rated_gen_T/(rated_gen_speed/97-synch_speed/97);   %coefficient linear approximation torque-speed generator
    rated_P=5.3;                    % Scaled power

    % State space initial conditions
    y_t_0 = 0;
    zeta_0 = 0;
    theta_s_0 = 0;
    y_t_dot_0 = 0;
    zeta_dot_0 = 0;
    Omega_r_0 = 1.2671;
    Omega_g_0 = 1.2671;
    beta_0 = 20;

    % Stiffness and damaping coefficients
    stiff_t=3.6475e+04;     
    damp_t=2.2515e+03;
    stiff_b=4.7916e+04;
    damp_b=  693.2792;
    stiff_s = 867637e3;
    damp_s = 6215e3;
    
   
    % clear vars we do not want in the workspace
    clear omega_n_s omega_n_b omega_n_t damp_ratio_s damp_ratio_b damp_ratio_t
    
    save("data")

end