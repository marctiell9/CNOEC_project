function f = WindTurbine_cost_constr(beta_d,z0,Yref,V,Ts,Q,R,A_tilde,B_tilde,fun_data,K,Uprec,coll)
%% Initialize simulation output 
N                   =    length(Yref);      % number of samples 
z_FFD               =       zeros(8,N);     % matrix with states
z_FFD(1:8,1)        =       z0;
F = zeros(3*N,1);

%% Downsamplng of control variables and collocation
if coll ~= 1
    beta_d1=zeros(N,1);
    for i=1:N-1
        beta_d1(i,1)=beta_d(floor((i-1)/K)+1,1);
    end
else
    if K > 1
        for i=1:N-1
            per = rem(i, K)/K;
            if floor(i/K) == 0
                pre = Uprec;
            else
                pre = beta_d(floor(i/K),1);
            end
            post = beta_d(ceil(i/K),1);
            beta_d1(i,1) = (post-pre)*(3*(per^2)-2*(per^3))+pre;
        end
    else
        beta_d1 = beta_d;
    end
end


%% Run simulation with FFD
for ind=2:N
    [zdot]               =   real(WindTurbine(z_FFD(:,ind-1),beta_d1(ind-1,1),V(ind-1,1),A_tilde,B_tilde,fun_data));
    z_FFD(:,ind)         =   z_FFD(:,ind-1)+Ts*zdot;
end

r_beta = z_FFD(8,:)';


% Constructing Beta_dot
beta_dot = zeros(N, 1);
beta_dot(1,1)=0; 
beta_dot(2:N,1)=(r_beta(2:N,1)-r_beta(1:N-1,1))/Ts;


%% Compute path constraints h(x)                            
h=[z_FFD(7,:)'+2.8*ones(N,1);  % omega_g constraints, limits are omega_gmax and omega_g min
   -z_FFD(7,:)'+2.8*ones(N,1);
   z_FFD(1,:)'+5*ones(N,1);  % y_t constraints 
   -z_FFD(1,:)'+5*ones(N,1);
   beta_dot+10*ones(N,1);    % beta_dot constraints
   -beta_dot+10*ones(N,1)];


%% F(x) computation
F(1:N,1)=Yref(1:N,1)-z_FFD(7,1:N)';        % term related to omega_g
F(N+1:2*N,1)=z_FFD(4,:)';                  % term related to yt_dot
F(1:2*N,1)=Q*F(1:2*N,1);
F(2*N+1:3*N,1)=R*beta_dot;
f=[F'*F;h];
end

