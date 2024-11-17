function [z_dot]=WindTurbine(z,beta_d,V,A_tilde,B_tilde,fun_data)
%% SYSTEM DYNAMIC EQUATIONS

Ve=V-(z(4)+z(5)*fun_data(1));
lambda=z(6)*fun_data(2)/Ve;
la_if = 1/(lambda+0.08*z(8))-0.035/(1+z(8)^3);
Cp = 0.5*(116*la_if-0.4*z(8)-5)*exp(-(21*la_if));
Tr=Cp/lambda*0.5*pi*fun_data(2)^3*fun_data(3)*Ve^2;
phi=acos(27*Cp/8-1);
a=2/3*cos(phi/3+2*pi/3)+2/3; 
Ct=real(4*a*(1-a)*(2*a*(1-a)/lambda^2+1)); 
Ft=Ct*0.5*pi*fun_data(2)^2*fun_data(3)*Ve^2;
Tg=fun_data(4)*(z(7)-fun_data(5)/97);
z_dot= A_tilde*z + B_tilde*[Ft;Tr;Tg;beta_d];
end





