clear variables
close all
clc

syms phi theta_1 theta_2 r l m g k r_0

c_phi = cos(phi);
s_phi = sin(phi);
c_theta_1 = cos(theta_1);
s_theta_1 = sin(theta_1);
c_theta_2 = cos(theta_2);
s_theta_2 = sin(theta_2);

A = [c_phi -s_phi 0 0; s_phi c_phi 0 0 ; 0 0 1 0; 0 0 0 1];
B = [c_theta_1 0 -s_theta_1 0; 0 1 0 0 ; s_theta_1 0 c_theta_1 0 ;0 0 0 1];
C = [1 0 0 0; 0 1 0 0 ; 0 0 1 -r; 0 0 0 1];

%Amneh says this should rotate like this
D = [c_theta_2 0 -s_theta_2 0; 0 1 0 0 ; s_theta_2 0 c_theta_2 0 ;0 0 0 1];
E = [1 0 0 0; 0 1 0 0 ; 0 0 1 -l; 0 0 0 1];
T = A*B*C*D*E;
xyz = simplify(T(1:3,4), 'Steps', 10);
J = simplify(jacobian(xyz, [phi r theta_1 theta_2]), 'Steps', 10);

%% Conservative energies  
T_OP = A*B*C;
xyz_OP = simplify(T_OP(1:3,4), 'Steps', 10); 
u_OP = xyz_OP./(xyz_OP(1)^2+xyz_OP(2)^2+xyz_OP(3)^2)^(1/2);
u_OP = simplify(u_OP, 'Steps', 10); 

G_f = m*g*xyz(3);
K_f = 0.5*k*(r-r_0)^2;

U = G_f + K_f;
U_F_phi = diff(U,phi);
U_F_theta_1 = diff(U,theta_1);
U_F_theta_2 = diff(U,theta_2);
U_F_r = diff(U,r); 

%% Kinetic energy
xyz_dot = simplify(J*[phi;theta_1;theta_2;r], 'Steps',10);
T = 0.5*m*sum(xyz_dot.^2);

%% Lagrangian
L = T - U;

