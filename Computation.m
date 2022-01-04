clear variables
close all
clc

syms phi(t) theta_1(t) theta_2(t) r(t) l m g k r_0 t

c_phi = cos(phi(t));
s_phi = sin(phi(t));
c_theta_1 = cos(theta_1(t));
s_theta_1 = sin(theta_1(t));
c_theta_2 = cos(theta_2(t));
s_theta_2 = sin(theta_2(t));

A = [c_phi -s_phi 0 0; s_phi c_phi 0 0 ; 0 0 1 0; 0 0 0 1];
B = [c_theta_1 0 -s_theta_1 0; 0 1 0 0 ; s_theta_1 0 c_theta_1 0 ;0 0 0 1];
C = [1 0 0 0; 0 1 0 0 ; 0 0 1 -r(t); 0 0 0 1];

%Amneh says this should rotate like this
D = [c_theta_2 0 -s_theta_2 0; 0 1 0 0 ; s_theta_2 0 c_theta_2 0 ;0 0 0 1];
E = [1 0 0 0; 0 1 0 0 ; 0 0 1 -l; 0 0 0 1];
T = A*B*C*D*E;
xyz = simplify(T(1:3,4), 'Steps', 10);
J = simplify(jacobian(xyz, [phi(t) r(t) theta_1(t) theta_2(t)]), 'Steps', 10);

%% Unit vector op
T_OP = A*B*C;
xyz_OP = simplify(T_OP(1:3,4), 'Steps', 10); 
u_OP = xyz_OP./(xyz_OP(1)^2+xyz_OP(2)^2+xyz_OP(3)^2)^(1/2);
u_OP = simplify(u_OP, 'Steps', 10); 

%% Conservative energies   
G_f = m*g*xyz(3);
K_f = 0.5*k*(r(t)-r_0)^2;

U = G_f + K_f; %not sure about the sign here
U_F_phi = diff(U,phi(t));
U_F_theta_1 = diff(U,theta_1(t));
U_F_theta_2 = diff(U,theta_2(t));
U_F_r = diff(U,r(t)); 

%% Kinetic energy
xyz_dot = simplify(J*[phi(t);theta_1(t);theta_2(t);r(t)], 'Steps',10);
T = 0.5*m*sum(xyz_dot.^2);

%% Lagrangian (using Jacobian)
L_a = T - U;

%% Lagrangian (not using Jacobian)
x_dot = diff(xyz(1),t);
y_dot = diff(xyz(2),t);
z_dot = diff(xyz(3),t);

T_b = 0.5*m*(x_dot^2 + y_dot^2 + z_dot^2);
L_b = T_b  - U;

%% diferential equations for Lagrangian
dL_dPhi = diff(L_b, phi(t));
dL_dPhi_dot = diff(L_b, diff(phi(t),t));
dt_dL_dPhi_dot = diff(dL_dPhi_dot,t);
eqPhi = dL_dPhi - dt_dL_dPhi_dot;
eqPhi = simplify(eqPhi, 'Steps', 10);

dL_dTheta_1 = diff(L_b, theta_1(t));
dL_dTheta_1_dot = diff(L_b, diff(theta_1(t),t));
dt_dL_dTheta_1_dot = diff(dL_dTheta_1_dot,t);
eqTheta_1 = dL_dTheta_1 - dt_dL_dTheta_1_dot;
eqTheta_1 = simplify(eqTheta_1, 'Steps', 10);

dL_dTheta_2 = diff(L_b, theta_2(t));
dL_dTheta_2_dot = diff(L_b, diff(theta_2(t),t));
dt_dL_dTheta_2_dot = diff(dL_dTheta_2_dot,t);
eqTheta_2 = dL_dTheta_2 - dt_dL_dTheta_2_dot;
eqTheta_2 = simplify(eqTheta_2, 'Steps', 10);

dL_dr = diff(L_b, r(t));
dL_dr_dot = diff(L_b, diff(r(t),t));
dt_dL_dr_dot = diff(dL_dr_dot,t);
eqR = dL_dr - dt_dL_dr_dot;
eqR = simplify(eqR, 'Steps', 10);

eqs = [eqPhi; eqTheta_1; eqTheta_2; eqR];

