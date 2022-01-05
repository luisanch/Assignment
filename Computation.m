clear variables
close all
clc

syms phi(t) theta_1(t) theta_2(t) r(t) l m g k r_0 t w lambda

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
x_dot = diff(xyz(1),t);
y_dot = diff(xyz(2),t);
z_dot = diff(xyz(3),t);

%c(p)*(s(t1)*r' + l*c12*(t1' + t2') + c(t1)*r*t1') - s(p)*(l*s12 + s(t1)*r)*p'
%s(p)*(s(t1)*r' + l*c12*(t1' + t2') + c(t1)*r*t1') + c(p)*(l*s12 + s(t1)*r)*p'
%l*s12*(t1' + t2') - c(t1)*r' + s(t1)*r*t1'

xyz_dot = [x_dot;
    y_dot;
    z_dot]; 
v_sq = simplify(x_dot^2 + y_dot^2 + z_dot^2, 'Steps', 10);

%(sp*(s1*r' + l*c12*(t1' + t2') + c1*r*t1') + cp*(l*s12 + s1*r)*p')^2 + (cp*(s1*r' + l*c12*(t1' + t2') + c1*r*t1') - sp*(l*s12 + s1*r)*p')^2 + (l*s12*(t1' + t2') - c1*r' + s1*r*t1')^2

T = 0.5*m*(x_dot^2 + y_dot^2 + z_dot^2);

%% Lagrangian
L = T - U;

%% diferential equations for Lagrangian
dL_dPhi = diff(L, phi(t));
dL_dPhi_dot = diff(L, diff(phi(t),t));
dt_dL_dPhi_dot = diff(dL_dPhi_dot,t);
Q_Phi = dL_dPhi - dt_dL_dPhi_dot;
Q_Phi = simplify( Q_Phi, 'Steps', 10);

dL_dTheta_1 = diff(L, theta_1(t));
dL_dTheta_1_dot = diff(L, diff(theta_1(t),t));
dt_dL_dTheta_1_dot = diff(dL_dTheta_1_dot,t);
Q_Theta_1 = dL_dTheta_1 - dt_dL_dTheta_1_dot;
Q_Theta_1 = simplify(Q_Theta_1, 'Steps', 10);

dL_dTheta_2 = diff(L, theta_2(t));
dL_dTheta_2_dot = diff(L, diff(theta_2(t),t));
dt_dL_dTheta_2_dot = diff(dL_dTheta_2_dot,t);
Q_Theta_2 = dL_dTheta_2 - dt_dL_dTheta_2_dot;
Q_Theta_2 = simplify(Q_Theta_2, 'Steps', 10);

dL_dr = diff(L, r(t));
dL_dr_dot = diff(L, diff(r(t),t));
dt_dL_dr_dot = diff(dL_dr_dot,t);
Q_r = dL_dr - dt_dL_dr_dot;
Q_r = simplify(Q_r, 'Steps', 10);

eqs = [ Q_Phi; Q_Theta_1; Q_Theta_2; Q_r];
%% rheonomic holonomic Constraint


