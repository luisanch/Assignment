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
xyz = T(1:3,4);

J = jacobian(xyz, [phi(t) r(t) theta_1(t) theta_2(t)]);

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

syms u(t)
%diff(u(t),t)

x_dot_u = cos(phi(t))*diff(u(t),t) - sin(phi(t))*u(t)*diff(phi(t),t);
y_dot_u = sin(phi(t))*diff(u(t),t) - cos(phi(t))*u(t)*diff(phi(t),t);
z_dot_u = l*sin(theta_1(t) + theta_2(t)) * diff(theta_2(t),t) - cos(theta_1(t))*diff(r(t),t) + diff(theta_1(t),t)*u(t);

%c(p)*(s(t1)*r' + l*c12*(t1' + t2') + c(t1)*r*t1') - s(p)*(l*s12 + s(t1)*r)*p'
%s(p)*(s(t1)*r' + l*c12*(t1' + t2') + c(t1)*r*t1') + c(p)*(l*s12 + s(t1)*r)*p'
%l*s12*(t1' + t2') - c(t1)*r' + s(t1)*r*t1'

xyz_dot = [x_dot;
    y_dot;
    z_dot]; 
v_sq = x_dot^2 + y_dot^2 + z_dot^2; 
T = 0.5*m*(x_dot^2 + y_dot^2 + z_dot^2);

T_u = 0.5*m*(x_dot_u^2 + y_dot_u^2 + z_dot_u^2);

%% Lagrangian
L = T - U;
L_u = T_u - U; 

%% diferential equations for Lagrangian
dL_dPhi = diff(L, phi(t));
dL_dPhi_dot = diff(L, diff(phi(t),t));
dt_dL_dPhi_dot = diff(dL_dPhi_dot,t);
Q_Phi = dL_dPhi - dt_dL_dPhi_dot; 

u_dL_dPhi = diff(L_u, phi(t));
u_dL_dPhi_dot = diff(L_u, diff(phi(t),t));
u_dt_dL_dPhi_dot = diff(dL_dPhi_dot,t);
u_Q_Phi = dL_dPhi - dt_dL_dPhi_dot; 

dL_dTheta_1 = diff(L, theta_1(t));
dL_dTheta_1_dot = diff(L, diff(theta_1(t),t));
dt_dL_dTheta_1_dot = diff(dL_dTheta_1_dot,t);
Q_Theta_1 = dL_dTheta_1 - dt_dL_dTheta_1_dot; 

dL_dTheta_2 = diff(L, theta_2(t));
dL_dTheta_2_dot = diff(L, diff(theta_2(t),t));
dt_dL_dTheta_2_dot = diff(dL_dTheta_2_dot,t);
Q_Theta_2 = dL_dTheta_2 - dt_dL_dTheta_2_dot; 

dL_dr = diff(L, r(t));
dL_dr_dot = diff(L, diff(r(t),t));
dt_dL_dr_dot = diff(dL_dr_dot,t);
Q_r = dL_dr - dt_dL_dr_dot; 

eqs = [ Q_Phi; Q_Theta_1; Q_Theta_2; Q_r];

%% rheonomic holonomic Constraint
lambda = dt_dL_dPhi_dot; %
lambda = simplify(lambda, 'Steps', 1000);
% m*(l*s12 + s1*r)*(2*s1*p!*r! + l*s12*p!! + s1*r*p!! + 2*c1*r*p!*t1! + 2*l*c12*p!*t1! + 2*l*c12*p!*t2!)
% 2*p!*m*(l*s12+r*s1)(s1*r!+c1*r*t1!+l*c12*t!+lc12t2!)
% u: (l*s12+r*s1), u!: (s1*r!+c1*r*t1!+l*c12*t!+lc12t2!)
% 2*p!*m*u*u!

test_r = k*r_0 - k*r(t) + g*m*cos(theta_1(t)) + (m*r(t)*diff(phi(t), t)^2)/2 - (l*m*cos(2*theta_1(t) + theta_2(t))*diff(phi(t), t)^2)/2 - (m*cos(2*theta_1(t))*r(t)*diff(phi(t), t)^2)/2 + (l*m*cos(theta_2(t))*diff(phi(t), t)^2)/2;
test_th2 =-l*m*(g*sin(theta_1(t) + theta_2(t)) - (l*sin(2*theta_1(t) + 2*theta_2(t))*diff(phi(t), t)^2)/2 - (sin(2*theta_1(t) + theta_2(t))*r(t)*diff(phi(t), t)^2)/2 + (sin(theta_2(t))*r(t)*diff(phi(t), t)^2)/2);
test_th1 = m*((sin(2*theta_1(t))*r(t)^2*diff(phi(t), t)^2)/2 - g*sin(theta_1(t))*r(t) + (l^2*sin(2*theta_1(t) + 2*theta_2(t))*diff(phi(t), t)^2)/2 - g*l*sin(theta_1(t) + theta_2(t)) + l*sin(2*theta_1(t) + theta_2(t))*r(t)*diff(phi(t), t)^2);

eqns = [test_r, test_th1, test_th2];
vars = [theta_2(t) theta_1(t)];

hmm = (m*r(t))/2 - (l*m*cos(2*theta_1(t) + theta_2(t)))/2 - (m*cos(2*theta_1(t))*r(t))/2 + (l*m*cos(theta_2(t)))/2;
centripetal = diff(phi(t), t)^2 * (l*sin(theta_1(t)+theta_2(t)) + r*sin(theta_1(t))) * m;
