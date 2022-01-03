clear variables
clc

syms phi theta_1 theta_2 r l

c_phi = cos(phi);
s_phi = sin(phi);
c_theta_1 = cos(theta_1);
s_theta_1 = sin(theta_1);
c_theta_2 = cos(theta_2);
s_theta_2 = sin(theta_2);

A = [c_phi -s_phi 0 0; s_phi c_phi 0 0 ; 0 0 1 0; 0 0 0 1];
B = [c_theta_1 0 -s_theta_1 0; 0 1 0 0 ; s_theta_1 0 c_theta_1 0 ;0 0 0 1];
C = [1 0 0 0; 0 1 0 0 ; 0 0 1 -r; 0 0 0 1];
D = [c_theta_2 0 s_theta_2 0; 0 1 0 0 ; -s_theta_2 0 c_theta_2 0 ;0 0 0 1];
E = [1 0 0 0; 0 1 0 0 ; 0 0 1 -l; 0 0 0 1];
T = A*B*C*D*E;
xyz = simplify(T(1:3,4), 'Steps', 10);
J = simplify(jacobian(xyz, [phi r theta_1 theta_2]), 'Steps', 10);