clear variables
clc

syms phi theta_1 theta_2 r l

c_phi = cosd(phi);
s_phi = sind(phi);
c_theta_1 = cosd(theta_1);
s_theta_1 = sind(theta_1);
c_theta_2 = cosd(theta_2);
s_theta_2 = sind(theta_2);

A = [c_phi -s_phi 0 0; s_phi c_phi 0 0 ; 0 0 1 0; 0 0 0 1];
B = [c_theta_1 0 -s_theta_1 0; 0 1 0 0 ; s_theta_1 0 c_theta_1 0 ;0 0 0 1];
C = [1 0 0 0; 0 1 0 0 ; 0 0 1 -r; 0 0 0 1];
D = [c_theta_2 0 s_theta_2 0; 0 1 0 0 ; -s_theta_2 0 c_theta_2 0 ;0 0 0 1];
E = [1 0 0 0; 0 1 0 0 ; 0 0 1 -l; 0 0 0 1];
T = A*B*C*D*E;
xyz = T(1:3,4); 