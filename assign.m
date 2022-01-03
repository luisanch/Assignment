clc
clear variables

phiRange = (0:1:360);
theta_1Range = (0:1:180);
theta_2Range = (0:1:180); 

phiSpinValidation = phiSpin(phiRange, 90, 0, 1,1);

function y = phiSpin(range, theta_1, theta_2, r,l)
y = zeros(size(range,2),3);
for index = range
y(index+1,:) = getCoordinatesM(index, theta_1, theta_2, r,l); 
end
end

function xyz = getCoordinatesM(phi, theta_1, theta_2, r,l)  
c_phi = cosd(phi);
s_phi = sind(phi);
c_theta_1 = cosd(theta_1);
s_theta_1 = sind(theta_1);
c_theta_2 = cosd(theta_2);
s_theta_2 = sind(theta_2);

A = [c_phi -s_phi 0 0; s_phi c_phi 0 0 ; 0 0 1 0; 0 0 0 1];
B = [1 0 0 0; 0 c_theta_1 s_theta_1 0 ; 0 -s_theta_1 c_theta_1 0 ;0 0 0 1];
C = [1 0 0 0; 0 1 0 0 ; 0 0 1 -r; 0 0 0 1];
D = [1 0 0 0; 0 c_theta_2 -s_theta_2 0; 0 s_theta_2 c_theta_2 0; 0 0 0 1];
E = [1 0 0 0; 0 1 0 0 ; 0 0 1 -l; 0 0 0 1];  
T = A*B*C*D*E;
xyz = T(1:3,4)';
end 





