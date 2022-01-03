clc
clear variables

validateKinematics()
function y = validateKinematics()
phiRange = (0:1:360);
theta_1Range = (0:1:180);
theta_2Range = (0:1:90); 

phiSpinValidation = phiSpin(phiRange, 90, 0, 1,1);
xyPlot(phiSpinValidation, 'phi')

th1SpinValidation = th1Spin(0, theta_1Range, 0, 1, 1);
xzPlot(th1SpinValidation, 'theta 1')
yzPlot(th1SpinValidation, 'theta 1')

th2SpinValidation = th2Spin(0, 90, theta_2Range, 1, 1);
xzPlot(th2SpinValidation, 'theta 2')
yzPlot(th2SpinValidation, 'theta 2')
end

function y = xyPlot(xyz, title)
figure('Name',title,'NumberTitle','off')
plot3(xyz(:,2),xyz(:,3),xyz(:,1)) 
ylabel('Y')
xlabel('X')
zlabel('deg') 
end

function y = xzPlot(xyz, title)
figure('Name',title,'NumberTitle','off')
plot3(xyz(:,2),xyz(:,4),xyz(:,1)) 
ylabel('Z')
xlabel('X')
zlabel('deg') 
end

function y = yzPlot(xyz, title)
figure('Name',title,'NumberTitle','off')
plot3(xyz(:,4),xyz(:,3),xyz(:,1)) 
ylabel('Y')
xlabel('Z')
zlabel('deg') 
end

function y = phiSpin(range, theta_1, theta_2, r,l)
y = zeros(size(range,2),3);
for index = range
y(index+1,:) = getCoordinatesM(index, theta_1, theta_2, r,l); 
end
y = [range',y];
end

function y = th1Spin(phi, range, theta_2, r,l)
y = zeros(size(range,2),3);
for index = range
y(index+1,:) = getCoordinatesM(phi, index, theta_2, r,l); 
end
y = [range',y];
end

function y = th2Spin(phi, theta_1, range, r,l)
y = zeros(size(range,2),3);
for index = range
y(index+1,:) = getCoordinatesM(phi, theta_1, index, r,l); 
end
y = [range',y];
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





