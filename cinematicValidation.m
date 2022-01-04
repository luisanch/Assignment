clc
close all
clear variables
 
validateKinematics()

function y = validateKinematics()
phiRange = (0:1:360);
theta_1Range = (0:1:180);
theta_2Range = (0:1:180);

phiSpinValidation = phiSpin(phiRange, 90, 0, 1,1);
xyzPlot(phiSpinValidation, 'phiRange')
xyzPlot3(phiSpinValidation, 'phiRange')

th1SpinValidation = th1Spin(0, theta_1Range, 0, 1, 1);
xyzPlot(th1SpinValidation, 'theta_1Range')
xyzPlot3(th1SpinValidation, 'theta_1Range')

th2SpinValidation = th2Spin(0, 90, theta_2Range, 1, 1);
xyzPlot(th2SpinValidation, 'theta_2Range')
xyzPlot3(th2SpinValidation, 'theta_2Range')
end

function y = xyzPlot3(xyz, title)
figure('Name',title,'NumberTitle','off')
plot3(xyz(:,2),xyz(:,3),xyz(:,4))
ylabel('Y')
xlabel('X')
zlabel('Z')
end

function y = xyzPlot(xyz, title)
figure('Name',title,'NumberTitle','off')

plot(xyz(:,2))
hold on
plot(xyz(:,3))
hold on
plot(xyz(:,4))
hold off
legend({'x','y','z'},'Location','southwest')
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
B = [c_theta_1 0 -s_theta_1 0; 0 1 0 0 ; s_theta_1 0 c_theta_1 0 ;0 0 0 1];
C = [1 0 0 0; 0 1 0 0 ; 0 0 1 -r; 0 0 0 1];

%Amneh says this should rotate like this
D = [c_theta_2 0 -s_theta_2 0; 0 1 0 0 ; s_theta_2 0 c_theta_2 0 ;0 0 0 1];
E = [1 0 0 0; 0 1 0 0 ; 0 0 1 -l; 0 0 0 1];
T = A*B*C*D*E;
xyz = T(1:3,4).';
end  