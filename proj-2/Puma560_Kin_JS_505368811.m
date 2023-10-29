%% Jacob Sayono

% 505368811
% MAE C163B
% Midterm

%% Part 1 - Kinematic Analysis

close; clear; clc

% DH PARAMETERS
syms th1 th2 th3 th4 th5 th6

% Modified DH parameters
alpha = [0, -pi/2, 0, -pi/2, pi/2, -pi/2];
a = [0, 0, 0.4318, 0.019, 0, 0];
d = [0, 0, 0.1254, 0.4318, 0, 0];
th = [th1, th2, th3, th4, th5, th6];

th1 = 0;
th2 = 0;
th3 = 0;
th4 = 0;
th5 = 0;
th6 = 0;


% FORWARD KINEMATICS
syms theta1 theta2 theta3 theta4 theta5 theta6 theta_t

T01 = [cos(theta1) -sin(theta1) 0 0;
sin(theta1) cos(theta1) 0 0;
0 0 1 0;
0 0 0 1];

T12 = [cos(theta2) -sin(theta2) 0 0;
0 0 1 0;
-sin(theta2) -cos(theta2) 0 0;
0 0 0 1];

T23 = [cos(theta3) -sin(theta3) 0 0.4318;
sin(theta3) cos(theta3) 0 0;
0 0 1 0.1254;
0 0 0 1];

T34 = [cos(theta4) -sin(theta4) 0 0.019;
0 0 1 0.4318;
-sin(theta4) -cos(theta4) 0 0;
0 0 0 1];

T45 = [cos(theta5) -sin(theta5) 0 0;
0 0 -1 0;
sin(theta5) cos(theta5) 0 0;
0 0 0 1];

T56 = [cos(theta6) -sin(theta6) 0 0;
0 0 1 0;
-sin(theta6) -cos(theta6) 0 0;
0 0 0 1];

T6t = [cos(theta_t) -sin(theta_t) 0 -0.1;
sin(theta_t) cos(theta_t) 0 0;
0 0 1 0.08;
0 0 0 1];

% find forward kinematics
T = T01 * T12 * T23 * T34 * T45 * T56;
T0t = T * T6t;
Tinv6t = inv(T6t);

% set some random joint angles in radians
q = [0.2, 0.3, -0.5, 0.4, 0.1, 0.8];

% check forward kinematics function
T = fk(q)

% separate parameters from matrix
position = T(1:3, 4)
orientation = tform2eul(T, 'XYZ')

% INVERSE KINEMATICS
[T1, T2, T3, T4, T5, T6, Tt] = IK(T)

% create table
PUMA_ik = SerialLink(L, 'name', 'PUMA_INVERSE_KINEMATICS')

%% Part 2 - Simulation

% FORWARD KINEMATICS

% DH parameters for Puma560
L1 = Link('revolute', 'alpha', 0, 'a', 0, 'd', 0, 'modified');
L2 = Link('revolute', 'alpha', -pi/2, 'a', 0, 'd', 0, 'modified');
L3 = Link('revolute', 'alpha', 0, 'a', 0.4318, 'd', 0.1254, 'modified');
L4 = Link('revolute', 'alpha', -pi/2, 'a', 0.019, 'd', 0.4318, 'modified');
L5 = Link('revolute', 'alpha', pi/2, 'a', 0, 'd', 0, 'modified');
L6 = Link('revolute', 'alpha', -pi/2, 'a', 0, 'd', 0, 'modified');
% L7 = Link('revolute', 'alpha', 0, 'a', -0.1, 'd', .008, 'modified');
tool = transl(0, -0.1, .008);

% create robot using the DH parameters
Puma560 = SerialLink([L1 L2 L3 L4 L5 L6], 'name', 'Puma 560', 'tool', tool);

% calculate the forward kinematics for the given joint angles
T = Puma560.fkine(q)


% INVERSE KINEMATICS

% set the desired end-effector position and orientation
position = [0.5, -0.3, 0.4];
orientation = [pi/2, 0, pi/2];

% convert the orientation to a rotation matrix
R = eul2r(orientation, 'XYZ');

% calculate the inverse kinematics for the desired pose
q = Puma560.ikine(transl(position) * rpy2tr(orientation), 'mask', [1 1 1 0 0 0])

% calculate the forward kinematics of the resulting joint angles to verify the solution
T = Puma560.fkine(q)

% visualize robot
figure(1)
Puma560.plot([T1, T2, T3, T4, T5, T6, Tt],'workspace',[-500,500,-500,500,0,500])

%% Function Definitions

% matrix function given DH parameters
function matrix = mat_from_DH(alpha_iminus1, a_iminus1, d_i, theta_i)
    matrix = [cosd(theta_i) -sind(theta_i) 0 a_iminus1;
              sind(theta_i)*cosd(alpha_iminus1) cosd(theta_i)*cosd(alpha_iminus1) -sind(alpha_iminus1) -sind(alpha_iminus1)*d_i;
              sind(theta_i)*sind(alpha_iminus1) cosd(theta_i)*sind(alpha_iminus1) cosd(alpha_iminus1) cosd(alpha_iminus1)*d_i;
              0 0 0 1];
end

% forward kinematics
function T = fk(theta1, theta2, theta3, theta4, theta5, theta6, thetat)
    % joint lengths
    a2 = 0.4318;
    a3 = 0.0191;
    d3 = 0.1254;
    d4 = 0.4318;
    
    T01 = mat_from_DH(0, 0, 0, theta1);
    T12 = mat_from_DH(-90, 0, 0, theta2);
    T23 = mat_from_DH(0, a2, d3, theta3);
    T34 = mat_from_DH(-90, a3, d4, theta4);
    T45 = mat_from_DH(90, 0, 0, theta5);
    T56 = mat_from_DH(-90, 0, 0, theta6);
    T6t = mat_from_DH(0, -0.1, 0.08, thetat);
    
    T = T01*T12*T23*T34*T45*T56*T6t;

end

% inverse kinematics
function [t1, t2, t3, t4, t5, t6] = ik(T_location)
    % joint lengths
    a2 = 0.4318;
    a3 = 0.0191;
    d3 = 0.1254;
    d4 = 0.4318;

    % transformation matrices
    T_Gto6 = [1 0 0 0; 0 1 0 0; 0 0 1 0.05625; 0 0 0 1];
    T_TtoG = [1 0 0 -0.1; 0 1 0 0; 0 0 1 0.08; 0 0 0 1];
    T_input = T_location * inverse(T_TtoG) * inverse(T_Gto6);

    % extract variables
    r = T_input(1:3,1:3);
    p = T_input(1:3,4);
    px = p(1);
    py = p(2);
    pz = p(3);

    % compute inverse kinematics
    t1 = atan2(py, px);
    t3 = atan2(a3, -d4) - atan2(sqrt(1 - ((px^2 + py^2 + pz^2 - a2^2 - a3^2 - d3^2 - d4^2)/(2*a2))^2), (px^2 + py^2 + pz^2 - a2^2 - a3^2 - d3^2 - d4^2)/(2*a2));
    t23 = atan2((a3 + a2*cos(t3))*pz - (cos(t1)*px + sin(t1)*py)*(a2*sin(t3) - d4), (a3 + a2*cos(t3))*(cos(t1)*px + sin(t1)*py) + (a2*sin(t3) - d4)*pz);
    t2 = t23 - t3;
    t4 = atan2(r(3,1)*sin(t1) - r(3,2)*cos(t1), -r(1,3)*sin(t1)*cos(t2+t3) + r(2,3)*cos(t1)*cos(t2+t3) - r(3,3)*sin(t2+t3));
    t5 = atan2(sqrt(1 - (r(1,3)*sin(t1)*cos(t2+t3) - r(2,3)*cos(t1)*cos(t2+t3) + r(3,3)*sin(t2+t3))^2), r(1,3)*cos(t1)*cos(t2+t3) + r(2,3)*sin(t1)*cos(t2+t3) + r(3,3)*cos(t2+t3));
    t6 = atan2(-r(1,2)*sin(t1)*cos(t4) - r(2,2)*cos(t1)*cos(t4) + r(3,2)*sin(t4), r(1,1)*sin(t1)*cos(t4) + r(2,1)*cos(t1)*cos(t4) - r(3,1)*sin(t4)*cos(t2+t3));
end
