%% MAE C163B Project 3
% Jacob Sayono
% 505368811

%% Jacobian Derivation

% VELOCITY PROPAGATION
close all; clear all; clc;

syms theta1 theta2 theta3 l1 l2 d1 d4

% Modified Denavit-Hartenberg (DH) parameters
alpha = [0, 0, 0, 0];
a = [0, l1, l2, 0];
d = [d1, 0, 0, -d4];
theta = [theta1, theta2, theta3, 0];
L(1) = Link('revolute', 'alpha', alpha(1), 'a', a(1), 'd', d(1), 'modified');
L(2) = Link('revolute', 'alpha', alpha(2), 'a', a(2), 'd', d(2), 'modified');
L(3) = Link('revolute', 'alpha', alpha(3), 'a', a(3), 'd', d(3), 'modified');
L(4) = Link('prismatic', 'alpha', alpha(4), 'a', a(4), 'theta', 0, 'modified');
scara_robot = SerialLink(L, 'name', 'scara_robot')

% Joint angles
syms theta1 theta2 theta3 d4
q = [theta1, theta2, theta3, d4];

% Transformation matrices
T01 = scara_robot.A(1, q);
T12 = scara_robot.A(2, q);
T23 = scara_robot.A(3, q);
T34 = scara_robot.A(4, q);
T04 = scara_robot.A([1 2 3 4], q);
T04 = simplify(T04)

% Calculate rotation matrices and positions for each link
[R01, P01] = tr2rt(T01); 
[R12, P12] = tr2rt(T12); 
[R23, P23] = tr2rt(T23); 
[R34, P34] = tr2rt(T34); 
[R10, ~] = tr2rt(inv(T01)); 
[R21, ~] = tr2rt(inv(T12)); 
[R32, ~] = tr2rt(inv(T23)); 
[R43, ~] = tr2rt(inv(T34)); 
[R04, P04] = tr2rt(T04);
[R02, P02] = tr2rt(simplify(T01*T12));
[R03, P03] = tr2rt(simplify(T01*T12*T23));

% Define joint velocities
syms dq1 dq2 dq3 dd4
dq = [dq1; dq2; dq3; dd4];

% Calculate angular velocities for each frame
w00 = [0; 0; 0];
w11 = simplify(R10*w00 + [0; 0; dq1]);
w22 = simplify(R21*w11 + [0; 0; dq2]);
w33 = simplify(R32*w22 + [0; 0; dq3]);
w44 = simplify(R43*w33 + [0; 0; 0]);

% Define velocities for frame {0}
v0 = [0; 0; 0];

% Calculate velocities for each link
v01 = cross(w00, P01) + v0; 
v11 = simplify(R10*v01);
v12 = cross(w11, P12) + v11; 
v22 = simplify(R21*v12);
v23 = cross(w22, P23) + v22; 
v33 = simplify(R32*v23);
v34 = cross(w33, P34) + v33; 
v44 = simplify(R43*v34 + [0; 0; dd4]);

% Calculate the Jacobian matrix
JV4 = simplify(jacobian(v44, dq)); % Linear velocity
JW4 = simplify(jacobian(w44, dq)); % Angular velocity
J4_VP = simplify([JV4; JW4]) % Combined velocity

% Jacobian expressed in frame {0}
JT = [R04, zeros(3,3);zeros(3,3), R04];
J0_VP = simplify(JT * J4_VP)

% FORCE PROPAGATION
syms f_x f_y f_z n_x n_y n_z
f_end = [f_x f_y f_z n_x n_y n_z];
f4 = [f_x;f_y;f_z];
f3 = R34 * f4;
f2 = R23 * f3;
f1 = R12 * f2; f1 = simplify(collect(f1, f_end));
n4 = [n_x; n_y; n_z];
n3 = R34*n4 + cross(P34,f3); n3 = simplify(collect(n3, f_end));
n2 = R23*n3 + cross(P23,f2); n2 = simplify(collect(n2, f_end));
n1 = R12*n2 + cross(P12,f1); n1 = simplify(collect(n1, f_end));

% Torque Matrix
tq1 = collect(transpose(n1) * [0;0;1], [f_x f_y f_z n_x n_y n_z]);
tq2 = collect(transpose(n2) * [0;0;1], [f_x f_y f_z n_x n_y n_z]);
tq3 = collect(transpose(n3) * [0;0;1], [f_x f_y f_z n_x n_y n_z]);
tq4 = collect(transpose(f4) * [0;0;1], [f_x f_y f_z n_x n_y n_z]);
TQ1 = equationsToMatrix(tq1, [f_x f_y f_z n_x n_y n_z]);
TQ2 = equationsToMatrix(tq2, [f_x f_y f_z n_x n_y n_z]);
TQ3 = equationsToMatrix(tq3, [f_x f_y f_z n_x n_y n_z]);
TQ4 = equationsToMatrix(tq4, [f_x f_y f_z n_x n_y n_z]);
JT = [R04, zeros(3,3); zeros(3,3), R04];
J4_FP_T = [TQ1; TQ2; TQ3; TQ4];
J4_FP = transpose(J4_FP_T);
J4_FP = simplify(JT * J4_FP)

% EXPLICIT METHOD
Jxyz1 = diff(P04, theta1);
Jxyz2 = diff(P04, theta2);
Jxyz3 = diff(P04, theta3);
Jxyz4 = diff(P04, d4);
Jxyz = [Jxyz1 Jxyz2 Jxyz3 Jxyz4];
Orn = [0; 0; theta1 + theta2 + theta3]; % roll pitch yaw
JOrn1 = diff(Orn, theta1); % derivative of rpy wrp theta 1
JOrn2 = diff(Orn, theta2);
JOrn3 = diff(Orn, theta3);
JOrn4 = diff(Orn, d4);
JOrn = [JOrn1 JOrn2 JOrn3 JOrn4];
Jdiff = [Jxyz; JOrn]

% Using Siciliano
Jxyzs1 = cross(R01(:,3), P04 - P01);
Jxyzs2 = cross(R02(:,3), P04 - P02);
Jxyzs3 = cross(R03(:,3), P04 - P03);
Jxyzs4 = R04(:,3);
Jxyzs = [Jxyzs1 Jxyzs2 Jxyzs3 Jxyzs4];
JOrns1 = R01(:,3);
JOrns2 = R02(:,3);
JOrns3 = R03(:,3);
JOrns4 = zeros(3,1);
JOrns = [JOrns1 JOrns2 JOrns3 JOrns4];
Jdiffs = [Jxyzs; JOrns]

%% Singularities

% POSE DEFINITION
J_reduced = [Jdiff(1:3,: ); Jdiff(6,:)];
find_singularity = simplify(det(J_reduced));
sol1 = solve(find_singularity == 0, theta2)
sol2 = solve([find_singularity == 0, theta2 ~= sol1], theta2)
sol3 = solve([find_singularity == 0, theta2 ~= sol1, theta2 ~= sol2], theta2)

% IMPLICATIONS
% The following are the implications of the singular configuration:
% 1) The manipulator loses one degree of freedom.
%    As a result, the end effector can only move along the tangent direction of the arm.
%    Movement along the radial direction is not possible.
% 2) In this singular configuration,
%    it is possible to apply a finite force to the end effector that does not produce any torque at the robot's joints.
%    Consequently, the manipulator can "lock up."
% 3) The determinant of the Jacobian matrix is zero.
% 4) theta2 = -pi so the manipulator folds back from the other direction

%% Jacobian Ellipsoid

% From TA posted code:
close all; clear all; clc;

syms l1 l2 theta1 theta2
alpha = [0, 0];
a = [l1, l2];
d = [0, 0];
L(1) = Link('revolute', 'alpha', alpha(1), 'a', a(1), 'd', d(1), 'standard');
L(2) = Link('revolute', 'alpha', alpha(2), 'a', a(2), 'd', d(2), 'standard');
RR_robot = SerialLink(L, 'name', 'RR_robot');

% Define Jacobian Matrix
q_sym = [theta1 theta2];
T01 = RR_robot.A(1, q_sym);
T12 = RR_robot.A(2, q_sym);
T02 = simplify(T01*T12);
[R01, P01] = tr2rt(T01); R10 = transpose(R01);
[R12, P12] = tr2rt(T12); R21 = transpose(R12);
[R02, P02] = tr2rt(simplify(T01*T12));
Jxyz1 = diff(P02, theta1);
Jxyz2 = diff(P02, theta2);
Jxyz = [Jxyz1 Jxyz2];
roll = 0;
pitch = 0;
yaw = theta1 + theta2;
Orn = [roll; pitch; yaw];
JOrn1 = diff(Orn, theta1);
JOrn2 = diff(Orn, theta2);
JOrn = [JOrn1 JOrn2];
Jdiff = [Jxyz; JOrn]

% Define eigenvalues and eigenvectors for each config

% Config 1: L1 = 0.25, L2 = 0.75
config1 = subs(Jdiff, [l1,l2,theta1,theta2], [0.25, 0.75, 0, 0]);
% Get reduced Jacobian
Jconfig1 = config1(1:2, 1:2)
Jconfig1_T = transpose(Jconfig1);
% Get V1 and D1
[V1,D1] = eig(Jconfig1*Jconfig1_T)

% Config 2: L1 = 0.5, L2 = 0.5
config2 = subs(Jdiff, [l1,l2,theta1,theta2], [0.5, 0.5, 0, 0]);
% Get reduced Jacobian
Jconfig2 = config2(1:2, 1:2)
Jconfig2_T = transpose(Jconfig2);
% Get V2 and D2
[V2,D2] = eig(Jconfig2*Jconfig2_T)

% Config 3: L1 = 0.75, L2 = 0.25
config3 = subs(Jdiff, [l1,l2,theta1,theta2], [0.75, 0.25, 0, 0]);
% Get reduced Jacobian
Jconfig3 = config3(1:2, 1:2)
Jconfig3_T = transpose(Jconfig3);
% Get V3 and D3
[V3,D3] = eig(Jconfig3*Jconfig3_T)

% Plot

% Config 1
alpha = [0, 0];
a = [0.25, 0.75];
d = [0, 0];

L(1) = Link('revolute', 'alpha', alpha(1), 'a', a(1), 'd', d(1), 'standard');
L(2) = Link('revolute', 'alpha', alpha(2), 'a', a(2), 'd', d(2), 'standard');

RR_robot = SerialLink(L, 'name', 'RR_robot');

% can reach from 0.5 to 1
q1 = zeros(6,2);
% to reach x = 0.5, theta1 = pi, theta2 = 0
q1(1,:) = [pi, 0];
% to reach x = 1, theta1 = 0, theta2 = 0
q1(6,:) = [0, 0];

for i = 2:5
q1(i,:) = ik(0.25, 0.75, 0.5+(i-1)*0.1);
figure(i-1);
xlim([-inf inf]);
ylim([-inf inf]);
zlim([-inf inf]);
RR_robot.plot(q1(i,:))
RR_robot.fellipse(q1(i,:))
RR_robot.vellipse(q1(i,:))
% saveas(figure(i-1),sprintf('config1_%d.jpg',i-1))
end

% Config 2
alpha = [0, 0];
a = [0.5, 0.5];
d = [0, 0];

L(1) = Link('revolute', 'alpha', alpha(1), 'a', a(1), 'd', d(1), 'standard');
L(2) = Link('revolute', 'alpha', alpha(2), 'a', a(2), 'd', d(2), 'standard');

RR_robot2 = SerialLink(L, 'name', 'RR_robot');

% can reach from 0 to 1
q2 = zeros(11,2);
% to reach x = 0, theta1 = 0, theta2 = pi
q2(1,:) = [0, pi];
% to reach x = 1, theta1 = 0, theta2 = 0
q2(11,:) = [0, 0];

for i = 2:10
q2(i,:) = ik(0.5, 0.5, (i-1)*0.1);
figure(i-1);
xlim([-inf inf]);
ylim([-inf inf]);
zlim([-inf inf]);
RR_robot2.plot(q2(i,:))
RR_robot2.fellipse(q2(i,:))
RR_robot2.vellipse(q2(i,:))
% saveas(figure(i-1),sprintf('config2_%d.jpg',i-1))
end

% Config 3
alpha = [0, 0];
a = [0.75, 0.25];
d = [0, 0];
L(1) = Link('revolute', 'alpha', alpha(1), 'a', a(1), 'd', d(1), 'standard');
L(2) = Link('revolute', 'alpha', alpha(2), 'a', a(2), 'd', d(2), 'standard');

RR_robot3 = SerialLink(L, 'name', 'RR_robot');

% can reach from 0.5 to 1
q3 = zeros(6,2);
% to reach x = 0, theta1 = 0, theta2 = pi
q3(1,:) = [0, pi];
% to reach x = 1, theta1 = 0, theta2 = 0
q3(6,:) = [0, 0];

for i = 2:5
 q3(i,:) = ik(0.75, 0.25, 0.5+(i-1)*0.1);
 figure(i-1);
 xlim([-inf inf]);
 ylim([-inf inf]);
 zlim([-inf inf]);
 RR_robot3.plot(q3(i,:))
 RR_robot3.fellipse(q3(i,:))
 RR_robot3.vellipse(q3(i,:))
%  saveas(figure(i-1),sprintf('config3_%d.jpg',i-1))
end

% Best position for configurations

% Config 1
w1 = zeros(1, 4);
for i = 6:9
    w1(1, i-5) = w_func(0.25, 0.75, i*0.1);
end

max = w1(1,1);
for i = 1:4
if w1(1,i) > max
max = w1(1,i);
end
end
pos = find(w1 == max)

% ANSWER: best position for configuration 1 is x = 0.8

% Config 2
w2 = zeros(1, 9);
for i = 1:9
w2(1, i) = w_func(0.5, 0.5, i*0.1);
end

max = w2(1,1);
for i = 1:9
if w2(1,i) > max
max = w2(1,i);
end
end
pos = find(w2 == max)

% ANSWER: best position for configuration 2 is x = 0.7

% Config 3
w3 = zeros(1, 4);
for i = 6:9
 w3(1, i-5) = w_func(0.75, 0.25, i*0.1);
end
max = w3(1,1);
for i = 1:4
 if w3(1,i) > max
 max = w3(1,i);
 end
end
pos = find(w3 == max)

% ANSWER: best position for configuration 3 is x = 0.8

%% Arm Optimization
close all; clear all; clc;

% Brute Force
syms th1 th2 th3 L1 L2
% Modified DH parameters
alpha = [0, 0, 0, 0];
a = [0, L1, L2, 0];
d = [0.4, 0, 0, -0.15];
th = [th1, th2, th3, 0];
L(1) = Link('revolute', 'alpha', alpha(1), 'a', a(1), 'd', d(1), 'modified');
L(2) = Link('revolute', 'alpha', alpha(2), 'a', a(2), 'd', d(2), 'modified');
L(3) = Link('revolute', 'alpha', alpha(3), 'a', a(3), 'd', d(3), 'modified');
L(4) = Link('prismatic', 'alpha', alpha(4), 'a', a(4), 'theta', 0, 'modified');
SCARA = SerialLink(L, 'name', 'SCARA')

% Joint angle
syms q1 q2 q3 d4
q = [q1 q2 q3 0.15];
% Transformations, rotations, translations
T01 = SCARA.A(1,q);
T12 = SCARA.A(2,q);
T23 = SCARA.A(3,q);
T34 = SCARA.A(4,q);
T04 = SCARA.A([1 2 3 4],q);
T04 = simplify(T04)

[R01, P01] = tr2rt(T01); R10 = transpose(R01);
[R12, P12] = tr2rt(T12); R21 = transpose(R12);
[R23, P23] = tr2rt(T23); R32 = transpose(R23);
[R34, P34] = tr2rt(T34); R43 = transpose(R34);
[R04, P04] = tr2rt(T04);
[R02, P02] = tr2rt(simplify(T01*T12));
[R03, P03] = tr2rt(simplify(T01*T12*T23));
Jxyzs1 = cross(R01(:,3), P04 - P01);
Jxyzs2 = cross(R02(:,3), P04 - P02);
Jxyzs3 = cross(R03(:,3), P04 - P03);
Jxyzs4 = R04(:,3);
Jxyzs = [Jxyzs1 Jxyzs2 Jxyzs3 Jxyzs4]

JOrns1 = R01(:,3);
JOrns2 = R02(:,3);
JOrns3 = R03(:,3);
JOrns4 = zeros(3,1);
JOrns = [JOrns1 JOrns2 JOrns3 JOrns4];
Jdiffs = [Jxyzs; JOrns];
J0 = [Jdiffs(1:2,1:3); Jdiffs(6,1:3)]

J0T = transpose(J0);
JJT = simplify(J0*J0T)

l1 = [0.01:0.01:0.5];
l2 = [0.01:0.01:0.5];
h = 0.1;
w = 0.14;
d = 0.2;
count = 0;
kappa_array = [];
C= [];
l1_plot = [];
l2_plot = [];
X= linspace(d, d+w, 15);
Y= linspace(-h/2, h/2, 11);
for i=1:length(l1)
    for j=1:length(l2)
        kappa_array = [];
        if (abs(l1(i)-l2(j)) < d) && (l1(i)+l2(j) > ((d+w)^2 + (h/2)^2)^0.5)
            count = count+1;
            for m=1:length(X)
                for n=1:length(Y)
                    [theta] = inverseKinematicsScara3(X(m), Y(n), l1(i), l2(j));
                    JJTsubs = subs(JJT, [q1, q2, L1, L2], [theta(1), theta(2), l1(i), l2(j)]);
                    JJTsubs = vpa(JJTsubs);
                    lambda = eig(JJTsubs);
                    lambda_max = max(lambda);
                    lambda_min = min(lambda);
                    kappa = vpa((lambda_min/lambda_max)^0.5);
                    kappa_array = [kappa_array; kappa];
                end
            end
            C(count) = (sum(kappa_array)*min(kappa_array))/(l1(i)^3 + l2(j)^3);
            l1_plot(count) = l1(i);
            l2_plot(count) = l2(j);
            fprintf('Count = %d \t C-Value = %f \n',count, C(count));
        end
    end
end


% Plot values of Goal function
scatter3(l1_plot,l2_plot,C, 'filled', 'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75])
xlabel('L1', 'FontWeight','bold')
ylabel('L2', 'FontWeight','bold')
zlabel('C', 'FontWeight','bold')

% Find optimal values
[Cmax,Imax] = max(C)
L1_max = l1_plot(Imax)
L2_max = l2_plot(Imax)

% Plot value of Ki for entire workspace for best and worst combinations
kappa_max = [];
for m=1:length(X)
    for n=1:length(Y)
        [theta] = inverseKinematicsScara3(X(m), Y(n), L1_max, L2_max);
        JJTsubs = subs(JJT, [q1, q2, L1, L2], [theta(1), theta(2), L1_max, L2_max]);
        JJTsubs = vpa(JJTsubs);
        lambda = eig(JJTsubs);
        lambda_max = max(lambda);
        lambda_min = min(lambda);
        kappa = vpa((lambda_min/lambda_max)^0.5);
        kappa_max(m,n) = kappa;
    end
end

surf(X,Y, kappa_max')
xlabel('X(m)', 'FontWeight','bold')
ylabel('Y(m)', 'FontWeight','bold')
zlabel('K', 'FontWeight','bold', 'FontSize', 14)

% Here, we can see that the plot is symmetric
% Based on the C_max, the best manipulability point is at x = 0.24
% while the worst manipulability point is at x = 0.20


%% Functions

% IK function
function q = ik(l1, l2, x)
 c1 = (l1^2 + x^2 - l2^2)/(2*l1*x);
 a = (l1^2 - l2^2 + x^2) / (2*x);
 h = sqrt(l1^2 - a^2);
 s1 = h/l1;
 theta1 = atan2(s1, c1);
 ca = (x-a)/l2;
 sa = h/l2;
 a = atan2(sa, ca);
 theta2 = -(theta1+a);
 q = [theta1 theta2];
end

% w function
function w = w_func(L1, L2, x, Jdiff)
    syms l1 l2 theta1 theta2
    
    alpha = [0, 0];
    a = [l1, l2];
    d = [0, 0];
    
    L(1) = Link('revolute', 'alpha', alpha(1), 'a', a(1), 'd', d(1), 'standard');
    L(2) = Link('revolute', 'alpha', alpha(2), 'a', a(2), 'd', d(2), 'standard');
    
    RR_robot = SerialLink(L, 'name', 'RR_robot');
    q_sym = [theta1 theta2];
    
    T01 = RR_robot.A(1, q_sym);
    T12 = RR_robot.A(2, q_sym);
    T02 = simplify(T01*T12);

    [R01, P01] = tr2rt(T01); R10 = transpose(R01);
    [R12, P12] = tr2rt(T12); R21 = transpose(R12);
    [R02, P02] = tr2rt(simplify(T01*T12));

    Jxyz1 = diff(P02, theta1);
    Jxyz2 = diff(P02, theta2);
    Jxyz = [Jxyz1 Jxyz2];

    roll = 0;
    pitch = 0;
    yaw = theta1 + theta2;
    Orn = [roll; pitch; yaw];
    JOrn1 = diff(Orn, theta1);
    JOrn2 = diff(Orn, theta2);
    JOrn = [JOrn1 JOrn2];
    Jdiff = [Jxyz; JOrn];
    
    theta = ik(L1, L2, x);
    q1 = theta(1);
    q2 = theta(2);
    
    config = subs(Jdiff, [l1,l2,theta1,theta2], [L1, L2, q1, q2]);
    Jconfig = config(1:2, 1:2); % reduced Jacobian
    Jconfig_T = transpose(Jconfig);
    
    w = sqrt(det(Jconfig*Jconfig_T));
end

function [theta] = inverseKinematicsScara3(X, Y, a1, a2)
    c2 = (X^2 + Y^2 - a1^2 - a2^2)/(2*a1*a2);
    s2 = (1 - c2^2)^(0.5);

    t2s1 = atan2(s2, c2);
    t2s2 = atan2(-s2, c2);

    if(t2s1 >= 0)
        theta2 = t2s1;
    else
        theta2 = t2s2;
    end

    beta = atan2(X, Y);
    psi = acos((X^2 + Y^2 + a1^2 - a2^2)/(2*a1*(X^2 + Y^2)^(0.5)));
    theta1 = -(beta + psi - pi/2);
    theta = [theta1; theta2];
end
