%% Jacob Sayono

% MAE C163B

% Midterm Report

% Part 3


%% Forward Kinematics

% Symbolic Expressions

syms t1 t2 t3 t4 t5 t6 a2 a3 d2 d3 d4

DH = [      0        0       0       t1;     %alpha, a, d, theta
            -pi/2     0       d2      t2;
            0       a2      d3      t3;
            pi/2    a3      d4      t4;
            -pi/2   0       0       t5;
            pi/2    0       0       t6
        ]

T_01 = transformationMatrix(DH(1,:));
T_12 = transformationMatrix(DH(2,:));
T_23 = transformationMatrix(DH(3,:));
T_34 = transformationMatrix(DH(4,:));
T_45 = transformationMatrix(DH(5,:));
T_56 = transformationMatrix(DH(6,:));

T_06 = T_01*T_12*T_23*T_34*T_45*T_56; 
T_06 = simplify(T_06)



% Numercial Expressions (Modified DH Parameters)


DH = [      0        0       0       t1;     %alpha, a, d, theta
            -pi/2     0       .2435      t2;
            0       0.4318      -0.0934      t3;
            pi/2    -0.0203      0.4331      t4;
            -pi/2   0       0       t5;
            pi/2    0       0       t6
        ]

T_01 = transformationMatrix(DH(1,:));
T_12 = transformationMatrix(DH(2,:));
T_23 = transformationMatrix(DH(3,:));
T_34 = transformationMatrix(DH(4,:));
T_45 = transformationMatrix(DH(5,:));
T_56 = transformationMatrix(DH(6,:));

% tool
T_6T = [    1       0       0       -0.1;
            0       1       0       0;
            0       0       1       0.13625;
            0       0       0       1];

% combine all expressions
T_06 = T_01*T_12*T_23*T_34*T_45*T_56;
T_06 = vpa(simplify(T_06))
T_0T = T_06*T_6T;
T_0T = vpa(simplify(T_0T))
T_T6 = (T_6T)^(-1)



%% Inverse Kinematics

% Numerical Expressions

% DH = [      0        0       0       t1;     %alpha, a, d, theta
%             -pi/2     0       0.2435      t2;
%             0       0.4318      -0.0934      t3;
%             pi/2    -0.0203      0.4331      t4;
%             -pi/2   0       0       t5;
%             pi/2    0       0       t6
%         ];

syms t1 t2 t3 t4 t5 t6 a2 a3 d2 d3 d4

DH = [      0        0       0       t1;     %alpha, a, d, theta
            -pi/2     0       d2      t2;
            0       a2      d3      t3;
            pi/2    a3      d4      t4;
            -pi/2   0       0       t5;
            pi/2    0       0       t6
        ]

T_01 = transformationMatrix(DH(1,:));
T_12 = transformationMatrix(DH(2,:));
T_23 = transformationMatrix(DH(3,:));
T_34 = transformationMatrix(DH(4,:));
T_45 = transformationMatrix(DH(5,:));
T_56 = transformationMatrix(DH(6,:));

T_06 = T_01*T_12*T_23*T_34*T_45*T_56;
vpa(simplify(combine(T_06)), 5)


% T_04 = T_01*T_12*T_23*T_34;
% T_04 = vpa(simplify(T_04), 5)

syms R11 R12 R13 R21 R22 R23 R31 R32 R33 Px Py Pz

Goal = [R11 R12 R13 Px;
        R21 R22 R23 Py;
        R31 R32 R33 Pz;
        0   0   0   1];



T_01_inverse = simplify(T_01^-1);

Step_1_LHS = T_01_inverse * Goal 
Step_1_RHS = vpa(simplify(T_12*T_23*T_34),5)


T_02_inverse = simplify(T_12^-1)*T_01_inverse;
Step_2_LHS = (T_02_inverse * Goal) 
Step_2_RHS = (simplify(T_23*T_34))

T_03 = T_01*T_12*T_23
T_03_inverse = simplify(T_03^-1)
Step_3_LHS = (T_03_inverse * Goal)
Step_3_RHS = (simplify(T_34*T_45*T_56))

T_04 = T_01*T_12*T_23*T_34
T_04_inverse = simplify(T_04^-1)
Step_4_LHS = (T_04_inverse * Goal)
Step_4_RHS = (simplify(T_45*T_56))

T_05 = T_01*T_12*T_23*T_34*T_45
T_05_inverse = simplify(T_05^-1)
Step_5_LHS = (T_05_inverse * Goal)
Step_5_RHS = (simplify(T_56))

% test IK
DH = [      0        0       0       pi/2;     %alpha, a, d, theta
            -pi/2     0       0.2435        -pi/2;
            0       0.4318      -0.0934      0;
            pi/2    -0.0203      0.4331      0;
            -pi/2   0       0      pi/3;
            pi/2    0       0       0
        ]

T_01 = transformationMatrix(DH(1,:));
T_12 = transformationMatrix(DH(2,:));
T_23 = transformationMatrix(DH(3,:));
T_34 = transformationMatrix(DH(4,:));
T_45 = transformationMatrix(DH(5,:));
T_56 = transformationMatrix(DH(6,:));


T_06 = T_01*T_12*T_23*T_34*T_45*T_56

theta = IKPuma(T_06, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331)


%% Trajectory Generation (Joint Space)
clear all
close all

% Numercial Expressions

DH = [      0        0       0       0;     %alpha, a, d, theta
            -pi/2     0       0.2435      0;
            0       0.4318      -0.0934      pi;
            pi/2    -0.0203      0.4331      0;
            -pi/2   0       0       0;
            pi/2    0       0       0
        ]

T_01 = transformationMatrix(DH(1,:));
T_12 = transformationMatrix(DH(2,:));
T_23 = transformationMatrix(DH(3,:));
T_34 = transformationMatrix(DH(4,:));
T_45 = transformationMatrix(DH(5,:));
T_56 = transformationMatrix(DH(6,:));
T_6T = [    1       0       0       -0.1;
            0       1       0       0;
            0       0       1       0.13625;
            0       0       0       1]

T_06 = T_01*T_12*T_23*T_34*T_45*T_56
T_0T = T_06*T_6T
T_T6 = (T_6T)^(-1)

% Intial Configuration
T_initial = FKPuma([0 0 pi 0 0 0])

 FKPuma([0 0 0 0 pi/2 0])

T0_temp = FKPuma([0 -pi/2 3.5*pi/4 0 0 0])
[R0_temp, P0_temp] = tr2rt(T0_temp);
R0 = [  0 1 0;
        0 0 1;
        1 0 0];
T0 = rt2tr(R0, P0_temp)

IKPuma(T0, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331)

% This the common corner point of the cube. This will be considered
% origin to find other points
T0_corner = T0*T_6T


[R0_corner, P0_corner] = tr2rt(T0_corner);
P_now = [0; 0; 0];
theta1 = IKPuma2(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 1)

P = [];
positionError = [];
orientationError = [];

% First Letter (S)
Sxy = [ 0.075, 0.075;
        0.05, 0.1;
        0.025, 0.075;
        0.05, 0.05;
        0.075, 0.025;
        0.05, 0;
        0.025, 0.025];

ax_lim_s = [-0.02 0.12 -0.02 0.12];

[sx, sy] = plotLetter(Sxy, 200, ax_lim_s);

% Phase 1: Corner Point to Start of First Letter Transition 1
Tr1xy = [0, 0;
         sx(1), sy(1)];
[tr1x, tr1y] = plotLines(Tr1xy, 50);

theta1 = [];
theta0 = [0 0 pi 0 0 0];
theta_corner = IKPuma3([0; 0; 0], P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 1)
theta1 = [theta1; JSTrajectory2(theta0, theta_corner, 50)];

for i=1:length(tr1x)
    Px = tr1x(i); Py = tr1y(i); Pz = 0;
    P_now = [Px; Py; Pz];
    theta_now = IKPuma3(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 1);
    theta1 = [theta1; theta_now];
    P = [P P_now];
    [positionErrorNow, orientationErrorNow] = calculatePositionError(theta_now, P_now, P0_corner, 1);
    positionError = [positionError, positionErrorNow];
    orientationError = [orientationError, orientationErrorNow];
end

figure
time1 = 10;
time = plotThetas(theta1, time1, 1);

figure
plotxyz(theta1, time1, 1);

% Phase 2 : Trace the first letter
theta2 = [];
for i=1:length(sx)
    Px = sx(i); Py = sy(i); Pz = 0;
    P_now = [Px; Py; Pz];
    theta_now = IKPuma3(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 1);
    theta2 = [theta2; theta_now];
    P = [P P_now];
    [positionErrorNow, orientationErrorNow] = calculatePositionError(theta_now, P_now, P0_corner, 1);
    positionError = [positionError, positionErrorNow];
    orientationError = [orientationError, orientationErrorNow];
end

% Plot all joint angles of phase 2
figure
time2 = 20;
time_prev = time(end);
time = [time, plotThetas(theta2, time2, 2) + time_prev];

figure
plotxyz(theta2, time2, 2);

% Second Letter: B
Bxy_straight_1 = [ 0.1, -0.05; 
                    0.1, 0;
                    0.05, 0;
                    0.025, 0];
Bxy_spline = [   0.025, 0
                 0, -0.025;
                 0.025, -0.05;
                 0, -0.075;
                 0.025, -0.1];
Bxy_straight_2 = [ 0.025, -0.1;
                    0.1, -0.1;
                    0.1, -0.05];


figure
ax_lim_b = [-0.02 0.12 -0.12 0.02];

[bx_1, by_1] = plotLines(Bxy_straight_1, 50);
[bx_2, by_2] = plotLetter(Bxy_spline, 100, ax_lim_b);
[bx_3, by_3] = plotLines(Bxy_straight_2, 50);
bx = [bx_1, bx_2, bx_3];
by = [by_1, by_2, by_3];
figure
plotScript(bx, by, ax_lim_b);

%Phase 3: Transition from letter S to B
Tr3xy1 = [sx(end), sy(end);
         0, 0];
[tr3x1, tr3y1] = plotLines(Tr3xy1, 50);

theta3 = [];
for i=1:length(tr3x1)
    Px = tr3x1(i); Py = tr3y1(i); Pz = 0;
    P_now = [Px; Py; Pz];
    theta_now = IKPuma3(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 1);
    theta3 = [theta3; theta_now];
    P = [P P_now];
    [positionErrorNow, orientationErrorNow] = calculatePositionError(theta_now, P_now, P0_corner, 1);
    positionError = [positionError, positionErrorNow];
    orientationError = [orientationError, orientationErrorNow];
end

theta3 = [theta3; JSTrajectory1(1, 2, P0_corner, 50)];

figure
plotEulerAngles(1,2, 10/3); 
% Euler angles are plotted only where euler angles are changing. Rest
% everywhere they are 0.

Tr3xy2 = [0, 0;
         bx(1), by(1)];
[tr3x2, tr3y2] = plotLines(Tr3xy2, 50);

for i=1:length(tr3x2)
    Px = 0; Py = tr3x2(i) ; Pz = tr3y2(i);
    P_now = [Px; Py; Pz];
    theta_now = IKPuma3(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 2);
    theta3 = [theta3; theta_now];
    P = [P P_now];
    [positionErrorNow, orientationErrorNow] = calculatePositionError(theta_now, P_now, P0_corner, 2);
    positionError = [positionError, positionErrorNow];
    orientationError = [orientationError, orientationErrorNow];
end

figure
time_prev = time(end);
time = [time, plotThetas(theta3, 10, 3)+time_prev];

figure
plotxyz(theta3, 10, 3);

% Phase 4 : Trace the second letter B
theta4 = [];
for i=1:length(bx)
    Px = 0; Py = bx(i); Pz = by(i);
    P_now = [Px; Py; Pz];
    theta_now = IKPuma3(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 2);
    theta4 = [theta4; theta_now];
    P = [P P_now];
    [positionErrorNow, orientationErrorNow] = calculatePositionError(theta_now, P_now, P0_corner, 2);
    positionError = [positionError, positionErrorNow];
    orientationError = [orientationError, orientationErrorNow];
end

figure
time_prev = time(end);
time = [time, plotThetas(theta4, 20, 4)+time_prev];

figure
plotxyz(theta4, 20, 4);

% Third letter: J
Jxy_spline = [0, -0.05;
                0.025, -0.1;
                0.05, -0.05];

Jxy_straight_1 = [0.05, -0.05;
                0.05, 0];

Jxy_straight_2 = [0.1, 0;
                    0, 0];


figure
ax_lim_p = [-0.02 0.12 -0.12 0.02];

[px_1, py_1] = plotLetter(Jxy_spline, 50, ax_lim_p);
[px_2, py_2] = plotLines(Jxy_straight_1, 50);
[px_3, py_3] = plotLines(Jxy_straight_2, 50);
px = [px_1, px_2, px_3];
py = [py_1, py_2, py_3];
figure
plotScript(px, py, ax_lim_p);

%Phase 5: Transition from letter B to J
Tr5xy1 = [bx(end), by(end);
         0, 0];
[tr5x1, tr5y1] = plotLines(Tr5xy1, 50);

theta5 = [];
for i=1:length(tr5x1)
    Px = 0; Py = tr5x1(i) ; Pz = tr5y1(i);
    P_now = [Px; Py; Pz];
    theta_now = IKPuma3(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 2);
    theta5 = [theta5; theta_now];
    P = [P P_now];
    [positionErrorNow, orientationErrorNow] = calculatePositionError(theta_now, P_now, P0_corner, 2);
    positionError = [positionError, positionErrorNow];
    orientationError = [orientationError, orientationErrorNow];
end

theta5 = [theta5; JSTrajectory1(2, 3, P0_corner, 50)];

figure
plotEulerAngles(2,3, 10/3);

Tr5xy2 = [0, 0;
         px(1), py(1)];
[tr5x2, tr5y2] = plotLines(Tr5xy2, 50);

for i=1:length(tr5x2)
    Px = tr5x2(i); Py = 0; Pz = tr5y2(i);
    P_now = [Px; Py; Pz];
    theta_now = IKPuma3(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 3);
    theta5 = [theta5; theta_now];
    P = [P P_now];
    [positionErrorNow, orientationErrorNow] = calculatePositionError(theta_now, P_now, P0_corner, 3);
    positionError = [positionError, positionErrorNow];
    orientationError = [orientationError, orientationErrorNow];
end

figure
time_prev = time(end);
time = [time, plotThetas(theta5, 10, 5)+time_prev];

figure
plotxyz(theta5, 10, 5);

% Phase 6 : Trace the third letter J
theta6 = [];
for i=1:length(px)
    Px = px(i); Py = 0; Pz = py(i);
    P_now = [Px; Py; Pz];
    theta_now = IKPuma3(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 3);
    theta6 = [theta6; theta_now];
    P = [P P_now];
    [positionErrorNow, orientationErrorNow] = calculatePositionError(theta_now, P_now, P0_corner, 3);
    positionError = [positionError, positionErrorNow];
    orientationError = [orientationError, orientationErrorNow];
end

figure
time_prev = time(end);
time = [time, plotThetas(theta6, 20, 6)+time_prev];

figure
plotxyz(theta6, 20, 6);

% Phase 7: Letter J to the Corner
Tr7xy1 = [px(end), py(end);
         0, 0];
[tr7x1, tr7y1] = plotLines(Tr7xy1, 50);

theta7 = [];
for i=1:length(tr7x1)
    Px = tr7x1(i); Py = 0; Pz = tr7y1(i);
    P_now = [Px; Py; Pz];
    theta_now = IKPuma3(P_now, P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, 3);
    theta7 = [theta7; theta_now];
    P = [P P_now];
    [positionErrorNow, orientationErrorNow] = calculatePositionError(theta_now, P_now, P0_corner, 3);
    positionError = [positionError, positionErrorNow];
    orientationError = [orientationError, orientationErrorNow];
end

figure
theta7 = [theta7; JSTrajectory1(3, 1, P0_corner, 50)];

theta7 = [theta7; JSTrajectory2(theta_corner, theta0, 50)];

time_prev = time(end);
time = [time, plotThetas(theta7, 10, 7)+time_prev];

figure
plotxyz(theta7, 10, 7);

figure
P(1,:) = P(1,:)+P0_corner(1); 
P(2,:) = P(2,:)+P0_corner(2);
P(3,:) = P(3,:)+P0_corner(3);
scatter3(P(1,:),P(2,:),P(3,:))
axis([-0.02+P0_corner(1) 0.12+P0_corner(1) -0.02+P0_corner(2) 0.12+P0_corner(2) -0.12+P0_corner(3) 0.02+P0_corner(3)])
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
title('XYZ Trajectory in 3D')

figure
plot(positionError)
ylabel('Position Error')

figure
plot(orientationError)
ylabel('Orientation Error')



%% User Functions
function [T_06] = FKPuma(theta)
    
    % Numercial Expressions

        DH = [      0        0       0       theta(1);     %alpha, a, d, theta
                    -pi/2     0       0.2435      theta(2);
                    0       0.4318      -0.0934      theta(3);
                    pi/2    -0.0203      0.4331      theta(4);
                    -pi/2   0       0       theta(5);
                    pi/2    0       0       theta(6)
                ];
        
        T_01 = transformationMatrix(DH(1,:));
        T_12 = transformationMatrix(DH(2,:));
        T_23 = transformationMatrix(DH(3,:));
        T_34 = transformationMatrix(DH(4,:));
        T_45 = transformationMatrix(DH(5,:));
        T_56 = transformationMatrix(DH(6,:));
        
        T_06 = T_01*T_12*T_23*T_34*T_45*T_56;
    

end

function [theta] = IKPuma(M, a2, a3, d2, d3, d4 )
    
    Px = M(1,4); Py = M(2,4); Pz = M(3,4);

    R11 = M(1,1); R12 = M(1,2); R13 = M(1,3);
    R21 = M(2,1); R22 = M(2,2); R23 = M(2,3);
    R31 = M(3,1); R32 = M(3,2); R33 = M(3,3);
    % theta1
    t1s1 = atan2((Px^2 + Py^2 - (d2 + d3)^2)^0.5, d2 + d3) + atan2(-Px,Py);
    t1s2 = atan2(-(Px^2 + Py^2 - (d2 + d3)^2)^0.5, d2 + d3) + atan2(-Px,Py);

    %theta3 with t1 as t1s1
    t1 = t1s1;
    K3 = ((Px*cos(t1)+Py*sin(t1))^2 + Pz^2 - (a2^2 + a3^2 + d4^2))/(2*a2);
    t3s1 = atan2((a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s2 = atan2(-(a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s1 = wrapToPi(t3s1);
    t3s2 = wrapToPi(t3s2);

    %theta2 with t1s1 and t3s1
    t1 = t1s1;
    t3 = t3s1;
    A = Px*cos(t1)+Py*sin(t1);  C = -Pz;    D = a2 + a3*cos(t3) + d4*sin(t3);
    E = -Pz;                F = -A;         G = a3*sin(t3) - d4*cos(t3);
    t2s1 = atan2(D*E-A*G, C*G-D*F);

    %theta2 with t1s1 and t3s2
    t1 = t1s1;
    t3 = t3s2;
    A = Px*cos(t1)+Py*sin(t1);  C = -Pz;    D = a2 + a3*cos(t3) + d4*sin(t3);
    E = -Pz;                F = -A;         G = a3*sin(t3) - d4*cos(t3);
    t2s2 = atan2(D*E-A*G, C*G-D*F);

    %theta3 with t1 as t1s2
    t1 = t1s2;
    K3 = ((Px*cos(t1)+Py*sin(t1))^2 + Pz^2 - (a2^2 + a3^2 + d4^2))/(2*a2);
    t3s3 = atan2((a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s4 = atan2(-(a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s3 = wrapToPi(t3s3);
    t3s4 = wrapToPi(t3s4);

    %theta2 with t1s2 and t3s3
    t1 = t1s2;
    t3 = t3s3;
    A = Px*cos(t1)+Py*sin(t1);  C = -Pz;    D = a2 + a3*cos(t3) + d4*sin(t3);
    E = -Pz;                F = -A;         G = a3*sin(t3) - d4*cos(t3);
    t2s3 = atan2(D*E-A*G, C*G-D*F);
    
    %theta2 with t1s2 and t3s4
    t1 = t1s2;
    t3 = t3s4;
    A = Px*cos(t1)+Py*sin(t1);  C = -Pz;    D = a2 + a3*cos(t3) + d4*sin(t3);
    E = -Pz;                F = -A;         G = a3*sin(t3) - d4*cos(t3);
    t2s4 = atan2(D*E-A*G, C*G-D*F);
    
    theta123 = [t1s1, t2s1, t3s1;
            t1s1, t2s2, t3s2;
            t1s2, t2s3, t3s3;
            t1s2, t2s4, t3s4];

    %theta4 with t1s1, t2s1, t3s1 
    t1 = t1s1;  t2 = t2s1;  t3 = t3s1;
    t4s1 = atan2(R23*cos(t1)-R13*sin(t1), R23*sin(t1)*cos(t2+t3) - R13*cos(t1)*cos(t2+t3) - R33*sin(t2+t3));
    t4s1 = wrapToHalfPi(t4s1);
    t4 = t4s1;
    t5s1 = atan2((R31*cos(t2+t3) + R21*sin(t1)*sin(t2+t3) + R11*cos(t1)*sin(t2+t3)), R31*sin(t2+t3)*cos(t4) - R21*(cos(t1)*sin(t4)...
        + sin(t1)*cos(t4)*cos(t2+t3)) - R11*(-sin(t1)*sin(t4) + cos(t1)*cos(t4)*cos(t2+t3)) );
    t6s1 = atan2((R21*cos(t1)*cos(t4) - R21*sin(t1)*sin(t4)*cos(t2+t3) - R11*sin(t1)*cos(t4) - R11*cos(t1)*sin(t4)*cos(t2+t3) + R31*sin(t2+t3)*sin(t4)),...
        (R22*cos(t1)*cos(t4) - R22*sin(t1)*sin(t4)*cos(t2+t3) - R12*sin(t1)*cos(t4) - R12*cos(t1)*sin(t4)*cos(t2+t3) + R32*sin(t2+t3)*sin(t4))  );
    
    
    t5s1 = wrapToHalfPi(t5s1);
    t6s1 = wrapToHalfPi(t6s1);

    %theta4 with t1s1, t2s2, t3s2
    t1 = t1s1; t2 = t2s2; t3 = t3s2;
    t4s2 = atan2(R23*cos(t1)-R13*sin(t1), R23*sin(t1)*cos(t2+t3) - R13*cos(t1)*cos(t2+t3) - R33*sin(t2+t3));
    t4s2 = wrapToHalfPi(t4s2);
    t4 = t4s2;
    t5s2 = atan2((R31*cos(t2+t3) + R21*sin(t1)*sin(t2+t3) + R11*cos(t1)*sin(t2+t3)), R31*sin(t2+t3)*cos(t4) - R21*(cos(t1)*sin(t4)...
        + sin(t1)*cos(t4)*cos(t2+t3)) - R11*(-sin(t1)*sin(t4) + cos(t1)*cos(t4)*cos(t2+t3)) );
    t6s2 = atan2((R21*cos(t1)*cos(t4) - R21*sin(t1)*sin(t4)*cos(t2+t3) - R11*sin(t1)*cos(t4) - R11*cos(t1)*sin(t4)*cos(t2+t3) + R31*sin(t2+t3)*sin(t4)),...
        (R22*cos(t1)*cos(t4) - R22*sin(t1)*sin(t4)*cos(t2+t3) - R12*sin(t1)*cos(t4) - R12*cos(t1)*sin(t4)*cos(t2+t3) + R32*sin(t2+t3)*sin(t4))  );
    
    
    t5s2 = wrapToHalfPi(t5s2);
    t6s2 = wrapToHalfPi(t6s2);
    

    %theta4 with t1s2, t2s3, t3s3
    t1 = t1s2; t2 = t2s3; t3 = t3s3;
    t4s3 = atan2(R23*cos(t1)-R13*sin(t1), R23*sin(t1)*cos(t2+t3) - R13*cos(t1)*cos(t2+t3) - R33*sin(t2+t3));
    t4s3 = wrapToHalfPi(t4s3);
    t4 = t4s3;
    t5s3 = atan2((R31*cos(t2+t3) + R21*sin(t1)*sin(t2+t3) + R11*cos(t1)*sin(t2+t3)), R31*sin(t2+t3)*cos(t4) - R21*(cos(t1)*sin(t4)...
        + sin(t1)*cos(t4)*cos(t2+t3)) - R11*(-sin(t1)*sin(t4) + cos(t1)*cos(t4)*cos(t2+t3)) );
    t6s3 = atan2((R21*cos(t1)*cos(t4) - R21*sin(t1)*sin(t4)*cos(t2+t3) - R11*sin(t1)*cos(t4) - R11*cos(t1)*sin(t4)*cos(t2+t3) + R31*sin(t2+t3)*sin(t4)),...
        (R22*cos(t1)*cos(t4) - R22*sin(t1)*sin(t4)*cos(t2+t3) - R12*sin(t1)*cos(t4) - R12*cos(t1)*sin(t4)*cos(t2+t3) + R32*sin(t2+t3)*sin(t4))  );
    
    
    t5s3 = wrapToHalfPi(t5s3);
    t6s3 = wrapToHalfPi(t6s3);


    %theta4 with t1s2, t2s4, t3s4
    t1 = t1s2; t2 = t2s4; t3 = t3s4;
    t4s4 = atan2(R23*cos(t1)-R13*sin(t1), R23*sin(t1)*cos(t2+t3) - R13*cos(t1)*cos(t2+t3) - R33*sin(t2+t3));
    t4s4 = wrapToHalfPi(t4s4);
    t4 = t4s4;
    t5s4 = atan2((R31*cos(t2+t3) + R21*sin(t1)*sin(t2+t3) + R11*cos(t1)*sin(t2+t3)), (R31*sin(t2+t3)*cos(t4) - R21*(cos(t1)*sin(t4)...
        + sin(t1)*cos(t4)*cos(t2+t3)) - R11*(-sin(t1)*sin(t4) + cos(t1)*cos(t4)*cos(t2+t3))) );
    t6s4 = atan2((R21*cos(t1)*cos(t4) - R21*sin(t1)*sin(t4)*cos(t2+t3) - R11*sin(t1)*cos(t4) - R11*cos(t1)*sin(t4)*cos(t2+t3) + R31*sin(t2+t3)*sin(t4)),...
        (R22*cos(t1)*cos(t4) - R22*sin(t1)*sin(t4)*cos(t2+t3) - R12*sin(t1)*cos(t4) - R12*cos(t1)*sin(t4)*cos(t2+t3) + R32*sin(t2+t3)*sin(t4))  );
    
    
    t5s4 = wrapToHalfPi(t5s4);
    t6s4 = wrapToHalfPi(t6s4);

    theta = [t1s1, t2s1, t3s1, t4s1, t5s1, t6s1;
            t1s1, t2s2, t3s2, t4s2, t5s2, t6s2;
            t1s2, t2s3, t3s3, t4s3, t5s3, t6s3;
            t1s2, t2s4, t3s4, t4s4, t5s4, t6s4];

end

function [theta] = IKPuma2(P, P0_corner, a2, a3, d2, d3, d4, O )
    
    P = P0_corner + P;

    if (O==1)
        M = [   0   1   0   P(1);
                0   0   1   P(2);
                1   0   0   P(3)
                0   0   0   1];
    elseif (O==2)
        M = [   -1   0   0   P(1);
                0   1   0   P(2);
                0   0   -1   P(3)
                0   0   0   1];
    elseif (O==3)
        M = [   0   0   -1   P(1);
                -1   0   0   P(2);
                0   1   0   P(3)
                0   0   0   1];
    end

    T_6T = [    1       0       0       -0.1;
            0       1       0       0;
            0       0       1       0.13625;
            0       0       0       1];
    
    
    M = M*(T_6T^(-1));

    Px = M(1,4); Py = M(2,4); Pz = M(3,4);

    R11 = M(1,1); R12 = M(1,2); R13 = M(1,3);
    R21 = M(2,1); R22 = M(2,2); R23 = M(2,3);
    R31 = M(3,1); R32 = M(3,2); R33 = M(3,3);



    % theta1
    t1s1 = atan2((Px^2 + Py^2 - (d2 + d3)^2)^0.5, d2 + d3) + atan2(-Px,Py);
    t1s2 = atan2(-(Px^2 + Py^2 - (d2 + d3)^2)^0.5, d2 + d3) + atan2(-Px,Py);
    if(abs(t1s1) <= abs(t1s2))
        t1 = t1s1;
    else
        t1 = t1s2;
    end
    
    % theta3
    K3 = ((Px*cos(t1)+Py*sin(t1))^2 + Pz^2 - (a2^2 + a3^2 + d4^2))/(2*a2);
    t3s1 = atan2((a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s2 = atan2(-(a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s1 = wrapToPi(t3s1);
    t3s2 = wrapToPi(t3s2);
    if (abs(t3s1)>=pi/2)
        t3 = t3s1;
    else
        t3 = t3s2;
    end

    % theta2
    A = Px*cos(t1)+Py*sin(t1);  C = -Pz;    D = a2 + a3*cos(t3) + d4*sin(t3);
    E = -Pz;                F = -A;         G = a3*sin(t3) - d4*cos(t3);
    t2 = atan2(D*E-A*G, C*G-D*F);

   

    %theta4, theta5, theta6
    t4 = atan2(R23*cos(t1)-R13*sin(t1), R23*sin(t1)*cos(t2+t3) - R13*cos(t1)*cos(t2+t3) - R33*sin(t2+t3));
    t4 = wrapToHalfPi(t4);
    t5 = atan2((R31*cos(t2+t3) + R21*sin(t1)*sin(t2+t3) + R11*cos(t1)*sin(t2+t3)), R31*sin(t2+t3)*cos(t4) - R21*(cos(t1)*sin(t4)...
        + sin(t1)*cos(t4)*cos(t2+t3)) - R11*(-sin(t1)*sin(t4) + cos(t1)*cos(t4)*cos(t2+t3)) );
    t6 = atan2((R21*cos(t1)*cos(t4) - R21*sin(t1)*sin(t4)*cos(t2+t3) - R11*sin(t1)*cos(t4) - R11*cos(t1)*sin(t4)*cos(t2+t3) + R31*sin(t2+t3)*sin(t4)),...
        (R22*cos(t1)*cos(t4) - R22*sin(t1)*sin(t4)*cos(t2+t3) - R12*sin(t1)*cos(t4) - R12*cos(t1)*sin(t4)*cos(t2+t3) + R32*sin(t2+t3)*sin(t4))  );
    
    
    t5 = wrapToHalfPi(t5);
    t6 = wrapToHalfPi(t6);


    theta = [t1, t2, t3, t4, t5, t6];

end

function [theta] = IKPuma3(P, P0_corner, a2, a3, d2, d3, d4, O )
    
    P = P0_corner + P;

    if (O==1)
        M = [   0   1   0   P(1);
                0   0   1   P(2);
                1   0   0   P(3)
                0   0   0   1];
    elseif (O==2)
        M = [   -1   0   0   P(1);
                0   1   0   P(2);
                0   0   -1   P(3)
                0   0   0   1];
    elseif (O==3)
        M = [   0   0   -1   P(1);
                -1   0   0   P(2);
                0   1   0   P(3)
                0   0   0   1];
    end

    T_6T = [    1       0       0       -0.1;
            0       1       0       0;
            0       0       1       0.13625;
            0       0       0       1];
    
    
    M = M*(T_6T^(-1));

    Px = M(1,4); Py = M(2,4); Pz = M(3,4);

    R11 = M(1,1); R12 = M(1,2); R13 = M(1,3);
    R21 = M(2,1); R22 = M(2,2); R23 = M(2,3);
    R31 = M(3,1); R32 = M(3,2); R33 = M(3,3);



    % theta1
    t1s1 = atan2((Px^2 + Py^2 - (d2 + d3)^2)^0.5, d2 + d3) + atan2(-Px,Py);
    t1s2 = atan2(-(Px^2 + Py^2 - (d2 + d3)^2)^0.5, d2 + d3) + atan2(-Px,Py);
    if(abs(t1s1) <= abs(t1s2))
        t1 = t1s1;
    else
        t1 = t1s2;
    end
    
    % theta3
    K3 = ((Px*cos(t1)+Py*sin(t1))^2 + Pz^2 - (a2^2 + a3^2 + d4^2))/(2*a2);
    t3s1 = atan2((a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s2 = atan2(-(a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s1 = wrapToPi(t3s1);
    t3s2 = wrapToPi(t3s2);
    if (abs(t3s1)>=pi/2)
        t3 = t3s1;
    else
        t3 = t3s2;
    end

    % theta2
    A = Px*cos(t1)+Py*sin(t1);  C = -Pz;    D = a2 + a3*cos(t3) + d4*sin(t3);
    E = -Pz;                F = -A;         G = a3*sin(t3) - d4*cos(t3);
    t2 = atan2(D*E-A*G, C*G-D*F);

   

    %theta4, theta5, theta6
    t4 = atan2(R23*cos(t1)-R13*sin(t1), R23*sin(t1)*cos(t2+t3) - R13*cos(t1)*cos(t2+t3) - R33*sin(t2+t3));
%     t4 = wrapToHalfPi(t4);
    t5 = atan2((R31*cos(t2+t3) + R21*sin(t1)*sin(t2+t3) + R11*cos(t1)*sin(t2+t3)), R31*sin(t2+t3)*cos(t4) - R21*(cos(t1)*sin(t4)...
        + sin(t1)*cos(t4)*cos(t2+t3)) - R11*(-sin(t1)*sin(t4) + cos(t1)*cos(t4)*cos(t2+t3)) );
    t6 = atan2((R21*cos(t1)*cos(t4) - R21*sin(t1)*sin(t4)*cos(t2+t3) - R11*sin(t1)*cos(t4) - R11*cos(t1)*sin(t4)*cos(t2+t3) + R31*sin(t2+t3)*sin(t4)),...
        (R22*cos(t1)*cos(t4) - R22*sin(t1)*sin(t4)*cos(t2+t3) - R12*sin(t1)*cos(t4) - R12*cos(t1)*sin(t4)*cos(t2+t3) + R32*sin(t2+t3)*sin(t4))  );
    
    
%     t5 = wrapToHalfPi(t5);
%     t6 = wrapToHalfPi(t6);


    theta = [t1, t2, t3, t4, t5, t6];

end

function [theta] = IKPuma4(M, a2, a3, d2, d3, d4 )
    

    T_6T = [    1       0       0       -0.1;
            0       1       0       0;
            0       0       1       0.13625;
            0       0       0       1];
    
    
    M = M*(T_6T^(-1));

    Px = M(1,4); Py = M(2,4); Pz = M(3,4);

    R11 = M(1,1); R12 = M(1,2); R13 = M(1,3);
    R21 = M(2,1); R22 = M(2,2); R23 = M(2,3);
    R31 = M(3,1); R32 = M(3,2); R33 = M(3,3);



    % theta1
    t1s1 = atan2((Px^2 + Py^2 - (d2 + d3)^2)^0.5, d2 + d3) + atan2(-Px,Py);
    t1s2 = atan2(-(Px^2 + Py^2 - (d2 + d3)^2)^0.5, d2 + d3) + atan2(-Px,Py);
    if(abs(t1s1) <= abs(t1s2))
        t1 = t1s1;
    else
        t1 = t1s2;
    end
    
    % theta3
    K3 = ((Px*cos(t1)+Py*sin(t1))^2 + Pz^2 - (a2^2 + a3^2 + d4^2))/(2*a2);
    t3s1 = atan2((a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s2 = atan2(-(a3^2 + d4^2 - K3^2)^0.5, K3) + atan2(d4, a3);
    t3s1 = wrapToPi(t3s1);
    t3s2 = wrapToPi(t3s2);
    if (abs(t3s1)>=pi/2)
        t3 = t3s1;
    else
        t3 = t3s2;
    end

    % theta2
    A = Px*cos(t1)+Py*sin(t1);  C = -Pz;    D = a2 + a3*cos(t3) + d4*sin(t3);
    E = -Pz;                F = -A;         G = a3*sin(t3) - d4*cos(t3);
    t2 = atan2(D*E-A*G, C*G-D*F);

   

    %theta4, theta5, theta6
    t4 = atan2(R23*cos(t1)-R13*sin(t1), R23*sin(t1)*cos(t2+t3) - R13*cos(t1)*cos(t2+t3) - R33*sin(t2+t3));
%     t4 = wrapToHalfPi(t4);
    t5 = atan2((R31*cos(t2+t3) + R21*sin(t1)*sin(t2+t3) + R11*cos(t1)*sin(t2+t3)), R31*sin(t2+t3)*cos(t4) - R21*(cos(t1)*sin(t4)...
        + sin(t1)*cos(t4)*cos(t2+t3)) - R11*(-sin(t1)*sin(t4) + cos(t1)*cos(t4)*cos(t2+t3)) );
    t6 = atan2((R21*cos(t1)*cos(t4) - R21*sin(t1)*sin(t4)*cos(t2+t3) - R11*sin(t1)*cos(t4) - R11*cos(t1)*sin(t4)*cos(t2+t3) + R31*sin(t2+t3)*sin(t4)),...
        (R22*cos(t1)*cos(t4) - R22*sin(t1)*sin(t4)*cos(t2+t3) - R12*sin(t1)*cos(t4) - R12*cos(t1)*sin(t4)*cos(t2+t3) + R32*sin(t2+t3)*sin(t4))  );
    
    
%     t5 = wrapToHalfPi(t5);
%     t6 = wrapToHalfPi(t6);


    theta = [t1, t2, t3, t4, t5, t6];

end

function [theta] = JSTrajectory1(O1, O2, P0_corner, n)
    
     thetaA = IKPuma3([0 0 0], P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, O1);
     thetaB = IKPuma3([0 0 0], P0_corner, 0.4318, -0.0203, 0.2435, -0.0934, 0.4331, O2);
     

%      for i=1:6
%          theta(:,i) = linspace(thetaA(i), thetaB(i), n);
%      end
%     vmax = [8; 10; 10; 5; 5; 5];
    [m10, m11, m12, m13] = cubicTrajectory(thetaA', thetaB', n);
    t1 = linspace(0, n, n);
    theta = m10 + m11*t1 + m12*t1.^2 + m13*t1.^3;
    theta = theta';
end

function [theta] = JSTrajectory2(thetaA, thetaB, n)
%     for i=1:6
%          theta(:,i) = linspace(thetaA(i), thetaB(i), n);
%     end
    [m10, m11, m12, m13] = cubicTrajectory(thetaA', thetaB', n);
    t1 = linspace(0, n, n);
    theta = m10 + m11*t1 + m12*t1.^2 + m13*t1.^3;
    theta = theta';
end

function [t] = plotEulerAngles(O1,O2, time)
    t = linspace(0, time, 50);
    if O1 == 1 && O2 ==2
        alpha = (pi/50)*ones(1, 50);
        beta = 0;
        gamma = ((pi/2)/50)*ones(1, 50);
    end
    if O1 ==2 && O2 == 3
        alpha = ((-pi/2)/50)*ones(1, 50);
        beta = 0;
        gamma = ((-pi/2)/50)*ones(1, 50);
    end
    
    figure
    hold on
    plot(t, alpha);
    plot(t, beta);
    plot(t, gamma);
    hold off
    legend('$\alpha$','$\beta$','$\gamma$', 'interpreter', 'latex')
    xlabel('Time(s)')
    ylabel('Euler Angles (rad)')
end

function [lx, ly] = plotLetter(Lxy, n, ax_lim )
    xxa=Lxy(:,1);
    yya=Lxy(:,2);
    distF=[0 ;sqrt(sum(diff([xxa,yya]).^2,2))];
    distFSum=cumsum(distF);
    t = linspace(min(distFSum),max(distFSum),n);
    lx=spline(distFSum,xxa,t);
    ly=spline(distFSum,yya,t);
    plot(lx,ly,'y','LineWidth',2)
    hold on
    plot(xxa,yya,'ko')
    axis(ax_lim)
    hold off
end

function [lx, ly] = plotLines(Lxy, n)
    xxa=Lxy(:,1);
    yya=Lxy(:,2);
    lx = [];
    ly = [];
    for i=1:length(xxa)-1
        xxa_intermediate = linspace(xxa(i), xxa(i+1), n);
        yya_intermediate = linspace(yya(i), yya(i+1), n);
        lx = [lx, xxa_intermediate];
        ly = [ly, yya_intermediate];
    end
end

function [done] = plotScript(lx, ly, ax_lim)
    plot(lx,ly,'y','LineWidth',2)
    axis(ax_lim)
    done= 1;
end

function [t] = plotThetas(theta, time, phase)
    
t = linspace(0, time, length(theta(:,1)));
figure
hold on
plot(t, theta(:, 1));
plot(t, theta(:, 2));
plot(t, theta(: ,3));
plot(t, theta(: ,4));
plot(t, theta(:, 5));
plot(t, theta(:, 6));
hold off
legend('$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','$\theta_5$','$\theta_6$', 'interpreter', 'latex')
xlabel('Time(s)')
ylabel('Joint Angles (rad)')
title(sprintf('Joint Angles of Phase: %d', phase) )
fprintf('The amount of via points used in this Phase is: %d', length(theta(:,1)))
end

function [t] = plotxyz(theta, time, phase)
    t = linspace(0, time, length(theta(:,1)));
     T_6T = [    1       0       0       -0.1;
            0       1       0       0;
            0       0       1       0.13625;
            0       0       0       1];
     x = []; y = []; z = [];
    for i=1:length(theta(:,1))
        T_actual = FKPuma(theta(i,:));
%         T_tool = T_actual*T_6T;
%         T_tool = T_actual;
        T_tool = T_actual;
        T_tool(1,4) = T_tool(1,4) + T_6T(1,4);
        T_tool(2,4) = T_tool(2,4) + T_6T(2,4);
        T_tool(3,4) = T_tool(3,4) + T_6T(3,4);
        [R_tool, P_tool] = tr2rt(T_tool);
        x = [x, P_tool(1)];
        y = [y, P_tool(2)];
        z = [z, P_tool(3)];
    end

    figure
    hold on
    plot(t, x);
    plot(t, y);
    plot(t, z);
    hold off
    legend('x','y','z')
    xlabel('Time(s)')
    ylabel('x, y, z (m)')
    title(sprintf('X,Y, Z Positions of Tool Tip of Phase: %d', phase) )
end

function [posError, orientError] = calculatePositionError(theta, P_required, P0_corner, O)
    T_actual = FKPuma(theta);
    T_6T = [    1       0       0       -0.1;
            0       1       0       0;
            0       0       1       0.13625;
            0       0       0       1];
    
    P = P0_corner + P_required;

    if (O==1)
        M = [   0   1   0   P(1);
                0   0   1   P(2);
                1   0   0   P(3)
                0   0   0   1];
    elseif (O==2)
        M = [   -1   0   0   P(1);
                0   1   0   P(2);
                0   0   -1   P(3)
                0   0   0   1];
    elseif (O==3)
        M = [   0   0   -1   P(1);
                -1   0   0   P(2);
                0   1   0   P(3)
                0   0   0   1];
    end

    
    
    T_required = M*(T_6T^(-1));
    
    [R_required, P_required] = tr2rt(T_required);
    [R_actual, P_actual] = tr2rt(T_actual);
    posError = sqrt(sum((P_actual - P_required) .^ 2));
    orientError = sqrt(sum((R_actual - R_required) .^ 2, 'all'));
end

function [m0, m1, m2, m3] = cubicTrajectory(theta0, theta1, n)

    tf = n;
    m0 = theta0;
    m1 = 0;
    m2 = (3/tf^2)*(theta1 - theta0);
    m3 = (-2/tf^3)*(theta1 - theta0);
end

function [T] = transformationMatrix(DH_row)

 T = [cos(DH_row(4))                     -sin(DH_row(4))                    0                   DH_row(2);
      sin(DH_row(4))*cos(DH_row(1))      cos(DH_row(4))*cos(DH_row(1))      -sin(DH_row(1))     -sin(DH_row(1))*DH_row(3);
      sin(DH_row(4))*sin(DH_row(1))      cos(DH_row(4))*sin(DH_row(1))      cos(DH_row(1))      cos(DH_row(1))*DH_row(3);
      0                                  0                                  0                   1]; 
end

function [wrapped_angle] = wrapToHalfPi(lambda)
    
    tmp = mod(lambda+pi/2,pi);
    wrapped_angle = tmp+pi*(lambda>0&tmp==0)-pi/2;
end

