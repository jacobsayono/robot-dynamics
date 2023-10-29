%% Jacob Sayono
% 505368811
% MAE C163B Project 4

%% Equations of Motion -- Derivation
clear all; close all; clc

syms l1 l2 t1 t2 t3 m1 m2 dt1 dt2 ddt1 ddt2 g f3x f3y n3z

L(1) = Link('revolute','d', 0, 'a', 0, 'alpha', 0 ,'modified');
L(2) = Link('revolute','d', 0, 'a', l1, 'alpha', 0 ,'modified');
L(3) = Link('revolute','d', 0, 'a', l2, 'alpha', 0 ,'modified');

RR = SerialLink(L, 'name', 'RR-Manipulator');

%% Newton-Euler Formulation

th = [t1 t2 0]

T_01 = RR.A([1], th);
T_12 = RR.A([2], th);
T_2T = RR.A([3], th);
T_0T = RR.A([1 2 3], th);
T_0T = simplify(T_0T)

[R_01, P_01] = tr2rt(T_01); R_10 = transpose(R_01);
[R_12, P_12] = tr2rt(T_12); R_21 = transpose(R_12);
[R_2T, P_2T] = tr2rt(T_2T); R_32 = transpose(R_2T);
[R_0T, P_0T] = tr2rt(T_0T);
R_0T = simplify(R_0T)
P_0T = simplify(P_0T)

PC1 = [l1/2; 0 ; 0];
PC2 = [l2/2; 0 ; 0];

IC1 = (1/12) * m1 * l1^2 * [0 0 0; 0 1 0; 0 0 1];
IC2 = (1/12) * m2 * l2^2 * [0 0 0; 0 1 0; 0 0 1];

f3 = [f3x; f3y; 0]; 
n3 = [0;0;n3z]; 

w0 = zeros(3,1);
wd0 = zeros(3,1); 

v0 = zeros(3,1); 
vd0 = [0 ; 0 ; -g];

% Inward Iteration

% i = 0
w1 = R_10 * w0 + dt1*R_01(1:3,3)
wd1 = R_10 * wd0 + R_10 * cross(w0, dt1*R_01(1:3,3)) + ddt1*R_01(1:3,3)

vd1 = R_10 * (cross(wd0, P_01) + cross(w0, cross(w0, P_01)) + vd0)
vcd1 = cross(wd1,PC1) + cross(w1,cross(w1,PC1)) + vd1

F1 = m1 * vcd1 
N1 = IC1 * wd1 + cross(w1,IC1*w1)

% i = 1
w2 = R_21 * w1 + dt2*R_12(1:3,3)
wd2 = R_21 * wd1 + R_21 * cross(w1, dt2*R_12(1:3,3)) + ddt2*R_12(1:3,3)

vd2 = R_21 * (cross(wd1, P_12) + cross(w1, cross(w1, P_12)) + vd1)
vcd2 = cross(wd2,PC2) + cross(w2,cross(w2,PC2)) + vd2

F2 = m2 * vcd2 
N2 = IC2 * wd2 + cross(w2,IC2*w2)

% Outward Iteration

% i = 2
f2 = R_2T * f3 + F2;
n2 = N2 + R_2T*n3 + cross(PC2, F2) + cross(P_2T, R_2T*f3);
f2 = simplify(f2)
n2 = simplify(n2)

% i = 1
f1 = R_12 * f2 + F1;
n1 = N1 + R_12*n2 + cross(PC1, F1) + cross(P_12, R_12*f2);
f1 = simplify(f1)
n1 = simplify(n1)

%% Design Trajectory
clear all; close all; clc

l1 = 0.5;
l2 = 0.5;

L(1) = Link('revolute','d', 0, 'a', 0, 'alpha', 0 ,'modified');
L(2) = Link('revolute','d', 0, 'a', l1, 'alpha', 0 ,'modified');
L(3) = Link('revolute','d', 0, 'a', l2, 'alpha', 0 ,'modified');

RR = SerialLink(L, 'name', 'RR-Manipulator');

g = 0
rho = 1000;
r_outer = 0.1;
r_inner = 0.005;
m1 = rho*l1*pi*(r_outer^2 - r_inner^2);
m2 = rho*l2*pi*(r_outer^2 - r_inner^2);

PC1 = [l1/2; 0 ; 0];
PC2 = [l2/2; 0 ; 0];

Ix = 0.5*m1*(r_outer^2 + r_inner^2);
Iy = Ix/2 + (1/12)*m1*l1^2; Iz = Ix/2 + + (1/12)*m2*l2^2;

IC1 = [Ix 0 0; 0 Iy 0; 0 0 Iz]
IC2 = [Ix 0 0; 0 Iy 0; 0 0 Iz]

f3 = [-10;0;0]; 
n3 = [0;0;10]; 

w0 = zeros(3,1);
wd0 = zeros(3,1); 

vd0 = [0 ; 0 ; 0];

t1_initial = -acos(0.45/0.5)
t2_initial = 2*acos(0.45/0.5)
t1_final = -acos(0.05/0.5)
t2_final = 2*acos(0.05/0.5)

N = 100;
t1 = linspace(t1_initial, t1_final, N+2)
t2 = linspace(t2_initial, t2_final, N+2)
totalTime = 4;
dt = 4/N;
dt1 = (diff(t1))/dt
dt2 = (diff(t2))/dt
ddt1 = (diff(t1,2))/dt^2
ddt2 = (diff(t2,2))/dt^2

for j=1:N

    th = [t1(j) t2(j) 0];

    T_01 = RR.A([1], th);
    T_12 = RR.A([2], th);
    T_2T = RR.A([3], th);
    T_0T = RR.A([1 2 3], th);

    [R_01, P_01] = tr2rt(T_01); R_10 = transpose(R_01);
    [R_12, P_12] = tr2rt(T_12); R_21 = transpose(R_12);
    [R_2T, P_2T] = tr2rt(T_2T); R_32 = transpose(R_2T);
    [R_0T, P_0T] = tr2rt(T_0T);

    x(j) = P_0T(1); y(j) = P_0T(2);


    % i = 0
    w1 = R_10 * w0 + dt1(j)*R_01(1:3,3);
    wd1 = R_10 * wd0 + R_10 * cross(w0, dt1(j)*R_01(1:3,3)) + ddt1(j)*R_01(1:3,3);
    
    vd1 = R_10 * (cross(wd0, P_01) + cross(w0, cross(w0, P_01)) + vd0);
    vcd1 = cross(wd1,PC1) + cross(w1,cross(w1,PC1)) + vd1;
    
    F1 = m1 * vcd1 ;
    N1 = IC1 * wd1 + cross(w1,IC1*w1);
    
    % i = 1
    w2 = R_21 * w1 + dt2(j)*R_12(1:3,3);
    wd2 = R_21 * wd1 + R_21 * cross(w1, dt2(j)*R_12(1:3,3)) + ddt2(j)*R_12(1:3,3);
    
    vd2 = R_21 * (cross(wd1, P_12) + cross(w1, cross(w1, P_12)) + vd1);
    vcd2 = cross(wd2,PC2) + cross(w2,cross(w2,PC2)) + vd2;
    
    F2 = m2 * vcd2 ;
    N2 = IC2 * wd2 + cross(w2,IC2*w2);
    
    % i = 2
    f2 = R_2T * f3 + F2;
    n2(:,j) = N2 + R_2T*n3 + cross(PC2, F2) + cross(P_2T, R_2T*f3);
    
    % i = 1
    f1 = R_12 * f2 + F1;
    n1(:,j) = N1 + R_12*n2(:,j) + cross(PC1, F1) + cross(P_12, R_12*f2);
end

time = linspace(0, totalTime, N);

% Plot Trajectory
figure(1)
title('Trajectory Animation')
xlabel('X-Direction (m)')
ylabel('Y-Direction (m)')
h = animatedline;
axis([0 1 -0.5 0.5])
for j=1:N
    addpoints(h,x(j),y(j));
    drawnow
    pause(0.05)
end

% End Effector Position
figure(2)
plot(time, x);
title('EE Position vs. Time')
ylabel('X and Y Plane (m)')
xlabel('Time (s)')
hold on
plot(time, y);
hold off

% Joint Torque
figure(3)
plot(time, n1);
title('Joint Torque vs. Time')
ylabel('Joint Torque (Nm)')
xlabel('Time (s)')
hold on
plot(time, n2);
legend('Joint 1','Joint 2')

%% Inertia Tensor
clear all; close all; clc;

% Body A
h = 0.1; l = 0.1; w = 0.1; r = 0.1; d=0.4;
I_cube = [h^2+l^2 0 0; 0 w^2+h^2 0; 0 0 l^2+h^2];
I_cyl = [(1/12)*(3*r^2+h^2) 0 0; 0 (1/12)*(3*r^2+h^2) 0; 0 0 0.5^r^2];
I_Acm = I_cube - I_cyl;
I_A = I_Acm + ([-d 0 0]*[-d; 0; 0]*eye(3) - [d^2 0 0; 0 0 0; 0 0 0]) 

% Body B
rB = 0.05; lB = 0.8;
I_B = [0.5*rB^2 0 0; 0 (1/12)*(3*rB^2+lB^2) 0; 0 0 (1/12)*(3*rB^2+lB^2)]

% Body C
I_Ccm = I_Acm;
I_C1 = I_Ccm + ([d 0 0]*[d; 0; 0]*eye(3) - [d^2 0 0; 0 0 0; 0 0 0])
Rotx = [1 0 0; 0 1/2^0.5 1/2^0.5; 0 -1/2^0.5 1/2^0.5];
I_C = Rotx*I_C1*Rotx'

% Total Inertia
syms mA mB mC
I = mA*I_A + mB*I_B + mC*I_C;
I = vpa(I, 4)