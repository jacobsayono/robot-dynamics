%% MAE C163B Final Exam

% Jacob Sayono

% 505368811

%% Jacobian

clear all; close all; clc;

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
simplify(T_06)

[R_01, P_01] = tr2rt(T_01); R_10 = transpose(R_01);
[R_12, P_12] = tr2rt(T_12); R_21 = transpose(R_12);
[R_23, P_23] = tr2rt(T_23); R_32 = transpose(R_23);
[R_34, P_34] = tr2rt(T_34); R_43 = transpose(R_34);
[R_45, P_45] = tr2rt(T_45); R_54 = transpose(R_45);
[R_56, P_56] = tr2rt(T_56); R_65 = transpose(R_56);
[R_06, P_06] = tr2rt(T_06); R_60 = transpose(R_06);

T_02 = T_01*T_12;
T_03 = T_02*T_23;
T_04 = T_03*T_34;
T_05 = T_04*T_45;
T_06 = T_05*T_56; 
[R_02, P_02] = tr2rt(T_02); R_20 = transpose(R_02);
[R_03, P_03] = tr2rt(T_03); R_30 = transpose(R_03);
[R_04, P_04] = tr2rt(T_04); R_40 = transpose(R_04);
[R_05, P_05] = tr2rt(T_05); R_50 = transpose(R_05);

J0_DD = sym(zeros(6,3));
J0_DD(:,1) = [cross(R_01(:, 3), (P_06 - P_01));R_01(:,3)];
J0_DD(:,2) = [cross(R_02(:, 3), (P_06 - P_02));R_02(:,3)];
J0_DD(:,3) = [cross(R_03(:, 3), (P_06 - P_03));R_03(:,3)];
J0_DD(:,4) = [cross(R_04(:, 3), (P_06 - P_04));R_04(:,3)];
J0_DD(:,5) = [cross(R_05(:, 3), (P_06 - P_05));R_05(:,3)];
J0_DD(:,6) = [cross(R_06(:, 3), (P_06 - P_06));R_06(:,3)];
J0_DD = simplify(J0_DD)

J0_det = det(J0_DD)

simplify(J0_det)

%% Dynamics

% Symbolic Expressions
syms t1 t2 t3 t4 t5 t6 a2 a3 d2 d3 d4 
syms c1 c2 c3 Ix1 Ix2 Ix3 Iy1 Iy2 Iy3 Iz1 Iz2 Iz3 m1 m2 m3 m4;

DH = [      0        0       0       t1;     %alpha, a, d, theta
            -pi/2     0       d2      t2;
            0       a2      d3      t3;
            pi/2    a3      d4      0;
        ]

T_01 = transformationMatrix(DH(1,:));
T_12 = transformationMatrix(DH(2,:));
T_23 = transformationMatrix(DH(3,:));
T_34 = transformationMatrix(DH(4,:));

T_04 = T_01*T_12*T_23*T_34; 
T_04 = simplify(T_04);

[R_01, P_01] = tr2rt(T_01); R_10 = transpose(R_01);
[R_12, P_12] = tr2rt(T_12); R_21 = transpose(R_12);
[R_23, P_23] = tr2rt(T_23); R_32 = transpose(R_23);
[R_34, P_34] = tr2rt(T_34); R_43 = transpose(R_34);
[R_04, P_04] = tr2rt(T_04); R_40 = transpose(R_04);

PC1 = [0; d2/2 ; 0];
PC2 = [a2/2; 0 ; 0];
PC3 = [0 ; -d4/2; 0];

IC1 = (1/12)*m1*(d2^2)*[1 0 0; 0 0 0; 0 0 1];
IC2 = (1/12)*m2*(a2^2)*[0 0 0; 0 1 0; 0 0 1];
IC3 = (1/12)*m3*(d4^2)*[1 0 0; 0 0 0; 0 0 1]...
    +  (m4*(d4/2)^2)* [1 0 0; 0 0 0; 0 0 1];


syms f4x f4y f4z n4x n4y n4z g  dt1 dt2 dt3 ddt1 ddt2 ddt3  ;

f4 = [f4x; f4y; f4z]; 
n4 = [n4x; n4y; n4z]; 

w0 = zeros(3,1);
wd0 = zeros(3,1); 

v0 = zeros(3,1); 
vd0 = [0 ; 0 ; -g];

% Inward Iteration

% i = 0
w1 = R_10 * w0 + dt1*R_01(1:3,3);
wd1 = R_10 * wd0 + R_10 * cross(w0, dt1*R_01(1:3,3)) + ddt1*R_01(1:3,3);

vd1 = R_10 * (cross(wd0, P_01) + cross(w0, cross(w0, P_01)) + vd0);
vcd1 = cross(wd1,PC1) + cross(w1,cross(w1,PC1)) + vd1;

F1 = m1 * vcd1 ;
N1 = IC1 * wd1 + cross(w1,IC1*w1);


% i = 1
w2 = R_21 * w1 + dt2*R_12(1:3,3);
wd2 = R_21 * wd1 + R_21 * cross(w1, dt2*R_12(1:3,3)) + ddt2*R_12(1:3,3);

vd2 = R_21 * (cross(wd1, P_12) + cross(w1, cross(w1, P_12)) + vd1);
vcd2 = cross(wd2,PC2) + cross(w2,cross(w2,PC2)) + vd2;

F2 = m2 * vcd2 ;
N2 = IC2 * wd2 + cross(w2,IC2*w2);


% i = 3
w3 = R_32 * w2 + dt3*R_23(1:3,3);
wd3 = R_32 * wd2 + R_32 * cross(w2, dt3*R_23(1:3,3)) + ddt3*R_23(1:3,3);

vd3 = R_32 * (cross(wd2, P_23) + cross(w2, cross(w2, P_23)) + vd2);
vcd3 = cross(wd3,PC3) + cross(w3,cross(w3,PC3)) + vd3;

F3 = (m3+m4) * vcd2; 
N3 = IC3 * wd3 + cross(w3,IC3*w3);



% Outward iteration 

% i = 3
f3 = R_34 * f4 + F3;
n3 = N3 + R_34*n4 + cross(PC3, F3) + cross(P_34, R_34*f4);
f3 = simplify(f3)
n3 = simplify(n3)

% i = 2
f2 = R_23 * f3 + F2;
n2 = N2 + R_23*n3 + cross(PC2, F2) + cross(P_23, R_23*f3);
f2 = simplify(f2)
n2 = simplify(n2)

% i = 1
f1 = R_12 * f2 + F1;
n1 = N1 + R_12*n2 + cross(PC1, F1) + cross(P_12, R_12*f2);
f1 = simplify(f1)
n1 = simplify(n1)

tau1 = n1(3);
tau2 = n2(3);
tau3 = n3(3);

TAU = [tau1; tau2;tau3]

M11 = subs(tau1, [ddt1, ddt2, ddt3, dt1, dt2, dt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
M12 = subs(tau1, [ddt1, ddt2, ddt3, dt1, dt2, dt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
M13 = subs(tau1, [ddt1, ddt2, ddt3, dt1, dt2, dt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

M21 = subs(tau2, [ddt1, ddt2, ddt3, dt1, dt2, dt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
M22 = subs(tau2, [ddt1, ddt2, ddt3, dt1, dt2, dt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
M23 = subs(tau2, [ddt1, ddt2, ddt3, dt1, dt2, dt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

M31 = subs(tau3, [ddt1, ddt2, ddt3, dt1, dt2, dt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
M32 = subs(tau3, [ddt1, ddt2, ddt3, dt1, dt2, dt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
M33 = subs(tau3, [ddt1, ddt2, ddt3, dt1, dt2, dt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

V1F1 = subs(tau1, [ddt1, ddt2, ddt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
V2F2 = subs(tau2, [ddt1, ddt2, ddt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
V3F3 = subs(tau3, [ddt1, ddt2, ddt3, g, f4x, f4y, f4z, n4x, n4y, n4z], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

G1 = subs(tau1, [ddt1, ddt2, ddt3, dt1, dt2, dt3], [0, 0, 0, 0, 0, 0])
G2 = subs(tau2, [ddt1, ddt2, ddt3, dt1, dt2, dt3], [0, 0, 0, 0, 0, 0])
G3 = subs(tau3, [ddt1, ddt2, ddt3, dt1, dt2, dt3], [0, 0, 0, 0, 0, 0])

%% Functions

function [T] = transformationMatrix(DH_row)

 T = [cos(DH_row(4))                     -sin(DH_row(4))                    0                   DH_row(2);
      sin(DH_row(4))*cos(DH_row(1))      cos(DH_row(4))*cos(DH_row(1))      -sin(DH_row(1))     -sin(DH_row(1))*DH_row(3);
      sin(DH_row(4))*sin(DH_row(1))      cos(DH_row(4))*sin(DH_row(1))      cos(DH_row(1))      cos(DH_row(1))*DH_row(3);
      0                                  0                                  0                   1]; 
end