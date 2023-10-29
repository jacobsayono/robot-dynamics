% Main function for MAE 263B Project 2
close all; clear all; clc;

% Fixed lengths [m]
a1 = 0.325;
a2 = 0.225;
d1 = 0.416;
d4 = 0.093;

% Define key points:
default = [0, -1, 0, .325;
           1, 0, 0, .225;
           0, 0, 1, .203;
           0, 0, 0, 1];
feeder =  [0, -1, 0, .325;
           1, 0, 0, .225;
           0, 0, 1, .180;
           0, 0, 0, 1];
goal1 =   [1, 0, 0, .280;
           0, 1, 0, .240;
           0, 0, 1, .180;
           0, 0, 0, 1];
goal2 =   [0, -1, 0, .280;
           1, 0, 0, .330;
           0, 0, 1, .180;
           0, 0, 0, 1];
goal3 =   [-1, 0, 0, .370;
           -1, 0, 0, .240;
           0, 0, 1, .180;
           0, 0, 0, 1];
goal4 =   [0, 1, 0, .370;
           -1, 0, 0, .240;
           0, 0, 1, .180;
           0, 0, 0, 1];

% Define via points
a = [0, -1, 0, .325;
     1, 0, 0, .225;
     0, 0, 1, .200;
     0, 0, 0, 1];
b = [1, 0, 0, .280;
     0, 1, 0, .240;
     0, 0, 1, .200;
     0, 0, 0, 1];
c = [0, -1, 0, .280;
     1, 0, 0, .330;
     0, 0, 1, .200;
     0, 0, 0, 1];
d = [-1, 0, 0, .370;
     0, -1, 0, .330;
     0, 0, 1, .200;
     0, 0, 0, 1];
e = [0, 1, 0, .370;
     -1, 0, 0, .240;
     0, 0, 1, .200;
     0, 0, 0, 1];

% Example of both Forward and Inverse Kinematics functions
T = FK(0,0,.120,0)
[t1, t2, d3, t4] = IK(a)


% Build the joints of the arm
deg = pi/180;
L1 = Revolute('alpha',0,'a',0,'d',d1,'qlim',[-170 170]*deg);
L2 = Revolute('alpha',0,'a',a1,'d',0,'qlim',[-145 145]*deg);
L3 = Prismatic('alpha',pi,'a',a2,'theta',0,'qlim',[0 .15]);
L4 = Revolute('alpha',0,'a',0,'d',d4,'qlim',[-360 360]*deg);

% Tool length already incorporated into parameters
tool = transl(0,0,0);

% Build robot
PCB_Builder = SerialLink([L1 L2 L3 L4], 'name', 'PCB_Builder', 'tool', tool);


% For joint space trajectories, convert waypoints into joint space values,
% then use interpolation to get between points




