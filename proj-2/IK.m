% Inverse Kinematics for MAE 263B Project 2

function [theta1, theta2, d3, theta4] = IK(T)

[x, y, z, theta] = mat_to_pos(T);

% Fixed lengths [m]
a1 = 0.325;
a2 = 0.225;
d1 = 0.416;
d4 = 0.093;

% Joint limits
t1_min = -170;
t1_max = 170;       %[deg]
t2_min = -145;
t2_max = 145;       %[deg]
d3_min = 0;
d3_max = 0.150;     %[m]
t4_min = -360;
t4_max = 360;       %[deg]

% z is independant with d3
d3_val = .323-z;

% Solve for theta 2 using x,y, law of cosines
c_theta_2 = (x^2+y^2-a1^2-a2^2)/(2*a1*a2);
s_theta_2_pos = sqrt(1-c_theta_2^2);
s_theta_2_neg = -sqrt(1-c_theta_2^2);

theta_2_val1 = atan2(s_theta_2_pos,c_theta_2);
theta_2_val2 = atan2(s_theta_2_neg,c_theta_2);

% Solve for theta 1 using the leftover angle
L3_1 = a1+cos(theta_2_val1)*a2;
L3_2 = a1+cos(theta_2_val2)*a2;
L4_1 = sin(theta_2_val1)*a2;
L4_2 = sin(theta_2_val2)*a2;

theta_1_val1 = atan2(y,x) - atan2(L4_1,L3_1);
theta_1_val2 = atan2(y,x) - atan2(L4_2,L3_2);

% Find theta 4 by making it match the desired orientation
theta_4_val1 = theta - pi/4 - theta_1_val1 - theta_2_val1;
theta_4_val2 = theta - pi/4 - theta_1_val2 - theta_2_val2;

% Check if any of the values are outside the joint ranges
if (d3_val < d3_min || d3_val > d3_max)
    d3 = -1;
else
    d3 = d3_val;
end

% Prefer the value 1, if it doesn't work check for val2
if (check_angles(theta_1_val1, theta_2_val1, theta_4_val1) == 1)
    theta1 = theta_1_val1;
    theta2 = theta_2_val1;
    theta4 = theta_4_val1;
    return
elseif (check_angles(theta_1_val2, theta_2_val2, theta_4_val2) == 1)
    theta1 = theta_1_val2;
    theta2 = theta_2_val2;
    theta4 = theta_4_val2;
    return
else
    theta1 = -1;
    theta2 = -1;
    theta4 = -1;
    return
end
end



function check = check_angles(t1_test, t2_test, t4_test)
    % Joint limits
    t1_min = -170;
    t1_max = 170;       %[deg]
    t2_min = -145;
    t2_max = 145;       %[deg]
    d3_min = 0;
    d3_max = 0.150;     %[m]
    t4_min = -360;
    t4_max = 360;       %[deg]

    check = -1;
    if (t1_test < t1_min || t1_test > t1_max)
        return
    elseif (t2_test < t2_min || t2_test > t2_max)
        return
    elseif (t4_test < t4_min || t4_test > t4_max)
        return
    else
        check = 1;
        return
    end
end

function [x, y, z, theta] = mat_to_pos(T)

x = T(1,4);
y = T(2,4);
z = T(3,4);
eul_angles = rotm2eul(T(1:3,1:3));
theta = eul_angles(1);

end