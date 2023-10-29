% Forward Kinematics for MAE 263B Project 2

function T = FK(theta1, theta2, d3, theta4)

% Fixed lengths [m]
a1 = 0.325;
a2 = 0.225;
d1 = 0.416;
d4 = 0.093;

% Forward Kinematics matrices
T01 = mat_from_DH(0,0,d1,theta1);
T12 = mat_from_DH(0,a1,0,90+theta2);
T23 = mat_from_DH(0,a2,-d3,0);
T34 = mat_from_DH(0,0,-d4,theta4);

T = T01*T12*T23*T34;

end

function matrix = mat_from_DH(alpha_im1, a_im1, d_i, theta_i)
    %Returns a matrix given the DH parameter inputs
    %angles are in degrees

    matrix = [cosd(theta_i) -sind(theta_i) 0 a_im1;
              sind(theta_i)*cosd(alpha_im1) cosd(theta_i)*cosd(alpha_im1) -sind(alpha_im1) -sind(alpha_im1)*d_i;
              sind(theta_i)*sind(alpha_im1) cosd(theta_i)*sind(alpha_im1) cosd(alpha_im1) cosd(alpha_im1)*d_i;
              0 0 0 1];
end