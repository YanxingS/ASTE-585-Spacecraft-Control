clc; clear; close all; 
ImportData;
%problem 1 hw 2 compute satellite attitude wrt ECI frame over 5 orbital
%period
%% Define parameters
global e a miu_earth i omega wp period t_span alpha epsilon big_delta roll yaw pitch
miu_earth = 3.986E5; % celestial constant of earth
a = 9304; % semi-major axis in km
e = 0.25; % ecentricity
i = deg2rad(28.5); % orbit inclination 
omega = deg2rad(30); % longitude of ascending node
wp = deg2rad(45); % argument of perigee 
true_anomaly = []; % angle of spacecraft from perigee
period = 8931; % second
epsilon = deg2rad(11);
big_delta = deg2rad(20);
t_span = period*5; %span over five orbital period

E = [];
r = [];
v = [];
index = 1; % index for E
x_axis = [];
%% solve for ecentric anomaly

for t = 0:60:t_span 
    f = @(E)[(a^3/miu_earth)^0.5*(E-e*sin(E))-t]; % puts everything on one side
    E(index) = fsolve(f,0);
    x_axis(index) = t; % store x_axis value 
    % r(t+1) = a*(1-e*cos(E(t+1))); 
    % v(t+1) = sqrt(miu_earth*(2/r(t+1)-1/a));
    index = index+1;
end
%% solve for true anomaly
true_anomaly = 2*atan(((1+e)/(1-e))^0.5*tan(E./2));

true_anomaly = true_anomaly*180/pi;
%plot true_anomoly over one period

figure(1);
plot(x_axis,true_anomaly);
title('true anomoly vs one period')
xlabel('time (seconds)')
ylabel('true anomoly in degrees')

true_anomaly = true_anomaly*pi/180; % convert true_anomaly back to radian

%% this section computes attitude of spacecraft wrt ECI frame

roll = deg2rad(20);
pitch = deg2rad(15);
yaw = deg2rad(30);

for j = 1:length(true_anomaly)

% generate full rotation matrix going from I to B and then taking it's
% euler angles
alpha = wp+true_anomaly(j);
R_I2N = [cos(omega) sin(omega) 0 ;-cos(i)*sin(omega) cos(i)*cos(omega) sin(i); sin(i)*sin(omega) -sin(i)*cos(omega) cos(i)];
R_N2O = [cos(alpha) sin(alpha) 0 ; -sin(alpha) cos(alpha) 0; 0 0 1];
R_O2R = [0 1 0;0 0 -1;-1 0 0];
R_R2B = [cos(yaw)*cos(pitch) sin(yaw)*cos(roll)+cos(yaw)*sin(pitch)* ...
    sin(roll) sin(yaw)*sin(roll)-cos(yaw)*sin(pitch)*cos(roll); ...
    -sin(yaw)*cos(pitch) cos(yaw)*cos(roll)-sin(yaw)*sin(pitch)*sin(roll) cos(yaw)*sin(roll)+sin(yaw)*sin(pitch)*cos(roll); ...
    sin(pitch) -cos(pitch)*sin(roll) cos(pitch)*cos(roll)];

Final_R = R_R2B*R_O2R*R_N2O*R_I2N;

attitude_ECI(1:3,j) = [atan2(Final_R(2,3),Final_R(3,3));asin(-Final_R(1,3));atan2(Final_R(1,2),Final_R(1,1))];

end
%% plotting attitude of satellite wrt ECI 

figure(2)
plot(x_axis,attitude_ECI(1,:),x_axis,attitude_ECI(2,:),x_axis,attitude_ECI(3,:),'LineWidth',3);
legend('roll','pitch','yaw')
xlabel('time (seconds)')
ylabel('attitude (radian)')
grid on

%% problem 2a, generate commanded attitude of satellite following ZPOP-YPSL

% extract sun movement data from NASA Horizon System, between July-9th to
% July 10th

% Note, use import to get longitudinal data

S = horizonsresults1{:,1};
S = deg2rad(S);
ecliptic_i = deg2rad(23.44);

for j = 1:length(true_anomaly)
alpha = wp+true_anomaly(j);

R_I2N = [cos(omega) sin(omega) 0 ;-cos(i)*sin(omega) cos(i)*cos(omega) sin(i); sin(i)*sin(omega) -sin(i)*cos(omega) cos(i)];
R_N2O = [cos(alpha) sin(alpha) 0 ; -sin(alpha) cos(alpha) 0; 0 0 1];
R_O2R = [0 1 0;0 0 -1;-1 0 0];
R_R2B = [cos(yaw)*cos(pitch) sin(yaw)*cos(roll)+cos(yaw)*sin(pitch)* ...
    sin(roll) sin(yaw)*sin(roll)-cos(yaw)*sin(pitch)*cos(roll); ...
    -sin(yaw)*cos(pitch) cos(yaw)*cos(roll)-sin(yaw)*sin(pitch)*sin(roll) cos(yaw)*sin(roll)+sin(yaw)*sin(pitch)*cos(roll); ...
    sin(pitch) -cos(pitch)*sin(roll) cos(pitch)*cos(roll)];
R_E2I = transpose([cos(S(j)) sin(S(j)) 0; -sin(S(j)) cos(S(j)) 0 ; 0 0 1]*[1 0 0; 0 cos(ecliptic_i) sin(ecliptic_i); 0 -sin(ecliptic_i) cos(ecliptic_i)]);

R_E2R = R_O2R*R_N2O*R_I2N*R_E2I;

commanded_attitude(:,j) = [pi/2;atan2(-R_E2R(3,1),R_E2R(1,1));0];

end
commanded_attitude = commanded_attitude*180/pi;

figure(3)
plot(x_axis,commanded_attitude(1,:),x_axis,commanded_attitude(2,:),x_axis,commanded_attitude(3,:),'LineWidth',3);
legend('roll','pitch','yaw')
xlabel('time (seconds)')
ylabel('commanded attitude (degrees)')
grid on

%% problem 2b, assume that we now have a solar arrary and can rotate about satellite's y body axis , plot and compute the gimbal angle of solar arrary
% the solar arrary's z axis must be pointing in sun line for optimal
% effeciency

for j = 1:length(true_anomaly)

alpha = wp+true_anomaly(j);

R_I2N = [cos(omega) sin(omega) 0 ;-cos(i)*sin(omega) cos(i)*cos(omega) sin(i); sin(i)*sin(omega) -sin(i)*cos(omega) cos(i)];
R_N2O = [cos(alpha) sin(alpha) 0 ; -sin(alpha) cos(alpha) 0; 0 0 1];
R_O2R = [0 1 0;0 0 -1;-1 0 0];
R_R2B = [cos(yaw)*cos(pitch) sin(yaw)*cos(roll)+cos(yaw)*sin(pitch)* ...
    sin(roll) sin(yaw)*sin(roll)-cos(yaw)*sin(pitch)*cos(roll); ...
    -sin(yaw)*cos(pitch) cos(yaw)*cos(roll)-sin(yaw)*sin(pitch)*sin(roll) cos(yaw)*sin(roll)+sin(yaw)*sin(pitch)*cos(roll); ...
    sin(pitch) -cos(pitch)*sin(roll) cos(pitch)*cos(roll)];
R_E2I = transpose([cos(S(j)) sin(S(j)) 0; -sin(S(j)) cos(S(j)) 0 ; 0 0 1]*[1 0 0; 0 cos(ecliptic_i) sin(ecliptic_i); 0 -sin(ecliptic_i) cos(ecliptic_i)]);

R_B2E = transpose(R_R2B*R_O2R*R_N2O*R_I2N*R_E2I);

gimbal_angle(j) = atan2(R_B2E(1,1),R_B2E(1,3));

end

gimbal_angle = gimbal_angle*180/pi;

figure(4)
plot(x_axis,gimbal_angle);
xlabel('time (seconds)')
ylabel('delta gimbal angle (degrees)')
grid on

%% problem 3 , find Earth's magnetic field wrt body frame

%first we need to find the time which the prime meridian passes the vernal
%equinox, then add 5 hours to that.

alpha_e_test = hour_angle(9763500); % base on 2024 march 20, and to 5 hours after first intersection near July 9th 4pm

%after some trial and error I found that first intersection of
%prime meridian and vernal equinox on July 9th is approximately 9763500
%seconds after march 20th, adding five hours to it would be 9781500, after
%inserting this time to time calculator we found that our start time would
%be July 11th 8:11 AM, which is 144660 seconds after July 9th 4 pm

start_time = 144660; % wrt July 9th 4 pm

magnetic_vector = [0 0 1]';

%here we need to solve for true anomaly one more time because our start
%time is different

x_axis_p3 = [];
index_p3 = 1;
E_p3 = [];
alpha_e = [];
for t = start_time:60:start_time+t_span 
    f_p3 = @(E_p3)[(a^3/miu_earth)^0.5*(E_p3-e*sin(E_p3))-t]; % puts everything on one side
    E_p3(index_p3) = fsolve(f_p3,0);
    x_axis_p3(index_p3) = t; % store x_axis value 
    alpha_e(index_p3) = hour_angle(9723240+t);
    % r(t+1) = a*(1-e*cos(E(t+1))); 
    % v(t+1) = sqrt(miu_earth*(2/r(t+1)-1/a));
    index_p3 = index_p3+1;
end

true_anomaly_p3 = 2*atan(((1+e)/(1-e))^0.5*tan(E_p3./2));

magnetic_vector_B = [];

for j = 1:length(true_anomaly_p3)

alpha = wp+true_anomaly_p3(j);

R_I2N = [cos(omega) sin(omega) 0 ;-cos(i)*sin(omega) cos(i)*cos(omega) sin(i); sin(i)*sin(omega) -sin(i)*cos(omega) cos(i)];
R_N2O = [cos(alpha) sin(alpha) 0 ; -sin(alpha) cos(alpha) 0; 0 0 1];
R_O2R = [0 1 0;0 0 -1;-1 0 0];
R_R2B = [cos(yaw)*cos(pitch) sin(yaw)*cos(roll)+cos(yaw)*sin(pitch)* ...
    sin(roll) sin(yaw)*sin(roll)-cos(yaw)*sin(pitch)*cos(roll); ...
    -sin(yaw)*cos(pitch) cos(yaw)*cos(roll)-sin(yaw)*sin(pitch)*sin(roll) cos(yaw)*sin(roll)+sin(yaw)*sin(pitch)*cos(roll); ...
    sin(pitch) -cos(pitch)*sin(roll) cos(pitch)*cos(roll)];
% R_E2I = transpose([cos(S(j)) sin(S(j)) 0; -sin(S(j)) cos(S(j)) 0 ; 0 0 1]*[1 0 0; 0 cos(ecliptic_i) sin(ecliptic_i); 0 -sin(ecliptic_i) cos(ecliptic_i)]);

R_B2E = transpose(R_R2B*R_O2R*R_N2O*R_I2N*R_E2I);

R_M2G = transpose([1 0 0; 0 cos(epsilon) sin(epsilon) ; 0 -sin(epsilon) cos(epsilon)]*[cos(big_delta) sin(big_delta) 0; -sin(big_delta) cos(big_delta) 0; 0 0 1]);

R_G2I = transpose([cos(alpha_e(j)) sin(alpha_e(j)) 0; -sin(alpha_e(j)) cos(alpha_e(j)) 0; 0 0 1]);

R_M2B = R_R2B*R_O2R*R_N2O*R_I2N*R_G2I*R_M2G;

magnetic_vector_B(:,j) = R_M2B*magnetic_vector;

end


figure(5)
plot(x_axis_p3,magnetic_vector_B(1,:),x_axis_p3,magnetic_vector_B(2,:),x_axis_p3,magnetic_vector_B(3,:),'LineWidth',3);
legend('x_B','y_B','z_B')
xlabel('time (seconds)')
ylabel('magnetic field in B frame')
title('Magnetic field vector in B frame')
grid on


function angle = hour_angle(time)
    rate = 7.272e-5; % rad/s
    %time to July 9th 4pm is 9723240
    current_angle = rate*time*180/pi;
    %get remainder
    angle = rem(current_angle,360);
end