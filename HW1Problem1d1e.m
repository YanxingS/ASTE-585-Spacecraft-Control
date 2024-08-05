clc; clear; close all;
%% Comments

% this program solves for the true anomoly of the satellite over time using
% two equations that is given from lecture 2 page 21. We use fsolve to
% compute E
%% setting up parameters

global e a miu_earth

e = 0.25; % ecentricity of orbit
a = 9294; % semi-major axis in km
miu_earth = 3.986E5; % celestial parameter of earth
period = 8917; % one period for orbit

E = [];
true_anomoly = [];
r = []; % this is for next question
v = []; % this is for next question
%% first compute ecentric anomoly
% refer to lecture 2 ASTE 585 page 21

for t = 0:1:period
    f = @(E)[(a^3/miu_earth)^0.5*(E-e*sin(E))-t]; % puts everything on one side
    E(t+1) = fsolve(f,0);

    r(t+1) = a*(1-e*cos(E(t+1))); 
    v(t+1) = sqrt(miu_earth*(2/r(t+1)-1/a));
end

%% solve for true anomoly and plot

true_anomoly = 2*atan(((1+e)/(1-e))^0.5*tan(E./2));
%covert to angels
true_anomoly = true_anomoly.*180/pi;
%plot true_anomoly over one period
t_span = 0:1:period;

figure(1);
plot(t_span,true_anomoly);
title('true anomoly vs one period')
xlabel('time (seconds)')
ylabel('true anomoly in degrees')

clc;

%
%
%% this is part e, compute the velocity of the satelite as a function of time
% we know E from previous calculation and thus we also got r 
% and from r we can get v from vis-viva equation, now we can also plot it
% here

figure(2)
plot(t_span,v);
title('velocity of satellite over one period')
xlabel('time (seconds)')
ylabel('velocity in km/s')
drawnow