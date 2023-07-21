clc
clear

%% Parameters
air_speed = 11; % m/s

wind_speed = 8; % m/s
wind_from = 90; % degrees, wind from

startLocalDate = datenum(2020,11,12,12,0,0); % Local Time
timezone = -5;
Latitude = 43.838; %
Longitude = -80.442; %

%% Variables
loiter_radius = [50:500]';
direc = {'CW', 'CCW'};

%% objective function
f_obj = @(x,y)loiter_score(x, y, 36, air_speed, ...
    wind_speed, wind_from, startLocalDate, timezone, Latitude, Longitude);

for i = 1:length(loiter_radius)
    for j = 1:2
        sc(i,j) = f_obj(direc{j}, loiter_radius(i));
    end
end

[ia,ib] = min(sc,[],'all','linear');
[idx_r, idx_d] = ind2sub(size(sc), ib);

disp([direc{idx_d}, ' circuit at a loiter radius of ' num2str(loiter_radius(idx_r)), ' meters -> Score: ', num2str(ia)]);




