function [Azimuth, Zenith, Heading, Elevation, uVec] = solarPosition(GMTdatenum, Latitude, Longitude)

GMTdatenum = reshape(GMTdatenum,[],1); % ensures the input time vector is nx1
theta = deg2rad(Longitude); % convert Longitude input to radians
psi = deg2rad(Latitude); % convert Latitude input to radians


% Grena 2012. Sec.3.1 - Starting point: time scale computation
% t = datenum(y,m,d,h,0,0)-datenum(2060,0,0)-1; % identical to Grena 2012. p.1329 Eq.2
t = GMTdatenum - datenum(2060,0,0) - 1; % identical to Grena 2012. p.1329 Eq.2

dtau = 96.4 + 0.00158.*t; %Grena 2012. p.1325 Eq.1
te = t + 1.1574e-5 .* dtau; %Grena 2012. p.1329 Eq.3

% Grena 2012. Sec.3.2 - Algorithm 1
omega = 0.017202786; %[day^-1]
alpha = -1.38880 + 1.72027920e-2.*te ...
    + 3.199e-2.*sin(omega.*te) - 2.65e-3.*cos(omega.*te) ...
    + 4.050e-2.*sin(2.*omega.*te) + 1.525e-2.*cos(2.*omega.*te);  %Grena 2012. p.1329 Eq.4 [0,2pi]
%alpha = mod(alpha, 2.*pi); %not needed %%%% whyy?
delta = 6.57e-3 ...
    + 7.347e-2.*sin(omega.*te) - 3.9919e-1.*cos(omega.*te) ...
    + 7.3e-4.*sin(2.*omega.*te) - 6.6e-3.*cos(2.*omega.*te);  %Grena 2012. p.1329 Eq.5 [-pi/2,pi/2]
H = 1.75283 + 6.3003881.*t + theta - alpha;
H = mod(H + pi, 2.*pi) - pi; % [-pi,pi]

% Grena 2012. Sec.3.7 - Final steps
e0 = asin(sin(psi).*sin(delta) + cos(psi).*cos(delta).*cos(H)); %Grena 2012. p.1332 Eq.21
dpe = -4.26e-5.*cos(e0); %Grena 2012. p.1332 Eq.22
ep = e0 + dpe; % elevation of the sun
Gamma = atan2(sin(H),(cos(H).*sin(psi))-(tan(delta).*cos(psi))); %Grena 2012. p.1332 Eq.23 [-pi,pi]
Z = pi/2 - ep; %Grena 2012. p.1332 Eq.25

% constant right now
P = 0.96; % atm @ 350m
T = 15; % deg C

dre = ((0.08422*P)/(273+T)).*(1./(tan(ep+(0.003138./(ep+0.08919))))); %Grena 2012. p.1332 Eq.24 refraction correction
Z = Z-dre; %Grena 2012. p.1332 Eq.25 [0,pi]

Zenith = rad2deg(Z);
Elevation = rad2deg(ep); % Elevation Angle in degrees
Azimuth = rad2deg(Gamma); % Azimuthal angle on the horizontal plane, 0 towards south, +ve towards west, and -ve towards east. [-180,180]
Heading = rem(Azimuth+180,360); % Calculate heading bearing by adding 180 degrees to Azimuthal angle [0,360]

[x,y,z] = sph2cart(-pi/2-Gamma, ep, 1); % Convert elevation and azimuth angles to cartesian coordinates with radius of 1
uVec = [-x,-y,-z]; % Create unit vector from position of the sun to origin

