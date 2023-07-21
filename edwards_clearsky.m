function [I, I_direct, I_diffuse]= edwards_clearsky(Time, Elevation, altitude)

% INPUTS
%   altitude = flight altitude in m
%   Elevation = sun elevation (from solar models)
%
%   Time is a struct with the following elements, note that all elements
%     can be column vectors, but they must all be the same length. Time is
%     entered as a struct which must include a value for the offset from
%     UTC. Time is absolutely specified by the date, time, and the number
%     of hours of offset from that date and time to UTC. For example, if
%     you live in Boston, USA, and your data is timestamped in local standard
%     time, your UTC offset should be -5.
%
%   Time.year = The year in the gregorian calendar
%   Time.month = the month of the year (January = 1 to December = 12)
%   Time.day = the day of the month
%   Time.hour = the hour of the day
%   Time.minute = the minute of the hour
%   Time.second = the second of the minute
%   Time.UTCOffset = the UTC offset code, using the convention
%     that a positive UTC offset is for time zones east of the prime meridian
%     (e.g. EST = -5)


% References
% [1] Edwards, 2016, Maximizing Net Power in Circular Turns for Solar and 
%       Autonomous Soaring Aircraft, section III(c)
% [2] Hottel, 1976, A simple model for stimating the transmittance of
%       direct solar radation through claer atmospheres

h = altitude/1000; % convert to km
Elevation(Elevation<0.0) = 0; % model doesn't do well at low sun elevation angles
z = 90 - Elevation;

% Require the field Location.altitude
p = inputParser;
p.addRequired('Time',@isstruct);
p.addRequired('Location',@(x) all(isstruct(x) & isfield(x,'altitude')));

% Determine day of year and extraterrestrial normal radiation for each time
% instant in Time
DayOfYear = pvl_date2doy(Time.year, Time.month, Time.day);
n = DayOfYear;

% Solar constant
I_sc = 1367; % W/m2

% Extraterrestrial direct normal insolation
I_on = I_sc.*(1+(0.033.*(cosd(360.*n/365))));

% Correction factors for climate type and season [2]
% Assumed mid-latitude summer constants, 23 km visibility
r0 = 0.97;
r1 = 0.99;
rk = 1.02;

% Constants for standard atmmosphere for altitudes less than 2.5 km and
% visibility greater than 23 km
a0 = 0.4237 - 0.00821.*((6.0-h).^2);
a1 = 0.5055 + 0.00595.*((6.5-h).^2);
k = 0.2711 + 0.01858.*((2.5-h).^2);

% Direct atmospheric transmittance
tau_direct = r0.*a0 + r1.*a1.*exp((-rk.*k)./(cosd(z)));

% Clear-sky direct normal radiation
I_direct = I_on.*tau_direct;

% Diffuse transmittance
tau_diffuse = 0.271  -0.294.*tau_direct;

% Diffuse insolation
I_diffuse = I_on.*tau_diffuse;

% Total normal solar insolation
I = I_direct + I_diffuse;

end
