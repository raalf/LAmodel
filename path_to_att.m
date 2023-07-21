function [roll, hdg, time_datenum, trk, ground_speed, mid_p, avec, gvec, wvec, dt, fp_angle,mid_zp] = path_to_att(direc, x_p, y_p, air_speed, wind_speed, wind_from, z_p)

%% sectional midpoints
% mid_p - midpoint of each section
% slen - section length
% dt - time on this section

% gvec - ground vector
% wvec - wind vector
% avec - airspeed vector

% ground_speed - ground speed magnitude
% wind_speed - winds speed magnitude
% air_speed - airspeed magnitude

min_gspd = 3.0; % display warning if gspd drops below min_gspd

% remove repeated points
idx = find(diff(x_p)==0 & diff(y_p)==0 & diff(z_p)==0);
x_p(idx) = [];
y_p(idx) = [];
z_p(idx) = [];

mid_p = [(x_p(1:end-1) + x_p(2:end))./2 (y_p(1:end-1) + y_p(2:end))./2];
mid_zp = (z_p(1:end-1) + z_p(2:end))./2;

% Tangential and normal vectors for each sectional midpoint (tangent to ground path)
gvec = [(x_p(2:end) - x_p(1:end-1)), (y_p(2:end) - y_p(1:end-1))];
slen = sqrt(sum(gvec.^2,2)); % section length;
gvec = gvec./slen;

gvec_z = [(x_p(2:end) - x_p(1:end-1)), (y_p(2:end) - y_p(1:end-1)), (z_p(2:end) - z_p(1:end-1))];
slen_z = sqrt(sum(gvec_z.^2,2));
fp_angle = asind((z_p(2:end)-z_p(1:end-1))./slen_z); % flight path angle

air_speed = air_speed.*(cosd(fp_angle)); % or assume cosd(fp_angle)=1 since angle is small

wvec = repmat([sind(wind_from-180) cosd(wind_from-180)], size(mid_p,1), 1);

theta = acosd(dot(gvec, wvec, 2)); % angle between the wind vector and tangential vector of ground vector (difference between the wind direction and our ground track)
% wta=acosd(dot(gvec, wvec, 2)./(sqrt(sum(gvec.^2,2)).*sqrt(sum(wvec.^2,2))));
alpha = asind((wind_speed.*sind(theta))./air_speed); % drift angle: takes into account the way the wind will change our direction of travel (difference between the heading and our ground track)
beta = 180 - theta - alpha; % angle to correct the effects of wind on our heading, (the angle between wind and heading)

ground_speed = (air_speed.*sind(beta))./sind(theta); % sine rule
% ground_speed = (sqrt((air_speed^2)-((wind_speed.^2).*((sind(theta)).^2))))-(wind_speed.*(cosd(theta)));
% V_g = (sqrt((air_speed^2)-((wind_speed.^2).*((sind(theta)).^2))))-(wind_speed.*(cosd(theta)));

% corrections
ground_speed(alpha == 180) = wind_speed - air_speed(alpha == 180);
ground_speed(theta == 180) = air_speed(theta == 180) - wind_speed;
ground_speed(beta == 180) = air_speed(beta == 180) + wind_speed;

dt = slen./ground_speed;

avec = (gvec.*ground_speed - wvec.*wind_speed)./air_speed;
avec = avec./sqrt(sum(avec.^2,2)); % avec should already be a unit vector, but its sometimes a little off. Problem???

% Calculating roll angle at each segment
[~,R,~] = curvature([mid_p(end,:); mid_p; mid_p(1,:)]);
R = R(2:end-1);

% Bank angle is solely a function of cruise speed and turn radius
roll = atan2d((ground_speed.^2),(R.*9.81)); % atan2d(y,x)
if strcmpi(direc, 'CCW')
    roll = roll.*-1;
end

roll(isnan(roll)) = 0;

% Track heading, assuming distance between two points is small
trk = rem(rad2deg(cart2pol(gvec(:,2), gvec(:,1)))+360, 360);
hdg = rem(rad2deg(cart2pol(avec(:,2), avec(:,1)))+360, 360);

time_datenum = cumsum(dt)/3600/24;

if min(ground_speed) < min_gspd
    warning("Ground speed drops below %.1f m/s.\n", min_gspd);
end

% 
% 
% hFig1 = figure(1);
% clf(1);
% subplot(1,2,1);
% plot(x_p, y_p, '-k');
% hold on
% quiver(mid_p(:,1), mid_p(:,2), air_speed.*avec(:,1), air_speed.*avec(:,2), 'm')
% 
% % quiver(mid_p(:,1), mid_p(:,2), wind_vec(:,1), wind_vec(:,2), 1, 'b')
% quiver(mid_p(:,1), mid_p(:,2), gvec(:,1), gvec(:,2), 'r')
% hold off
% axis equal
% grid minor
% box on
% axis tight
% 
% subplot(1,2,2);
% plot(cumsum(dt), roll, '-ok');
% grid minor
% box on
% axis tight
% 
% yyaxis right
% tmp = acosd(dot(avec, gvec, 2));
% plot(cumsum(dt), tmp, '--bs');

% 
% 
% figure(2)
% plot(cumsum(dt), trk)
% hold on
% plot(cumsum(dt), hdg)
% plot(cumsum(dt), ground_speed)
% hold off
% grid minor
end

function [L,R,k] = curvature(X)
% Radius of curvature and curvature vector for 2D or 3D curve
%  [L,R,Kappa] = curvature(X)
%   X:   2 or 3 column array of x, y (and possibly z) coordiates
%   L:   Cumulative arc length
%   R:   Radius of curvature
%   k:   Curvature vector
  N = size(X,1);
  dims = size(X,2);
  if dims == 2
    X = [X,zeros(N,1)];  % Do all calculations in 3D
  end
  L = zeros(N,1);
  R = NaN(N,1);
  k = NaN(N,3);
  for i = 2:N-1
    [R(i),~,k(i,:)] = circumcenter(X(i,:)',X(i-1,:)',X(i+1,:)');
    L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
  end
  i = N;
  L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
  if dims == 2
    k = k(:,1:2);
  end
end

function [R,M,k] = circumcenter(A,B,C)
% Center and radius of the circumscribed circle for the triangle ABC
%  A,B,C  3D coordinate vectors for the triangle corners
%  R      Radius
%  M      3D coordinate vector for the center
%  k      Vector of length 1/R in the direction from A towards M
%         (Curvature vector)
  D = cross(B-A,C-A);
  b = norm(A-C);
  c = norm(A-B);
  if nargout == 1
    a = norm(B-C);     % slightly faster if only R is required
    R = a*b*c/2/norm(D);
    return
  end
  E = cross(D,B-A);
  F = cross(D,C-A);  
  G = (b^2*E-c^2*F)/norm(D)^2/2;
  M = A + G;
  R = norm(G);  % Radius of curvature
  if R == 0
    k = G;
  else
    k = G'/R^2;   % Curvature vector
  end
end







