function [roll, hdg, time_datenum, trk, ground_speed, mid_p, avec, gvec, wvec, dt, Preq, aoa,...
    V_inf, CL_tf, fp_angle, mid_zp] = path_to_att_mc4(direc, x_p, y_p, wind_speed, ...
    wind_from, C_L, mass,powerdrawdata, eta, z_p)

% if nargin==5
CreateVData_FLT = powerdrawdata.data;
% else
% load('simpleModel.mat')
% end
velocities = CreateVData_FLT(:,1);
CD_data = CreateVData_FLT(:,2);
CL_data = CreateVData_FLT(:,3);
AoA_data = CreateVData_FLT(:,4);

w = mass * 9.81;
b = 6.28;
c = 0.395;
S = b * c;

rho = 1.1;

V_w = wind_speed;
psi_w = wind_from;

min_gspd = 3.0; % display warning if gspd drops below min_gspd

% remove repeated points
idx = find(diff(x_p)==0 & diff(y_p)==0);
x_p(idx) = [];
y_p(idx) = [];

% Midpoint of section
mid_p = [(x_p(1:end-1) + x_p(2:end))./2 (y_p(1:end-1) + y_p(2:end))./2];
mid_zp = (z_p(1:end-1) + z_p(2:end))./2;

% 2. Find phi (roll)
[~,R,~] = curvature([mid_p(end,:); mid_p; mid_p(1,:)]);
R = R(2:end-1);

% Tangential ground path vector for each section (unit vector)
gvec = [(x_p(2:end) - x_p(1:end-1)), (y_p(2:end) - y_p(1:end-1))];
slen = sqrt(sum(gvec.^2,2)); % section length;
gvec = gvec./slen;

gvec_z = [(x_p(2:end) - x_p(1:end-1)), (y_p(2:end) - y_p(1:end-1)), (z_p(2:end) - z_p(1:end-1))];
slen_z = sqrt(sum(gvec_z.^2,2));
fp_angle = asind((z_p(2:end)-z_p(1:end-1))./slen_z); % flight path angle

% Wind vector (unit vector)
wvec = repmat([-sind(psi_w), -cosd(psi_w)], size(mid_p,1), 1);

% Wind-track angle (angle between gvec and wvec)
phi_w = 180 - acosd(dot(gvec, wvec, 2)); % is this correct now?

% syms V_inf w C_L rho S V_w R g phi_w V_g
% eq0 = (V_inf^2) - ((V_w^2)*((sind(phi_w))^2)) - (2*(V_w*(cosd(phi_w)))*(sqrt((V_inf^2) - (V_w^2)*(sind(phi_w))^2))) + (V_w^2)*((cosd(phi_w))^2) == V_g^2
% eq = sqrt((2*w)/(C_L*rho*S*(cosd(atan2d((V_inf^2) - ((V_w^2)*((sind(phi_w))^2)) - (2*(V_w*(cosd(phi_w)))*(sqrt((V_inf^2) - (V_w^2)*(sind(phi_w))^2))) + (V_w^2)*((cosd(phi_w))^2),R*g))))) == V_inf
% eq = sqrt((2*w)./(C_L*rho*S*(cosd(atan2d(((V_inf^2) - ((V_w^2)*((sind(phi_w))^2)) - (2*(V_w*(cosd(phi_w)))*(sqrt((V_inf^2) - (V_w^2)*(sind(phi_w))^2))) + (V_w^2)*((cosd(phi_w))^2)),(R*g)))))) == V_inf;
% clear phi_old V_inf_old V_g_old
% syms phi_old V_inf_old V_g_old phi_w_old R_old
% assume(phi_old>-180 & phi_old<180 & V_g_old>0 & V_inf_old>0) % try putting in ranges for variables and using above method outside loop
%
% % 1. Find V_g
% % Groundspeed (from Mechanics of Flight, Phillips, pg. 295)
% eq1 = (((V_inf_old.^2)-((V_w.^2).*((sin(phi_w_old)).^2)))^0.5)-(V_w.*(cos(phi_w_old))) == V_g_old;
%
% eq2 = atan((V_g_old.^2)/ (R_old.*9.81)) == phi_old;
%
% eq3 = ((2.*w)./(C_L.*rho.*S.*(cos(phi_old))))^0.5 == V_inf_old;
%
% sol = solve([eq1,eq2,eq3],[phi_old,V_inf_old,V_g_old])
% sol.phi_old = subs(sol.phi_old,R_old,R)
% sol.phi_old = subs(sol.phi_old,phi_w_old,phi_w)
% sol.V_inf_old = subs(sol.V_inf_old,R_old,R)
% sol.V_inf_old = subs(sol.V_inf_old,phi_w_old,phi_w)
% sol.V_g_old = subs(sol.V_inf_old,R_old,R)
% sol.V_g_old = subs(sol.V_inf_old,phi_w_old,phi_w)
% double(sol.phi_old)

syms phi_old V_inf_old V_g_old

for i = 1:length(slen)
    if isnan(R(i)) == 1
        phi_old1 = 0;
        eq1 = (sqrt((V_inf_old.^2)-((V_w.^2).*((sind(phi_w(i))).^2))))-(V_w.*(cosd(phi_w(i)))) == V_g_old;
        eq3 = (sqrt((2.*w)./(C_L.*rho.*S))) == V_inf_old;
        sol = vpasolve([eq1,eq3],[V_inf_old,V_g_old]);
        sol.phi_old = phi_old1;
    elseif isnan(R(i)) == 0
        % 1. Find V_g
        % Groundspeed (from Mechanics of Flight, Phillips, pg. 295)
        eq1 = (sqrt((V_inf_old.^2)-((V_w.^2).*((sind(phi_w(i))).^2))))-(V_w.*(cosd(phi_w(i)))) == V_g_old;

        % 2. Find phi (roll)
        eq2 = atan2d((V_g_old.^2), (R(i).*9.81)) == phi_old;

        % 3. Find V_inf (airspeed)
        eq3 = (sqrt((2.*w)./(C_L.*rho.*S.*(cosd(phi_old))))) == V_inf_old;

        sol = vpasolve([eq1,eq2,eq3],[phi_old,V_inf_old,V_g_old]);
    end

    phi(i,1) = sol.phi_old;
    V_inf(i,1) = sol.V_inf_old;
    V_g(i,1) = sol.V_g_old;

end

phi = double(phi(:,1));
V_inf = double(V_inf(:,1)).*(cosd(fp_angle()));
V_g = double(V_g(:,1));

% Time spent in section
dt = slen./V_g;

% Airspeed vector (V_g = V_inf + V_w)
avec = ((V_g.*gvec) - (V_w.*wvec))./V_inf;
avec = real(avec);

% if CCW, roll angle is opposite
roll = phi;
if strcmpi(direc, 'CCW')
    roll = roll.*-1;
end

roll(isnan(roll)) = 0;
roll = real(roll);

% following steps are in powerdraw_mc4.m
% 3. Find C_L
% C_L = 2*W/V_inf*rho*S*cosd(phi);
% % 4. Find C_D (powerdraw data)
% C_D = interp1()
%
% % 5. Find D
% D = (W/cosd(phi))*(C_D/C_L);
%
% % 6. Find P = V_inf.*D
% P = V_inf*D;

time_datenum = cumsum(dt)/3600/24;

if min(V_g) < min_gspd
    warning("Ground speed drops below %.1f m/s.\n", min_gspd);
end

ground_speed = real(V_g);

% Heading (from airspeed vector) (angle of airspeed vector, direction aircraft is pointing in the air out of 360 deg)
hdg = rem(rad2deg(cart2pol(avec(:,2), avec(:,1)))+360, 360);

% Track (from groundspeed vector) (angle of groundspeed vector, direction aircraft is pointing on the ground out of 360 deg)
trk = rem(rad2deg(cart2pol(gvec(:,2), gvec(:,1)))+360, 360);

%% powerdraw_mc4.m
CL_tf = ones(size(phi_w,1),1).*C_L;

%% Drag
% Fuselage drag is already included in the drag polar (CD, CL)
CD_tot = CD_data;

%% CD_tf
CD_tf = interp1(CL_data, CD_tot, CL_tf); % the CL_tf is greater than CL_data from cell 94-106 so our A/C stalls

%% Power
%Ptf = ((0.5*rho*S.*(Vel.^3).*CD_tf) ./ eta)  +6; % 0.63 is for the new propeller, +5W to avionics and losses
aoa = interp1(CL_data, AoA_data, CL_tf);

Pclimb = w.*V_inf.*sind(fp_angle);
thrust = (0.5*rho*S.*(V_inf.^2).*CD_tf) + (w.*(sind(fp_angle)));
P_tf = 0.5*rho*S.*(V_inf.^3).*CD_tf;

for i = 1:length(thrust)
    [pt_eff(i), motor_eta(i), J(i)] = powertrain_eff(thrust(i), V_inf(i), '20x8');
    i = i+1;
end

pt_eff = pt_eff'; % set to zero when CP is less than zero (and very close)
J = J';
motor_eta = motor_eta';

Preq = ((P_tf+Pclimb)./(pt_eff.*motor_eta)) + 6;

end

%% Curvature Function (Old)
function [L,R,k] = curvature(X)
% Radius of curvature and curvature vector for 2D or 3D curve
%  [L,R,Kappa] = curvature(X)
%   X:   2 or 3 column array of x, y (and possibly z) coordinates
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

%% Circumcenter Function (Old)
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

