function [pt_eta, motor_eta, J, n] = powertrain_eff(thrust, V_tas, prop)
%% From Propeller Performance code (Will)
% 1. Estimate thrust required
% 2. Use thrust to find RPM using getRPM
% 3. Use RPM to find J
% 4. Plug J back in to get efficiency

% CT = Thrust Coefficient
% CP = Power Coeficient
% eff = overall efficiency for 20x8 system from CP and CT fits

if strcmpi(prop, '20x8')
    diameter = 0.512; % for aeronaut 20x8
elseif strcmpi(prop,'185x12')
    diameter = 0.4699;
end

rho = 1.1483;

[n] = getRPM(thrust, rho, V_tas, diameter, prop);
[Q, J] = getTorque(rho, n, V_tas, diameter, prop);
[motor_eta, curr] = motor_eff(n, Q);

CT = thrust_coeff(J, prop);
CP = power_coeff(J, prop);
pt_eta = J.*(CT./CP);

end

% Overall thrust coefficient with polynomial fit
function [CT] = thrust_coeff(J, prop)
if strcmpi(prop, '20x8')
    CT = (0.1298.*(J.^3)) - (0.2679.*(J.^2)) - (0.02553.*J) + 0.07525;
elseif strcmpi(prop,'185x12')
    CT = (-1.636.*(J.^5)) + (3.933.*(J.^4)) - (3.246.*(J.^3)) + (0.8995.*(J.^2)) - (0.09467.*J) + 0.08651;
end
end

% Overall power coefficient with polynomial fit
function [CP] = power_coeff(J, prop)
if strcmpi(prop, '20x8')
    CP = (-0.1618.*(J.^3)) + (0.0292.*(J.^2)) - (0.002126.*J) + 0.02489;
elseif strcmpi(prop,'185x12')
    CP = (0.2741.*(J.^4)) - (0.5853.*(J.^3)) + (0.3012.*(J.^2)) - (0.05987.*J) + 0.04802;
end
end

%% Finding drag of freewheeling propeller
function [CT_fw, J_new] = freewheel_tcoeff(J_old)
%     J_old = 0.5;
dJ = 0.01;
error_dem = 1e-5;
cp = power_coeff(J_old);

while cp > error_dem
    CP_slope1 = power_coeff(J_old+(0.5.*dJ));
    CP_slope2 = power_coeff(J_old-(0.5.*dJ));
    slope = (CP_slope1-CP_slope2)./dJ;

    CP_J_new = power_coeff(J_old);
    J_new = J_old - (CP_J_new./slope);

    cp = power_coeff(J_new);
    J_old = J_new;
end

CT_fw = thrust_coeff(J_old);
end

%% Finding RPS required for certain thrust @ true airspeed
function [n] = getRPM(thrust, rho, V_tas, diameter, prop)
n_old = 40;
J_old = V_tas./(n_old.*diameter);
dn = 0.2;

% error_dem = ones(length(thrust),1).*0.001; % 1e-5
error_dem = 1e-5;

[CT_T_act] = thrust_coeff(J_old, prop);
T_act = CT_T_act.*rho.*(n_old.^2).*(diameter.^4);
error = abs(thrust-T_act);

%     for i = 1:(length(error))
while error > error_dem
    J_f_old = V_tas./(n_old.*diameter); % Old guess at AR
    [CT_f_old] = thrust_coeff(J_f_old, prop); % Old thrust coefficient
    f_old = CT_f_old.*rho.*(n_old.^2).*(diameter.^4); % Old thrust

    J_slope1 = V_tas./((n_old-(0.5.*dn)).*diameter);
    J_slope2 = V_tas./((n_old+(0.5.*dn)).*diameter);
    [CT_slope1] = thrust_coeff(J_slope1, prop);
    [CT_slope2] = thrust_coeff(J_slope2, prop);
    slope = ((CT_slope1.*rho.*((n_old-(0.5.*dn)).^2).*(diameter.^4))-...
        (CT_slope2.*rho.*((n_old+(0.5.*dn)).^2).*(diameter.^4))).*(dn.^(-1));

    n = n_old - ((thrust-f_old)./slope);

    J_T_act = V_tas./(n.*diameter);
    [CT_T_act] = thrust_coeff(J_T_act, prop);
    T_act = CT_T_act.*rho.*(n.^2).*(diameter.^4); % keeps getting very smol, no convergence

    error = abs(thrust-T_act);
    n_old = n;
end

%     end
end



% Finding torque from propeller RPS and airspeed
function [Q, J] = getTorque(rho, n, V_tas, diameter, prop)
J = V_tas./(n.*diameter);
[CP_CQ] = power_coeff(J, prop);
CQ = CP_CQ.*(sqrt(1/(2*pi)));
Q = CQ.*rho.*(n.^2).*(diameter.^5);
end

%% From Motor Functions code (Will)
function [motor_eta, curr] = motor_eff(n, Q)
% n = rotations per second
% i = current draw (measured by ESC)

curr = (17.98.*(Q.^2)) + (8.493.*Q) + 0.363;

% For Myxa ESC (from fitting)
i00 = 0.1035;
i01 = 0.0003938;
i02 = 1.128e-6;

Kv = 280; % For U7-V2
Kq = 304;
tau = -0.0002069;
R = 0.01541;

omega = 2.*pi.*n;
i0 = i00 + (i01.*omega) + (i02.*(omega.^2));

motor_eta = (1-(i0./curr)).*(Kv/Kq).*((1+(tau.*omega)+((curr.*R.*Kv).*(omega.^(-1)))).^(-1));

end
