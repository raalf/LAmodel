%% Power Prediction for turning flight (To integrate into optimization code)

function [Preq, aoa, CL_tf, fp_a] = powerdraw_mc4(roll, air_speed, mass,powerdrawdata, eta, fp_angle, V_g)
% Inputs:
%       * Roll Angle (in degrees)
%       * Airspeed: Always use true airspeed (m/s)
%       * Altitude: in meters
% Loads data and assigns names to each coloumn needed
% if nargin==6
    CreateVData_FLT = powerdrawdata.data;
% else
%     load('simpleModel_lineaerCLa.mat')
% end

velocities = CreateVData_FLT(:,1);
CD_data = CreateVData_FLT(:,2);
CL_data = CreateVData_FLT(:,3);
AoA_data = CreateVData_FLT(:,4);

% velocities = CreateVData_FLT(:,1);
% CD_data = CreateVData_FLT.CD;
% CL_data = CreateVData_FLT.CL;
% AoA_data = CreateVData_FLT.aoa;

w = mass * 9.81; 
b = 6.28;
c = 0.395;
S = b * c;

rho = 1.1483;
Vel = air_speed;

%% CLtf for Roll/Bank Angle 
phi = (roll); 
cl_temp = w ./(0.5*rho*S*(Vel.^2)); 

% CL_tf = cl_temp ./ cosd(phi); 

for i = 1:1:size(phi, 1)
     CL_tf(i,:) = (cl_temp)./ cosd(phi(i)); % Calculated CL_tf data for each roll/bank angle using the CL data
end

%% Drag
% Fuselage drag is already included in the drag polar (CD, CL)
CD_tot = CD_data;
%% CD_tf 
CD_tf = interp1(CL_data, CD_tot, CL_tf); % the CL_tf is greater than CL_data from cell 94-106 so our A/C stalls
%% Power
% Ptf = ((0.5*rho*S.*(Vel.^3).*CD_tf) ./ eta)  +6; % 0.63 is for the new propeller, +5W to avionics and losses
% account for accelerations? mgh + 0.5mv^2

aoa = interp1(CL_data, AoA_data, CL_tf);
Pclimb = w.*V_g.*tand(fp_angle);
sinkrate = V_g.*tand(fp_angle);
fp_a = asind(sinkrate./Vel);
thrust = (0.5*rho*S.*(Vel.^2).*CD_tf) + (Pclimb./Vel);
P_tf = 0.5*rho*S.*(Vel.^3).*CD_tf;

for i = 1:length(thrust)
    [pt_eff(i), motor_eta(i), J(i), n(i)] = powertrain_eff(thrust(i), Vel, '20x8');
    i = i+1;
end

pt_eff = pt_eff'; % set to zero when CP is less than zero (and very close)
J = J';
motor_eta = motor_eta';
n = n';

% pt_eff = 0.59;
% motor_eta = 1;

Preq = ((P_tf+Pclimb)./(pt_eff.*motor_eta)) + 6;
% Preq = ((P_tf+Pclimb)./(eta)) + 6;

end

