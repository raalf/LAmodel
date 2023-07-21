%% Power Prediction for turning flight (To integrate into optimization code)

function[Ptf] = powerdraw(roll, air_speed)

% Loads data and assigns names to each coloumn needed
load('CreateVDataNotTrimmed');
alpha = CreateVDataNotTrimmed(:,1);
CL_data = CreateVDataNotTrimmed(:,3);
% CDi_data = CreateVDataNotTrimmed(:,4);
% CDv_data = CreateVDataNotTrimmed(:,5);
CD_data = CreateVDataNotTrimmed(:,6);
w =(12.6 + (0.05*0*6)) * 9.81;  %add extra batt here!
b = 6.28;
c = 0.395;
S = b * c;
rho = 1.2;
Vel = air_speed;
%% CLtf for Roll/Bank Angle 


phi = (roll); 
cl_temp = w /( 0.5*rho*S.*(Vel^2));

for i = 1:1:size(phi, 1)
    CL_tf(i, :) = (cl_temp)./ cosd(phi(i)); % Calculated CL_tf data for each roll/bank angle using the CL data
end

%% Drag

Df = 0.01;

% Drag to find K
D = 0.5*rho*S*CD_data.*(Vel^2); 
D_tot = Df + D;

CD_tot = 2*D_tot/(rho*S*(Vel^2));

%% CD_tf 

CD_tf = interp1(CL_data, CD_tot, CL_tf);


%% Power


Ptf = ( 0.5*rho*S*(Vel^3).*CD_tf ) ./ 0.61 ;



