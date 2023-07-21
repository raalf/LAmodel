function [out_score, LOCALdatenum, eta_solar, eta_solar_straightup, mid_p, wvec, sunVec, p_eta, p_req, p_res,...
    net_energy,chargerate,roll, hdg, gspd, Elevation, CL_tf,aoa,mid_zp] = score_path(aircraft, direc, x_p, y_p, air_speed, ...
    wind_speed, wind_from, startLocalDate, timezone, Latitude, Longitude, Altitude, irradiance,mass,TL,...
    powerdrawdata,solarmodel,powertrain_eta,C_L, z_p)

addpath(genpath('Include'))
addpath(genpath('PVLib 1.4 Release'))

if strcmpi(aircraft, 'CREATeV')
    % For CREATeV (includes 5 deg wing tilt + panel angles on airfoil)
    paneltiltA = 11.966;
    paneltiltB = 18.378;

    % paneltiltA = 0;
    % paneltiltB = 0;

    area_array = 1.48743; % 96 * (125 * 125) - (8.09 * 8.09 * 2)
elseif strcmpi(aircraft, 'George')
    paneltiltA = 10;
    paneltiltB = 12;
    paneltiltC = 15;

    area_array = 2.2;
end

if isnan(air_speed) % calculate changing airspeed
    [roll, hdg, time_datenum, trk, gspd, mid_p, avec, gvec, wvec, dt, p_req, aoa,...
        air_speed, CL_tf, fp_angle, mid_zp] = path_to_att_mc4(direc, x_p, y_p, wind_speed, ...
        wind_from, C_L, mass,powerdrawdata, powertrain_eta, z_p);

elseif isnan(C_L) % assuming constant airspeed
    % Convert Flight Path to Roll/Yaw based on airspeed and wind
    [roll, hdg, time_datenum, trk, gspd, mid_p, avec, gvec, wvec, dt, fp_angle,mid_zp] = ...
        path_to_att(direc, x_p, y_p, air_speed, wind_speed, wind_from, z_p);

    % Determine AoA and power required for turning flight
    [p_req, aoa, CL_tf] = powerdraw_mc4(roll, air_speed, mass,powerdrawdata,powertrain_eta, fp_angle);

end

if any(isnan(p_req))
    warning('Aircraft stalled')
end

% aoa = [zeros(size(roll,1),1)];
% aoa = aoa-2;

if strcmpi(aircraft,'CREATeV')
    % Calculate PV Panel vector based on roll and yaw, assume zero pitch
    [panelVecA] = attitude_to_panel_vector(paneltiltA, roll, fp_angle+aoa, hdg);
    [panelVecB] = attitude_to_panel_vector(paneltiltB, roll, fp_angle+aoa, hdg);
elseif strcmpi(aircraft,'George')
    % Calculate PV Panel vector based on roll and yaw, assume zero pitch
    [panelVecA] = attitude_to_panel_vector(paneltiltA, roll, fp_angle+aoa, hdg);
    [panelVecB] = attitude_to_panel_vector(paneltiltB, roll, fp_angle+aoa, hdg);
    [panelVecC] = attitude_to_panel_vector(paneltiltC, roll, fp_angle+aoa, hdg);
end

% Prepare Time Vectors
LOCALdatenum = startLocalDate+(time_datenum); % need to check
% GMTdatenum = LOCALdatenum-(timezone/24); % changing sun elevation for loiter
GMTdatenum = (ones(length(LOCALdatenum),1)).*(startLocalDate-(timezone/24)); % fixed sun elevation for loiter

% Calculate Solar Position (Alton's Model)
% [~, ~, SunAz, Elevation, sunVec] = solarPosition(GMTdatenum, Latitude, Longitude);
Time = pvl_maketimestruct(GMTdatenum, 0);
Location = pvl_makelocationstruct(Latitude, Longitude, Altitude);

% Calculate Solar Position (Ephemeris)
% [~, ~, ApparentSunEl, ~] = pvl_ephemeris(Time, Location); %%%%%
[SunAz, Elevation, ApparentSunEl, ~]= pvl_ephemeris(Time, Location);
% [SunAz, Elevation, ApparentSunEl] = pvl_spa(Time, Location);
SunZen = 90-Elevation;
[x,y,z] = sph2cart(pi/2-(deg2rad(SunAz)), deg2rad(Elevation), 1); % Convert elevation and azimuth angle (pi/2-azimuth)
sunVec = [-x,-y,-z]; % Create unit vector from position of the sun to origin

wvec = [wvec, zeros(size(wvec,1),1)];

% Calculate Solar Eff.
[ ~, eta_solar_straightup ] = panelTiltAngle( sunVec, [0,0,1] );
if strcmpi(aircraft,'CREATeV')
    [ AOI_A, eta_solarA ] = panelTiltAngle( sunVec, panelVecA );
    [ AOI_B, eta_solarB ] = panelTiltAngle( sunVec, panelVecB ); % used in model 2, should we be using spherical trig?
    [eta_solar] = min([eta_solarA, eta_solarB],[],2); % try using average
    % [eta_solar] = mean([eta_solarA, eta_solarB],2);
    eta_solar(Elevation<0.0) = 0.0; % can't charge before sunrise or after sunset, even if we are rolled. Update this when we add diffusion
    tiltangle = max(AOI_A,AOI_B);
    % tiltangle = AOI_A;
    AOI = tiltangle;
elseif strcmpi(aircraft,'George')
    [ AOI_A, eta_solarA ] = panelTiltAngle( sunVec, panelVecA );
    [ AOI_B, eta_solarB ] = panelTiltAngle( sunVec, panelVecB ); % used in model 2, should we be using spherical trig?
    [ AOI_C, eta_solarC ] = panelTiltAngle( sunVec, panelVecA );
    eta_solarA(Elevation<0.0) = 0.0; % can't charge before sunrise or after sunset, even if we are rolled. Update this when we add diffusion
    eta_solarB(Elevation<0.0) = 0.0;
    eta_solarC(Elevation<0.0) = 0.0;
    eta_solar = min([eta_solarA, eta_solarB, eta_solarC],[],2);
end

%%
if solarmodel == 1 % using ground station measurements assuming no diffused component
    % p_eta =  eta_solar .* area_array .* irradiance .* 0.195; % maybe (no diffused)
    p_eta =  (cosd(AOI)).* area_array .* irradiance .* 0.195;
    p_eta(Elevation<0.0)=0.0;
    p_eta(AOI>90.0)=0.0;

elseif solarmodel == 2 % using sandia predictions
    pressure = 101325; %could change

    DOY = pvl_date2doy(Time.year, Time.month, Time.day);
    HExtra = pvl_extraradiation(DOY);

    % [SurfAz, SurfEl, ~] = cart2sph(panelVecB(:,2),panelVecB(:,1),panelVecB(:,3));
    % [SurfAz, SurfEl, ~] = cart2sph(panelVecA(:,2),panelVecA(:,1),panelVecA(:,3));
    if strcmpi(aircraft,'CREATeV')
        panelVecavg_1 = mean([panelVecA(:,2),panelVecB(:,2)],2);
        panelVecavg_2 = mean([panelVecA(:,1),panelVecB(:,1)],2);
        panelVecavg_3 = mean([panelVecA(:,3),panelVecB(:,3)],2);
        [SurfAz, SurfEl, ~] = cart2sph(panelVecavg_1,panelVecavg_2,panelVecavg_3);
    
        SurfAz = wrapTo360(rad2deg(SurfAz));
        SurfEl = rad2deg(SurfEl);
        SurfTilt = 90-SurfEl;
    elseif strcmpi(aircraft,'George')
        [SurfAzA, SurfElA, ~] = cart2sph(panelVecA(:,2),panelVecA(:,1),panelVecA(:,3));
        [SurfAzB, SurfElB, ~] = cart2sph(panelVecB(:,2),panelVecB(:,1),panelVecB(:,3));
        [SurfAzC, SurfElC, ~] = cart2sph(panelVecC(:,2),panelVecC(:,1),panelVecC(:,3));
    
        SurfAzA = wrapTo360(rad2deg(SurfAzA));
        SurfAzB = wrapTo360(rad2deg(SurfAzB));
        SurfAzC = wrapTo360(rad2deg(SurfAzC));

        SurfElA = rad2deg(SurfElA);
        SurfElB = rad2deg(SurfElB);
        SurfElC = rad2deg(SurfElC);

        SurfTiltA = 90-SurfElA;
        SurfTiltB = 90-SurfElB;
        SurfTiltC = 90-SurfElC;
    end

    % Compute Clear Sky Irrandances using Ineichen model
    if isnan(TL)
        [GHI, ClearSkyDNI, ClearSkyDHI]= pvl_clearsky_ineichen(Time, Location);
    else
        load(['Required Data' filesep 'LinkeTurbidities.mat']);
        LatitudeIndex = round(LinearlyScale(Latitude, 90, -90, 1, 2160));
        LongitudeIndex = round(LinearlyScale(Longitude, -180, 180, 1, 4320));
        Lookup3D = @(array,a,b,c) array(((c-1).*numel(array(:,:,1))+(b-1).*numel(array(:,1,1))+(a-1)+1));
        L1 = Lookup3D(LinkeTurbidity, LatitudeIndex, LongitudeIndex, Time.month);
        TL = double(L1)./20;
        [~, ClearSkyDNI, ClearSkyDHI]= pvl_clearsky_ineichen(Time, Location, TL);
    end

    %      DHI = irradiance ./ 1050 .* ClearSkyDHI; %apply mu and sigma
    %      DNI = irradiance ./ 1050 .* ClearSkyDNI; %(irradiance is supposed to be based on location values and run multiple times)
    DHI = ClearSkyDHI;
    DNI = ClearSkyDNI;

    AppSunZen = 90-ApparentSunEl;

    AMrelative = pvl_relativeairmass(AppSunZen); % doesn't result in significant changes to model results
    AMa = pvl_absoluteairmass(AMrelative, pressure);
    
    if strcmpi(aircraft,'CREATeV')
    [SkyDiffuse,SkyDiffuse_Iso,SkyDiffuse_Cir,SkyDiffuse_Hor] = pvl_perez(SurfTilt, SurfAz, DHI, DNI, HExtra, AppSunZen, SunAz, AMa); % where morning and afternoon differ

    Eb = DNI.*(cosd(AOI)); % DNI from apparent sun and AOI from Alton
    Eb(isnan(Eb)) = 0;
    Eb(Elevation<0.0)=0.0;
    Eb(AOI>90.0)=0.0;

    %     Flux = Eb;
    Flux = Eb + SkyDiffuse;
    % Flux = SkyDiffuse;
    % Flux = SkyDiffuse_Iso;
    % Flux = Eb + SkyDiffuse_Cir;
    % Flux = SkyDiffuse_Hor;
    %     cell_eta_data = load('cell_eta_data.mat');
    %     cell_eta = interp1(cell_eta_data.cell_eta.x,cell_eta_data.cell_eta.eff,Eb);
    cell_eta = (ones(length(Flux),1)).*0.195;
    %     cell_eta(AOI>mean(AOI)) = 0.19;
    p_eta = Flux.*cell_eta.*area_array; % 0.195

    elseif strcmpi(aircraft,'George')
    [SkyDiffuseA,~,~,~] = pvl_perez(SurfTiltA, SurfAzA, DHI, DNI, HExtra, AppSunZen, SunAz, AMa); % where morning and afternoon differ
    [SkyDiffuseB,~,~,~] = pvl_perez(SurfTiltB, SurfAzB, DHI, DNI, HExtra, AppSunZen, SunAz, AMa); % where morning and afternoon differ
    [SkyDiffuseC,~,~,~] = pvl_perez(SurfTiltC, SurfAzC, DHI, DNI, HExtra, AppSunZen, SunAz, AMa); % where morning and afternoon differ


    EbA = DNI.*(cosd(AOI_A)); % DNI from apparent sun and AOI from Alton
    EbA(isnan(EbA)) = 0;
    EbA(Elevation<0.0)=0.0;
    EbA(AOI_A>90.0)=0.0;

    EbB = DNI.*(cosd(AOI_B)); % DNI from apparent sun and AOI from Alton
    EbB(isnan(EbB)) = 0;
    EbB(Elevation<0.0)=0.0;
    EbB(AOI_B>90.0)=0.0;

    EbC = DNI.*(cosd(AOI_C)); % DNI from apparent sun and AOI from Alton
    EbC(isnan(EbC)) = 0;
    EbC(Elevation<0.0)=0.0;
    EbC(AOI_C>90.0)=0.0;

    FluxA = EbA + SkyDiffuseA;
    FluxB = EbB + SkyDiffuseB;
    FluxC = EbC + SkyDiffuseC;

    cell_eta = (ones(length(FluxA),1)).*0.21;

    p_etaA = FluxA.*cell_eta.*(1/3)*area_array;
    p_etaB = FluxB.*cell_eta.*(1/3)*area_array;
    p_etaC = FluxC.*cell_eta.*(1/3)*area_array;
    
    p_eta = p_etaA + p_etaB + p_etaC;
    end

    

elseif solarmodel == 3 %Using ground station data

    Time = pvl_maketimestruct(GMTdatenum, 0);
    % Location = pvl_makelocationstruct(Latitude, Longitude, FMT.ORGN.Alt(end));
    [~, Zenith, Heading, Elevation, sunVec] = solarPosition(GMTdatenum, Latitude, Longitude);


    [SurfAz, SurfEl, ~] = cart2sph(panelVecB(:,2),panelVecB(:,1),panelVecB(:,3));

    SurfAz = wrapTo360(rad2deg(SurfAz));
    SurfEl = rad2deg(SurfEl);
    SurfTilt = 90-SurfEl;

    % Extraterrestrial radiation based on day of the year
    DOY = pvl_date2doy(Time.year, Time.month, Time.day);
    HExtra = pvl_extraradiation(DOY);
    GHI = ones(size(Zenith,1),1).*irradiance;
    GHI(Elevation<10.0) = 0; %ground station sensor data is bad when sun below 10deg
    [DNI2, DHI2, ~] = pvl_reindl_2(GHI, Zenith, DOY);
    % pressure = interp1(FMT.BARO.TimeS,FMT.BARO.Press,FMT.ATT.TimeS);
    pressure = 101325;
    % temppressure = rmmissing(pressure);
    % pressure(isnan(pressure)) = temppressure(1); %overwrite nan pressures with start of flight?

    AMrelative = pvl_relativeairmass(90-Elevation); %should be apparent sun
    AMa = pvl_absoluteairmass(AMrelative, pressure);

    DHI2(Elevation<0.0) = NaN;
    DNI2(Elevation<0.0 | DNI2<0.0) = NaN;
    [SkyDiffuse,~,~,~] = pvl_perez(SurfTilt, SurfAz, DHI2, DNI2, HExtra, Zenith, Heading, AMa);

    % total solar power income
    p_eta = (eta_solar .* area_array .* DNI2 .*0.195) + (SkyDiffuse .* area_array .* 0.195);

elseif solarmodel == 4 % edwards

    [I, I_direct, I_diffuse]= edwards_clearsky(Time, Elevation, Altitude);

    p_eta =  (cosd(AOI)).* area_array .* I .* 0.195;
    p_eta(Elevation<0.0)=0.0;
    p_eta(AOI>90.0)=0.0;

else

    error('Solar model must be 1 or 2 or 3 or 4');

end

% p_req(:) = ptf(:) + pclimb(:);
p_req(p_req<0) = 0; % power required shouldn't be below zero
p_res(:) = p_eta(:) - p_req(:);  % For each occurrance of bank angle wrt alpha calculate resulting power of p_eta - ptf


% This calculates the difference between the plane and the baseline solar
% efficiencies, integrated over the loiter circle
% % out_score = trapz(LOCALdatenum, eta_solar_straightup - eta_solar)/(max(LOCALdatenum)-min(LOCALdatenum));
net_energy = (trapz(LOCALdatenum, p_res)) * 24;    %time at each point multiplied by power at that point, cummulative

chargerate = net_energy / ((LOCALdatenum(end) - startLocalDate) *24);
out_score = chargerate;
end

function OutputMatrix = LinearlyScale(inputmatrix, inputmin, inputmax, outputmin, outputmax)
% OutputMatrix = LinearlyScale(inputmatrix, inputmin, inputmax, outputmin, outputmax)
% Linearly scales the inputmatrix. Maps all values from inputmin to
% outputmin, and from inputmax to outputmax. Linear mapping from one point
% to the other.
inputrange = inputmax-inputmin;
outputrange  = outputmax-outputmin;

OutputMatrix = (inputmatrix-inputmin)*outputrange/inputrange + outputmin;

end
