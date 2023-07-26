% Evaluate loiter score for one loiter given time & location
clc
close all
clear

aircraft = 'CREATeV';
% aircraft = 'George';
shape = 'circle';
% shape = 'ellipse';
% shape = 'racetrack'; % circle, elliptical, racetrack
direc = 'CW';
num_points = 1500; % 1500;
delta_alt = -25; % m
air_speed = 10;
C_L = NaN;
% air_speed = NaN;
% C_L = 1;

wind_speed = 8;
wind_from = 270;
loiter_radius = 500;

Altitude = 350;
mass = 12.685;
powerdrawdata = load('simpleModel_lineaerCLa.mat');

SolarModel = 2;
timezone = -4;
Latitude = 44.04548765; 
Longitude = -79.84323204;

irradiance = 1000;

% startLocalDate = datenum(2022,06,21,05,42,54); % sunrise 0 deg sun el
% startLocalDate = datenum(2022,06,21,20,59,34); % sunset 0 deg sun el
% startLocalDate = datenum(2022,06,21,06,40,00); % morning 9 deg sun el
% startLocalDate = datenum(2022,06,21,06,45,00); % morning 10 deg sun el
% startLocalDate = datenum(2022,06,21,19,57,00); % afternoon 10 deg sun el
% startLocalDate = datenum(2022,06,21,07,44,00); % morning 20 deg sun el
% startLocalDate = datenum(2022,06,21,18,58,00); % afternoon 20 deg sun el
% startLocalDate = datenum(2022,06,21,08,13,00); % morning sun 25 deg el (charging starts)
% startLocalDate = datenum(2022,06,21,18,29,25); % afternoon sun 25 deg el (chargerate ends)
% startLocalDate = datenum(2022,06,21,08,40,00); % morning 30 deg sun el
% startLocalDate = datenum(2022,06,21,18,02,00); % afternoon 30 deg sun el
% startLocalDate = datenum(2022,06,21,09,36,46); % morning 40 deg sun el
startLocalDate = datenum(2022,06,21,17,05,42); % afternoon 40 deg sun el
% startLocalDate = datenum(2022,06,21,10,33,00); % morning 50 deg sun el
% startLocalDate = datenum(2022,06,21,16,10,00); % afternoon 50 deg sun el
% startLocalDate = datenum(2022,06,21,11,35,00); % morning 60 deg sun el
% startLocalDate = datenum(2022,06,21,15,08,00); % afternoon 60 deg sun el
% startLocalDate = datenum(2022,06,21,13,21,15); % solar noon ~70 deg sun el

% startLocalDate = datenum(2022,06,21,17,05,42); 
% startLocalDate = datenum(2022,06,21,18,29,25);

%%
[out_score, LOCALdatenum, eta_solar, eta_solar_straightup, mid_p, wvec, sunVec, p_eta, p_req, p_res, net_energy,...
    chargerate, roll, hdg, gspd, Elevation, CL_tf, aoa, mid_zp] = loiter_score(aircraft, direc, loiter_radius, num_points,...
    air_speed, wind_speed, wind_from, startLocalDate, timezone, Latitude, Longitude, Altitude, irradiance, 1,...
    mass, NaN, powerdrawdata, SolarModel, 0.59,C_L,shape,delta_alt/2);

p_solar = ((trapz(LOCALdatenum, p_eta)) * 24)/ ((LOCALdatenum(end) - startLocalDate) *24)
p_flight = ((trapz(LOCALdatenum, p_req)) * 24)/ ((LOCALdatenum(end) - startLocalDate) *24)
chargerate
% chargerate = net_energy / ((LOCALdatenum(end) - LOCALdatenum(1)) *24); 
%% Plotting
str_title = sprintf("%s \n Airspeed: %.1f m/s \n Wind: %.1f m/s from %.f \n %s Loiter Radius: %.f m \n Charge Rate (Wh/h): %f", ...
              datestr(startLocalDate, 'dd-mmm-yyyy HH:MM'), air_speed, wind_speed, wind_from, direc, loiter_radius, chargerate);
str_title2 = sprintf("Elevation: %.1f deg \n Airspeed: %.1f m/s \n Wind: %.1f m/s from %.f \n %s Loiter Radius: %.f m \n Charge Rate (Wh/h): %f \n Altitude Delta: %.f (m)", ...
              Elevation(1), air_speed, wind_speed, wind_from, direc, loiter_radius, chargerate, delta_alt);
str_title3 = sprintf("Elevation: %.1f deg \n Airspeed: %.1f m/s \n Wind: %.1f m/s from %.f \n %s Loiter Radius: %.f m \n Solar Power Input (Wh/h): %f \n Altitude Delta: %.f (m)", ...
              Elevation(1), air_speed, wind_speed, wind_from, direc, loiter_radius, p_solar, delta_alt);
str_title4 = sprintf("Elevation: %.1f deg \n Airspeed: %.1f m/s \n Wind: %.1f m/s from %.f \n %s Loiter Radius: %.f m \n Power Req Output (Wh/h): %f \n Altitude Delta: %.f (m)", ...
              Elevation(1), air_speed, wind_speed, wind_from, direc, loiter_radius, p_flight, delta_alt);

%% Efficiency
% hFig1 = figure(4);
% clf(4);
% 
% plot(LOCALdatenum, eta_solar_straightup)
% hold on
% plot(LOCALdatenum, eta_solar)
% hold off
% datetick('x', 'HH:MM:SS')
% title(str_title)
% grid minor
% box on
% axis tight
% legend('Baseline', 'Aircraft', 'Location', 'Best')
% xlabel('Time');
% ylabel('Solar Efficiency')

%% 2D Efficiency
% hFig12 = figure(12);
% clf(12);
% hFig12.Units = "centimeters";
% hFig12.Position = [18 2 20 14];
% 
% plot(mid_p(:,1), mid_p(:,2), '-k');
% hold on
% q_scale = 0.1.*loiter_radius;
% q1 = quiver(0, 0, wvec(1,1).*wind_speed.*q_scale, wvec(1,2).*wind_speed.*q_scale, 1, 'b');
% q1.AutoScale = 'off';
% q1.MaxHeadSize = 1;
% 
% sunscale = loiter_radius.*mean(eta_solar_straightup);
% q2 = quiver(0, 0, sunVec(1,1).*sunscale, sunVec(1,2).*sunscale, 0, 'color',[255 69 0]./255);
% q2.AutoScale = 'off';
% q2.MaxHeadSize = 1;
% 
% scatter(mid_p(:,1), mid_p(:,2), 100, p_res, 'filled');
% 
% mk = 'v';
% if strcmpi(direc,'CW')
%     mk = '^';
% end
% scatter(mid_p(1,1), mid_p(1,2), 500, mk,'m','filled')
% colormap(jet)
% colorbar;
% h = colorbar;
% ylabel(h, 'Chargerate (W)', 'fontweight', 'bold','FontSize',14)
% 
% xlabel('\Deltax (m)')
% ylabel('\Deltay (m)')
% lgd = legend('','Wind Vector','Sun Vector','','Direction of Flight');
% lgd.Location = 'southoutside';
% lgd.FontSize = 14;
% lgd.NumColumns = 3;
% 
% hold off
% % grid on
% % grid minor
% % box on
% % axis tight
% axis equal
% axis  off
% 
% 
% % exportgraphics(hFig12,'bestCR.pdf','ContentType','vector')
% 
% % title(str_title);
% % title(str_title2);

%% Power Comparison
% hFig2 = figure(2);
% clf(2);
% hFig2.Units = "centimeters";
% hFig2.Position = [7 2 26 14];
% 
% plot(LOCALdatenum, p_res, Color='blue')
% hold on
% plot(LOCALdatenum, p_req, Color='black')
% plot(LOCALdatenum, p_eta, Color='red')
% hold off
% datetick('x', 'HH:MM:SS')
% % title(str_title)
% grid minor
% box on
% axis tight
% lgd = legend('Batt Power', 'Power Req OUT', 'Solar Power IN', 'Location', 'southoutside');
% % lgd = legend('Batt Power', 'Power Req OUT', 'Solar Power IN', 'Location', 'southoutside');
% % lgd.NumColumns = 3;
% xlabel('Time');
% ylabel('Power (W)')
% % exportgraphics(hFig2,'changingalt_3mps_west.pdf','ContentType','vector')

% % hFig2 = figure(2);
% % clf(2);
% % hFig2.Units = "centimeters";
% % hFig2.Position = [7 2 26 14];
% % 
% % plot(LOCALdatenum_plus, p_res_plus, Color='blue')
% % hold on
% % plot(LOCALdatenum_minus, p_res_minus, Color='blue',LineWidth=3)
% % plot(LOCALdatenum_plus, p_req_plus, Color='black')
% % plot(LOCALdatenum_minus, p_req_minus, Color='black',LineWidth=3)
% % plot(LOCALdatenum_plus, p_eta_plus, Color='red')
% % plot(LOCALdatenum_minus, p_eta_minus, Color='red',LineWidth=3)
% % hold off
% % datetick('x', 'HH:MM:SS')
% % % title(str_title)
% % grid minor
% % box on
% % axis tight
% % lgd = legend('\Delta+25 Batt Power', '\Delta-25 Batt Power', '\Delta+25 Power Req OUT', '\Delta-25 Power Req OUT', '\Delta+25 Solar Power IN', '\Delta-25 Solar Power IN', 'Location', 'southoutside');
% % % lgd = legend('Batt Power', 'Power Req OUT', 'Solar Power IN', 'Location', 'southoutside');
% % lgd.NumColumns = 3;
% % xlabel('Time');
% % ylabel('Power (W)')
% % % exportgraphics(hFig2,'changingalt_3mps_west.pdf','ContentType','vector')

%% Path Efficiency Visualization
% hFig4 = figure(1);
% clf(1);
% hFig4.Units = "centimeters";
% hFig4.Position = [18 2 20 14];
% 
% plot3(mid_p(:,1), mid_p(:,2), mid_zp(:,1), '-k');
% hold on
% q_scale = 0.1.*loiter_radius;
% q1 = quiver3(0, 0, Altitude, wvec(1,1).*wind_speed.*q_scale, wvec(1,2).*wind_speed.*q_scale, 1, 'b');
% q1.AutoScale = 'off';
% q1.MaxHeadSize = 1;
% 
% sunscale = loiter_radius.*mean(eta_solar_straightup);
% q2 = quiver3(0, 0, Altitude, sunVec(1,1).*sunscale, sunVec(1,2).*sunscale, 0, 'color',[255 69 0]./255);
% q2.AutoScale = 'off';
% q2.MaxHeadSize = 1;
% 
% scatter3(mid_p(:,1), mid_p(:,2), mid_zp(:,1), 100, p_res, 'filled');
% 
% mk = 'v';
% if strcmpi(direc,'CW')
%     mk = '^';
% end
% % scatter3(mid_p(1,1), mid_p(1,2),Altitude, 500, mk,'m','filled')
% colormap(jet)
% colorbar;
% h = colorbar;
% ylabel(h, 'Chargerate (W)', 'fontweight', 'bold','FontSize',14)
% 
% xlabel('\Deltax (m)')
% ylabel('\Deltay (m)')
% zlabel('Altitude (SL) + \Deltaz (m)')
% % lgd = legend('','Wind Vector','Sun Vector','','Direction of Flight');
% % lgd.Location = 'southoutside';
% % lgd.FontSize = 14;
% % lgd.NumColumns = 2;
% 
% hold off
% grid on
% % grid minor
% box on
% % axis tight
% axis equal
% % axis  off
% % axis off
% % box off
% view(0,90)
% 
% % exportgraphics(hFig4,'CR_25m_racetrack_9mps90deg.pdf','ContentType','vector')
% 
% % title(str_title);
% title(str_title2);

%% Solar Power Input Visualization
% hFig5 = figure(4);
% clf(4);
% hFig5.Units = "centimeters";
% hFig5.Position = [10 2 20 14];
% 
% plot3(mid_p(:,1), mid_p(:,2), mid_zp(:,1), '-k');
% hold on
% q_scale = 0.1.*loiter_radius;
% q1 = quiver3(0, 0, Altitude, wvec(1,1).*wind_speed.*q_scale, wvec(1,2).*wind_speed.*q_scale, 1, 'b');
% q1.AutoScale = 'off';
% q1.MaxHeadSize = 1;
% 
% sunscale = loiter_radius.*mean(eta_solar_straightup);
% q2 = quiver3(0, 0, Altitude, sunVec(1,1).*sunscale, sunVec(1,2).*sunscale, 0, 'color',[255 69 0]./255);
% q2.AutoScale = 'off';
% q2.MaxHeadSize = 1;
% 
% scatter3(mid_p(:,1), mid_p(:,2), mid_zp(:,1), 100, p_eta, 'filled');
% 
% mk = 'v';
% if strcmpi(direc,'CW')
%     mk = '^';
% end
% % scatter3(mid_p(1,1), mid_p(1,2),Altitude, 500, mk,'m','filled')
% colormap("hot")
% colorbar;
% h = colorbar;
% ylabel(h, 'Solar Power In (W)', 'fontweight', 'bold','FontSize',14)
% 
% xlabel('\Deltax (m)')
% ylabel('\Deltay (m)')
% zlabel('Altitude (SL) + \Deltaz (m)')
% 
% hold off
% grid on
% box on
% axis equal
% axis off
% box off
% view(0,90)

% exportgraphics(hFig5,'Sol_n25m_circ_xy.pdf','ContentType','vector')

% title(str_title3);

% figure(69)
% clf(69)
% plot(LOCALdatenum,hdg)
% xlabel('date')
% ylabel('heading')
% datetick('x', 'HH:MM:SS')
% box on
% axis tight
% grid on

%% Power Required Output Visualization
% hFig8 = figure(8);
% clf(8);
% hFig8.Units = "centimeters";
% hFig8.Position = [10 2 20 14];
% 
% plot3(mid_p(:,1), mid_p(:,2), mid_zp(:,1), '-k');
% hold on
% q_scale = 0.1.*loiter_radius;
% q1 = quiver3(0, 0, Altitude, wvec(1,1).*wind_speed.*q_scale, wvec(1,2).*wind_speed.*q_scale, 1, 'b');
% q1.AutoScale = 'off';
% q1.MaxHeadSize = 1;
% 
% sunscale = loiter_radius.*mean(eta_solar_straightup);
% q2 = quiver3(0, 0, Altitude, sunVec(1,1).*sunscale, sunVec(1,2).*sunscale, 0, 'color',[255 69 0]./255);
% q2.AutoScale = 'off';
% q2.MaxHeadSize = 1;
% 
% scatter3(mid_p(:,1), mid_p(:,2), mid_zp(:,1), 100, p_req, 'filled');
% 
% mk = 'v';
% if strcmpi(direc,'CW')
%     mk = '^';
% end
% % scatter3(mid_p(1,1), mid_p(1,2),Altitude, 500, mk,'m','filled')
% colormap("parula")
% colorbar;
% h = colorbar;
% ylabel(h, 'Power Req Out (W)', 'fontweight', 'bold','FontSize',14)
% 
% xlabel('\Deltax (m)')
% ylabel('\Deltay (m)')
% zlabel('Altitude (SL) + \Deltaz (m)')
% 
% hold off
% grid on
% box on
% axis equal
% axis off
% box off
% view(0,90)
% % exportgraphics(hFig8,'40deg_morn_wind9_from270_rad500_CW.png','ContentType','vector')
% 
% % title(str_title4);

%% Path
% hFig24 = figure(24);
% hold on
% % clf(24);
% hFig24.Units = "centimeters";
% hFig24.Position = [10 2 20 12];
% plot3(x_p,y_p,z_p,'-')
% xlabel('x (m)', 'fontweight', 'bold','FontSize',14)
% ylabel('y (m)', 'fontweight', 'bold','FontSize',14)
% zlabel('alt + \Delta z (m)', 'fontweight', 'bold','FontSize',14)
% % grid on
% % title('Circular Flight Path')
% axis tight
% % exportgraphics(figure(24),'2d_allpath.pdf','ContentType','vector')
% lgd = legend('Circle','Ellipse','Racetrack');
% lgd.FontSize = 12;
% lgd.Location = 'southoutside';
% lgd.NumColumns = 3;

%% 3D Bar for Power Result (Power vs. Bank Angle)
% figure(53);
% clf(53)
% 
% hold on
% plot(LOCALdatenum,ptf,'.');
% plot(LOCALdatenum,p_eta,'.');
% plot(LOCALdatenum,p_res,'.');
% 
% title(str_title)
% legend('Power Req','Solar Input','Charging Power');
% % legend('Solar Input');
% datetick('x', 'HH:MM:SS')
% %xlim([0 40]);
% axis tight;
% %ylim([0 70]);
% %zlim([-200 200]);
% ylabel('Power (Watts)', 'fontweight', 'bold');
% xlabel('Local Time', 'fontweight', 'bold');
% zlabel('Power (Watts)', 'fontweight', 'bold');
% grid on;
