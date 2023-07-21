function [x_p, y_p, z_p] = path_circle(direc, loiter_radius, num_points, repeat, shape, delta_alt,alt)
% PATH_CIRCLE Summary of this function goes here

% Discretize the loiter circle into x,y points
if strcmpi(shape, 'circle')
    x_p = -loiter_radius.*cos(linspace(0, pi, num_points)');
    y_p = [sqrt(loiter_radius.^2 - x_p.^2); -sqrt(loiter_radius.^2 - x_p(2:end).^2)];
    x_p = [x_p; flipud(x_p(1:end-1))];

    %     x_p = (x_p.*(cosd(45))) - (y_p.*(sind(45)));
    %     y_p = (y_p.*(cosd(45))) + (x_p.*(sind(45)));

    z_p = (delta_alt/loiter_radius).*x_p + alt;

elseif strcmpi(shape, 'racetrack') % no longer bwoken
    const_length = loiter_radius/2;

    x_p1 = linspace(0,const_length,num_points/4)';
    y_p1 = ones(num_points/4,1)*loiter_radius;

    y_p2 = loiter_radius.*cos(linspace(0, pi/2, num_points/2)');
    x_p2 = sqrt(loiter_radius.^2 - y_p2.^2); %sqrt(loiter_radius.^2 - x_p2.^2);

    y_p3 = -loiter_radius.*cos(linspace(0, pi/2, num_points/2)');
    x_p3 = sqrt(loiter_radius.^2 - y_p3.^2); %sqrt(loiter_radius.^2 - x_p2.^2);

    x_p = [x_p1(1:end-1); const_length+x_p2(1:end-1); const_length+flipud(x_p3); flipud(x_p1(1:end-1))];
    x_p = [x_p; -flipud(x_p(1:end-1))]; % prev x_p(2:end-1)
    y_p = [y_p1(1:end-1); y_p2(1:end-1); flipud(y_p3); -y_p1(1:end-1)];
    y_p = [y_p; flipud(y_p(1:end-1))]; % prev y_p(2:end-1)

    z_p = (delta_alt/loiter_radius).*x_p + alt;

elseif strcmpi(shape, 'ellipse')
    a = 1.5.*loiter_radius;
    b = loiter_radius;
    y_p1 = b.*cos(linspace(0, pi/2, num_points/2)');
    x_p1 = sqrt((1-((y_p1.^2)./(b.^2))).*(a^2));

    y_p2 = -b.*cos(linspace(0, pi/2, num_points/2)');
    x_p2 = sqrt((1-((y_p2.^2)./(b.^2))).*(a^2));

    x_p = [x_p1(1:end-1); flipud(x_p2)];
    x_p = [x_p; -flipud(x_p(1:end-1))]; % prev x_p(2:end-1)
    y_p = [y_p1(1:end-1); flipud(y_p2)];
    y_p = [y_p; flipud(y_p(1:end-1))]; % prev y_p(2:end-1)

    %     x_p = (x_p.*(cosd(45))) - (y_p.*(sind(45)));
    %     y_p = (y_p.*(cosd(45))) + (x_p.*(sind(45)));

    z_p = (delta_alt/loiter_radius).*x_p + alt;

else
    disp('Loiter shape does not exist')
end

x_p = repmat(x_p, repeat, 1);
y_p = repmat(y_p, repeat, 1);
z_p = repmat(z_p, repeat, 1);

if strcmpi(direc, 'CCW')
    x_p = flipud(x_p);
    y_p = flipud(y_p);
    z_p = flipud(z_p);
elseif strcmpi(direc, 'CW')

else
    disp('Circuit direction wrong')
end

end

% figure(1)
% clf(1)
% plot3(x_p,y_p,z_p,'-')
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('alt + \Delta z (m)')
