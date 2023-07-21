function [panelVec] = attitude_to_panel_vector(paneltilt, roll, pitch, yaw)
%ATTITUDE_TO_PANEL_VECTOR Summary of this function goes here
%   Detailed explanation goes here

% paneltilt: how many degrees do the pv panels lean back on the wing

qA = quaternion([-yaw-90 pitch+paneltilt roll*-1],'eulerd','zyx','frame');
panelRotMatA = rotmat(qA,'frame');

% panelVec2 = ([0;0;1].*panelRotMat);
panelVec = (repmat([0;0;1], 1,1,size(qA,1)).*panelRotMatA);
panelVec = squeeze(panelVec(3,:,:,:))';

end

