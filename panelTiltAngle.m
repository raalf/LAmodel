function [ tiltAngle, eta_solar ] = panelTiltAngle( sunVec, panelVec )
% PANELTILTANGLE Summary of this function goes here
% Calculates angle between vectors using acos of dot product
% Also calculates solar efficiency (loss) using angle between vectors (https://en.wikipedia.org/wiki/Solar_tracker)

% Compare sizes of both input vectors
[sr, sc] = size(sunVec);
[pr, pc] = size(panelVec);


if sc == 3 && pc == 3 && (sr == 1 || pr == 1) || (sc == 3 && pc == 3)
    if pr > 1 && sr == 1
        % single sunVec, multiple panelVec
        sunVec = repmat(sunVec,pr,1);
        
    elseif sr > 1 && pr == 1
        % multiple sunVec, single panelVec
        panelVec = repmat(panelVec,sr,1);
    end
    
    
    % flip panel vec to ensure angle is correct
    panelVec = panelVec .* -1;
    
    % normalize both vectors
    sunVecNorm = sunVec ./ sqrt(sum(abs(sunVec).^2,2));
    panelVecNorm = panelVec ./ sqrt(sum(abs(panelVec).^2,2));
    
    % angle between two vectors
    tiltAngle = acosd(dot(sunVecNorm,panelVecNorm,2));
    
    loss = 1-cosd(tiltAngle); % https://en.wikipedia.org/wiki/Solar_tracker
    % loss of solar due to misallignment of panel to sun
    eta_solar = max(1-loss,0); % clamp efficiency = 1 - loss between 0 to 1
    
else
    warning('Both sunVec and panelVec should have 3 columns. Only one of them can have more than one row');
    
end



end

