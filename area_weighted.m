function [S] = area_weighted(xx, yy)
% AREA_WEIGHTED
% compute area for each grid with a given resolution
R = 6378.1 * 1000;   % unit: m
S = nan(xx,yy);
for x=1:xx
    % angle 2 radian
    lat2 = 90 - (x-1)*(180/xx);
    lat1 = 90 - x*(180/xx);
    a2 = lat2 *pi /180;
    a1 = lat1 *pi /180;
    S(x,:) = repmat(R*R*((180/xx)*pi/180)*(sin(a2)-sin(a1)),1,yy);
    clear lat2 lat1 a2 a1
end
end

