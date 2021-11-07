function [ Dist ] = Dist( x0, y0, x, y )
% Distance between two points using simple coordinate geometry.
Dist = sqrt((x0-x)^2+(y0-y)^2); % in meters
Dist = Dist/1852; % in nautical miles (nm)
end

