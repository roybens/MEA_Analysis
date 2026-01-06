function cmap = customDivergingColorMap(nColors)
% FINALCUSTOMCOLORMAP builds black-blue-white-orange-red with correct anchors
% nColors: number of colors (default 256)

if nargin < 1
    nColors = 256;
end

% [Position (0-1), R, G, B]
controlPoints = [
   
    0.00, 0.0, 0.0, 0.0;  
    0.01,0.0,0.4,0.7;
    0.05, 0.0, 0.5, 1.0;    % Light Blue
    0.75, 1.0, 1.0, 1.0;    % White (center at 75%)
    0.95, 1.0, 0.5, 0.2;    % Light Orange
    1.00,1.0,0.3,0.0;
];

% Split into components
x = controlPoints(:,1);
r = controlPoints(:,2);
g = controlPoints(:,3);
b = controlPoints(:,4);

% Query points
xq = linspace(0,1,nColors);

% Interpolate each channel
rq = interp1(x, r, xq, 'linear');
gq = interp1(x, g, xq, 'linear');
bq = interp1(x, b, xq, 'linear');

% Build colormap
cmap = [rq(:), gq(:), bq(:)];
end