function [D, scalp] = electrode_distances(eloc, varargin)
%ELECTRODE_DISTANCES  Distance between every pair of EEG electrodes, in cm along the scalp
%
%   Usage:
%       D = electrode_distances(eloc)
%       [D, scalp] = electrode_distances(eloc, 'head_radius', 9.5, 'nonscalp', {'Chin'})
%
%   Inputs:
%       eloc : 1xC struct - channel locations as used by plot_topo -- required
%
%   Name-Value Pairs:
%       'head_radius' : scalar - radius of the head in cm; 8.9 is a 56 cm
%                       circumference. Pass 1 for the arc in radians on the
%                       unit sphere. (default: 8.9)
%       'nonscalp'    : cellstr - labels of further channels to leave out,
%                       beyond the ones scalp_channels recognizes (default: {})
%       'verbose'     : logical - print the nearest-neighbor summary (default: true)
%
%   Outputs:
%       D     : C x C double - great-circle distance between channels i and
%               j along the scalp, in cm; 0 on the diagonal; NaN in the row
%               and column of every non-EEG channel
%       scalp : C x 1 logical - true for the EEG channels, from scalp_channels
%
%   Notes:
%       The distance is the arc on a sphere of head_radius through the
%       positions plot_topo uses: the azimuthal equidistant projection of
%       xyz_to_eloc, where radius is the polar angle down from the vertex
%       with the equator at 0.5 and theta the azimuth. The head is taken as
%       a sphere, so the centimeters are nominal and scale with head_radius;
%       a head of 53 or 59 cm circumference is about 5% off the default.
%       Electrodes below the equator are handled like any other.
%
%       One line is printed: the median, smallest and largest distance from
%       each EEG electrode to its nearest neighbor. That spacing is the
%       smallest spatial scale the cap resolves and the smallest useful
%       width for smoothing across the scalp.
%
%   Example:
%       eloc = cap_montage('AC-64');
%       D = electrode_distances(eloc);              % median nearest neighbor about 2.5 cm
%       up = triu(true(size(D)), 1);
%       histogram(D(up & isfinite(D)), 0:0.5:20)    % all pair distances
%
%   See also: scalp_channels, plot_topo, xyz_to_eloc, cap_montage
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

p = inputParser;
p.FunctionName = 'electrode_distances';
addRequired(p, 'eloc', @(x) validateattributes(x, {'struct'}, {'vector'}));
addParameter(p, 'head_radius', 8.9, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'finite'}));
addParameter(p, 'nonscalp', {}, @(x) iscellstr(x) || isstring(x) || ischar(x));
addParameter(p, 'verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, eloc, varargin{:});
opt = p.Results;

scalp = scalp_channels(eloc, cellstr(opt.nonscalp));

% Back from the projection to points on the unit sphere, then the arc
% between every pair: acos of the dot product, clipped so rounding never
% makes it complex
polar = [eloc.radius]' * pi;                     % angle down from the vertex, radians (0.5 -> the equator)
azim  = [eloc.theta]' * pi / 180;
xyz = [sin(polar) .* cos(azim), sin(polar) .* sin(azim), cos(polar)];
D = opt.head_radius * acos(min(max(xyz * xyz', -1), 1));
D(1:size(D, 1) + 1:end) = 0;                     % exact zeros on the diagonal
D(~scalp, :) = NaN;
D(:, ~scalp) = NaN;

if ~opt.verbose, return, end
Dnn = D(scalp, scalp);
Dnn(1:size(Dnn, 1) + 1:end) = Inf;
nearest = min(Dnn, [], 2);
fprintf('electrode_distances: %d EEG channel(s) of %d; nearest neighbor median %.2f, smallest %.2f, largest %.2f (head radius %g)\n', ...
    nnz(scalp), numel(scalp), median(nearest), min(nearest), max(nearest), opt.head_radius);
end
