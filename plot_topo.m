function [hax, ex, ey] = plot_topo(values, eloc, varargin)
%PLOT_TOPO  Plot electrode locations (and optional scalp map) on a head model
%
%   Usage:
%       [hax, ex, ey] = plot_topo([], eloc)
%       [hax, ex, ey] = plot_topo([], eloc, 'electrodes', 'ptslabels')
%       [hax, ex, ey] = plot_topo(values, eloc, 'style', 'map', 'conv', 'on')
%
%   Inputs:
%       values : Cx1 double - one value per channel to render as an
%                interpolated scalp map. Pass [] for a head-and-electrodes
%                plot with no map -- required (may be empty)
%       eloc   : 1xC struct - EEGLAB-style channel locations. Must supply
%                either the polar pair (theta, radius) or the cartesian
%                triple (X, Y, Z); polar is used when present. An optional
%                'labels' field (char or numeric) is used by the label
%                display modes -- required
%
%   Name-Value Pairs:
%       'style'      : char - 'map' draws the interpolated scalp map,
%                      'blank' draws only the head and electrodes
%                      (default: 'blank' when values is empty, else 'map')
%       'electrodes' : char - 'off' | 'pts' | 'labels' | 'ptslabels'
%                      (default: 'pts')
%       'maplimits'  : 1x2 double - [lo hi] color limits, or 'absmax' for
%                      symmetric limits about zero (default: 'absmax')
%       'headrad'    : double - drawn radius of the head cartoon. Only 0.5
%                      is anatomically correct; 0 omits the head
%                      (default: 0.5)
%       'gridres'    : integer - side length of the interpolation grid.
%                      Drives the smoothness of the masked map edge
%                      (default: 200)
%       'shading'    : char - 'interp' or 'flat' (default: 'interp')
%       'conv'       : char - 'on' masks the map to the convex hull of the
%                      electrodes, 'off' masks it to the head rim
%                      (default: 'off')
%       'emarker'    : cell - {marker, color, size, linewidth}. Any element
%                      may be [] to keep its default
%                      (default: {'.', 'k', [], 1})
%       'headcolor'  : 1x3 double - color of the head/nose/ear lines
%                      (default: [0 0 0])
%       'parent'     : axes handle - axes to draw into (default: gca)
%
%   Outputs:
%       hax : axes handle - the axes drawn into
%       ex  : Cx1 double - electrode screen x coordinates, in channel order
%       ey  : Cx1 double - electrode screen y coordinates, in channel order
%
%   Notes:
%       This is a self-contained replacement for EEGLAB's topoplot family.
%       It reads an eloc struct directly and owns its head model, so the
%       repository carries no EEGLAB dependency.
%
%       Coordinate convention follows EEGLAB so that existing montages plot
%       unchanged: theta is in degrees, radius is a normalized arc length
%       where 0.5 is the head rim (the ears/eyes equator), and the nose
%       points toward +y on screen. Channel k lands at
%
%           ex(k) = radius(k)*sin(theta(k)*pi/180)
%           ey(k) = radius(k)*cos(theta(k)*pi/180)
%
%       Returning ex/ey is the point of the interface: callers that need
%       electrode positions read them from the second and third outputs
%       rather than scraping XData back off the plotted line objects.
%
%       Coordinates are plotted at their projected radius with no rescaling.
%       Electrodes below the equator have radius > 0.5 and therefore fall
%       outside the drawn head, in the conventional "skirt" -- FT9/FT10 and
%       TP9/TP10 on a 64-channel cap sit at ~0.63 and are genuinely below the
%       ears/eyes line. The axes are sized to keep them visible.
%
%       Channels whose values are NaN or Inf are dropped from the map
%       interpolation but still drawn as electrodes.
%
%   Example:
%       load eloc64;
%       [~, ex, ey] = plot_topo(randn(64,1), eloc, 'style', 'map');
%
%   See also: channelbrowse, channelbrowse_hist, channelcheck
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%************************************************************
%                      INPUT HANDLING
%************************************************************

p = inputParser;
addRequired(p,  'values',      @(x) validateattributes(x, {'numeric'}, {}));
addRequired(p,  'eloc',        @(x) validateattributes(x, {'struct'}, {'nonempty'}));
addParameter(p, 'style',      '',      @(x) any(validatestring(x, {'map','blank'})) || isempty(x));
addParameter(p, 'electrodes', 'pts',   @(x) any(validatestring(x, {'off','on','pts','labels','ptslabels'})));
addParameter(p, 'maplimits',  'absmax',@(x) (ischar(x) && any(validatestring(x, {'absmax','maxmin'}))) || ...
                                            (isnumeric(x) && numel(x) == 2));
addParameter(p, 'headrad',    0.5,     @(x) validateattributes(x, {'numeric'}, {'scalar','nonnegative','<=',1}));
% The map mask is applied per grid cell, so gridres sets how smooth the masked edge looks.
% EEGLAB defaulted to 67 (and topoplotFast to 32) for speed on 1990s hardware; a 200 grid
% costs well under a tenth of a second here and removes the visible staircase.
addParameter(p, 'gridres',    200,     @(x) validateattributes(x, {'numeric'}, {'scalar','integer','>=',32}));
addParameter(p, 'shading',    'interp',@(x) any(validatestring(x, {'interp','flat'})));
addParameter(p, 'conv',       'off',   @(x) any(validatestring(x, {'on','off'})));
addParameter(p, 'emarker',    {'.','k',[],1}, @(x) validateattributes(x, {'cell'}, {'vector'}));
addParameter(p, 'headcolor',  [0 0 0], @(x) validateattributes(x, {'numeric'}, {'vector','numel',3,'>=',0,'<=',1}));
addParameter(p, 'parent',     [],      @(x) isempty(x) || isgraphics(x, 'axes'));
parse(p, values, eloc, varargin{:});

values     = p.Results.values(:);
headrad    = p.Results.headrad;
gridres    = p.Results.gridres;
elec_mode  = lower(p.Results.electrodes);
headcolor  = p.Results.headcolor(:)';
use_hull   = strcmpi(p.Results.conv, 'on');
do_interp  = strcmpi(p.Results.shading, 'interp');

% 'on' is accepted as a synonym for 'pts' so that EEGLAB-style callsites port over verbatim
if strcmp(elec_mode, 'on')
    elec_mode = 'pts';
end

% Style defaults to whichever mode the caller implied by passing (or not passing) data
style = p.Results.style;
if isempty(style)
    if isempty(values)
        style = 'blank';
    else
        style = 'map';
    end
end
draw_map = strcmpi(style, 'map') && ~isempty(values);

hax = p.Results.parent;
if isempty(hax)
    hax = gca;
end

%************************************************************
%                   ELECTRODE PROJECTION
%************************************************************

[ex, ey, labels] = eloc_to_screen(eloc);
num_chans = numel(ex);

assert(isempty(values) || numel(values) == num_chans, ...
    'plot_topo:sizeMismatch', ...
    'values must have one entry per channel (%d channels, %d values).', num_chans, numel(values));

% The head rim sits at 0.5 by anatomical convention; this is the model's fixed scale, not a tunable
rmax = 0.5;

% Electrode coordinates are used exactly as projected -- deliberately NOT rescaled to
% pull the outermost ring onto the rim. Radius 0.5 already IS the rim (the ears/eyes
% equator), so Fp1 at 0.5 must land ON the circle. Scaling to fit the outermost
% electrode would drag every electrode inward and open a gap between the frontal
% ring and the head outline, which is wrong: it is the below-equator electrodes
% (FT9/TP9 on a 64 cap, at radius ~0.63) that belong OUTSIDE the head, in the
% conventional "skirt". The axes below are sized to include that skirt.
max_rad = max(hypot(ex, ey));

%************************************************************
%                    INTERPOLATED MAP
%************************************************************

set(hax, 'NextPlot', 'add');

if draw_map
    % NaN/Inf channels cannot participate in the interpolation, but they are still real
    % electrodes and are drawn as markers below.
    good = isfinite(values);
    assert(sum(good) >= 3, 'plot_topo:tooFewValues', ...
        'Need at least 3 finite values to interpolate a scalp map (%d given).', sum(good));

    gx = linspace(-rmax, rmax, gridres);
    gy = linspace(-rmax, rmax, gridres);
    [Xi, Yi] = meshgrid(gx, gy);

    % 'v4' is the biharmonic spline used by EEGLAB. It extrapolates smoothly past the
    % electrode ring, which is what makes the map fill to the head rim rather than
    % leaving a NaN halo the way a linear Delaunay fit would.
    Zi = griddata(double(ex(good)), double(ey(good)), double(values(good)), Xi, Yi, 'v4');

    % Never paint past the drawn head. The grid is a square, so its corners fall
    % outside the circle, and on a cap whose hull reaches the skirt electrodes
    % (FT9/TP9/PO9 at radius ~0.63) a hull-only mask would spill colour beyond the
    % outline. Clip to the head first, then narrow further if asked.
    outside = hypot(Xi, Yi) > rmax;

    if use_hull
        % Also mask to the electrode convex hull: past it the spline is pure
        % extrapolation, so showing it would imply coverage the montage lacks.
        hull = convhull(double(ex), double(ey));
        outside = outside | ~inpolygon(Xi, Yi, ex(hull), ey(hull));
    end
    Zi(outside) = nan;

    if do_interp
        face = 'interp';
    else
        face = 'flat';
    end
    surface(Xi, Yi, zeros(size(Zi)), Zi, ...
        'EdgeColor', 'none', 'FaceColor', face, 'Parent', hax);

    set(hax, 'CLim', resolve_maplimits(p.Results.maplimits, values(good)));
end

%************************************************************
%                        HEAD MODEL
%************************************************************

if headrad > 0
    draw_head(hax, headrad, headcolor);
end

%************************************************************
%                       ELECTRODES
%************************************************************

[marker, ecolor, msize, mwidth] = unpack_emarker(p.Results.emarker, num_chans);

% Electrodes and labels are lifted above the map surface so they are never z-hidden by it
if ~strcmp(elec_mode, 'off')
    if any(strcmp(elec_mode, {'pts','ptslabels'}))
        plot3(ex, ey, 2 * ones(num_chans, 1), marker, ...
            'Color', ecolor, 'MarkerSize', msize, 'LineWidth', mwidth, ...
            'LineStyle', 'none', 'Parent', hax);
    end

    if any(strcmp(elec_mode, {'labels','ptslabels'}))
        for ii = 1:num_chans
            text(ex(ii), ey(ii), 2.1, ['  ' labels{ii}], ...
                'Parent', hax, 'FontSize', 9, 'Color', ecolor, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
    end
end

%************************************************************
%                       AXES SETUP
%************************************************************

% 1.35x leaves room for the nose and ears; max_rad keeps any below-equator skirt
% electrodes (radius > 0.5) inside the axes rather than clipped at the edge
lim = 1.35 * max([headrad, rmax, max_rad]);
set(hax, 'XLim', [-lim lim], 'YLim', [-lim lim], ...
    'DataAspectRatio', [1 1 1], 'XTick', [], 'YTick', [], ...
    'Visible', 'off', 'NextPlot', 'replacechildren');
view(hax, 2);

if nargout == 0
    clear hax
end

end

%************************************************************
%                 PROJECT ELOC ONTO SCREEN
%************************************************************
function [ex, ey, labels] = eloc_to_screen(eloc)
%ELOC_TO_SCREEN  Convert an eloc struct to screen coordinates and label strings
%
%   Inputs:
%       eloc : 1xC struct - EEGLAB-style channel locations -- required
%
%   Outputs:
%       ex     : Cx1 double - screen x (positive toward the right ear)
%       ey     : Cx1 double - screen y (positive toward the nose)
%       labels : Cx1 cell - per-channel label strings

num_chans = numel(eloc);

has_polar = all(isfield(eloc, {'theta','radius'})) && ...
    ~any(cellfun(@isempty, {eloc.theta})) && ~any(cellfun(@isempty, {eloc.radius}));
has_cart = all(isfield(eloc, {'X','Y','Z'})) && ...
    ~any(cellfun(@isempty, {eloc.X}));

if has_polar
    theta = [eloc.theta];
    radius = [eloc.radius];
elseif has_cart
    % Fall back to the cartesian triple, matching EEGLAB's convertlocs cart2topo:
    % head X points at the nose, Y at the left ear, Z at the vertex. Radius is the
    % normalized polar angle from the vertex, so the equator (elevation 0) lands at 0.5.
    X = [eloc.X];
    Y = [eloc.Y];
    Z = [eloc.Z];
    theta = -atan2(Y, X) * 180/pi;
    radius = 0.5 - atan2(Z, hypot(X, Y)) / pi;
else
    error('plot_topo:noCoords', ...
        'eloc must supply either (theta, radius) or (X, Y, Z) coordinates.');
end

assert(numel(theta) == num_chans && numel(radius) == num_chans, ...
    'plot_topo:badCoords', ...
    'Every channel must have a coordinate; found %d theta and %d radius for %d channels.', ...
    numel(theta), numel(radius), num_chans);

% EEGLAB convention: theta is measured clockwise from the nose, so sin drives screen x
ex = radius(:) .* sin(theta(:) * pi/180);
ey = radius(:) .* cos(theta(:) * pi/180);

% Labels may be char, numeric, or absent depending on how the montage was built
labels = cell(num_chans, 1);
for ii = 1:num_chans
    if isfield(eloc, 'labels') && ~isempty(eloc(ii).labels)
        lab = eloc(ii).labels;
        if ischar(lab)
            labels{ii} = lab;
        elseif isstring(lab)
            labels{ii} = char(lab);
        else
            labels{ii} = num2str(lab);
        end
    else
        labels{ii} = num2str(ii);
    end
end

end

%************************************************************
%                     DRAW THE HEAD MODEL
%************************************************************
function draw_head(hax, headrad, headcolor)
%DRAW_HEAD  Draw the cartoon head outline, nose, and ears
%
%   Inputs:
%       hax       : axes handle - axes to draw into -- required
%       headrad   : double - drawn head radius -- required
%       headcolor : 1x3 double - line color -- required
%
%   Outputs:
%       none (side effects only)

linewidth = 1.7;
zh = 1.5;  % above the map surface (z=0), below the electrode markers (z=2)

% --- Head outline ---------------------------------------------------
t = linspace(0, 2*pi, 200);
hx = headrad * cos(t);
hy = headrad * sin(t);
plot3(hx, hy, zh * ones(size(hx)), '-', ...
    'Color', headcolor, 'LineWidth', linewidth, 'Parent', hax);

% --- Nose -----------------------------------------------------------
% A triangle seated on the rim at +y. Its base half-angle is chosen so the nose meets
% the circle exactly, which keeps the joint clean at any headrad.
nose_halfwidth = 0.18 * headrad;
nose_len = 0.18 * headrad;
base_ang = asin(nose_halfwidth / headrad);
nose_x = [-nose_halfwidth, 0, nose_halfwidth];
nose_y = [headrad*cos(base_ang), headrad + nose_len, headrad*cos(base_ang)];
plot3(nose_x, nose_y, zh * ones(size(nose_x)), '-', ...
    'Color', headcolor, 'LineWidth', linewidth, 'Parent', hax);

% --- Ears -----------------------------------------------------------
% Each ear is the outer half of an ellipse hung off the side of the head, so it reads as
% an ear without depending on EEGLAB's hard-coded outline.
ear_a = 0.09 * headrad;   % how far the ear sticks out
ear_b = 0.22 * headrad;   % ear height
et = linspace(-pi/2, pi/2, 60);
ear_x = headrad + ear_a * cos(et);
ear_y = ear_b * sin(et);

% Close each ear back onto the head rim so it does not float free of the outline
seat = asin(max(min(ear_y([1 end]) / headrad, 1), -1));
ear_x = [headrad*cos(seat(1)), ear_x, headrad*cos(seat(2))];
ear_y = [ear_y(1), ear_y, ear_y(end)];

plot3(ear_x, ear_y, zh * ones(size(ear_x)), '-', ...
    'Color', headcolor, 'LineWidth', linewidth, 'Parent', hax);   % right ear
plot3(-ear_x, ear_y, zh * ones(size(ear_x)), '-', ...
    'Color', headcolor, 'LineWidth', linewidth, 'Parent', hax);   % left ear

end

%************************************************************
%                    UNPACK EMARKER SPEC
%************************************************************
function [marker, ecolor, msize, mwidth] = unpack_emarker(spec, num_chans)
%UNPACK_EMARKER  Resolve the {marker color size linewidth} cell into plot arguments
%
%   Inputs:
%       spec      : cell - {marker, color, size, linewidth}, any element may be []
%                   -- required
%       num_chans : integer - channel count, used to scale the default marker size
%                   -- required
%
%   Outputs:
%       marker : char - marker symbol
%       ecolor : char or 1x3 double - marker color
%       msize  : double - marker size in points
%       mwidth : double - marker line width

defaults = {'.', 'k', [], 1};
spec = [spec(:)', defaults(numel(spec)+1:end)];

marker = spec{1};
ecolor = spec{2};
msize  = spec{3};
mwidth = spec{4};

if isempty(marker), marker = defaults{1}; end
if isempty(ecolor), ecolor = defaults{2}; end
if isempty(mwidth), mwidth = defaults{4}; end

% Dense montages need smaller dots to stay legible; scale the default the way EEGLAB does
if isempty(msize)
    if num_chans > 100
        msize = 3;
    elseif num_chans > 64
        msize = 4;
    else
        msize = 6;
    end
    % '.' renders visually smaller than the outline markers at the same point size
    if strcmp(marker, '.')
        msize = msize * 2;
    end
end

end

%************************************************************
%                    RESOLVE COLOR LIMITS
%************************************************************
function clim = resolve_maplimits(maplimits, values)
%RESOLVE_MAPLIMITS  Turn a maplimits option into a valid [lo hi] color limit pair
%
%   Inputs:
%       maplimits : char or 1x2 double - 'absmax', 'maxmin', or explicit limits
%                   -- required
%       values    : Nx1 double - finite channel values -- required
%
%   Outputs:
%       clim : 1x2 double - strictly increasing color limits

if ischar(maplimits)
    switch lower(maplimits)
        case 'absmax'
            % Symmetric about zero so that diverging colormaps put zero at the midpoint
            amax = max(abs(values));
            clim = [-amax amax];
        case 'maxmin'
            clim = [min(values) max(values)];
    end
else
    clim = sort(maplimits(:))';
end

% A flat map (all values equal, or a single finite channel) yields a degenerate range,
% which set(hax,'CLim',...) rejects. Widen it rather than erroring on valid data.
if ~(clim(2) > clim(1))
    clim = clim(1) + [-1 1] * max(abs(clim(1)) * 1e-3, eps);
end

end
