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
%       'interp'     : char - interpolant for the scalp map. 'cubic' is the
%                      Clough-Tocher scheme MNE-Python uses, ported from scipy
%                      rather than taken from griddata, which implements a
%                      different cubic. 'linear' and 'nearest' are MNE's other
%                      two modes, and 'v4' is the biharmonic spline that EEGLAB
%                      and FieldTrip use (default: 'cubic')
%       'extrapolate': char - how the map is carried past the electrodes,
%                      following MNE. 'head' rings the map with synthetic
%                      points and paints out to the clip circle, 'local'
%                      pushes the electrode convex hull out by one electrode
%                      spacing and masks to it, 'box' places four far corner
%                      points (default: 'head')
%       'border'     : char or double - value given to the synthetic points.
%                      'mean' gives each the average of its neighbouring
%                      electrodes; a number pins them all to that value
%                      (default: 'mean')
%       'cliprad'    : double - hard radius cap on the painted map, applied on
%                      top of whatever mask 'extrapolate' produces. This is a
%                      deliberate departure from MNE, which has no such mode:
%                      it can only narrow the map, never widen it, and it does
%                      not touch the interpolation, so values inside the cap
%                      are unchanged. Pass 0.5 to confine the map to the head
%                      rim the way this function did before it was matched to
%                      MNE. [] leaves MNE's behaviour alone (default: [])
%       'shading'    : char - 'interp' or 'flat' (default: 'interp')
%       'conv'       : char - kept for backwards compatibility. 'on' is a
%                      synonym for 'extrapolate','local' (default: 'off')
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
%       The scalp map reproduces MNE-Python's plot_topomap: the same azimuthal
%       equidistant projection, the same Clough-Tocher cubic interpolant, the
%       same synthetic boundary ring carrying a 'mean' border condition, and
%       the same clip circle, which is expanded past the head rim to just
%       beyond the outermost electrode rather than cutting the map at the rim.
%       The head, nose, and ear outlines use MNE's geometry as well.
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
%       Coordinates are plotted at their projected radius with no rescaling,
%       which is what MNE does too -- it keeps the head outline at the sphere
%       radius and widens the clip circle instead. (EEGLAB reaches the same
%       geometry the other way round, shrinking electrodes and head cartoon
%       together by rmax/plotrad.)
%       Electrodes below the equator have radius > 0.5 and therefore fall
%       outside the drawn head, in the conventional "skirt" -- FT9/FT10 and
%       TP9/TP10 on a 64-channel cap sit at ~0.63 and are genuinely below the
%       ears/eyes line. The axes are sized to keep them visible.
%
%       Channels whose values are NaN or Inf are dropped from the map
%       interpolation but still drawn as electrodes.
%
%   Example:
%       eloc = cap_montage('AC-64');
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
% It only samples the interpolant more finely -- it does not change the underlying surface.
% MNE defaults to 64 and EEGLAB to 67 (topoplotFast to 32) for speed on 1990s hardware;
% a 200 grid costs well under a tenth of a second here and removes the visible staircase.
addParameter(p, 'gridres',    200,     @(x) validateattributes(x, {'numeric'}, {'scalar','integer','>=',32}));
addParameter(p, 'interp',     'cubic', @(x) any(validatestring(x, {'cubic','linear','nearest','v4'})));
addParameter(p, 'extrapolate','',      @(x) isempty(x) || any(validatestring(x, {'head','local','box'})));
addParameter(p, 'border',     'mean',  @(x) (ischar(x) && strcmpi(x, 'mean')) || ...
                                            (isnumeric(x) && isscalar(x) && isfinite(x)));
% Genuine "empty or a positive scalar" disjunction, which attributes cannot express,
% so this is a predicate rather than a validateattributes call -- same as 'parent' below.
addParameter(p, 'cliprad',    [],      @(x) isempty(x) || ...
                                            (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
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
do_interp  = strcmpi(p.Results.shading, 'interp');
interp_fcn = lower(p.Results.interp);
cliprad    = p.Results.cliprad;
border     = p.Results.border;
if ischar(border)
    border = lower(border);
end

% 'conv','on' predates the MNE-shaped interface and means the same thing as
% 'extrapolate','local', so it is folded in here rather than kept as a second code path.
extrapolate = lower(p.Results.extrapolate);
if isempty(extrapolate)
    if strcmpi(p.Results.conv, 'on')
        extrapolate = 'local';
    else
        extrapolate = 'head';
    end
elseif strcmpi(p.Results.conv, 'on') && ~strcmp(extrapolate, 'local')
    error('plot_topo:conflictingMask', ...
        '''conv'',''on'' means ''extrapolate'',''local''; it conflicts with ''extrapolate'',''%s''.', ...
        extrapolate);
end

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

    [Xi, Yi, Zi] = mne_scalp_map(ex, ey, values, good, rmax, gridres, ...
        extrapolate, border, interp_fcn, cliprad);

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
%                    MNE-STYLE SCALP MAP
%************************************************************
function [Xi, Yi, Zi] = mne_scalp_map(ex, ey, values, good, rmax, gridres, extrapolate, border, method, cliprad)
%MNE_SCALP_MAP  Interpolate a scalp map the way MNE-Python's plot_topomap does
%
%   Inputs:
%       ex          : Cx1 double - electrode screen x, all channels -- required
%       ey          : Cx1 double - electrode screen y, all channels -- required
%       values      : Cx1 double - one value per channel -- required
%       good        : Cx1 logical - channels whose value is finite -- required
%       rmax        : double - head rim radius -- required
%       gridres     : integer - side length of the output grid -- required
%       extrapolate : char - 'head', 'local', or 'box' -- required
%       border      : char or double - 'mean', or a fixed value -- required
%       method      : char - griddata method: 'cubic', 'linear', 'nearest', 'v4'
%                     -- required
%       cliprad     : double or empty - extra hard radius cap on the mask, or []
%                     for MNE's behaviour -- required (may be empty)
%
%   Outputs:
%       Xi : gridres x gridres double - grid x coordinates
%       Yi : gridres x gridres double - grid y coordinates
%       Zi : gridres x gridres double - interpolated values, NaN outside the mask
%
%   Notes:
%       Mirrors MNE's _setup_interp / _GridData / _get_extra_points. MNE never lets the
%       interpolant extrapolate: it rings the electrodes with synthetic points carrying
%       a boundary condition and interpolates inside that ring. 'v4' is the exception,
%       since the biharmonic spline extrapolates by construction, which is exactly how
%       EEGLAB and FieldTrip use it -- so no ring is added for 'v4'.
%
%       The clip circle is widened past the head rim to just beyond the outermost
%       electrode (MNE's mask_scale), so below-equator channels get painted instead of
%       cut off. Every channel counts toward that radius, including ones dropped from
%       the fit for a non-finite value, so all drawn electrodes sit inside the map.
%
%       cliprad, when set, masks the finished grid to that radius. It is applied after
%       the interpolation and does not move the grid or the boundary ring, so the map
%       inside the cap is bit-identical to the uncapped one -- only less of it shows.

ex     = double(ex(:));
ey     = double(ey(:));
values = double(values(:));

% MNE: mask_scale = max(1, max||pos|| * 1.01 / radius). The 1.01 is MNE's; EEGLAB
% independently settled on 1.02 for the same job.
mask_scale  = max(1, max(hypot(ex, ey)) * 1.01 / rmax);
clip_radius = rmax * mask_scale;

pos = [ex(good), ey(good)];
v   = values(good);

% Triangulation-based interpolants collapse on coincident points, and every sleep cap
% stacks its EOG/EMG/ECG channels on a single placeholder position -- see cap_montage.
assert(size(unique(pos, 'rows'), 1) == size(pos, 1), 'plot_topo:duplicatePositions', ...
    ['Two or more channels with finite values share the same position. Exclude the ' ...
     'non-scalp channels (EOG/EMG/ECG) before interpolating a scalp map.']);

% --- Synthetic boundary points --------------------------------------
[extra, mask_xy] = boundary_points(pos, extrapolate, clip_radius, method);

if isempty(extra)
    fit_pos = pos;
    fit_val = v;
else
    fit_pos = [pos; extra];
    fit_val = [v; border_values(pos, extra, v, border)];
end

% --- Interpolate ----------------------------------------------------
gx = linspace(-clip_radius, clip_radius, gridres);
[Xi, Yi] = meshgrid(gx, gx);
if strcmp(method, 'cubic')
    % Deliberately not griddata's own 'cubic': that is a different triangulation-based
    % cubic scheme and lands several percent of the data range away from MNE. See
    % clough_tocher_2d.
    Zi = clough_tocher_2d(fit_pos, fit_val, Xi, Yi);
else
    Zi = griddata(fit_pos(:,1), fit_pos(:,2), fit_val, Xi, Yi, method);
end

% --- Mask -----------------------------------------------------------
% 'local' clips to its own polygon alone -- MNE does not also intersect it with the
% circle, so a montage that reaches past the rim keeps its full hull.
if strcmp(extrapolate, 'local') && ~isempty(mask_xy)
    outside = ~inpolygon(Xi, Yi, mask_xy(:,1), mask_xy(:,2));
else
    outside = hypot(Xi, Yi) > clip_radius;
end

% Applied last and to the mask only, so it never perturbs the interpolation. The grid
% still spans MNE's clip radius, so a capped map is a strict subset of the uncapped one.
if ~isempty(cliprad)
    outside = outside | (hypot(Xi, Yi) > cliprad);
end

Zi(outside) = nan;

end

%************************************************************
%                 SYNTHETIC BOUNDARY POINTS
%************************************************************
function [extra, mask_xy] = boundary_points(pos, extrapolate, clip_radius, method)
%BOUNDARY_POINTS  Place the synthetic points that carry MNE's boundary condition
%
%   Inputs:
%       pos         : Nx2 double - electrode positions used in the fit -- required
%       extrapolate : char - 'head', 'local', or 'box' -- required
%       clip_radius : double - radius of the painted disc -- required
%       method      : char - griddata method, used only to detect 'v4' -- required
%
%   Outputs:
%       extra   : Mx2 double - synthetic point positions, empty for 'v4'
%       mask_xy : Kx2 double - clip polygon for 'local', empty otherwise

extra   = zeros(0, 2);
mask_xy = [];

% The biharmonic spline extrapolates on its own -- that is precisely why EEGLAB and
% FieldTrip use it bare. Ringing it with synthetic points would change the surface
% rather than stabilize it, so 'v4' opts out and stays EEGLAB-comparable.
skip_ring = strcmp(method, 'v4');

switch extrapolate
    case 'head'
        if skip_ring
            return
        end
        distance = median_electrode_spacing(pos);
        % MNE spaces the ring so successive points subtend roughly one electrode
        % spacing, and pushes it out to 1.1*clip_radius + distance. That is outside the
        % painted disc, so every visible pixel is interpolated, never extrapolated.
        ang    = asin(min(distance / clip_radius, 1));
        n_pnts = max(12, round(2*pi / ang));
        t      = (0:n_pnts-1)' * (2*pi / n_pnts);
        r_ring = clip_radius * 1.1 + distance;
        extra  = [cos(t) * r_ring, sin(t) * r_ring];

    case 'local'
        distance = median_electrode_spacing(pos);
        [hull_extra, mask_xy] = local_hull_points(pos, distance);
        if ~skip_ring
            extra = hull_extra;
        end

    case 'box'
        if skip_ring
            return
        end
        % Four corners of a box three times the data range on a side
        lo = min(pos, [], 1);
        hi = max(pos, [], 1);
        d  = hi - lo;
        lo = lo - d;
        hi = hi + d;
        extra = [lo; lo(1) hi(2); hi(1) lo(2); hi];
end

end

%************************************************************
%                VALUES AT THE BOUNDARY POINTS
%************************************************************
function v_extra = border_values(pos, extra, v, border)
%BORDER_VALUES  Value assigned to each synthetic boundary point
%
%   Inputs:
%       pos    : Nx2 double - real electrode positions -- required
%       extra  : Mx2 double - synthetic point positions -- required
%       v      : Nx1 double - values at the real electrodes -- required
%       border : char or double - 'mean', or a fixed value -- required
%
%   Outputs:
%       v_extra : Mx1 double - value for each synthetic point

n_real  = size(pos, 1);
n_extra = size(extra, 1);

if ~ischar(border)
    v_extra = repmat(double(border), n_extra, 1);
    return
end

% border == 'mean': each synthetic point takes the average of the real electrodes it is
% adjacent to in the Delaunay graph of the combined point set.
DT = delaunayTriangulation([pos; extra]);
assert(size(DT.Points, 1) == n_real + n_extra, 'plot_topo:degenerateTriangulation', ...
    'The electrode layout and its boundary ring could not be triangulated cleanly.');

E       = edges(DT);
v_extra = zeros(n_extra, 1);
used    = false(n_extra, 1);
for ii = 1:n_extra
    jj  = n_real + ii;
    ngb = [E(E(:,1) == jj, 2); E(E(:,2) == jj, 1)];
    ngb = ngb(ngb <= n_real);
    if ~isempty(ngb)
        used(ii)    = true;
        v_extra(ii) = mean(v(ngb));
    end
end

% A ring point with no real neighbour falls back to the mean of the ones that had some,
% which is what MNE does rather than leaving it at zero.
if any(used) && ~all(used)
    v_extra(~used) = mean(v_extra(used));
end

end

%************************************************************
%               MEDIAN INTER-ELECTRODE SPACING
%************************************************************
function d = median_electrode_spacing(pos)
%MEDIAN_ELECTRODE_SPACING  Median edge length of the electrode Delaunay graph
%
%   Inputs:
%       pos : Nx2 double - electrode positions -- required
%
%   Outputs:
%       d : double - median inter-electrode distance, MNE's scale for the ring

if size(pos, 1) >= 4
    try
        DT = delaunayTriangulation(pos);
        T  = DT.ConnectivityList;
        if ~isempty(T)
            % Two of the three edges per triangle, matching MNE's _get_extra_points.
            % Which two is arbitrary in both libraries -- the median barely moves.
            e1 = sqrt(sum((pos(T(:,1),:) - pos(T(:,2),:)).^2, 2));
            e2 = sqrt(sum((pos(T(:,2),:) - pos(T(:,3),:)).^2, 2));
            d  = median([e1; e2]);
            return
        end
    catch
        % Collinear layouts make delaunayTriangulation throw; fall through.
    end
end

% Collinear, or too few electrodes to triangulate: median nearest-neighbour spacing
D = hypot(pos(:,1) - pos(:,1)', pos(:,2) - pos(:,2)');
D(1:size(D,1)+1:end) = inf;
d = median(min(D, [], 2));

end

%************************************************************
%             LOCAL (CONVEX HULL) BOUNDARY POINTS
%************************************************************
function [extra, mask_xy] = local_hull_points(pos, distance)
%LOCAL_HULL_POINTS  Push the electrode convex hull outward, MNE's 'local' extrapolation
%
%   Inputs:
%       pos      : Nx2 double - electrode positions -- required
%       distance : double - median inter-electrode spacing -- required
%
%   Outputs:
%       extra   : Mx2 double - synthetic point positions, one spacing outside the hull
%       mask_xy : Kx2 double - clip polygon, half a spacing outside the hull

DT   = delaunayTriangulation(pos);
loop = convexHull(DT);            % closed loop of point indices
p1   = pos(loop(1:end-1), :);     % one row per hull edge
p2   = pos(loop(2:end), :);

% Each hull vertex is pushed radially away from the electrode centroid: the fit points by
% a full electrode spacing, the mask by half of one, so the mask always sits inside the
% ring of fit points and no painted pixel is extrapolated.
c  = mean(pos, 1);
u1 = (p1 - c) ./ sqrt(sum((p1 - c).^2, 2));
u2 = (p2 - c) ./ sqrt(sum((p2 - c).^2, 2));

e1 = p1 + u1 * distance;
e2 = p2 + u2 * distance;

mask_xy    = unique([p1 + u1 * distance * 0.5; p2 + u2 * distance * 0.5], 'rows');
mc         = mean(mask_xy, 1);
[~, order] = sort(atan2(mask_xy(:,2) - mc(2), mask_xy(:,1) - mc(1)));
mask_xy    = mask_xy(order, :);

% Long hull edges get intermediate points so the ring stays about one spacing apart
hull_diff = p2 - p1;
n_seg     = round(0.25 * sqrt(sum(hull_diff.^2, 2)) / distance);
add       = zeros(0, 2);
for nn = 2:max([n_seg; 0])
    sel = (n_seg == nn);
    if ~any(sel)
        continue
    end
    for kk = 1:(nn - 1)
        add = [add; e1(sel,:) + hull_diff(sel,:) * (kk / nn)]; %#ok<AGROW>
    end
end

extra = [unique([e1; e2], 'rows'); add];

end

%************************************************************
%              CLOUGH-TOCHER CUBIC INTERPOLANT
%************************************************************
function Zq = clough_tocher_2d(P, v, Xq, Yq)
%CLOUGH_TOCHER_2D  Piecewise cubic C1 interpolant, matching scipy's CloughTocher2DInterpolator
%
%   Inputs:
%       P  : Nx2 double - interpolation point coordinates -- required
%       v  : Nx1 double - values at those points -- required
%       Xq : MxK double - query x coordinates -- required
%       Yq : MxK double - query y coordinates -- required
%
%   Outputs:
%       Zq : MxK double - interpolated values, NaN outside the convex hull of P
%
%   Notes:
%       MATLAB's own griddata 'cubic' is also triangulation-based, but it is a
%       different cubic scheme: on a 64-channel cap it lands about 6.5% of the data
%       range (40% at the worst pixel) away from what MNE draws. This is a direct
%       port of the interpolant scipy actually uses -- a cubic Bezier patch per
%       triangle under the Clough-Tocher split, with vertex gradients from the
%       approximate-curvature-minimization network of Nielson (1983) and Renka &
%       Cline (1984), and the affine-invariant choice of cross-boundary direction
%       that scipy's interpnd.pyx documents.
%
%       On a generic point set this reproduces scipy to machine precision (~4e-15).
%       On a real montage the synthetic boundary ring is exactly cocircular, which is
%       a degenerate Delaunay configuration that MATLAB and Qhull resolve differently;
%       the resulting disagreement is confined to the annulus outside the outermost
%       electrodes and stays under about 1.2% of the data range anywhere a real
%       electrode contributes.

P = double(P);
v = double(v(:));

DT = delaunayTriangulation(P);
T  = DT.ConnectivityList;
NB = neighbors(DT);            % NB(t,k) is the triangle opposite vertex k, NaN at a hull edge

grad = estimate_gradients_2d(DT, P, v);

% --- Bezier control net, one triangle at a time ----------------------
% Vertex values and edge-projected gradients give the boundary control points; the
% interior ones follow from the C1 conditions. Names follow scipy's c<ijkl> convention.
i1 = T(:,1);
i2 = T(:,2);
i3 = T(:,3);
e12 = P(i2,:) - P(i1,:);
e23 = P(i3,:) - P(i2,:);
e31 = P(i1,:) - P(i3,:);

df12 =  sum(grad(i1,:) .* e12, 2);
df21 = -sum(grad(i2,:) .* e12, 2);
df23 =  sum(grad(i2,:) .* e23, 2);
df32 = -sum(grad(i3,:) .* e23, 2);
df31 =  sum(grad(i3,:) .* e31, 2);
df13 = -sum(grad(i1,:) .* e31, 2);

c3000 = v(i1);
c2100 = (df12 + 3*c3000)/3;
c2010 = (df13 + 3*c3000)/3;
c0300 = v(i2);
c1200 = (df21 + 3*c0300)/3;
c0210 = (df23 + 3*c0300)/3;
c0030 = v(i3);
c1020 = (df31 + 3*c0030)/3;
c0120 = (df32 + 3*c0030)/3;

c2001 = (c2100 + c2010 + c3000)/3;
c0201 = (c1200 + c0300 + c0210)/3;
c0021 = (c1020 + c0120 + c0030)/3;

% --- Cross-boundary direction weights --------------------------------
% Picking the edge normal here would make the interpolant non-affine-invariant and let
% sliver triangles blow up, so scipy instead points at the neighbouring triangle's
% centroid, expressed in barycentric coordinates. A hull edge has no neighbour and
% falls back to -1/2, the direction of this triangle's own centroid.
gw = zeros(size(T,1), 3);
for kk = 1:3
    nbt = NB(:,kk);
    has = ~isnan(nbt);
    gw(~has, kk) = -1/2;
    if any(has)
        tri_idx = find(has);
        cen = (P(T(nbt(tri_idx),1),:) + P(T(nbt(tri_idx),2),:) + P(T(nbt(tri_idx),3),:)) / 3;
        c = cartesianToBarycentric(DT, tri_idx, cen);
        switch kk
            case 1
                gw(tri_idx,kk) = (2*c(:,3) + c(:,2) - 1) ./ (2 - 3*c(:,3) - 3*c(:,2));
            case 2
                gw(tri_idx,kk) = (2*c(:,1) + c(:,3) - 1) ./ (2 - 3*c(:,1) - 3*c(:,3));
            case 3
                gw(tri_idx,kk) = (2*c(:,2) + c(:,1) - 1) ./ (2 - 3*c(:,2) - 3*c(:,1));
        end
    end
end

c0111 = (gw(:,1).*(-c0300 + 3*c0210 - 3*c0120 + c0030) + (-c0300 + 2*c0210 - c0120 + c0021 + c0201))/2;
c1011 = (gw(:,2).*(-c0030 + 3*c1020 - 3*c2010 + c3000) + (-c0030 + 2*c1020 - c2010 + c2001 + c0021))/2;
c1101 = (gw(:,3).*(-c3000 + 3*c2100 - 3*c1200 + c0300) + (-c3000 + 2*c2100 - c1200 + c2001 + c0201))/2;

c1002 = (c1101 + c1011 + c2001)/3;
c0102 = (c1101 + c0111 + c0201)/3;
c0012 = (c1011 + c0111 + c0021)/3;
c0003 = (c1002 + c0102 + c0012)/3;

% --- Evaluate --------------------------------------------------------
qxy = [Xq(:), Yq(:)];
Zq  = nan(size(qxy,1), 1);
tq  = pointLocation(DT, qxy);
in  = ~isnan(tq);          % pointLocation returns NaN outside the convex hull
if any(in)
    t = tq(in);
    b = cartesianToBarycentric(DT, t, qxy(in,:));

    % Extended (4-coordinate) barycentrics for the Clough-Tocher split: dropping the
    % smallest coordinate picks the sub-triangle, and b4 is the weight on the centroid.
    bmin = min(b, [], 2);
    b1 = b(:,1) - bmin;
    b2 = b(:,2) - bmin;
    b3 = b(:,3) - bmin;
    b4 = 3*bmin;

    % One of b1..b4 is zero by construction, so the terms scipy omits here are zero too
    Zq(in) = b1.^3.*c3000(t)        + 3*b1.^2.*b2.*c2100(t) + 3*b1.^2.*b3.*c2010(t) + ...
             3*b1.^2.*b4.*c2001(t)  + 3*b1.*b2.^2.*c1200(t) + 6*b1.*b2.*b4.*c1101(t) + ...
             3*b1.*b3.^2.*c1020(t)  + 6*b1.*b3.*b4.*c1011(t) + 3*b1.*b4.^2.*c1002(t) + ...
             b2.^3.*c0300(t)        + 3*b2.^2.*b3.*c0210(t) + 3*b2.^2.*b4.*c0201(t) + ...
             3*b2.*b3.^2.*c0120(t)  + 6*b2.*b3.*b4.*c0111(t) + 3*b2.*b4.^2.*c0102(t) + ...
             b3.^3.*c0030(t)        + 3*b3.^2.*b4.*c0021(t) + 3*b3.*b4.^2.*c0012(t) + ...
             b4.^3.*c0003(t);
end

Zq = reshape(Zq, size(Xq));

end

%************************************************************
%              VERTEX GRADIENTS BY CURVATURE MIN
%************************************************************
function grad = estimate_gradients_2d(DT, P, v)
%ESTIMATE_GRADIENTS_2D  Vertex gradients that approximately minimize surface curvature
%
%   Inputs:
%       DT : delaunayTriangulation - triangulation of P -- required
%       P  : Nx2 double - point coordinates -- required
%       v  : Nx1 double - values at those points -- required
%
%   Outputs:
%       grad : Nx2 double - [dF/dx, dF/dy] at each point
%
%   Notes:
%       Port of scipy's _estimate_gradients_2d_global. Restricted to one edge the
%       Clough-Tocher interpolant is a cubic in the edge parameter, so the bending
%       energy over that edge is a quadratic form in the two end gradients. Summing
%       over the edges at a vertex leaves a 2x2 system for that vertex's gradient,
%       and the whole network is relaxed by Gauss-Seidel -- each vertex is solved
%       against its neighbours' current values, in index order, until the largest
%       relative change falls below tol.
%
%       scipy's defaults (tol 1e-6, maxiter 400) are kept. A 64-channel montage plus
%       its boundary ring converges in roughly a dozen sweeps.

tol     = 1e-6;
maxiter = 400;

num_pts = size(P, 1);

% Adjacency in the Delaunay graph, which is what the edge sum below runs over
E = edges(DT);
adj = cell(num_pts, 1);
for ii = 1:size(E,1)
    adj{E(ii,1)}(end+1) = E(ii,2);
    adj{E(ii,2)}(end+1) = E(ii,1);
end

grad = zeros(num_pts, 2);

for iter_num = 1:maxiter
    err = 0;
    for ipt = 1:num_pts
        nbr = adj{ipt}(:);
        ev  = P(nbr,:) - P(ipt,:);          % edge vectors away from this vertex
        L3  = sum(ev.^2, 2).^1.5;

        % Neighbour gradients projected onto the edge, held fixed for this sweep
        df2 = -sum(ev .* grad(nbr,:), 2);

        Qm   = 4 * (ev' * (ev ./ L3));
        coef = (6*(v(ipt) - v(nbr)) - 2*df2) ./ L3;
        sv   = (coef' * ev)';

        detQ = Qm(1,1)*Qm(2,2) - Qm(1,2)*Qm(2,1);
        r    = [( Qm(2,2)*sv(1) - Qm(1,2)*sv(2)) / detQ; ...
                (-Qm(2,1)*sv(1) + Qm(1,1)*sv(2)) / detQ];

        change = max(abs(grad(ipt,1) + r(1)), abs(grad(ipt,2) + r(2)));
        grad(ipt,:) = -r';

        % Relative where the gradient is large, absolute where it is small
        change = change / max(1, max(abs(r(1)), abs(r(2))));
        err = max(err, change);
    end
    if err < tol
        break
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
t = linspace(0, 2*pi, 101);   % 101 points, as MNE's _make_head_outlines uses
hx = headrad * cos(t);
hy = headrad * sin(t);
plot3(hx, hy, zh * ones(size(hx)), '-', ...
    'Color', headcolor, 'LineWidth', linewidth, 'Parent', hax);

% --- Nose -----------------------------------------------------------
% MNE writes the nose seat as exp(1i*acos(deg2rad(12))), which is a roundabout way of
% getting the pair (0.2094, 0.9778): half width 0.2094*r, base at 0.9778*r, tip at
% 1.15*r. Reproduced as-is so the cartoon matches plot_topomap rather than approximating it.
nose_seat = exp(1i * acos(deg2rad(12)));
nose_dx = real(nose_seat);
nose_dy = imag(nose_seat);
nose_x = [-nose_dx, 0, nose_dx] * headrad;
nose_y = [nose_dy, 1.15, nose_dy] * headrad;
plot3(nose_x, nose_y, zh * ones(size(nose_x)), '-', ...
    'Color', headcolor, 'LineWidth', linewidth, 'Parent', hax);

% --- Ears -----------------------------------------------------------
% MNE's ear outline, which it in turn inherited from EEGLAB's EarX/EarY. The coefficients
% are fractions of the head diameter, hence the 2*headrad scaling. The polyline is open at
% both ends by design: its first and last points sit just inside the rim at 0.994*headrad.
ear_x = [0.497 0.510 0.518 0.5299 0.5419 0.54 0.547 0.532 0.510 0.489] * (2 * headrad);
ear_y = [0.0555 0.0775 0.0783 0.0746 0.0555 -0.0055 -0.0932 -0.1313 -0.1384 -0.1199] * (2 * headrad);

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
