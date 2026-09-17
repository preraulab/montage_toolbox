function [hax, ex, ey, hmap] = plot_topo(values, eloc, varargin)
%PLOT_TOPO  Plot electrode locations (and optional scalp map) on a head model
%
%   Usage:
%       [hax, ex, ey] = plot_topo([], eloc)
%       [hax, ex, ey] = plot_topo([], eloc, 'electrodes', 'ptslabels')
%       [hax, ex, ey, hmap] = plot_topo(values, eloc, 'style', 'map')
%       plot_topo(newvalues, [], 'update', hmap)   % fast redraw, same montage
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
%       'update'     : image handle - a map returned earlier as hmap. Repaints it
%                      from new values, reusing the cached geometry, and returns
%                      immediately. eloc may be [] and every other option is
%                      ignored, since they are fixed by the original call
%                      (default: [])
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
%       hax  : axes handle - the axes drawn into
%       ex   : Cx1 double - electrode screen x coordinates, in channel order
%       ey   : Cx1 double - electrode screen y coordinates, in channel order
%       hmap : image handle - the scalp map object, or [] when no map was drawn.
%              Pass it back as 'update' to repaint with new values; the
%              interpolation geometry is cached on it.
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
%   Credit:
%       The scalp map is a deliberate reimplementation of MNE-Python's
%       mne.viz.plot_topomap (mne/viz/topomap.py): its projection, synthetic
%       boundary ring and 'mean' border condition, clip radius, extrapolation
%       modes, and head/nose/ear geometry. The Clough-Tocher interpolant in
%       cache_clough_tocher and topo_evaluate is a port of SciPy's
%       CloughTocher2DInterpolator (scipy/interpolate/interpnd.pyx). Both
%       projects are BSD-3-Clause, as is this file; neither endorses this port.
%       Please cite Gramfort et al. (2013) and Virtanen et al. (2020) alongside
%       this toolbox. See README.md for the full reference list.
%
%   See also: cap_montage, read_montage, default_montage, xyz_to_eloc
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%************************************************************
%                      INPUT HANDLING
%************************************************************

p = inputParser;
addRequired(p,  'values',      @(x) validateattributes(x, {'numeric'}, {}));
% eloc may be empty when 'update' is given -- the geometry is already cached there.
% "empty or a nonempty struct" is a disjunction across types, which attributes cannot
% express, so this is a predicate rather than a validateattributes call.
addRequired(p,  'eloc',        @(x) isempty(x) || (isstruct(x) && ~isempty(x)));
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
addParameter(p, 'update',     [],      @(x) isempty(x) || isgraphics(x, 'image'));
addParameter(p, 'parent',     [],      @(x) isempty(x) || isgraphics(x, 'axes'));
parse(p, values, eloc, varargin{:});

%************************************************************
%                      FAST UPDATE PATH
%************************************************************

% Everything below this block is setup that an update does not need: the geometry is
% already cached on the map handle, and the head, markers and axes are already drawn.
hmap = p.Results.update;
if ~isempty(hmap)
    S = getappdata(hmap, 'plot_topo_cache');
    assert(~isempty(S), 'plot_topo:noCache', ...
        ['That image carries no interpolation cache. ''update'' takes the hmap output ' ...
         'of an earlier plot_topo call, not an arbitrary image handle.']);

    vals = double(p.Results.values(:));
    assert(numel(vals) == S.num_chans, 'plot_topo:sizeMismatch', ...
        'values must have one entry per channel (%d channels, %d values).', ...
        S.num_chans, numel(vals));
    assert(isequal(isfinite(vals), S.good), 'plot_topo:goodChanged', ...
        ['Which channels are finite has changed since the map was drawn, so the cached ' ...
         'triangulation no longer applies. Call plot_topo without ''update'' to rebuild.']);

    [Zi, S] = topo_evaluate(S, vals);
    setappdata(hmap, 'plot_topo_cache', S);
    set(hmap, 'CData', Zi);

    hax = ancestor(hmap, 'axes');
    ex  = S.ex;
    ey  = S.ey;
    if nargout == 0
        clear hax
    end
    return
end

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
hmap = [];

hax = p.Results.parent;
if isempty(hax)
    hax = gca;
end

%************************************************************
%                   ELECTRODE PROJECTION
%************************************************************

assert(~isempty(eloc), 'plot_topo:noEloc', ...
    'eloc is required unless ''update'' is given.');

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

    % Everything that does not depend on the values is built once and kept on the map
    % handle, so an 'update' call can skip straight to topo_evaluate.
    S = topo_cache(ex, ey, good, rmax, gridres, extrapolate, border, interp_fcn, cliprad);
    [Zi, S] = topo_evaluate(S, double(values));

    % An image, not a surface: the grid is regular, so this is the right primitive and it
    % renders about 3x faster. AlphaData carries the mask, since an image has no NaN
    % transparency of its own.
    hmap = image('XData', S.gx, 'YData', S.gx, 'CData', Zi, ...
        'CDataMapping', 'scaled', 'AlphaData', double(isfinite(Zi)), 'Parent', hax);

    % 'shading' maps onto the image's own interpolation. The property is recent, so fall
    % back to flat rather than erroring on older releases.
    if isprop(hmap, 'Interpolation')
        if do_interp
            set(hmap, 'Interpolation', 'bilinear');
        else
            set(hmap, 'Interpolation', 'nearest');
        end
    end

    setappdata(hmap, 'plot_topo_cache', S);

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

% Drawn after the map, so child order alone puts these on top; uistack below makes that
% explicit rather than relying on it. This is what MNE expresses with zorder.
if ~strcmp(elec_mode, 'off')
    if any(strcmp(elec_mode, {'pts','ptslabels'}))
        plot(ex, ey, marker, ...
            'Color', ecolor, 'MarkerSize', msize, 'LineWidth', mwidth, ...
            'LineStyle', 'none', 'Parent', hax);
    end

    if any(strcmp(elec_mode, {'labels','ptslabels'}))
        for ii = 1:num_chans
            text(ex(ii), ey(ii), ['  ' labels{ii}], ...
                'Parent', hax, 'FontSize', 9, 'Color', ecolor, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
    end
end

if ~isempty(hmap)
    uistack(hmap, 'bottom');
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
%                  CACHED MAP GEOMETRY
%************************************************************
function S = topo_cache(ex, ey, good, rmax, gridres, extrapolate, border, method, cliprad)
%TOPO_CACHE  Precompute everything about a scalp map that does not depend on the values
%
%   Inputs:
%       ex          : Cx1 double - electrode screen x, all channels -- required
%       ey          : Cx1 double - electrode screen y, all channels -- required
%       good        : Cx1 logical - channels whose value is finite -- required
%       rmax        : double - head rim radius -- required
%       gridres     : integer - side length of the output grid -- required
%       extrapolate : char - 'head', 'local', or 'box' -- required
%       border      : char or double - 'mean', or a fixed value -- required
%       method      : char - 'cubic', 'linear', 'nearest', or 'v4' -- required
%       cliprad     : double or empty - extra hard radius cap -- required (may be empty)
%
%   Outputs:
%       S : struct - geometry cache, consumed by topo_evaluate
%
%   Notes:
%       This is the expensive half of a redraw and none of it depends on the data, which
%       is why it is worth keeping. Splitting the interpolator this way is what MNE's
%       _GridData does: "computing parameters for a fixed set of true points, and
%       allowing the values at those points to be set independently".
%
%       The cubic branch is a port of SciPy's CloughTocher2DInterpolator
%       (scipy/interpolate/interpnd.pyx, BSD-3-Clause), split into this cache and
%       topo_evaluate. MATLAB's own griddata 'cubic' is a different triangulation-based
%       cubic scheme and does not reproduce MNE's output.
%
%       For 'cubic' the cache also holds the inverse of each vertex's 2x2
%       curvature-minimizing system. Q = 4*sum(e*e'/L^3) is built from edge vectors
%       alone, so it is pure geometry even though it lives inside the value solve.

S = struct();
S.ex        = double(ex(:));
S.ey        = double(ey(:));
S.good      = logical(good(:));
S.num_chans = numel(S.ex);
S.method    = method;
S.extrap    = extrapolate;
S.border    = border;

% MNE: mask_scale = max(1, max||pos|| * 1.01 / radius)
mask_scale    = max(1, max(hypot(S.ex, S.ey)) * 1.01 / rmax);
S.clip_radius = rmax * mask_scale;

pos = [S.ex(S.good), S.ey(S.good)];

assert(size(unique(pos, 'rows'), 1) == size(pos, 1), 'plot_topo:duplicatePositions', ...
    ['Two or more channels with finite values share the same position. Exclude the ' ...
     'non-scalp channels (EOG/EMG/ECG) before interpolating a scalp map.']);

[extra, mask_xy] = boundary_points(pos, extrapolate, S.clip_radius, method);
S.pos     = pos;
S.extra   = extra;
S.n_real  = size(pos, 1);
S.n_extra = size(extra, 1);
S.fit_pos = [pos; extra];

S.gx = linspace(-S.clip_radius, S.clip_radius, gridres);
[Xi, Yi] = meshgrid(S.gx, S.gx);
S.size = size(Xi);

% --- Mask ------------------------------------------------------------
if strcmp(extrapolate, 'local') && ~isempty(mask_xy)
    outside = ~inpolygon(Xi, Yi, mask_xy(:,1), mask_xy(:,2));
else
    outside = hypot(Xi, Yi) > S.clip_radius;
end
if ~isempty(cliprad)
    outside = outside | (hypot(Xi, Yi) > cliprad);
end

% The Delaunay graph is needed by the cubic interpolant, and also by border == 'mean'
% whichever interpolant is in use, so it is built once here for both.
need_tri = strcmp(method, 'cubic') || (S.n_extra > 0 && ischar(border));
S.ring_nbr = cell(S.n_extra, 1);
if need_tri
    DT = delaunayTriangulation(S.fit_pos);
    assert(size(DT.Points, 1) == size(S.fit_pos, 1), 'plot_topo:degenerateTriangulation', ...
        'The electrode layout and its boundary ring could not be triangulated cleanly.');
    adj = delaunay_adjacency(DT, size(S.fit_pos, 1));
    for kk = 1:S.n_extra
        nb = adj{S.n_real + kk}(:);
        S.ring_nbr{kk} = nb(nb <= S.n_real);
    end
end

if strcmp(method, 'cubic')
    S = cache_clough_tocher(S, DT, adj, Xi, Yi, outside);
else
    % griddata rebuilds its own triangulation each call, so only the grid is cacheable
    S.Xi = Xi;
    S.Yi = Yi;
    S.outside = outside;
end

end

%************************************************************
%              CACHED CLOUGH-TOCHER GEOMETRY
%************************************************************
function S = cache_clough_tocher(S, DT, adj, Xi, Yi, outside)
%CACHE_CLOUGH_TOCHER  Triangulation, point location and the per-vertex 2x2 systems
%
%   Inputs:
%       S       : struct - partially built cache -- required
%       DT      : delaunayTriangulation - of the electrodes plus boundary ring -- required
%       adj     : cell - Delaunay adjacency list, one entry per point -- required
%       Xi      : MxM double - grid x -- required
%       Yi      : MxM double - grid y -- required
%       outside : MxM logical - mask, true where nothing is painted -- required
%
%   Outputs:
%       S : struct - cache with the triangulation and evaluation tables added

P = S.fit_pos;
S.DT = DT;
T    = DT.ConnectivityList;
S.T  = T;
NB   = neighbors(DT);
n_pts = size(P, 1);

% Per-vertex geometry of the curvature-minimizing system. Q holds no values, so its
% inverse is cached here and the per-frame solve becomes a single 2x2 multiply.
S.nbr   = cell(n_pts, 1);
S.edge  = cell(n_pts, 1);
S.L3    = cell(n_pts, 1);
S.Qinv  = cell(n_pts, 1);
for ip = 1:n_pts
    nb = adj{ip}(:);
    ev = P(nb,:) - P(ip,:);
    L3 = sum(ev.^2, 2).^1.5;
    Q  = 4 * (ev' * (ev ./ L3));
    S.nbr{ip}  = nb;
    S.edge{ip} = ev;
    S.L3{ip}   = L3;
    S.Qinv{ip} = inv(Q);
end

% Edge vectors and the affine-invariant cross-boundary weights, both value-free
S.i1 = T(:,1);
S.i2 = T(:,2);
S.i3 = T(:,3);
S.e12 = P(S.i2,:) - P(S.i1,:);
S.e23 = P(S.i3,:) - P(S.i2,:);
S.e31 = P(S.i1,:) - P(S.i3,:);

S.gw = zeros(size(T,1), 3);
for kk = 1:3
    nbt = NB(:,kk);
    has = ~isnan(nbt);
    S.gw(~has, kk) = -1/2;
    if any(has)
        tri_idx = find(has);
        cen = (P(T(nbt(tri_idx),1),:) + P(T(nbt(tri_idx),2),:) + P(T(nbt(tri_idx),3),:)) / 3;
        c = cartesianToBarycentric(DT, tri_idx, cen);
        switch kk
            case 1
                S.gw(tri_idx,kk) = (2*c(:,3) + c(:,2) - 1) ./ (2 - 3*c(:,3) - 3*c(:,2));
            case 2
                S.gw(tri_idx,kk) = (2*c(:,1) + c(:,3) - 1) ./ (2 - 3*c(:,1) - 3*c(:,3));
            case 3
                S.gw(tri_idx,kk) = (2*c(:,2) + c(:,1) - 1) ./ (2 - 3*c(:,2) - 3*c(:,1));
        end
    end
end

% Point location and extended barycentrics: the single most expensive cacheable step
tq = pointLocation(DT, [Xi(:), Yi(:)]);
S.in = ~isnan(tq) & ~outside(:);
S.tq = tq(S.in);
b = cartesianToBarycentric(DT, S.tq, [Xi(S.in), Yi(S.in)]);
bmin = min(b, [], 2);
S.b1 = b(:,1) - bmin;
S.b2 = b(:,2) - bmin;
S.b3 = b(:,3) - bmin;
S.b4 = 3 * bmin;

S.grad = zeros(n_pts, 2);   % warm-start seed for the first solve

end

%************************************************************
%                   DELAUNAY ADJACENCY LIST
%************************************************************
function adj = delaunay_adjacency(DT, n_pts)
%DELAUNAY_ADJACENCY  Neighbour list for every point of a triangulation
%
%   Inputs:
%       DT    : delaunayTriangulation - the triangulation -- required
%       n_pts : integer - number of points in it -- required
%
%   Outputs:
%       adj : n_ptsx1 cell - indices of each point's Delaunay neighbours

E = edges(DT);
adj = cell(n_pts, 1);
for ii = 1:size(E, 1)
    adj{E(ii,1)}(end+1) = E(ii,2);
    adj{E(ii,2)}(end+1) = E(ii,1);
end

end

%************************************************************
%                EVALUATE A CACHED MAP
%************************************************************
function [Zi, S] = topo_evaluate(S, values)
%TOPO_EVALUATE  Interpolate one frame's values onto a cached grid
%
%   Inputs:
%       S      : struct - cache from topo_cache -- required
%       values : Cx1 double - one value per channel -- required
%
%   Outputs:
%       Zi : gridres x gridres double - interpolated values, NaN outside the mask
%       S  : struct - cache with the gradient solution carried forward
%
%   Notes:
%       The gradient sweep and the Bezier net are ported from SciPy's
%       _estimate_gradients_2d_global and _clough_tocher_2d_single (interpnd.pyx,
%       BSD-3-Clause), the curvature-minimizing network of Nielson (1983) and
%       Renka & Cline (1984).
%
%       S is returned so the vertex gradients survive to the next frame. They are the
%       initial guess for the next solve, which roughly halves the Gauss-Seidel cost on
%       smoothly varying data and changes nothing about the converged answer beyond the
%       1e-6 tolerance that stops it.

v = double(values(:));
v = v(S.good);

if ~strcmp(S.method, 'cubic')
    fit_val = fit_values(S, v);
    Zi = griddata(S.fit_pos(:,1), S.fit_pos(:,2), fit_val, S.Xi, S.Yi, S.method);
    Zi(S.outside) = nan;
    return
end

V = fit_values(S, v);

% --- Vertex gradients, Gauss-Seidel over the Delaunay graph ----------
g   = S.grad;
tol = 1e-6;
for iter_num = 1:400
    err = 0;
    for ip = 1:numel(V)
        nb = S.nbr{ip};
        ev = S.edge{ip};
        df2 = -sum(ev .* g(nb,:), 2);
        sv  = (((6*(V(ip) - V(nb)) - 2*df2) ./ S.L3{ip})' * ev)';
        r   = S.Qinv{ip} * sv;
        change = max(abs(g(ip,1) + r(1)), abs(g(ip,2) + r(2)));
        g(ip,:) = -r';
        err = max(err, change / max(1, max(abs(r))));
    end
    if err < tol
        break
    end
end
S.grad = g;

% --- Bezier control net ----------------------------------------------
df12 =  sum(g(S.i1,:) .* S.e12, 2);
df21 = -sum(g(S.i2,:) .* S.e12, 2);
df23 =  sum(g(S.i2,:) .* S.e23, 2);
df32 = -sum(g(S.i3,:) .* S.e23, 2);
df31 =  sum(g(S.i3,:) .* S.e31, 2);
df13 = -sum(g(S.i1,:) .* S.e31, 2);

c3000 = V(S.i1);
c2100 = (df12 + 3*c3000)/3;
c2010 = (df13 + 3*c3000)/3;
c0300 = V(S.i2);
c1200 = (df21 + 3*c0300)/3;
c0210 = (df23 + 3*c0300)/3;
c0030 = V(S.i3);
c1020 = (df31 + 3*c0030)/3;
c0120 = (df32 + 3*c0030)/3;
c2001 = (c2100 + c2010 + c3000)/3;
c0201 = (c1200 + c0300 + c0210)/3;
c0021 = (c1020 + c0120 + c0030)/3;

c0111 = (S.gw(:,1).*(-c0300 + 3*c0210 - 3*c0120 + c0030) + (-c0300 + 2*c0210 - c0120 + c0021 + c0201))/2;
c1011 = (S.gw(:,2).*(-c0030 + 3*c1020 - 3*c2010 + c3000) + (-c0030 + 2*c1020 - c2010 + c2001 + c0021))/2;
c1101 = (S.gw(:,3).*(-c3000 + 3*c2100 - 3*c1200 + c0300) + (-c3000 + 2*c2100 - c1200 + c2001 + c0201))/2;
c1002 = (c1101 + c1011 + c2001)/3;
c0102 = (c1101 + c0111 + c0201)/3;
c0012 = (c1011 + c0111 + c0021)/3;
c0003 = (c1002 + c0102 + c0012)/3;

% --- Evaluate --------------------------------------------------------
t  = S.tq;
b1 = S.b1; b2 = S.b2; b3 = S.b3; b4 = S.b4;
z = b1.^3.*c3000(t)        + 3*b1.^2.*b2.*c2100(t) + 3*b1.^2.*b3.*c2010(t) + ...
    3*b1.^2.*b4.*c2001(t)  + 3*b1.*b2.^2.*c1200(t) + 6*b1.*b2.*b4.*c1101(t) + ...
    3*b1.*b3.^2.*c1020(t)  + 6*b1.*b3.*b4.*c1011(t) + 3*b1.*b4.^2.*c1002(t) + ...
    b2.^3.*c0300(t)        + 3*b2.^2.*b3.*c0210(t) + 3*b2.^2.*b4.*c0201(t) + ...
    3*b2.*b3.^2.*c0120(t)  + 6*b2.*b3.*b4.*c0111(t) + 3*b2.*b4.^2.*c0102(t) + ...
    b3.^3.*c0030(t)        + 3*b3.^2.*b4.*c0021(t) + 3*b3.*b4.^2.*c0012(t) + ...
    b4.^3.*c0003(t);

Zi = nan(S.size);
Zi(S.in) = z;

end

%************************************************************
%              VALUES INCLUDING THE BOUNDARY RING
%************************************************************
function V = fit_values(S, v)
%FIT_VALUES  Concatenate electrode values with the synthetic points' border values
%
%   Inputs:
%       S : struct - cache from topo_cache -- required
%       v : Nx1 double - values at the finite-valued electrodes -- required
%
%   Outputs:
%       V : (N+M)x1 double - electrode values followed by the boundary values

if S.n_extra == 0
    V = v;
    return
end

if ~ischar(S.border)
    V = [v; repmat(double(S.border), S.n_extra, 1)];
    return
end

% border == 'mean', using the cached real-electrode neighbour lists
v_extra = zeros(S.n_extra, 1);
used    = false(S.n_extra, 1);
for ii = 1:S.n_extra
    nb = S.ring_nbr{ii};
    if ~isempty(nb)
        used(ii)    = true;
        v_extra(ii) = mean(v(nb));
    end
end
if any(used) && ~all(used)
    v_extra(~used) = mean(v_extra(used));
end

V = [v; v_extra];

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
%
%   Notes:
%       Drawn in 2-D. The map is an image pushed to the bottom of the child list, so
%       these lines sit above it by draw order rather than by a z offset.

linewidth = 1.7;

% --- Head outline ---------------------------------------------------
t = linspace(0, 2*pi, 101);   % 101 points, as MNE's _make_head_outlines uses
hx = headrad * cos(t);
hy = headrad * sin(t);
plot(hx, hy, '-', ...
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
plot(nose_x, nose_y, '-', ...
    'Color', headcolor, 'LineWidth', linewidth, 'Parent', hax);

% --- Ears -----------------------------------------------------------
% MNE's ear outline, which it in turn inherited from EEGLAB's EarX/EarY. The coefficients
% are fractions of the head diameter, hence the 2*headrad scaling. The polyline is open at
% both ends by design: its first and last points sit just inside the rim at 0.994*headrad.
ear_x = [0.497 0.510 0.518 0.5299 0.5419 0.54 0.547 0.532 0.510 0.489] * (2 * headrad);
ear_y = [0.0555 0.0775 0.0783 0.0746 0.0555 -0.0055 -0.0932 -0.1313 -0.1384 -0.1199] * (2 * headrad);

plot(ear_x, ear_y, '-', ...
    'Color', headcolor, 'LineWidth', linewidth, 'Parent', hax);   % right ear
plot(-ear_x, ear_y, '-', ...
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
