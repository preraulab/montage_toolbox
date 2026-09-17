function eloc = read_montage(filename, varargin)
%READ_MONTAGE  Read a manufacturer electrode-location file into an eloc struct
%
%   Usage:
%       eloc = read_montage('biosemi64.sfp')
%       eloc = read_montage('waveguard64.elc', 'orient', 'ras')
%       eloc = read_montage(file, 'center', 'fit')
%
%   Inputs:
%       filename : char - path to a vendor electrode file. The format is taken
%                  from the extension: '.sfp' (EGI geodesic nets), '.elc' (ASA
%                  electrode file, used by ANT Neuro), or '.bvef'
%                  (BrainVision / EasyCap) -- required
%
%   Name-Value Pairs:
%       'orient'   : char - axis convention OF THE FILE. 'als' is +x anterior,
%                    +y left, +z up (the EEGLAB/eloc convention). 'ras' is
%                    +x right, +y anterior, +z up. 'auto' infers it from 10-20
%                    landmark labels and errors if it cannot (default: 'auto')
%       'center'   : char - 'origin' projects about the file's own origin,
%                    which is correct for an idealized cap defined on a sphere.
%                    'fit' least-squares fits a sphere to the electrodes and
%                    projects about that centre, which is what a subject
%                    digitization needs (default: 'origin')
%       'validate' : logical - sanity-check the montage against 10-20 landmarks
%                    and error if it is implausible (default: true)
%
%   Outputs:
%       eloc : 1xC struct - channel locations in the convention plot_topo
%              expects, with fields labels, X, Y, Z, theta, radius,
%              sph_theta, sph_phi, sph_radius
%
%   Notes:
%       Vendors disagree on axis order and units, and a file that is silently
%       mis-oriented still plots -- it just plots wrong, which is exactly the
%       failure this function exists to prevent. So orientation is resolved
%       from the data rather than assumed: 'auto' locates 10-20 landmarks by
%       label and picks the convention that puts them where anatomy says they
%       go. See ORIENTATION INFERENCE below.
%
%       Units are normalized away entirely -- the projection depends only on
%       angles, so mm vs cm vs m does not matter and no unit conversion is done.
%
%       'center' is the knob for the distortion that motivated this function.
%       radius = 0.5 - elevation/pi measures elevation from the coordinate
%       origin, so it is only meaningful when the origin IS the head centre.
%       Canonical vendor caps are defined on a sphere about the origin, so
%       'origin' is right for them. Digitized montages are not, and want 'fit'.
%
%       GND, REF, CMS and DRL are dropped alongside the fiducials: they are
%       real electrode positions but not recorded data channels, so keeping
%       them would put a marker on the map for a channel the data does not have.
%
%       EOG/EMG/ECG channels are NOT dropped -- they are recorded channels, and
%       removing them would break the correspondence between eloc(k) and the
%       data's channel k. Be aware that vendors give them placeholder positions
%       (BrainVision stacks EMG1/EMG2/EMG3/ECG on top of Fpz), which is why a
%       coincidentElectrodes warning fires on those montages. Exclude them
%       before interpolating a scalp map.
%
%   -----------------------------------------------------------------------
%   ORIENTATION INFERENCE  (how 'auto' decides)
%   -----------------------------------------------------------------------
%   Both candidate conventions are tried; each is scored on how well the
%   landmarks it can find land where they should, in the ALS frame:
%
%       Cz          -> near the vertex (small radius)
%       Fpz/Fp1/Fp2 -> anterior  (+x)
%       Oz/O1/O2    -> posterior (-x)
%       T7/T3       -> left      (+y)
%       T8/T4       -> right     (-y)
%
%   The higher-scoring convention wins. A tie, or too few landmarks to score,
%   is an error rather than a guess -- pass 'orient' explicitly in that case.
%
%   Example:
%       eloc = read_montage('standard_64.sfp');
%       plot_topo([], eloc, 'electrodes', 'ptslabels');
%
%   See also: plot_topo, cap_montage, xyz_to_eloc, build_montage_library
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%************************************************************
%                      INPUT HANDLING
%************************************************************

p = inputParser;
addRequired(p,  'filename',   @(x) validateattributes(x, {'char','string'}, {'scalartext'}));
addParameter(p, 'orient',     'auto',   @(x) any(validatestring(x, {'auto','als','ras'})));
addParameter(p, 'center',     'origin', @(x) any(validatestring(x, {'origin','fit'})));
addParameter(p, 'validate',   true,     @(x) validateattributes(x, {'logical'}, {'scalar'}));
parse(p, filename, varargin{:});

filename = char(p.Results.filename);
do_validate = p.Results.validate;

assert(exist(filename, 'file') == 2, 'read_montage:fileNotFound', ...
    'Montage file not found: %s', filename);

[~, ~, ext] = fileparts(filename);

switch lower(ext)
    case '.sfp'
        [labels, xyz] = parse_sfp(filename);
    case '.elc'
        [labels, xyz] = parse_elc(filename);
    case '.bvef'
        [labels, xyz] = parse_bvef(filename);
    otherwise
        error('read_montage:unsupportedFormat', ...
            'Unsupported montage format ''%s''. Supported: .sfp, .elc, .bvef', ext);
end

assert(~isempty(labels), 'read_montage:noChannels', ...
    'No electrodes were parsed from %s.', filename);

%************************************************************
%                   DROP NON-ELECTRODES
%************************************************************

% Fiducials and headshape points are not channels, and they sit far off the cap
% (nasion/LPA/RPA are at the equator or below). Leaving them among the
% electrodes would drag the sphere fit and the plot radius around.
%
% They are kept for orientation first, though, and that ordering matters: an EGI
% HydroCel .sfp labels its electrodes E1..E64 and carries exactly one 10-20 name
% (Cz), which sits on the z axis and so cannot distinguish ALS from RAS. Its
% FidNz/FidT9/FidT10 are the only usable orientation signal in the file.
is_fid = is_fiducial(labels);
labels_orient = labels;
xyz_orient    = xyz;

if any(is_fid)
    labels = labels(~is_fid);
    xyz    = xyz(~is_fid, :);
end

assert(~isempty(labels), 'read_montage:onlyFiducials', ...
    'File %s contained only fiducials/headshape points, no electrodes.', filename);

%************************************************************
%                       ORIENTATION
%************************************************************

switch p.Results.orient
    case 'auto'
        convention = infer_convention(labels_orient, xyz_orient, filename);
    otherwise
        convention = p.Results.orient;
end

if strcmpi(convention, 'ras')
    xyz = ras_to_als(xyz);
    xyz_orient = ras_to_als(xyz_orient);
end

%************************************************************
%                        PROJECTION
%************************************************************

% Fit the sphere regardless of the centring mode: under 'fit' it defines the
% projection, and under 'origin' it is still the diagnostic that says whether
% trusting the origin was reasonable.
[c, R_fit] = fit_sphere_centre(xyz);

if strcmpi(p.Results.center, 'fit')
    xyz = xyz - c;
    xyz_orient = xyz_orient - c;   % keep the landmark frame in step with the electrodes
    centre_offset = 0;
else
    centre_offset = norm(c) / R_fit;
end

eloc = xyz_to_eloc(labels, xyz);

if do_validate
    % Direction checking uses the landmark set (electrodes + fiducials), since a
    % geodesic net's only orientation evidence lives in its fiducials.
    validate_montage(eloc, xyz_orient, labels_orient, centre_offset, filename);
end

end

%************************************************************
%                      PARSE .SFP FILES
%************************************************************
function [labels, xyz] = parse_sfp(filename)
%PARSE_SFP  Read an EGI/BioSemi simple electrode file
%
%   Inputs:
%       filename : char - path to a .sfp file -- required
%
%   Outputs:
%       labels : Cx1 cell - electrode label strings
%       xyz    : Cx3 double - raw coordinates, in the file's own convention

% Format is one electrode per line: "<label> <x> <y> <z>", whitespace separated.
% Blank lines and '#' comments are tolerated.
fid = fopen(filename, 'r');
assert(fid > 0, 'read_montage:cannotOpen', 'Could not open %s', filename);
cleaner = onCleanup(@() fclose(fid));

labels = {};
xyz = [];
line_num = 0;

while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    line_num = line_num + 1;

    tline = strtrim(tline);
    if isempty(tline) || tline(1) == '#'
        continue;
    end

    parts = strsplit(tline);
    if numel(parts) < 4
        continue;   % not an electrode row
    end

    coords = str2double(parts(2:4));
    if any(isnan(coords))
        error('read_montage:badSfpLine', ...
            '%s line %d: expected "<label> <x> <y> <z>", got "%s".', ...
            filename, line_num, tline);
    end

    labels{end+1,1} = parts{1};  %#ok<AGROW>
    xyz(end+1,:) = coords;       %#ok<AGROW>
end

end

%************************************************************
%                      PARSE .ELC FILES
%************************************************************
function [labels, xyz] = parse_elc(filename)
%PARSE_ELC  Read an ASA electrode file (ANT Neuro waveguard and others)
%
%   Inputs:
%       filename : char - path to a .elc file -- required
%
%   Outputs:
%       labels : Cx1 cell - electrode label strings
%       xyz    : Cx3 double - raw coordinates, in the file's own convention

% Layout is a "Positions" block of "x y z" rows followed by a "Labels" block.
% The two blocks are parallel, which is why they are read separately and then
% length-checked against each other rather than zipped as they are read.
txt = fileread(filename);
lines = strsplit(txt, {'\r\n', '\n', '\r'});

mode = '';
xyz = [];
labels = {};

for ii = 1:numel(lines)
    tline = strtrim(lines{ii});
    if isempty(tline) || tline(1) == '#'
        continue;
    end

    if strncmpi(tline, 'Positions', 9)
        mode = 'pos';
        continue;
    elseif strncmpi(tline, 'Labels', 6)
        mode = 'lab';
        continue;
    elseif ~isempty(regexpi(tline, '^(ReferenceLabel|UnitPosition|NumberPositions|Polygons|HeadShape)', 'once'))
        mode = '';   % header/other block
        continue;
    end

    switch mode
        case 'pos'
            % Rows may optionally carry a leading "Label:" prefix
            tline = regexprep(tline, '^\s*\S+\s*:\s*', '');
            vals = sscanf(tline, '%f');
            if numel(vals) >= 3
                xyz(end+1,:) = vals(1:3)';  %#ok<AGROW>
            end
        case 'lab'
            parts = strsplit(tline);
            parts = parts(~cellfun(@isempty, parts));
            labels = [labels; parts(:)];  %#ok<AGROW>
    end
end

assert(~isempty(xyz), 'read_montage:noPositions', ...
    'No "Positions" block found in %s.', filename);

if isempty(labels)
    % A positions-only file is usable, but the labels are what make orientation
    % inference and validation possible, so say so rather than silently guessing.
    error('read_montage:noLabels', ...
        ['No "Labels" block found in %s. Labels are required to verify the ' ...
         'montage orientation; pass ''orient'' explicitly and ''validate'',false ' ...
         'if you really want to load it unchecked.'], filename);
end

assert(numel(labels) == size(xyz,1), 'read_montage:blockMismatch', ...
    '%s: %d labels but %d positions -- the blocks do not correspond.', ...
    filename, numel(labels), size(xyz,1));

end

%************************************************************
%                     PARSE .BVEF FILES
%************************************************************
function [labels, xyz] = parse_bvef(filename)
%PARSE_BVEF  Read a BrainVision / EasyCap electrode file
%
%   Inputs:
%       filename : char - path to a .bvef file -- required
%
%   Outputs:
%       labels : Cx1 cell - electrode label strings
%       xyz    : Cx3 double - coordinates in the RAS convention
%
%   Notes:
%       BVEF is XML, but the schema is a flat list of <Electrode> blocks, so it
%       is read with regexp rather than an XML parser -- no toolbox or Java
%       dependency, and nothing here needs the DOM.

txt = fileread(filename);

blocks = regexp(txt, '<Electrode>(.*?)</Electrode>', 'tokens');
assert(~isempty(blocks), 'read_montage:noElectrodes', ...
    'No <Electrode> blocks found in %s.', filename);

labels = cell(numel(blocks),1);
xyz = zeros(numel(blocks),3);
keep = true(numel(blocks),1);

for ii = 1:numel(blocks)
    b = blocks{ii}{1};

    name  = regexp(b, '<Name>(.*?)</Name>',     'tokens', 'once');
    theta = regexp(b, '<Theta>(.*?)</Theta>',   'tokens', 'once');
    phi   = regexp(b, '<Phi>(.*?)</Phi>',       'tokens', 'once');
    rad   = regexp(b, '<Radius>(.*?)</Radius>', 'tokens', 'once');

    % An electrode with no angles is a placeholder for an unpopulated holder;
    % vendors do emit these, and they are not positions.
    if isempty(name) || isempty(theta) || isempty(phi)
        keep(ii) = false;
        continue;
    end

    th = str2double(theta{1});
    ph = str2double(phi{1});
    if isempty(rad)
        r = 1;
    else
        r = str2double(rad{1});
    end
    if isnan(th) || isnan(ph) || isnan(r)
        keep(ii) = false;
        continue;
    end

    % EasyCap spherical convention, verified against a real actiCAP AC-21 file:
    % Theta is the signed polar angle down from +z (the vertex), negative on the
    % left hemisphere; Phi is the azimuth of the xy projection measured from +x.
    % The result is RAS -- Cz(0,0) -> (0,0,r); T7(-90,0) -> (-r,0,0) = left;
    % Fz(45,90) -> (0,+.71r,+.71r) = anterior and up; TP9(-113,18) lands below
    % the equator at the left mastoid, as it should.
    th = th * pi/180;
    ph = ph * pi/180;
    xyz(ii,:) = r * [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)];
    labels{ii} = strtrim(name{1});
end

labels = labels(keep);
xyz = xyz(keep,:);

end

%************************************************************
%                   IDENTIFY NON-ELECTRODES
%************************************************************
function tf = is_fiducial(labels)
%IS_FIDUCIAL  Flag fiducial / headshape entries that are not real channels
%
%   Inputs:
%       labels : Cx1 cell - electrode label strings -- required
%
%   Outputs:
%       tf : Cx1 logical - true where the label names a fiducial

% Vendors spell these several ways; match the common set case-insensitively
fid_names = {'nz','nas','nasion','fidnz','lpa','rpa','fidt9','fidt10', ...
             'al','ar','a1','a2','m1','m2','cms','drl','ref','gnd','com'};
tf = false(numel(labels),1);
for ii = 1:numel(labels)
    lab = lower(strtrim(labels{ii}));
    tf(ii) = any(strcmp(lab, fid_names)) || strncmp(lab, 'fid', 3);
end

end

%************************************************************
%                  LANDMARK LOOKUP TABLE
%************************************************************
function idx = find_label(labels, names)
%FIND_LABEL  Index of the first label matching any of the given names
%
%   Inputs:
%       labels : Cx1 cell - electrode label strings -- required
%       names  : 1xN cell - acceptable names, tried in order -- required
%
%   Outputs:
%       idx : scalar double - index into labels, or [] if none matched

idx = [];
lab = strtrim(labels);
for ii = 1:numel(names)
    hit = find(strcmpi(lab, names{ii}), 1);
    if ~isempty(hit)
        idx = hit;
        return;
    end
end

end

%************************************************************
%                  SCORE AN ORIENTATION
%************************************************************
function [score, n_used] = score_orientation(labels, xyz)
%SCORE_ORIENTATION  Rate how well landmarks point the way anatomy says they should
%
%   Inputs:
%       labels : Cx1 cell - electrode label strings -- required
%       xyz    : Cx3 double - coordinates interpreted as ALS -- required
%
%   Outputs:
%       score  : double - mean directional cosine in [-1 1], 1 is perfect
%       n_used : integer - how many landmarks were found and scored

% In the ALS frame: +x anterior, +y left, +z up. Each landmark contributes the
% cosine between its own position vector and the axis it ought to point along.
% Normalizing by the electrode's OWN radius (not a global scale) makes each term
% a true direction cosine, so the score is a mean of comparable quantities and a
% single distant electrode cannot dominate.
%
% Cz is deliberately excluded: it sits on the +z axis, so it is invariant under
% the rotation about z that separates ALS from RAS and contributes nothing to
% telling them apart. It is checked separately in validate_montage.
%
% Fiducials are included because they are often the ONLY usable signal: a
% geodesic net labels its electrodes E1..E256 and carries no 10-20 names at all,
% but always carries nasion and the two preauricular points.
tests = {
    {'Fpz','Fp1','Fp2'},            1, +1   % anterior
    {'Oz','O1','O2'},               1, -1   % posterior
    {'T7','T3'},                    2, +1   % left
    {'T8','T4'},                    2, -1   % right
    {'FidNz','Nz','Nas','Nasion'},  1, +1   % nasion  -> anterior
    {'FidT9','LPA','T9'},           2, +1   % left preauricular
    {'FidT10','RPA','T10'},         2, -1   % right preauricular
    };

total = 0;
n_used = 0;
for ii = 1:size(tests,1)
    idx = find_label(labels, tests{ii,1});
    if isempty(idx), continue; end
    r = norm(xyz(idx,:));
    if r <= 0, continue; end
    axis_i = tests{ii,2};
    want   = tests{ii,3};
    total = total + want * xyz(idx, axis_i) / r;
    n_used = n_used + 1;
end

if n_used == 0
    score = -inf;
else
    score = total / n_used;
end

end

%************************************************************
%                  INFER THE AXIS CONVENTION
%************************************************************
function convention = infer_convention(labels, xyz_raw, filename)
%INFER_CONVENTION  Choose between the ALS and RAS conventions from landmarks
%
%   Inputs:
%       labels   : Cx1 cell - labels INCLUDING fiducials -- required
%       xyz_raw  : Cx3 double - coordinates in the file's own convention
%                  -- required
%       filename : char - source file, for error messages -- required
%
%   Outputs:
%       convention : char - 'als' or 'ras'

[s_als, n_als] = score_orientation(labels, xyz_raw);
[s_ras, n_ras] = score_orientation(labels, ras_to_als(xyz_raw));

n_used = max(n_als, n_ras);
if n_used < 2
    error('read_montage:cannotInferOrient', ...
        ['%s: cannot infer the axis convention -- found only %d usable landmark(s). ' ...
         'Need at least 2 of: Fpz/Fp1/Fp2, Oz/O1/O2, T7/T3, T8/T4, or the ' ...
         'fiducials FidNz/Nz, FidT9/LPA, FidT10/RPA. ' ...
         'Pass ''orient'',''als'' or ''orient'',''ras'' explicitly.'], filename, n_used);
end

% Require a clear winner. A near-tie means the landmarks did not actually
% discriminate, and picking one at random is how a montage plots wrong.
if abs(s_als - s_ras) < 0.25
    error('read_montage:ambiguousOrient', ...
        ['%s: axis convention is ambiguous (ALS score %.3f vs RAS score %.3f over ' ...
         '%d landmarks). Pass ''orient'' explicitly.'], ...
        filename, s_als, s_ras, n_used);
end

if s_als >= s_ras
    convention = 'als';
else
    convention = 'ras';
end

end

%************************************************************
%                    RAS -> ALS CONVERSION
%************************************************************
function xyz = ras_to_als(xyz_ras)
%RAS_TO_ALS  Rotate right/anterior/superior coordinates into the eloc frame
%
%   Inputs:
%       xyz_ras : Cx3 double - +x right, +y anterior, +z up -- required
%
%   Outputs:
%       xyz : Cx3 double - +x anterior, +y left, +z up

% anterior <- +y_ras;  left <- -x_ras;  up is shared
xyz = [xyz_ras(:,2), -xyz_ras(:,1), xyz_ras(:,3)];

end

%************************************************************
%                  FIT A SPHERE TO THE CAP
%************************************************************
function [c, R] = fit_sphere_centre(xyz)
%FIT_SPHERE_CENTRE  Least-squares sphere centre for a set of electrodes
%
%   Inputs:
%       xyz : Cx3 double - electrode coordinates -- required
%
%   Outputs:
%       c : 1x3 double - fitted sphere centre
%       R : double - fitted sphere radius

% Algebraic (Kasa) fit: |p-c|^2 = R^2 expands to a linear system in (c, R^2-|c|^2),
% which avoids an iterative solve. Good enough here -- the electrodes cover most
% of a hemisphere, which is the regime where this fit is well conditioned.
assert(size(xyz,1) >= 4, 'read_montage:tooFewForFit', ...
    'Need at least 4 electrodes to fit a sphere centre (%d given).', size(xyz,1));

A = [2*xyz, ones(size(xyz,1),1)];
b = sum(xyz.^2, 2);
w = A\b;
c = w(1:3)';
R = sqrt(max(w(4) + c*c', 0));

end

%************************************************************
%                   SANITY-CHECK THE RESULT
%************************************************************
function validate_montage(eloc, xyz_lm, labels_lm, centre_offset, filename)
%VALIDATE_MONTAGE  Error if the projected montage contradicts 10-20 anatomy
%
%   Inputs:
%       eloc          : 1xC struct - projected channel locations -- required
%       xyz_lm        : Mx3 double - ALS coordinates of the LANDMARK set
%                       (electrodes plus fiducials) -- required
%       labels_lm     : Mx1 cell - labels of the landmark set -- required
%       centre_offset : double - fitted sphere centre distance from the
%                       projection origin, as a fraction of the fitted radius
%                       -- required
%       filename      : char - source file, for error messages -- required
%
%   Outputs:
%       none (throws or warns on an implausible montage)

% Two label sets are in play and they are NOT interchangeable: the landmark set
% carries fiducials and is only for direction scoring, while radius checks must
% index the electrode set that eloc was built from.
labels = {eloc.labels};
radius = [eloc.radius];

% --- Landmark directions ---------------------------------------------
% This is the check that catches a mis-declared axis convention. Reading RAS as
% ALS is a 90-degree rotation about z, which leaves every radius and the overall
% aspect ratio untouched -- the cap looks perfectly plausible and is simply
% turned. Only the landmark DIRECTIONS reveal it.
[dir_score, n_dir] = score_orientation(labels_lm, xyz_lm);
if n_dir >= 2 && dir_score < 0.5
    error('read_montage:landmarksMisplaced', ...
        ['%s: 10-20 landmarks do not point where they should ' ...
         '(direction score %.2f over %d landmarks; 1.0 is perfect). ' ...
         'Fpz/Oz/T7/T8 should lie anterior/posterior/left/right. ' ...
         'The axis convention is probably wrong.'], ...
        filename, dir_score, n_dir);
end

% Every electrode should be on or above the head, i.e. inside the skirt. Beyond
% radius 1 is below the equator by more than a hemisphere, which means the
% projection is wrong rather than the cap being unusual.
if any(radius > 1)
    bad = find(radius > 1);
    error('read_montage:implausibleRadius', ...
        ['%s: %d electrode(s) project past radius 1.0 (max %.2f, e.g. ''%s''). ' ...
         'The axis convention or the sphere centre is probably wrong.'], ...
        filename, numel(bad), max(radius), labels{bad(1)});
end

% Cz, if present, pins the vertex. If it is not near the centre the projection
% is off regardless of what everything else does.
idx_cz = find_label(labels, {'Cz'});
if ~isempty(idx_cz) && radius(idx_cz) > 0.2
    error('read_montage:czOffCentre', ...
        ['%s: Cz projects to radius %.3f, but it should sit near the vertex ' ...
         '(radius ~0). The sphere centre is probably wrong -- try ''center'',''fit''.'], ...
        filename, radius(idx_cz));
end

% --- Coincident electrodes -------------------------------------------
% Vendors park non-scalp channels (EOG, EMG, ECG, IO) on a placeholder angle
% rather than a real position, so several land on exactly the same spot -- on a
% BrainVision sleep cap, EMG1/EMG2/EMG3/ECG all sit on top of Fpz. They cannot
% be dropped, because they ARE recorded data channels and removing them would
% break the one-to-one match between eloc(k) and the data's channel k. So warn
% instead: interpolating a map over them silently averages contradictory values
% at one location and returns a plausible-looking wrong answer rather than
% failing. 25 of the 109 montages in the library are affected.
ex_v = radius(:) .* sin([eloc.theta]'*pi/180);
ey_v = radius(:) .* cos([eloc.theta]'*pi/180);
Dsq = hypot(ex_v - ex_v', ey_v - ey_v');
Dsq(1:numel(eloc)+1:end) = inf;
[ra, ca] = find(Dsq < 1e-9);
if ~isempty(ra)
    pairs = unique(sort([ra ca], 2), 'rows');
    shown = min(4, size(pairs,1));
    desc = arrayfun(@(k) sprintf('%s=%s', labels{pairs(k,1)}, labels{pairs(k,2)}), ...
        1:shown, 'UniformOutput', false);
    warning('read_montage:coincidentElectrodes', ...
        ['%s: %d electrode pair(s) share the same position (%s%s). These are ' ...
         'usually non-scalp channels (EOG/EMG/ECG) given placeholder coordinates ' ...
         'by the vendor. They are kept so channel order still matches the data, ' ...
         'but exclude them before interpolating a scalp map -- otherwise the map ' ...
         'averages several different values at one point.'], ...
        filename, size(pairs,1), strjoin(desc, ', '), ...
        repmat(', ...', 1, size(pairs,1) > shown));
end

% --- Off-centre coordinates ------------------------------------------
% The decisive test for the distortion this function exists to prevent. An
% offset along z shrinks every electrode's radius by the SAME amount, so the cap
% stays perfectly circular and merely clusters inward -- no aspect-ratio or
% per-landmark radius check can see it. Comparing the fitted sphere centre to
% the projection origin catches it directly, and works even for a cap with no
% recognizable 10-20 labels. 5% is well below the level that is visible in a
% plot but well above the noise of a real digitization.
% Threshold calibrated against real vendor files rather than an idealized cap.
% Genuine EGI HydroCel nets sit 12-14% off the electrode set's best-fit sphere
% centre and are perfectly good -- their origin is the fiducial-defined head
% centre, and the net's coverage extends asymmetrically down over the face,
% which drags a sphere fit. The old picture-digitized eloc64 starburst (since
% removed from the repo) sat at 50%. 25% separates them with room to spare.
%
% Deliberately NOT checked: sphericity. A real head is not a sphere, so the
% residual cannot discriminate -- a legitimate EGI 256 net scores 8.1% and that
% same junk starburst scored 8.5%. A check that cannot tell good from bad is worse
% than no check, because people learn to ignore it.
if centre_offset > 0.25
    warning('read_montage:offCentre', ...
        ['%s: electrodes are not centred on the projection origin -- the ' ...
         'best-fit sphere centre is %.1f%% of its radius away. Projecting about ' ...
         'the origin distorts the layout (electrodes clustered inward, or ' ...
         'squashed along one axis). Try ''center'',''fit'', and check the result: ' ...
         'if the layout is still not circular, the coordinates themselves are suspect.'], ...
        filename, 100*centre_offset);
end

% A cap much wider than it is tall is the other signature of an off-centre
% projection -- specifically an offset in x or y rather than z.
ex = radius(:) .* sin([eloc.theta]'*pi/180);
ey = radius(:) .* cos([eloc.theta]'*pi/180);
span_x = max(ex) - min(ex);
span_y = max(ey) - min(ey);
if span_x > 0 && span_y > 0
    aspect = max(span_x/span_y, span_y/span_x);
    if aspect > 1.35
        warning('read_montage:anisotropic', ...
            ['%s: projected cap is %.2f:1 anisotropic (x span %.2f, y span %.2f). ' ...
             'A canonical cap should be close to circular -- real vendor nets ' ...
             'run about 1.1-1.2:1. See any offCentre warning above for the ' ...
             'likely cause.'], ...
            filename, aspect, span_x, span_y);
    end
end

end
