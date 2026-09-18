function [eloc, dig] = fif_to_eloc(filename, varargin)
%FIF_TO_ELOC  Read a Polhemus/Neuromag .fif digitization into an eloc struct
%
%   Usage:
%       eloc = fif_to_eloc('NSS05-2.fif')
%       eloc = fif_to_eloc(file, 'labels', {'1','2',...})
%       [eloc, dig] = fif_to_eloc(file, 'skipfirst', true)
%
%   Inputs:
%       filename : char - path to a .fif file holding digitizer points
%                  -- required
%
%   Name-Value Pairs:
%       'labels'    : cell - one label per EEG point, in order. Defaults to
%                     '1'..'N', which is what a cap with numbered rather than
%                     named electrodes wants (default: {})
%       'skipfirst' : logical - drop the first EEG point before building the
%                     eloc. Many montages digitize the reference electrode
%                     first, and it is not a recorded data channel
%                     (default: false)
%       'center'    : 1x3 double - subtract this from every point before
%                     projecting. [] uses the file's own origin, which for
%                     Neuromag head coordinates is already the correct head
%                     centre -- see Notes (default: [])
%
%   Outputs:
%       eloc : 1xC struct - channel locations ready for plot_topo
%       dig  : struct - everything read from the file, with fields
%              fiducials (struct with lpa/nasion/rpa), hpi, eeg, headshape,
%              and raw (Nx5: [kind ident x y z])
%
%   Notes:
%       Reads the FIF tag stream directly, so it needs no MNE, FieldTrip or
%       Python. FIF is big-endian, and every tag is a 16-byte header -- kind,
%       type, size, next -- followed by size bytes of payload. Digitizer
%       points are FIFF_DIG_POINT (kind 213) with a 20-byte payload holding
%       int32 kind, int32 ident and three float32 coordinates in metres.
%
%       Coordinate convention. Neuromag/MNE head coordinates are RAS: +x
%       toward the right preauricular point, +y toward the nasion, +z up. The
%       origin is the midpoint of the left/right preauricular points, so all
%       three fiducials lie in the z = 0 plane with LPA and RPA on the y axis.
%       That origin IS the interaural midpoint, which is the anatomically
%       right centre for the azimuthal projection -- so no re-centring is
%       needed and 'center' should normally be left empty. A least-squares
%       sphere fit is a poor substitute: an EEG cap covers the top of the head
%       rather than a sphere, so the fit is biased upward and pushes far more
%       electrodes past the head rim than a real cap has.
%
%       This function checks the fiducials against that convention and warns
%       if they do not match, since a file in some other frame would otherwise
%       project silently and wrongly.
%
%       xyz_to_eloc wants ALS (+x anterior, +y left, +z up), so the RAS input
%       is permuted to [y, -x, z] on the way in. Getting that wrong does not
%       error -- it mirrors or rotates the head, and a mirrored topography
%       looks entirely plausible.
%
%   Example:
%       eloc = fif_to_eloc('NSS05-2.fif', 'skipfirst', true);
%       plot_topo(values, eloc, 'style', 'map');
%
%   See also: xyz_to_eloc, read_montage, cap_montage, plot_topo
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%************************************************************
%                      INPUT HANDLING
%************************************************************

p = inputParser;
addRequired(p,  'filename',  @(x) validateattributes(x, {'char','string'}, {'scalartext'}));
addParameter(p, 'labels',    {},    @(x) validateattributes(x, {'cell'}, {}));
addParameter(p, 'skipfirst', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(p, 'center',    [],    @(x) isempty(x) || ...
                                         (isnumeric(x) && numel(x) == 3));
parse(p, filename, varargin{:});

filename  = char(p.Results.filename);
labels    = p.Results.labels;
skipfirst = p.Results.skipfirst;
center    = p.Results.center;

assert(exist(filename, 'file') == 2, 'fif_to_eloc:fileNotFound', ...
    'File not found: %s', filename);

%************************************************************
%                    READ THE TAG STREAM
%************************************************************

raw = read_dig_tags(filename);
assert(~isempty(raw), 'fif_to_eloc:noDigPoints', ...
    '%s contains no FIFF_DIG_POINT tags. Is it a digitization file?', filename);

% FIFFV_POINT_* : 1 cardinal (fiducial), 2 HPI, 3 EEG, 4 extra (head shape)
dig.raw       = raw;
dig.hpi       = raw(raw(:,1) == 2, 3:5);
dig.eeg       = raw(raw(:,1) == 3, 3:5);
dig.headshape = raw(raw(:,1) == 4, 3:5);

card = raw(raw(:,1) == 1, :);
dig.fiducials = struct('lpa', [], 'nasion', [], 'rpa', []);
for ii = 1:size(card,1)
    switch card(ii,2)
        case 1, dig.fiducials.lpa    = card(ii,3:5);
        case 2, dig.fiducials.nasion = card(ii,3:5);
        case 3, dig.fiducials.rpa    = card(ii,3:5);
    end
end

check_head_frame(dig.fiducials);

assert(~isempty(dig.eeg), 'fif_to_eloc:noEegPoints', ...
    '%s has digitizer points but none of kind EEG.', filename);

%************************************************************
%                      BUILD THE ELOC
%************************************************************

xyz = dig.eeg;
if skipfirst
    xyz = xyz(2:end, :);
end

if ~isempty(center)
    xyz = xyz - center(:)';
end

num_chans = size(xyz,1);
if isempty(labels)
    labels = cellstr(string(1:num_chans));
end
assert(numel(labels) == num_chans, 'fif_to_eloc:labelMismatch', ...
    '%d labels for %d EEG points%s.', numel(labels), num_chans, ...
    char("" + string(repmat(' (after skipfirst)', 1, skipfirst))));

% RAS in, ALS out: +x anterior is RAS y, +y left is minus RAS x, +z up unchanged
eloc = xyz_to_eloc(labels(:), [xyz(:,2), -xyz(:,1), xyz(:,3)]);

end

%************************************************************
%                  READ FIFF_DIG_POINT TAGS
%************************************************************
function raw = read_dig_tags(filename)
%READ_DIG_TAGS  Pull every digitizer point out of a FIF tag stream
%
%   Inputs:
%       filename : char - path to the .fif file -- required
%
%   Outputs:
%       raw : Nx5 double - one row per point, [kind ident x y z], metres
%
%   Notes:
%       Walks the file sequentially rather than following the tag 'next'
%       pointers. Digitization files are written as a flat stream, and a
%       sequential walk cannot loop or miss a block the way a corrupted
%       pointer chain can.

FIFF_DIG_POINT = 213;
DIG_PAYLOAD    = 20;     % int32 kind + int32 ident + 3 x float32

fid = fopen(filename, 'r', 'ieee-be');   % FIF is big-endian
assert(fid > 0, 'fif_to_eloc:cannotOpen', 'Could not open %s', filename);
closer = onCleanup(@() fclose(fid));

raw = zeros(0,5);
while true
    hdr = fread(fid, 4, 'int32');        % kind, type, size, next
    if numel(hdr) < 4
        break
    end
    kind = hdr(1);
    sz   = hdr(3);
    if sz < 0
        break                            % malformed length, stop rather than seek wildly
    end
    if kind == FIFF_DIG_POINT && sz == DIG_PAYLOAD
        pkind = fread(fid, 1, 'int32');
        pid   = fread(fid, 1, 'int32');
        r     = fread(fid, 3, 'float32');
        if numel(r) < 3
            break
        end
        raw(end+1, :) = [pkind pid r(:)'];  %#ok<AGROW>
    else
        if fseek(fid, sz, 'cof') ~= 0
            break
        end
    end
end

end

%************************************************************
%                  VERIFY THE HEAD FRAME
%************************************************************
function check_head_frame(f)
%CHECK_HEAD_FRAME  Warn if the fiducials do not look like Neuromag head coordinates
%
%   Inputs:
%       f : struct - fiducials with fields lpa, nasion, rpa -- required
%
%   Outputs:
%       none (warns only)
%
%   Notes:
%       In head coordinates LPA and RPA sit on the x axis either side of the
%       origin and the nasion sits on +y, all at z = 0. If that does not hold
%       the file is in some other frame, and the ALS permutation below would
%       rotate or mirror the montage without erroring.

if isempty(f.lpa) || isempty(f.nasion) || isempty(f.rpa)
    warning('fif_to_eloc:noFiducials', ...
        ['No cardinal points in this file, so the coordinate frame could not be ' ...
         'checked. Positions are assumed to be RAS head coordinates.']);
    return
end

tol = 5e-3;    % 5 mm, loose enough for digitizer noise, tight enough to catch a wrong frame
ok = abs(f.lpa(3)) < tol && abs(f.rpa(3)) < tol && abs(f.nasion(3)) < tol && ...
     abs(f.lpa(2)) < tol && abs(f.rpa(2)) < tol && ...
     f.lpa(1) < 0 && f.rpa(1) > 0 && f.nasion(2) > 0;

if ~ok
    warning('fif_to_eloc:unexpectedFrame', ...
        ['Fiducials do not match Neuromag head coordinates (LPA/RPA on -x/+x, ' ...
         'nasion on +y, all at z=0):\n   LPA %s\n   Nasion %s\n   RPA %s\n' ...
         'The montage may come out rotated or mirrored.'], ...
        mat2str(round(f.lpa,4)), mat2str(round(f.nasion,4)), mat2str(round(f.rpa,4)));
end

end
