function lib = build_montage_library(varargin)
%BUILD_MONTAGE_LIBRARY  Compile vendor montage files into the shipped library
%
%   Usage:
%       build_montage_library
%       lib = build_montage_library('outfile', '/tmp/test_lib.mat')
%
%   Inputs:
%       none
%
%   Name-Value Pairs:
%       'vendordir' : char - directory tree of raw vendor files to compile
%                     (default: channel_locs/vendor)
%       'outfile'   : char - where to write the library
%                     (default: channel_locs/montage_library.mat)
%       'verbose'   : logical - report progress (default: true)
%
%   Outputs:
%       lib : struct - the library that was written. See LIBRARY FORMAT below
%
%   Notes:
%       Run this once, when vendor files are added or updated. Everyday use goes
%       through cap_montage, which reads the built library and never touches the
%       raw files -- so the ~1.4 MB of vendor XML does not need to be present,
%       or distributed, for the montages to work.
%
%       Two montages are derived rather than read from a file, because they are
%       clinical standards rather than products and no vendor ships them:
%       '10-20' (the 19-channel system) and 'AASM-6' (F3 F4 C3 C4 O1 O2). Both
%       are cut from the 21-channel actiCAP, so their electrode positions are
%       the vendor's; only the selection and the conventional ordering are ours.
%
%       Deduplication is by ordered channel list AND coordinates. Vendors ship
%       the same arrangement under many product names -- eleven 32-channel codes
%       share one layout -- so the library stores each distinct layout once and
%       maps every name onto it. 109 names collapse to 65 layouts.
%
%       Both halves of that signature are load-bearing. Order matters because
%       eloc(k) must line up with the data's channel k. Geometry matters because
%       an identical channel list is NOT always the same cap: AP-32 and
%       AS-32_NO_REF list the same channels in the same order yet place FC5
%       1.5e-02 apart. Deduping on labels alone merges them and silently moves
%       electrodes between cap families.
%
%       For the same reason, positions are not factored into a shared name table
%       even though that would be smaller: 19 scalp names -- FC5, FCz, and the
%       F11/FT11/TP11 series among them -- genuinely sit in different places on
%       different cap families. FC5 alone has three distinct positions across
%       90 caps.
%
%   -----------------------------------------------------------------------
%   LIBRARY FORMAT
%   -----------------------------------------------------------------------
%   lib.layouts : 1xL struct, one per distinct ordered channel list
%       .labels : 1xC cell   - channel labels, in the vendor's channel order
%       .xyz    : Cx3 double - ALS unit-sphere coordinates
%       .vendor : char       - 'brainvision' | 'egi' | ...
%       .flags  : struct     - what read_montage's validation found:
%                              .coincident (logical), .n_coincident (double)
%   lib.names   : 1xM cell   - every selectable montage name
%   lib.map     : 1xM double - index into lib.layouts for each name
%   lib.built   : char       - ISO date the library was compiled
%   lib.source  : char       - where the raw files came from
%
%   Only labels and xyz are stored; theta/radius and the spherical fields are
%   derived by xyz_to_eloc on load. Storing them would let the library drift
%   out of step with read_montage's projection.
%
%   Example:
%       build_montage_library;
%       eloc = cap_montage('AC-64');
%
%   See also: cap_montage, read_montage, xyz_to_eloc
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%************************************************************
%                      INPUT HANDLING
%************************************************************

here = fileparts(mfilename('fullpath'));

p = inputParser;
addParameter(p, 'vendordir', fullfile(here,'channel_locs','vendor'), ...
    @(x) validateattributes(x, {'char','string'}, {'scalartext'}));
addParameter(p, 'outfile', fullfile(here,'channel_locs','montage_library.mat'), ...
    @(x) validateattributes(x, {'char','string'}, {'scalartext'}));
addParameter(p, 'verbose', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
parse(p, varargin{:});

vendordir = char(p.Results.vendordir);
outfile   = char(p.Results.outfile);
verbose   = p.Results.verbose;

assert(exist(vendordir,'dir') == 7, 'build_montage_library:noVendorDir', ...
    ['Vendor directory not found: %s\n' ...
     'The raw files are only needed to build the library -- see ' ...
     'channel_locs/vendor/README.md for how to fetch them.'], vendordir);

%************************************************************
%                     FIND VENDOR FILES
%************************************************************

files = {};
for ext = {'*.bvef','*.sfp','*.elc'}
    found = dir(fullfile(vendordir, '**', ext{1}));
    for ii = 1:numel(found)
        files{end+1,1} = fullfile(found(ii).folder, found(ii).name); %#ok<AGROW>
    end
end

assert(~isempty(files), 'build_montage_library:noFiles', ...
    'No .bvef/.sfp/.elc files found under %s', vendordir);

if verbose
    fprintf('compiling %d vendor files from %s\n', numel(files), vendordir);
end

%************************************************************
%                       READ THEM ALL
%************************************************************

names   = cell(numel(files),1);
vendors = cell(numel(files),1);
labels  = cell(numel(files),1);
xyzs    = cell(numel(files),1);
flags   = cell(numel(files),1);
keep    = true(numel(files),1);

for ii = 1:numel(files)
    [folder, base] = fileparts(files{ii});
    parts = strsplit(folder, filesep);

    % read_montage's warnings are the point of the build, not noise: capture
    % them per montage so callers see them without re-parsing the raw file
    lastwarn('', '');
    ws = warning('off','all');
    try
        e = read_montage(files{ii});
    catch ME
        warning(ws);
        keep(ii) = false;
        fprintf(2, '  SKIP %-22s [%s] %s\n', base, ME.identifier, ME.message);
        continue;
    end
    warning(ws);
    [~, wid] = lastwarn;

    names{ii}   = base;
    vendors{ii} = parts{end};
    labels{ii}  = {e.labels};
    xyzs{ii}    = [[e.X]', [e.Y]', [e.Z]'];

    f = struct('coincident', strcmp(wid,'read_montage:coincidentElectrodes'), ...
               'n_coincident', 0);
    if f.coincident
        f.n_coincident = count_coincident(xyzs{ii});
    end
    flags{ii} = f;
end

names = names(keep); vendors = vendors(keep);
labels = labels(keep); xyzs = xyzs(keep); flags = flags(keep);

%************************************************************
%                   DERIVED STANDARD MONTAGES
%************************************************************

% The clinical standards are not vendor products, so no cap file carries them,
% but they are what most people actually want by name. Cut them from the
% 21-channel actiCAP, which was verified to contain the full 10-20 (19/19) and
% the AASM set (6/6). Positions are the vendor's; only the selection and order
% are ours, and the order is the conventional clinical listing rather than the
% cap's own numbering.
STANDARDS = {
    '10-20', {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8', ...
              'P7','P3','Pz','P4','P8','O1','O2'}
    'AASM-6', {'F3','F4','C3','C4','O1','O2'}
    };

src = find(strcmp(names, 'AC-21'), 1);
if isempty(src)
    warning('build_montage_library:noStandardSource', ...
        ['AC-21 not found, so the derived standard montages (%s) were not built. ' ...
         'They are cut from that cap.'], strjoin(STANDARDS(:,1)', ', '));
else
    for ii = 1:size(STANDARDS,1)
        want = STANDARDS{ii,2};
        [tf, loc] = ismember(lower(want), lower(labels{src}));
        assert(all(tf), 'build_montage_library:standardMissing', ...
            'AC-21 lacks %s, needed for the %s montage.', ...
            strjoin(want(~tf), ' '), STANDARDS{ii,1});
        names{end+1,1}   = STANDARDS{ii,1};           %#ok<AGROW>
        vendors{end+1,1} = 'standard';                %#ok<AGROW>
        labels{end+1,1}  = want;                      %#ok<AGROW>
        xyzs{end+1,1}    = xyzs{src}(loc, :);         %#ok<AGROW>
        flags{end+1,1}   = struct('coincident', false, 'n_coincident', 0); %#ok<AGROW>
    end
end

%************************************************************
%                   DEDUPE BY ORDERED LIST
%************************************************************

% The signature is the ordered label list AND the coordinates. Both are needed:
% order matters because eloc(k) must line up with the data's channel k, and
% geometry matters because the same channel list is not always the same cap.
% AP-32 and AS-32_NO_REF list identical channels in identical order yet place
% FC5 1.5e-02 apart -- deduping on labels alone silently merges them and moves
% electrodes between cap families.
sig = cell(numel(labels),1);
for ii = 1:numel(labels)
    sig{ii} = sprintf('%s#%s', strjoin(labels{ii},'|'), ...
        sprintf('%.9f,', round(xyzs{ii}(:), 9)));
end
[~, first, map] = unique(sig, 'stable');

layouts = struct('labels', {}, 'xyz', {}, 'vendor', {}, 'flags', {});
for jj = 1:numel(first)
    k = first(jj);
    layouts(jj).labels = labels{k};
    layouts(jj).xyz    = xyzs{k};  
    layouts(jj).vendor = vendors{k};
    layouts(jj).flags  = flags{k}; 
end

% Cross-check the dedupe rather than trusting it: identical signatures must also
% have identical coordinates, or the signature is not capturing what matters.
for jj = 1:numel(first)
    same = find(map == jj);
    for k = same(:)'
        d = max(abs(xyzs{k}(:) - layouts(jj).xyz(:)));
        assert(d < 1e-12, 'build_montage_library:aliasMismatch', ...
            ['%s and %s share a channel list but their coordinates differ by %.2e. ' ...
             'They are not the same layout.'], names{k}, names{first(jj)}, d);
    end
end

lib = struct();
lib.layouts = layouts;
lib.names   = names(:)';
lib.map     = map(:)';
lib.built   = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
lib.source  = 'vendor files; see channel_locs/vendor/README.md';

%************************************************************
%                          WRITE
%************************************************************

save(outfile, 'lib', '-v7');   % -v7 compresses; the label lists are repetitive

if verbose
    d = dir(outfile);
    raw = sum(cellfun(@(f) subsref(dir(f), substruct('.','bytes')), files));
    fprintf('  %d names -> %d distinct layouts\n', numel(lib.names), numel(lib.layouts));
    fprintf('  channels: %d..%d\n', min(cellfun(@numel,{layouts.labels})), ...
        max(cellfun(@numel,{layouts.labels})));
    nco = sum(arrayfun(@(L) L.flags.coincident, layouts));
    fprintf('  layouts with coincident (non-scalp) electrodes: %d\n', nco);
    fprintf('  wrote %s  (%.0f KB, from %.1f MB of raw files -- %.0fx smaller)\n', ...
        outfile, d.bytes/1024, raw/1e6, raw/d.bytes);
end

end

%************************************************************
%                  COUNT COINCIDENT PAIRS
%************************************************************
function n = count_coincident(xyz)
%COUNT_COINCIDENT  Number of electrode pairs sharing a position
%
%   Inputs:
%       xyz : Cx3 double - electrode coordinates -- required
%
%   Outputs:
%       n : double - count of coincident pairs

D = squeeze(sqrt(sum((permute(xyz,[1 3 2]) - permute(xyz,[3 1 2])).^2, 3)));
D(1:size(xyz,1)+1:end) = inf;
n = nnz(triu(D < 1e-9));

end
