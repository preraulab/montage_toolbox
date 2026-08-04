function varargout = cap_montage(name)
%CAP_MONTAGE  Look up an EEG cap montage by name
%
%   Usage:
%       eloc  = cap_montage('AC-64')
%       cap_montage list
%       names = cap_montage('list')
%       cap_montage layouts
%
%   Inputs:
%       name : char - the montage name (e.g. 'AC-64', 'BC-SL-64'), case
%              insensitive. 'list' enumerates every montage; 'layouts' shows
%              them grouped by distinct layout, so it is clear which names are
%              the same arrangement under different product codes. An unknown
%              name is reported with near-miss suggestions, and a partial name
%              matching several montages lists them (default: 'list')
%
%   Outputs:
%       varargout : with a montage name, the 1xC eloc struct ready for
%                   plot_topo or channelbrowse. With 'list', a cell array of
%                   names, or nothing (prints a table) when called without an
%                   output
%
%   Notes:
%       Reads channel_locs/montage_library.mat, a compiled store built from the
%       cap manufacturers' own files by build_montage_library. The raw vendor
%       files are not needed at run time and are not distributed.
%
%       Names are the vendors' product codes. The suffix conventions are theirs:
%
%           AC-64            actiCAP, 64 channels
%           BC-SL-64         BrainCap Sleep, 64 channels
%           AS-64_REF        variant including the reference electrode
%           AS-64_NO_REF     the same cap without it
%
%       Many names share one layout -- eleven 32-channel product codes are the
%       same arrangement. "cap_montage layouts" shows the groups.
%
%       Some montages stack non-scalp channels (EOG/EMG/ECG) on a placeholder
%       position; every sleep cap does. Those channels are kept, because
%       dropping them would break the match between eloc(k) and the data's
%       channel k, but a warning fires and they must be excluded before
%       interpolating a scalp map -- see the README.
%
%   Example:
%       eloc = cap_montage('BC-SL-64');
%       scalp = eloc(~ismember({eloc.labels}, {'EOG1','EOG2','EMG1','EMG2','EMG3','ECG'}));
%       channelbrowse(mdata, 'eloc', scalp);
%
%   See also: build_montage_library, read_montage, default_montage, plot_topo
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%************************************************************
%                      INPUT HANDLING
%************************************************************

if nargin < 1
    name = 'list';
end
validateattributes(name, {'char','string'}, {'scalartext'}, 'cap_montage', 'name');
name = char(name);

lib = library();

%************************************************************
%                     LIST / LAYOUTS
%************************************************************

if strcmpi(name, 'list')
    if nargout > 0
        varargout{1} = lib.names;
    else
        print_list(lib);
        varargout = {};
    end
    return;
end

if strcmpi(name, 'layouts')
    print_layouts(lib);
    varargout = {};
    return;
end

%************************************************************
%                        RESOLVE
%************************************************************

hit = find(strcmpi(lib.names, name), 1);

if isempty(hit)
    % No exact hit. A substring usually means a half-remembered name
    % ('BC-SL', 'sleep'), so show the candidates rather than only saying no.
    near = find(contains(lower(lib.names), lower(name)));
    if isscalar(near)
        hit = near;
    elseif ~isempty(near)
        error('cap_montage:ambiguousName', ...
            '''%s'' matches %d montages: %s\nUse the full name.', ...
            name, numel(near), strjoin(lib.names(near), ', '));
    else
        error('cap_montage:unknownName', ...
            'Unknown montage ''%s''.%s\nRun "cap_montage list" to see all %d.', ...
            name, suggestion_text(lib.names, name), numel(lib.names));
    end
end

L = lib.layouts(lib.map(hit));
varargout{1} = xyz_to_eloc(L.labels(:), L.xyz);

% Re-raise what the build recorded, so a caller sees the same caveat they would
% have seen reading the vendor file directly
if L.flags.coincident
    warning('cap_montage:coincidentElectrodes', ...
        ['%s: %d electrode pair(s) share a position -- non-scalp channels ' ...
         '(EOG/EMG/ECG) with placeholder coordinates. They are kept so channel ' ...
         'order still matches the data, but exclude them before interpolating a ' ...
         'scalp map, or it will average different values at one point.'], ...
        lib.names{hit}, L.flags.n_coincident);
end

end

%************************************************************
%                     LOAD THE LIBRARY
%************************************************************
function lib = library()
%LIBRARY  Load and cache the compiled montage store
%
%   Inputs:
%       none
%
%   Outputs:
%       lib : struct - the montage library (see build_montage_library)

% Cached across calls: channelbrowse resolves a montage on every launch, and
% re-reading the .mat each time is pure overhead. clear('cap_montage') to reload
% after a rebuild.
persistent cached
if ~isempty(cached)
    lib = cached;
    return;
end

here = fileparts(mfilename('fullpath'));
f = fullfile(here, 'channel_locs', 'montage_library.mat');

assert(exist(f,'file') == 2, 'cap_montage:noLibrary', ...
    ['Montage library not found at %s\n' ...
     'Build it with build_montage_library (needs the raw vendor files; see ' ...
     'channel_locs/vendor/README.md).'], f);

S = load(f, 'lib');
cached = S.lib;
lib = cached;

end

%************************************************************
%                       PRINT THE LIST
%************************************************************
function print_list(lib)
%PRINT_LIST  Show every montage, grouped by vendor
%
%   Inputs:
%       lib : struct - the montage library -- required
%
%   Outputs:
%       none (prints)

vend = arrayfun(@(k) lib.layouts(lib.map(k)).vendor, 1:numel(lib.names), ...
    'UniformOutput', false);
n = arrayfun(@(k) numel(lib.layouts(lib.map(k)).labels), 1:numel(lib.names));

fprintf('\n%d montages, %d distinct layouts (built %s)\n', ...
    numel(lib.names), numel(lib.layouts), lib.built);
fprintf('  cap_montage(''<name>'') to load one; cap_montage layouts to see equivalences\n');

uv = unique(vend);
for vv = 1:numel(uv)
    sel = find(strcmp(vend, uv{vv}));
    [~, o] = sort(n(sel));
    sel = sel(o);
    fprintf('\n  %s (%d)\n', uv{vv}, numel(sel));
    for ii = 1:numel(sel)
        if mod(ii-1,3) == 0, fprintf('    '); end
        fprintf('%-24s', sprintf('%s (%d)', lib.names{sel(ii)}, n(sel(ii))));
        if mod(ii,3) == 0 || ii == numel(sel), fprintf('\n'); end
    end
end
fprintf('\n');

end

%************************************************************
%                     PRINT THE LAYOUTS
%************************************************************
function print_layouts(lib)
%PRINT_LAYOUTS  Show distinct layouts and the names that map onto each
%
%   Inputs:
%       lib : struct - the montage library -- required
%
%   Outputs:
%       none (prints)

n = arrayfun(@(L) numel(L.labels), lib.layouts);
[~, order] = sort(n);

fprintf('\n%d distinct layouts across %d names:\n\n', ...
    numel(lib.layouts), numel(lib.names));
fprintf('  %-4s %-6s %s\n', 'N', 'vendor', 'names sharing this layout');
for jj = order(:)'
    mem = lib.names(lib.map == jj);
    fprintf('  %-4d %-6s %s\n', n(jj), lib.layouts(jj).vendor(1:min(6,end)), ...
        strjoin(mem, ', '));
end
fprintf('\n');

end

%************************************************************
%                   SUGGEST NEAR MISSES
%************************************************************
function txt = suggestion_text(names, name)
%SUGGESTION_TEXT  Build a "did you mean" fragment for an unknown name
%
%   Inputs:
%       names : 1xM cell - all known montage names -- required
%       name  : char - the name the caller asked for -- required
%
%   Outputs:
%       txt : char - a suggestion fragment, or '' if nothing is close

d = cellfun(@(s) edit_distance(lower(s), lower(name)), names);
[ds, order] = sort(d);
keep = order(ds <= max(3, ceil(numel(name)/3)));
if isempty(keep)
    txt = '';
else
    keep = keep(1:min(4,numel(keep)));
    txt = sprintf(' Did you mean: %s?', strjoin(names(keep), ', '));
end

end

%************************************************************
%                      EDIT DISTANCE
%************************************************************
function d = edit_distance(a, b)
%EDIT_DISTANCE  Levenshtein distance between two strings
%
%   Inputs:
%       a : char - first string -- required
%       b : char - second string -- required
%
%   Outputs:
%       d : double - number of single-character edits between a and b

na = numel(a); nb = numel(b);
D = zeros(na+1, nb+1);
D(:,1) = (0:na)';
D(1,:) = 0:nb;
for ii = 2:na+1
    for jj = 2:nb+1
        cost = ~isequal(a(ii-1), b(jj-1));
        D(ii,jj) = min([D(ii-1,jj)+1, D(ii,jj-1)+1, D(ii-1,jj-1)+cost]);
    end
end
d = D(end,end);

end
