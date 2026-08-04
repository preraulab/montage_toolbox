function [eloc, source] = default_montage(n_chans)
%DEFAULT_MONTAGE  Pick the built-in sensor montage for a given channel count
%
%   Usage:
%       eloc = default_montage(64)
%       [eloc, source] = default_montage(6)
%
%   Inputs:
%       n_chans : integer - number of channels in the data -- required
%
%   Outputs:
%       eloc   : 1xC struct - channel locations, with real 10-20 labels when
%                the montage came from a vendor file
%       source : char - description of where the montage came from, for
%                titles and error messages
%
%   Notes:
%       Montages come from the compiled library (channel_locs/montage_library.mat)
%       via cap_montage, and are Brain Products actiCAP layouts: named electrodes
%       (Cz, F3) with correct geometry. The legacy eloc*.mat files are the
%       fallback and are known to be distorted -- see the README -- so falling
%       back warns rather than passing silently.
%
%       Any cap in the library can be used instead by name:
%       channelbrowse(mdata, 'eloc', cap_montage('BC-SL-64')).
%
%       6 and 19 channels resolve to the clinical standards -- 'AASM-6'
%       (F3 F4 C3 C4 O1 O2) and '10-20' -- both of which are also addressable
%       by name through cap_montage.
%
%   Example:
%       [eloc, src] = default_montage(64);
%       plot_topo([], eloc, 'electrodes', 'ptslabels');
%
%   See also: cap_montage, read_montage, plot_topo, channelbrowse
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%************************************************************
%                      INPUT HANDLING
%************************************************************

validateattributes(n_chans, {'numeric'}, {'scalar','integer','positive'}, ...
    'default_montage', 'n_chans');

%************************************************************
%                    RESOLVE THE MONTAGE
%************************************************************

switch n_chans
    case 6
        [eloc, source] = whole_cap('AASM-6', 'AASM clinical six');
    case 19
        [eloc, source] = whole_cap('10-20', 'standard 10-20');
    case 21
        [eloc, source] = whole_cap('AC-21', 'actiCAP 21 (10-20 + mastoids)');
    case 32
        [eloc, source] = whole_cap('AC-32', 'actiCAP 32');
    case 64
        [eloc, source] = whole_cap('AC-64', 'actiCAP 64');
    case 96
        [eloc, source] = whole_cap('AC-96', 'actiCAP 96');
    case 128
        [eloc, source] = whole_cap('AC-128', 'actiCAP 128');
    otherwise
        % No cap matches. Fall back to the 64 rather than guessing, and say so:
        % the electrode count will not line up with the data, and the caller's
        % own channel-count assert is what should stop the run.
        warning('default_montage:noMatch', ...
            ['No built-in montage for %d channels; falling back to the 64-channel ' ...
             'cap. Pass an explicit ''eloc'' -- cap_montage(''<name>'') will load ' ...
             'any cap in the library.'], n_chans);
        [eloc, source] = whole_cap('AC-64', 'actiCAP 64');
end

if isempty(eloc)
    [eloc, source] = legacy_fallback(n_chans);
end

end

%************************************************************
%                     LOAD A WHOLE CAP
%************************************************************
function [eloc, source] = whole_cap(cap, name)
%WHOLE_CAP  Load a library cap in full
%
%   Inputs:
%       cap  : char - the montage name in the library (e.g. 'AC-64') -- required
%       name : char - human-readable montage name -- required
%
%   Outputs:
%       eloc   : 1xC struct - channel locations, or [] if the library is absent
%       source : char - montage name, or '' if the library is absent

eloc = [];
source = '';
try
    eloc = cap_montage(cap);
    source = name;
catch ME
    % A missing library is the expected case on a fresh clone and the caller
    % falls back; anything else is a real fault and should not be swallowed.
    if ~strcmp(ME.identifier, 'cap_montage:noLibrary')
        rethrow(ME);
    end
end

end

%************************************************************
%                LEGACY .MAT FALLBACK
%************************************************************
function [eloc, source] = legacy_fallback(n_chans)
%LEGACY_FALLBACK  Load a bundled eloc*.mat when no vendor file is available
%
%   Inputs:
%       n_chans : integer - number of channels in the data -- required
%
%   Outputs:
%       eloc   : 1xC struct - channel locations
%       source : char - montage name

switch n_chans
    case 6
        file = 'eloc6';
    otherwise
        file = 'eloc64';
end

warning('default_montage:legacyMontage', ...
    ['Montage library not found; falling back to the bundled %s.mat, whose ' ...
     'geometry is known to be distorted and whose channels are numbered rather ' ...
     'than named. Build the library with build_montage_library, or pass an ' ...
     'explicit ''eloc''.'], file);

montage = load(file, 'eloc');
eloc = montage.eloc;
source = sprintf('%s.mat (legacy)', file);

end
