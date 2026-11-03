function scalp = scalp_channels(eloc, nonscalp)
%SCALP_CHANNELS  Which channels of an eloc are EEG: true unless the label marks a non-EEG channel
%
%   Usage:
%       scalp = scalp_channels(eloc)
%       scalp = scalp_channels(eloc, {'Chin'})
%       scalp = scalp_channels(labels)
%
%   Inputs:
%       eloc     : 1xC struct with a labels field, or a cellstr of C labels
%                  -- required
%       nonscalp : cellstr - further labels to treat as non-EEG, matched
%                  without regard to case (default: {})
%
%   Outputs:
%       scalp : Cx1 logical - true where the channel is EEG
%
%   Notes:
%       The rule is by name because montages carry no channel type: a label
%       of EOG, EMG, ECG or EKG, with or without a leading H, V, L or R and
%       a trailing number, or IO (infraorbital), is non-EEG. The library's
%       110 montages use ECG, EMG1-3, EOG1-2, IO and VEOG, all of which this
%       catches; a scalp label (Iz, O1, FC5, ...) never matches. plot_topo
%       and electrode_distances apply this one rule, so a channel left out
%       of the map is also left out of the distances.
%
%   Example:
%       eloc = cap_montage('BC-SL-64');
%       nnz(scalp_channels(eloc))                 % 62: EOG, EMG and ECG are out
%
%   See also: plot_topo, electrode_distances, cap_montage
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

if nargin < 2 || isempty(nonscalp), nonscalp = {}; end
if isstruct(eloc), labels = {eloc.labels}; else, labels = cellstr(eloc); end
non_eeg = '^(?:[HVLR]?(?:EOG|EMG|ECG|EKG)|IO)\d*$';
scalp = cellfun(@isempty, regexpi(labels(:), non_eeg, 'once')) & ~ismember(lower(labels(:)), lower(cellstr(nonscalp)));
end
