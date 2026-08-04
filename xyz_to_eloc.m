function eloc = xyz_to_eloc(labels, xyz)
%XYZ_TO_ELOC  Build an eloc struct from ALS cartesian electrode coordinates
%
%   Usage:
%       eloc = xyz_to_eloc(labels, xyz)
%
%   Inputs:
%       labels : Cx1 cell - electrode label strings -- required
%       xyz    : Cx3 double - coordinates in the ALS convention: +x anterior,
%                +y left, +z up. Scale is irrelevant, only direction matters
%                -- required
%
%   Outputs:
%       eloc : 1xC struct - channel locations with fields labels, X, Y, Z,
%              theta, radius, sph_theta, sph_phi, sph_radius
%
%   Notes:
%       The azimuthal equidistant projection used throughout this toolbox, in
%       one place. Both read_montage (parsing a vendor file) and cap_montage
%       (reading the built library) route through here, so a montage cannot
%       project differently depending on which door it came in by.
%
%       radius is linear in the polar angle down from the vertex, so the
%       equator -- the ears/eyes line -- lands on the head rim at 0.5, and
%       electrodes below it exceed 0.5 and plot outside the drawn head.
%
%   Example:
%       eloc = xyz_to_eloc({'Cz'}, [0 0 1]);   % eloc.radius == 0
%
%   See also: read_montage, cap_montage, plot_topo
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

validateattributes(labels, {'cell'}, {'vector'}, 'xyz_to_eloc', 'labels');
validateattributes(xyz, {'numeric'}, {'2d','ncols',3}, 'xyz_to_eloc', 'xyz');
assert(numel(labels) == size(xyz,1), 'xyz_to_eloc:sizeMismatch', ...
    '%d labels but %d coordinate rows.', numel(labels), size(xyz,1));

X = xyz(:,1); Y = xyz(:,2); Z = xyz(:,3);

[az, elev, r] = cart2sph(X, Y, Z);

sph_theta = az * 180/pi;
sph_phi   = elev * 180/pi;
theta     = -sph_theta;
radius    = 0.5 - elev/pi;

eloc = struct('labels', labels(:)', ...
    'X', num2cell(X'), 'Y', num2cell(Y'), 'Z', num2cell(Z'), ...
    'theta', num2cell(theta'), 'radius', num2cell(radius'), ...
    'sph_theta', num2cell(sph_theta'), 'sph_phi', num2cell(sph_phi'), ...
    'sph_radius', num2cell(r'));

end
