# montage_toolbox

EEG electrode montages and scalp topography rendering for MATLAB. Reads vendor
electrode files into EEGLAB-style `eloc` structs, ships a compiled library of 109 caps,
and renders head models and interpolated scalp maps that match MNE-Python's
`plot_topomap`.

**Base MATLAB only** — no toolboxes, no EEGLAB, no other labcode submodules.

```matlab
eloc = cap_montage('AC-64');                    % actiCAP 64 from the library
plot_topo(values, eloc, 'style', 'map');        % interpolated scalp map
[~, ex, ey] = plot_topo([], eloc);              % head + electrodes, screen coords out
```

## Entry points

| Function | Purpose |
|---|---|
| `plot_topo` | Head model and interpolated scalp map; returns electrode screen coordinates |
| `cap_montage` | Look up a cap montage by vendor product code; `cap_montage list` to enumerate |
| `read_montage` | Read a `.sfp` / `.bvef` / `.elc` file into an `eloc` struct |
| `default_montage` | Pick a sensible montage from a channel count |
| `xyz_to_eloc` | Cartesian triple to a full `eloc` struct (polar + spherical fields) |
| `build_montage_library` | Recompile `channel_locs/montage_library.mat` from raw vendor files (rarely run) |

An `eloc` is an EEGLAB-style struct array. `plot_topo` needs either the polar pair
(`theta`, `radius`) or the cartesian triple (`X`, `Y`, `Z`); polar is used when both are
present. An optional `labels` field (char or numeric) drives the label display modes.

## The montage library — `cap_montage.m`

Pick a cap **by name**:

```matlab
cap_montage list                    % everything in the library
cap_montage layouts                 % grouped by layout, showing equivalent names
eloc = cap_montage('AC-64');        % actiCAP 64
eloc = cap_montage('BC-SL-64');     % BrainCap Sleep 64
```

**109 montages, 6–257 channels** — the full Brain Products catalogue (actiCAP, BrainCap, BrainCap Sleep/MR/MEG/TMS, LiveCap, R-Net, Xpress Twist) plus the EGI HydroCel nets. All 109 load clean, and Cz lands at radius 0.000 in all 106 that have it.

Names are the vendors' product codes: `AC-64` is actiCAP 64, `BC-SL-64` is BrainCap Sleep 64, and a `_REF` / `_NO_REF` suffix marks variants with and without the reference electrode. Lookup is case-insensitive; an unknown name suggests near misses (`AC64` → *did you mean AC-64, AP-64, BC-64?*) and a partial name matching several lists them.

### How it's stored

Everything lives in `channel_locs/montage_library.mat` — **77 KB**, built from the vendors' own files by `build_montage_library`. The raw files (1.2 MB of XML) are **not needed at run time and not distributed**; the library is self-contained and verified to reproduce `read_montage` bit-exactly for all 109 montages.

The store is deduplicated: **109 names → 65 distinct layouts**, since vendors ship one arrangement under many product codes (eleven 32-channel codes are the same layout). Each distinct layout is stored once as an ordered label list plus unit-sphere coordinates; `theta`/`radius` and the spherical fields are derived on load by `xyz_to_eloc`, so the library cannot drift from the projection `read_montage` uses.

Deduplication keys on the channel list **and** the coordinates, and both halves matter. Order matters because `eloc(k)` must line up with the data's channel *k*. Geometry matters because an identical channel list is not always the same cap — `AP-32` and `AS-32_NO_REF` list the same channels in the same order yet place `FC5` 1.5e-02 apart. For the same reason positions are not factored into a shared name→position table, tempting as that is: 19 scalp names (`FC5`, `FCz`, the `F11`/`FT11`/`TP11` series) genuinely sit in different places on different cap families, and `FC5` alone has three distinct positions across 90 caps.

To add caps: drop vendor files into `channel_locs/vendor/` (see its README) and re-run `build_montage_library`.

`read_montage` is the lower level, for a file outside the library:

```matlab
eloc = read_montage('/path/to/my_cap.bvef');   % .bvef / .sfp / .elc
```

### Non-scalp channels are stacked — exclude them from maps

**25 of the 109 montages contain electrodes at identical positions**, including *every* sleep cap. Vendors give non-scalp channels a placeholder angle rather than a real position: on a BrainVision sleep cap, `EMG1`, `EMG2`, `EMG3` and `ECG` all sit exactly on top of `Fpz`; several BrainCaps stack `IO` or `ECG` there too.

They are deliberately **not** dropped — they are recorded channels, so removing them would break the correspondence between `eloc(k)` and the data's channel *k*. `read_montage` warns (`read_montage:coincidentElectrodes`) instead.

This matters because it fails silently: interpolating a map over stacked points does not error, it averages several contradictory values at one location and returns a plausible-looking wrong answer. Drop them before mapping:

```matlab
eloc = cap_montage('BC-SL-64');
scalp = eloc(~ismember({eloc.labels}, {'EOG1','EOG2','EMG1','EMG2','EMG3','ECG'}));
plot_topo(vals(1:numel(scalp)), scalp, 'style', 'map');
```

Supported: `.sfp` (EGI geodesic nets), `.bvef` (BrainVision / EasyCap), `.elc` (ASA / ANT Neuro). Vendors disagree on axis order and units, so the loader resolves orientation *from the data* — it locates 10-20 landmarks and fiducials by label and picks the convention that puts them where anatomy says they go, rather than trusting an assumed axis order. Units are irrelevant (the projection depends only on angles). See `channel_locs/vendor/README.md` for where to obtain each vendor's files and why they are not committed.

The standard 10-20 and the AASM clinical six are subsets of the 21-channel actiCAP:

```matlab
eloc = read_montage('channel_locs/vendor/brainvision/AC-21.bvef');
six  = eloc(ismember({eloc.labels}, {'F3','F4','C3','C4','O1','O2'}));
```

### Migrating off the old `eloc*.mat` built-ins

The legacy `eloc6.mat`, `eloc64.mat`, `eloc64_2d.mat` and `eloc_10_20.mat` files have been **removed**. They were digitized off pictures, their channels were numbered rather than named, and their geometry did not describe a head — `eloc64` in particular was a radial starburst whose electrode radii varied 30-fold and whose best-fit sphere centre sat 50% of a radius off the origin. Every channel count they covered is now served by a real vendor layout from the library.

**The caveat that matters for old data:** those montages numbered their channels rather than naming them, so nothing recorded which physical electrode channel *k* was. The 64-channel default is now `AC-64`, which assumes the data's channel order matches actiCAP's (Fp1, Fp2, F7, F3, Fz, ...). For data that really came from a BrainVision 64 that is right; for anything else it silently mislabels channels. Verify the order against the recording header, and if the data is not a BrainVision 64, pass an explicit `eloc` from `read_montage` or `cap_montage`.
## Scalp rendering — `plot_topo.m`

`plot_topo` reads an `eloc` struct, projects the electrodes onto a 2-D head model, and optionally renders an interpolated scalp map.

```matlab
[hax, ex, ey] = plot_topo(values, eloc, 'style', 'map', 'conv', 'on');
```

- `values` is one number per channel, or `[]` for a head-and-electrodes plot with no map.
- `ex`/`ey` are the electrode screen coordinates in channel order — read them from the outputs rather than scraping `XData` back off the plotted lines. `channelbrowse` uses these for nearest-electrode hover lookup.

Coordinate convention follows EEGLAB, so existing montages plot unchanged: `theta` is in degrees, `radius` is a normalized arc length where **0.5 is the head rim** (the ears/eyes equator), and the nose points toward +y on screen. Channel *k* lands at

```
ex(k) = radius(k) * sin(theta(k) * pi/180)
ey(k) = radius(k) * cos(theta(k) * pi/180)
```

Electrode coordinates are used exactly as projected — nothing is rescaled to pull the outermost ring onto the rim. Radius 0.5 already *is* the rim, so Fp1 lands on the circle, and below-equator channels (`radius > 0.5`, e.g. FT9/TP9/PO9 at ~0.63 on a 64-channel cap) sit outside it in the conventional skirt. This is what MNE does: it keeps the head outline at the sphere radius and widens the clip circle instead. EEGLAB reaches the same electrode-to-head geometry from the other direction, shrinking electrodes *and* head cartoon together by `rmax/plotrad`.

### Interpolation

The scalp map reproduces MNE-Python's `plot_topomap` rather than EEGLAB's `topoplot`:

| | `plot_topo` | MNE `plot_topomap` | EEGLAB `topoplot` |
|---|---|---|---|
| projection | azimuthal equidistant | azimuthal equidistant | azimuthal equidistant |
| interpolant | Clough-Tocher cubic | Clough-Tocher cubic | `griddata 'v4'` biharmonic |
| edge behavior | synthetic ring, `border` mean | synthetic ring, `border` mean | free spline extrapolation |
| clip radius | `max(rim, 1.01 × outermost)` | `max(rim, 1.01 × outermost)` | `max(rim, 1.02 × outermost)` |
| grid | 200 | 64 | 67 |

MNE never lets the interpolant extrapolate. It rings the electrodes with synthetic points well outside the painted disc, gives each one the mean of its neighbouring electrodes (`'border','mean'`), and interpolates inside that ring — so every visible pixel is interpolated, not extrapolated. `plot_topo` does the same.

`'interp','cubic'` is a port of scipy's `CloughTocher2DInterpolator`, **not** `griddata`'s own `'cubic'`. MATLAB's is a different triangulation-based cubic scheme and lands several percent of the data range away from what MNE draws. Validated against scipy directly: identical to machine precision (~4e-15) on a generic point set, and on a real 64-channel montage the finished map agrees with MNE's to a maximum of 0.19% of the data range inside the head rim (rms 0.03%, r = 0.999999). The residual comes from the synthetic ring being exactly cocircular, a degenerate Delaunay configuration that MATLAB and Qhull break differently; it is confined to the skirt annulus outside the outermost electrodes.

`'cliprad'` is the one deliberate departure from MNE, and it is off by default. Set it to a radius to cap the painted map — `'cliprad', 0.5` confines the map to the head rim, the way this function behaved before it was matched to MNE. It is applied to the mask after interpolation and never moves the grid or the boundary ring, so values inside the cap are bit-identical to the uncapped map; it can only narrow the map, never widen it.

`'extrapolate'` selects MNE's boundary mode — `'head'` (default), `'local'`, or `'box'`. `'conv','on'` is kept as a synonym for `'extrapolate','local'`, which pushes the electrode convex hull out by one electrode spacing and masks the map to it. Set `'interp','v4'` for the EEGLAB/FieldTrip biharmonic spline instead; in that mode no synthetic ring is added, since the biharmonic spline extrapolates by construction, which is precisely how those toolboxes use it.

Note that the display interpolation is a separate question from *data* interpolation — repairing a bad channel's time series, where the field standard is unambiguously the Perrin et al. (1989) spherical splines used by MNE's `interpolate_bads`, EEGLAB's `eeg_interp`, and FieldTrip's `ft_channelrepair`. `plot_topo` does display only.

The head model (outline, nose, ears) uses MNE's geometry, scaled by `headrad`.

## Credits and prior work

`plot_topo`'s scalp map is a deliberate reimplementation of **MNE-Python's**
[`mne.viz.plot_topomap`](https://mne.tools/stable/generated/mne.viz.plot_topomap.html).
The projection, the synthetic boundary ring and its `'mean'` border condition, the clip
radius, the extrapolation modes, and the head/nose/ear outline geometry all follow MNE's
`mne/viz/topomap.py`. Where this toolbox departs from MNE it says so — currently only the
opt-in `'cliprad'` option and the denser default grid.

The Clough-Tocher interpolant in `plot_topo`'s local functions `clough_tocher_2d` and
`estimate_gradients_2d` is a MATLAB port of **SciPy's**
[`CloughTocher2DInterpolator`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CloughTocher2DInterpolator.html)
(`scipy/interpolate/interpnd.pyx`) — the curvature-minimizing gradient network and the
Bézier control net including SciPy's affine-invariant choice of cross-boundary direction.
MATLAB's own `griddata 'cubic'` is a different scheme and does not reproduce MNE's output.

MNE-Python and SciPy are both BSD-3-Clause, as is this toolbox. Neither project endorses
this port, and any bug here is this toolbox's, not theirs. If you use the scalp maps in
published work, please cite MNE and SciPy alongside this repository.

### References

- Gramfort A, Luessi M, Larson E, et al. (2013). MEG and EEG data analysis with MNE-Python. *Frontiers in Neuroscience* 7:267. doi:10.3389/fnins.2013.00267
- Virtanen P, Gommers R, Oliphant TE, et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17:261–272. doi:10.1038/s41592-019-0686-2
- Nielson GM (1983). A method for interpolating scattered data based upon a minimum norm network. *Mathematics of Computation* 40:253–271. doi:10.1090/S0025-5718-1983-0679444-7
- Renka RJ, Cline AK (1984). A triangle-based C¹ interpolation method. *Rocky Mountain Journal of Mathematics* 14:223–237. doi:10.1216/RMJ-1984-14-1-223
- Alfeld P (1984). A trivariate Clough–Tocher scheme for tetrahedral data. *Computer Aided Geometric Design* 1:169–181. doi:10.1016/0167-8396(84)90029-3
- Farin G (1986). Triangular Bernstein–Bézier patches. *Computer Aided Geometric Design* 3:83–127. doi:10.1016/0167-8396(86)90016-6
- Sandwell DT (1987). Biharmonic spline interpolation of GEOS-3 and SEASAT altimeter data. *Geophysical Research Letters* 14:139–142. doi:10.1029/GL014i002p00139 — the `'interp','v4'` mode
- Perrin F, Pernier J, Bertrand O, Echallier JF (1989). Spherical splines for scalp potential and current density mapping. *Electroencephalography and Clinical Neurophysiology* 72:184–187. doi:10.1016/0013-4694(89)90180-6 — the standard for *data* interpolation, which this toolbox does not do
- Delorme A, Makeig S (2004). EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics. *Journal of Neuroscience Methods* 134:9–21. doi:10.1016/j.jneumeth.2003.10.009 — the `eloc` struct convention and the `'v4'` topoplot lineage

### History

This code began inside
[`channelbrowse`](https://github.com/preraulab/channelbrowse), which vendored
`topoplotFast.m` (a trimmed EEGLAB `topoplot`) plus eight EEGLAB-derived helpers.
Replacing them removed ~110 KB of dependency and resolved a license conflict: the EEGLAB
files are GPL, which sits badly in a BSD-3 repository. The montage and rendering code was
then split out here so `channelbrowse`, `eeg_analysis`, and anything else can share one
copy. File history is preserved from `channelbrowse`.

## Files

- `plot_topo.m` — head model and scalp-map plotter
- `cap_montage.m` — montage lookup by product code
- `read_montage.m` — vendor electrode file reader (`.sfp` / `.bvef` / `.elc`)
- `default_montage.m` — montage from a channel count
- `xyz_to_eloc.m` — cartesian triple to full `eloc`
- `build_montage_library.m` — recompile the library from raw vendor files
- `channel_locs/montage_library.mat` — the compiled store (77 KB, 109 caps)
- `channel_locs/vendor/` — where raw vendor files go; see its README

## License

BSD 3-Clause. See `LICENSE`.
