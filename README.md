# ZStack GUI

MATLAB UI for browsing ND2 stacks with Bio-Formats, interactive ROI selection, frame navigation, and ROI-only histograms.

## Requirements
- MATLAB with `uifigure` (R2020a+ recommended).
- Bio-Formats `bfmatlab` toolbox: download from https://www.openmicroscopy.org/bio-formats/downloads/ and place it inside `bfmatlab/` so it contains `bfopen.m`, `bfGetReader.m`, `bfCheckJavaPath.m`, etc.

## Quick start
1. Place Bio-Formats in `bfmatlab/` (see above).
2. In MATLAB, `cd` into this folder and run:
   ```matlab
   ZStackGUI
   ```
3. Use the right panel:
   - `Load ND2` / `Load folder` to import files (folder load is parallel when the Parallel Toolbox is available).
   - Selector: dropdown, Prev/Next buttons, frame slider, and `Save`.
   - Define ROI: draw rectangle or edit X/Y/Width/Height, set histogram bins.

## Features
- Left view shows the current Z slice; right view shows the intensity histogram restricted to the ROI with adjustable bins.
- Prev/Next video buttons, frame slider, and dropdown for multi-video browsing.
- ROI is editable (drag/resize) and fields stay synced.
- Save button crops all frames of all loaded videos to the current ROI and writes a single `roiData` struct array (includes stack name/path, ROI, frames, bins, per-frame ΔT/Z when present, acquisition date).

## Notes
- Default start directory for file dialogs: `\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS` (falls back to `pwd` if unavailable).
- Histogram requires an active ROI; otherwise a prompt is shown.
