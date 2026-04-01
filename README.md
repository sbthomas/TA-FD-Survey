# Telescope Array FD Survey: ITRF00 Reference Frame Corrections

GPS survey data, total station measurements, and processing code for the
Telescope Array fluorescence detector (FD) telescope mirror positions at
Black Rock Mesa, Long Ridge, and Middle Drum.

This repository accompanies the technical note:

> S.B. Thomas, "Reference Frame Discrepancy in TA Fluorescence Detector
> Survey Marker Coordinates: Impact on Opt-copter Calibration and
> Implications for Xmax Reconstruction," TA Internal Technical Note,
> April 2026 (revised).

## The Problem

The eight GPS survey markers used for FD telescope geometry were stored in
the analysis database using **NAD 83** ellipsoidal heights. The Opt-copter
RTK-GPS reports positions in **WGS 84 / ITRF**. The height difference
between these frames at the TA site is **~0.73 m**, producing a systematic
**+0.17 deg** elevation angle bias that matches the offsets observed by
Matsuzawa (2025) at Middle Drum and Nakazawa (2020) at Black Rock Mesa.

## Repository Contents

```
.
├── BR1/ BR2/ CLF1/ CLF2/        GPS RINEX observations (.11o, .11n)
│   LR1/ LR2/ MD1/ MD2/          per monument, plus BR1 OPUS report
├── drawings/                     FD station site drawings
├── progs/
│   ├── survey_processor.c        Survey reduction code (C)
│   ├── survey_averager.c         Multi-processor coordinate averager (C)
│   ├── xyz2csv.py                ECEF XYZ to geodetic converter (Python 3)
│   ├── Survey_monuments.csv      Monument coordinates in 5 reference frames
│   ├── Survey_Record.csv         Total station field measurements
│   └── ELS.txt                   ELS reference station data
├── Survey_Record.csv             Total station measurements (master copy)
├── Survey_Monuments.txt          Monument coordinate summary (wiki format)
├── antenna_heights.txt           GPS antenna heights above mark
├── TA_FD_Survey_Marker_Reference_Frame_Correction.tex   Technical note (LaTeX)
└── TA.pdf                        Technical note (compiled PDF)
```

## Reproducing the Analysis

### Requirements

- A C compiler (gcc, clang, or cc)
- Python 3 (standard library only, no external packages)

### Build

```bash
cd progs

# Summary output (one line per measurement per reference frame):
cc -DSUMMARY -o survey_summary survey_processor.c -lm

# Averaged mirror positions per station (wiki table format):
cc -DDISTANCES -o survey_distances survey_processor.c -lm

# Full detail (all measurement parameters and computed positions):
cc -DFULL -o survey_full survey_processor.c -lm
```

### Run

All programs read `Survey_Record.csv` and `Survey_monuments.csv` from the
current working directory and write to stdout.

```bash
cd progs

# Produce corrected telescope positions in all 5 reference frames:
./survey_summary

# Produce averaged mirror position tables (as in the report):
./survey_distances

# Verify the XYZ-to-geodetic conversion of Table 5 coordinates:
python3 xyz2csv.py
```

### Interpreting the Output

The processing code runs all five reference frames in sequence:

| Index | Frame   | Source |
|-------|---------|--------|
| 0     | OPUS    | Single OPUS ITRF00 solution per marker |
| 1     | CSRS    | NRCan CSRS-PPP (ITRF05) |
| 2     | AusPos  | Geoscience Australia (ITRF2005) |
| 3     | APPS    | JPL APPS (WGS 84) |
| 4     | ITRF00  | Inverse-variance weighted average of multiple OPUS solutions |

**The ITRF00 frame (index 4) contains the corrected coordinates** that
should replace the NAD 83 values in the analysis database. The other four
frames are retained for cross-validation.

### Validation Checks

After building and running, verify:

1. **Monument heights match Table 5 of the report:**
   The ITRF00 base heights in the FULL output (column 23) should be:
   BR1=1395.053, BR2=1394.659, LR1=1537.847, LR2=1538.340, MD1=1587.222 m.

2. **Inter-mirror distances are frame-independent:**
   Distances between adjacent mirrors should agree across all five frames
   to within 1 mm.

3. **XYZ conversion self-consistency:**
   Running `python3 xyz2csv.py` should show height residuals < 3 mm for
   all seven markers.

4. **CLF1 has no ITRF00 solution:**
   CLF1-based measurements produce invalid output in the ITRF00 frame.
   This is expected (the OPUS solution failed; CLF1 requires resurvey).

## Key Results

| Marker | ITRF00 h (m) | NAD 83 h (m) | Delta h (m) |
|--------|-------------|-------------|------------|
| BR1    | 1395.053    | 1395.791    | -0.738     |
| BR2    | 1394.659    | 1395.405    | -0.746     |
| CLF2   | 1369.979    | 1370.716    | -0.737     |
| LR1    | 1537.847    | 1538.577    | -0.730     |
| LR2    | 1538.340    | 1539.073    | -0.733     |
| MD1    | 1587.222    | 1587.964    | -0.742     |
| MD2    | 1582.449    | 1583.179    | -0.730     |

Weighted mean offset: **-0.736 m** (NAD 83 heights are systematically
higher than ITRF).

## Data Provenance

- **GPS observations:** March 2011 campaign using Trimble TRM41249.00 antenna.
- **OPUS processing:** NGS Online Positioning User Service, April 2011.
- **Total station measurements:** March-April 2011, recorded in `Survey_Record.csv`.
- **ITRF00 weighted averages:** Computed from multiple OPUS static solutions
  (n = 1-7 per marker), inverse-variance weighted for height, simple averaged
  for XYZ. See Table 1 of the report for per-marker solution counts.
- **CLF1:** OPUS solution failed (3% observations used, 130 m peak-to-peak
  vertical uncertainty). Resurvey required.

## License

This data is provided for use by the Telescope Array collaboration.
Contact S.B. Thomas (ORCID: 0000-0002-8828-7856) for questions.
