#!/usr/bin/env python3
"""
Convert ITRF00 ECEF XYZ coordinates from Table 5 of the TA FD Survey Marker
Reference Frame Correction paper to geodetic lat/lon/h for Survey_monuments.csv.

Uses WGS84 ellipsoid parameters matching the xyz2geo() function in survey_processor.c.
Inverse-variance weighted average heights from multiple OPUS solutions.
"""

import math

# WGS84 ellipsoid parameters
a = 6378137.0           # semi-major axis (m)
f = 1.0 / 298.257223563 # flattening
b = a * (1.0 - f)       # semi-minor axis
e2 = 2.0 * f - f * f    # first eccentricity squared
ep2 = e2 / (1.0 - e2)   # second eccentricity squared


def xyz2geo(x, y, z):
    """Convert ECEF XYZ to geodetic lat (rad), lon (rad), h (m).
    Bowring's iterative method, matching survey_processor.c xyz2geo()."""
    p = math.sqrt(x * x + y * y)
    lon = math.atan2(y, x)

    # Bowring's initial approximation
    theta = math.atan2(z * a, p * b)
    lat = math.atan2(
        z + ep2 * b * math.sin(theta) ** 3,
        p - e2 * a * math.cos(theta) ** 3
    )

    # One iteration for convergence
    for _ in range(10):
        sin_lat = math.sin(lat)
        N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
        lat_new = math.atan2(z + e2 * N * sin_lat, p)
        if abs(lat_new - lat) < 1e-14:
            break
        lat = lat_new

    sin_lat = math.sin(lat)
    N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    h = p / math.cos(lat) - N

    return lat, lon, h


def deg2dms(decimal_degrees):
    """Convert decimal degrees to (dd, mm, ss.sssss)."""
    d = abs(decimal_degrees)
    dd = int(d)
    d = (d - dd) * 60.0
    mm = int(d)
    ss = (d - mm) * 60.0
    return dd, mm, ss


# Table 5: Corrected ITRF00 coordinates (inverse-variance weighted averages)
# From TA_FD_Survey_Marker_Reference_Frame_Correction.tex, lines 308-314
markers = [
    ("BR1",  -1911712.757, -4567269.863,  4009427.954, 1395.053, 0.004),
    ("BR2",  -1911674.654, -4567248.389,  4009469.677, 1394.659, 0.009),
    ("CLF2", -1924405.289, -4553771.178,  4018597.785, 1369.979, 0.016),
    ("LR1",  -1943693.195, -4552383.580,  4011197.010, 1537.847, 0.009),
    ("LR2",  -1943699.086, -4552343.485,  4011240.156, 1538.340, 0.010),
    ("MD1",  -1926374.126, -4539528.317,  4033980.925, 1587.222, 0.004),
    ("MD2",  -1926193.139, -4539634.752,  4033940.294, 1582.449, 0.020),
]

print("# ITRF00 inverse-variance weighted averages from OPUS solutions")
print("# Converted from Table 5 of TA FD Survey Marker Reference Frame Correction")
print()
print("# Self-consistency check: computed h vs paper h")

max_h_residual = 0.0
for name, x, y, z, h_paper, sigma_h in markers:
    lat_rad, lon_rad, h_computed = xyz2geo(x, y, z)
    residual = abs(h_computed - h_paper)
    max_h_residual = max(max_h_residual, residual)
    print(f"#   {name}: h_computed={h_computed:.3f}  h_paper={h_paper:.3f}  "
          f"residual={residual:.4f} m")

print(f"# Max height residual: {max_h_residual:.4f} m")
if max_h_residual > 0.050:
    print("# WARNING: height residual exceeds 50 mm -- check input coordinates!")
print()

# Generate CSV lines
print("# CSV lines for Survey_monuments.csv:")
print()
for name, x, y, z, h_paper, sigma_h in markers:
    lat_rad, lon_rad, h_computed = xyz2geo(x, y, z)

    lat_deg = math.degrees(lat_rad)
    lon_deg = math.degrees(lon_rad)  # negative for west

    lat_dd, lat_mm, lat_ss = deg2dms(lat_deg)
    lon_dd, lon_mm, lon_ss = deg2dms(lon_deg)  # positive output (west convention)

    # Use paper's h (weighted average), not computed h (which may differ at mm level
    # due to rounding of XYZ to mm precision in the paper)
    h = h_paper

    # Lat/lon uncertainties: 0.010 m (conservative; not used by survey_processor.c)
    sigma_lat = 0.010
    sigma_lon = 0.010

    print(f"{name},ITRF00,{lat_dd},{lat_mm},{lat_ss:.5f},{sigma_lat},{lon_dd},{lon_mm},{lon_ss:.5f},{sigma_lon},{h},{sigma_h}")
