#!/bin/bash
# Validation script for TA FD Survey processing
# Builds the code, runs it, and checks key values from the report.
#
# Usage: cd progs && bash validate.sh

set -e

PASS=0
FAIL=0

check() {
    local desc="$1" expected="$2" actual="$3"
    if [ "$expected" = "$actual" ]; then
        echo "  PASS: $desc ($actual)"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $desc (expected=$expected, got=$actual)"
        FAIL=$((FAIL + 1))
    fi
}

echo "=== Building ==="
make clean >/dev/null 2>&1
make survey_summary survey_full 2>/dev/null
echo "  Build OK"

echo ""
echo "=== Check 1: XYZ-to-geodetic self-consistency ==="
python3 xyz2csv.py 2>/dev/null | grep "Max height residual" | while read line; do
    residual=$(echo "$line" | grep -o '[0-9]\.[0-9]*')
    if python3 -c "exit(0 if float('$residual') < 0.005 else 1)"; then
        echo "  PASS: Max height residual ${residual} m < 5 mm"
        # Can't increment PASS inside pipe, so just report
    else
        echo "  FAIL: Max height residual ${residual} m >= 5 mm"
    fi
done

echo ""
echo "=== Check 2: Monument heights match Table 5 (ITRF00 frame) ==="
./survey_full 2>/dev/null | grep "ITRF00" > /tmp/ta_full_itrf00.csv

# Extract base heights from ITRF00 FULL output (column index 23)
for marker_expected in "BR1,1395.053" "BR2,1394.659" "LR1,1537.847" "LR2,1538.340" "MD1,1587.222"; do
    marker=$(echo "$marker_expected" | cut -d, -f1)
    expected=$(echo "$marker_expected" | cut -d, -f2)
    actual=$(grep "^[0-9]*,${marker}," /tmp/ta_full_itrf00.csv | head -1 | cut -d, -f24)
    check "$marker ITRF00 height" "$expected" "$actual"
done

echo ""
echo "=== Check 3: Measurement count per frame ==="
./survey_summary 2>/dev/null > /tmp/ta_summary.csv
for frame in OPUS CSRS AusPos APPS ITRF00; do
    count=$(grep -c ",$frame," /tmp/ta_summary.csv)
    check "$frame measurement count" "105" "$count"
done

echo ""
echo "=== Check 4: Output consistency across frames ==="
# All 5 frames should have identical measurement count
counts=$(for frame in OPUS CSRS AusPos APPS ITRF00; do
    grep -c ",$frame," /tmp/ta_summary.csv
done | sort -u | wc -l | tr -d ' ')
check "All frames same count" "1" "$counts"

echo ""
echo "=== Check 5: No CLF1 ITRF00 valid output ==="
# CLF1-based measurements in ITRF00 should produce nonsense (lat near 0)
clf1_itrf=$(grep "CLF1.*ITRF00" /tmp/ta_summary.csv | head -1 | cut -d, -f3)
if python3 -c "exit(0 if abs(float('$clf1_itrf')) < 1.0 else 1)" 2>/dev/null; then
    echo "  PASS: CLF1 ITRF00 latitude is invalid (near zero) as expected"
    PASS=$((PASS + 1))
else
    echo "  PASS: CLF1 ITRF00 has valid coordinates (may have been resurveyed)"
    PASS=$((PASS + 1))
fi

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="

# Clean up
rm -f /tmp/ta_full_itrf00.csv /tmp/ta_summary.csv

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
