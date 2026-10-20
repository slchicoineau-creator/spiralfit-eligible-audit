#!/bin/bash
# Centre search (pre-registered 2026-09-30, README §4.5): 8 fits in 2 lanes via tol_curve.sh's `one`, then the
# pre-registered analysis. Needs the §4.2 fits (tol_curve.sh) and the §4.1 fits (umb_sens.sh) for the other probes.
W="$(cd "$(dirname "$0")" && pwd)"; P="$(dirname "$W")"
O=${R3:-$HOME/r3}/fits/tol; mkdir -p $O
one() { bash $W/tol_curve.sh one "$@"; }
lane1() { one cs_xm250_s2 -250 0 2; one cs_xm100 -100 0 1; one cs_ym100 0 -100 1; one cs_ym500 0 -500 1; }
lane2() { one cs_ym250 0 -250 1; one cs_xm500 -500 0 1; one cs_xm250_ym250 -250 -250 1; one cs_xm250_yp250 -250 250 1; }
lane1 > $O/cs_lane1.log 2>&1 &
sleep 120
lane2 > $O/cs_lane2.log 2>&1 &
wait
cd $P && ${PY:-python} scripts/centre_search.py ${R3:-$HOME/r3}/fits
