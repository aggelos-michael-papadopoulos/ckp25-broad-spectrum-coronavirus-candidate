#!/usr/bin/env python3
"""
_mdrun_progress.py — tqdm progress bar for gmx mdrun by polling the log file.

Usage (called internally by run_simulation_gmx.sh):
    python _mdrun_progress.py <log_file> <total_steps> <label> [dt_ps]
"""
import sys
import re
import time
import os
from tqdm import tqdm

log_file    = sys.argv[1]            if len(sys.argv) > 1 else "md.log"
total_steps = int(sys.argv[2])       if len(sys.argv) > 2 else 50000
label       = sys.argv[3]            if len(sys.argv) > 3 else "MD"
dt_ps       = float(sys.argv[4])     if len(sys.argv) > 4 else 0.002

total_ns = total_steps * dt_ps / 1000.0
scale    = dt_ps / 1000.0   # steps → ns

bar = tqdm(total=round(total_ns, 4), desc=label, unit="ns",
           bar_format="{l_bar}{bar}| {n:.3f}/{total:.3f} ns  [{elapsed}<{remaining}, {rate_fmt}]",
           dynamic_ncols=True)

last_step = 0

# Wait for the log file to appear (mdrun creates it shortly after start)
for _ in range(60):
    if os.path.exists(log_file):
        break
    time.sleep(1)

# Poll the log file until the run completes
step_pattern = re.compile(r'^\s+Step\s+Time\s*$')
done_patterns = ("Finished mdrun", "gmx mdrun", "GROMACS reminds you")

with open(log_file, "r", errors="replace") as fh:
    want_next = False
    while True:
        line = fh.readline()
        if not line:
            # No new data — check if mdrun finished
            fh.seek(0, 2)       # seek to end to get fresh reads
            if last_step >= total_steps:
                break
            time.sleep(0.5)
            continue

        # MD format: "Step  Time" header, value on next line
        if step_pattern.match(line):
            want_next = True
            continue
        if want_next:
            want_next = False
            m = re.match(r'^\s+(\d+)\s+[\d.]+', line)
            if m:
                step = int(m.group(1))
                if step > last_step:
                    delta_ns = (step - last_step) * scale
                    bar.update(round(delta_ns, 6))
                    last_step = step
                if last_step >= total_steps:
                    break

        # EM format: "Step=    1, Dmax= ..."  or  "Stepsize too small..."
        m_em = re.match(r'^Step=\s*(\d+),', line)
        if m_em:
            step = int(m_em.group(1))
            if step > last_step:
                delta_ns = (step - last_step) * scale
                bar.update(round(delta_ns, 6))
                last_step = step
            if last_step >= total_steps:
                break

        # EM converged line signals completion
        if "Potential Energy  =" in line or "converged to Fmax" in line:
            break

# Fill any remaining gap (e.g. last batch of steps)
remaining = (total_steps - last_step) * scale
if remaining > 0:
    bar.update(round(remaining, 6))
bar.close()
