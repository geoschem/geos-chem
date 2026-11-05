#!/usr/bin/env python3

# Author: Dandan Zhang (Harvard)

import sys
import re

if len(sys.argv) < 3:
    print("Usage: expand_carbon_eqn.py carbon.eqn.template Njac", file=sys.stderr)
    sys.exit(1)

infile = sys.argv[1]
Njac = int(sys.argv[2])

with open(infile, "r") as f:
    lines = f.readlines()

out_lines = []

in_defvar = False
in_equations = False

for line in lines:
    stripped = line.lstrip()
    tokens = stripped.split()

    # Section markers
    if stripped.startswith("#DEFVAR"):
        in_defvar = True
        in_equations = False
        out_lines.append(line)
        continue

    if stripped.startswith("#DEFFIX"):
        in_defvar = False
        in_equations = False
        out_lines.append(line)
        continue

    if stripped.startswith("#EQUATIONS"):
        in_defvar = False
        in_equations = True
        out_lines.append(line)
        continue

    # Always keep the original line
    out_lines.append(line)

    # =====================================================
    # 1) DEFVAR expansions: CH4_jacXXXX and L* dummy species
    # =====================================================
    if in_defvar and tokens:
        key = tokens[0]

        # CH4 main species -> CH4_jac000N
        if key == "CH4" and len(tokens) > 1 and tokens[1] == "=":
            for k in range(1, Njac + 1):
                idx = f"{k:04d}"
                out_lines.append(
                    f"CH4_jac{idx}   = IGNORE;  {{ Active methane Jacobian tracer }}\n"
                )

    # =====================================================
    # 2) EQUATIONS expansions: CH4 reactions -> CH4_jac000N
    # =====================================================
    if in_equations and tokens and tokens[0] == "CH4":

        eq_pos = line.find("=")
        if eq_pos == -1:
            continue  # not an equation line we care about

        for k in range(1, Njac + 1):
            idx = f"{k:04d}"
            new_line = line

            # Replace CH4 with CH4_jacXXXX as a whole word (avoid touching LCH4...)
            new_line = re.sub(r"\bCH4\b", f"CH4_jac{idx}", new_line)

            out_lines.append(new_line)

# Write result to stdout
sys.stdout.writelines(out_lines)
