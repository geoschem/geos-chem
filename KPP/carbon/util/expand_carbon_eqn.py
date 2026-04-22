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

        # LCH4byOH -> LCH4byOH_jac000N
        elif key == "LCH4byOH":
            for k in range(1, Njac + 1):
                idx = f"{k:04d}"
                out_lines.append(
                    f"LCH4byOH_jac{idx} = IGNORE;  "
                    f"{{ Dummy spc to track loss of CH4_jac{idx} by OH     }}\n"
                )

        # LCH4byCl -> LCH4byCl_jac000N
        elif key == "LCH4byCl":
            for k in range(1, Njac + 1):
                idx = f"{k:04d}"
                out_lines.append(
                    f"LCH4byCl_jac{idx} = IGNORE;  "
                    f"{{ Dummy spc to track loss of CH4_jac{idx} by Cl     }}\n"
                )

        # LCH4inStrat -> LCH4inStrat_jac000N
        elif key == "LCH4inStrat":
            for k in range(1, Njac + 1):
                idx = f"{k:04d}"
                out_lines.append(
                    f"LCH4inStrat_jac{idx} = IGNORE;  "
                    f"{{ Dummy spc to track loss of CH4_jac{idx} in strat }}\n"
                )

    # =====================================================
    # 2) EQUATIONS expansions: CH4 reactions -> CH4_jac000N
    # =====================================================
    if in_equations and tokens and tokens[0] == "CH4":
        # We need to find the product species (right after "=")
        eq_pos = line.find("=")
        if eq_pos == -1:
            continue  # not an equation line we care about

        # substring after "=" up to ":" (if present)
        rest = line[eq_pos + 1 :]
        colon_pos = rest.find(":")
        if colon_pos != -1:
            lhs_products = rest[:colon_pos]
        else:
            lhs_products = rest

        lhs_products_stripped = lhs_products.strip()
        if not lhs_products_stripped:
            continue

        prod = lhs_products_stripped.split()[0]  # e.g. LCH4byOH, LCH4byCl, LCH4inStrat

        for k in range(1, Njac + 1):
            idx = f"{k:04d}"
            new_line = line

            # Replace CH4 with CH4_jacXXXX as a whole word (avoid touching LCH4...)
            new_line = re.sub(r"\bCH4\b", f"CH4_jac{idx}", new_line)

            # Replace the product species once with product_jacXXXX
            new_line = re.sub(
                rf"\b{re.escape(prod)}\b",
                f"{prod}_jac{idx}",
                new_line,
                count=1,
            )

            out_lines.append(new_line)

# Write result to stdout
sys.stdout.writelines(out_lines)
