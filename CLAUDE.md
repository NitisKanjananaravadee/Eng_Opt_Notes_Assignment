# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Coursework for **Engineering Optimisation (M-EOP)** at SDU: lecture notes as PDFs
plus small, self-contained MATLAB demo scripts that accompany them.

This is not an application. There is no build, no package manager, no test
framework, and no shared library — every script runs standalone and exists to
make one pedagogical point from its lecture. Treat the scripts as *teaching
artefacts whose numbers and vocabulary are pinned to the lecture notes*, not as
code to be refactored toward reuse.

## Layout

```
L01_formulation/           L01.pdf, L01_slides.pdf + demos/Matlab/*.m
L02_linear_programming/    L02_slides.pdf, Handwritten_note_L02.pdf + demo/*.m
```

One directory per lecture, `L<nn>_<topic>/`, each holding the lecture PDFs and a
demo folder. **The demo folder name is inconsistent between lectures**
(`demos/Matlab/` in L01, `demo/` in L02) — this comes from the upstream course
zips, so don't normalise it unprompted. The `.zip` files at the root are those
upstream distributions, kept as-received.

Course material arrives as new `L<nn>_*.zip` files; extracting one into its own
folder is how a lecture gets added.

## Running a demo

```bash
cd L02_linear_programming/demo && matlab -batch l02_demo_ex32
```

**Always `cd` into the demo folder first.** Scripts call
`exportgraphics(gcf, 'name.pdf')` with a bare filename and probe for
`../../figures/` — both resolve against the *current working directory*, not the
script location. Running from the repo root is why `l02_ex32_graphical.pdf` and
`l02_parallel_edge.pdf` are sitting at the root instead of beside their scripts.

MATLAB's cold start is very slow here (tens of seconds at best; one measured
attempt exceeded 200s without finishing). Batch a whole run rather than
iterating one script at a time, and never assume a hung-looking MATLAB has
failed.

`-batch` is headless and non-interactive: figure creation and `exportgraphics`
still work, but the script must not wait for input.

## Environment check

```bash
cd L01_formulation/demos/Matlab && matlab -batch check_toolboxes
```

Verifies toolboxes, resolves every solver the course uses, and finishes with a
real `fmincon` solve (a product can be listed but unlicensed — only the solve
catches that). Its header notes it was written without MATLAB available and
never run by its author; it is built only from `ver` and `exist`, so verify by
hand if it misbehaves.

Two toolboxes beyond base MATLAB, with the lecture that first needs each:

| Solver | Lecture | Toolbox |
|---|---|---|
| `linprog` | L2 | Optimization |
| `fminunc` | L3 | Optimization |
| `fmincon` | L4 | Optimization |
| `ga`, `particleswarm` | L7 | Global Optimization |
| `fgoalattain`, `paretosearch` | L8 | Global Optimization |
| `intlinprog` | L9 | Optimization |

## How a demo script is built

Every script follows the same shape, and new ones should:

1. **Header block** — the problem in mathematical notation with units, the
   reasoning behind each coefficient, any modelling assumption stated as an
   assumption, a `Requires:` line, and a provenance line
   (`RUN AND VERIFIED in MATLAB R2025b on 2026-08-19`). Update that line when
   you re-verify; don't inherit a stale one.
2. `clear; clc; close all`
3. **Problem data as explicit `c` / `A` / `b` / `lb` / `ub`**, laid out so the
   matrix rows visually match the constraints written in the header, with a
   trailing comment naming each row's resource.
4. **Solve** with `optimoptions(..., 'Display', 'none')`, requesting the
   `lambda` output whenever shadow prices are part of the lesson.
5. **Report block** — `fprintf` that *interprets* the answer (which resource
   binds, what a zero price means, why a fractional answer is honest rather than
   an error), not just a dump of numbers.
6. **Figure**, then `exportgraphics` (see the dual-write rule below).
7. **Local functions** at the bottom where needed.

## Conventions that are deliberate — preserve them

**Minimisation is stated directly.** Problems are written as
`minimise f = -profit` from the outset and the profit reported as `-fval`. There
is no negate-then-negate-back step; don't reintroduce one.

**Vocabulary is staged by lecture.** L01 says `dfstar/db` and deliberately does
*not* call it a Lagrange multiplier, because the sign convention needs the
Lagrangian from L04. Introducing later terminology into an earlier demo breaks
the notes. Script comments flag this explicitly, including references to
numbered "Course Plan" open items — that plan is not in this repository, so
raise such a reference rather than acting on it.

**Hand-derived values are hardcoded and cross-checked — this is the test
suite.** `l02_demo_ex32.m` carries `hand_x`, `hand_p`, `hand_s`, `hand_y` from
the note's final simplex tableau and asserts agreement at `tol = 1e-6`;
`l01_demo_1D_fmincon.m` computes an exact answer by enumeration and prints
`NO <-- only a LOCAL minimum` when the solver disagrees. **A script that prints
disagreement has failed**, even though nothing throws. Keep these checks and
their tolerances when editing, and read the printed output rather than just the
exit status.

**Figures dual-write.** Each figure goes beside the script *and* into
`../../figures/` when that folder exists, because the lecture notes embed it.
The folder does not currently exist; the `isfolder` guard means the second write
silently no-ops, which is intended.

**Comment density is a feature.** Comments explain *why*, including plotting
subtleties — keeping `yline`/`xline` handles out of legends, positioning an
annotation below the axis because `df*/db` is negative. Match this register;
terse code is the wrong style here.

**Some demos are built to fail instructively.** `l01_demo_1D_fmincon.m` is
designed so changing `x0` from 2.5 to 0.5 makes `fmincon` return the wrong
answer for three of four bounds *while reporting success every time*. Don't
"fix" a demo into always converging — that destroys the lesson.

## Git

The repository has no commits yet and everything is untracked. Generated PDFs
(`l02_*.pdf` at the root) are script output, not source.
