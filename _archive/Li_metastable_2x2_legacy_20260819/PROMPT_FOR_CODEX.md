# Codex Master Prompt — Li Σ5 [001](310) Metastable Grain Boundary

## 0. Task objective

Create a **new, isolated production/research workflow** for studying a metastable grain-boundary structure in BCC Li.

The system is fixed as:

* Material: **BCC Li**
* Grain boundary: **symmetric tilt Σ5 [001](310)**
* Misorientation: **36.86989765°**
* Macroscopic GB type: **must remain exactly the same as the validated old model**
* Potential:
  `SNAP_ZuoChenLi_2019_Li__MO_732106099012_000`

The scientific goal is NOT to search for a different GB type.

The goal is:

> Starting from the validated Li Σ5 [001](310) bicrystal, thermally reconstruct the GB at approximately 80% of the Li melting temperature, rapidly quench it, and determine whether the resulting disordered/reconstructed GB can remain stable for more than 1 ps at 10 K and 0 bar while the grain interiors retain BCC structure.

The advisor's criterion is:

> If the system can exist stably for more than 1 ps at 10 K and zero pressure, it can be regarded as worth continuing to study. Phase-transition classification is not required at this stage.

---

# 1. Use the old validated project as the starting point

The authoritative old project is:

`C:\Users\1\Documents\Codex\2026-08-01\bcc-li-5-001-sigma5-310`

Inspect the actual directory before doing anything.

Do NOT guess file contents from filenames.

Reuse the validated parts of the old workflow, including where applicable:

```text
generator/
lammps/
tests/
validated input structures
parameter/provenance files
README / final status / environment information
```

In particular inspect and reuse where appropriate:

```text
lammps/in.0_check
lammps/in.1_minimize
lammps/in.2_anneal
lammps/in.2_anneal_resume
lammps/in.3_quench
lammps/in.3_quench_resume
lammps/in.4_final
lammps/in.dump2data
lammps/in.bulk_energy
lammps/run_all.sh
lammps/resume.sh
```

Also inspect the validated generator and its tests.

The old package is a **reference/provenance source**. Do not overwrite it.

Create a completely new experiment directory.

Suggested name:

```text
Li_Sigma5_310_metastable_2x2
```

---

# 2. Critical scientific interpretation

The new experiment is based on the **same macroscopic Σ5 [001](310) GB** as the old model.

Do NOT:

* change Σ5 to another Σ value;
* change the [001] rotation axis;
* change the (310) macroscopic GB plane;
* search for another GB geometry;
* interpret atomic reconstruction as a change of macroscopic GB type.

Distinguish explicitly between:

### Macroscopic GB character

```text
Σ5 [001](310)
```

and

### Microscopic atomic configuration

The latter is allowed to reconstruct strongly during high-temperature MD.

The scientific question concerns the **microscopic GB structure**, not changing the macroscopic crystallographic identity.

---

# 3. Model size: fixed decision

The old validated model contains approximately:

```text
3600 atoms
```

For the new experiment:

```text
GB-plane direction 1 → 2×
GB-plane direction 2 → 2×
GB-normal direction → unchanged
```

Therefore the GB area becomes approximately:

```text
4×
```

and the atom count should be approximately:

```text
3600 × 4 ≈ 14400 atoms
```

The exact atom count must be obtained from the generated structure, not assumed.

## VERY IMPORTANT

Do not blindly multiply an old generator tuple by 2.

First determine from the authoritative old generator/output which physical directions lie in the GB plane.

The old production documentation indicates that the GB plane is the **yz plane** and x is the GB-normal direction.

However, because earlier generator versions may use different coordinate conventions, verify this against:

1. authoritative generator;
2. actual generated structure;
3. validated old output.

The physical requirement is:

> **double the two dimensions parallel to the GB, not merely two arbitrary generator parameters.**

Do not change the GB-normal thickness unless there is a scientifically justified reason.

---

# 4. Stage 0 — Build and validate the 2×2 model

## Objective

Generate a 2×2 in-plane version of the validated Σ5 [001](310) Li bicrystal.

Prefer modifying the generator's supercell parameters and regenerating the model rather than blindly scaling an already generated structure.

## Validate

Check:

* atom count;
* box dimensions;
* PBC;
* GB plane;
* GB-normal direction;
* grain orientation;
* Σ5 relationship;
* misorientation;
* grain IDs;
* minimum interatomic distance;
* absence of accidental overlaps;
* bulk BCC structure;
* correct Li element/mass;
* correct potential compatibility.

Expected scale:

```text
~14400 atoms
```

but report the actual number.

## Outputs

Create:

```text
stage0/
    structure.data
    structure.xyz
    grain_map
    validation_report
    parameters.json
```

Do not mix Stage 0 outputs with old completed outputs.

---

# 5. Stage 1 — Initial relaxation/equilibration

Start from the new 2×2 structure.

Reuse the validated old minimization/equilibration logic as much as possible.

The purpose of Stage 1 is:

> obtain a clean and mechanically reasonable initial Σ5(310) GB.

It is NOT the metastable-state generation stage.

Perform:

1. energy minimization;
2. validated equilibration procedure from the old workflow;
3. structural sanity checks.

Save:

```text
stage1/
    minimized.data
    trajectory
    log
    restart
    thermodynamic_history
```

Record:

* potential energy;
* temperature;
* pressure;
* box dimensions;
* structural indicators.

---

# 6. Stage 2 — High-temperature NPT reconstruction

This is the key new experiment.

The advisor's proposed procedure is:

> After the GB has been established, directly use an NPT ensemble to give the system a high temperature, for example 80% of the Li melting temperature, hold for 2 ps, and look for a disordered/reconstructed GB while the grain interiors remain BCC.

## Target temperature

Li melting temperature is approximately:

```text
454 K
```

80% gives approximately:

```text
0.8 × 454 ≈ 363 K
```

Therefore use:

```text
T_high = 363 K
```

unless the validated implementation gives a scientifically justified reason to use a nearby value.

## Ensemble

```text
NPT
T = 363 K
P = 0 bar
duration = 2 ps
```

Use the validated timestep from the old workflow unless there is a clear reason to change it.

Target timestep:

```text
0.001 ps
```

if consistent with the existing workflow.

## Purpose

The purpose is NOT to melt the system.

The desired condition is:

```text
GB:
    strongly reconstructed / disordered

Bulk grains:
    still recognizably BCC
```

The target is therefore:

```text
GB disorder
+
bulk BCC retention
```

---

# 7. Stage 2 must be evaluated before quenching

Do not automatically proceed to Stage 3 just because the 2 ps run finished.

Inspect the complete trajectory.

At minimum inspect:

```text
t = 0
t ≈ 0.5 ps
t ≈ 1 ps
t ≈ 1.5 ps
t = 2 ps
```

The important question is:

> Did the GB undergo a genuine thermally activated reconstruction while the bulk crystal remained ordered?

## Stage 2 classification

### Case A — Insufficient reconstruction

GB remains essentially unchanged.

Then:

```text
Stage 2 = insufficient
```

Do not claim success.

Possible future parameter studies may vary temperature or hold time, but do not silently change the protocol.

### Case B — Desired reconstruction

GB becomes clearly disordered/reconstructed and bulk remains BCC.

Then:

```text
Stage 2 = PASS
```

and proceed to Stage 3.

### Case C — Global melting/disordering

If the grain interiors lose BCC structure significantly:

```text
Stage 2 = FAIL
```

Do not treat this as the desired metastable GB.

---

# 8. Stage 3 — Rapid quench

Starting from the accepted Stage-2 high-temperature configuration:

```text
T ≈ 363 K
↓
10 K
```

with:

```text
cooling rate = 10 K/ps
```

Temperature interval:

```text
363 - 10 = 353 K
```

Quench duration:

```text
353 / 10 = 35.3 ps
```

Therefore:

```text
quench ≈ 35.3 ps
```

Use a clearly defined temperature-control protocol.

Do not silently alter the cooling rate.

Save the complete quench trajectory.

Also save restart/checkpoint files.

---

# 9. Stage 3 analysis

Do NOT only inspect the final frame.

Track the GB structure throughout:

```text
high-T reconstructed GB
        ↓
cooling
        ↓
partial relaxation
        ↓
frozen/reconstructed low-T GB
```

It is acceptable for the GB to relax somewhat during cooling.

The question is whether the essential reconstructed/disordered character survives.

Track:

* temperature;
* pressure;
* potential energy;
* local structure;
* GB width;
* bulk BCC fraction;
* atomic motion.

The fact that atomic motion becomes very small during quenching is NOT by itself evidence of metastability.

---

# 10. Stage 4 — 10 K / 0 bar stability test

After the quench:

```text
T = 10 K
P = 0 bar
```

Hold the system.

Advisor's minimum criterion:

```text
> 1 ps
```

Recommended actual production hold:

```text
5 ps
```

if computationally practical.

10 ps may be used if inexpensive.

## Minimum success criterion

The candidate passes if:

1. GB disorder/reconstruction remains;
2. grain interiors remain BCC;
3. there is no obvious structural collapse;
4. the system survives at 10 K and 0 bar for more than 1 ps.

Formally:

$$
T=10K
$$

$$
P=0
$$

$$
t>1ps
$$

with persistent GB disorder.

---

# 11. Important distinction: frozen does not automatically mean metastable

This is a major lesson from the previous attempts.

At 10 K, atomic motion can become extremely small.

Therefore:

```text
low atomic displacement
≠
proof of metastability
```

The structure must be shown to remain a distinct GB configuration rather than simply being kinetically frozen by the low temperature and finite simulation time.

Therefore compare:

```text
Stage 1 initial GB
Stage 2 high-T GB
Stage 3 quenched GB
Stage 4 final 10-K GB
```

and determine whether the high-T reconstruction survives.

---

# 12. Stage 5 — Structural and physical audit

Only candidates passing Stage 4 should proceed here.

Do not begin sophisticated SOAP analysis before the basic MD survival test is passed.

## 12.1 Bulk BCC analysis

Use a reproducible local structure method such as:

* CNA / adaptive CNA;
* PTM;
* another validated BCC classifier.

Quantify the BCC fraction in the grain interiors.

The goal is to distinguish:

```text
localized GB disorder
```

from:

```text
global crystal disorder
```

---

# 13. Stage 5 — GB disorder profile

Construct a spatial profile across the GB.

For example:

$$
Q_\mathrm{BCC}(x)
$$

or another local-order parameter.

Because the GB normal is x in the authoritative production convention, the profile should normally be evaluated along x, after verifying the actual coordinate convention.

Also calculate a local-energy profile:

$$
E_\mathrm{local}(x)
$$

and, where appropriate, a disorder profile:

$$
D(x)
$$

Use these to estimate:

$$
w_\mathrm{GB}
$$

the effective GB disorder width.

Compare:

```text
initial
high-T
quenched
10-K final
```

---

# 14. Stage 5 — Energy analysis

Track:

$$
E_\mathrm{pot}(t)
$$

but do not use total energy as the sole success criterion.

Energy should be treated as supporting evidence.

A stable-looking energy curve can coexist with an undesirable structural state.

---

# 15. Stage 5 — SOAP analysis

Only after a candidate passes the MD criteria.

Use SOAP as a **structural descriptor**, not as the primary proof of metastability.

The previous Ni SOAP work demonstrated an important methodological point:

> A difference in SOAP representation does not by itself establish physical metastability or an independent energetic mechanism.

Therefore the logical order must be:

```text
MD stability
      ↓
structural verification
      ↓
energy / GB profile
      ↓
SOAP characterization
```

NOT:

```text
SOAP difference
      ↓
claim metastability
```

For accepted Li candidates, compare SOAP environments among:

* initial GB;
* high-T reconstructed GB;
* quenched GB;
* final 10-K GB;
* bulk reference.

---

# 16. Stage 6 — Size dependence

Do NOT immediately run a 3×3 model.

First establish whether the 2×2 model can produce a successful candidate.

If successful:

```text
2×2
↓
3×3
```

Repeat the same protocol.

Compare:

* GB disorder;
* GB width;
* bulk BCC fraction;
* local energy;
* final structure;
* 10-K survival;
* SOAP characteristics if needed.

The purpose is to determine whether the result is strongly dependent on finite in-plane dimensions.

The first 2×2 success should therefore be treated as:

```text
candidate metastable state
```

not as proof of a size-independent physical phenomenon.

---

# 17. Lessons from previous failed work

These rules are mandatory.

## 17.1 Do not modify the old authoritative project

Never overwrite the old validated package.

Create a new experiment directory.

---

## 17.2 Do not reuse old completed outputs as new outputs

Old files such as:

```text
stage1_minimized.data
stage2_annealed.data
stage3_quenched.data
stage4_final.data
```

may be used as reference/provenance.

They must NOT be silently presented as results from this new experiment.

---

## 17.3 Do not mix lattice constants or construction conventions

The old workflow contains validated lattice/provenance information.

Do not blindly replace it with the earlier nominal:

```text
a = 3.51 Å
```

if the validated production model uses a relaxed value.

Inspect the authoritative old parameters first.

---

## 17.4 Do not confuse finite-size freezing with metastability

This is one of the central lessons.

Small periodic GB cells can artificially constrain atomic rearrangement and can make a rapidly quenched structure appear more persistent than it would in a larger in-plane model.

That is why the new experiment uses:

```text
2× GB-plane dimensions
```

and later a 3×3 size check.

---

## 17.5 Do not use one final snapshot as proof

Always retain trajectories for:

```text
Stage 2
Stage 3
Stage 4
```

The formation and freezing history matters.

---

## 17.6 Do not use total energy alone

Energy is necessary but insufficient.

Always combine:

```text
energy
+
local structure
+
spatial localization of disorder
+
time evolution
```

---

## 17.7 Do not use SOAP as a substitute for dynamics

SOAP describes local atomic environments.

It does not directly establish dynamical stability.

---

## 17.8 Do not make phase-transition claims yet

The advisor explicitly says the phase-transition issue can temporarily be ignored.

Therefore the current task is simply:

> Does a structurally disordered Σ5(310) GB survive the specified thermal protocol at 10 K and 0 bar?

Do not spend computational effort trying to classify the result as a particular thermodynamic phase unless later requested.

---

## 17.9 Remember that the periodic bicrystal contains two GBs

When calculating GB energy or interpreting the structure, account for the fact that the periodic bicrystal contains two interfaces.

Do not accidentally interpret a two-GB cell as a one-GB system.

---

## 17.10 One successful trajectory is not an ensemble statement

A single successful thermal trajectory demonstrates:

```text
a candidate pathway exists
```

It does NOT demonstrate:

```text
all Σ5(310) structures are metastable
```

or:

```text
this is the unique metastable state
```

If later required, multiple independent seeds/trajectories should be added.

---

# 18. File organization

Create:

```text
Li_Sigma5_310_metastable_2x2/
│
├── README.md
├── PROTOCOL.md
├── PARAMETERS.json
├── RUN_STATUS.md
├── PROVENANCE.md
│
├── legacy/
│   └── references only
│
├── generator/
│
├── stage0_build_validate/
│
├── stage1_relax/
│
├── stage2_highT_NPT/
│
├── stage3_quench/
│
├── stage4_10K_hold/
│
├── stage5_analysis/
│
├── stage6_size_effect/
│
└── scripts/
```

Keep outputs stage-specific.

---

# 19. Required provenance

Record:

```text
old project path
old generator version
new generator changes
old lattice parameters
new supercell parameters
actual atom count
box dimensions
GB orientation
potential name
LAMMPS version
OpenKIM version
timestep
thermostat parameters
barostat parameters
random seeds
Stage 2 temperature
Stage 2 duration
Stage 3 cooling rate
Stage 3 duration
Stage 4 duration
```

Every scientifically relevant parameter change must be recorded.

---

# 20. Required validation before expensive runs

Before starting a long simulation:

1. inspect old project;
2. generate 2×2 model;
3. validate geometry;
4. validate atom count;
5. validate PBC;
6. validate minimum distance;
7. validate grain IDs;
8. validate BCC bulk;
9. validate LAMMPS input syntax;
10. validate SNAP potential loading;
11. run a short smoke test;
12. only then run the production trajectory.

Do not launch all expensive stages blindly.

---

# 21. Required final report

At the end, report:

## Model

```text
GB:
Σ5 [001](310)

in-plane scaling:
2× × 2×

actual atom count:
...

box:
...

GB-normal:
...

potential:
...
```

## Stage 0

* generation success/failure;
* geometry validation;
* BCC validation.

## Stage 1

* minimization/equilibration result;
* energy;
* structure.

## Stage 2

* T;
* P;
* duration;
* whether GB reconstructed;
* whether bulk remained BCC.

## Stage 3

* initial temperature;
* final temperature;
* cooling rate;
* actual duration;
* trajectory result.

## Stage 4

* T = 10 K;
* P = 0;
* hold duration;
* GB disorder retained?;
* bulk BCC retained?;
* success/failure.

## Stage 5

* structural metrics;
* GB width;
* energy;
* SOAP characterization if applicable.

## Stage 6

* only if performed;
* size-dependence result.

---

# 22. Final scientific decision rule

The experiment is considered a **successful candidate** only if:

$$
\boxed{
\begin{aligned}
&\Sigma5[001](310)\ \text{macroscopic GB preserved}\\
&\text{GB becomes disordered/reconstructed at high T}\\
&\text{bulk grains remain BCC}\\
&\text{quench at }10K/ps\\
&10K,\ 0bar\\
&\text{GB disorder persists for } >1ps
\end{aligned}
}
$$

If all conditions are met:

```text
SUCCESS — candidate metastable disordered GB
```

If not:

```text
FAIL / INSUFFICIENT — do not overinterpret
```

Do not call a candidate successful merely because:

* the temperature reached 10 K;
* energy became flat;
* atoms barely moved;
* the final snapshot looks visually unchanged;
* SOAP distance is large.

The structural and dynamical criteria must support the conclusion.

---

# 23. Execution philosophy

Work conservatively.

First understand the validated old project.

Then make the smallest necessary modification to create the new 2×2 experiment.

Do not rewrite the entire old workflow unnecessarily.

The new scientific change is primarily:

```text
old validated Σ5(310) GB
        ↓
2×2 in-plane enlargement
        ↓
363 K NPT / 2 ps
        ↓
10 K/ps quench
        ↓
10 K / 0 bar hold
        ↓
structural stability test
```

Everything else should remain as close as reasonably possible to the validated old workflow.

If an ambiguity is discovered between the old files and this protocol, **do not silently choose an interpretation**. Inspect the authoritative source, document the discrepancy, and choose the physically correct implementation.

Do not claim a result that has not been demonstrated by the actual trajectory.
