# ELMxx

ELMxx is a second, alternate LND component (`COMP_LND=elmxx`), intended to become a
Kokkos/C++ port of Fortran ELM. It was added by mirroring how RDycore was added as
an alternate ROF, minus PETSc.

## Where things are — read this first

Three kinds of document, three homes, one rule each. **Do not create a fourth
kind.**

| Kind | Lives at | Rule |
|---|---|---|
| **Always-on context** | `components/elmxx/CLAUDE.md` (this file) | Facts that stay true across sessions: layout, invariants, build gotchas. **Never** holds status or plans. |
| **The plan** | `plans/active/` | **Exactly one file.** If there are two, one is wrong -- fix that before doing anything else. Long-lived; edited deliberately. |
| **Where we are now** | `plans/STATUS.md` | **One file, rewritten in place** each session. Never dated, never appended to. **Keep it under ~200 lines.** It carries a physics ledger and an outstanding-bug list, not narrative. |
| **Retired plans** | `plans/archive/` | Dated on retirement. Read-only thereafter. |

**Start any session by reading `plans/STATUS.md`, then the single file in
`plans/active/`.** STATUS says where the cursor is; the plan says where we are
going; this file says what is always true.

Why the split: a plan and a status snapshot have different lifecycles. The plan
changes rarely and on purpose; status changes every session. Keeping them as
peers in one directory -- which is how this started -- makes it impossible to tell
which is authoritative. If a finding is durable (a build gotcha, an invariant) it
belongs **here**; if it changes what we will do, it belongs in the **plan**; if it
is just where we stopped, it belongs in **STATUS.md**.

Do not add dated handoff files. That is what STATUS.md replaced.

**And do not let STATUS.md grow into one.** It reached 2063 lines once, 56% of
it a single section of accreted session narrative, and had to be retired to
`plans/archive/2026_08_16_status_stages_1_4.md`. The rule that prevents a
repeat: **a session's story goes in the commit message** -- searchable, attached
to the change, and never needing to be pruned. STATUS.md answers only four
questions: what is ported, what is validated, what is broken, what is next.

**Current state is NOT recorded here** -- it changes every session and this file
does not. See `plans/STATUS.md`, which carries the physics ledger (what is
ported, what is validated, with what evidence) and the outstanding-bug list.

## Running a case

Two helper scripts (untracked, in the `cime` submodule, alongside the pre-existing
`brazil.sh`):

```bash
cd cime/scripts
./brazil_elmxx.sh    # I1850ELMXX @ 1x1_brazil,  NTASKS=1
./f19_f19_elmxx.sh   # I1850ELMXX @ f19_f19,     NTASKS=4
```

Then `./case.build`, and run. Three gotchas, all of which will bite:

1. **`--driver mct` is required.** This repo defaults `COMP_INTERFACE=moab` (from
   `driver-moab/cime_config/config_component.xml`), for ELM cases too. ELMxx only
   ships an MCT layer, so `cime_config/buildnml` asserts `COMP_INTERFACE == mct`
   and fails at `case.setup` if not.
2. **`./case.submit` fails on this laptop** with an hwloc error -- the machine's
   `mpirun` args include `--bind-to hwthread`, which macOS does not support. Run
   the exe directly instead:
   ```bash
   cd run && mkdir -p timing/checkpoints
   mpirun -np 4 ../bld/e3sm.exe > e3sm.manual.log 2>&1
   ```
   `run/timing/checkpoints` must exist or the run aborts.
3. **`1x1_brazil` must run on exactly 1 rank.** It is a single grid cell; DATM
   aborts (`shr_tInterp_getCosz: lon lat cosz sizes disagree`) on ranks that get
   zero points, and `elmxx_init` aborts if `npes > numg`.

Check `run/lnd.log.*` for the ELMxx banner and per-rank cell counts, and
`run/cpl.log.*` for `lnd model present = T` and `--- checking atm/land domains ---`.

## Source layout

| File | Role |
|---|---|
| `src/cpl/lnd_comp_mct.F90` | MCT layer: `lnd_{init,run,final}_mct`, gsMap, domain |
| `src/main/elmxxMod.F90`    | Decomposition + namelist + `elmxx_{init,run,final}` |
| `src/main/elmxxIO.F90`     | PIO read of the land domain file |
| `src/main/elmxxSpmdMod.F90`| `mpicom_lnd`, `iam`, `npes`, `masterproc`, `LNDID` |

`elmxxMod.F90` deliberately mirrors `components/rdycore/src/main/rdycoreMod.F90`
(same public names: `num_cells_owned`, `num_cells_global`, `natural_id_cells_owned`)
but has no PETSc.

There is no `elmxx_cpl_indices.F90` yet -- add it (from
`components/elm/src/cpl/elm_cpl_indices.F90`) in the same change that first reads
field values out of `x2l_l`.

## Naming — standing rules for ELMxx C++ (owner-approved 2026-10-04)

Goal: an ELM developer finds the ELMxx counterpart of any ELM subroutine by
grepping its ELM name. Details and the full rename inventory: development
plan, Part 5.

- **R1 Files** are named after the ELM module they port, minus `Mod`
  (`SoilMoistStressMod` → `SoilMoistStressImpl.h` / `…Data.h`).
- **Directories** mirror ELM's `src/`: `src/{biogeophys,biogeochem,data_types,main,utils}`
  in the submodule, each module's header and source together; `include/` holds
  only the public API (`ELMxx.h`, `ELMxxMacros.h`, `elmxx_mod.F90`, `finclude/`).
  The include root is `src/`: `#include "biogeophys/SoilTemperatureImpl.h"`.
  Containers live in `data_types/`, the C API and history/restart in `main/`.
- **R2 Kernel functions** carry the ELM subroutine name, plus at most one element
  suffix: `Column`, `Patch`, `Landunit`.
- **R3** Rename only where ELMxx uses *different words* from ELM. Pure style
  variants stay (`PhaseChangeBeta` for `Phasechange_beta`, `SnicarAdRt` for
  `SNICAR_AD_RT`) and are covered by R6.
- **R4 Data structs** for one kernel: ELM subroutine name + `Data`.
- **R5 Field names are not changed to ELM's spelling, and ELMxx does not add
  underscores.** ELMxx's camelCase names are its convention. **Snow and soil are
  separate arrays, the snow index reversed relative to ELM's** (`SnowIndexing.h`),
  and the urban per-surface members stay — one kernel serves all five surfaces.
- **R6** Every ported function's header names its ELM origin by subroutine, not
  line number: `// ELM: SoilTemperatureMod :: Phasechange_beta`. Every field
  declaration names its ELM variable: `// ELM: col_ws%h2osoi_liq`,
  `// ELM: SoilFluxes local egsmax`. Check an `inst%member` against ELM's type
  definitions before writing it: 68 earlier references named owners that do not exist.
- **R7** ELMxx-only constructs say so: `// ELMxx-only: block inside
  HydrologyNoDrainage, no ELM subroutine`.
- Validation test names and namelist kernel tokens keep their current names
  (tests match instrumented ELM's diagnostic labels; tokens are user-facing).

## Coupler constraints -- do not "fix" these

Read out of `driver-mct/main/seq_domain_mct.F90`; each one aborts the run if broken.

- **gsMap `gsize` is `nlon_g*nlat_g`, not `num_cells_global`.** The segments cover
  only active land cells (`mask == 1`, ELM's `numg`), but the global size must match
  the atm grid or `seq_domain_mct.F90:211` aborts on `gatmsize /= glndsize`. This is
  exactly what ELM does in `lnd_setgsmap_mct`.
- `natural_id_cells_owned` holds global **grid** IDs (indices into the `ni*nj` arrays),
  not land-cell ordinals, so it indexes the domain arrays directly.
- **`area` passes through unconverted.** The domain file stores radians^2, which is
  what the coupler wants. (ELM multiplies by `re**2` on read and divides by `re*re`
  on export -- a round trip.)
- **`aream` is left at `-9999`.** The atm-lnd mapper fills it in.
- `mask` and `frac` come straight from the domain file, which is what satisfies
  `seq_domain_check_fracmask` (it aborts on `frac /= 0` where `mask == 0`).
- **`lnd_domain` is a MOAB-only `seq_infodata_PutData` argument.** It exists in
  `driver-moab/shr/seq_infodata_mod.F90` but not in driver-mct's; passing it will
  not compile under mct.

## Namelist

The file is **`lnd_in`** (namelist files are named after the component *class*, so
ELMxx uses the same name as ELM), group `elmxx_inparm`, holding `do_elmxx` and
`fatmlndfrc`. It is produced by the perl `bld/build-namelist` (adapted from
`components/rdycore/bld/build-namelist`), driven by `cime_config/buildnml`.
`fatmlndfrc` is set from `LND_DOMAIN_PATH`/`LND_DOMAIN_FILE` -- change those with
`./xmlchange`, not by overriding the namelist, or the domain will stop matching the
one the coupler checks against.

## Registration lives outside this directory

Adding or renaming compsets also touches, in the repo root:

- `cime_config/config_files.xml` -- `COMP_ROOT_DIR_LND`, `COMPSETS_SPEC_FILE`,
  `PES_SPEC_FILE` entries for `elmxx`
- `cime_config/config_archive.xml` -- the `elmxx`/`lnd` spec
- `cime_config/allactive/config_compsets.xml` -- the `LND = [...]` help line

Note the regex `_ELM` also matches `_ELMXX`, so ELMXX inherits ELM's defaults in
`driver-mct/cime_config/config_component_e3sm.xml` and ELM's grid aliases in
`cime_config/config_grids.xml`. That is currently desirable, but it means
`ELM_USRDAT` grids appear eligible while `ELM_USRDAT_NAME` is undefined for elmxx
cases -- use a standard grid.

## Not supported yet

Multi-instance (`NINST_LND > 1` is asserted against in `buildnml`), the moab
driver, `RUN_TYPE=hybrid` (it cold-starts ELMxx: `finidat` reads ELM files,
not ELMxx restarts), and restarting urban-column / urban-patch state
(never allocated on the coupled path; `elmxxRestMod` aborts naming the field
if one ever is). Lake state restarts (`LAKECOL`/`LAKEPATCH`, since Stage 6.5).

## The lake surface

Lake is its own surface, a peer of natural and urban: `lakeCol`/`lakePatch`
(`src/data_types/Lake{Column,Patch}Data.h`), driven from
`src/main/elmxxLakeMod.F90`. Durable facts:

- **Seeded and read by NAME** through the restart registry
  (`ELMxxRestartFieldFind` + `ELMxxRestartFieldGet/Set`, e.g.
  `"lakecol:t_lake"`), not per-field setters. The registry crosses a flat
  ROW-major buffer; `elmxxLakeMod`'s helpers transpose Fortran arrays.
- **ELM's combined snow+ground index space**, not naturalCol's reversed
  snow index: slot j <-> ELM layer j-4; `zi` slot i <-> ELM `zi(i-5)`.
  LakeTemperature reads `zi` ONE SLOT EARLIER than LakeHydrology, so it gets
  `zi_lt`, refreshed from `zi` each step. Never hand it `zi` directly.
- **One patch per lake column, packed alike** (LakeHydrology indexes
  p = c); `ELMxxAllocateLakeSurface` refuses anything else.
- **Forcing comes from the cell's natural column** (`col_natcol`), which
  every cell has: ELM allocates the natural landunit whatever its weight.
- **History is the natural+lake landunit mean**, per field by `LakeRule`
  (`HistData.h`), each rule read off ELM's own h0 over a 100% lake cell
  (`brazil_lake100`), not guessed. The h0 guard is natural + lake weight = 1.
- **ELM's lake column keeps `frac_sno` = 0 under snow** (lake code never
  sets it), so its lake ground albedo is the snow-free one, and SNICAR enters
  only through the per-layer factors, weighted `(1-albsod)*flx/(1-albsnd)`
  (ELM's "subgridflag == 0 .or. is_lake" branch) -- not the natural path's
  `(1-albsnd)` weighting. `SnowAge_grain` on lake reads `frac_sno_eff` (1).
- **Lake snow carries aerosol, and it dominates its albedo** (0.2 by March
  on the glc lake twin): deposition, meltwater scavenging, combine/divide,
  dead-slot zeroing all run in `LakeHydrologyImpl.h`. Optional on
  `LakeHydrologyData` (the replay structs carry none).
- **Lake negative controls need >= 1e-6 K.** Convective mixing averages a
  1e-12 K perturbation of `t_lake` below a 300 K ULP; it vanishes.

## Restarts

ELMxx writes ALL persistent state (`src/main/elmxxRestMod.F90`): every field of
the C++ restart registry (`ELMxxRestartField*`, generated from the state-container
headers by `external_models/elmxx/tools/gen_restart_registry.py`), the active
history accumulators, and the Fortran-side `nstep` and history interval. Written
when the coupler's restart alarm rings, dated with the step's END time
(`<case>.elmxx.r.YYYY-MM-DD-SSSSS.nc`), plus `rpointer.lnd`. Read on
`CONTINUE_RUN=TRUE` (from `rpointer.lnd`) or a branch (`elmxx_nrevsn`).

- **The read happens at the top of the first `elmxx_run`**, after forcing the
  first-step lazy inits and BEFORE the forcing push and root stress. Moving it
  later breaks bit-for-bit: forcing views are restart fields, and `btran` would
  come from the cold start.
- **A continue run takes no duplicate nstep-0 pass** and its clock starts at the
  coupler's current time, as ELM's does.
- **Layout-independent, as ELM's are:** a restart written on N ranks reads on
  M. Each registry field carries an entity KIND (`natcol`, `natpatch`, `cell`,
  `urblun`, ...; emitted by the generator from the view's first allocation
  extent -- an unclassified extent fails generation), the file's leading dim is
  that kind over the whole domain, and an entity's position is (global cell
  ordinal, occurrence within the cell), per-cell counts from an allreduce.
  Packed urban order is (density class, cell), so this is a real reorder.
  `TOPO` fields (`col_gridcell`, `patch_column`, `filter_*`, ...) hold
  rank-local indices: rebuilt at init, never written or read. The
  `KindMatchesLeadingExtent` CTest and a per-field extent check at write/read
  guard the classification.
- **Reading is strict:** every allocated registry field must be on the file with
  the same shape, or the run aborts. A new view must be added to the registry
  (re-run the generator; the `RestartRegistrySync` CTest fails otherwise) and
  then old restart files no longer read.
- **Closure test (passes bit-for-bit, 2026-10-04):** run A 45 days from cold
  start with `REST_OPTION=ndays REST_N=15`; set A's day-45 restart and Jan h0
  aside; point `rpointer.lnd`, `rpointer.drv` AND `rpointer.atm` at the day-15
  date; `CONTINUE_RUN=TRUE STOP_N=30`; compare with
  `ncdump -p 17,17 | tail -n +2` (drop `date_written`/`time_written`). A mid-month
  restart must reproduce that month's h0 exactly. Negative control: perturb a
  PROGNOSTIC field (soil T) in the restart by 1e-12 -- not `t_grnd`, which
  CanopyTemperature recomputes at the start of every step.
- **Layout test (passes, 2026-10-04):** a 4-cell grid at four coordinates
  (so cells differ; `tools/make_4cell_rsttest.py` builds it from the 2x1
  files: urban region 0,3,3,3, zero urban weight -- h0 needs natural weight
  1), the 2x1 case with the domain/fsurdat swapped, `REST_N=5`, 10 d,
  on 1 and on 2 ranks: h0 and both restarts identical across layouts; a 2-rank
  continue from the 1-rank day-5 restart, and vice versa, reproduce the day-10
  restart and the Jan h0. The driver rejects another case's `cpl.r` on a
  continue, so each run keeps its own `cpl.r`/`datm.rs1` and only
  `rpointer.lnd` points at the other layout's file. A grid of identical cells
  cannot catch a permutation -- every reorder compares equal.
- **A branch starts history fresh**, as ELM's does (`hist_restart_ncd` zeroes
  `ntimes` on `nsrBranch`): the `hist:*` accumulators and the interval
  start/count are not read, so the first h0 holds only the branch's own steps
  and `elmxx_hist_fincl` may change. `buildnml` sets `elmxx_nrevsn` to
  `$RUNDIR/$RUN_REFCASE.elmxx.r.$RUN_REFDATE-$RUN_REFTOD.nc` unless
  `user_nl_elmxx` sets it. **Branch test (passes, 2026-10-04):** A = 45 d
  startup, `REST_N=15`; B = clone with `RUN_TYPE=branch RUN_REFCASE=A
  RUN_REFDATE=0001-01-16 GET_REFCASE=FALSE STOP_N=30`, A's 01-16 `cpl.r`,
  `datm.rs1` and `elmxx.r` copied into B's run dir and `rpointer.drv`/`.atm`
  written by hand. B's day-45 restart and Feb h0 equal A's; B's Jan h0 has
  768 samples (Jan 16–31). **`create_clone --keepexe` keeps the source's
  `RUNDIR`** -- `./xmlchange RUNDIR=$PWD/run` before `preview_namelists`, or
  the clone overwrites the source's namelists.
- `ncdump -v` cannot select the registry's `natcol:...` names (the colon); dump
  the whole file.


## Validating a Fortran-side port against ELM

ELMxx's Fortran-side ports (`SurfaceAlbedo`, the ground heat flux, `btran`) are
validated by **replaying ELM's own inputs through the port**, not by comparing
two runs. The two models' states differ and are allowed to; a replay removes
that entirely.

**The shape that makes it possible** — copy it for any new port:

| File | Holds | Depends on |
|---|---|---|
| `elmxx<X>KernelMod.F90` | the physics | **its arguments only** — no ELMxx object, no subgrid maps, no setters, no `shr_sys_abort` |
| `elmxx<X>Mod.F90` | gather, call, push | everything else |

The kernel then links into a standalone `gfortran` program, so the test drives
**the same source the coupled run calls** rather than a copy of it. Two worked
examples, and they differ in where the reference comes from:

| Port | Reference | Tools |
|---|---|---|
| `SurfaceAlbedo` | an ELM **restart** — the fields are model state | `surface_albedo_replay.F90`, `validate_surface_albedo.py` |
| ground heat flux | ELM's **diagnostic binary** — `hs_soil` and friends are intermediates, computed and consumed inside one call to SoilTemperature, so no restart ever holds them | `ground_heat_flux_replay.F90`, `validate_ground_heat_flux.py` |

Prefer the diagnostic binary where it exists: it gives one sample per timestep
rather than one per restart, and it needs no twin re-run. See `plans/STATUS.md`
§J.2 for the commands.

**Constants must come from the same file the model reads.** The checker that
hardcoded the textbook `5.67e-8` for Stefan-Boltzmann instead of E3SM's
`5.670374419e-8` reported a 0.02 W/m² disagreement on `hs_soil` — small enough
to shrug at, large enough to hide a real one. Parse `SHR_CONST_*` out of
`share/util/shr_const_mod.F90`.

**A validation is only as good as what the twin exercises.** Measure the
coverage and say what was NOT tested: on the brazil twins `snl` is 0 at every
step, so any snow branch is exercised only at its no-snow values, and
`hs_top_snow` is identically `hs_soil`.

**Never compare a temperature relatively.** A soil temperature is ~274 K, so a
relative tolerance of 1e-7 permits 2.7e-5 K of error — while the deep layers
of a 5-day run sit ~1e-6 K above their initial condition. Those layers can be
*entirely* wrong and every relative check still passes; that is exactly how the
`nlevbed` solve truncation survived. Assert on **absolute** error in kelvin,
and report it **per layer**: the shape of a disagreement identifies its cause,
where the maximum alone does not. A layer whose disagreement equals its
movement is a layer the port never updated.

### Kokkos kernels are validated by the CTest suite, not by a replay driver

Where the physics is already C++ there is nothing to extract, and a data-driven
test usually exists in `tests/validation/`. Build and run it — do not build a
third harness:

```bash
E=~/projects/e3sm/elmpy/e3sm/cime/scripts E2=components/elmxx/external_models/elmxx
cmake -B <build> -S $E2 -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_ENABLE_OPENMP=OFF -DKokkos_ENABLE_SERIAL=ON \
  -DELMXX_DIAG_BINARY_PATH=$E/2x1_brazil.I1850ELM.*/run/elm_diagnostics.bin \
  -DELMXX_LAKE_DIAG_BINARY_PATH=$E/brazil_with_lakes.*/run/elm_diagnostics.bin
cmake --build <build> -j 8 && ctest --test-dir <build>
```

- **`-DKokkos_ENABLE_OPENMP=OFF` is required on this machine.** Kokkos defaults
  to the OpenMP backend and Apple clang has no OpenMP, so configure fails in
  `kokkos_tpls.cmake`. This matches what EKAT does inside the E3SM build.
- Both diag paths must be set or the data-driven tests **skip and report
  `Passed`** (see the note further down).
- **The suite is currently RED: seven pre-existing failures** at submodule
  commit `2dcf189`, the one E3SM builds against — `CanopyTemperatureTests`,
  `test_cantemp_val`, `test_bgflux_val`, `test_canflx_val`, `test_soilflx_val`,
  `test_surfruninf_val`, `test_integration`. Baseline before blaming your
  change: `git stash`, rebuild, re-run.

Two things bite every time:

- **The ELM twins' only restart is at night.** Day 6 00Z is 20:00 local at
  60°W, so every albedo field in it is zero. Re-run the twin with
  `STOP_OPTION=nsteps STOP_N=272 REST_OPTION=nsteps REST_N=8` for 34 restarts
  across the diurnal cycle, then restore `ndays`/5. Back up
  `elm_diagnostics.bin` first if that case feeds the C++ validation tests.
- **`ncdump` writes an unwritten element as a bare `_`.** A regex that scans
  for numbers does not skip it, it **shifts every value after it** — and an ELM
  restart is full of them (`H2OSOI_LIQ` is `(column, levtot)`, and a snow-free
  column never writes its five snow slots). Split the CDL value list on commas
  and map `_` to NaN.

## Standing traps — check these against any new work

Every one of these was found the expensive way, more than once in some cases.

**A zero is a legal number.** An unseeded view reads zero and produces
plausible garbage several kernels downstream — `z0mr`, `displar`, `dleaf`,
`qflx_top_soil`, `qflx_snow_melt`, and PAR in a validation test that then
silently solved the night problem. Ask of every parameter a kernel divides by
whether anything sets it.

**A single-rank, single-cell run hides a whole class of defect.** Two were
found the first time ELMxx ran on >1 cell or >1 rank (2026-10-04): four
namelist variables (`elmxx_do_albedo`, `elmxx_do_photosynthesis`,
`elmxx_stomata_closed`, `elmxx_co2_ppmv`) were never `mpi_bcast`, so every
rank but 0 ran with albedo and photosynthesis off; and there was ONE topounit
per rank, so every cell got the rank's first cell's precipitation in its water
balance. Every new namelist member needs a broadcast line, and any per-rank
scalar is a per-cell array in disguise. Check new work on the 4-cell, 2-rank
grid ("Restarts"): 1 vs 2 ranks must be bit-identical.

**`minval`/`maxval` skip NaN.** Every range report carries a non-finite count
for this reason.

**A constant a port writes for itself is a constant it can silently disagree
about.** Three times: Stefan-Boltzmann (`5.670374419e-8`, not `5.67e-8`), the
gas constant (`AVOGAD*BOLTZ = 8314.467591`, not `8314.4598`), and `fwet`'s
exponent — where ELM's own literal `0.666666666666_r8` is *truncated* and
writing the more accurate `2.0/3.0` broke the comparison. Derive constants from
the same file the model reads, and match ELM's literals even when they are
worse.

**Relative tolerance is wrong for offset or near-zero quantities.** A soil
temperature is ~274 K, so 1e-7 relative permits 2.7e-5 K while the deep-layer
signal is ~1e-6 K. `errsoi` spans O(1) to 1e-10 within one timestep; `h2ocan`
comes out at 5e-18 against an exact 0. Assert absolutely, in the quantity's own
units, and report the measured worst so the bound stays honest.

**Some expressions are ill-conditioned and no tolerance will fix them.**
`int_snow` divides by `0.5*(cos(rpi*x)+1)` with `cos -> -1`; measured, that
denominator is 2.06e-7 and amplifies argument error by 4.85e+06. One ULP
between Fortran's `**` and `Kokkos::pow` becomes 5e-10. Write the arithmetic
down next to the loosened bound.

**`snl == 0` does not mean "no snow".** A pack too thin for a layer is carried
as `h2osno` on a layerless column and handled by a separate thin-snow branch.
Splitting those from genuinely bare columns is what made the cold SoilTemperature
bug legible.

**ELM never clears inactive layers.** It recomputes `tk`, `fn` and friends only
for `j >= snl+1`; slots above hold whatever they had when that layer last
existed — 0.86 W/m² of "flux" at a step with no snow layers. Start every
layer comparison at `NLEVSNO + snl` or you are grading ELM's bookkeeping.

**The layer below is not the slot below.** ELMxx inserts a standing-surface-water
node between snow and soil where ELM has none. Five separate places took `j+1`
to mean "the next layer"; `SoilTempCon::below()` now exists so they share one
rule. Suspect it anywhere snow meets soil.

**`ncdump` prints an unwritten element as a bare `_`.** A number-scanning regex
does not skip it, it SHIFTS every value after it — soil layer 1 reads as layer
6 with a plausible number. Split on the comma and map `_` to NaN.

**A diagnostic binary can be stale in a way that looks like working.** Count
its labels: a current one has ~1006. `1x1_glc`'s has 956 and predates the
Photosynthesis instrumentation; `1x1_brazil`'s once had 18.

**`cmake -B` with a relative `-S .` silently reconfigures nothing** if the
shell's directory is not what you assume. Check `CMakeCache.txt` actually
changed before believing a result — a whole "suite is green against the snow
case" run turned out to be the wrong binary.

**`case.build` needs `export SDKROOT=$(xcrun --show-sdk-path)` on this laptop.**
Without it the PIO configure fails with `ld: library 'System' not found` —
before any ELMxx code is compiled, so it looks unrelated to your change.

**A new source file needs `./case.build --clean-all`**; `--clean lnd` is not
enough and fails on a missing `.mod`. So does any `env_build` change.

**Running `e3sm.exe` directly does not regenerate namelists.** After
`./xmlchange STOP_N/STOP_OPTION/CONTINUE_RUN/REST_*`, run
`./preview_namelists` -- otherwise `drv_in` keeps the old values and the
"continue" silently runs from the cold start (happened 2026-10-04: a
two-month diagnostics run became a full year).

**A data-driven test can skip because its binary STARTS late.** The lake
replays read topology at "step 1"; a binary recorded from a restart has no
step 1, the suite skipped, and gtest counted it passed. Read topology from
`reader.timesteps().front()`, and assert the step you validate has records.
Fixed-step replays (steps 10/50/100/200) also never reach a seasonal branch:
sweep a year (`ELMXX_LAKE_SWEEP=all`) before calling a kernel validated.

**`cime` can silently leave its branch.** Found detached at upstream master with
`uses_kokkos()`'s elmxx case absent. Parent `git status` shows `? cime` instead
of ` M cime`; confirm with `grep -n elmxx cime/CIME/build.py`.

## Build gotchas — durable, hard-won

### ELMxx's C++ under E3SM

- **Kokkos comes from EKAT, not ELMxx's submodule.** ELMxx's `CMakeLists.txt`
  detects an existing `Kokkos::kokkos` target and skips its bundled copy
  (`ELMXX_STANDALONE` guard). Both are the same `E3SM-Project/kokkos` fork at
  4.5.1, two commits apart, so either satisfies ELMxx identically. `USE_KOKKOS`
  is switched on for elmxx cases by `uses_kokkos()` in `cime/CIME/build.py`.
- **Keep the `externals/kokkos` submodule.** It is what makes the standalone
  build and the CTest suite work (googletest is vendored at
  `externals/kokkos/tpls/gtest`), and it costs nothing when unused.
- **Never use `CMAKE_SOURCE_DIR` inside ELMxx's CMake.** Under E3SM's
  `add_subdirectory` it resolves to *E3SM's* root. Use `PROJECT_SOURCE_DIR` /
  `PROJECT_BINARY_DIR`, which are correct in both standalone and nested builds
  because ELMxx calls `project()`.
- **Do not add ELMxx's Fortran bindings to `Filepath`.** They are compiled into
  `libelmxx_fortran.a`; listing them as well double-compiles them and produces
  duplicate `.mod` files.

### macOS + Homebrew GCC: the `_Static_assert` trap

macOS SDK 26's `mach/port.h` defines
`xnu_static_assert_struct_size(...)` as `_Static_assert(...)`, and
`mach/message.h` invokes it at file scope. **`_Static_assert` is C11**: clang
accepts it in C++ as an extension, GCC's C++ frontend does not. Any C++
translation unit reaching a mach header fails to parse -- here that is
googletest's `gtest-all.cc` (it includes `mach/mach_init.h` for timing), which
alone made the entire test suite unbuildable.

Fixed in ELMxx's `if(APPLE)` block with
`add_compile_definitions($<$<COMPILE_LANGUAGE:CXX>:_Static_assert=static_assert>)`,
guarded on `CMAKE_CXX_COMPILER_ID STREQUAL "GNU"`. C is deliberately untouched --
`_Static_assert` is standard there.

Forcing an older SDK is **not** a workaround: `-isysroot` to MacOSX14/15 moves the
failure into GCC's own fixed headers, since Homebrew GCC 14.3.0 was built against
the newer SDK. `urbanxx/build_gcc.sh` does not solve this either -- its `mpicxx`
resolves to the same `g++-14`, and its own `gtest-all.cc.o` was never built.

### Running ELMxx's test suite — always set the diag paths

With `ELMXX_DIAG_BINARY_PATH` unset, the 17 data-driven validation tests **skip at
runtime and report `Passed`**, so the suite looks green when nothing was
validated. Always configure both paths:

```bash
E=~/projects/e3sm/elmpy/e3sm/cime/scripts
cmake -B <build> -S . \
  -DELMXX_DIAG_BINARY_PATH=$E/2x1_brazil.I1850ELM.PNNL-L07D666226.gnu11/run/elm_diagnostics.bin \
  -DELMXX_LAKE_DIAG_BINARY_PATH=$E/brazil_with_lakes.1x1_brazil.I1850ELM.PNNL-L07D666226.gnu11/run/elm_diagnostics.bin
```

Configure log should say "Validation tests enabled: <path>", not "enabled with
empty path". Test count goes 28 -> 30 once the paths are set.
