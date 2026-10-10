# ATHAD_PERID — the atmosphere a peridotitic magma ocean left behind

An atmospheric general-circulation model of the Hadean epoch that follows
[ATHAD_COND](../ATHAD_COND): the magma ocean has solidified, the steam has long since rained
out, and the volatiles a **peridotitic** (ultramafic, mantle-composition) melt outgassed have
cooled into a **CO₂-dominated atmosphere of 5–20 bar over a 20–150 °C ocean**. Nominal run
point: **12.5 bar / 358.15 K (85 °C)**, the midpoint of both ranges.

**"Peridotitic" names where the volatiles came from, not the present surface.** The melt is
gone; what this model carries is its outgassed inventory, cooled. So the thermodynamic regime
is ATHAD_COND's — subcritical water, a liquid sea, condensation live from the surface up — and
NOT ATHAD's supercritical column. The fork is cut from ATHAD_COND for that reason and keeps its
invariants rather than inverting them.

**Forked from `ATHAD_COND` @ `3f0b05d` (2026-08-20), carrying its full history** — 63 commits
reaching back to ATHAD's own fork-from-ATOM commit — so fixes cherry-pick in every direction
(`git remote cond` and `git remote athad` point at the two sibling working copies). ATHAD_COND
came from `ATHAD` @ `29ca2f9`, ATHAD from `ATOM_Precipitation` @ `1e3f319`. There is no
paleogeography and no time-slice series — ATHAD_PERID is one epoch. There is an **ocean**, but
as a boundary condition only: no prognostic hydrosphere.

**Where it differs from ATHAD_COND:** the same physics at a fifth of the pressure and 155 K
colder, with a genuinely CO₂-dominated composition — 74 mole-% CO₂ against COND's 35 %, and
water down from 56 mole-% to 4.6 %. The sea still sets the surface humidity, but a much drier
atmosphere sits on it, and that difference is why several of ATHAD_COND's measured nulls
(the surface albedo especially) must be re-measured here rather than assumed.

**EVERYTHING IN THIS FILE BELOW THE COMPOSITION SECTION IS INHERITED FROM ATHAD_COND AND WAS
MEASURED ON IT** — a 60 bar / 513 K world — unless it says otherwise. The code is shared and
the reasoning travels; the numbers do not.

## Build, run, test

```bash
make perid                                      # -> cli/perid
make test                                       # IAPWS + column self-tests, run these first
cd python && OMP_NUM_THREADS=8 ../cli/perid config_perid.xml
```

**`atom_log.txt` is reserved by the model — never redirect a run into it.**
`lib/Utils.cpp:45` opens `atom_log.txt` in the run directory with `std::ofstream::out`,
which truncates it. A shell redirect to the same name loses the whole console output when
the model's handle closes; a 100-iteration run's printouts were destroyed this way. Use any
other filename.

Output lands in `python/output_Hadean_peridotite/` — the ATOM line's convention
(`output_<name>/`, underscore, relative to the run directory; the giant-planet
siblings hyphenate instead).

`make` regenerates the parameter bindings from `param.py` first. The generated files
(`atmosphere/*.inc`, `python/atmosphere_pxd.pxi`, `python/pyperid.pyx`,
`cli/config_perid.xml`, `python/config_perid.xml`) are **tracked on purpose**, so a
`param.py` change shows its full effect in the diff. Regenerate and commit them together.

Removing a parameter from `param.py` deletes a C++ member, so it must be removed
together with its uses or the build breaks.

## Atmospheric composition

Mole fractions are the input; the model works in mass fractions. **The configured mole fractions
are the SEA-SURFACE values.** Stated inputs: p = 5–20 bar, T = 20–150 °C, CO₂ 60–80 %,
N₂ 10–20 %, and SO₂, CH₄, NH₃, CO, H₂ each ~1 % by mole.

| Species | Mole frac. | Mass frac. |
|---|---|---|
| H₂O | 0.0463 | 0.021235 |
| CO₂ | 0.7417 | 0.831045 |
| N₂ | 0.1590 | — |
| SO₂ / CH₄ / NH₃ / CO / H₂ | 0.0106 each | — |
| background (N₂ + the five traces) | 0.2120 | 0.147719 |

- **M_mean = 39.278 g/mol**, **R_mix = 211.68 J/(kg·K)**
- **Background** (everything but H₂O and CO₂): M_bg = 27.369 g/mol, **R_bg = 303.79 J/(kg·K)**.
  This is what the `R_Air` parameter means — it is not air. `M_nonwater` = 40.311 g/mol is the
  different number that matters for saturation.
- **p_surf = 12.5 bar**, **T_surf = 358.15 K** (prescribed), **ρ_surf = 16.49 kg/m³**
- cp = 947.9, cv = 736.2 J/(kg·K); **scale height 7.73 km**, dry lapse **10.35 K/km**

**`x_H₂O` is not an input.** It is `p_sat(358.15 K)/p_0` = 0.5787/12.5 = 0.0463, set by the sea.
CO₂, N₂ and the five traces fill the remaining 95.37 % in the stated dry ratio
70 : 15 : 1 : 1 : 1 : 1 : 1, with CO₂ absorbing the rounding so the eight sum to 1.0000000000.
Verified: `q_sat` at the sea = 0.021233 against `q_H₂O` = 0.021235, so **the surface is saturated
by construction**; water boils at 462.96 K, leaving the sea 105 K below boiling.

**The five traces are ~1 % each and are NOT the inherited values.** The fork initially carried
ATHAD_COND's CH₄/NH₃/H₂/CO at 1.4 % simply because they came with the clone, which left the mole
fractions summing to **1.056** and `M_bg` 13 % adrift until the startup guard refused to run.
By mass the background is N₂ 0.768, SO₂ 0.117, CO 0.051, NH₃ 0.031, CH₄ 0.029, H₂ 0.004, so
**Σf_i·κ_i = 6.39e-04, 639× the lumped `kappa_bg`** — item 60's per-species split doing its work
on a third composition.

**The cold-trap floor is derived from this column, not inherited**: `c_h2o_dry_top` = 4.06e-04,
the `q_sat` minimum at level 24 (14.59 km, 254.0 K, 1.233 bar) where the profile goes isothermal.
ATHAD_COND's 0.00452 is **11.1× wetter than this column can hold**.

Only **H₂O (`c`) and CO₂ (`co2`) are prognostic**, both as mass fractions. The rest is a fixed
well-mixed background entering `R_mix`, `cp_mix` and the opacity.

## Where the physics lives

| File | What it owns |
|---|---|
| `MixtureAtm.h` | Composition → mass fractions, `R_of`, `cp_of` (Shomate), `M_of`, `M_nonwater`, water's critical point |
| `SaturationH2O.h` | IAPWS saturation + sublimation, Watson `latentHeat`, exact `saturationMassFraction`, `dewPoint`, `satDerivFactor`, exact `dqSatdT`, **`moistLapse`** |
| `MultiLayerRadiation.h` | Grey optical depth from column mass with pressure broadening; surface energy balance |
| `ThermoAtm.h` | The column integration (`densities`), the sea saturation boundary (`waterVapourEvaporation`), diagnostics |
| `test/saturation_selftest.cpp` | IAPWS reference points |
| `test/cond_column_selftest.cpp` | The ATHAD_COND regime: both mixtures, the config's own arithmetic, the moist lapse (with an Earth control) |

## Four invariants — do not silently break these

1. **There is no topography, and the planet is hemispherically symmetric.** `h ≡ 0`,
   `i_topography ≡ 0`, `is_land()` false everywhere — now a global *ocean*, the same
   statement in code. `LandOceanFraction()` throws if a land point appears. Any
   north–south asymmetry is a defect: symmetric surface parabola, mirrored insolation, no
   obliquity, no seasons.

2. **Water is subcritical everywhere and condensation is LIVE.** *This is ATHAD's
   invariant 2, inverted, and it is why this fork exists.* Every path ATHAD made a genuine
   no-op is a path that now runs. Here the sea is at **358.15 K and 12.5 bar**, water boils
   at **462.96 K**, so the surface is 105 K below boiling — a wider margin than
   ATHAD_COND's 35 K, and the critical point (647.096 K / 220.64 bar) is far away in both
   coordinates. What carries over is which formulas are allowed: **never Magnus**
   (calibrated to ~320 K) and **never the dilute `q_sat = ep·E/(p−E)`**. The dilute form is
   *less* wrong here than in ATHAD_COND — water is 4.6 mole-% at the sea rather than 56 % —
   which makes it more dangerous, not less: it now looks right in the printouts almost
   everywhere. Use the exact mass-fraction form, as the rest of the tree does.

3. **Radiation runs in mode 2** (direct σT⁴). Modes 0/1/3/4/5 lean on the Scotese snapshot
   or the 280 ppm CO₂ reference; neither exists in a 74 mole-% CO₂ Hadean atmosphere.
   Radiation must *set* the profile rather than be handed one — **and it does not yet**:
   `densities()` re-imposes the adiabat every iteration, and the OLR currently equals
   `σT_lid⁴` exactly (README item 6).

4. **The column is on its own adiabat, integrated not fitted, and the adiabat is MOIST
   where the air is saturated.** `SaturationH2O::moistLapse`, hydrostatic on the layer-mean
   T, isothermal above `t_skin`. Do **not** restore the COSMO `T = T₀√(1−coeff·h)` form.
   Do **not** revert to `dT/dz = −g/cp`: `cosmo_lapse_fraction = 1.0` is justified in ATHAD
   *because nothing condenses there*, and that justification is void here.
   And do not use `g/(cp + L·dq_sat/dT)` — that is only the denominator of the moist-static-
   energy balance. The dry lapse here is **10.35 K/km**; the moist one is what the model
   integrates where the air is saturated.

## Assumptions vs. results — read this before quoting any number

| Parameter | Value | Status |
|---|---|---|
| `geothermal_flux` | 150 W/m² | **Inherited from ATHAD's magma ocean and almost certainly wrong.** 55 % of the whole energy input, against modern Earth's 0.087 W/m². Biggest unexamined input |
| `kappa_CO2` / `kappa_H2O` | 0.001 / 0.01 m²/kg | Factor-of-2 uncertain. `kappa_CO2` dominates (CO₂ column 5.8e5 kg/m² vs water 2.6e3) |
| `kappa_N2/CH4/NH3/H2/CO/SO2` | 1e-6 / 3e-3 / 1e-2 / 1e-5 / 1e-4 / 2e-3 m²/kg | ATHAD's values, and they now matter here: Σf_i·κ_i = **1.944e-3**, 1944× the lumped `kappa_bg` = 1e-6 that was exact while the background was pure N₂. Six assumptions of the same standing as the two above |
| trace-gas mole fractions | 0.014 each (CH₄/NH₃/H₂/CO/SO₂) | ATHAD's numbers, adopted 2026-08-20. **The fractions match, the inventory does not** (1.4 % of 250 bar ≠ 1.4 % of 60 bar), and they push dry CO₂ to 77.4 %, below the quoted 89–95 % |
| `t_surf_equator` / `t_surf_pole` | 513.15 / 503.15 K | **Prescribed, not solved** |
| `t_skin` | fixed-point iterate | Currently *sets* the OLR rather than following it. **Its fixed point is `sigma*T^4 = F/2` since 2026-08-21** (ATHAD item 67); it used to be `= F`, which made the budget close by construction. `ATM_SKIN_GREY=0` restores that. Effect size here UNMEASURED |
| `albedo_cloud` | 0.50 | IS the planetary albedo (0.4997 measured); saturates wherever condensate exists, which is now everywhere |
| `albedo_surface` | **0.06** | A **water** value since 2026-08-20; it was 0.08 with ATHAD's basaltic-melt justification attached. Clear-sky only, and the cloud bump overwrites it almost everywhere, so its reach is small — measured, not assumed |
| `omega` | 3.17e-4 (5.5 h day) | Inherited; this epoch is later and slower |
| `mc_M_max` | **100 kg/(m²s)** | Sized 2026-08-20 from the density ratio (Earth's 3.0 implies σw = 2.5 m/s; ×36–42 kg/m³ here) and confirmed against the uncapped scheme's own 13.0–13.7. Was Earth's 3.0 and **pinned from 2.2 to 21.9 km in 177 of 181 columns** |
| `delta_i_c` | 500 s | Bechtold convective timescale, Earth-calibrated. A time, not a pressure — no unit error, and no evidence here to replace it |

## What the model currently says

> **THE GREY SKIN FACTOR IS IN, AND IT IS ON BY DEFAULT (ported from ATHAD item 67, 2026-08-21).**
> `t_skin` was solving `sigma*T^4 = F` where `F` = absorbed SW + geothermal. That is the planet's
> EFFECTIVE EMISSION temperature; the top layer of a grey atmosphere is at `T_eff/2^(1/4)`, because
> it sees no downward flux and re-emits half of what passes up through it. **The precondition was
> checked in THIS tree, not assumed**: `SW_abs` appears in exactly one place, `MultiLayerRadiation`'s
> surface energy balance, and nowhere in the column, so the classical skin result `sigma*T^4 = F/2`
> is exact here rather than a fit. `MultiLayerRadiation` already contained the right relation twice.
>
> **Consequence: `t_skin` is now ~2^(1/4) colder, and the OLR follows it down.** Because the reported
> OLR descends onto `sigma*t_skin^4`, setting `sigma*t_skin^4 = F` made "OLR -> absorbed" an
> identity. In ATHAD, correcting it turned a -1.41 W/m2 imbalance at 200 iterations into **+123.76**.
>
> **EVERY OLR AND IMBALANCE FIGURE RECORDED BELOW PREDATES THIS AND WAS MEASURED WITH THE FACTOR
> OFF.** None is retracted -- they are correct measurements of a branch still reachable with
> `ATM_SKIN_GREY=0` -- but none describes the shipped model.
>
> **THE SIZE OF THE EFFECT IN THIS TREE IS NOT MEASURED.** ATHAD is 250 bar and this is a different
> column; the defect and its repair transfer, the magnitude does not. Run the A/B before quoting a
> number here.


> **THE 2026-08-24 PORT FROM ATHAD (its README items 68-80). FOUR DEFAULTS CHANGED, AND THE
> A/B IS RUN.** Every number recorded in this file and the README predates the port; they are
> not retracted -- each of the four is restorable with one environment variable -- but none of
> them describes the shipped model. **Measured, 40 iterations, 24 threads, env-only arms: OLR
> 136.98 -> 136.28 W/m2 (-0.5 %), `Psi_max` -27 %, `div(rho u)/rho` rms +62 %, `P_snow` 0 ->
> 9.46e4 mm/a.** Three things the README's port entry says at length and this table cannot:
> **`ATM_SAT_SUPERHEAT` is INERT here** (reached 1.67 million times, rejects zero, zero
> forbidden cells at all four census stages -- ATHAD's largest effect is a no-op in a
> subcritical column); **the precipitation flux is the `P_max_flux` cap by a factor of 1.6e5**,
> so no precip number here means anything yet; and **this fork is the control that says the
> divergence rise is not radiative** -- its OLR moved 0.5 % while `div` rose 62 %.
>
> | knob | was | is | ATHAD's finding |
> |---|---|---|---|
> | `ATM_PROJ_SWEEPS` | 1 | **10** | the initial pressure projection was under-converged: `Psi(ground)`, which must be 0, was **2.09x the interior circulation**, and 10 sweeps take 52.5 % of that off. A knee, not a cure (item 68) |
> | `ATM_RAD_DIRECT` | absent | **on** | the closed-form two-stream equilibrium solve, exact in two O(N) passes, replaces the 4-sweep Lambda iteration. `n_lambda = 4` is an Earth constant: **15 % high** on ATHAD's column, and the exact answer costs nothing (items 30, 71). `ATM_N_LAMBDA` is new alongside it |
> | `ATM_SAT_SUPERHEAT` | absent | **on** | `clampAndFade` tested whether a cell could hold a condensed phase, condensed, added the latent heat that makes the answer FALSE, and never re-tested. The ice scheme deleted the result every iteration. **The largest single effect ATHAD has recorded: OLR -49.4 %** (item 75) |
> | `ATM_PRECIP_BANDS` | absent | **on** | the precipitation phase bands were Earth's, and they did two jobs at once: `(band) ? (inherited + produced) : 0` DESTROYS a flux arriving from above rather than only withholding production. Snow's -20 C floor goes; the inherited flux always passes (item 76) |
>
> **Why they ship on rather than off.** Each is a repair whose ARGUMENT is universal -- a
> projection that has not converged, a relaxation that has not converged, a routine disagreeing
> with `IceSchemeCommon::canCondense` about where a phase can exist, a flux destroyed by the
> temperature of the air it falls through. **The SIZES were ATHAD's and did not transfer, which
> is now measured rather than warned about**: one of the four turned out to be inert here and
> another turned out to be capped. Restore the old branch with
> `ATM_PROJ_SWEEPS=1 ATM_RAD_DIRECT=0 ATM_SAT_SUPERHEAT=0 ATM_PRECIP_BANDS=0`.
> **The one-knob decomposition was run in ATHAD_COND (2026-08-24)**: the OLR there is
> `ATM_PRECIP_BANDS` (+14.4 %), `ATM_RAD_DIRECT` +1.6 %, `ATM_PROJ_SWEEPS` and
> `ATM_SAT_SUPERHEAT` 0.00, and the divergence rise is entirely `ATM_PROJ_SWEEPS`.
> **It has NOT been repeated in this tree**, whose all-knob OLR moved only -0.5 %.
>
> **Default-off knobs ported with it, off-branch bit-identical:** `ATM_CELL_ALTERNATE` (item 69 --
> four of five prescribed cells turn the same way, and `edgeRadialCoeff` already assumes the
> parity it does not impose), `ATM_METRIC_EXACT` (item 80 -- `exp_rm = 1/(rm+1)` is a QUADRATIC
> stretch's Jacobian applied to an exponential grid; the repair also adds the `-(J'/J)f'`
> curvature term the Laplacian omits in both metrics. **Default off because in ATHAD the correct
> Jacobian made the projection WORSE**), `ATM_PRECIP_CAP` (item 77). Print-only: the `ATM_ICE_CENSUS`
> stages, the `canCondense`-violation counter, the `P_rain` cap probe and the `S_r` decomposition
> (items 74, 77, 78), and the zonal ParaView writer's **true-height vertical axis** (item 70 --
> it was level index, which on this stretched grid distorts by ~18x and varies with altitude, so
> no aspect setting could undo it; no physics reads those fields).
>
> **Deliberately NOT ported.** `ATM_SKIN_TAU` (item 73): measured in ATHAD and **withdrawn** --
> `T_rad(tau)` exceeds the prescribed profile at all 41 levels there, so there is no
> radiative-convective crossing to switch at. It also presumes ATHAD's isothermal lid over a
> prescribed dry adiabat and its `ATM_PROGNOSTIC_T` knob, and this tree has neither: its adiabat
> is moist and state-dependent. `moist_phys_start_iter = 0` (item 74) is already the default here
> and has been since the fork.

**One 40-iteration baseline, and its headline is a limitation.** OLR 270.91 W/m² = absorbed
SW + geothermal 270.91 exactly; photosphere at 19.2 km and **262.91 K = `t_skin`**; **emission
from the isothermal skin 100.0 % of columns**. The OLR carries no radiative information at all
in this configuration — it is arithmetic. Across the family the collapse is monotone in
pressure: ATHAD 1.1 % of columns, ATHAD_COND 58.6 %, this fork **100 %**, because at 12.5 bar
with 2 % water there is no optical depth to put the radiating level below the isothermal top.
**Do not test opacity, composition or clouds against the OLR here until the profile is freed.**

`Psi_max` = 9389 (1e9 kg/s), a quarter of ATHAD_COND's and stable to −0.08 % over 20 iterations;
mean albedo 0.4994 = `albedo_cloud`, so the deck saturates here as in both siblings.

The run that follows is `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1` — ATHAD item 66 showed that pairing
converges rather than diverging, and this fork has the most to gain from it.
**HALF OF THAT PAIR NOW EXISTS HERE: `ATM_RAD_DIRECT` came in with the 2026-08-24 port and is
the default. `ATM_PROGNOSTIC_T` STILL DOES NOT EXIST IN THIS TREE** — it is ATHAD's knob on
ATHAD's prescribed dry adiabat, and this fork's adiabat is moist and state-dependent, so the
run named above cannot be launched as written until that half is ported or rewritten.

## Relationship to the family

Siblings live beside this directory. **Three Hadean models now, in time order:**
`ATHAD` (~4.4 Ga magma ocean, 250 bar of supercritical steam over 1500 K), `ATHAD_COND`
(the steam rained out: 60 bar over a 513 K sea), and this one (**12.5 bar over a 358 K
ocean**, the outgassed inventory cooled further). Also `ATOM_Precipitation` (modern Earth),
`ATJUP`, `ATSAT`, `ATURAN`, `ATNEPT` (giants) and `ASTIM` (impacts).

**The shell is renamed and nothing else:** `libperid.a`, `cli/perid`, `config_perid.xml`,
`pyperid`, `test/perid_column_selftest`. C++ class, file and function names stay identical
across the family — that is what makes `git cherry-pick` work between the three Hadean
trees, and `git remote cond` / `git remote athad` point at the sibling working copies.

Preserve that. The git history is ATHAD's, so a fix made in any of the three can travel to
the other two — use it rather than re-implementing something that already exists next door.

**The pattern to expect, and why it is worse here than in ATHAD.** ATHAD found seventeen
defects of one shape: *Earth's numbers as bare literals inside physics kernels.* At 250 bar
those are wrong by orders of magnitude and obvious. **At 60 bar they are wrong by a factor,
and plausible.** `p_stat <= 1000 hPa` fails at every level of ATHAD's column and so does
nothing; here it selects everything above 53 km and fires the deep-convection scheme in the
stratosphere. When something behaves oddly, look for a constant that was true at 1 bar and
288 K — and do not assume that a branch which now *runs* is therefore correct.

**⚠️ DO NOT CHERRY-PICK ATHAD'S BALANCED INITIAL STATE AS IT WAS FIRST WRITTEN.** ATHAD
items 26-28 build a `initBalancedState` that this model does not have. Item 27's version
balances the **θ-momentum equation only**, and porting it here would reproduce a failure
ATHAD measured in full: it drove `max |u|` from 0.114 to **11.17 m/s** over 200 iterations,
reversed the tropics to sinking over the hottest surface on the planet, and the meridional
streamfunction — **which is built from `v` alone and therefore cannot see a radial failure**
— reported the circulation as healthy throughout. The streamfunction here was just repaired
(`a85122c`, the density inside the integral); that fixes what it measures, not what it is
blind to.

The reason it fails is a **switch**, and this model has the identical one. At `u = v = 0`:

| | ATHAD | here |
|---|---|---|
| `rhs_u` | `RHS_Atm_Turb.cpp:1001` | line 1001, identical |
| `AtomUtils::coriolis_nontraditional()` (gates `coriolis_rad`) | default **false** | `lib/Utils.h:52`, default **false** |
| `AtomUtils::metric_curvature()` (gates the `−(v²+w²)/r` term) | default **false** | `lib/Utils.h:90`, default **false** |
| `buoyancy_ramp` at iteration 0 | **0** (`buoyancy_ramp_iters` = 300) | **0**, same 300 |
| `ATM_METRIC_RADIUS` | default on | default on |

So the radial momentum equation at iteration 0 is `rhs_u = −dp_dyn/dr·exp_rm` **and nothing
else**, and any radial gradient a balance writes is a pure unopposed vertical force. The
same switch governs the `w²cotθ` term item 27's θ-balance is *built on*, so with the
defaults that balance is also computed against a force the model does not apply — though
that part is worth only 0.3 % of `F_θ` (both balances are ~99.7 % Coriolis).

**And turning the switch on does not fix it**, which was worth testing before assuming:
with `ATOM_METRIC_CURVATURE=1` the radial equation gains a force of 0.075 rms against the
~293 a θ-only balance creates — **~4000× too small**. The θ-force's shape varies with height
because `w` does, and no radial force of comparable size exists to pay for that variation.
Completing the metric does not rescue a one-component balance; only balancing both does.

**The lesson generalises past the balance, and is the reason this note is here rather than
in a commit message: a term written in `RHS_Atm_Turb.cpp` is not necessarily a term the
model applies.** Read the switch, not just the source. This is the family's constant trap
one level up — not an Earth number in a kernel, an Earth-tested code path that never runs.

If the balance is ported later, port **ATHAD item 28's two-component version** (it minimises
the unbalanced acceleration in both equations at once and reads `F_r`/`F_θ` switch by switch,
so it follows whatever this model is configured to do) **together with its residual
diagnostic**, which prints what each component is left holding. There is no version of "port
item 27 now and fix it after". Whether a balanced initial state helps here **at all** is a
separate question needing its own A/B in this regime — 120 km shell, saturated sea — and
ATHAD's own 200-iteration confirmation was still running when this was written (2026-08-15).

Fixes made here that belong upstream in ATHAD: the exact `dqSatdT` (ATHAD's dilute form is
correct only in the `x → 0` limit, and ATHAD's water is not dilute either — it is the bulk
gas, so the error there is larger, not smaller); `moistLapse`; and the `M_bg`/`M_nonwater`
distinction at the saturation call sites.

Traps already solved elsewhere in the family — check before re-deriving:
Coriolis/centrifugal signs (ATURAN `8b284cb`, `4201957`; ATNEPT `024c37f`, `e412b1b`);
mass- not mole-weighted mixture properties (ATNEPT `c116d71`); in-place Gauss–Seidel as a
threading defect (ATURAN `ffd0e0e`); report failures and limits in the README (ATURAN
`74b4ded`, ATNEPT `34286b8`).

**THE BUILD CAN SHIP A STALE OBJECT, AND THE `-MMD` REPAIR DOES NOT COVER IT** (ported
2026-08-26 from ATOM_Precipitation, where it actually bit). `-MMD` tracks headers only for
objects ALREADY compiled with it; one built before the flag has no `.d`, so make sees only its
`.cpp` and never rebuilds it. Upstream's `cli/atm.o` sat four weeks stale through a large port,
and because it holds the model as a **stack local**, adding members to the class made `main`
reserve the OLD `sizeof` while the library constructed the NEW one -- the constructor ran off
the end of main's frame onto the stack canary. Every run aborted with
`*** stack smashing detected ***` AFTER completing successfully, with no compiler warning, and
**AddressSanitizer saw nothing** (an ODR size mismatch is not an out-of-bounds access to
anything it tracks). Found by disassembling `main`. The cure is one line -- **every object
depends on the Makefile itself** -- and it is PROPHYLACTIC here: every C++ object in this tree
already had its `.d`, checked before applying.

**The constant-density streamfunction claim in ATHAD's CLAUDE.md was already out of date about
this tree** -- the density-weighted version is present here, and ATOM_Precipitation got it on
2026-08-26 (2.28x overweight at its cell core, `Psi_max` 851.68 -> 373.84 (1e9 kg/s), location
unmoved). Nothing to do here; recorded so it is not re-derived.

**The longitudinal ParaView slice was on level index here too**, fixed the same day: item 70
put the ZONAL writer on true height and fixed only that one, so every longitude-height figure
kept the distortion. Plot-only; any conclusion drawn by eye from an older longitudinal plot
should be re-examined.

## Open risks

- **THE CONVECTION SCHEME IS A MOIST CONVECTIVE ADJUSTMENT SINCE 2026-10-10; TIEDTKE IS NOT CALLED**
  (README, first section; the user's instruction). `ATM_CONV_SCHEME=tiedtke` restores it, and
  every number recorded before that date needs it. `adjust` adds no routine: `densities()`
  already rebuilds `t` on the moist adiabat every moist step (invariant 4) and
  `SaturationAdjustment` + the ice scheme condense and rain; `adjust` stops calling
  `MoistConvection` and zeroes `MC_t`, `MC_q`, `MC_v`, `MC_w`, `P_conv` and the updraft / downdraft
  fields. Measured here, 20 iterations, 8 threads: 127 of 180 printed extrema identical, the
  scheme's own fields go to zero, everything else moves in the 6th-9th digit. **ATHAD_COND's
  findings about the mass-flux scheme (no rain, the `iter_prec` parity cycle) were NOT repeated
  in this tree.** **The long run, 400 iterations, three arms (README; noise floor zero): this is
  the tree where the scheme mattered most** — max T −0.107 K, up to 1.9 K in 109 cells near
  1.9 bar, max vapour −0.41 %, mean-water drift +0.8855 % (`tiedtke`) against +1.0141 %
  (`adjust`); winds unchanged to 1e-6; the vapour difference is still growing. (This bullet is this tree's own; the rest of the section
  is ATHAD_COND's, as the note below says.)

> **THIS SECTION WAS COPIED FROM `ATHAD_COND` AND IS BYTE-IDENTICAL TO IT** (diffed 2026-08-29:
> `diff <(sed -n '344,$p' ATHAD_COND/CLAUDE.md) <(sed -n '353,$p' ATHAD_PERID/CLAUDE.md)` is
> empty). Every measured number below was taken in ATHAD_COND. **None of them has been
> re-measured here**, and this fork is a different atmosphere — 12.5 bar CO2-dominated at
> 358.15 K, `M_mean` 39.278 g/mol, `R_mix` 211.68 — so they are a PRIOR inherited from the
> parent fork, not results.
>
> **What DOES carry, and why.** The grid is identical: `im` = 61, `zeta` = 3.0,
> `L_atm` = 6287.5 in both. So the two GRID-DERIVED facts transfer exactly — the
> `checkRadialMetric` spread of **12.29x** (it depends only on `zeta` and `im`) and
> `drag_n_layers = 5` meaning **1 786 m** of physical depth. Those two can be read as this
> fork's own.
>
> **What does NOT carry**: every physics measurement. The photosphere at 59.7 km / 278.6 K,
> `tau_above` 2.55e6, `tau_layer` 1.11e5, `N2` 4.30e-5, `Psi` 4.2212365e13, the whole `ubud_*`
> table (`ubud_pgf` 0.100816 and the rest), the OLR 283.81 against sigma*T_lid^4 = 270.62 and
> its -13.03 W/m2 imbalance, and `ATM_MC_GEOPOTENTIAL`'s `max c_u` 3.116e-03 at 24 854 m with
> its 20 km convective column. All of those are ATHAD_COND's atmosphere.
>
> **Three sentences below are self-referential and read as nonsense here**: "Unverified in this
> fork — no COND run has been made against the renamed writers" (this is not COND); "the reason
> it works here and not in ATHAD is geometry: **this fork's** convective column is 20 km deep";
> and "treat **COND's own** 0.004-0.03 % OLR deltas as inside the envelope until the same null
> pair is run here".
>
> The text is kept rather than deleted because a fork cut from ATHAD_COND inherits its structure
> and these remain the right things to look for. **Re-measure before quoting any of it.**


- **Field names changed 2026-08-20** (ported from ATHAD README item 61). The Brunt-Vaisala
  frequency squared is `brunt_N2` in the code and `BruntVaisala_N2` in the VTK files — it used to
  be `N2`, in the same file as nitrogen's `q_N2`. The eight species are written under bare names
  (`H2O CO2 N2 CH4 NH3 H2 CO SO2`), all mass fractions from `AtmMixture::split()`, so what is
  plotted is what is integrated; the raw transported CO2 array is `CO2_tracer`. **Unverified in
  this fork — no COND run has been made against the renamed writers.**
- **`ATM_MC_GEOPOTENTIAL` PAYS HERE — measured 2026-08-20, and it is the first updraft
  condensation in either fork.** `max c_u` 0 → **3.116e-03 g/kg/s at 24 854 m** (the 7.19e-04 first measured was at the old `M_max` = 3.0, which suppressed the mass flux 4.3× — superseded), with `max s_u`
  moving from the cloud base (661 m) to the top of the column (24 854 m). The reason it works here
  and not in ATHAD is geometry: this fork's convective column is **20 km deep** (base 1786–2200 m,
  LFS 20.5–28.1 km, 47 % of columns convecting) because its triggers have been fractions of surface
  pressure since the fork was cut, while ATHAD's updraft is one grid level deep and has no ascent
  for the geopotential to act on. Nothing integrated moves — OLR, photosphere, Ψ and the water
  column are identical — which is the albedo wall, not a null result about the physics.
  **ON BY DEFAULT here since 2026-08-20** (`ATM_MC_GEOPOTENTIAL=0` restores the old behaviour);
  ATHAD keeps it off, because its updraft is one grid level deep. `M_max` is sized now (100, see
  the table above), so the mass flux is no longer a cap.
- **The 40-iteration OLR is reproducible to only ~0.5 % across a thread-count change in ATHAD**
  (item 61), and this fork shares the convective trigger set that amplifies it. Ψ, albedo and the
  photosphere are unmoved at 1e-6. **State the thread count with any OLR number**, and treat COND's
  own 0.004-0.03 % OLR deltas as inside the envelope until the same null pair is run here.

- **The dynamics and the radiation do not agree where the levels are — measured here now,
  not inherited as an estimate** (ported from ATHAD README item 39, 2026-08-18).
  `exp_rm = 1/(rm+1)` is documented in two places as the Jacobian of the radial coordinate
  transformation (`TurbulenceAtm.h`; `PressureSolverAtm.h` writes it as
  `dp/dr_physical = exp_rm * dp/d(rad.z)`). It is the Jacobian of a **quadratic** stretch,
  `z ~ (rm+1)²/2`, while `init_layer_heights` builds an **exponential** one. `checkRadialMetric()`
  now prints the spread at every startup, and it is **12.29×** here: the core's radial length
  unit is 9 837 m at the surface and 120 850 m at the top, while `get_layer_height()` — which
  the radiation, `ConvectiveAdjustment` and every diagnostic use — has the true heights.
  CLAUDE.md previously carried "~12×" as an inference from `zeta = 3.0`; that inference was
  right, and it is now a measurement. The test is unit-free
  (`ratio(i) = [inv_2dr·exp_rm]/[1/(z[i+1]−z[i−1])]` must be constant in `i`), so it cannot be
  argued away as a units convention. **`im` is not the lever**: the spread depends only on
  `zeta` and `im`, and moving `im` 61 → 41 changes it 12.29 → 11.77. `zeta = ln(1.5) = 0.405`
  would make `exp_rm` correct by construction. `ATM_METRIC_CHECK=1` adds the per-level table.
  The diagnostic is print-only and self-silencing below a spread of 1.05.
- **The surface drag is Earth's, twice over, and its depth moves with the grid**
  (ported from ATHAD README item 44, 2026-08-18). `rayleigh_kf` and `drag_n_layers` were
  `constexpr` in `RHS_Atm_Turb.cpp` with comments justifying them by Earth's geography — a jet
  off the west coast of South America for the strength, an Andes orographic feature for the
  depth. Invariant 1 makes `is_land()` false everywhere here too. **`drag_n_layers = 5` is a
  count of AIR CELLS, not a length**, so its physical depth is whatever the grid makes it:
  236 m on `ATOM_Precipitation`'s grid, **1 786 m here**, 7 152 m in ATHAD. Both are config
  parameters now, bit-identically. **The part ATHAD's writeup did not have to state**: the
  depth moves with `im` as well as `zeta`, and ATHAD holds `im` fixed at 41 while this fork is
  at 61 — the same constant means 2 861 m at `im` 41, 1 786 m at 61 and 1 297 m at 81, a 2.2×
  swing from the level count alone. **A grid-refinement study here silently rescales the drag**,
  so `im` and `drag_n_layers` cannot be varied independently until the depth is reformulated as
  a length in metres. Not yet scanned here; ATHAD's four-arm scan is the template.
- **`Psi`, `tau_above`, `tau_layer` and `N2` are fields now, and `N2` VALIDATES the moist
  adiabat** (ported from ATHAD's item 42, 2026-08-18). All four were absent: `Psi` existed only
  as a CSV of zonal means plus the scalar `Psi_max`, and the other three were computed **nowhere
  in this model**. Measured at 20 iterations:

  | field | value |
  |---|---|
  | `tau_above` | 0 at the lid → **2.55e6** at the surface |
  | `tau_layer` | up to **1.11e5** per layer, at 7.7 km |
  | `N2` | **4.30e-5** 1/s² at the surface, 4.25e-4 at 114 km |
  | `Psi` | 4.2212365e13 kg/s, symmetric about the equator |

  **`N2` is the result worth having.** A column on the MOIST adiabat is DRY-stable, so the dry
  N² must be positive and equal to `(g/T)(Γ_dry − Γ_moist)`. With the self-test's own lapse
  rates that predicts **4.270e-5**, against a measured **4.30e-5** — agreement to **0.7 %**.
  So invariant 4 in its moist form is now *measured* rather than asserted, and the saturated
  sweep in `densities()` is integrating the lapse rate it claims to. Note this differs from
  ATHAD by design: there N² ≈ 0 through the column because the adiabat is dry.
  **`Psi` also caught a defect** — it read identically zero until the `write_meridional_
  streamfunction` / `print_min_max_atm` order was swapped (ATHAD's item 2, which had not been
  ported with the rest). The CSV was always right; the printed extrema and the field were not.
  **Not ported deliberately**: the `.vts` panorama carries none of these in either fork.
  **The photosphere is measured now too** (block ported 2026-08-18, 20 iterations, 24
  threads): **59.7 km, T there = 278.6 K, and 0.0 % of columns emitting from within 1 K of
  `t_skin`.** Against ATHAD's 239.4 km / 325.2 K with a skin fraction rising 1.1 → 39.2 %.
  **That difference matters more than the numbers.** ATHAD's `t_skin` pinning works by the
  radiating level migrating into the isothermal top, and here it has not: the photosphere sits
  at half the shell depth, 16 K above `t_skin`, in **zero** columns of the skin. So this fork's
  "the OLR is not independent of `t_skin`" risk cannot be the *same* mechanism as ATHAD's, and
  the OLR is correspondingly not sitting on it — 283.81 W/m² against σT_lid⁴ = 270.62 and
  absorbed = 270.78, a **−13.03 W/m² imbalance that is not identically the distance to a
  constant**. Whether it closes onto one is now checkable rather than assumable, which is what
  the instrument was for. Note the lid emissivity is still 0.0000, so the *lid* is transparent
  while the photosphere is deep — the two are decoupled here, unlike ATHAD at 230 km.
- **The radial momentum balance is the pressure gradient and nothing else, and the buoyancy
  body force is absent** (`ubud_*` ported from ATHAD's item 42, 2026-08-18, measured here at
  20 iterations). Ψ is built from `v` alone, so it cannot see a radial failure; `vbud_*`/`wbud_*`
  covered two components of three and this closes the set.

  | term | value | vs `ubud_pgf` |
  |---|---|---|
  | `ubud_pgf` | 0.100816 | — |
  | `ubud_cor` | **0.000000** | exactly zero |
  | `ubud_advv` | 0.002092 | 2.1 % |
  | `ubud_advh` | 0.002844 | 2.8 % |
  | `ubud_diff` | 0.000445 | 0.44 % |
  | `ubud_buoy` | **0.000000** | < 5e-7, so >2e5 down |

  Same structure as ATHAD's item 42 (pgf 1.15, cor exactly 0, advh 0.0085, buoy 1e-6) at ~10×
  smaller magnitude. Two of these are self-checks that now pass rather than being asserted:
  **`ubud_cor` is identically zero**, so `coriolis_nontraditional()` is genuinely off and not
  merely documented off — the CLAUDE.md warning that "a term written in `RHS_Atm_Turb.cpp` is
  not necessarily a term the model applies" is now instrumented rather than trusted. And
  **`ubud_buoy` is effectively zero**, which measures ATHAD item 34's *other* extra-`*dt` term —
  the half **not** repaired by the `surf_drag` fix in `a609d2b`. The buoyancy body force in a
  model whose density spans ~3 orders of magnitude is currently not acting.
- **The OLR is not independent of `t_skin`.** First thing to fix, and the same task as
  invariant 3: the profile is prescribed, so radiation cannot set it.
- **`geothermal_flux` is over half the energy budget** and is ATHAD's magma-ocean number.
  Nothing quantitative survives it being wrong.
- **The surface temperature is prescribed, not solved.** Every result is conditional on it.
- **Rain has nowhere to go.** The sea is a boundary condition, not a reservoir, so the
  water budget is not a closed-system test: `moist_phys_start_iter` is 0 here, not ATHAD’s 300, because condensation is the subject.
- **Boussinesq.** `ATM_ANELASTIC` is ported from ATHAD and ships default-off; it has not
  been measured in this regime. The density span here is ~3 orders of magnitude instead of
  ~5, which makes this the better testbed for whether it matters at all.
- **The column air mass is not conserved** (inherited): `p_stat.x[0]` is re-anchored every
  iteration to `r_air·R_mix·T_surf`.
- `time_start/end/step` remain because the time-slice loop is still structural, though only
  one slice ever runs.
