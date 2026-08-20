// Standalone self-test for the ATHAD_COND regime — the post-condensation Hadean column.
//
// ATHAD's numbers were checked once and then lived in CLAUDE.md; here they are asserted,
// because ATHAD_COND inverts the thermodynamic regime and almost every one of them moves.
// This test owns the composition arithmetic that the config values are derived FROM, so a
// change to either has to be made in both places and the disagreement shows up here.
//
// Build and run:
//   g++ -std=c++17 -Iatmosphere -Ilib -o test/cond_column_selftest test/cond_column_selftest.cpp
//   ./test/cond_column_selftest
//
// or from the project root:  make test-cond-column

#include "MixtureAtm.h"
#include "SaturationH2O.h"

#include <cstdio>
#include <cmath>

static int failures = 0;

static void check(const char* what, double got, double want, double rel_tol)
{
    const double rel = std::fabs(got - want) / std::fabs(want);
    const bool ok = (rel <= rel_tol);
    if (!ok) failures++;
    std::printf("  %-54s got %12.6g   want %12.6g   %s\n",
                what, got, want, ok ? "ok" : "FAIL");
}

static void check_true(const char* what, bool cond)
{
    if (!cond) failures++;
    std::printf("  %-54s %s\n", what, cond ? "ok" : "FAIL");
}

// The nominal run point. These four numbers are the model's inputs; everything else in
// this file is derived from them, and the config must agree.
static const double P_SURF_BAR = 12.5;
static const double T_SURF     = 358.15;      // 85 C, the midpoint of 20-150 C
static const double T_POLE     = 348.15;      // 75 C

// THE STATED DRY RATIO, which is what the surface mixture inherits: CO2 60-80 %, N2 10-20 %,
// and SO2 / CH4 / NH3 / CO / H2 each ~1 % by mole. Midpoints, as parts of a whole that the
// water then displaces from the top.
//
// Unlike ATHAD_COND's, these are NOT a literature dry-atmosphere composition quoted at some
// other height -- they are the outgassed inventory of a peridotitic melt, stated directly.
// So there is no "the composition and the temperature describe different heights" trap here;
// what the cold trap does to the water is a separate, computed thing (see c_h2o_dry_top).
static const double R_CO2      = 70.0;
static const double R_N2       = 15.0;
static const double R_TRACE    = 1.0;         // each of SO2, CH4, NH3, CO, H2
static const double R_SUM      = R_CO2 + R_N2 + 5.0 * R_TRACE;   // 90

// The SURFACE composition, copied from param.py. The checks below re-derive them from
// T_SURF and P_SURF_BAR and fail if the config has drifted from the physics it encodes.
static const double X_H2O_CFG  = 0.0463;
static const double X_CO2_CFG  = 0.7417;
static const double X_N2_CFG   = 0.1590;
static const double X_TRACE    = 0.0106;      // each of the five
static const double R_AIR_CFG  = 16.49;       // [kg/m3]
static const double T_SKIN     = 262.91;      // config t_skin — the isothermal top
static const double C_DRY_TOP  = 7.59e-4;     // c_h2o_dry_top: q_sat at the cold trap

int main()
{
    std::printf("\nATHAD_PERID column self-test  (%.1f bar, %.2f K)\n\n", P_SURF_BAR, T_SURF);

    // ------------------------------------------------------------------
    // 1. The regime itself. ATHAD's invariant 2 is that water is supercritical; this
    //    model exists because it is not. If these fail, nothing below means anything.
    // ------------------------------------------------------------------
    std::printf(" the regime is subcritical and the sea is liquid\n");
    check_true("T_surf below the critical temperature",
               T_SURF < AtmMixture::T_CRIT_H2O);
    check_true("p_surf below the critical pressure",
               P_SURF_BAR * 1000.0 < AtmMixture::P_CRIT_H2O);

    // Boiling point at the surface pressure: the sea must not be boiling, or there is no
    // liquid ocean to hold the surface layer at saturation.
    const double T_boil = SaturationH2O::dewPoint(P_SURF_BAR * 1000.0);
    check("boiling temperature at 12.5 bar",          T_boil, 462.96, 2e-3);
    check_true("the sea is below boiling at the equator",  T_SURF < T_boil);
    check_true("the sea is below boiling at the pole",     T_POLE < T_boil);
    // For the record, and because the 27 bar corner of the input range fails this:
    std::printf("      (the sea is %.1f K below boiling; ATHAD_COND's margin was 35 K)\n",
                T_boil - T_SURF);

    // ------------------------------------------------------------------
    // 2. The surface is SATURATED. This is the modelling decision that separates
    //    ATHAD_COND from a dry-surface variant: there is an ocean, so the air in contact
    //    with it carries water at its own vapour pressure — not at the dry-atmosphere
    //    0.4-2 %, which is what survives above the cold trap.
    // ------------------------------------------------------------------
    std::printf("\n the surface mole fraction is set by the sea, not by the dry composition\n");
    const double p_sat_surf = SaturationH2O::saturationPressure(T_SURF);        // [hPa]
    check("p_sat(T_surf) [bar]",                      p_sat_surf / 1000.0, 0.5787, 5e-3);

    const double x_H2O_surf = p_sat_surf / (P_SURF_BAR * 1000.0);
    check("saturated x_H2O at the surface",           x_H2O_surf, 0.0463, 5e-3);
    check_true("water is a minor constituent here, unlike in ATHAD_COND",
               x_H2O_surf < 0.10);

    // The rest of the surface column is the five trace gases at their fixed 1.4 % each,
    // and then CO2 and N2 in the dry ratio within whatever is left. The order matters and
    // is the physics: x_H2O is not free — it is p_sat(T_surf)/p_0, set by the sea — so
    // anything added to the mixture displaces CO2 and N2, never water.
    const double x_rest      = 1.0 - x_H2O_surf;
    const double x_CO2_surf   = x_rest * R_CO2   / R_SUM;
    const double x_N2_surf    = x_rest * R_N2    / R_SUM;
    const double x_trace_surf = x_rest * R_TRACE / R_SUM;
    check("x_CO2 at the surface",                     x_CO2_surf,   0.7417, 1e-2);
    check("x_N2  at the surface",                     x_N2_surf,    0.1590, 2e-2);
    check("each trace at the surface",                x_trace_surf, 0.0106, 2e-2);
    check_true("CO2 is inside the stated 60-80 %",    x_CO2_surf > 0.60 && x_CO2_surf < 0.80);
    check_true("N2 is inside the stated 10-20 %",     x_N2_surf  > 0.10 && x_N2_surf  < 0.20);

    // ...and the config carries exactly those, to the four decimals it stores them in.
    check("config x_H2O matches the saturated value", X_H2O_CFG, x_H2O_surf, 1e-3);
    check("config x_CO2 matches the derived value",   X_CO2_CFG, x_CO2_surf, 2e-3);
    check("config x_N2  matches the derived value",   X_N2_CFG,  x_N2_surf,  2e-3);
    check_true("config mole fractions sum to 1",
               std::fabs(X_H2O_CFG + X_CO2_CFG + X_N2_CFG + 5.0 * X_TRACE - 1.0) < 1e-9);

    // ------------------------------------------------------------------
    // 3. The two mixtures, from the model's own resolve(). The config's x_* and r_air are
    //    the SURFACE values; the dry mixture is what the column tends to aloft.
    // ------------------------------------------------------------------
    // Resolved from the CONFIG values, not the exact ones, so what this test reports is
    // what the model will actually run with.
    std::printf("\n mixture properties at the sea surface\n");
    const AtmMixture::Composition S =
        AtmMixture::resolve(X_H2O_CFG, X_CO2_CFG, X_N2_CFG,
                            X_TRACE, X_TRACE, X_TRACE, X_TRACE, X_TRACE);
    check_true("mole fractions sum to 1",             S.valid);
    check("M_mean [g/mol]",     S.M_mean * 1000.0,    39.278, 1e-2);
    check("R_mix [J/(kg K)]",   S.R_mix,              211.68, 1e-2);
    check("q_H2O [kg/kg]",      S.q_H2O,              0.0212, 2e-2);
    check_true("CO2 is the bulk of the mass here",    S.q_CO2 > 0.80);

    // The background is N2 plus the five ~1 % traces. N2 still dominates it by mass here --
    // 77 % -- unlike ATHAD and ATHAD_COND, where the traces carry ~71 % of the background,
    // because there the traces are 1.4 % against a much smaller N2 share.
    check("M_bg [g/mol]",       S.M_bg * 1000.0,      27.369, 1e-3);
    check("R_bg [J/(kg K)]",    S.R_bg,               303.79, 1e-3);
    check_true("N2 is the majority of the background by mass", S.f_bg[0] > 0.70);
    check_true("the six background fractions sum to 1",
               std::fabs(S.f_bg[0]+S.f_bg[1]+S.f_bg[2]+S.f_bg[3]+S.f_bg[4]+S.f_bg[5]-1.0) < 1e-9);

    const double cp_surf = AtmMixture::cp_of(S.q_H2O, S.q_CO2, T_SURF, S.M_bg);
    check("cp(T_surf) [J/(kg K)]", cp_surf,           947.9, 3e-2);

    // r_air is the column's anchor: ThermoAtm::densities() sets p_stat[0] = 1e-2*r_air*
    // R_mix*T_surf, so r_air and R_mix together ARE the 60 bar. Getting this wrong does not
    // fail loudly — it quietly runs a different atmosphere.
    const double rho_surf = P_SURF_BAR * 1.0e5 / (S.R_mix * T_SURF);
    check("rho_surf [kg/m3]",                         rho_surf, 16.49, 2e-2);
    check("config r_air reproduces it",               R_AIR_CFG, rho_surf, 1e-3);
    check("...and so the anchor gives back p_0 [hPa]",
          1e-2 * R_AIR_CFG * S.R_mix * T_SURF, P_SURF_BAR * 1000.0, 1e-3);
    const double H_surf = S.R_mix * T_SURF / 9.81;
    check("scale height [m]",   H_surf,               7728.0, 3e-2);

    const double M_nw_surf = AtmMixture::M_nonwater(S.q_H2O, S.q_CO2, S.M_bg);

    std::printf("\n the cold trap, and the water that survives it\n");
    // c_h2o_dry_top is NOT a literature figure here and NOT inherited: it is q_sat at the
    // coldest point of this model's own column. printColumnProfile puts the temperature
    // minimum at level 24 -- 14.59 km, 254.0 K, 1.2330 bar -- above which the profile is
    // isothermal, so that is where the water content is set for everything above.
    //
    // ATHAD_COND's value was 0.00452 and belongs to a 60 bar / 513 K column; carried across
    // unchanged it is 11.1x wetter than this column can hold. That fork's own README lists
    // "the cold trap is eight times too wet" as an open item, so the fork imported a known
    // defect -- which is exactly what a self-test is for.
    // The trap is level 23 -- 13.57 km, 263.7 K, 1.6082 bar -- one level BELOW where the
    // profile goes isothermal at t_skin, because q_sat keeps falling until the temperature
    // stops falling and then rises again as pressure drops. Its position moves with t_skin:
    // at the inherited 254.0 K it was level 24 and 4.06e-04, at the measured 262.91 it is
    // level 23 and 7.59e-04. That coupling is why both numbers are checked here.
    const double q_cold = SaturationH2O::saturationMassFractionAt(263.7, 1608.19, M_nw_surf);
    check("q_sat at the cold trap [kg/kg]",           q_cold, C_DRY_TOP, 2e-2);
    check_true("the trap is colder than the isothermal top is warm",
               263.7 - T_SKIN < 1.0 && 263.7 > T_SKIN);
    check_true("q_sat rises again above the trap",
               SaturationH2O::saturationMassFractionAt(T_SKIN, 1103.05, M_nw_surf) > q_cold);
    check("config c_h2o_dry_top matches it",          C_DRY_TOP, q_cold, 2e-2);
    check_true("and it is far drier than the sea",    q_cold < 0.05 * S.q_H2O);
    check_true("ATHAD_COND's 0.00452 would be far too wet here", 0.00452 > 5.0 * q_cold);
    check_true("and the inherited 254 K trap value would be too dry", 4.06e-4 < 0.7 * q_cold);

    // ------------------------------------------------------------------
    // 4. Why the dry adiabat cannot be used here.
    //
    //    ATHAD sets cosmo_lapse_fraction = 1.0 and integrates dT/dz = -g/cp, and that is
    //    correct there BECAUSE nothing condenses. The saturated lapse carries the extra
    //    heat capacity L*dq_sat/dT, and at 35 % water by mass that term is an order of
    //    magnitude larger than cp itself. This check is the justification for the moist
    //    adiabat in ThermoAtm::densities(); if it ever stops holding, that change is no
    //    longer needed.
    // ------------------------------------------------------------------
    std::printf("\n the saturated lapse rate is nothing like the dry one\n");
    // The "other" molar mass in the saturation formula is everything that is NOT water,
    // i.e. CO2 AND the background — M_nonwater, not M_bg. Passing M_bg here (N2 only,
    // 28.014) makes the sea surface come out at q_sat = 0.448 instead of 0.346, because
    // it weighs the non-water gas as if the 41 % CO2 were not there. That is a 29 % error
    // in the single quantity this whole model turns on, and it is the same confusion
    // CLAUDE.md records M_nonwater itself having had.
    check("M_nonwater at the surface [g/mol]",        M_nw_surf * 1000.0, 40.311, 1e-2);

    const double q_sat_surf = SaturationH2O::saturationMassFractionAt(
                                  T_SURF, P_SURF_BAR * 1000.0, M_nw_surf);
    check("q_sat at the sea surface [kg/kg]",         q_sat_surf, S.q_H2O, 5e-2);

    const double L_surf = SaturationH2O::latentHeat(T_SURF);
    check("latent heat at T_surf [J/kg]",             L_surf, 2.268e6, 5e-2);

    // NOTE: dqSatdT is the DILUTE Clausius-Clapeyron form, q_sat*L/(R_v T^2). At q_sat =
    // 0.35 it is not the derivative of the exact mass fraction — see the exact form below,
    // which differs by M_other/(x*M_H2O + (1-x)*M_other) ~ 1.45 here. Both are printed so
    // the gap is visible; Phase 2 makes the exact one the one the model uses.
    const double dqdT_dilute = q_sat_surf * L_surf / (AtmMixture::R_H2O * T_SURF * T_SURF);
    const double dqdT_exact  = SaturationH2O::dqSatdT(T_SURF, P_SURF_BAR * 1000.0, M_nw_surf);
    std::printf("      dq_sat/dT  dilute %.5g /K   exact %.5g /K   ratio %.3f\n",
                dqdT_dilute, dqdT_exact, dqdT_exact / dqdT_dilute);
    check("the exact/dilute ratio is M_other/den",
          dqdT_exact / dqdT_dilute,
          M_nw_surf / (x_H2O_surf * AtmMixture::M_H2O + (1.0 - x_H2O_surf) * M_nw_surf), 1e-2);

    // The lapse rate itself. NOTE the trap this check was written wrong the first time:
    // the familiar g/(cp + L dq/dT) keeps only the denominator of the moist-static-energy
    // balance and gives 0.7 K/km here, which is seven times too small. The pressure
    // dependence of q_sat is not negligible when water is 56 % by mole.
    const double gamma_dry     = 9.81 / cp_surf;
    const double gamma_trunc   = 9.81 / (cp_surf + L_surf * dqdT_exact);
    const double gamma_moist   = SaturationH2O::moistLapse(T_SURF, P_SURF_BAR * 1000.0,
                                     M_nw_surf, cp_surf, S.R_mix, 9.81);
    std::printf("      dry %.3f K/km   saturated %.3f K/km   (truncated form would give %.3f)\n",
                gamma_dry * 1000.0, gamma_moist * 1000.0, gamma_trunc * 1000.0);
    check("dry lapse [K/km]",                         gamma_dry * 1000.0,   10.35, 1e-2);
    check("saturated lapse [K/km]",                   gamma_moist * 1000.0, 5.70, 2e-2);
    check_true("the saturated lapse is well below the dry one",
               gamma_moist < 0.75 * gamma_dry);
    // The margin over the truncated form is 1.65x here against ATHAD_COND's 6.8x, and that
    // is PHYSICS, not a loosened tolerance. g/(cp + L*dq_sat/dT) drops the pressure dependence
    // of q_sat, which matters in proportion to how much of the mixture is water: 56 mole-% in
    // ATHAD_COND, 4.6 % here. So the wrong formula looks nearly right in this fork -- which
    // makes it more dangerous, not less, and is why the check stays rather than being deleted.
    check_true("...and still above what the truncated form claims",
               gamma_moist > 1.4 * gamma_trunc);

    // Earth, as a control on the formula itself: the standard saturated adiabatic lapse
    // rate at 300 K / 1000 hPa is ~3.6-3.9 K/km. A formula that only ever runs at 60 bar
    // has nothing to be checked against; this is the one point where it does.
    const double gamma_earth = SaturationH2O::moistLapse(300.0, 1000.0, 0.02896,
                                                          1005.0, 287.0, 9.81);
    check("Earth control: saturated lapse at 300 K / 1 bar [K/km]",
          gamma_earth * 1000.0, 3.79, 3e-2);

    // And the dry limit. Above the critical point there is no saturation curve at all, so
    // moistLapse must return g/cp EXACTLY — this is the path every unsaturated layer of the
    // column takes, and a formula that only approximately reduces to the dry adiabat would
    // put a slow bias through all of it.
    check("dry limit: moistLapse above T_crit returns g/cp",
          SaturationH2O::moistLapse(700.0, 60000.0, M_nw_surf, cp_surf, S.R_mix, 9.81),
          9.81 / cp_surf, 1e-15);
    // At 200 K there IS a saturation curve (over ice) but the vapour pressure is ~1e-4 hPa
    // of a 60 bar column, so the correction is real and utterly negligible. Asserted as
    // "negligible" rather than "zero" because it is not zero.
    check("cold limit: the moist correction at 200 K is below 1e-5 relative",
          SaturationH2O::moistLapse(200.0, 60000.0, M_nw_surf, cp_surf, S.R_mix, 9.81),
          9.81 / cp_surf, 1e-5);

    // ------------------------------------------------------------------
    // 5. Known limits of the imported fits at this regime, asserted so they cannot be
    //    forgotten rather than silently relied on.
    // ------------------------------------------------------------------
    std::printf("\n limits of the imported fits\n");
    check_true("Shomate cp fits cover the sea surface",
               T_SURF > AtmMixture::T_FIT_MIN && T_SURF < AtmMixture::T_FIT_MAX);
    check_true("but NOT the cold trap: cp below 298 K is clamped, not extrapolated",
               AtmMixture::T_FIT_MIN > 250.0);

    std::printf("\n%s  (%d failure%s)\n\n",
                failures ? "self-test FAILED" : "self-test passed",
                failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
