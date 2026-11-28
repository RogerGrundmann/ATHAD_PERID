/*
 * TraceGasAtm.h -- SO2 and NH3 as prognostic trace gases (2026-11-28, register row PERID-TRACE-CHEM).
 *
 * WHY. SO2 and NH3 were a FIXED well-mixed share of the background (f_bg[5], f_bg[2] of q_b = 1 - q_v - q_CO2):
 * no transport, no source, no sink. They are not inert here -- the background opacity is 639x the lumped
 * kappa_bg, SO2 is 11.7 % and NH3 3.1 % of the background by mass -- and over a liquid 85 C sea with live rain
 * both are highly soluble, so rain and ocean uptake are their dominant sinks. The model had neither.
 *
 * WHAT, three steps, each behind its own knob (all need ATM_TRACE_PROG=1; with it unset nothing here runs and
 * the output is byte-identical to the model before):
 *
 *  (a) ATM_TRACE_PROG=1  -- `so2`, `nh3` are prognostic mass fractions [kg/kg of the full mixture] on CO2's path:
 *      the same RHS (advection u.grad q with CO2's central differences, diffusion with CO2's Schmidt number
 *      sc_CO2 and its density term), the same in-place RK4 stages (clamped at 0, which CO2 is not), the same
 *      radial / theta / phi boundary passes. Initial condition = the fixed background they replace,
 *      q_x = f_bg[x] * (1 - q_v - q_CO2), so at iteration 0 every derived quantity equals the old one.
 *      R_mix and cp still use the fixed background composition: the trace mass that leaves is not taken out
 *      of R_bg / M_bg (a < 0.5 % effect on R at the shipped 1 mole-%).
 *
 *  (b) ATM_TRACE_WET (default 1 with PROG) -- in-cloud removal by precipitation, operator-split after the moist
 *      physics of each iteration. The precipitation PRODUCTION rate of a cell, R_p = max(0, S_r + S_s + S_g)
 *      [kg/kg/s] (ThreeCatIceScheme: P(i) = P(i+1) + rho S dz), removes the fraction R_p dt / q_cond of the
 *      condensate; the gas goes with it in the share dissolved in the LIQUID cloud water,
 *          f_aq = H* R T L / (1 + H* R T L),   L = q_cloud rho / rho_w  [m3 water / m3 air],
 *      with the effective Henry constant (Seinfeld & Pandis, Atmospheric Chemistry and Physics, Tables 7.2-7.4):
 *          SO2  H = 1.23 M/atm (3145 K), K1 = 1.3e-2 M (1960 K), K2 = 6.6e-8 M (1500 K)
 *               H* = H (1 + K1/[H+] + K1 K2/[H+]^2)
 *          NH3  H = 62 M/atm (4110 K), Kb = 1.7e-5 M (-450 K), Kw = 1.0e-14 M2 (-6710 K)
 *               H* = H (1 + Kb [H+] / Kw)
 *      each constant K(T) = K(298.15) exp(B (1/T - 1/298.15)), [H+] = 10^-pH, pH = ATM_TRACE_PH (default 5.0,
 *      an assumption: the cloud-water pH is not computed). Ice and graupel take up nothing. The removed gas
 *      leaves the atmosphere (it is not re-released where the rain evaporates).
 *      ocean uptake, ATM_TRACE_OCEAN (default 1 with PROG): a deposition velocity at the first air cell above the
 *      ground row, q_1 *= exp(-v_dep dt / dz_1), v_dep = ATM_TRACE_VDEP [m/s] (default 0.01; both gases are
 *      highly soluble, so the uptake is limited by the air side). The ocean is an infinite sink (no back-flux),
 *      and only where FluidPhase says the surface is liquid water (code 10); elsewhere nothing is deposited.
 *      Time: dt = m.dt * L_atm / u_0 [s] per iteration (as ATM_SFC_COUPLE), times ATM_TRACE_ACCEL (default 1).
 *      One iteration is ~0.03 s here. The OCEAN sink is invisible at ACCEL 1 (4e-4 kg/m2 SO2 in 20 it); the RAIN sink
 *      is NOT: the model's microphysics converts condensate within seconds (rates like ice / dt_snow_dim), so the
 *      wet removal inherits that speed and takes ~0.1 % of the SO2 column per iteration even at ACCEL 1 (measured
 *      11-28). ACCEL is part of what has to be decided, and so is whether that precipitation rate is physical.
 *
 *  (c) ATM_TRACE_OPAC (default 1 with PROG) -- the grey gas optical depth takes the LOCAL so2, nh3:
 *          kappa_bg q_b  ->  kappa_rest (q_b - q_SO2 - q_NH3) + kappa_SO2 q_SO2 + kappa_NH3 q_NH3,
 *      kappa_rest = (kappa_bg_eff - f_SO2 kappa_SO2 - f_NH3 kappa_NH3) / (1 - f_SO2 - f_NH3), the other four.
 *      With the initial fields this equals the old term up to rounding. Not applied with ATM_BG_LUMPED=1.
 *
 * Output: the species fields `NH3` and `SO2` (mass fractions, in the panorama, sphere and slices; ATM_SPECIES_3D) are
 * the prognostic fields when the knob is on, the fixed background share otherwise.
 *
 * Monitor: traceColumn() prints the global mass-weighted mean of each and the cumulative mass removed by rain and
 * by the ocean [kg/m2, global mean], at the diagnostic stride.
 */
#pragma once

#include "cAtmosphereModel.h"
#include "MixtureAtm.h"
#include "SaturationH2O.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace TraceGas {

    inline double env_d(const char* name, double def){
        const char* e = getenv(name);
        return (e && *e) ? atof(e) : def;
    }
    inline bool on(){                                   // (a)
        static const bool v = env_d("ATM_TRACE_PROG", 0.0) != 0.0;
        return v;
    }
    inline bool wet(){                                  // (b) rain
        static const bool v = on() && env_d("ATM_TRACE_WET", 1.0) != 0.0;
        return v;
    }
    inline bool ocean(){                                // (b) ocean
        static const bool v = on() && env_d("ATM_TRACE_OCEAN", 1.0) != 0.0;
        return v;
    }
    inline bool opac(){                                 // (c)
        static const bool v = on() && env_d("ATM_TRACE_OPAC", 1.0) != 0.0
                              && env_d("ATM_BG_LUMPED", 0.0) == 0.0;      // the lumped kappa_bg has no split
        return v;
    }
    inline double pH()    { static const double v = env_d("ATM_TRACE_PH", 5.0);   return v; }
    inline double v_dep() { static const double v = env_d("ATM_TRACE_VDEP", 0.01); return v > 0.0 ? v : 0.0; }
    inline double accel() { static const double v = env_d("ATM_TRACE_ACCEL", 1.0); return v > 0.0 ? v : 1.0; }

    constexpr int I_NH3 = 2, I_SO2 = 5;                 // positions in f_bg (MixtureAtm.h BG_NAMES)

    // K(T) = K298 exp(B (1/T - 1/298.15)), B = -dH/R [K]
    inline double vant_hoff(double K298, double B, double T){
        return K298 * std::exp(B * (1.0 / T - 1.0 / 298.15));
    }
    // Effective Henry constants [M/atm]
    inline double henrySO2(double T, double h){
        const double H  = vant_hoff(1.23,   3145.0, T);
        const double K1 = vant_hoff(1.3e-2, 1960.0, T);
        const double K2 = vant_hoff(6.6e-8, 1500.0, T);
        return H * (1.0 + K1 / h + K1 * K2 / (h * h));
    }
    inline double henryNH3(double T, double h){
        const double H  = vant_hoff(62.0,    4110.0, T);
        const double Kb = vant_hoff(1.7e-5,  -450.0, T);
        const double Kw = vant_hoff(1.0e-14, -6710.0, T);
        return H * (1.0 + Kb * h / Kw);
    }
    // Aqueous share in a cell with liquid volume fraction L
    inline double f_aq(double Hstar, double T, double L){
        constexpr double R_atm = 0.082057;              // L atm / (mol K)
        const double x = Hstar * R_atm * T * L;
        return (x > 0.0) ? x / (1.0 + x) : 0.0;
    }

    // Cumulative removal, global mean [kg/m2]
    inline double& removed_rain_so2()  { static double v = 0.0; return v; }
    inline double& removed_rain_nh3()  { static double v = 0.0; return v; }
    inline double& removed_ocean_so2() { static double v = 0.0; return v; }
    inline double& removed_ocean_nh3() { static double v = 0.0; return v; }

    // (a) initial condition: the fixed background share the fields replace.
    inline void init(cAtmosphereModel& m){
        if (!on()) return;
        const double f_s = m.m_comp.f_bg[I_SO2], f_n = m.m_comp.f_bg[I_NH3];
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++){
                    const double q_v = std::max(0.0, m.c.x[i][j][k]);
                    const double q_c = std::max(0.0, AtmMixture::q_CO2_of(m.c.x[i][j][k], m.co2.x[i][j][k]));
                    const double q_b = std::max(0.0, 1.0 - q_v - q_c);
                    m.so2.x[i][j][k]  = f_s * q_b;  m.so2n.x[i][j][k] = m.so2.x[i][j][k];
                    m.nh3.x[i][j][k]  = f_n * q_b;  m.nh3n.x[i][j][k] = m.nh3.x[i][j][k];
                }
        printf("      AGCM: ATM_TRACE_PROG = 1: SO2 / NH3 prognostic (f_bg %.4f / %.4f); wet %d, ocean %d (v_dep %.3g m/s),"
               " opacity %d, pH %.2f, ACCEL %.3g\n", f_s, f_n, (int)wet(), (int)ocean(), v_dep(), (int)opac(), pH(), accel());
    }

    // (c) opacity of the background in a cell, replacing kappa_bg * q_b.
    inline double kappaTimesBackground(const cAtmosphereModel& m, double kappa_bg, double q_b, int i, int j, int k){
        if (!opac()) return kappa_bg * q_b;
        const double f_s = m.m_comp.f_bg[I_SO2], f_n = m.m_comp.f_bg[I_NH3];
        const double f_r = 1.0 - f_s - f_n;
        if (f_r <= 0.0) return kappa_bg * q_b;
        const double kappa_rest = (kappa_bg - f_s * m.kappa_SO2 - f_n * m.kappa_NH3) / f_r;
        const double q_s = std::min(std::max(0.0, m.so2.x[i][j][k]), q_b);
        const double q_n = std::min(std::max(0.0, m.nh3.x[i][j][k]), q_b - q_s);
        return kappa_rest * (q_b - q_s - q_n) + m.kappa_SO2 * q_s + m.kappa_NH3 * q_n;
    }

    // (b) sinks, operator-split, once per iteration after the moist physics.
    // THE RK BASE STATE. The iteration runs moist block (here) -> RK (so2 = so2n + k dt) -> storeIntermediateData3D
    // (so2n = so2), so so2n is the state from BEFORE this call: a change to so2 alone is overwritten by the RK step
    // in every cell it integrates (measured 11-28, ACCEL 1, 20 it: the bookkeeping removed 25 kg/m2 between
    // iterations 10 and 20 while the mean moved 4e-8). Every removal is therefore applied to so2n as well.
    inline void sinks(cAtmosphereModel& m){
        if (!wet() && !ocean()) return;
        const double dt_s = m.dt * m.L_atm / m.u_0 * accel();       // [s]
        const double h    = std::pow(10.0, -pH());
        constexpr double rho_w = 1000.0, q_thr = 1.0e-12;
        double rr_s = 0.0, rr_n = 0.0, ro_s = 0.0, ro_n = 0.0, w_sum = 0.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:rr_s, rr_n, ro_s, ro_n, w_sum)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++){
                const double wgt = std::sin(m.the.z[j]);               // area weight of the row
                w_sum += wgt;
                const int ig = m.i_topography[j][k];
                for (int i = ig + 1; i < m.im - 1; i++){
                    const double dz  = 0.5 * (m.get_layer_height(i + 1) - m.get_layer_height(i - 1));
                    const double rho = m.r_humid.x[i][j][k];
                    if (!(rho > 0.0) || !(dz > 0.0)) continue;
                    if (wet()){
                        const double R_p  = std::max(0.0, m.S_r.x[i][j][k] + m.S_s.x[i][j][k] + m.S_g.x[i][j][k]);
                        const double q_l  = std::max(0.0, m.cloud.x[i][j][k]);
                        const double q_cd = q_l + std::max(0.0, m.ice.x[i][j][k]) + std::max(0.0, m.gr.x[i][j][k]);
                        if (R_p > 0.0 && q_l > q_thr && q_cd > q_thr){
                            const double T    = m.t.x[i][j][k] * m.t_0;
                            const double L    = q_l * rho / rho_w;
                            const double frac = std::min(1.0, R_p * dt_s / q_cd);   // condensate share precipitated
                            const double ds = m.so2.x[i][j][k] * f_aq(henrySO2(T, h), T, L) * frac;
                            const double dn = m.nh3.x[i][j][k] * f_aq(henryNH3(T, h), T, L) * frac;
                            m.so2.x[i][j][k] -= ds;  m.nh3.x[i][j][k] -= dn;
                            m.so2n.x[i][j][k] -= ds; m.nh3n.x[i][j][k] -= dn;   // the RK base state too (see below)
                            rr_s += wgt * ds * rho * dz;  rr_n += wgt * dn * rho * dz;
                        }
                    }
                    if (ocean() && i == ig + 1){
                        // only over liquid water: the FluidPhase surface test (code 10)
                        const double T_s = m.t.x[ig][j][k] * m.t_0, p_s = m.p_stat.x[ig][j][k];
                        bool liquid = false;
                        if (T_s >= 273.16 && T_s < AtmMixture::T_CRIT_H2O) {
                            const double e = SaturationH2O::saturationPressure(T_s);
                            liquid = (e != SaturationH2O::NO_SATURATION && e < p_s);
                        }
                        if (liquid){
                            const double dz1 = 0.5 * (m.get_layer_height(ig + 2) - m.get_layer_height(ig));
                            const double fac = 1.0 - std::exp(-v_dep() * dt_s / dz1);
                            const double ds = m.so2.x[i][j][k] * fac, dn = m.nh3.x[i][j][k] * fac;
                            m.so2.x[i][j][k] -= ds;  m.nh3.x[i][j][k] -= dn;
                            m.so2n.x[i][j][k] -= ds; m.nh3n.x[i][j][k] -= dn;
                            ro_s += wgt * ds * rho * dz;  ro_n += wgt * dn * rho * dz;
                        }
                    }
                }
            }
        if (w_sum > 0.0){
            removed_rain_so2()  += rr_s / w_sum;  removed_rain_nh3()  += rr_n / w_sum;
            removed_ocean_so2() += ro_s / w_sum;  removed_ocean_nh3() += ro_n / w_sum;
        }
    }

    // Monitor: global (area x mass weighted) mean mass fraction and the cumulative removal.
    inline void traceColumn(cAtmosphereModel& m, bool report){
        if (!on() || !report) return;
        double ms = 0.0, mn = 0.0, mm = 0.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:ms, mn, mm)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++){
                const double wgt = std::sin(m.the.z[j]);
                for (int i = m.i_topography[j][k] + 1; i < m.im - 1; i++){
                    const double dz = 0.5 * (m.get_layer_height(i + 1) - m.get_layer_height(i - 1));
                    const double dm = wgt * m.r_humid.x[i][j][k] * dz;
                    ms += dm * m.so2.x[i][j][k];  mn += dm * m.nh3.x[i][j][k];  mm += dm;
                }
            }
        double w_sum = 0.0;
        for (int j = 0; j < m.jm; j++) w_sum += m.km * std::sin(m.the.z[j]);
        if (mm > 0.0 && w_sum > 0.0)
            printf("      AGCM: trace gases: mean q SO2 %.6e  NH3 %.6e kg/kg;  column [kg/m2] SO2 %.6e  NH3 %.6e;"
                   "  removed so far [kg/m2] rain SO2 %.4e NH3 %.4e, ocean SO2 %.4e NH3 %.4e\n", ms / mm, mn / mm,
                   ms / w_sum, mn / w_sum, removed_rain_so2(), removed_rain_nh3(), removed_ocean_so2(), removed_ocean_nh3());
    }
}
