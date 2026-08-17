#!/usr/bin/env python3
"""Reproduce the numerical benchmarks in the spinor-lift kappa note.

This script does NOT claim that the Planck or reduced-Planck scale is a
physical spinor partner. It evaluates a conditional hypothesis:

    P'^mu = kappa^2 Lambda^mu_nu P^nu

For timelike states, if m -> M under the scalar-lift map, then

    kappa = sqrt(M/m).

All masses are handled as energy equivalents in eV (c = 1 convention for ratios).
"""

from decimal import Decimal, getcontext
import csv
import json
from pathlib import Path

getcontext().prec = 50

PI = Decimal("3.1415926535897932384626433832795028841971693993751")
E9 = Decimal(10) ** 9

# External/benchmark inputs
PLANCK_MASS_GEV = Decimal("1.220890e19")  # 2022 CODATA central value
PLANCK_MASS_EV = PLANCK_MASS_GEV * E9
REDUCED_PLANCK_MASS_GEV = PLANCK_MASS_GEV / (Decimal(8) * PI).sqrt()
REDUCED_PLANCK_MASS_EV = REDUCED_PLANCK_MASS_GEV * E9

DELTA_M2_ATM_EV2 = Decimal("2.513e-3")  # NuFIT 6.0 normal-ordering benchmark
M_ATM_EV = DELTA_M2_ATM_EV2.sqrt()      # assumes negligible lightest mass for scale estimate
M_BENCHMARK_EV = Decimal("0.05")
KATRIN_UPPER_EV = Decimal("0.45")       # 90% C.L. upper limit on m_beta

KAPPA_MODULAR = Decimal("2.81005508341944e14")


def kappa_from_mass_ratio(high_ev: Decimal, low_ev: Decimal) -> Decimal:
    return (high_ev / low_ev).sqrt()


def low_mass_from_kappa(high_ev: Decimal, kappa: Decimal) -> Decimal:
    return high_ev / (kappa * kappa)


def high_mass_from_kappa(low_ev: Decimal, kappa: Decimal) -> Decimal:
    return low_ev * kappa * kappa


def sci(x: Decimal, digits: int = 12) -> str:
    return f"{x:.{digits}E}"


def main() -> None:
    rows = []
    for label, low_mass in [
        ("Rounded neutrino benchmark", M_BENCHMARK_EV),
        ("Atmospheric-splitting benchmark", M_ATM_EV),
        ("KATRIN upper-limit benchmark", KATRIN_UPPER_EV),
    ]:
        rows.append({
            "case": label,
            "m_low_eV": str(low_mass),
            "kappa_Planck": str(kappa_from_mass_ratio(PLANCK_MASS_EV, low_mass)),
            "kappa_reduced_Planck": str(kappa_from_mass_ratio(REDUCED_PLANCK_MASS_EV, low_mass)),
        })

    m_pred_planck = low_mass_from_kappa(PLANCK_MASS_EV, KAPPA_MODULAR)
    m_pred_reduced = low_mass_from_kappa(REDUCED_PLANCK_MASS_EV, KAPPA_MODULAR)
    partner_for_005 = high_mass_from_kappa(M_BENCHMARK_EV, KAPPA_MODULAR) / E9

    seesaw = []
    for heavy_gev in [Decimal("1e14"), Decimal("1e15")]:
        heavy_ev = heavy_gev * E9
        seesaw.append({
            "M_R_GeV": str(heavy_gev),
            "m_light_eV": str(M_BENCHMARK_EV),
            "kappa_if_scalar_lift": str(kappa_from_mass_ratio(heavy_ev, M_BENCHMARK_EV)),
        })

    result = {
        "inputs": {
            "Planck_mass_GeV": str(PLANCK_MASS_GEV),
            "Reduced_Planck_mass_GeV": str(REDUCED_PLANCK_MASS_GEV),
            "Delta_m2_atmospheric_eV2": str(DELTA_M2_ATM_EV2),
            "m_atmospheric_scale_eV": str(M_ATM_EV),
            "KATRIN_mbeta_upper_eV_90CL": str(KATRIN_UPPER_EV),
            "kappa_modular_benchmark": str(KAPPA_MODULAR),
        },
        "mass_ratio_estimates": rows,
        "inverse_predictions_using_kappa_modular": {
            "m_from_Planck_eV": str(m_pred_planck),
            "m_from_reduced_Planck_eV": str(m_pred_reduced),
            "partner_mass_for_m0p05_GeV": str(partner_for_005),
        },
        "seesaw_comparison": seesaw,
    }

    out_dir = Path(__file__).resolve().parent
    with (out_dir / "spinor_lift_results.json").open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    with (out_dir / "spinor_lift_results.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["quantity", "value"])
        writer.writerow(["Planck mass energy equivalent (GeV)", PLANCK_MASS_GEV])
        writer.writerow(["Reduced Planck mass (GeV)", REDUCED_PLANCK_MASS_GEV])
        writer.writerow(["Atmospheric mass-squared splitting (eV^2)", DELTA_M2_ATM_EV2])
        writer.writerow(["Derived atmospheric mass scale (eV)", M_ATM_EV])
        writer.writerow(["KATRIN m_beta upper limit, 90% C.L. (eV)", KATRIN_UPPER_EV])
        writer.writerow(["kappa modular benchmark", KAPPA_MODULAR])
        writer.writerow(["kappa(P, m=0.05 eV)", kappa_from_mass_ratio(PLANCK_MASS_EV, M_BENCHMARK_EV)])
        writer.writerow(["kappa(Mbar_P, m=0.05 eV)", kappa_from_mass_ratio(REDUCED_PLANCK_MASS_EV, M_BENCHMARK_EV)])
        writer.writerow(["kappa(P, m=sqrt(Delta m^2))", kappa_from_mass_ratio(PLANCK_MASS_EV, M_ATM_EV)])
        writer.writerow(["kappa(Mbar_P, m=sqrt(Delta m^2))", kappa_from_mass_ratio(REDUCED_PLANCK_MASS_EV, M_ATM_EV)])
        writer.writerow(["KATRIN-implied lower kappa bound with Planck partner", kappa_from_mass_ratio(PLANCK_MASS_EV, KATRIN_UPPER_EV)])
        writer.writerow(["KATRIN-implied lower kappa bound with reduced-Planck partner", kappa_from_mass_ratio(REDUCED_PLANCK_MASS_EV, KATRIN_UPPER_EV)])
        writer.writerow(["Predicted m using Planck scale and modular kappa (eV)", m_pred_planck])
        writer.writerow(["Predicted m using reduced Planck scale and modular kappa (eV)", m_pred_reduced])
        writer.writerow(["Required partner scale for m=0.05 eV and modular kappa (GeV)", partner_for_005])
        for entry in seesaw:
            writer.writerow([f"Scalar-lift kappa for M_R={entry['M_R_GeV']} GeV, m=0.05 eV", entry["kappa_if_scalar_lift"]])

    print("Spinor-lift kappa numerical benchmarks")
    print(f"Planck mass:            {sci(PLANCK_MASS_GEV)} GeV")
    print(f"Reduced Planck mass:    {sci(REDUCED_PLANCK_MASS_GEV)} GeV")
    print(f"sqrt(Delta m_atm^2):    {sci(M_ATM_EV)} eV")
    print(f"kappa_P (0.05 eV):      {sci(kappa_from_mass_ratio(PLANCK_MASS_EV, M_BENCHMARK_EV))}")
    print(f"kappa_Mbar (0.05 eV):   {sci(kappa_from_mass_ratio(REDUCED_PLANCK_MASS_EV, M_BENCHMARK_EV))}")
    print(f"m_pred(P):              {sci(m_pred_planck)} eV")
    print(f"m_pred(Mbar):           {sci(m_pred_reduced)} eV")
    print(f"partner for 0.05 eV:    {sci(partner_for_005)} GeV")


if __name__ == "__main__":
    main()
