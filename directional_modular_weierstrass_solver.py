#!/usr/bin/env python3
"""High-precision solver for the directional Modular–Weierstrass closure.

The implementation follows the equations and normalization used in the accompanying
Universe manuscript.  It supports:
  * regional and residual physical branches;
  * symmetric finite lattice sums for g2 and g3;
  * an independent Eisenstein q-series cross-check;
  * real near-square and optional complex-deformation solving;
  * complex scale coordinate K and magnitude kappa = |K|;
  * truncation scans and JSON/CSV outputs.

Only Python's standard library and mpmath are required.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import mpmath as mp

MPC_M = mp.mpf("3.0856775814913673e22")
C_M_S = mp.mpf("299792458")


class ConfigurationError(ValueError):
    """Raised when the input configuration is invalid."""


def mpf(value: Any) -> mp.mpf:
    """Convert JSON-friendly numeric input to an arbitrary-precision real."""
    if isinstance(value, bool):
        raise ConfigurationError("Boolean values are not valid numeric parameters.")
    return mp.mpf(str(value))


def mpc_to_dict(z: mp.mpc, digits: int = 30) -> Dict[str, str]:
    return {
        "real": mp.nstr(mp.re(z), digits),
        "imag": mp.nstr(mp.im(z), digits),
    }


def real_to_str(x: mp.mpf, digits: int = 30) -> str:
    return mp.nstr(x, digits)


def km_s_mpc_to_s_inverse(value: mp.mpf) -> mp.mpf:
    return value * mp.mpf("1000") / MPC_M


@dataclass(frozen=True)
class PhysicalBranch:
    name: str
    h_average_km_s_mpc: mp.mpf
    h_reference_km_s_mpc: mp.mpf
    eta: mp.mpf
    t0_s: mp.mpf

    def calculate(self) -> Dict[str, mp.mpf]:
        if self.h_average_km_s_mpc <= 0:
            raise ConfigurationError(f"{self.name}: average coefficient must be positive.")
        if self.h_reference_km_s_mpc <= 0:
            raise ConfigurationError(f"{self.name}: reference coefficient must be positive.")
        if self.eta < 1:
            raise ConfigurationError(f"{self.name}: eta must be >= 1.")
        if self.t0_s <= 0:
            raise ConfigurationError(f"{self.name}: primitive response time must be positive.")

        h_avg = km_s_mpc_to_s_inverse(self.h_average_km_s_mpc)
        h_ref = km_s_mpc_to_s_inverse(self.h_reference_km_s_mpc)
        tau_avg = 1 / h_ref
        l_eff = C_M_S / h_ref
        tau_s = self.eta * self.t0_s
        duty_factor = tau_s / tau_avg
        h_inst = h_avg / duty_factor
        beta = mp.sqrt(mp.pi * h_inst * l_eff / C_M_S)
        return {
            "h_average_s_inverse": h_avg,
            "h_reference_s_inverse": h_ref,
            "tau_average_s": tau_avg,
            "l_effective_m": l_eff,
            "tau_s": tau_s,
            "duty_factor": duty_factor,
            "h_instantaneous_s_inverse": h_inst,
            "beta": beta,
        }


def lattice_invariants(phi: mp.mpc, truncation_n: int) -> Tuple[mp.mpc, mp.mpc]:
    """Evaluate g2 and g3 using the symmetric square lattice truncation."""
    if truncation_n < 1:
        raise ConfigurationError("lattice_truncation must be >= 1.")
    g2_sum = mp.mpc(0)
    g3_sum = mp.mpc(0)
    for m in range(-truncation_n, truncation_n + 1):
        for n in range(-truncation_n, truncation_n + 1):
            if m == 0 and n == 0:
                continue
            omega = mp.mpf(m) + mp.mpf(n) * phi
            g2_sum += omega ** -4
            g3_sum += omega ** -6
    return 60 * g2_sum, 140 * g3_sum


def divisor_power_sums(max_n: int, power: int) -> List[int]:
    """Sieve for sigma_power(n), n=0..max_n."""
    sigma = [0] * (max_n + 1)
    for divisor in range(1, max_n + 1):
        contribution = divisor ** power
        for multiple in range(divisor, max_n + 1, divisor):
            sigma[multiple] += contribution
    return sigma


def q_series_invariants(phi: mp.mpc, terms: int) -> Tuple[mp.mpc, mp.mpc]:
    """Independent Eisenstein-series evaluation in the lattice normalization.

    G_4 = 2*zeta(4) E_4, G_6 = 2*zeta(6) E_6,
    g2 = 60 G_4 = (4*pi^4/3) E_4,
    g3 = 140 G_6 = (8*pi^6/27) E_6.
    """
    if terms < 1:
        raise ConfigurationError("q_series_terms must be >= 1.")
    q = mp.e ** (2 * mp.pi * 1j * phi)
    sigma3 = divisor_power_sums(terms, 3)
    sigma5 = divisor_power_sums(terms, 5)
    e4 = mp.mpc(1)
    e6 = mp.mpc(1)
    qn = q
    for n in range(1, terms + 1):
        e4 += 240 * sigma3[n] * qn
        e6 -= 504 * sigma5[n] * qn
        qn *= q
    g2 = (4 * mp.pi ** 4 / 3) * e4
    g3 = (8 * mp.pi ** 6 / 27) * e6
    return g2, g3


def modular_quantities(g2: mp.mpc, g3: mp.mpc) -> Dict[str, mp.mpc]:
    if abs(g3) == 0:
        raise ZeroDivisionError("g3 is zero at the exact square lattice; use a nonzero deformation.")
    k_complex = g2 / (mp.pi * g3)
    b_mod = (k_complex ** 2) * g2 / 3
    b_mod_alt = mp.pi * (k_complex ** 3) * g3 / 3
    discriminant = g2 ** 3 - 27 * g3 ** 2
    j_invariant = 1728 * g2 ** 3 / discriminant
    return {
        "K": k_complex,
        "kappa": abs(k_complex),
        "phase_rad": mp.arg(k_complex),
        "B_mod": b_mod,
        "B_mod_alt": b_mod_alt,
        "discriminant": discriminant,
        "j_invariant": j_invariant,
    }


def phi_from_delta(delta: mp.mpc) -> mp.mpc:
    return 1j * (1 + delta)


def solve_real_delta(
    target_beta: mp.mpf,
    truncation_n: int,
    delta_initial: mp.mpf,
) -> Dict[str, Any]:
    if target_beta <= 0 or delta_initial <= 0:
        raise ConfigurationError("target beta and initial delta must be positive.")

    def residual(log_delta: mp.mpf) -> mp.mpf:
        delta = mp.e ** log_delta
        g2, g3 = lattice_invariants(phi_from_delta(delta), truncation_n)
        return mp.re(modular_quantities(g2, g3)["B_mod"]) - target_beta

    x0 = mp.log(delta_initial)
    x1 = mp.log(delta_initial * mp.mpf("1.05"))
    log_root = mp.findroot(residual, (x0, x1), solver="secant", tol=mp.eps * 100)
    delta = mp.e ** log_root
    phi = phi_from_delta(delta)
    g2, g3 = lattice_invariants(phi, truncation_n)
    modular = modular_quantities(g2, g3)
    return {
        "delta": mp.mpc(delta, 0),
        "phi": phi,
        "g2": g2,
        "g3": g3,
        **modular,
        "target": mp.mpc(target_beta, 0),
        "closure_residual": modular["B_mod"] - target_beta,
        "alternate_form_residual": modular["B_mod"] - modular["B_mod_alt"],
    }


def solve_complex_delta(
    target_beta: mp.mpf,
    response_phase_rad: mp.mpf,
    truncation_n: int,
    delta_initial_real: mp.mpf,
    delta_initial_imag: mp.mpf,
) -> Dict[str, Any]:
    target = target_beta * mp.e ** (1j * response_phase_rad)

    def equations(delta_r: mp.mpf, delta_i: mp.mpf) -> Tuple[mp.mpf, mp.mpf]:
        delta = mp.mpc(delta_r, delta_i)
        g2, g3 = lattice_invariants(phi_from_delta(delta), truncation_n)
        residual = modular_quantities(g2, g3)["B_mod"] - target
        return mp.re(residual), mp.im(residual)

    delta_r, delta_i = mp.findroot(
        equations,
        (delta_initial_real, delta_initial_imag),
        tol=mp.eps * 100,
        maxsteps=100,
    )
    delta = mp.mpc(delta_r, delta_i)
    phi = phi_from_delta(delta)
    g2, g3 = lattice_invariants(phi, truncation_n)
    modular = modular_quantities(g2, g3)
    return {
        "delta": delta,
        "phi": phi,
        "g2": g2,
        "g3": g3,
        **modular,
        "target": target,
        "closure_residual": modular["B_mod"] - target,
        "alternate_form_residual": modular["B_mod"] - modular["B_mod_alt"],
    }


def solve_q_series_real_delta(
    target_beta: mp.mpf,
    q_terms: int,
    delta_initial: mp.mpf,
) -> Dict[str, Any]:
    def residual(log_delta: mp.mpf) -> mp.mpf:
        delta = mp.e ** log_delta
        g2, g3 = q_series_invariants(phi_from_delta(delta), q_terms)
        return mp.re(modular_quantities(g2, g3)["B_mod"]) - target_beta

    x0 = mp.log(delta_initial)
    x1 = mp.log(delta_initial * mp.mpf("1.05"))
    log_root = mp.findroot(residual, (x0, x1), solver="secant", tol=mp.eps * 100)
    delta = mp.e ** log_root
    phi = phi_from_delta(delta)
    g2, g3 = q_series_invariants(phi, q_terms)
    modular = modular_quantities(g2, g3)
    return {
        "delta": mp.mpc(delta, 0),
        "phi": phi,
        "g2": g2,
        "g3": g3,
        **modular,
        "target": mp.mpc(target_beta, 0),
        "closure_residual": modular["B_mod"] - target_beta,
        "alternate_form_residual": modular["B_mod"] - modular["B_mod_alt"],
    }


def serialize_solution(solution: Mapping[str, Any], digits: int) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in solution.items():
        if isinstance(value, mp.mpc):
            result[key] = mpc_to_dict(value, digits)
        elif isinstance(value, mp.mpf):
            result[key] = real_to_str(value, digits)
        else:
            result[key] = value
    return result


def serialize_physical(values: Mapping[str, mp.mpf], digits: int) -> Dict[str, str]:
    return {key: real_to_str(value, digits) for key, value in values.items()}


def get_required(config: Mapping[str, Any], key: str) -> Any:
    if key not in config:
        raise ConfigurationError(f"Missing required input field: {key}")
    return config[key]


def build_branch(config: Mapping[str, Any], name: str, h_average: mp.mpf, eta: mp.mpf) -> PhysicalBranch:
    return PhysicalBranch(
        name=name,
        h_average_km_s_mpc=h_average,
        h_reference_km_s_mpc=mpf(get_required(config, "H_reference_km_s_Mpc")),
        eta=eta,
        t0_s=mpf(get_required(config, "primitive_time_s")),
    )


def run(config: Mapping[str, Any]) -> Dict[str, Any]:
    precision = int(config.get("working_precision_digits", 100))
    if precision < 50:
        raise ConfigurationError("At least 50 decimal digits are required for the near-square branch.")
    mp.mp.dps = precision

    truncation_n = int(config.get("lattice_truncation", 8))
    q_terms = int(config.get("q_series_terms", 80))
    delta_initial = mpf(config.get("delta_initial_real", "1.13e-16"))
    delta_initial_imag = mpf(config.get("delta_initial_imag", "0"))
    response_phase = mpf(config.get("response_phase_rad", "0"))
    solve_mode = str(config.get("solve_mode", "real")).lower()

    h_region = mpf(get_required(config, "H_region_km_s_Mpc"))
    h_background = mpf(config.get("H_background_km_s_Mpc", "0"))
    principal_h_average = h_region - h_background
    eta_in = mpf(config.get("eta_in", "1"))
    eta_out = mpf(config.get("eta_out", "1"))

    branches: Dict[str, Any] = {}
    for direction, eta in (("radiation_to_matter", eta_in), ("matter_to_radiation", eta_out)):
        physical_branch = build_branch(config, f"principal_{direction}", principal_h_average, eta)
        physical = physical_branch.calculate()
        if solve_mode == "complex" or response_phase != 0 or delta_initial_imag != 0:
            solution = solve_complex_delta(
                physical["beta"], response_phase, truncation_n, delta_initial, delta_initial_imag
            )
        else:
            solution = solve_real_delta(physical["beta"], truncation_n, delta_initial)
        branches[direction] = {
            "physical": serialize_physical(physical, precision),
            "finite_lattice_solution": serialize_solution(solution, precision),
        }

    principal_solution_raw = solve_real_delta(
        build_branch(config, "principal_reference", principal_h_average, eta_in).calculate()["beta"],
        truncation_n,
        delta_initial,
    )
    q_solution_raw = solve_q_series_real_delta(
        build_branch(config, "principal_reference", principal_h_average, eta_in).calculate()["beta"],
        q_terms,
        mp.re(principal_solution_raw["delta"]),
    )
    relative_q_shift = (
        q_solution_raw["kappa"] - principal_solution_raw["kappa"]
    ) / principal_solution_raw["kappa"]

    scan_values: List[Dict[str, str]] = []
    for n in config.get("truncation_scan", [4, 6, 8, 10, 12, 16, 20]):
        n_int = int(n)
        solution_n = solve_real_delta(
            build_branch(config, "scan", principal_h_average, eta_in).calculate()["beta"],
            n_int,
            mp.re(principal_solution_raw["delta"]),
        )
        scan_values.append(
            {
                "N": str(n_int),
                "delta_real": real_to_str(mp.re(solution_n["delta"]), precision),
                "kappa": real_to_str(solution_n["kappa"], precision),
                "K_real": real_to_str(mp.re(solution_n["K"]), precision),
                "K_imag": real_to_str(mp.im(solution_n["K"]), precision),
            }
        )

    residual_output: Dict[str, Any] | None = None
    if bool(config.get("include_residual_branch", True)):
        h_local = mpf(get_required(config, "H_local_km_s_Mpc"))
        h_baseline = mpf(get_required(config, "H_baseline_km_s_Mpc"))
        h_residual = h_local - h_baseline
        residual_branch = build_branch(config, "residual", h_residual, eta_in)
        residual_physical = residual_branch.calculate()
        residual_solution = solve_real_delta(
            residual_physical["beta"], truncation_n, mp.re(principal_solution_raw["delta"]) * 2
        )
        residual_output = {
            "H_residual_km_s_Mpc": real_to_str(h_residual, precision),
            "physical": serialize_physical(residual_physical, precision),
            "finite_lattice_solution": serialize_solution(residual_solution, precision),
            "kappa_ratio_to_principal": real_to_str(
                residual_solution["kappa"] / principal_solution_raw["kappa"], precision
            ),
        }

    return {
        "solver": {
            "name": "directional_modular_weierstrass_solver.py",
            "working_precision_digits": precision,
            "lattice_truncation": truncation_n,
            "q_series_terms": q_terms,
            "solve_mode": solve_mode,
        },
        "input_echo": dict(config),
        "principal_branches": branches,
        "independent_q_series_check": {
            "solution": serialize_solution(q_solution_raw, precision),
            "relative_kappa_shift_vs_finite_N": real_to_str(relative_q_shift, precision),
            "relative_kappa_shift_percent": real_to_str(100 * relative_q_shift, precision),
        },
        "truncation_scan": scan_values,
        "residual_branch": residual_output,
    }


def write_csv(output: Mapping[str, Any], path: Path) -> None:
    rows: List[Tuple[str, str]] = []

    def flatten(prefix: str, value: Any) -> None:
        if isinstance(value, Mapping):
            for key, subvalue in value.items():
                flatten(f"{prefix}.{key}" if prefix else str(key), subvalue)
        elif isinstance(value, list):
            for index, subvalue in enumerate(value):
                flatten(f"{prefix}[{index}]", subvalue)
        elif value is not None:
            rows.append((prefix, str(value)))

    flatten("", output)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["field", "value"])
        writer.writerows(rows)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_json", type=Path, help="Input configuration JSON.")
    parser.add_argument("--output-json", type=Path, default=None, help="Output JSON path.")
    parser.add_argument("--output-csv", type=Path, default=None, help="Flattened output CSV path.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        with args.input_json.open("r", encoding="utf-8") as handle:
            config = json.load(handle)
        output = run(config)
        output_json = args.output_json or args.input_json.with_name(args.input_json.stem + "_output.json")
        output_csv = args.output_csv or args.input_json.with_name(args.input_json.stem + "_output.csv")
        output_json.write_text(json.dumps(output, indent=2, ensure_ascii=False), encoding="utf-8")
        write_csv(output, output_csv)
        print(f"Wrote {output_json}")
        print(f"Wrote {output_csv}")
        principal = output["principal_branches"]["radiation_to_matter"]["finite_lattice_solution"]
        print(f"delta = {principal['delta']['real']} + {principal['delta']['imag']} i")
        print(f"K = {principal['K']['real']} + {principal['K']['imag']} i")
        print(f"kappa = {principal['kappa']}")
        return 0
    except (OSError, json.JSONDecodeError, ConfigurationError, ValueError, ZeroDivisionError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
