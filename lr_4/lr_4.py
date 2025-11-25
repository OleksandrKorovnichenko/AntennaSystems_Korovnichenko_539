import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from dataclasses import dataclass


@dataclass(frozen=True)
class AntennaParams:
    l: float = 23.7
    lam: float = 4.0
    d: float = 2.3
    xi: float = 1.08
    h: float = 15.0


THETA = np.arange(0, 91, 1)
THETA_RAD = np.deg2rad(THETA)
SIN = np.sin(THETA_RAD)
COS = np.cos(THETA_RAD)
TABLE_THETA = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
TABLE_H = np.array([1.0, 0.975, 0.906, 0.806, 0.693, 0.569, 0.457, 0.389, 0.349, 0.321])
TABLE_E = np.array([1.0, 0.961, 0.859, 0.771, 0.653, 0.365, 0.237, 0.131, 0.061, 0.0])


def base_pattern(cfg: AntennaParams) -> np.ndarray:
    ratio = cfg.l / cfg.lam
    numerator = (cfg.xi - 1) * np.sin(np.pi * ratio * (cfg.xi - COS))
    denominator = np.sin(np.pi * ratio * (cfg.xi - 1)) * (cfg.xi - COS)
    return np.abs(numerator / denominator)


def subarray_pattern(cfg: AntennaParams) -> np.ndarray:
    term = np.cos((np.pi * cfg.h / cfg.lam) * SIN)
    return np.abs(term)


def lookup_tables() -> tuple[np.ndarray, np.ndarray]:
    return (
        np.interp(THETA, TABLE_THETA, TABLE_H),
        np.interp(THETA, TABLE_THETA, TABLE_E),
    )


def normalize(data: np.ndarray) -> np.ndarray:
    peak = np.max(data)
    return data / peak if peak else data


def beamwidth(theta: np.ndarray, values: np.ndarray, level: float = 0.707) -> float:
    below = np.where(values < level)[0]
    if len(below) == 0:
        return 90.0
    idx_low = below[0]
    idx_high = max(idx_low - 1, 0)
    if values[idx_low] == values[idx_high]:
        return float(theta[idx_low])
    interp = theta[idx_high] + (level - values[idx_high]) * (theta[idx_low] - theta[idx_high]) / (values[idx_low] - values[idx_high])
    return round(interp, 1)


def extrema(theta: np.ndarray, values: np.ndarray):
    minima_idx, _ = find_peaks(-values, prominence=1e-4)
    peaks_idx, _ = find_peaks(values, prominence=1e-4)
    zeros = theta[minima_idx][theta[minima_idx] > 0]
    lobes = theta[peaks_idx][theta[peaks_idx] > 0]
    return zeros, lobes, values[peaks_idx][peaks_idx >= 0]


def render_plot(pattern: np.ndarray, title: str, filename: str, cfg: AntennaParams, overlay=None):
    overlay = overlay or {}
    bw = beamwidth(THETA, pattern)
    zeros, lobes, lobes_level = extrema(THETA, pattern)
    rbp_db = 20 * np.log10(pattern[30]) if pattern[30] > 0 else -np.inf

    fig, ax = plt.subplots(figsize=(13, 7))
    ax.plot(THETA, base_pattern(cfg), '#D2691E', linewidth=2.3, label='|F_b(θ)|')
    for label, data, color in overlay.values():
        ax.plot(THETA, data, color, linewidth=2.3, label=label)
    ax.plot(THETA, pattern, '#1E90FF', linewidth=3, label='F(θ)')
    ax.axhline(0.707, color='red', linestyle='--', linewidth=1.2, alpha=0.8)
    ax.axvline(bw, color='green', linestyle=':', linewidth=2)
    ax.text(bw + 1, 0.65, f'{bw}°', color='green', fontsize=10)

    ax.scatter(zeros, np.zeros_like(zeros), color='black', s=40, zorder=5)
    for z in zeros:
        ax.text(z, 0.06, f'{z}°', color='black', ha='center', fontsize=9)

    ax.scatter(lobes, pattern[np.isin(THETA, lobes)], color='gold', s=60, zorder=5, edgecolors='black')
    for ang, lvl in zip(lobes, pattern[np.isin(THETA, lobes)]):
        ax.text(ang, lvl + 0.05, f'{ang}°', color='gold', ha='center', fontsize=9)

    summary = (
        f"θ₀ = 0°\n"
        f"θ₀.₇₀₇ = {bw}°\n"
        f"ШГП = {2 * bw}°\n"
        f"РБП ≈ 30° / {rbp_db:.1f} дБ\n"
        f"Нулі: {', '.join([f'{int(z)}°' for z in zeros]) or 'немає'}\n"
        f"Бічні пелюстки: {', '.join([f'{int(p)}°' for p in lobes]) or 'немає'}"
    )

    ax.text(45, 0.35, summary, bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.95), fontsize=8.5, ha='center', va='center', linespacing=1.3)
    ax.set_title(f"{title}\nλ = {cfg.lam} см, l = {cfg.l} см, d = {cfg.d} см, ξ = {cfg.xi}, h = {cfg.h} см", fontsize=13, pad=15)
    ax.set_xlabel('θ, °')
    ax.set_ylabel('F(θ), |F(θ)|')
    ax.grid(True, alpha=0.4)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_xlim(0, 90)
    ax.set_ylim(0, 1.1)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Збережено: {filename}")


def main():
    cfg = AntennaParams()
    fb = base_pattern(cfg)
    f1h, f1e = lookup_tables()
    fa = subarray_pattern(cfg)

    h_pattern = normalize(fb * f1h)
    render_plot(
        h_pattern,
        "ДС однострижневої антени — площина H",
        "ДС_однострижнева_площина_H_v2.png",
        cfg,
        overlay={"table": (f'F₁H(θ)', f1h, '#228B22')},
    )

    e_pattern = normalize(fb * f1e)
    render_plot(
        e_pattern,
        "ДС однострижневої антени — площина E",
        "ДС_однострижнева_площина_E_v2.png",
        cfg,
        overlay={"table": (f'F₁E(θ)', f1e, '#228B22')},
    )

    h2_pattern = normalize(fb * fa)
    render_plot(
        h2_pattern,
        "ДС двохстрижневої антени — площина H",
        "ДС_двохстрижнева_площина_H_v2.png",
        cfg,
        overlay={"array": ('|cos(πh/λ·sinθ)|', fa, '#8B4513')},
    )

    e2_pattern = normalize(fb * fa * COS)
    render_plot(
        e2_pattern,
        "ДС двохстрижневої антени — площина E",
        "ДС_двохстрижнева_площина_E_v2.png",
        cfg,
        overlay={"array": ('|cos(πh/λ·sinθ)|', fa, '#8B4513')},
    )


if __name__ == "__main__":
    main()
