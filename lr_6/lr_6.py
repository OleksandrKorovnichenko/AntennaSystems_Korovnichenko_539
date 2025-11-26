import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from scipy.special import jv
from scipy.signal import find_peaks


@dataclass
class VariantParams:
    wavelength: float
    aperture: float
    focus: float

    @property
    def radius(self):
        return self.aperture / 2

    @property
    def k(self):
        return 2 * np.pi / self.wavelength

    @property
    def p(self):
        return 2 * self.focus

    @property
    def v(self):
        return 3.5 * self.radius / self.p

    @property
    def v_ratio(self):
        return 1.5 * self.v


def safe_div(numerator, denominator, eps=1e-12):
    out = np.zeros_like(numerator)
    mask = np.abs(denominator) > eps
    out[mask] = numerator[mask] / denominator[mask]
    return out


def angle_grid(max_deg=90, samples=50001):
    theta_deg = np.linspace(0, max_deg, samples)
    return theta_deg, np.deg2rad(theta_deg)


def bessel_terms(u, v, v_ratio):
    j0_u = jv(0, u)
    j1_u = jv(1, u)
    j2_u = jv(2, u)
    j0_v = jv(0, v)
    j1_v = jv(1, v)
    j1_ratio = jv(1, v_ratio)
    j2_ratio = jv(2, v_ratio)
    term_a = safe_div(v * j1_v * j0_u - u * j1_u * j0_v, v**2 - u**2)
    term_b = np.empty_like(u)
    zero_mask = np.isclose(u, 0.0)
    term_b[zero_mask] = 0.5
    term_b[~zero_mask] = j1_u[~zero_mask] / u[~zero_mask]
    term_c = safe_div(u * j1_u * j2_ratio - v_ratio * j1_ratio * j2_u, v_ratio**2 - u**2)
    denom = 0.74 * (j1_v / v) + 0.13
    numer_h = 0.74 * term_a + 0.26 * term_b - 0.25 * term_c
    numer_e = 0.74 * term_a + 0.26 * term_b + 0.25 * term_c
    return numer_h, numer_e, denom


def build_fields(params):
    theta_deg, theta = angle_grid()
    u_vals = params.k * params.radius * np.sin(theta)
    numer_h, numer_e, denom = bessel_terms(u_vals, params.v, params.v_ratio)
    cos_sq = np.cos(theta / 2) ** 2
    fh = cos_sq * (numer_h / denom)
    fe = cos_sq * (numer_e / denom)
    return theta_deg, u_vals, fh, fe


def normalize_field(field):
    return np.abs(field) / np.max(np.abs(field))


def find_hpbw(pattern, angles_deg, level=0.707):
    idx = np.argmax(pattern)
    left = right = None
    for i in range(idx, 0, -1):
        if pattern[i] >= level and pattern[i - 1] < level:
            left = np.interp(level, [pattern[i - 1], pattern[i]], [angles_deg[i - 1], angles_deg[i]])
            break
    for i in range(idx, len(pattern) - 1):
        if pattern[i] >= level and pattern[i + 1] < level:
            right = np.interp(level, [pattern[i], pattern[i + 1]], [angles_deg[i], angles_deg[i + 1]])
            break
    if left is None or right is None:
        below = np.where(pattern < level)[0]
        if below.size == 0:
            return None, None, None
        left = angles_deg[0]
        right = angles_deg[below[0]]
    return (right - left) * 2, left, right


def select_even_peaks(indices):
    return [indices[i] for i in range(len(indices)) if i % 2 == 1]


def axes_value(angle_deg, params):
    if angle_deg is None:
        return None
    return params.k * params.radius * np.sin(np.deg2rad(angle_deg))


def plot_patterns(u_vals, fh, fe, mins_h, mins_e, peaks_h, peaks_e, h_u_h, h_u_e):
    plt.figure(figsize=(12, 6))
    plt.plot(u_vals, fh, label='F_H(θ)', color='blue')
    plt.plot(u_vals, fe, label='F_E(θ)', color='orange', linestyle='--')
    plt.axhline(0.707, color='black', linestyle=':', label='0.707')
    if h_u_h is not None:
        plt.axvline(h_u_h, color='blue', linestyle=':')
    if h_u_e is not None:
        plt.axvline(h_u_e, color='orange', linestyle=':')
    plt.scatter(u_vals[mins_h], fh[mins_h], color='red', s=25, label='Мінімум H')
    plt.scatter(u_vals[mins_e], fe[mins_e], color='red', s=25, marker='x', label='Мінімум E')
    plt.scatter(u_vals[peaks_h], fh[peaks_h], color='yellow', s=50, label='Бічні H (парні)')
    plt.scatter(u_vals[peaks_e], fe[peaks_e], color='green', s=50, label='Бічні E (парні)')
    plt.xlabel('k * R0 * sin(θ)')
    plt.ylabel('Амплітуда поля')
    plt.title('Діаграма спрямованості H та E')
    plt.grid(True)
    plt.legend()
    plt.show()


def report(label, peaks, theta_deg, u_vals, field):
    print(f'\n{label}')
    for idx in peaks:
        print(f'θ = {theta_deg[idx]:.4f}°, kR0sinθ = {u_vals[idx]:.4f}, F = {field[idx]:.4f}')


def main():
    params = VariantParams(0.027, 0.6, 0.25)
    theta_deg, u_vals, fh, fe = build_fields(params)
    fh_norm = normalize_field(fh)
    fe_norm = normalize_field(fe)
    mins_h, _ = find_peaks(-fh_norm)
    mins_e, _ = find_peaks(-fe_norm)
    peaks_h, _ = find_peaks(fh_norm)
    peaks_e, _ = find_peaks(fe_norm)
    peaks_h_sel = select_even_peaks(peaks_h)
    peaks_e_sel = select_even_peaks(peaks_e)
    hpbw_h, left_h, right_h = find_hpbw(fh_norm, theta_deg)
    hpbw_e, left_e, right_e = find_hpbw(fe_norm, theta_deg)
    left_h_u = axes_value(left_h, params)
    left_e_u = axes_value(left_e, params)
    plot_patterns(u_vals, fh, fe, mins_h, mins_e, peaks_h_sel, peaks_e_sel, left_h_u, left_e_u)
    print('Результати:')
    print(f'ШГП H-площини (подвоєна): {hpbw_h:.4f}°' if hpbw_h else 'HPBW H не знайдено')
    print(f'ШГП E-площини (подвоєна): {hpbw_e:.4f}°' if hpbw_e else 'HPBW E не знайдено')
    report('Бічні максимумі H', peaks_h, theta_deg, u_vals, fh)
    report('Бічні максимумі E', peaks_e, theta_deg, u_vals, fe)


if __name__ == '__main__':
    main()
