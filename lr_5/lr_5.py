import numpy as np
import matplotlib.pyplot as plt


def slot_pattern(angle):
    cos_val = np.cos(angle)
    safe = np.where(np.abs(cos_val) < 1e-12, 1e-12, cos_val)
    return np.cos(0.5 * np.pi * np.sin(angle)) / safe


def array_factor(angle, slots, spacing, lambda0, lambda_wg):
    k0 = 2 * np.pi / lambda0
    phase = k0 * spacing * np.sin(angle) - 2 * np.pi * spacing / lambda_wg
    result = np.ones_like(angle, dtype=float)
    denom = np.sin(phase / 2)
    mask = np.abs(denom) > 1e-12
    result[mask] = np.sin(0.5 * slots * phase[mask]) / (slots * denom[mask])
    return result


def normalize(values):
    peak = np.max(np.abs(values))
    return np.abs(values) / peak if peak else values


def find_extrema(values):
    idx = np.arange(1, len(values) - 1)
    minima = idx[(values[idx] < values[idx - 1]) & (values[idx] < values[idx + 1])]
    maxima = idx[(values[idx] > values[idx - 1]) & (values[idx] > values[idx + 1])]
    return minima, maxima


params = dict(slots=5, spacing=4.0, slot_len=2.0, lambda0=2.7, lambda_wg=3.335)
angles = np.linspace(-90.0, 90.0, 4001)
theta = np.radians(angles)

single_slot = slot_pattern(theta)
grating = array_factor(theta, params["slots"], params["spacing"], params["lambda0"], params["lambda_wg"])
pattern = normalize(single_slot * grating)
slot_norm = normalize(single_slot)
grating_norm = normalize(grating)

mins, maxs = find_extrema(pattern)
zero_angles = mins[pattern[mins] <= 1e-3]
main = np.argmax(pattern)
side_maxs = maxs[maxs != main]

plt.figure(figsize=(12, 7))
plt.plot(angles, pattern, color="navy", linewidth=2, label="|F_H(θ)|")
plt.plot(angles, slot_norm, "--", color="darkorange", linewidth=1.5, label="|F₁ᴴ(θ)|")
plt.plot(angles, grating_norm, ":", color="seagreen", linewidth=1.5, label="|FᴴC(θ)|")

for t in angles[zero_angles]:
    plt.axvline(t, color="gray", linestyle=":", alpha=0.6)

plt.scatter(angles[mins], pattern[mins], color="crimson", s=30, zorder=5, label="Мінімуми")
plt.scatter(angles[side_maxs], pattern[side_maxs], color="gold", edgecolor="black", s=40, zorder=5, label="Бічні максимуми")

plt.xlabel("θ, градуси", fontsize=12)
plt.ylabel("Нормована ДС", fontsize=12)
plt.title("ДС у площині H, варіант 17", fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend(loc="upper right", fontsize=10)
plt.xlim(-90, 90)
plt.ylim(0, 1.1)
plt.tight_layout()
plt.show()

print("Кути нульових рівнів:", np.round(angles[zero_angles], 2))
print("Локальні мінімуми:", np.round(angles[mins], 2))
print("Бічні максимуми:", np.round(angles[side_maxs], 2))