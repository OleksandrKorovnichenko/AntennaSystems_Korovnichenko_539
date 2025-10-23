import numpy as np
import matplotlib.pyplot as plt
import math

# === Вхідні дані ===
lam = 3.2e-2          # довжина хвилі, м
a = 13e-2              # розмір апертури по Е, м
b = 13e-2              # розмір апертури по Н, м
k = 2 * np.pi / lam    # хвильове число

print(f"λ = {lam*100:.2f} см")
print(f"a = {a*100:.2f} см, b = {b*100:.2f} см")
print(f"k = {k:.2f} рад/м")

# === Масив кутів ===
theta = np.arange(-np.pi/2, np.pi/2, 0.0001)
theta_deg = np.degrees(theta)

# === Допоміжна функція sinc ===
def sinc(x):
    return np.where(np.abs(x) < 1e-12, 1, np.sin(x)/x)

# === Формули нормованих ДС ===
FE = np.abs(sinc(k * a / 2 * np.sin(theta)))**2
FH = np.abs(sinc(k * b / 2 * np.sin(theta)))**2

# === Нормування ===
FE /= np.max(FE)
FH /= np.max(FH)

# === Визначення ширини головної пелюстки на рівні 0.707 ===
def find_hpbw(angle_deg, pattern):
    indices = np.where(pattern >= 0.707)[0]
    if len(indices) == 0:
        return 0, 0, 0
    left, right = angle_deg[indices[0]], angle_deg[indices[-1]]
    width = right - left
    return width, left, right

E_hpbw, E_left, E_right = find_hpbw(theta_deg, FE)
H_hpbw, H_left, H_right = find_hpbw(theta_deg, FH)

print(f"\nШирина головної пелюстки в площині E ≈ {E_hpbw:.2f}° (від {E_left:.2f}° до {E_right:.2f}°)")
print(f"Ширина головної пелюстки в площині H ≈ {H_hpbw:.2f}° (від {H_left:.2f}° до {H_right:.2f}°)")

# === Знаходження максимумів та мінімумів ===
def find_extrema(values, angles):
    max_x, max_y, min_x, min_y = [], [], [], []
    for i in range(1, len(values) - 1):
        if values[i] > values[i - 1] and values[i] > values[i + 1]:
            max_x.append(angles[i])
            max_y.append(values[i])
        if values[i] < values[i - 1] and values[i] < values[i + 1]:
            min_x.append(angles[i])
            min_y.append(values[i])
    return max_x, max_y, min_x, min_y

max_x_FE, max_y_FE, min_x_FE, min_y_FE = find_extrema(FE, theta_deg)
max_x_FH, max_y_FH, min_x_FH, min_y_FH = find_extrema(FH, theta_deg)

# === Таблиця екстремумів ===
print("\nТаблиця 1 – Екстремуми ДС рупорної антени")
print("-------------------------------------------")
print("| № | θmin(H) | θmin(E) | θmax(H) | FH(θ) | θmax(E) | FE(θ) |")
for i in range(max(len(max_x_FE), len(max_x_FH))):
    n = i + 1
    tminH = f"{min_x_FH[i]:.2f}" if i < len(min_x_FH) else " - "
    tminE = f"{min_x_FE[i]:.2f}" if i < len(min_x_FE) else " - "
    tmaxH = f"{max_x_FH[i]:.2f}" if i < len(max_x_FH) else " - "
    FmaxH = f"{max_y_FH[i]:.3f}" if i < len(max_y_FH) else " - "
    tmaxE = f"{max_x_FE[i]:.2f}" if i < len(max_x_FE) else " - "
    FmaxE = f"{max_y_FE[i]:.3f}" if i < len(max_y_FE) else " - "
    print(f"| {n:<2}|{tminH:<8}|{tminE:<8}|{tmaxH:<8}|{FmaxH:<7}|{tmaxE:<8}|{FmaxE:<7}|")
print("-------------------------------------------")

# === Побудова графіка ===
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(theta_deg, FE, label="E-площина", color="tab:blue")
ax.plot(theta_deg, FH, label="H-площина", color="tab:orange")
ax.axhline(0.707, color='r', linestyle='--', linewidth=0.8, label='Рівень 0.707 (½ потужності)')

# Маркери ШГП
ax.plot([E_left, E_right], [0.707, 0.707], 'bo')
ax.plot([H_left, H_right], [0.707, 0.707], 'ro')

# Анотації
ax.annotate(f"ШГП E = {E_hpbw:.2f}°", xy=(E_right, 0.707), xytext=(E_right + 3, 0.73),
             arrowprops=dict(arrowstyle="->", color="blue"), fontsize=8, color="blue")
ax.annotate(f"ШГП H = {H_hpbw:.2f}°", xy=(H_right, 0.707), xytext=(H_right + 3, 0.63),
             arrowprops=dict(arrowstyle="->", color="red"), fontsize=8, color="red")

ax.set_title("Нормовані діаграми спрямованості рупорної антени (варіант 8)")
ax.set_xlabel("Кут θ, градуси")
ax.set_ylabel("Нормована амплітуда")
ax.set_xlim(-90, 90)
ax.set_ylim(0, 1.05)
ax.legend()
ax.grid(True, linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.savefig("ДС_рупорна_антена_вар8.png", dpi=400)
plt.show()
