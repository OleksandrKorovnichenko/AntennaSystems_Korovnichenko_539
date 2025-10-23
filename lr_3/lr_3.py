import numpy as np
import matplotlib.pyplot as plt

# Константи для розрахунків
wavelength = 3.2
dimension_h = 13  # розмір для H-площини в см
dimension_e = 13  # розмір для E-площини в см

# Обчислення нормованої функції
def calculate_normalized_function(F_values, angle_array=None):
    if angle_array is None:
        angle_array = theta
    max_value = np.abs(F_values[angle_array >= 0][0])  # максимальне значення при theta=0
    return np.abs(F_values) / max_value

# Обчислення коефіцієнта Гюйгенса
def calculate_huygens_coefficient(angle_deg):
    angle_rad = np.deg2rad(angle_deg)
    return (1 + np.cos(angle_rad)) / 2

# Обчислення sinc функції з обробкою особливого випадку
def calculate_sinc_function(x):
    return np.sinc(x / np.pi)

# Розрахунок діаграми спрямованості для H-площини
def calculate_h_plane_pattern(angle_deg):
    angle_rad = np.deg2rad(angle_deg)
    cosine_argument = np.pi * dimension_h * np.sin(angle_rad) / wavelength
    denominator_argument = 2 * dimension_h * np.sin(angle_rad) / wavelength
    cosine_component = np.cos(cosine_argument)
    denominator = 1 - denominator_argument**2
    return np.where(np.abs(denominator) > 1e-10, cosine_component / denominator, 0.0)

# Розрахунок діаграми спрямованості для E-площини
def calculate_e_plane_pattern(angle_deg):
    angle_rad = np.deg2rad(angle_deg)
    argument = np.pi * dimension_e * np.sin(angle_rad) / wavelength
    return np.where(np.abs(argument) > 1e-10, calculate_sinc_function(argument), 1.0)

# Повна діаграма спрямованості для H-площини
def calculate_full_h_pattern(angle_deg):
    return calculate_h_plane_pattern(angle_deg) * calculate_huygens_coefficient(angle_deg)

# Повна діаграма спрямованості для E-площини
def calculate_full_e_pattern(angle_deg):
    return calculate_e_plane_pattern(angle_deg) * calculate_huygens_coefficient(angle_deg)

# Створення масиву кутів для розрахунків
theta = np.arange(0, 90.1, 0.1)

# Обчислення значень для E-площини
e_plane_base_values = calculate_e_plane_pattern(theta)
huygens_coefficients = calculate_huygens_coefficient(theta)
e_plane_full_values = calculate_full_e_pattern(theta)
e_plane_normalized = calculate_normalized_function(e_plane_full_values)

# Обчислення значень для H-площини
h_plane_base_values = calculate_h_plane_pattern(theta)
h_plane_full_values = calculate_full_h_pattern(theta)
h_plane_normalized = calculate_normalized_function(h_plane_full_values)

# Визначення характеристичних точок діаграми спрямованості
# Обчислення нульових значень для E-площини: theta_0n = arcsin(n * lambda / dimension_e)
e_plane_zeros = np.arcsin(np.arange(1, 5) * wavelength / dimension_e)
e_plane_zeros_deg = np.rad2deg(e_plane_zeros[np.sin(e_plane_zeros) <= 1])

# Обчислення нульових значень для H-площини: theta_0n = arcsin((2n+1) * lambda / (2 * dimension_h))
h_plane_zeros = np.arcsin((2 * np.arange(1, 4) + 1) * wavelength / (2 * dimension_h))
h_plane_zeros_deg = np.rad2deg(h_plane_zeros[np.sin(h_plane_zeros) <= 1])

# Обчислення максимальних значень для E-площини: theta_mn = arcsin((2n+1) * lambda / (2 * dimension_e))
e_plane_maxima = np.arcsin((2 * np.arange(1, 5) + 1) * wavelength / (2 * dimension_e))
e_plane_maxima_deg = np.rad2deg(e_plane_maxima[np.sin(e_plane_maxima) <= 1])

# Обчислення максимальних значень для H-площини: theta_mn = arcsin((n+1) * lambda / dimension_h)
h_plane_maxima = np.arcsin((np.arange(1, 4) + 1) * wavelength / dimension_h)
h_plane_maxima_deg = np.rad2deg(h_plane_maxima[np.sin(h_plane_maxima) <= 1])

# Рівні бічних пелюсток для E-площини: 2 / ((2n+1) * pi)
e_plane_sidelobe_levels = [2 / ((2 * n + 1) * np.pi) for n in range(1, len(e_plane_maxima_deg) + 1)]

# Рівні бічних пелюсток для H-площини: 1 / (1 - [2(n+1)]^2)
h_plane_sidelobe_levels = [1 / (1 - (2 * (n + 1))**2) for n in range(1, len(h_plane_maxima_deg) + 1)]

# Створення графічного представлення для E-площини
e_plane_figure, e_plane_axes = plt.subplots(figsize=(8, 6))
e_plane_axes.plot(theta, e_plane_normalized, 'k-', label='|F_E(θ)|')
e_plane_axes.plot(theta, e_plane_base_values / np.max(e_plane_base_values), 'k--', label='|F(θ)|')
e_plane_axes.plot(theta, huygens_coefficients / np.max(huygens_coefficients), 'k-.', label='F_e(θ)')
e_plane_axes.set_xlabel('θ°')
e_plane_axes.set_ylabel('|F_E(θ)|')
e_plane_axes.set_title('Розрахункова ДС пірамідального рупора в площині E')
e_plane_axes.set_xlim(0, 90)
e_plane_axes.set_ylim(0, 1.1)
e_plane_axes.grid(False)

# Додавання допоміжних горизонтальних ліній
e_plane_axes.axhline(0.5, color='k', linestyle='--', linewidth=0.5)
for level in e_plane_sidelobe_levels[:len(e_plane_maxima_deg)]:
    e_plane_axes.axhline(level, color='k', linestyle='--', linewidth=0.5)

# Позначення нульових точок
for zero_point in e_plane_zeros_deg:
    e_plane_axes.plot(zero_point, 0, 'ro')  # Червона точка для нулів
    e_plane_axes.text(zero_point, -0.05, f'θ_{int(np.round(zero_point)):02d}', ha='center', va='top', fontsize=8)

# Позначення максимальних точок
for max_point in e_plane_maxima_deg:
    e_plane_axes.plot(max_point, np.interp(max_point, theta, e_plane_normalized), 'bo')  # Синя точка для максимумів

# Створення легенди з позначеннями
legend_labels = ['|F_E(θ)|', '|F(θ)|', 'F_e(θ)']
legend_points = [f'θ_{int(np.round(z)):02d}={z:.2f}°' for z in e_plane_zeros_deg] + [f'θ_m{int(np.round(m)):02d}={m:.2f}°' for m in e_plane_maxima_deg]
e_plane_axes.legend(legend_labels + legend_points, loc='upper right')

e_plane_figure.savefig('E_plane_plot.png')

# Створення графічного представлення для H-площини
h_plane_figure, h_plane_axes = plt.subplots(figsize=(8, 6))
h_plane_axes.plot(theta, h_plane_normalized, 'k-', label='|F_H(θ)|')
h_plane_axes.plot(theta, h_plane_base_values / np.max(h_plane_base_values), 'k--', label='|F(θ)|')
h_plane_axes.plot(theta, huygens_coefficients / np.max(huygens_coefficients), 'k-.', label='F_h(θ)')
h_plane_axes.set_xlabel('θ°')
h_plane_axes.set_ylabel('|F_H(θ)|')
h_plane_axes.set_title('Розрахункова ДС пірамідального рупора в площині H')
h_plane_axes.set_xlim(0, 90)
h_plane_axes.set_ylim(0, 1.1)
h_plane_axes.grid(False)

# Додавання допоміжних горизонтальних ліній
h_plane_axes.axhline(0.5, color='k', linestyle='--', linewidth=0.5)
for level in np.abs(h_plane_sidelobe_levels[:len(h_plane_maxima_deg)]):
    h_plane_axes.axhline(level, color='k', linestyle='--', linewidth=0.5)

# Позначення нульових точок
for zero_point in h_plane_zeros_deg:
    h_plane_axes.plot(zero_point, 0, 'ro')  # Червона точка для нулів
    h_plane_axes.text(zero_point, -0.05, f'θ_{int(np.round(zero_point)):02d}', ha='center', va='top', fontsize=8)

# Позначення максимальних точок
for max_point in h_plane_maxima_deg:
    h_plane_axes.plot(max_point, np.interp(max_point, theta, h_plane_normalized), 'bo')  # Синя точка для максимумів

# Створення легенди з позначеннями
legend_labels = ['|F_H(θ)|', '|F(θ)|', 'F_h(θ)']
legend_points = [f'θ_{int(np.round(z)):02d}={z:.2f}°' for z in h_plane_zeros_deg] + [f'θ_m{int(np.round(m)):02d}={m:.2f}°' for m in h_plane_maxima_deg]
h_plane_axes.legend(legend_labels + legend_points, loc='upper right')

h_plane_figure.savefig('H_plane_plot.png')

# Створення полярних діаграм спрямованості
# Полярна діаграма для E-площини
e_polar_figure, e_polar_axes = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
theta_rad = np.deg2rad(theta)

# Обчислення значень для полярної діаграми (тільки для кутів 0-180°)
theta_polar = np.arange(0, 180.1, 0.1)
theta_polar_rad = np.deg2rad(theta_polar)

# Розрахунок значень для повного діапазону кутів
e_polar_values = calculate_full_e_pattern(theta_polar)
e_polar_normalized = calculate_normalized_function(e_polar_values, theta_polar)

# Побудова полярної діаграми для E-площини
e_polar_axes.plot(theta_polar_rad, e_polar_normalized, 'b-', linewidth=2, label='E-площина')
e_polar_axes.set_ylim(0, 1)
e_polar_axes.set_title('Полярна діаграма спрямованості пірамідального рупора\nE-площина', pad=20, fontsize=14)
e_polar_axes.grid(True)
e_polar_axes.legend(loc='upper right')

# Додавання позначень для нульових і максимальних точок
for zero_point in e_plane_zeros_deg:
    e_polar_axes.plot(np.deg2rad(zero_point), 0, 'ro', markersize=8)
    e_polar_axes.text(np.deg2rad(zero_point), 0.1, f'{zero_point:.1f}°', ha='center', va='bottom', fontsize=10)

for max_point in e_plane_maxima_deg:
    max_value = np.interp(max_point, theta_polar, e_polar_normalized)
    e_polar_axes.plot(np.deg2rad(max_point), max_value, 'bo', markersize=8)
    e_polar_axes.text(np.deg2rad(max_point), max_value + 0.1, f'{max_point:.1f}°', ha='center', va='bottom', fontsize=10)

e_polar_figure.savefig('E_plane_polar_plot.png')

# Полярна діаграма для H-площини
h_polar_figure, h_polar_axes = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))

# Розрахунок значень для H-площини
h_polar_values = calculate_full_h_pattern(theta_polar)
h_polar_normalized = calculate_normalized_function(h_polar_values, theta_polar)

# Побудова полярної діаграми для H-площини
h_polar_axes.plot(theta_polar_rad, h_polar_normalized, 'r-', linewidth=2, label='H-площина')
h_polar_axes.set_ylim(0, 1)
h_polar_axes.set_title('Полярна діаграма спрямованості пірамідального рупора\nH-площина', pad=20, fontsize=14)
h_polar_axes.grid(True)
h_polar_axes.legend(loc='upper right')

# Додавання позначень для нульових і максимальних точок
for zero_point in h_plane_zeros_deg:
    h_polar_axes.plot(np.deg2rad(zero_point), 0, 'ro', markersize=8)
    h_polar_axes.text(np.deg2rad(zero_point), 0.1, f'{zero_point:.1f}°', ha='center', va='bottom', fontsize=10)

for max_point in h_plane_maxima_deg:
    max_value = np.interp(max_point, theta_polar, h_polar_normalized)
    h_polar_axes.plot(np.deg2rad(max_point), max_value, 'bo', markersize=8)
    h_polar_axes.text(np.deg2rad(max_point), max_value + 0.1, f'{max_point:.1f}°', ha='center', va='bottom', fontsize=10)

h_polar_figure.savefig('H_plane_polar_plot.png')

# Комбінована полярна діаграма для обох площин
combined_polar_figure, combined_polar_axes = plt.subplots(figsize=(12, 12), subplot_kw=dict(projection='polar'))

# Побудова обох діаграм на одному графіку
combined_polar_axes.plot(theta_polar_rad, e_polar_normalized, 'b-', linewidth=2, label='E-площина')
combined_polar_axes.plot(theta_polar_rad, h_polar_normalized, 'r-', linewidth=2, label='H-площина')
combined_polar_axes.set_ylim(0, 1)
combined_polar_axes.set_title('Полярні діаграми спрямованості пірамідального рупора\nE-площина та H-площина', pad=20, fontsize=14)
combined_polar_axes.grid(True)
combined_polar_axes.legend(loc='upper right')

# Додавання позначень для E-площини (сині точки)
for zero_point in e_plane_zeros_deg:
    combined_polar_axes.plot(np.deg2rad(zero_point), 0, 'bo', markersize=6, alpha=0.7)

for max_point in e_plane_maxima_deg:
    max_value = np.interp(max_point, theta_polar, e_polar_normalized)
    combined_polar_axes.plot(np.deg2rad(max_point), max_value, 'bs', markersize=6, alpha=0.7)

# Додавання позначень для H-площини (червоні точки)
for zero_point in h_plane_zeros_deg:
    combined_polar_axes.plot(np.deg2rad(zero_point), 0, 'ro', markersize=6, alpha=0.7)

for max_point in h_plane_maxima_deg:
    max_value = np.interp(max_point, theta_polar, h_polar_normalized)
    combined_polar_axes.plot(np.deg2rad(max_point), max_value, 'rs', markersize=6, alpha=0.7)

combined_polar_figure.savefig('Combined_polar_plot.png')

plt.show()
