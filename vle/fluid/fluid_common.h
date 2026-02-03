#pragma once
#include <vector>
#include <limits>
#include <Eigen/Dense>

namespace vlelib {
;


/// @brief Выбор типа алгоритма
enum class phase_equilibrium_algorithm_t {
    RaoultDalton = 0, PengRobinson = 1
};

/// @brief Фаза флюида
enum class fluid_phase_t {
    Vapor = 0b01, Liquid = 0b10
};

/// @brief Состояние флюида
enum class flash_type_t {
    Vapor = 0b01, Liquid = 0b10, TwoPhases = 0b11
};

/// @brief Распределение величины по фазам и смеси
struct amounts_per_phase {
    /// @brief Величина для жидкой фазы
    double liquid{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Величина для паровой фазы
    double vapor{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Величина для смеси
    double mix{ std::numeric_limits<double>::quiet_NaN() };
};

/// @brief Величины в мольном и массовом выражении
struct amounts_molar_and_mass {
    /// @brief Величина в мольном выражении
    amounts_per_phase molar;
    /// @brief Величина в массовом выражении
    amounts_per_phase mass;
};

/// @brief Плотность идеального газа по Менделееву-Клапейрону
double density_ideal_gas(double pressure, double temperature, double molar_mass);


/// @brief Расчет линейной температурурной поправки по ГОСТ
/// На скорую руку нашел тут: https://lektsia.com/5x7fab.html
/// @param temperature
/// @param density_20
/// @return
double density_liquid_gost(double temperature, double density_20);

/// @brief Расчет плотности нефтепродукта по формуле А.К. Мановяна
/// @param temperature Температура (формула работает до 300 град. цельсия)
/// @param density_20 Плотность вещества при 20 градусах цельсия
/// @return Плотность с учетом температурной поправки
double density_liquid_manovyan1(double temperature, double density_20);

class fluid_t;


/// @brief Заполняет вектор концентраций
/// @param fluid
/// @param vector_concentration
/// @param components_count
void fill_concentration_from_fluid(fluid_t* fluid, std::vector<double>* vector_concentration, size_t components_count);

/// @brief Нормировка концентраций
/// @param molar_fraction Вектор концентраций
void normalize_concentration(Eigen::VectorXd& molar_fraction);


Eigen::VectorXd get_fracs_as_VectorXd(const std::vector<double>& amounts);


}
