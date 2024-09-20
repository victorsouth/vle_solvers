#pragma once
#ifndef __COMPONENTS_DB__
#define __COMPONENTS_DB__

//#include "vlelib_cmn.h"

#include <fixed/fixed.h>
#include "common_db.h"
#include <limits>
#include <array>
#include <map>
#include <unordered_map>

namespace vle_solvers
{
;

/// @brief Используемые единицы количества вещества (мольные, массовые)
enum class AmountType { Molar, Mass };

/// @brief Коэффициенты уравнения Антуана. Нулевой коэффициент - А и т.д.
typedef std::array<double, 6> antoine_coefficients_t;

/// @brief Аппроксимационный полином теплоемкости
typedef std::vector<double> heat_capacity_coefficients_t;

/// @brief Тип формулы давления насыщенных паров 
/// (номера нужны для импорта из БД!!!)
enum class antoine_formula {
    /// @brief десятичный логарифм, бар, кельвины
    LgBarKelvin = 0,
    /// @brief десятичный логарифм, бар, цельсии
    LgBarCelcium = 2,
    /// @brief десятичный логарифм, мм. рт. ст.цельсии
    LgMmHgCelcium = 3,
    /// @brief натуральный логарифм, мм. рт. ст., кельвины
    LnMmHgKelvin = 1,
    /// @brief натуаальный логарифм, кПа, кельвины
    ExtLnKPaKelvin = 4
};

/// @brief Модель Антуана для давления насыщенных паров
struct antoine_model_t {
    /// @brief размерности уравнения Антуана
    antoine_formula formula{ antoine_formula::LnMmHgKelvin};
    /// @brief коэффициенты в уравнении Антуана
    antoine_coefficients_t antoine_coefficients{ std::numeric_limits<double>::quiet_NaN() , std::numeric_limits<double>::quiet_NaN() , std::numeric_limits<double>::quiet_NaN() };
    /// @brief Нижняя граница применимости урав. Антуана
    double min_bound{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Верхняя граница применимости урав. Антуана
    double max_bound{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief экстраполирующий коэффициент в формуле Антуана
    double extrapolation_coefficient{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief расчет давления насыщенных паров
    double get_saturated_pressure(double temperature) const;
    /// @brief расчет температуры, при которой давление насыщенных паров равно заданному
    double get_temperature_for_saturated_pressure(double saturated_pressure) const;
};

/// @brief Коэффициенты термодинамических функций (теплоемкость Cp, энтропия, энтальпия)
struct thermodynamic_functions_coefficients_t {
    heat_capacity_coefficients_t heat_capacity; /// @brief полиномиальная зависимость теплоемкости Cp от температуры
    double enthalpy; /// @brief для формулы энтальпии
    double entropy; /// @brief для формулы энтропии
};

struct component_properties_t;

/// @brief термодинамические функции (теплоемкость Cp, энтропия, энтальпия)
class thermodynamic_functions_t : public ranged_function_t<thermodynamic_functions_coefficients_t>
{
    friend component_properties_t;
public:
    double get_Cp_molar(double temperature) const;
    double get_enthalpy_gas_molar(double temperature) const;
    double get_entropy_molar(double temperature) const;
    thermodynamic_functions_t() = default; // для инициализации, там где нет RAII
    thermodynamic_functions_t(const std::vector<function_range_t<thermodynamic_functions_coefficients_t>>& ranges);
};

/// @brief Корреляция Ван-Вельцена для вязкости
struct van_velzen_viscosity_correlation
{
    double B; /// @brief коэффициент эмпирической зависимости am[20]
    double T0; /// @brief коэффициент эмпирической зависимости am[21]
};

/// @brief Параметры чистого вещества
struct component_properties_t {
    /// @brief название (формула)
    std::wstring name; 
    /// @brief молярная масса
    double molar_mass; 
    /// @brief плотность жидкости при 20 град
    double density_liquid_20; 
    /// @brief коэффициент сжимаемости жидкости
    double elastic_modulus;
    /// @brief температура кипения при нормальных условиях
    double normal_boiling_temperature; 
    /// @brief критическая температура
    double critical_temperature; 
    /// @brief критическое давление
    double critical_pressure; 
    /// @brief критический молярный объем
    double critical_molarvolume; 
    /// @brief фактор ацентричности Питцера
    double acentric_factor; 
    /// @brief теплота конденсации
    double condensation_heat_molar; 
    /// @brief газокинетический диаметр am[18]
    double gas_kinetic_diameter; 
    /// @brief равновесная энергия      am[19]
    double equilibrium_energy; 

    /// @brief Корреляция Ван-Вельцена для вязкости am[20,21]
    van_velzen_viscosity_correlation viscosity_correlation;

    /// @brief модель давления насыщенных паров Антуана
    antoine_model_t antoine_model;
    /// @brief термодинамические функции компонента в газовом состоянии (теплоемкость Cp, энтропия, энтальпия)
    thermodynamic_functions_t functions;
    /// @brief коэффициенты теплоемкости жидкой фазы 
    /// (теплоемкость мольная, проверено по воде и википедии 27.09.2022)
    ranged_polynom_t<heat_capacity_coefficients_t> heat_capacity_liquid;

    /// @brief удельная энтальпия вещества в жидком состоянии
    template <AmountType amount_type>
    double get_enthalpy_liquid(double pressure, double temperature) const;

    /// @brief удельная массовая энтальпия вещества в газообразном состоянии
    template <AmountType amount_type>
    double get_enthalpy_gas(double temperature) const
    {
        if constexpr (amount_type == AmountType::Mass) {
            return functions.get_enthalpy_gas_molar(temperature) / molar_mass;
        }
        else {
            return functions.get_enthalpy_gas_molar(temperature);
        }
    }

    /// @brief удельная внутренняя энергия вещества в газовом фазовом состоянии
    template <AmountType amount_type>
    double get_inner_energy_gas(double temperature) const;

    /// @brief Удельная внутренняя энергия вещества в жидком состоянии
    template <AmountType amount_type>
    double get_inner_energy_liquid(double temperature) const;

        /// @brief удельная мольная теплоемкость вещества в газообразном состоянии
    double get_Cp_gas_molar(double temperature) const;
    /// @brief удельная массовая теплоемкость вещества в газообразном состоянии
    double get_Cp_gas_mass(double temperature) const;
    /// @brief удельная мольная энтропия вещества в газообразном состоянии
    double get_entropy_gas_molar(double temperature) const;
    /// @brief удельная массовая энтропия вещества в газообразном состоянии
    double get_entropy_gas_mass(double temperature) const;
    /// @brief расчет давления насыщенных паров по формуле экстраполяции
    double get_saturated_pressure(double temperature) const;
    /// @brief Расчет производной давления насыщенных паров, 
    /// численный расчет вызывает get_saturated_pressure
    double get_saturated_pressure_derivative(double temperature) const;

    /// @brief расчет давления насыщенных паров по формуле экстраполяции
    double get_saturated_pressure_extrapolation(double temperature) const;
    /// @brief Оценка коэффициента экстраполяции давления насыщенных паров
    double estimate_antoine_exptraploation_coeff() const;
};

typedef std::unordered_map<std::wstring, component_properties_t> components_database_t;
typedef std::map<std::wstring, component_properties_t> sorted_components_database_t;
extern const components_database_t components_database;
extern const char* thermo_db_serialized;

}

#endif

/// @brief удельная внутренняя энергия вещества в газовом фазовом состоянии

