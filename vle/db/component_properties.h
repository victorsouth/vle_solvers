#pragma once
#ifndef __COMPONENT_PROPERTIES__
#define __COMPONENT_PROPERTIES__

/// @brief Используемые единицы количества вещества (мольные, массовые)
enum class AmountType { Molar, Mass };

/// @brief Корреляция Ван-Вельцена для вязкости
struct van_velzen_viscosity_correlation
{
    /// @brief коэффициент эмпирической зависимости am[20] -1 при отсутствии в данных вместо мусора
    double B{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief коэффициент эмпирической зависимости am[21] -1 при отсутствии в данных вместо мусора
    double T0{ std::numeric_limits<double>::quiet_NaN() };
};

/// @brief Параметры чистого вещества
struct component_properties_t {
    /// @brief название (формула)
    std::wstring name;
    /// @brief название (название)
    std::wstring component_name;
    /// @brief название (формула)
    std::wstring CASno;
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
    /// @brief газокинетический диаметр am[18] -1 при отсутствии в данных вместо мусора
    double gas_kinetic_diameter{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief равновесная энергия      am[19] -1 при отсутствии в данных вместо мусора
    double equilibrium_energy{ std::numeric_limits<double>::quiet_NaN() };

    /// @brief Корреляция Ван-Вельцена для вязкости am[20,21]
    van_velzen_viscosity_correlation viscosity_correlation;

    /// @brief модель давления насыщенных паров Антуана
    antoine_model_t antoine_model;
    /// @brief термодинамические функции компонента в газовом состоянии (теплоемкость Cp, энтропия, энтальпия)
    thermodynamic_functions_t functions;
    /// @brief коэффициенты теплоемкости жидкой фазы
    /// (теплоемкость мольная, проверено по воде и википедии 27.09.2022)
    fixed_solvers::ranged_polynom_t<heat_capacity_coefficients_t> heat_capacity_liquid;

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
    double estimate_antoine_extrapolation_coeff() const;
};


/// @brief База данных компонентов, индексированная по CAS-номеру.
using components_database_t = std::unordered_map<std::wstring, component_properties_t>;


#endif
