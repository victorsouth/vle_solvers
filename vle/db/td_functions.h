#pragma once
#ifndef __TD_FUNCTIONS__
#define __TD_FUNCTIONS__

/// @brief Аппроксимационный полином теплоемкости
typedef std::vector<double> heat_capacity_coefficients_t;

/// @brief Коэффициенты термодинамических функций (теплоемкость Cp, энтропия, энтальпия)
struct thermodynamic_functions_coefficients_t {
    /// @brief полиномиальная зависимость теплоемкости Cp от температуры
    heat_capacity_coefficients_t heat_capacity;
    /// @brief для формулы энтальпии
    double enthalpy;
    /// @brief для формулы энтропии
    double entropy;
    // далее фикс для энтропии (что за "фикс"?)
    /// @brief Должен использоваться дефолтный конструктор.
    thermodynamic_functions_coefficients_t() = default;
    /// @brief Должен использоваться дефолтный конструктор копии.
    thermodynamic_functions_coefficients_t(const thermodynamic_functions_coefficients_t&) = default;
};

struct component_properties_t;

/// @brief термодинамические функции (теплоемкость Cp, энтропия, энтальпия)
class thermodynamic_functions_t
    : public fixed_solvers::ranged_function_t<thermodynamic_functions_coefficients_t>
{
    friend component_properties_t;
public:
    /// @brief Возвращает размерную идеальногазовую теплоемкость индивидуального компонента.
    double get_Cp_molar(double temperature) const;
    /// @brief Возвращает размерную идеальногазовую энтальпию индивидуального компонента.
    double get_enthalpy_gas_molar(double temperature) const;
    /// @brief Возвращает размерную идеальногазовую стандартную энтропию индивидуального компонента.
    double get_entropy_molar(double temperature) const;
    /// @brief Должен использоваться дефолтный конструктор для инициализации. (там где нет RAII)
    thermodynamic_functions_t() = default;
    /// @brief Должен использоваться дефолтный конструктор копии.
    thermodynamic_functions_t(const std::vector<fixed_solvers::function_range_t<thermodynamic_functions_coefficients_t>>& ranges);
};

#endif
