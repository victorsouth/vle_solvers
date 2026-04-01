#include "../vle_solvers.h"

thermodynamic_functions_t::thermodynamic_functions_t(
    const std::vector<fixed_solvers::function_range_t<thermodynamic_functions_coefficients_t>>& ranges)
    : fixed_solvers::ranged_function_t<thermodynamic_functions_coefficients_t>(ranges)
{
}

double thermodynamic_functions_t::get_Cp_molar(double temperature) const
{
/*
    Источник: функция _GetCpSpec_IG из Maple
    формула : A[1] + Temp * A[2] + Temp ^ 2 * A[3] + Temp ^ 3 * A[4] + A[5] * Temp ^ 4;
    формула возвращает безразмерную идеальногазовую теплоемкость индивидуального компонента(отнесенную к R)
    *** В данной функции сделано приведение к размерным величинам! ***
*/
    size_t index = get_range_index(temperature);
    double result = M_R * fixed_solvers::polyval(
        ranges[index].coefficients.heat_capacity, temperature);
    return result;
}

double thermodynamic_functions_t::get_entropy_molar(double temperature) const
{
/*
    Источник: функция _GetEntrSpec_IG из Maple
    формула : A[1]*ln(T)+T*(A[2]+T*(A[3]/2.0+T*(A[4]/3.0+A[5]/4.0*T)))+A[7];
    формула возвращает безразмерную идеальногазовую стандартную энтропию индивидуального компонента(отнесенную к R)
    *** В данной функции сделано приведение к размерным величинам! ***
*/

    size_t index = get_range_index(temperature);
    const auto& A = ranges[index].coefficients.heat_capacity;
    const double& T = temperature;

    double entropy_dimensionless = A[0] * std::log(T) + T * (A[1] + T * (A[2] / 2.0 + T * (A[3] / 3.0 + A[4] / 4.0 * T))) +
        ranges[index].coefficients.entropy;

    return entropy_dimensionless * M_R;
}

double thermodynamic_functions_t::get_enthalpy_gas_molar(double temperature) const
{
/*
    Источник: функция _GetEnthSpec_IG из Maple
    формула : A[1] + T * (A[2] / 2.0 + T * (A[3] / 3.0 + T * (A[4] / 4.0 + A[5] / 5.0 * T))) + A[6] / T;
    формула возвращает безразмерную идеальногазовую энтальпию индивидуального компонента(отнесенную к R * T)
    *** В данной функции сделано приведение к размерным величинам! ***
*/

    size_t index = get_range_index(temperature);
    const auto& A = ranges[index].coefficients.heat_capacity;
    const double& T = temperature;

    double enthalpy_dimensionless = A[0] + T * (A[1] / 2.0 + T * (A[2] / 3.0 + T * (A[3] / 4.0 + A[4] / 5.0 * T))) +
        ranges[index].coefficients.enthalpy / T;

    double enthalpy = enthalpy_dimensionless * T * M_R;

    //enthalpy *= 1000; // перевод из кДж/моль в Дж/моль

    return enthalpy;
}
