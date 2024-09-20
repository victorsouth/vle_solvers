#include "components_db.h"

namespace vle_solvers
{
;

double antoine_model_t::get_saturated_pressure(double temperature) const
{
    const double& A = antoine_coefficients[0];
    const double& B = antoine_coefficients[1];
    const double& C = antoine_coefficients[2];

    // в градусах цельсия или кельвина, в зависимости от формулы
    double min_denumerator = 1.0; 
    // в паскалях
    double min_PsatPa = 1e-5; 

    double Psat_Pa;
    switch (formula) {
    case antoine_formula::LnMmHgKelvin: {
        double lnPsat = A - B / std::max(min_denumerator, temperature + C);
        double Psat_mmHg = exp(lnPsat);
        Psat_Pa = (Psat_mmHg / 760) * 1e5;
        break;
    }
    case antoine_formula::LgMmHgCelcium: {
        double lgPsat = A - B / std::max(min_denumerator, kelvin2celcium(temperature) + C);
        double Psat_mmHg = std::pow(10.0, lgPsat);
        Psat_Pa = (Psat_mmHg / 760) * 1e5;
        break;
    }
    case antoine_formula::ExtLnKPaKelvin: {
        const double& D = antoine_coefficients[3];
        const double& E = antoine_coefficients[4];
        const double& F = antoine_coefficients[5];
        double lnPsat = A + B / std::max(min_denumerator, temperature + C) +
            D * log(temperature) + E * pow(temperature, F);
        double Psat_kPa = exp(lnPsat);
        Psat_Pa = Psat_kPa * 1000;
        break;
    }
    case antoine_formula::LgBarKelvin: {
        double lgPsat = A - B / std::max(min_denumerator, temperature + C);
        double Psat_bar = std::pow(10.0, lgPsat);
        Psat_Pa = Psat_bar * 1e5;
        break;
    }
    case antoine_formula::LgBarCelcium: {
        double lgPsat = A - B / std::max(min_denumerator, kelvin2celcium(temperature) + C);
        double Psat_bar = std::pow(10.0, lgPsat);
        Psat_Pa = Psat_bar * 1e5;
        break;
    }
    default:
        throw std::logic_error("Unsupported formula");
    }
    return std::max(min_PsatPa, Psat_Pa);
}

double antoine_model_t::get_temperature_for_saturated_pressure(double saturated_pressure) const
{
    const double& A = antoine_coefficients[0];
    const double& B = antoine_coefficients[1];
    const double& C = antoine_coefficients[2];
    const double& Psat_Pa = saturated_pressure;

    switch (formula) {
    case antoine_formula::LnMmHgKelvin: {
        double Psat_mmHg = Psat_Pa * 760 / 1e5;
        double lnPsat = log(Psat_mmHg);
        double temperature = B / (A - lnPsat) - C;
        return temperature;
    }
    case antoine_formula::LgBarKelvin: {
        double Psat_bar = Psat_Pa / 1e5;
        double lgPsat = log10(Psat_bar);
        double temperature = B / (A - lgPsat) - C;
        return temperature;
    }
    default:
        throw std::logic_error("Unsupported formula");
    }
}

thermodynamic_functions_t::thermodynamic_functions_t(
    const std::vector<function_range_t<thermodynamic_functions_coefficients_t>>& ranges)
    : ranged_function_t(ranges)
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
    double result = M_R * polyval(ranges[index].coefficients.heat_capacity, temperature);
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

    double entropy_dimensionless = A[0] * log(T) + T * (A[1] + T * (A[2] / 2.0 + T * (A[3] / 3.0 + A[4] / 4.0 * T))) +
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

template <AmountType amount_type>
double component_properties_t::get_enthalpy_liquid(double pressure, double temperature) const
{
    if constexpr (amount_type == AmountType::Mass) {
        double U = get_inner_energy_liquid<AmountType::Mass>(temperature);
        double H = U + pressure / density_liquid_20; // стандартная плотность??

        return H;
    }
    else {
        // [Дж/кг] [кг/моль] -> [Дж/моль]
        return get_enthalpy_liquid<AmountType::Mass>(pressure, temperature) * molar_mass;
    }

    // Файл "Расчет полной энтальпии жидкости.doc" - там все плохо? (27.09.2022)
    //double Cp = heat_capacity_liquid.get_polynom_value(temperature);
    //double enthalpy =
    //    functions.get_enthalpy_gas_molar(normal_boiling_temperature) - condensation_heat_molar
    //    + Cp * (temperature - normal_boiling_temperature);

    //return enthalpy;
}

template<AmountType amount_type>
double component_properties_t::get_inner_energy_gas(double temperature) const
{
    // Описано в документе "Задачи на парожидкостное...", раздел про расчет внутренней энергии идеального газа (3.01.2023)
    // Учтены замечания от Ю.П.
    double U = -M_R * temperature + get_enthalpy_gas<AmountType::Molar>(temperature);
    if constexpr (amount_type == AmountType::Molar)
        return U;
    else
        return U / molar_mass;
}

template <typename Function>
double numerical_differentiate(Function f, double x, double eps = 1e-8) {
    double dx = std::max(1.0, abs(x)) * eps;
    double df = f(x + dx) - f(x - dx);
    return df / (2 * dx);
}

template <AmountType amount_type>
double component_properties_t::get_inner_energy_liquid(double temperature) const
{
    double Cp = heat_capacity_liquid.get_polynom_value(temperature);

    double inner_energy_molar;

    double Uliq_Tb = get_inner_energy_gas<AmountType::Molar>(normal_boiling_temperature) - condensation_heat_molar;

    if (temperature < critical_temperature) {
        inner_energy_molar = Uliq_Tb + Cp * (temperature - normal_boiling_temperature);
    }
    else
    {
        auto U = [&](double T) {
            return get_inner_energy_gas<AmountType::Molar>(T);
        };
        double dUdT = numerical_differentiate(U, critical_temperature);

        double Uliq_crit = Uliq_Tb + Cp * (critical_temperature - normal_boiling_temperature);
        inner_energy_molar = Uliq_crit + dUdT * (temperature - critical_temperature);
    }

    // возможна другая реализация - через интеграл от полинома
    // capacity.get_polynom_value_integral()
    if constexpr (amount_type == AmountType::Molar)
        return inner_energy_molar;
    else
        return inner_energy_molar / molar_mass;
}

double component_properties_t::get_Cp_gas_molar(double temperature) const
{
    return functions.get_Cp_molar(temperature);
}

double component_properties_t::get_Cp_gas_mass(double temperature) const
{
    return functions.get_Cp_molar(temperature) / molar_mass;
}

double component_properties_t::get_entropy_gas_molar(double temperature) const
{
    return functions.get_entropy_molar(temperature);
}

double component_properties_t::get_entropy_gas_mass(double temperature) const
{
    return functions.get_entropy_molar(temperature) / molar_mass;
}

double component_properties_t::get_saturated_pressure(double temperature) const
{
    if (temperature > antoine_model.max_bound) {
        double pant_max = antoine_model.get_saturated_pressure(antoine_model.max_bound);
        double pext_max = get_saturated_pressure_extrapolation(antoine_model.max_bound);
        double dp = pext_max - pant_max; // смещение - насколько экстраполяция больше Антуана

        double pext = get_saturated_pressure_extrapolation(temperature);
        double psat = pext - dp; // вычитаем смещение из экстраполяции
        return psat;
    }
    //else if (temperature < antoine_model.min_bound) {
    //    double pant_min = antoine_model.get_saturated_pressure(antoine_model.min_bound);
    //    double pext_min = get_saturated_pressure_extrapolation(antoine_model.min_bound);
    //    double dp = pext_min - pant_min; // смещение - насколько экстраполяция больше Антуана

    //    double pext = get_saturated_pressure_extrapolation(temperature);
    //    double psat = pext - dp; // вычитаем смещение из экстраполяции
    //    return psat;
    //}
    else {
        return antoine_model.get_saturated_pressure(temperature);
    }
}

double component_properties_t::get_saturated_pressure_extrapolation(double temperature) const
{
    double alpha = antoine_model.extrapolation_coefficient;
    double omega = acentric_factor;

    double result =
        critical_pressure * exp(alpha * (1 + omega) * (1 - critical_temperature / temperature));
    return result;
}

double component_properties_t::estimate_antoine_exptraploation_coeff() const
{
    throw std::runtime_error("Not impl");
    //size_t point_count = 20;

    //VectorXd Y(point_count);
    //MatrixXd X(point_count, 1);

    //double dT = (antoine_model.max_bound - antoine_model.min_bound) / (point_count - 1);
    //for (size_t index = 0; index < point_count; ++index) {
    //    double T = antoine_model.min_bound + index * dT;

    //    double Psat = get_saturated_pressure(T);
    //    Y(index) = log(Psat / critical_pressure);

    //    X(index, 0) = (1 + acentric_factor) * (1 - critical_temperature / T);

    //}

    //VectorXd alpha = (X.transpose() * X).inverse() * X.transpose() * Y;
    //return alpha(0);
}

double component_properties_t::get_saturated_pressure_derivative(double temperature) const
{
    constexpr double eps = 1e-8;
    double dT = numeric_derivative_delta(temperature, eps);

    double Psat_plus = get_saturated_pressure(temperature + dT);
    double Psat_minus = get_saturated_pressure(temperature - dT);

    return (Psat_plus - Psat_minus) / (2 * dT);
}

template double component_properties_t::get_enthalpy_liquid<AmountType::Molar>(double, double) const;
template double component_properties_t::get_enthalpy_liquid<AmountType::Mass>(double, double) const;


template double component_properties_t::get_inner_energy_gas<AmountType::Molar>(double temperature) const;
template double component_properties_t::get_inner_energy_gas<AmountType::Mass>(double temperature) const;

template double component_properties_t::get_inner_energy_liquid<AmountType::Molar>(double temperature) const;
template double component_properties_t::get_inner_energy_liquid<AmountType::Mass>(double temperature) const;

}