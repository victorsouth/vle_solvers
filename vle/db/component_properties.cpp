#include "../vle_solvers.h"

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
    double dx = std::max(1.0, std::abs(x)) * eps;
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
        critical_pressure * std::exp(alpha * (1 + omega) * (1 - critical_temperature / temperature));
    return result;
}

double component_properties_t::estimate_antoine_extrapolation_coeff() const
{
    size_t point_count = 20;

    Eigen::VectorXd Y(point_count);
    Eigen::MatrixXd X(point_count, 1);

    double dT = (antoine_model.max_bound - antoine_model.min_bound) / (point_count - 1);
    for (size_t index = 0; index < point_count; ++index) {
        double T = antoine_model.min_bound + index * dT;

        double Psat = get_saturated_pressure(T);
        Y(index) = std::log(Psat / critical_pressure);

        X(index, 0) = (1 + acentric_factor) * (1 - critical_temperature / T);

    }

    Eigen::VectorXd alpha = (X.transpose() * X).inverse() * X.transpose() * Y;
    return alpha(0);
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
