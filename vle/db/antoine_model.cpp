#include "../components_db.h"

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
        double Psat_mmHg = std::exp(lnPsat);
        Psat_Pa = (Psat_mmHg / 760) * 1e5;
        break;
    }
    case antoine_formula::LgMmHgCelcium: {
        double lgPsat = A - B / std::max(min_denumerator, vle_solvers::kelvin2celcium(temperature) + C);
        double Psat_mmHg = std::pow(10.0, lgPsat);
        Psat_Pa = (Psat_mmHg / 760) * 1e5;
        break;
    }
    case antoine_formula::ExtLnKPaKelvin: {
        const double& D = antoine_coefficients[3];
        const double& E = antoine_coefficients[4];
        const double& F = antoine_coefficients[5];
        double lnPsat = A + B / std::max(min_denumerator, temperature + C) +
            D * std::log(temperature) + E * std::pow(temperature, F);
        double Psat_kPa = std::exp(lnPsat);
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
        double lgPsat = A - B / std::max(min_denumerator, vle_solvers::kelvin2celcium(temperature) + C);
        double Psat_bar = std::pow(10.0, lgPsat);
        Psat_Pa = Psat_bar * 1e5;
        break;
    }
    default:
        throw std::logic_error("get_saturated_pressure::Unsupported formula");
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
        double lnPsat = std::log(Psat_mmHg);
        double temperature = B / (A - lnPsat) - C;
        return temperature;
    }
    case antoine_formula::LgBarKelvin: {
        double Psat_bar = Psat_Pa / 1e5;
        double lgPsat = std::log10(Psat_bar);
        double temperature = B / (A - lgPsat) - C;
        return temperature;
    }
    default:
        throw std::logic_error("get_temperature_for_saturated_pressure::Unsupported formula");
    }
}
