#include "../vle_solvers.h"


namespace vlelib {
;

double density_ideal_gas(double pressure, double temperature, double molar_mass)
{
    return pressure * molar_mass / (M_R * temperature);
}

double density_liquid_gost(double temperature, double density_20)
{
    double gamma = 0.001828 - 0.00132 * density_20 * 0.001;

    return density_20 - gamma * (temperature - 20);
}

double density_liquid_manovyan1(double temperature, double density_20)
{
    double relative_density_20 = density_20 / 1000;

    double Tcelcium = vle_solvers::kelvin2celcium(temperature);

    double result = density_20
        - 0.58 / relative_density_20 * (Tcelcium - 20)
        - fabs(Tcelcium - 1200 * (relative_density_20 - 0.68)) / 1000 * (Tcelcium - 20);

    return result;
}

std::pair<std::vector<double>, std::vector<double>> plot_enthalpy(
    fluid_t* fluid, double pressure, double Tfrom, double Tto, double Tstep /*= 0.1*/)
{
    std::pair<std::vector<double>, std::vector<double>> result;
    auto& T = result.first;
    auto& H = result.second;
    //vector<double> H, vapor_frac, Temp;
    for (double t = Tfrom; t < Tto; t += Tstep)
    {
        const auto vle = fluid->flash(pressure, t);
        double h = vle.enthalpy.mass.mix;

        T.push_back(t);
        H.push_back(h);
        //vapor_frac.push_back(vle.flash);
    }
    return result;
}





void fill_concentration_from_fluid(fluid_t* fluid, std::vector<double>* vector_concentration, size_t components_count)
{
    if (fluid == nullptr) {
        if (components_count != 0) {
            (*vector_concentration) =
                std::vector<double>(components_count, std::numeric_limits<double>::quiet_NaN());
        }
        return;
    }

    const Eigen::VectorXd& concentration = fluid->get_molar_fraction();
    if (vector_concentration->size() != concentration.size()) {
        (*vector_concentration) =
            std::vector<double>(concentration.size(), std::numeric_limits<double>::quiet_NaN());
    }

    std::copy(&concentration(0), &concentration(0) + concentration.size(),
        vector_concentration->begin());
}

void normalize_concentration(Eigen::VectorXd& molar_fraction, double min_sum_threshold)
{
    double sum = molar_fraction.sum();
    if (std::abs(sum - 1.0) > std::numeric_limits<double>::epsilon()) {
        if (std::abs(sum) < min_sum_threshold) {
            throw std::logic_error("wrong component concentrations");
        }
        molar_fraction /= sum;
    }
}

Eigen::VectorXd get_fracs_as_VectorXd(const std::vector<double>& amounts)
{
    double total_flow = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    Eigen::VectorXd fracs = Eigen::Map<Eigen::VectorXd>(const_cast<double*>(amounts.data()), amounts.size());
    fracs /= total_flow;
    return fracs;
}

std::vector<std::tuple<double, double, double>> phase_diagram(
    const fluid_t* fluid, double Pfrom, double Pto, double Tfrom, double Tto, size_t step_count /*= 10*/)
{
    double DP = (Pto - Pfrom) / step_count;
    double DT = (Tto - Tfrom) / step_count;

    std::vector<std::tuple<double, double, double>> values;
    for (double P = Pfrom; P < Pto; P += DP)
    {
        for (double T = Tfrom; T < Tto; T += DT) {
            const auto vle = fluid->flash(P, T);
            values.emplace_back(P, T, vle.flash);
        }

    }
    return values;
}

}
