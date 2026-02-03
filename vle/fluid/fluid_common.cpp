#include <cmath>
#include "physical_constants.h"
#include "helpers/physical_helpers.h"
#include "fluid_common.h"
#include "fluid_base.h"
#include <numeric>


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

void normalize_concentration(Eigen::VectorXd& molar_fraction)
{
    double sum = molar_fraction.sum();
    if (std::abs(sum - 1.0) > std::numeric_limits<double>::epsilon()) {
        if (std::abs(sum) < 1e-8) {
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

}
