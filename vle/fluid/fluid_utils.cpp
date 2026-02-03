#include "fluid_utils.h"
#include "fluid_base.h"

namespace vlelib{
;

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


