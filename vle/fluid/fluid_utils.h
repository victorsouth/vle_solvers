#ifndef FLUID_UTILS_H
#define FLUID_UTILS_H

#include "fluid_common.h"
#include "components_db.h"
namespace vlelib{
;

/// @brief Строит энтальпию как функцию температуры при фиксированном давлении
/// @param fluid Флюид
/// @param pressure Давление
/// @param Tfrom Начало температурного диапазона
/// @param Tto Конец температурного диапазона
/// @param Tstep Шаг по температуре
/// @return first - температуры, second - энтальпии
std::pair<std::vector<double>, std::vector<double>> plot_enthalpy(
    fluid_t* fluid, double pressure,
    double Tfrom, double Tto, double Tstep = 0.1);


std::vector<std::tuple<double, double, double>> phase_diagram(
    const fluid_t* fluid,
    double Pfrom, double Pto,
    double Tfrom, double Tto,
    size_t step_count = 10);

//template <typename ValueType>
//inline void print_values(const char* filename, const vector<ValueType>& values)
//{
//    std::ofstream f(filename, std::ofstream::out);
//    for (const auto& point : values) {
//        auto cb = [&](size_t index, double value) {
//            f << value << ";";
//        };
//        for_each(point, cb);
//        f << std::endl;
//    }
//}

/// @brief Рассчитывает давление насыщенных паров для заданного чистого вещества
/// в заданном диапазоне по модели Антуана
/// @param component_name
/// @param Tfrom Начало диапазона
/// @param Tto Конец диапазона
/// @param Tstep Шаг по температурному диапазону
/// @return Рассчитанные значения
inline std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> plot_antoine(
    const std::wstring& component_name, double Tfrom, double Tto, double Tstep)
{
    const auto& component = components_database.at(component_name);

    double alpha = component.estimate_antoine_extrapolation_coeff();

    std::vector<double> T, Psat, Psat_extra;

    for (double t = Tfrom; t < Tto; t += Tstep)
    {
        double psat = component.antoine_model.get_saturated_pressure(t);
        double psat_extra = component.get_saturated_pressure(t);

        T.push_back(t);
        Psat.push_back(psat);
        Psat_extra.push_back(psat_extra);
    }
    return std::make_tuple(T, Psat, Psat_extra);
}

}

#endif // FLUID_UTILS_H
