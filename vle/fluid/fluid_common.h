#pragma once

//using hydraulics::celcium2kelvin;
//using hydraulics::kelvin2celcium;

using std::pair;
using std::vector;

namespace vlelib {
;


/// @brief Выбор типа алгоритма
enum class phase_equilibrium_algorithm_t {
    RaoultDalton = 0, PengRobinson = 1
};

/// @brief Фаза флюида
enum class fluid_phase_t {
    Vapor = 0b01, Liquid = 0b10
};

/// @brief Состояние флюида
enum class flash_type_t {
    Vapor = 0b01, Liquid = 0b10, TwoPhases = 0b11
};

/// @brief Распределение величины по фазам и смеси
struct amounts_per_phase {
    /// @brief Величина для жидкой фазы
    double liquid{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Величина для паровой фазы
    double vapor{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Величина для смеси
    double mix{ std::numeric_limits<double>::quiet_NaN() };
};

/// @brief Любые значения
template <typename T>
struct values_per_two_phases_t {
    /// @brief Величина для жидкой фазы.
    T liquid;
    /// @brief Величина для паровой фазы.
    T vapor;
};

/// @brief Скаляр по фазам (жидкость / пар): например Z или молярный объём, м^3/моль.
struct amounts_per_two_phases_t {
    /// @brief Величина для жидкой фазы.
    double liquid{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Величина для паровой фазы.
    double vapor{ std::numeric_limits<double>::quiet_NaN() };
};


/// @brief Псевдокритические свойства смеси (Kay): P, T, мольный объём.
struct fluid_pseudocritical_properties_t {
    /// @brief Псевдокритическое давление.
    double pressure = std::numeric_limits<double>::quiet_NaN();
    /// @brief Псевдокритическая температура.
    double temperature = std::numeric_limits<double>::quiet_NaN();
    /// @brief Псевдокритический мольный объём.
    double molar_volume = std::numeric_limits<double>::quiet_NaN();
};


/// @brief Величины в мольном и массовом выражении
struct amounts_molar_and_mass {
    /// @brief Величина в мольном выражении
    amounts_per_phase molar;
    /// @brief Величина в массовом выражении
    amounts_per_phase mass;
};

/// @brief Термодинамические функции PR (моль / масса; liquid / vapor / mix).
/// У теплоёмкостей mix всегда NaN.
struct td_functions_t {
    /// @brief Энтальпия пара, жидкости, смеси
    amounts_molar_and_mass enthalpy;
    /// @brief Внутренняя энергия пара, жидкости
    amounts_molar_and_mass inner_energy;
    /// @brief Энтропия пара, жидкости, смеси
    amounts_molar_and_mass entropy;
    /// @brief Изобарная теплоёмкость; mix (molar и mass) всегда NaN.
    amounts_molar_and_mass heat_capacity_pressure;
    /// @brief Изохорная теплоёмкость; mix (molar и mass) всегда NaN.
    amounts_molar_and_mass heat_capacity_volume;
};

/// @brief Термические коэффициенты, скорость звука, показатель адиабаты PR (Савельев разд. 1 / 6).
struct thermical_coefficients_t {
    /// @brief Изотермический модуль всестороннего сжатия K_T, Па. Формула (6.14).
    amounts_per_phase bulk_modulus_isothermal;
    /// @brief Адиабатический модуль всестороннего сжатия K_S, Па. Формулы (1.33), (6.16).
    amounts_per_phase bulk_modulus_isentropic;
    /// @brief Изотермический коэффициент сжатия β_T, 1/Па. Формулы (1.34), (6.15).
    amounts_per_phase compressibility_isothermal;
    /// @brief Адиабатический коэффициент сжатия β_S, 1/Па. Формулы (1.35), (6.18).
    amounts_per_phase compressibility_isentropic;
    /// @brief Изобарный коэффициент термического расширения α_P, 1/K. Формулы (1.36), (6.24).
    amounts_per_phase expansion_isobaric;
    /// @brief Изохорный коэффициент давления γ_v, 1/K. Формулы (1.37), (6.12).
    amounts_per_phase pressure_coefficient_isochoric;
    /// @brief Показатель адиабаты γ, безразмерный. Формулы (1.38), (1.39).
    amounts_per_phase isentropic_exponent;
    /// @brief Скорость звука a, м/с. Формулы (1.40), (6.17).
    amounts_per_phase sonic_speed;
};

/// @brief Плотность идеального газа по Менделееву-Клапейрону
double density_ideal_gas(double pressure, double temperature, double molar_mass);


/// @brief Расчет линейной температурурной поправки по ГОСТ
/// На скорую руку нашел тут: https://lektsia.com/5x7fab.html
/// @param temperature
/// @param density_20
/// @return
double density_liquid_gost(double temperature, double density_20);

/// @brief Расчет плотности нефтепродукта по формуле А.К. Мановяна
/// @param temperature Температура (формула работает до 300 град. цельсия)
/// @param density_20 Плотность вещества при 20 градусах цельсия
/// @return Плотность с учетом температурной поправки
double density_liquid_manovyan1(double temperature, double density_20);

class fluid_t;

/// @brief Строит энтальпию как функцию температуры при фиксированном давлении
/// @param fluid Флюид
/// @param pressure Давление
/// @param Tfrom Начало температурного диапазона
/// @param Tto Конец температурного диапазона
/// @param Tstep Шаг по температуре
/// @return first - температуры, second - энтальпии
pair<vector<double>, vector<double>> plot_enthalpy(
    fluid_t* fluid, double pressure,
    double Tfrom, double Tto, double Tstep = 0.1);


/// @brief Заполняет вектор концентраций
/// @param fluid
/// @param vector_concentration
/// @param components_count
void fill_concentration_from_fluid(fluid_t* fluid, vector<double>* vector_concentration, size_t components_count);

/// @brief Нормировка концентраций
/// @param molar_fraction Вектор концентраций
/// @param min_sum_threshold Порог минимальной суммы для нормировки
void normalize_concentration(Eigen::VectorXd& molar_fraction, double min_sum_threshold = 1e-8);


Eigen::VectorXd get_fracs_as_VectorXd(const vector<double>& amounts);


vector<std::tuple<double, double, double>> phase_diagram(
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
/// @param component_formula 
/// @param Tfrom Начало диапазона
/// @param Tto Конец диапазона
/// @param Tstep Шаг по температурному диапазону
/// @return Рассчитанные значения 
inline std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> plot_antoine(
    const std::wstring& component_formula, double Tfrom, double Tto, double Tstep)
{
    const auto& component = components_database.get_component_by_formula(component_formula);

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