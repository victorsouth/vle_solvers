#include "../vle_solvers.h"


namespace vlelib {
;


flash_type_t flash_calculation_result_t::get_flash_status() const
{
    unsigned char status =
            (fluid_vapor.get() != nullptr)
            | ((fluid_liquid.get() != nullptr) << 1);

    if (status == 0b00)
        throw std::logic_error("flash result is not calculated");

    return static_cast<flash_type_t>(status);
}

bool flash_calculation_result_t::is_liquid_only() const
{
    return fluid_vapor.get() == nullptr && fluid_liquid.get() != nullptr;
}

bool flash_calculation_result_t::is_gas_only() const
{
    return fluid_vapor.get() != nullptr && fluid_liquid.get() == nullptr;
}

bool flash_calculation_result_t::is_two_phase() const
{
    return fluid_vapor.get() != nullptr && fluid_liquid.get() != nullptr;
}

double flash_calculation_result_t::get_liquid_molar_fraction() const
{
    return 1.0 - flash;
}

double flash_calculation_result_t::get_liquid_volume_fraction() const
{
    return 1.0 - vapor_volumetric_fraction;
}

double flash_calculation_result_t::get_liquid_mass_fraction() const
{
    return 1.0 - vapor_mass_fraction;
}

bool flash_calculation_result_t::was_calculated(double _pressure, double _temperature) const
{
    if (!std::isfinite(pressure) || !std::isfinite(temperature))
        return false;
    return pressure == _pressure && temperature == _temperature;
}


void flash_calculation_result_t::invalidate_calculation()
{
    pressure = std::numeric_limits<double>::quiet_NaN();
    temperature = std::numeric_limits<double>::quiet_NaN();
    k_value.clear();
    liquid_volume_shift_mix = std::numeric_limits<double>::quiet_NaN();
    z_factor.liquid = std::numeric_limits<double>::quiet_NaN();
    z_factor.vapor = std::numeric_limits<double>::quiet_NaN();
    fluid_vapor.reset();
    fluid_liquid.reset();
}


/// @brief Возвращает последний выполненный flash-расчет
/*const flash_calculation_result_t& fluid_t::get_last_flash_result() const {
    throw logic_error("not impl");
}*/


fluid_t::fluid_t(const fluid_t& other)
    : fluid_t(other.get_components(), other.get_molar_fraction(), other.get_binary_coeffs_ref())
{ }

fluid_t::fluid_t(const std::vector<const component_properties_t*>& components)
    : fluid_fundamental_data_t(components)
    , fluid_components_functions_t(this->get_components())
    , fluid_composition_functions_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_components_functions_t&>(*this))
    , fluid_phase_criteria_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_composition_functions_t&>(*this))
    , fluid_flash_functions_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_components_functions_t&>(*this),
        static_cast<fluid_composition_functions_t&>(*this))
{

}

fluid_t::fluid_t(const std::vector<const component_properties_t*>& components, const Eigen::VectorXd& components_concentration)
    : fluid_fundamental_data_t(components, components_concentration)
    , fluid_components_functions_t(this->get_components())
    , fluid_composition_functions_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_components_functions_t&>(*this))
    , fluid_phase_criteria_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_composition_functions_t&>(*this))
    , fluid_flash_functions_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_components_functions_t&>(*this),
        static_cast<fluid_composition_functions_t&>(*this))
{

}

fluid_t::fluid_t(const std::vector<const component_properties_t*>& components, const std::vector<double>& components_concentration)
    : fluid_fundamental_data_t(components, components_concentration)
    , fluid_components_functions_t(this->get_components())
    , fluid_composition_functions_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_components_functions_t&>(*this))
    , fluid_phase_criteria_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_composition_functions_t&>(*this))
    , fluid_flash_functions_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_components_functions_t&>(*this),
        static_cast<fluid_composition_functions_t&>(*this))
{

}

fluid_t::fluid_t(const std::vector<const component_properties_t*>& components, 
    const Eigen::VectorXd& components_concentration,
    const Eigen::MatrixXd& binary_coeffs)
    : fluid_fundamental_data_t(components, components_concentration, binary_coeffs)
    , fluid_components_functions_t(this->get_components())
    , fluid_composition_functions_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_components_functions_t&>(*this))
    , fluid_phase_criteria_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_composition_functions_t&>(*this))
    , fluid_flash_functions_t(
        static_cast<fluid_fundamental_data_t&>(*this),
        static_cast<fluid_components_functions_t&>(*this),
        static_cast<fluid_composition_functions_t&>(*this))
{

}



fluid_fundamental_data_t::fluid_fundamental_data_t(const std::vector<const component_properties_t*>& components, 
                                                   const std::vector<double>& components_concentration)
    : components_(components)
    , binary_coeffs_(std::make_shared<Eigen::MatrixXd>())
{
    if (components_concentration.empty()) {
        concentration_ = Eigen::VectorXd::Ones(components.size()) / components.size();
    }
    else {
        concentration_ = Eigen::VectorXd::Map(&components_concentration[0], components_concentration.size());
        normalize_concentration(concentration_);
    }
}

fluid_fundamental_data_t::fluid_fundamental_data_t(
        const std::vector<const component_properties_t*>& components,
        const Eigen::VectorXd& components_concentration)
    : components_(components)
    , binary_coeffs_(std::make_shared<Eigen::MatrixXd>())
    , concentration_(components_concentration)
{
    if (components_concentration.size() == 0) {
        concentration_ = Eigen::VectorXd::Ones(components.size()) / components.size();
    }
    else {
        normalize_concentration(concentration_);
    }

}

fluid_fundamental_data_t::fluid_fundamental_data_t(const std::vector<const component_properties_t*>& components)
    : components_(components)
    , binary_coeffs_(std::make_shared<Eigen::MatrixXd>())

{
    concentration_ = Eigen::VectorXd::Ones(components.size()) / components.size();
}

fluid_fundamental_data_t::fluid_fundamental_data_t(const fluid_fundamental_data_t& other)
    : concentration_(other.concentration_)
    , components_(other.components_)
    , binary_coeffs_(other.binary_coeffs_)
    , last_flash_result_()
{

}

fluid_fundamental_data_t::fluid_fundamental_data_t(
    const std::vector<const component_properties_t*>& components,
    const Eigen::VectorXd& components_concentration,
    const Eigen::MatrixXd& binary_coeffs)
    : components_(components)
    , concentration_(components_concentration)
    , binary_coeffs_(std::make_shared<Eigen::MatrixXd>(binary_coeffs))
{
    if (components_concentration.size() == 0) {
        concentration_ = Eigen::VectorXd::Ones(components.size()) / components.size();
    }
    else {
        normalize_concentration(concentration_);
    }
}


void fluid_fundamental_data_t::fill_state(fluid_state_t* state) const
{
    auto molar_fraction = get_molar_fraction();
    auto last_flash_result = get_last_flash_result();

    state->flash = last_flash_result.flash;
    state->pressure = last_flash_result.pressure;
    state->temperature = last_flash_result.temperature;
    state->concentration = std::vector<double>(molar_fraction.data(), molar_fraction.data() + molar_fraction.size());
}

void fluid_fundamental_data_t::set_state(const fluid_state_t& state)
{
    const auto& components = get_components();

    if (state.concentration.size() != components.size()) {
        std::stringstream ss;
        ss << "В состоянии неверное количество концентраций. Нужно " << components.size()
           << ", пришло:" << state.concentration.size();
        std::cerr << ss.str() << std::endl;
        throw std::logic_error(ss.str());
    }
    set_molar_fraction(Eigen::VectorXd::Map(&state.concentration[0], state.concentration.size()));
    flash(state.pressure, state.temperature, state.flash);
}

size_t fluid_fundamental_data_t::get_components_count() const
{
    return components_.size();
}

const std::vector<const component_properties_t*>& fluid_fundamental_data_t::get_components() const
{
    return components_;
}

const Eigen::MatrixXd& fluid_fundamental_data_t::get_binary_coeffs_ref() const
{
    return *binary_coeffs_;
}

/// @brief Отладочный запуск мьютексов
constexpr bool debug_disable_locks = false;

const Eigen::VectorXd fluid_fundamental_data_t::get_molar_fraction() const
{
    if constexpr (debug_disable_locks == false)
    {
        std::lock_guard lk(concentration_mutex);
        return concentration_;
    }
    else {
        return concentration_;
    }
}

void fluid_fundamental_data_t::set_molar_fraction(const std::vector<double>& fraction)
{
    double* data = const_cast<double*>(fraction.data());
    Eigen::Map<Eigen::VectorXd> vector_map(data, fraction.size());
    set_molar_fraction(vector_map);
}

void fluid_fundamental_data_t::set_molar_fraction(const Eigen::VectorXd& fraction)
{
    auto func = [&]() {
        if (components_.size() != (size_t)fraction.size())
            throw std::logic_error("Wrong fraction vector size");
        concentration_ = fraction;
        normalize_concentration(concentration_);// тут видимо тоже нужно
        last_flash_result_.invalidate_calculation();
    };

    if constexpr (debug_disable_locks == false)
    {
        std::lock_guard concentration_lock(concentration_mutex);
        std::lock_guard flash_and_cache_lock(flash_and_cache_mutex);
        func();
    }
    else {
        func();
    }
}

flash_calculation_result_t fluid_fundamental_data_t::get_last_flash_result() const
{
    if constexpr (debug_disable_locks == false) {
        std::lock_guard flash_and_cache_lock(flash_and_cache_mutex);
        return last_flash_result_;
    }
    else {
        return last_flash_result_;
    }
}

const flash_calculation_result_t fluid_fundamental_data_t::flash(double pressure
    , double temperature, double initial_estimation /*= std::numeric_limits<double>::quiet_NaN()*/) const
{
    if constexpr (debug_disable_locks == false) {
        std::lock_guard lk(concentration_mutex);
        std::lock_guard flash_lock(flash_and_cache_mutex);
        flash_unsafe(pressure, temperature, initial_estimation, last_flash_result_);
        return last_flash_result_;
    }
    else {
        flash_unsafe(pressure, temperature, initial_estimation, last_flash_result_);
        return last_flash_result_;
    }

}

const fluid_stubdata_t fluid_fundamental_data_t::get_mock_data() const
{
    auto molar_conc = get_molar_fraction();

    fluid_stubdata_t result;
    for (size_t index = 0; index < components_.size(); ++index)
    {
        result.component_list.emplace_back(components_[index]->name);
        result.molar_fraction.push_back(molar_conc(index));
    }

    return result;
}

fluid_composition_functions_t::fluid_composition_functions_t(const fluid_fundamental_data_t& composition, 
                                                             const fluid_components_functions_t& components_getters)
    : composition(composition)
    , components_getters(components_getters)
{

}

Eigen::VectorXd fluid_composition_functions_t::get_saturated_pressures(double temperature) const
{
    auto molar_fraction = composition.get_molar_fraction();
    const auto& components = composition.get_components();

    Eigen::VectorXd result = Eigen::VectorXd::Zero(molar_fraction.size());
    for (size_t index = 0; index < static_cast<size_t>(molar_fraction.size()); ++index) {
        if (molar_fraction(index) != 0) {
            result(index) = components[index]->antoine_model.get_saturated_pressure(temperature);
        }
    }
    return result;
}

Eigen::VectorXd fluid_composition_functions_t::get_K_values(double pressure, double temperature) const
{
    Eigen::VectorXd result = get_saturated_pressures(temperature);
    result /= std::max(1.0, pressure);

    if (fixed_solvers::has_not_finite(result)) {
        throw std::logic_error("Infinite K-values");
    }

    return result;
}

Eigen::VectorXd fluid_composition_functions_t::get_mass_fraction() const
{
    auto molar_fraction = composition.get_molar_fraction();

    Eigen::VectorXd mass_fracs = molar_fraction.cwiseProduct(components_getters.get_molar_masses());
    mass_fracs /= mass_fracs.sum();

    return mass_fracs;
}

double fluid_composition_functions_t::get_molar_mass() const
{
    double molar_mass = components_getters.get_molar_masses().dot(composition.get_molar_fraction());
    return molar_mass;
}

double fluid_composition_functions_t::get_gas_constant() const
{
    double R = M_R / get_molar_mass();
    return R;
}

double fluid_composition_functions_t::get_pseudocritical_temperature() const
{
    size_t components_count = composition.get_components_count();
    auto molar_fraction = composition.get_molar_fraction();
    const auto& components = composition.get_components();

    double result = 0;
    for (size_t index = 0; index < components_count; ++index) {
        result += components[index]->critical_temperature * molar_fraction(index);
    }
    return result;
}

double fluid_composition_functions_t::get_pseudocritical_pressure() const
{
    size_t components_count = composition.get_components_count();
    auto molar_fraction = composition.get_molar_fraction();
    const auto& components = composition.get_components();

    double result = 0;
    for (size_t index = 0; index < components_count; ++index) {
        result += components[index]->critical_pressure * molar_fraction(index);
    }
    return result;
}

double fluid_composition_functions_t::get_pseudocritical_molar_volume() const
{
    size_t components_count = composition.get_components_count();
    auto molar_fraction = composition.get_molar_fraction();
    const auto& components = composition.get_components();

    double result = 0;
    for (size_t index = 0; index < components_count; ++index) {
        result += components[index]->critical_molarvolume * molar_fraction(index);
    }
    return result;
}

fluid_pseudocritical_properties_t fluid_composition_functions_t::get_pseudocritical_properties() const
{
    fluid_pseudocritical_properties_t result{
        get_pseudocritical_pressure(),
        get_pseudocritical_temperature(),
        get_pseudocritical_molar_volume(),
    };
    return result;
}

double fluid_composition_functions_t::get_molar_volume_vapor(double pressure, double temperature) const
{
    auto molar_fraction = composition.get_molar_fraction();
    Eigen::VectorXd density = components_getters.get_densities_vapor(pressure, temperature);
    Eigen::VectorXd M = components_getters.get_molar_masses();

    double result = M.cwiseProduct(molar_fraction).cwiseProduct(density.cwiseInverse()).sum();
    return result;
}

double fluid_composition_functions_t::get_molar_volume_liquid(double pressure, double temperature) const
{
    auto molar_fraction = composition.get_molar_fraction();
    Eigen::VectorXd density = components_getters.get_densities_liquid(pressure, temperature);
    Eigen::VectorXd M = components_getters.get_molar_masses();
    double result = M.cwiseProduct(molar_fraction).cwiseProduct(density.cwiseInverse()).sum();
    return result;
}

double fluid_composition_functions_t::get_enthalpy_td_mass_as_vapor(double /*pressure*/, double temperature) const
{
    const auto& components = composition.get_components();
    size_t components_count = composition.get_components_count();
    auto molar_fraction = composition.get_molar_fraction();

    Eigen::VectorXd result(components_count);
    for (int index = 0; index < result.size(); ++index) {
        result(index) = components[index]->get_enthalpy_gas<AmountType::Mass>(temperature);
    }
    return result.dot(molar_fraction);
}

double fluid_composition_functions_t::get_enthalpy_td_mass_as_liquid(double pressure, double temperature) const
{
    const auto& components = composition.get_components();
    size_t components_count = composition.get_components_count();
    auto molar_fraction = composition.get_molar_fraction();

    Eigen::VectorXd result(components_count);
    for (int index = 0; index < result.size(); ++index) {
        result(index) = components[index]->get_enthalpy_liquid<AmountType::Mass>(pressure, temperature);
    }
    return result.dot(molar_fraction);
}

/// @brief Возвращает удельную внутренню энергию газа на заданную единицу вещества
/// @tparam amount_type Используемые единицы количества вещества (мольные, массовые)
/// @param pressure Давление
/// @param temperature Температура
template <AmountType amount_type>
double fluid_composition_functions_t::get_inner_energy_as_vapor(double pressure, double temperature) const
{
    const auto& components = composition.get_components();
    size_t components_count = composition.get_components_count();
    auto molar_fraction = composition.get_molar_fraction();

    Eigen::VectorXd result(components_count);
    for (int index = 0; index < result.size(); ++index) {
        result(index) = components[index]->get_inner_energy_gas<amount_type>(temperature);
    }
    if constexpr (amount_type == AmountType::Mass)
            return result.dot(get_mass_fraction());
    else
    return result.dot(molar_fraction);
}


double fluid_composition_functions_t::get_max_antoine_bound() const
{
    const auto& comps = composition.get_components();

    std::vector<double> max_antoine_bound;
    std::transform(comps.begin(), comps.end(), std::back_inserter(max_antoine_bound),
                   [](const component_properties_t* c) { return c->antoine_model.max_bound; }
    );
    // учитываем вес компонентов
    double mean_max = composition.get_molar_fraction().dot(
                          Eigen::VectorXd::Map(max_antoine_bound.data(), max_antoine_bound.size()));
    // а тут не учитываем весь компонентов, просто берем минимум,
    // сколько бы ни было его в смеси
    double max = *std::min_element(max_antoine_bound.begin(), max_antoine_bound.end());
    return mean_max;
}

double fluid_composition_functions_t::get_min_antoine_bound() const
{
    const auto& comps = composition.get_components();

    std::vector<double> min_antoine_bound;
    std::transform(comps.begin(), comps.end(), std::back_inserter(min_antoine_bound),
                   [](const component_properties_t* c) { return c->antoine_model.min_bound; }
    );
    // учитываем вес компонентов
    double mean_min = composition.get_molar_fraction().dot(
                          Eigen::VectorXd::Map(min_antoine_bound.data(), min_antoine_bound.size()));
    // а тут не учитываем весь компонентов, просто берем минимум,
    // сколько бы ни было его в смеси
    double min = *std::min_element(min_antoine_bound.begin(), min_antoine_bound.end());
    return mean_min;
}

fluid_components_functions_t::fluid_components_functions_t(
        const std::vector<const component_properties_t*>& components)
    : components_(components)
{ }

Eigen::VectorXd fluid_components_functions_t::get_heat_vaporization_by_component(
        double pressure, double temperature) const
{
    size_t components_count = components_.size();
    const auto& components = components_;

    // По формуле из Chung et al. Modeling of asphaltene and wax precipation 1991
    Eigen::VectorXd result(components_count);
    for (size_t index = 0; index < components_count; ++index) {
        const double& Tb = components[index]->normal_boiling_temperature;
        const double& Tc = components[index]->critical_temperature;
        double dHtb = 1.014 * Tb * (8.75 + 4.571 * log10(Tb));
        double Ht = dHtb * pow(
                        std::max<double>(0.0, Tc - temperature) / (Tc - Tb), 0.38);
        result(index) = Ht;
    }

    return result;
}

Eigen::VectorXd fluid_components_functions_t::get_molar_masses() const
{
    size_t components_count = components_.size();
    const auto& components = components_;

    Eigen::VectorXd result(components_count);
    for (size_t index = 0; index < components_count; ++index) {
        result(index) = components[index]->molar_mass;
    }
    return result;
}

Eigen::VectorXd fluid_components_functions_t::get_Cp_molar_vapor_by_components(double pressure, double temperature) const
{
    size_t components_count = components_.size();
    const auto& components = components_;

    Eigen::VectorXd result(components_count);
    for (int index = 0; index < result.size(); ++index) {
        result(index) = components[index]->get_Cp_gas_molar(temperature);
    }
    return result;
}

Eigen::VectorXd fluid_components_functions_t::get_Cp_mass_vapor_by_components(double pressure, double temperature) const
{
    size_t components_count = components_.size();
    const auto& components = components_;

    Eigen::VectorXd result(components_count);
    for (int index = 0; index < result.size(); ++index) {
        result(index) = components[index]->get_Cp_gas_mass(temperature);
    }
    return result;
}

Eigen::VectorXd fluid_components_functions_t::get_Cp_molar_liquid_by_components(double pressure, double temperature) const
{
    size_t components_count = components_.size();
    const auto& components = components_;
    Eigen::VectorXd result(components_count);
    for (int index = 0; index < result.size(); ++index) {
        auto& capacity = components[index]->heat_capacity_liquid;
        result(index) = capacity.get_polynom_value(temperature);
    }
    return result;
}


template <AmountType amount_type>
Eigen::VectorXd fluid_components_functions_t::get_enthalpy_liquid_by_components(
        double pressure, double temperature) const
{
    size_t components_count = components_.size();
    const auto& components = components_;

    Eigen::VectorXd result(components_count);
    for (int index = 0; index < result.size(); ++index) {
        result(index) = components[index]->get_enthalpy_liquid<amount_type>(pressure, temperature);
    }
    return result;
}












bool fluid_phase_criteria_t::is_vapor_only(double pressure, double temperature) const
{
    // Допустимое отклонение расчета, при котором все равно будет приниматься однофазное газовое решение
    // сделано с целью решения проблемы ошибки машинного расчета
    double epsilon = 1e-8; // Снижено с 1e-10 при отладке 29.07.2024 BUGS-74. Если и дальше придется снижать, следует разобраться подробнее

    if (pressure <= 0)
        return true;

    double total_x_sum = vapor_only_criteria(pressure, temperature);

    // <= 1 в соответствии с расчетом этого критерия
    return total_x_sum - 1.0 <= epsilon;
}

double fluid_phase_criteria_t::vapor_only_criteria(double pressure, double temperature) const
{
    auto molar_fraction = composition.get_molar_fraction();
    // Запрещаем отрицательные давления, которые возникают при расчете производных при нулевом давлении
    pressure = std::max(0.0, pressure);
    auto invK = getters.get_K_values(pressure, temperature);
    for (size_t index = 0; index < static_cast<size_t>(molar_fraction.size()); ++index) {

        if (molar_fraction(index) != 0) {
            invK(index) = 1.0 / invK(index);
        }
        else {
            invK(index) = 0;
        }
    }
    double total_x_sum = molar_fraction.dot(invK);

    if (fixed_solvers::has_not_finite(total_x_sum)) {
        throw std::logic_error("Infinite vapor_only_criteria");
    }
    return total_x_sum;
}

bool fluid_phase_criteria_t::is_liquid_only(double pressure, double temperature) const
{
    if (pressure <= 0)
        return false;
    // Skogestad flash calculation 7.5.1
    double total_y_sum = liquid_only_criteria(pressure, temperature);
    return total_y_sum < 1.0;
}

double fluid_phase_criteria_t::liquid_only_criteria(double pressure, double temperature) const
{
    // Запрещаем отрицательные давления, которые возникают при расчете производных при нулевом давлении
    pressure = std::max(0.0, pressure);

    // Skogestad flash calculation 7.5.1
    auto molar_fraction = composition.get_molar_fraction();

    auto K = getters.get_K_values(pressure, temperature);
    double total_y_sum = K.dot(molar_fraction);
    return total_y_sum;
}

fluid_phase_criteria_t::fluid_phase_criteria_t(const fluid_fundamental_data_t& composition, 
                                               const fluid_composition_functions_t& getters)
    : composition(composition)
    , getters(getters)
{

}

fluid_flash_functions_t::fluid_flash_functions_t(const fluid_fundamental_data_t& fluid_fundamental,
                                                 const fluid_components_functions_t& components,
                                                 const fluid_composition_functions_t& composition)
    : fluid_fundamental(fluid_fundamental)
    , components(components)
    , composition(composition)
{

}

double fluid_flash_functions_t::get_heat_vaporization_mass(double pressure, double temperature) const
{
    const auto& vle = fluid_fundamental.flash(pressure, temperature);

    if (vle.is_gas_only()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    Eigen::VectorXd frac = vle.fluid_liquid->get_mass_fraction();
    Eigen::VectorXd dH = components.get_heat_vaporization_by_component(pressure, temperature);
    double result = frac.dot(dH);
    return result;
}

double fluid_flash_functions_t::get_heat_capacity_molar(double pressure, double temperature) const
{
    const flash_calculation_result_t calc = fluid_fundamental.flash(pressure, temperature);

    Eigen::VectorXd Cp_vapor_vector = components.get_Cp_molar_vapor_by_components(pressure, temperature);
    Eigen::VectorXd Cp_liquid_vector = components.get_Cp_molar_liquid_by_components(pressure, temperature);

    if (calc.fluid_vapor.get() != nullptr && calc.fluid_liquid.get() != nullptr) {
        double Cp_vapor = calc.fluid_vapor->get_molar_fraction().dot(Cp_vapor_vector);
        double Cp_liquid = calc.fluid_liquid->get_molar_fraction().dot(Cp_liquid_vector);

        double Cp = calc.flash * Cp_vapor + (1 - calc.flash) * Cp_liquid;
        return Cp;
    }
    else if (calc.fluid_vapor.get() != nullptr) {
        double Cp_vapor = calc.fluid_vapor->get_molar_fraction().dot(Cp_vapor_vector);
        return Cp_vapor;
    }
    else {
        double Cp_liquid = calc.fluid_liquid->get_molar_fraction().dot(Cp_liquid_vector);
        return Cp_liquid;
    }
}

double fluid_flash_functions_t::get_heat_capacity_mass(double pressure, double temperature) const
{
    double Cp_molar = get_heat_capacity_molar(pressure, temperature);
    double M = composition.get_molar_mass();
    return Cp_molar / M;
}

double fluid_flash_functions_t::get_heat_capacity_isochoric(double pressure, double temperature) const
{
    return get_heat_capacity_molar(pressure, temperature) - M_R;
}

double fluid_flash_functions_t::get_adiabatic_exponent(double pressure, double temperature) const
{
    return get_heat_capacity_molar(pressure, temperature) / get_heat_capacity_isochoric(pressure, temperature);
}

/// @brief Специализация get_inner_energy_as_vapor для AmountType::Mass
template double
fluid_composition_functions_t::get_inner_energy_as_vapor<AmountType::Mass>
(double pressure, double temperature) const;
/// @brief Специализация get_inner_energy_as_vapor для AmountType::Molar
template double
fluid_composition_functions_t::get_inner_energy_as_vapor<AmountType::Molar>
(double pressure, double temperature) const;

/// @brief Специализация fluid_components_functions_t::get_enthalpy_liquid_by_components для AmountType::Mass
template Eigen::VectorXd
fluid_components_functions_t::get_enthalpy_liquid_by_components<AmountType::Mass>(
double pressure, double temperature) const;
/// @brief Специализация fluid_components_functions_t::get_enthalpy_liquid_by_components для AmountType::Molar
template Eigen::VectorXd
fluid_components_functions_t::get_enthalpy_liquid_by_components<AmountType::Molar>(
double pressure, double temperature) const;

/// @brief Специализация fluid_components_functions_t::get_enthalpy_liquid_derivatives_by_components для AmountType::Molar
template std::pair<Eigen::VectorXd, Eigen::VectorXd>
fluid_components_functions_t::get_enthalpy_liquid_derivatives_by_components<AmountType::Molar>
(double pressure, double temperature) const;
/// @brief Специализация fluid_components_functions_t::get_enthalpy_liquid_derivatives_by_components для AmountType::Mass
template std::pair<Eigen::VectorXd, Eigen::VectorXd>
fluid_components_functions_t::get_enthalpy_liquid_derivatives_by_components<AmountType::Mass>
(double pressure, double temperature) const;

/// @brief Специализация fluid_components_functions_t::get_enthalpy_vapor_derivatives_by_components для AmountType::Mass 
template std::pair<Eigen::VectorXd, Eigen::VectorXd>
fluid_components_functions_t::get_enthalpy_vapor_derivatives_by_components<AmountType::Mass>(
        double pressure, double temperature) const;
/// @brief Специализация fluid_components_functions_t::get_enthalpy_vapor_derivatives_by_components для AmountType::Molar 
template std::pair<Eigen::VectorXd, Eigen::VectorXd>
fluid_components_functions_t::get_enthalpy_vapor_derivatives_by_components<AmountType::Molar>(
        double pressure, double temperature) const;


}
