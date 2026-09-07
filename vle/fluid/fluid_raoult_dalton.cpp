#include "../vle_solvers.h"

using std::endl;

namespace vlelib {

/// @brief Возвращает удельную внутреннюю энергию для газа
/// @tparam amount_type Тип используемых единиц количества вещества (мольные, массовые)
/// @param c Чистое вещество (его свойства)
/// @param pressure Давление
/// @param temperature Температура
template <AmountType amount_type>
double get_energy_gas(const component_properties_t* c, double pressure, double temperature) {
    return c->get_inner_energy_gas<amount_type>(temperature);
};

/// @brief Возвращает удельную внутреннюю энергию для жидкости
/// @tparam amount_type Тип используемых единиц количества вещества (мольные, массовые)
/// @param c Чистое вещество (его свойства)
/// @param pressure Давление
/// @param temperature Температура
template <AmountType amount_type>
double get_energy_liq(const component_properties_t* c, double pressure, double temperature) {
    return c->get_inner_energy_liquid<amount_type>(temperature);
};

/// @brief Возвращает удельную энтальпию для газа
/// @tparam amount_type Тип используемых единиц количества вещества (мольные, массовые)
/// @param c Чистое вещество (его свойства)
/// @param pressure Давление 
/// @param temperature Температура
template <AmountType amount_type>
double get_enthalpy_gas(const component_properties_t* c, double pressure, double temperature) {
    return c->get_enthalpy_gas<amount_type>(temperature);
};

/// @brief Возвращает удельную энтальпию для жикости
/// @tparam amount_type Тип используемых единиц количества вещества (мольные, массовые)
/// @param c Чистое вещество (его свойства)
/// @param pressure Давление 
/// @param temperature Температура
template <AmountType amount_type>
double get_enthalpy_liq(const component_properties_t* c, double pressure, double temperature) {
    return c->get_enthalpy_liquid<amount_type>(pressure, temperature);
};


//std::atomic_int fluid_rault_dalton_t::cached{0};
//std::atomic_int fluid_rault_dalton_t::notcached{0};

/// @brief Создает флюид типа fluid_rault_dalton_t
/// @param component_names Имена компонентов
/// @param molar_fractions Мольные доли
/// @return Уникальный указатель на флюид
//template std::unique_ptr<fluid_rault_dalton_t>
//create_fluid<fluid_rault_dalton_t>(const std::vector<std::wstring>& component_names, const std::vector<double>& molar_fractions);

double fluid_rault_dalton_t::get_density_liquid(const component_properties_t& component,
                                                double pressure, double temperature)
{
    double density = density_liquid_gost(temperature, component.density_liquid_20);
    // Добавка для учета сжимаемости нефти
    density *= (1 + (pressure - ATMOSPHERIC_PRESSURE) / component.elastic_modulus);
    return density;
}

Eigen::VectorXd fluid_rault_dalton_t::get_densities_liquid(double pressure, double temperature) const
{
    const auto& components = get_components();
    size_t components_count = get_components_count();

    Eigen::VectorXd result(components_count);
    for (size_t index = 0; index < components_count; ++index) {
        // расчет плотности по Мановяну и сжимаемости
        // Считает через density_ghost, которой нет в документе [Задачи на парожидкостное равновесие]
        result(index) = get_density_liquid(*components[index], pressure, temperature);
    }
    return result;
}

Eigen::VectorXd fluid_rault_dalton_t::get_liquid_volumetric_fracs(double pressure, double temperature) const {
    auto molar_fraction = get_molar_fraction();

    Eigen::VectorXd densities = get_densities_liquid(pressure, temperature);
    Eigen::VectorXd molar_masses = get_molar_masses();

    // объем (в куб.м) доли каждого компонента для одного моля смеси
    Eigen::VectorXd frac_volumes_liquid =
            molar_fraction.cwiseProduct(molar_masses).cwiseProduct(densities.cwiseInverse());
    // сумма всех объемов - мольный объем смеси
    double molar_volume = frac_volumes_liquid.sum();

    Eigen::VectorXd volume_fraction = frac_volumes_liquid / molar_volume;
    return volume_fraction;
}

double fluid_rault_dalton_t::get_density_as_liquid(double pressure, double temperature, double* molar_volume_ptr /*= nullptr*/) const
{
    auto molar_fraction = get_molar_fraction();

    // полная копия get_liquid_volumetric_fracs до момента расчета плотности
    Eigen::VectorXd densities = get_densities_liquid(pressure, temperature);
    Eigen::VectorXd molar_masses = get_molar_masses();

    // объем (в куб.м) доли каждого компонента для одного моля смеси
    Eigen::VectorXd frac_volumes_liquid =
            molar_fraction.cwiseProduct(molar_masses).cwiseProduct(densities.cwiseInverse());
    // сумма всех объемов - мольный объем смеси
    double molar_volume = frac_volumes_liquid.sum();
    if (molar_volume_ptr != nullptr) {
        *molar_volume_ptr = molar_volume;
    }

    Eigen::VectorXd volume_fraction = frac_volumes_liquid / molar_volume;

    // расчет плотности
    // взвешиваем плотности по объемным долям
    double result = volume_fraction.dot(densities);
    return result;
}

double fluid_rault_dalton_t::get_density_as_liquid_pressure_derivative(double pressure, double temperature) const {
    const auto& components = get_components();
    size_t components_count = get_components_count();
    auto molar_fraction = get_molar_fraction();

    Eigen::VectorXd density = get_densities_liquid(pressure, temperature);

    double V = 0;
    double der_inv_V = 0; // производная по давлению инверсии объема
    for (size_t index = 0; index < components_count; ++index) {
        double beta = 1 / components[index]->elastic_modulus;
        double densityMan = density_liquid_gost(temperature,
                                                components[index]->density_liquid_20);
        double xM = molar_fraction(index) * components[index]->molar_mass;
        V += xM / density(index);
        der_inv_V += xM * beta * densityMan / fixed_solvers::sqr(density(index));
    }

    der_inv_V /= fixed_solvers::sqr(V);

    double molar_mass = get_molar_mass();
    double drho_dp = der_inv_V * molar_mass;
    return drho_dp;
}

const vlelib::flash_calculation_result_t fluid_rault_dalton_t::flash_liquid_only(double pressure, double temperature, flash_calculation_result_t& result) const
{
    result.pressure = pressure;
    result.temperature = temperature;
    result.molar_mass.mix = get_molar_mass();
    result.flash = 0;
    result.critical_temperature = get_pseudocritical_temperature();
    result.critical_pressure = get_pseudocritical_pressure();
    if (result.temperature >= result.critical_temperature &&
            result.pressure >= result.critical_pressure) {
        result.state_of_matter = state_of_matter_t::Critical;
    }
    else {
        result.state_of_matter = state_of_matter_t::Liquid;
    }
    result.vapor_volumetric_fraction = 0;
    result.vapor_mass_fraction = 0;

    result.saturated_pressure_liquid = get_bubble_point_at_given_temperature(temperature);

    result.density.mix = result.density.liquid =
            get_density_as_liquid(pressure, temperature, 
                &result.molar_volume.mix);

    // Добавка для корректности (?)
    result.molar_mass.liquid = result.molar_mass.mix;
    result.molar_mass.vapor = 0;
    result.density.vapor = 0;
    result.molar_volume.liquid = result.molar_volume.mix;
    result.molar_volume.vapor = 0;
    // Конец добавки

    // Копируем жидкость без мемоизации, только состав.
    // Это потом попадет в мемоизацию fluid_liquid
    result.fluid_liquid = create_copy(false);
    result.fluid_vapor = nullptr;
    result.z_factor.vapor = 1.0;

    //flash_enthalpy(pressure, temperature);
    result.td_functions.enthalpy.molar = flash_function<AmountType::Molar>(
                                get_enthalpy_gas<AmountType::Molar>, get_enthalpy_liq<AmountType::Molar>, result);
    result.td_functions.enthalpy.mass = flash_function<AmountType::Mass>(
                               get_enthalpy_gas<AmountType::Mass>, get_enthalpy_liq<AmountType::Mass>, result);

    result.td_functions.inner_energy.molar = flash_function<AmountType::Molar>(
                                    get_energy_gas<AmountType::Molar>, get_energy_liq<AmountType::Molar>, result);
    result.td_functions.inner_energy.mass = flash_function<AmountType::Mass>(
                                   get_energy_gas<AmountType::Mass>, get_energy_liq<AmountType::Mass>, result);

    /*
    // Дальше три строчки ахтунга для однофазных случаев! В двухфазке проще

    // Копируем жидкость без мемоизации, только состав.
    // Это потом попадет в мемоизацию fluid_liquid
    result.fluid_liquid = create_copy(false);
    result.fluid_vapor = nullptr;
    // в мемоизации fluid_liquid есть свой fluid_liquid,
    // он будет копией текущего флюида без мемоизации
    // сам же fluid_liquid будет с мемоизацией
    result.fluid_liquid = create_copy();

    */
    return result;
}

const vlelib::flash_calculation_result_t fluid_rault_dalton_t::flash_vapor_only(
        double pressure, double temperature, flash_calculation_result_t& result) const
{
    result.pressure = pressure;
    result.temperature = temperature;
    result.molar_mass.mix = get_molar_mass();
    result.critical_temperature = get_pseudocritical_temperature();
    result.critical_pressure = get_pseudocritical_pressure();
    result.flash = 1;
    if (result.temperature >= result.critical_temperature &&
            result.pressure >= result.critical_pressure) {
        result.state_of_matter = state_of_matter_t::Critical;
    }
    else {
        result.state_of_matter = state_of_matter_t::Gas;
    }
    result.vapor_volumetric_fraction = 1;
    result.vapor_mass_fraction = 1;

    {
        size_t components_count = get_components_count();
        auto molar_fraction = get_molar_fraction();

        Eigen::VectorXd densities_vapor = get_densities_vapor(pressure, temperature);
        if (densities_vapor.sum() != 0)
        {
            Eigen::VectorXd molar_masses = get_molar_masses();
            Eigen::VectorXd molar_volume_gas = molar_fraction.cwiseProduct(molar_masses).cwiseProduct(densities_vapor.cwiseInverse());
            result.molar_volume.mix = molar_volume_gas.sum();

            Eigen::VectorXd volume_fraction(components_count);
            volume_fraction = molar_volume_gas / result.molar_volume.mix;

            result.density.mix = result.density.vapor = volume_fraction.dot(densities_vapor);
        }
        else //densities_vapor  тождественный 0
        {
            result.density.mix = result.density.vapor = 0;
        }
    }

    result.z_factor.vapor = 1.0;

    result.fluid_liquid = nullptr;
    result.fluid_vapor = create_copy(false);

    //flash_enthalpy(pressure, temperature);

    result.td_functions.enthalpy.molar = flash_function<AmountType::Molar>(
                                get_enthalpy_gas<AmountType::Molar>, get_enthalpy_liq<AmountType::Molar>, result);
    result.td_functions.enthalpy.mass = flash_function<AmountType::Mass>(
                               get_enthalpy_gas<AmountType::Mass>, get_enthalpy_liq<AmountType::Mass>, result);



    result.td_functions.inner_energy.molar = flash_function<AmountType::Molar>(
                                    get_energy_gas<AmountType::Molar>, get_energy_liq<AmountType::Molar>, result);
    result.td_functions.inner_energy.mass = flash_function<AmountType::Mass>(
                                   get_energy_gas<AmountType::Mass>, get_energy_liq<AmountType::Mass>, result);
    /*
    // ВОТ ТУТ САМОЕ БОЛЬШОЕ ЗАБЛУЖДЕНИЕ. РЕЗАЛТ ТУТ НЕ СФОРМИРОВАН!!!!!!!!!!!!!!!
    // а я искал у себя почему у меня inf
    // Обязательно копируем в самом конце, только по окончании формирования result

    // Дальше три строчки ахтунга для однофазных случаев (жидкость, газ)! 8-)))))) тут полный ахтунг
    // В двухфазке с мемоизацией проще

    // Копируем газ без мемоизации, только состав.
    // Это потом попадет в мемоизацию fluid_liquid
    result.fluid_liquid = nullptr;
    result.fluid_vapor = create_copy(false);
    // в мемоизации fluid_vapor есть свой fluid_vapor,
    // он будет копией текущего флюида без мемоизации
    // сам же fluid_vapor будет с мемоизацией
    result.fluid_vapor = create_copy();
    */


    return result;
}


void fluid_rault_dalton_t::build_twophase_result2(double pressure, double temperature, 
    const rachford_rice_result_t& rr_result, flash_calculation_result_t& result) const
{
    two_phase_result_builder temp_obj;

    temp_obj.build(pressure, temperature, rr_result, *this, result);
}


void two_phase_result_builder::build(double pressure, double temperature, 
    const rachford_rice_result_t& rr_result, const fluid_rault_dalton_t& flash_in, 
    flash_calculation_result_t& result) const
{
    // Расчёт в данной функции проводится для получаемых параметров состояния 
    // (давления и температуры) и прописанного до вызова этой функции состава.
    result.pressure = pressure;
    result.temperature = temperature;

    // На основании входных данных в rachford_rice_result_t создаются объекты с
    // составами жидкой и газовой фаз. Владение передается в shared_ptr<fluid_t>
    // объекта result.
    // Поскольку create_copy вызывается с изменным составом,
    // мемоизация копироваться не будет
    result.fluid_vapor = flash_in.create_copy(rr_result.y);
    result.fluid_liquid = flash_in.create_copy(rr_result.x);

    // давление насыщенных паров жидкости (функция использует состав из this) S. Skogestad flash calculations.pdf
    result.saturated_pressure_liquid = flash_in.get_saturated_pressures(temperature).dot(rr_result.x);
    set_critical_parameters(flash_in, result);

    // мольная доля отгона (принимается,что уже ограниченная сверху и снизу)
    const double omega = result.flash = rr_result.vapor_split;
    intermediate_twophase_data intermediate(rr_result
        , flash_in.get_densities_vapor(pressure, temperature)
        , flash_in.get_densities_liquid(pressure, temperature)
        , flash_in.get_molar_masses()
    );

    result.molar_mass = get_molar_masses(intermediate, flash_in);
    result.molar_volume = get_molar_volumes(rr_result.vapor_split, intermediate, flash_in);
    result.density = get_densities(result.molar_mass, result.molar_volume);
    std::tie(result.vapor_mass_fraction, result.vapor_volumetric_fraction) =
        get_vapor_fractions(rr_result.vapor_split, result.molar_mass, result.molar_volume);

    result.td_functions.enthalpy = get_enthalpies(result, flash_in);
    result.td_functions.inner_energy = get_inner_energies(result, flash_in);
    result.z_factor.vapor = 1.0;

    // по сравнению с fluid_rault_dalton_t::build_twophase_result
    // отсутствует обработка случая (при всей его спорности)  
    //      if (rr_result.vapor_split >= 1)
}

void two_phase_result_builder::set_critical_parameters(const fluid_rault_dalton_t& flash_in, flash_calculation_result_t& result) const
{
    // скалярные произведения для получения пвсевдокритических температуры и давления
    // внутри функций используются явные циклы по компонентам
    result.critical_temperature = flash_in.get_pseudocritical_temperature();
    result.critical_pressure = flash_in.get_pseudocritical_pressure();
    if (result.temperature >= result.critical_temperature
        && result.pressure >= result.critical_pressure) {
        result.state_of_matter = state_of_matter_t::Critical;
    }
    else {
        result.state_of_matter = state_of_matter_t::TwoPhase;
    }
}
amounts_per_phase two_phase_result_builder::get_molar_masses(const intermediate_twophase_data& intermediate
    , const fluid_rault_dalton_t& flash_in) const
{
    amounts_per_phase molar_mass;
    // определяется молярная масса (скрытое (MatrixBase::dot) скалярное произведение 
    // вектора молярных масс компонентов и вектора мольных долей). Цепочка владения вектора 
    // мольных долей - fluid_composition_functions_t::fluid_t::fluid_rault_dalton_t ->
    // const fluid_fundamental_data_t& composition -> Eigen::VectorXd concentration_.
    // TODO (BDI): привести в соответствие описания функций / названия переменных связанных с
    // мольным составом (доли) и молярными концентрациями
    const double& Mmix = molar_mass.mix = flash_in.get_molar_mass();
    molar_mass.vapor = intermediate.yiMi.sum();
    molar_mass.liquid = intermediate.xiMi.sum();
    return molar_mass;
}
std::pair<double, double> two_phase_result_builder::get_vapor_fractions(
    double omega, const amounts_per_phase& molar_mass, const amounts_per_phase& molar_volume) const
{
    double vapor_mass_fraction = omega * molar_mass.vapor / molar_mass.mix;
    double vapor_volumetric_fraction = omega * molar_volume.vapor / molar_volume.mix;
    return std::make_pair(vapor_mass_fraction, vapor_volumetric_fraction);
}

amounts_per_phase two_phase_result_builder::get_molar_volumes(double omega
    , const intermediate_twophase_data& intermediate
    , const fluid_rault_dalton_t& flash_in) const
{
    amounts_per_phase molar_volume;
    // молярный объем газовой фазы (данный член структуры results не заполнялся в build1, 
    // кроме того, название molar_volume_vapor использовалось для вектора объемов 
    // компонентов в газе на 1 моль смеси)
    const double& v_vap = molar_volume.vapor = intermediate.w_vap.sum();
    // молярный объем жидкой фазы (данный член структуры results не заполнялся в build1, 
    // кроме того, название molar_volume_liquid использовалось для вектора объемов 
    // компонентов в жидкости на 1 моль смеси)
    const double& v_liq = molar_volume.liquid = intermediate.w_liq.sum();

    // молярный объем смеси (используется мольная доля отгона без ограничений)
    const double& v_mix = molar_volume.mix = omega * v_vap + (1 - omega) * v_liq;
    return molar_volume;
}


amounts_per_phase two_phase_result_builder::get_densities(const amounts_per_phase& molar_mass
    , const amounts_per_phase& molar_volume) const
{
    amounts_per_phase density;
    // молярный объем смеси (используется мольная доля отгона без ограничений)
    density.mix = molar_mass.mix / molar_volume.mix;
    density.vapor = molar_mass.vapor / molar_volume.vapor;
    density.liquid = molar_mass.liquid / molar_volume.liquid;
    return density;
}

amounts_molar_and_mass two_phase_result_builder::get_enthalpies(const flash_calculation_result_t& flash_result
    , const fluid_rault_dalton_t& flash_in) const {
    amounts_molar_and_mass result;
    result.molar = flash_in.flash_function<AmountType::Molar>(
        get_enthalpy_gas<AmountType::Molar>, get_enthalpy_liq<AmountType::Molar>, flash_result);
    result.mass = flash_in.flash_function<AmountType::Mass>(
        get_enthalpy_gas<AmountType::Mass>, get_enthalpy_liq<AmountType::Mass>, flash_result);
    return result;
}

amounts_molar_and_mass two_phase_result_builder::get_inner_energies(const flash_calculation_result_t& flash_result
    , const fluid_rault_dalton_t& flash_in) const {
    amounts_molar_and_mass result;
    result.molar = flash_in.flash_function<AmountType::Molar>(
            get_energy_gas<AmountType::Molar>, get_energy_liq<AmountType::Molar>, flash_result);
    result.mass = flash_in.flash_function<AmountType::Mass>(
        get_energy_gas<AmountType::Mass>, get_energy_liq<AmountType::Mass>, flash_result);
    return result;
}

void fluid_rault_dalton_t::flash_unsafe(double pressure, double temperature, double initial_estimation, flash_calculation_result_t& result) const
{
    if (temperature < 0)
        throw std::runtime_error("temperature < 0" + std::to_string(temperature));

    if (result.was_calculated(pressure, temperature)) {
        if (!std::isfinite(result.density.mix))
            throw std::runtime_error("flash_result.density.mix is nan");
        //++cached;
        return;
    }

    if (is_liquid_only(pressure, temperature)) {
        flash_liquid_only(pressure, temperature, result);
    }
    else if (is_vapor_only(pressure, temperature)) {
        flash_vapor_only(pressure, temperature, result);
    }
    else {

        rachford_rice2_t rr(this, pressure, temperature);

        fixed_solver_result_t<1> solver_result = rr.solve();
        if (solver_result.result_code != numerical_result_code_t::Converged) {
            pt_flash_stub_data_t data;
            data.fluid = get_mock_data();
            data.pressure = pressure;
            data.temperature = temperature;
            crashdump_and_throw("fluid_rault_dalton_t flash ", data);
            //throw std::runtime_error("fluid_rault_dalton_t::flash() RR not convered");
        }
        auto rr_result = rr.build_result(solver_result.argument);
        build_twophase_result2(pressure, temperature, rr_result, result);
    }
}

double fluid_rault_dalton_t::get_dew_point_at_given_temperature(double temperature) const
{
    const auto& components = get_components();
    auto molar_fraction = get_molar_fraction();

    double Pinv = 0;
    for (size_t index = 0; index < components.size(); ++index) {
        if (molar_fraction(index) != 0) {
            double Psat = components[index]->antoine_model.get_saturated_pressure(temperature);
            Pinv += molar_fraction(index) / Psat;
        }
    }
    return 1 / Pinv;
}


/// @brief Функция расчета температуры конденсации при заданном давлении
/// @param pressure Давление
/// @return Температуру конденсации
double fluid_rault_dalton_t::get_dew_point_at_given_pressure(double pressure) const
{
    const auto& components = get_components();
    size_t components_count = get_components_count();

    dew_point_at_given_pressure dew_point_calcultaion(pressure, this);
    fixed_solver_result_t<1> solver_result;

    if (components_count != 1)
    {
        return dew_point_calcultaion.solve(&solver_result);
    }

    return components.front()->antoine_model.get_temperature_for_saturated_pressure(pressure);
}

double fluid_rault_dalton_t::get_water_dew_point_at_given_pressure(double pressure) const
{
    double ret=std::numeric_limits<double>::quiet_NaN();
    for(size_t i=0;i<get_components_count();++i){
        auto c=get_components()[i];
        auto mf = get_molar_fraction()(i);
        if(c->name==L"H2O" && mf>std::numeric_limits<double>::epsilon()){
            auto tsat = c->antoine_model.get_temperature_for_saturated_pressure(pressure*mf);
            ret = tsat;
        }
    }
    //throw std::logic_error("Не реализовано");
    return ret;
}

double fluid_rault_dalton_t::get_bubble_point_at_given_temperature(double temperature) const
{
    const auto& components = get_components();
    auto molar_fraction = get_molar_fraction();

    double P = 0;
    for (size_t index = 0; index < components.size(); ++index) {
        if (molar_fraction(index) != 0) {
            double Psat = components[index]->antoine_model.get_saturated_pressure(temperature);
            P += molar_fraction(index) * Psat;
        }
    }
    return P;
}

/// @brief Функция расчета температуры кипения при заданном давлении
/// @param pressure Давление
/// @return Температуру кипения
double fluid_rault_dalton_t::get_bubble_point_at_given_pressure(double pressure) const
{
    const auto& components = get_components();
    size_t components_count = get_components_count();

    bubble_point_at_given_pressure bubble_point_calcultaion(pressure, this);
    fixed_solver_result_t<1> solver_result;

    if (components_count != 1)
    {
        return bubble_point_calcultaion.solve(&solver_result);
    }

    return components.front()->antoine_model.get_temperature_for_saturated_pressure(pressure);
}

/// @brief Начальное приближение по температуре при известном давлении
/// В смеси из компонентов температура будет где-то посредине между максимальной и минимальной
double temperature_initial_estimate(double pressure, const std::vector<const component_properties_t*>& components)
{
    std::vector<double> saturated_temperatures;
    for (const auto& component : components) {
        saturated_temperatures.push_back(
                    component->antoine_model.get_temperature_for_saturated_pressure(pressure));
    }
    double min_temperature = *std::min_element(saturated_temperatures.begin(), saturated_temperatures.end());
    double max_temperature = *std::max_element(saturated_temperatures.begin(), saturated_temperatures.end());
    double result = 0.5 * (min_temperature + max_temperature);

    return result;
}
void fluid_rault_dalton_t::change_liquid_volume_fraction(
        double pressure, double temperature, double liquid_volume_fraction)
{
    auto& f = flash(pressure, temperature);
    if (f.is_gas_only() || f.is_liquid_only())
        throw std::logic_error("cannot change liquid volume fraction");

    // Количество вещества по фазам, занимающем объем в метрах, численно равный объемной доле
    double molar_amount_liquid =
            f.density.liquid * liquid_volume_fraction / f.fluid_liquid->get_molar_mass();
    double molar_amount_gas =
            f.density.vapor * (1 - liquid_volume_fraction) / f.fluid_vapor->get_molar_mass();

    // Общее количество вещества по компонентам
    const Eigen::VectorXd& x = f.fluid_liquid->get_molar_fraction();
    const Eigen::VectorXd& y = f.fluid_vapor->get_molar_fraction();
    Eigen::VectorXd amount_vector = x * molar_amount_liquid + y * molar_amount_gas;

    //molar_fraction = amount_vector / amount_vector.sum();
    set_molar_fraction(amount_vector / amount_vector.sum());
    // Invalidate last result;
    //last_flash_result.invalidate_calculation();
    const auto f2 = flash(pressure, temperature);// избыточно?

}

double fluid_rault_dalton_t::fill_and_check_liquid_accumulation(
        double volume, double temperature, double molar_amount) const
{
    // Факическая плотность смеси, по входным данным
    double liquid_mix = get_molar_mass() * molar_amount / volume;

    // Давления начала кипения и соответствующая плотность жидкой смеси
    double Pbubble = get_bubble_point_at_given_temperature(temperature);
    double liquid_density_border = get_density_as_liquid(Pbubble, temperature);

    if (liquid_mix >= liquid_density_border) {
        double eps = 1e-5;
        double dp = Pbubble * eps;
        double drho =
                get_density_as_liquid(Pbubble + dp, temperature)
                - get_density_as_liquid(Pbubble, temperature);
        double derivative = dp / drho; // именно так, чувствительность давления к плотности!
        double result = Pbubble + derivative * (liquid_mix - liquid_density_border);
        // Плотность смеси больше граничной, полностью жидкая фаза
        return result;
    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
}


double fluid_rault_dalton_t::fill_volume_with_liquid(
    double volume, double temperature, double liquid_volume) const
{
    // Старый вызов (закомментирован)
    //desired_liquid_vol_fraction_equation solver(volume, temperature, liquid_volume, this);
    //return solver.solve();

    // Новый вызов на основе fixed_system
    desired_liquid_vol_fraction_equation_fixed solver(volume, temperature, liquid_volume, this);
    return solver.solve();
}

double fluid_rault_dalton_t::fill_volume_with_pressure(
        double pressure, double temperature, double total_volume) const
{
    double Pbubble = get_bubble_point_at_given_temperature(temperature);

    double density_mix;
    if (pressure >= Pbubble) {
        double liquid_density_border = get_density_as_liquid(Pbubble, temperature);

        double eps = 1e-4;
        double dp = Pbubble * eps;
        double drho = get_density_as_liquid(Pbubble + dp, temperature) - get_density_as_liquid(Pbubble, temperature);
        double drho_dp = drho / dp; // именно так, чувствительность давления к плотности!

        density_mix = liquid_density_border + drho_dp * (pressure - Pbubble);
    }
    else {
        auto vle = flash(pressure, temperature);
        density_mix = vle.density.mix;
    }

    double M_mix = get_molar_mass();
    double result = density_mix * total_volume / M_mix;
    return result;


}





std::pair<double, double> fluid_rault_dalton_t::fill_volume_with_total_moles2(
        double volume, double inner_energy_molar, double molar_amount,
        double initial_pressure, double initial_temperature) const
{
    double Tliq = find_liquid_temperature_with_inner_energy<AmountType::Molar>(inner_energy_molar);

    //auto HH = get_enthalpy_td_mass_as_liquid(1e5, Tliq);

    // Факическая плотность смеси, по входным данным
    double target_density = get_molar_mass() * molar_amount / volume;

    double Pliq_bubble = get_bubble_point_at_given_temperature(Tliq);
    double density_liquid_border = get_density_as_liquid(Pliq_bubble, Tliq);

    if (target_density >= density_liquid_border) {
        // только жидкость
        //double eps = 1e-5;
        //double dp = Pliq_bubble * eps;
        //double drho = get_density_as_liquid(Pliq_bubble + dp, Tliq) -
        //    get_density_as_liquid(Pliq_bubble - dp, Tliq);
        //double derivative = 2 * dp / drho; // именно так, чувствительность давления к плотности!

        // именно так, чувствительность давления к плотности!
        double derivative = 1 / get_density_as_liquid_pressure_derivative(Pliq_bubble, Tliq);

        double Pliq = Pliq_bubble + derivative * (target_density - density_liquid_border);
        // Плотность смеси больше граничной, полностью жидкая фаза
        return { Pliq, Tliq };
    }
    else {
        // Используем Tliq как начальное приближение, если initial_temperature не задана
        double initial_temp = std::isfinite(initial_temperature) ? initial_temperature : Tliq;

        double min_temperature = vlelib::get_min_antoine_bound(this);
        if (get_ideal_gas_inner_energy<AmountType::Molar>(min_temperature) < inner_energy_molar) {
            /*double Tgas2 = estimate_temperature_for_inner_energy<AmountType::Molar>(
                this, inner_energy_molar, initial_temp, 1e5);*/

            double Tgas = estimate_temperature_for_inner_energy_fixed<AmountType::Molar>(
                this, inner_energy_molar, initial_temp);

            // todo: обработка сходимости
            if (!std::isfinite(Tgas))
                throw std::logic_error("vle_desired_ideal_gas_enthalpy_equation not converged");
            double Pgas = molar_amount / volume * M_R * Tgas;
            double Pgas_dew = get_dew_point_at_given_temperature(Tgas);

            if (Pgas < Pgas_dew) {
                return { Pgas, Tgas };
            }
        }
    }

    // Если не газ и не жидкость, т.е. двухфазный расчет
    {
        uv_flash_over_RR uv(volume, molar_amount, inner_energy_molar, this,
                            false, initial_pressure, initial_temperature);

        auto uv_result = uv.solve();
        if (uv_result.result_code == numerical_result_code_t::NotConverged) {
            // Создается структура под сериализацию
            crashdump_and_throw("UV-flash two phase not converged", uv.get_mock_data());
            //throw std::runtime_error("UV-flash two phase not converged");
        }

        auto [Pmix, Tmix] = uv.denorm(uv_result.argument[0], uv_result.argument[1]);
        return std::make_pair(Pmix, Tmix);
    }
}

template <AmountType amount_type>
double fluid_rault_dalton_t::find_liquid_temperature_with_inner_energy(double inner_energy) const
{
    const auto& components = get_components();
    size_t components_count = get_components_count();
    auto molar_fraction = get_molar_fraction();


    if constexpr (amount_type == AmountType::Mass) {
        inner_energy = inner_energy * get_molar_mass();
    }

    double Sum1 = 0;
    double Sum2 = 0;
    for (size_t index = 0; index < components_count; ++index) {
        double zi = molar_fraction(index);
        // реально теплоемкость не зависит от температуры!
        double Ci = components[index]->heat_capacity_liquid.get_polynom_value(300);
        double Tbi = components[index]->normal_boiling_temperature;
        double Qi = components[index]->condensation_heat_molar;
        double Hvi = components[index]->get_inner_energy_gas<AmountType::Molar>(Tbi);

        Sum1 += zi * (Hvi - Qi - Ci * Tbi);
        Sum2 += zi * Ci;
    }
    double num = (inner_energy - Sum1);
    double r = num / Sum2;
    double result = (inner_energy - Sum1) / Sum2;


    return result;
}

/// @brief Эмпирический коэффициент, на который домножается нижнее ограничение по приведенной
/// температуре. Его значение было подобрано так, чтобы сошелся состав из задачи MM-260.
/// дальнейшее снижение коэффициента нежелательно, требуется исследование обоснованного задания
/// границ по приведенной температуре.
constexpr double temperature_constraint_coeff = 0.8;

double get_min_antoine_bound(const fluid_t* fluid)
{
    const auto& comps = fluid->get_components();

    std::vector<double> min_antoine_bound;
    std::transform(comps.begin(), comps.end(), std::back_inserter(min_antoine_bound),
        [](const component_properties_t* c) { return c->antoine_model.min_bound; }
    );
    // учитываем вес компонентов
    double mean_min = fluid->get_molar_fraction().dot(
        Eigen::VectorXd::Map(min_antoine_bound.data(), min_antoine_bound.size()));
    // а тут не учитываем весь компонентов, просто берем минимум,
    // сколько бы ни было его в смеси
    double min = *std::min_element(min_antoine_bound.begin(), min_antoine_bound.end());
    
    return mean_min * temperature_constraint_coeff;
}

double get_max_antoine_bound(const fluid_rault_dalton_t* fluid)
{
    const auto& comps = fluid->get_components();

    std::vector<double> max_antoine_bound;
    std::transform(comps.begin(), comps.end(), std::back_inserter(max_antoine_bound),
        [](const component_properties_t* c) { return c->antoine_model.max_bound; }
    );
    // учитываем вес компонентов
    double mean_max = fluid->get_molar_fraction().dot(
        Eigen::VectorXd::Map(max_antoine_bound.data(), max_antoine_bound.size()));
    // а тут не учитываем весь компонентов, просто берем минимум,
    // сколько бы ни было его в смеси
    double max = *std::min_element(max_antoine_bound.begin(), max_antoine_bound.end());
    return mean_max;
}

//
//template<> double fluid_rault_dalton_t::find_liquid_temperature_with_inner_energy<AmountType::Molar>(double inner_energy) const;
//template<> double fluid_rault_dalton_t::find_liquid_temperature_with_inner_energy<AmountType::Mass>(double inner_energy) const;


intermediate_twophase_data::intermediate_twophase_data(
    const rachford_rice_result_t& rr_result, 
    const Eigen::VectorXd& densities_vapor, 
    const Eigen::VectorXd& densities_liquid, const Eigen::VectorXd& M)
    : densities_vapor(densities_vapor)
    , densities_liquid(densities_liquid)
    , M(M)
{
    // массы компонентов в газе на один моль газа
    yiMi = M.cwiseProduct(rr_result.y);
    // массы компонентов в жидкости на один моль жидкости
    xiMi = M.cwiseProduct(rr_result.x);
    // TODO (BDI): удостовериться, что rachford_rice_result_t передает составы по фазам,
    // после чего явно указать это в описании структуры.

    // объемы чистых компонентов по газу на один моль газа
    w_vap = yiMi.cwiseProduct(fixed_solvers::invVectorXd(densities_vapor));
    // объемы чистых компонентов по жидкости на один моль жидкости
    w_liq = xiMi.cwiseProduct(fixed_solvers::invVectorXd(densities_liquid));
}

}


