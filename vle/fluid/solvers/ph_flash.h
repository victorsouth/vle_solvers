#pragma once

#ifndef __vlelib_fluid_equations_h__
#error "Do not include ph_flash.h directly. Use fluid_equations.h instead."
#endif
namespace vlelib {
;


/// @brief Данные для изолированного вызова задачи PH-flash
struct ph_flash_stub_data_t
{
    /// @brief Давление задачи PH-flash
    double pressure;
    /// @brief Энтальпия
    double enthalpy_mass;
    /// @brief Данные для изолированного вызова по смеси
    fluid_stubdata_t fluid;

#ifdef VLELIB_SERIALIZATION_SUPPORT
    /// @brief Для сериализации
    friend class boost::serialization::access;

    /// @brief Запись в заданный файл с помощью сериализатора нужного формата
    template <typename Serializer = boost::archive::xml_oarchive>
    void to_file(const std::string& filename = "ph_flash_failed.xml") const {
        std::ofstream ofs(filename);
        if (ofs.is_open()) {
            boost::archive::xml_oarchive oa(ofs);
            const auto& mock_data = *this;
            oa << BOOST_SERIALIZATION_NVP(mock_data);
        }
        else {
            throw std::runtime_error("Failed to write PH-flash mockdata");
        }
        ofs.flush();
        ofs.close();
    }
    /// @brief Интрузивная сериализация/десериализация
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_NVP(pressure);
        ar& BOOST_SERIALIZATION_NVP(enthalpy_mass);
        ar& BOOST_SERIALIZATION_NVP(fluid);
    }
#endif
};


/// \brief PH flash методом бисекции
/// упрощённый вариант
/// todo: необходимо проработать ситуацию, когда
/// предельная искомая энтальпия ниже
/// желаемого значения
struct desired_enthalpy : public fixed_system_t<1> {
    /// \brief исследуемый флюид
    fluid_t* fluid;
    /// \brief энтальпия для которой необходимо найти температуру флюида
    double given_enthalpy;
    /// \brief заданное давление флюида
    double given_pressure;
    /// \brief начальное приближение по температуре,
    /// в общем случае может быть не задано
    double initial_temperature;
public:
    /// \brief функция вычисления невязки
    /// \param x - аргумент
    /// \return невязка
    virtual function_type residuals(const var_type& x) override
    {
        return given_enthalpy - fluid->flash(given_pressure, x).td_functions.enthalpy.mass.mix;
    }
    /// \brief
    /// конструктор desired_enthalpy
    /// \param f исследуемый флюид
    /// \param de энтальпия для которой необходимо найти температуру флюида
    /// \param p заданное давление флюида
    /// \param t начальное приближение по температуре, в общем случае может быть не задано
    desired_enthalpy(fluid_t* f, double de, double p, double t = std::numeric_limits<double>::quiet_NaN())
        :fluid(f), given_enthalpy(de), given_pressure(p), initial_temperature(t)
    {}
    /// \brief
    /// решить PH flash методом бисекции
    /// \param numerical_result - результат
    /// \param analysis - анализ
    /// \param desired_precision - желаемая точность по невязке
    /// \return
    double solve(
        double desired_precision = std::numeric_limits<float>::epsilon(),
        fixed_bisection_result_t<1>* numerical_result = nullptr,
        fixed_bisection_result_analysis_t<1>* analysis = nullptr
    )
    {
        fixed_bisection_result_t<1> res;
        if (numerical_result == nullptr)
            numerical_result = &res;

        fixed_bisectional_parameters_t p;
        p.argument_limit_min = 10;
        p.argument_limit_max = 5000;
        // correcting limits by estimation if applicable
        if (std::isfinite(initial_temperature)) {
            double enthalpy_at_initial_temperature = fluid->flash(given_pressure, initial_temperature).td_functions.enthalpy.mass.mix;
            if (enthalpy_at_initial_temperature < given_enthalpy)p.argument_limit_min = initial_temperature;
            if (enthalpy_at_initial_temperature > given_enthalpy)p.argument_limit_max = initial_temperature;
        }
        double T_critical = fluid->get_pseudocritical_temperature();
        {
            double enthalpy_at_critical_temperature = fluid->flash(given_pressure, T_critical).td_functions.enthalpy.mass.mix;
            if (enthalpy_at_critical_temperature < given_enthalpy)p.argument_limit_min = T_critical;
            if (enthalpy_at_critical_temperature > given_enthalpy)p.argument_limit_max = T_critical;
        }
        p.argument_precision = desired_precision;
        p.residual_precision = desired_precision;
        p.argument_history = true;
        p.residual_history = true;
        p.solution_type = fixed_bisectional_solution_type::Combined;
        p.secant_treshhold_iterations=0;
        //p.secant_treshhold_max = 10;
        //p.secant_treshhold_min = 0.001;
        p.verbose = false;
        fixed_bisectional<1>::solve(p, *this, numerical_result, analysis);

        return numerical_result->argument;
    }
    /// \brief
    /// решить PH flash методом бисекции
    /// \param numerical_result - результат
    /// \param analysis - анализ
    /// \param desired_precision - желаемая точность по невязке
    /// \return
    double solve(
        fixed_bisection_result_t<1>* numerical_result = nullptr,
        fixed_bisection_result_analysis_t<1>* analysis = nullptr
    ) {
        return solve(std::numeric_limits<double>::epsilon() * 10000., numerical_result, analysis);
    }
    /// @brief Возвращает данные для изолированного вызова
    ph_flash_stub_data_t get_stub_data() const {
        ph_flash_stub_data_t result;
        result.pressure = given_pressure;
        result.enthalpy_mass = given_enthalpy;
        result.fluid = fluid->get_mock_data();
        return result;
    }
};

/// @brief Уравнение для поиска температуры по заданной энтальпии 
/// с учетом фазовых переходов
struct desired_enthalpy_fixed : fixed_system_t<1>
{
    /// @brief Флюид, для которого выполняется расчет
    const fluid_t* fluid;
    /// @brief Целевая удеальная массова энтальпия
    double target_enthalpy;
    /// @brief Давление для расчетной энтальпии
    double pressure;
    /// @brief Начальное приближение по температуре
    double temperature_initial;
    /// @brief Нижняя граница по температуре
    double temperature_minimum{ 0.0 };
    /// @brief Верхняя граница по температуре
    double temperature_maximum{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Точность расчета
    double desired_precission;
public:
    /// @brief Энтальпия задаётся как целевой параметр
    desired_enthalpy_fixed(const fluid_t* fluid, double pressure, double target_enthalpy_mass,
        double _temperature_initial = std::numeric_limits<double>::quiet_NaN(),
        double _prec = std::numeric_limits<double>::epsilon() * 10000.)
        : fluid(fluid)
        , target_enthalpy(target_enthalpy_mass)
        , pressure(pressure)
        , temperature_initial(_temperature_initial)
        , desired_precission(_prec)
    {
        if (!std::isfinite(temperature_initial))
        {
            auto fluid_rd = dynamic_cast<const fluid_rault_dalton_t*>(fluid);
            if (fluid_rd == nullptr)
                throw std::runtime_error("Can estimate min temperature only for Raoult-fluid Dalton");
            double T_critical = fluid->get_pseudocritical_temperature();
            double T_initial = T_critical;
            double T_min = get_min_antoine_bound(fluid_rd);
            temperature_initial = std::max(T_min, 0.5 * T_initial);

        }
    }
    /// @brief Энтальпия расчитывается по начальной, подводимому теплу и вносу энтальпии расходом 
    desired_enthalpy_fixed(const fluid_t* fluid, double p,
        double initial_entalpy, double Q, double mass_flow,
        double temperature_initial)
        : desired_enthalpy_fixed(fluid, p, initial_entalpy + Q / mass_flow, temperature_initial)
    {

    }
    /// @brief Невязки уравнения
    virtual double residuals(const double& temperature) override
    {
        const auto vle = fluid->flash(pressure, temperature);
        double r = vle.td_functions.enthalpy.mass.mix - target_enthalpy;
        return r;
    }

    /// @brief Оптимизация расчета, учитывающая, что jacobian запускается после residuals
    virtual double jacobian_dense(const double& temperature) override
    {
        double e = epsilon * std::max(1.0, abs(temperature));
        const auto vle = fluid->flash(pressure, temperature);//= fluid->get_last_flash_result();//спорно
        double mass_enthalpy_0 = vle.td_functions.enthalpy.mass.mix;

        const auto vle2 = fluid->flash(pressure, temperature + e);
        double mass_enthalpy_eps = vle2.td_functions.enthalpy.mass.mix;

        function_type J = (mass_enthalpy_eps - mass_enthalpy_0) / e;
        return J;
    }

    /// @brief Возвращает начальное приближение по температуре
    var_type estimation() const {
        return temperature_initial;
    }

    /// @brief Расчет PH-flash
    /// @param numerical_result 
    /// @param analysis_result 
    /// @return Искомая температура или NaN, если расчет не сошелся
    double solve(fixed_solver_result_t<1>* numerical_result = nullptr,
        fixed_solver_result_analysis_t<1>* analysis_result = nullptr)
    {
        fixed_solver_result_t<1> solver_result_carrier;
        if (numerical_result == nullptr)
            numerical_result = &solver_result_carrier;


        double T_initial = estimation();
        fixed_solver_parameters_t<1, 0, golden_section_search> solver_parameters;
        solver_parameters.line_search.function_decrement_factor = 10;
        solver_parameters.line_search.iteration_count = 30;
        solver_parameters.line_search_fail_action = line_search_fail_action_t::TreatAsFail;

        solver_parameters.constraints.minimum = temperature_minimum;
        solver_parameters.constraints.maximum = temperature_maximum;

        solver_parameters.argument_increment_norm = desired_precission;
        if (analysis_result != nullptr) {
            solver_parameters.analysis.argument_history = true;
            solver_parameters.analysis.steps = true;
            solver_parameters.analysis.line_search_explore = true;
        }

        fixed_newton_raphson<1>::solve_dense(
            *this, T_initial, solver_parameters, numerical_result, analysis_result);

        if (numerical_result->result_code == numerical_result_code_t::Converged) {
            return numerical_result->argument;
        }
        else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    /// @brief Возвращает данные для изолированного вызова
    ph_flash_stub_data_t get_stub_data() const {
        ph_flash_stub_data_t result;
        result.pressure = pressure;
        result.enthalpy_mass = target_enthalpy;
        result.fluid = fluid->get_mock_data();
        return result;
    }
};

/* НЕИСПОЛЬЗУЕМЫЙ КОД: Заменен на vle_desired_fluid_function_equation_fixed (на основе fixed_system_t)
/// @brief Уравнение для поиска температуры по заданному значению функции флюида
/// @tparam Function Тип функции флюида, принимающей (const fluid_t*, double) и возвращающей double
template <typename Function>
struct vle_desired_fluid_function_equation : num_methods::equation_system_t {
    /// @brief Флюид, для которого выполняется расчет
    const fluid_t* fluid;
    /// @brief Целевая энтальпия
    double target_value;
    /// @brief Функция флюида
    Function& function;
public:
    /// @brief Энтальпия задаётся как целевой параметр
    vle_desired_fluid_function_equation(
        const fluid_t* _fluid, double _target_value, Function _function)
        : fluid(_fluid)
        , target_value(_target_value)
        , function(_function)
    {}

    /// @brief Возвращает нормализованную невязку из residuals
    virtual double operator() (const Eigen::VectorXd& value) override {
        return residuals(value).norm();
    }

    /// @brief Возвращает невязку между целевой энтальпией и расчёту по функции
    virtual Eigen::VectorXd residuals(const Eigen::VectorXd& value) override {
        const double& T = value(0);
        // здесь при расчете энтальпии фиктивно отдаем атмосферное давление, 
        // т.к. идеально-газовая энтальпия зависит только от температуры
        double function_value = function(fluid, T);

        double diff = function_value - target_value;
        Eigen::VectorXd r(1);
        r(0) = diff;
        return r;
    }
    /// @brief Возвращает размерность системы уравнений
    /// @return Размерность системы (1)
    virtual size_t dimension() const override { return 1; }
    /// @brief Возвращает размерность аргумента
    /// @return Размерность аргумента (1)
    virtual size_t argument_dimension() const override { return 1; }
};
*/



/* НЕИСПОЛЬЗУЕМЫЙ КОД: Заменен на desired_liquid_vol_fraction_equation_fixed (на основе fixed_system_t)
/// @brief Уравнение для поиска давления по заданной объемной доле жидкости
/// на основе equation_system_t (старая реализация)
struct desired_liquid_vol_fraction_equation : num_methods::equation_system_t {
    /// @brief Общий объем системы
    double volume;
    /// @brief Температура системы
    double temperature;
    /// @brief Объем жидкости, который должен быть достигнут
    double liquid_volume;
    /// @brief Флюид, для которого выполняется расчет
    const fluid_t* fluid;

    /// @brief Конструктор уравнения
    /// @param volume Общий объем системы
    /// @param temperature Температура системы
    /// @param liquid_volume Объем жидкости, который должен быть достигнут
    /// @param fluid Флюид, для которого выполняется расчет
    desired_liquid_vol_fraction_equation(double volume, double temperature, double liquid_volume, const fluid_t* fluid)
        : volume(volume)
        , temperature(temperature)
        , liquid_volume(liquid_volume)
        , fluid(fluid)
    {
        epsilon = 1e-4;

        if (liquid_volume > volume) {
            throw std::logic_error("liquid volume must be less than total volume");
        }
        if (liquid_volume < 0) {
            throw std::logic_error("liquid volume must be nonnegative");
        }

    }
    /// @brief Начальное приближение (не реализовано)
    /// @throws logic_error Всегда выбрасывает исключение, т.к. начальное приближение нетривиально
    Eigen::VectorXd estimation() const override {
        throw std::logic_error("nontrivial estimation");
    }
    /// @brief Невязки уравнения
    /// @param w Вектор аргументов, содержащий давление
    /// @return Вектор невязок (размерности 1)
    virtual Eigen::VectorXd residuals(const Eigen::VectorXd& w) override {
        const double& pressure = w(0);
        Eigen::VectorXd r(1);

        const auto flash = fluid->flash(pressure, temperature);

        double liquid_frac = 1. - flash.vapor_volumetric_fraction;
        //r(0) = liquid_frac - liquid_volume / volume;
        r(0) = liquid_frac * volume - liquid_volume;
        return r;
    }
    /// @brief Размерность системы уравнений
    /// @return Размерность (1)
    virtual size_t dimension() const override {
        return 1;
    }

    /// @brief Размерность аргумента
    /// @return Размерность аргумента (1)
    virtual size_t argument_dimension() const override {
        return 1;
    }
    /// @brief Костылим расчет производной.
    /// В областях, близких однофазным производная обнуляется, что не дает считать численный метод
    virtual Eigen::SparseMatrix<double> jacobian(const Eigen::VectorXd& argument) override
    {
        double Pdew = fluid->get_dew_point_at_given_temperature(temperature);
        double Pbubble = fluid->get_bubble_point_at_given_temperature(temperature);

        constexpr double eps = 1e-6;
        // Проверка на верхнюю границу
        if (argument(0) > Pbubble * (1 - eps)) {

            Eigen::VectorXd arg(1);
            arg(0) = Pbubble * (1 - eps);
            return equation_system_t::jacobian(arg);
        }
        // Проверка на нижнюю границу
        else if (argument(0) < Pdew * (1 + eps)) {
            Eigen::VectorXd arg(1);
            arg(0) = Pdew * (1 + eps);
            return equation_system_t::jacobian(arg);

        }
        // Нормальный, обычный расчет
        return equation_system_t::jacobian(argument);
    }

    /// @brief Возвращает веса для невязок
    /// @return Вектор весов (размерности 1)
    virtual vector<double> get_gains() const override {
        return { 1e3 };
    }
    /// @brief Решение уравнения методом Гаусса-Ньютона
    /// @return Давление, при котором достигается заданная объемная доля жидкости
    double solve() {
        double Pdew = fluid->get_dew_point_at_given_temperature(temperature);
        double Pbubble = fluid->get_bubble_point_at_given_temperature(temperature);
        if (liquid_volume / volume < 4.65e-06)//уровень в тесте ~0.05 мм
            return Pdew;
        if ((1. - liquid_volume / volume) < 6.52e-06)//уровень в тесте ~9899.93 мм при высоте 9.9м
            return Pbubble;


        num_methods::equation_solver_parameters_t solver_parameters;
        solver_parameters.constraints.minimum[0] = Pdew;
        solver_parameters.constraints.maximum[0] = Pbubble;
        solver_parameters.constraints.relative_boundaries[0] = (Pbubble - Pdew) * 0.05;

        solver_parameters.argument_increment_norm = 1e-7;

        Eigen::VectorXd estimation(1);
        estimation(0) = (Pdew + Pbubble) * liquid_volume / volume;// так поближе
        //estimation(0) = (Pdew + Pbubble) / 2;
        double dp = Pbubble - Pdew;
        estimation(0) = std::max(Pdew + 1e-3 * dp, estimation(0));
        estimation(0) = std::min(Pbubble - 1e-3 * dp, estimation(0));

        if (fabs(residuals(estimation)(0)) < 1e-4)
            return estimation(0);

        num_methods::equation_solver_result_t solver_result;

        // Южанин 2022-11-03: Почему 6, почему 23.5 ?
        // Тут переработка нужна

        num_methods::solve_gauss_newton_bounded(*this, estimation, solver_parameters, &solver_result);
        int i = 0;
        while (!solver_result.converged && ++i < 6) {// если не сошлось за раз на тестах i(max)=3, возможно надо поправить параметры солвера.
            solver_parameters.line_search_min_step /= 23.5;
            try {
                num_methods::solve_gauss_newton_bounded(*this, estimation, solver_parameters, &solver_result);
            }
            catch (...) {
                continue;
            }
        }
        if (!solver_result.converged)
            throw std::runtime_error("gauss failed");
        return solver_result.argument(0);
    }
};
*/

/// @brief Уравнение для поиска давления по заданной объемной доле жидкости
/// на основе fixed_system_t
struct desired_liquid_vol_fraction_equation_fixed : fixed_system_t<1> {
    /// @brief Общий объем системы
    double volume;
    /// @brief Температура системы
    double temperature;
    /// @brief Объем жидкости, который должен быть достигнут
    double liquid_volume;
    /// @brief Флюид, для которого выполняется расчет
    const fluid_t* fluid;

    /// @brief Конструктор уравнения
    /// @param volume Общий объем системы
    /// @param temperature Температура системы
    /// @param liquid_volume Объем жидкости, который должен быть достигнут
    /// @param fluid Флюид, для которого выполняется расчет
    desired_liquid_vol_fraction_equation_fixed(double volume, double temperature, double liquid_volume, const fluid_t* fluid)
        : volume(volume)
        , temperature(temperature)
        , liquid_volume(liquid_volume)
        , fluid(fluid)
    {
        if (liquid_volume > volume) {
            throw std::logic_error("liquid volume must be less than total volume");
        }
        if (liquid_volume < 0) {
            throw std::logic_error("liquid volume must be nonnegative");
        }
    }

    /// @brief Невязки уравнения
    virtual double residuals(const double& pressure) override {
        const auto flash = fluid->flash(pressure, temperature);
        double liquid_frac = 1. - flash.vapor_volumetric_fraction;
        return liquid_frac * volume - liquid_volume;
    }

    /// @brief Костылим расчет производной.
    /// В областях, близких однофазным производная обнуляется, что не дает считать численный метод
    virtual double jacobian_dense(const double& pressure) override
    {
        double Pdew = fluid->get_dew_point_at_given_temperature(temperature);
        double Pbubble = fluid->get_bubble_point_at_given_temperature(temperature);

        constexpr double eps = 1e-6;
        // Проверка на верхнюю границу
        if (pressure > Pbubble * (1 - eps)) {
            return fixed_system_t<1>::jacobian_dense(Pbubble * (1 - eps));
        }
        // Проверка на нижнюю границу
        else if (pressure < Pdew * (1 + eps)) {
            return fixed_system_t<1>::jacobian_dense(Pdew * (1 + eps));
        }
        // Нормальный, обычный расчет
        return fixed_system_t<1>::jacobian_dense(pressure);
    }

    /// @brief Решение уравнения методом Ньютона-Рафсона
    /// @param solver_result Результат работы численного метода (опционально)
    /// @return Давление, при котором достигается заданная объемная доля жидкости
    double solve(fixed_solver_result_t<1>* solver_result = nullptr) {
        double Pdew = fluid->get_dew_point_at_given_temperature(temperature);
        double Pbubble = fluid->get_bubble_point_at_given_temperature(temperature);

        if (liquid_volume / volume < 4.65e-06)//уровень в тесте ~0.05 мм
            return Pdew;
        if ((1. - liquid_volume / volume) < 6.52e-06)//уровень в тесте ~9899.93 мм при высоте 9.9м
            return Pbubble;

        fixed_solver_result_t<1> solver_result_carrier;
        if (solver_result == nullptr)
            solver_result = &solver_result_carrier;

        fixed_solver_parameters_t<1, 0, golden_section_search> solver_parameters;
        solver_parameters.constraints.minimum = Pdew;
        solver_parameters.constraints.maximum = Pbubble;
        solver_parameters.argument_increment_norm = 1e-7;

        double estimation = (Pdew + Pbubble) * liquid_volume / volume;// так поближе
        double dp = Pbubble - Pdew;
        estimation = std::max(Pdew + 1e-3 * dp, estimation);
        estimation = std::min(Pbubble - 1e-3 * dp, estimation);

        if (fabs(residuals(estimation)) < 1e-4)
            return estimation;

        fixed_newton_raphson<1>::solve_dense(
            *this, estimation, solver_parameters, solver_result);

        if (solver_result->result_code != numerical_result_code_t::Converged) {
            throw std::runtime_error("desired_liquid_vol_fraction_equation_fixed not converged");
        }
        return solver_result->argument;
    }
};


/// @brief Уравнение для поиска температуры по заданному значению функции флюида (фиксированная размерность)
/// @tparam Function Тип функции флюида, принимающей (const fluid_t*, double) и возвращающей double
template <typename Function>
struct vle_desired_fluid_function_equation_fixed : fixed_system_t<1> {
    /// @brief Флюид, для которого выполняется расчет
    const fluid_t* fluid;
    /// @brief Целевое значение функции
    double target_value;
    /// @brief Функция флюида
    Function function;
public:
    /// @brief Конструктор
    vle_desired_fluid_function_equation_fixed(
        const fluid_t* _fluid, double _target_value, Function _function)
        : fluid(_fluid)
        , target_value(_target_value)
        , function(_function)
    {}
    /// @brief Невязки уравнения
    virtual function_type residuals(const var_type& temperature) override {
        // здесь при расчете энтальпии фиктивно отдаем атмосферное давление, 
        // т.к. идеально-газовая энтальпия зависит только от температуры
        double function_value = function(fluid, temperature);
        double diff = function_value - target_value;
        return diff;
    }
    /// @brief Оптимизация расчета, учитывающая, что target_value не влияет на якобиан
    virtual function_type jacobian_dense(const var_type& temperature) override {
        double e = epsilon * std::max(1.0, abs(temperature));
        double function_value_0 = function(fluid, temperature);
        double function_value_eps = function(fluid, temperature + e);
        function_type J = (function_value_eps - function_value_0) / e;
        return J;
    }
};

/* НЕИСПОЛЬЗУЕМЫЙ КОД: Заменен на estimate_temperature_for_inner_energy_fixed (использует fixed_system_t)
/// @brief Оценка температуры по заданной внутренней энергии
/// @tparam AmountType Тип используемых единиц количества вещества (мольные, массовые)
/// @param fluid Флюид, для которого выполняется расчет
/// @param target_inner_energy Целевая удельная внутренняя энергия
/// @param initial_temperature Начальное приближение по температуре (опционально)
/// @param pressure Давление для расчета внутренней энергии (по умолчанию 1e5 Па)
/// @return Найденная температура или NaN, если расчет не сошелся
template <AmountType amount_type>
double estimate_temperature_for_inner_energy(
    const fluid_t* fluid,
    double target_inner_energy,
    double initial_temperature = std::numeric_limits<double>::quiet_NaN(),
    double pressure = 1e5)
{
    auto get_inner_energy = [](const fluid_t* fluid, double T) -> double {
        return fluid->get_ideal_gas_inner_energy<amount_type>(T);
    };
    
    vle_desired_fluid_function_equation eq(fluid, target_inner_energy, get_inner_energy);

    Eigen::VectorXd estimation(1);
    if (std::isfinite(initial_temperature)) {
        estimation(0) = initial_temperature;
    }
    else {
        // Задаем начальное приближение по температуре
        // Используем псевдокритическую температуру как начальное приближение
        double T_critical = fluid->get_pseudocritical_temperature();
        estimation(0) = T_critical;
    }

    num_methods::equation_solver_parameters_t solver_parameters;
    solver_parameters.constraints.minimum[0] = 100;
    solver_parameters.constraints.relative_boundaries[0] = 100;
    solver_parameters.argument_increment_norm = 1e-8;
    solver_parameters.argument_increment_relative = true;
    
    num_methods::equation_solver_result_t solver_result;
    num_methods::solve_gauss_newton_bounded(eq, estimation, solver_parameters, &solver_result);

    if (!solver_result.converged) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return solver_result.argument(0);
}
*/



/// @brief Оценка температуры по заданной внутренней энергии (версия с fixed_system_t)
/// @tparam AmountType Тип используемых единиц количества вещества (мольные, массовые)
/// @param fluid Флюид, для которого выполняется расчет
/// @param target_inner_energy Целевая удельная внутренняя энергия
/// @param initial_temperature Начальное приближение по температуре (опционально)
/// @param pressure Давление для расчета внутренней энергии (по умолчанию 1e5 Па)
/// @return Найденная температура или NaN, если расчет не сошелся
template <AmountType amount_type>
double estimate_temperature_for_inner_energy_fixed(
    const fluid_t* fluid,
    double target_inner_energy,
    double initial_temperature = std::numeric_limits<double>::quiet_NaN())
{
    auto get_inner_energy = [](const fluid_t* fluid, double T) -> double {
        return fluid->get_ideal_gas_inner_energy<amount_type>(T);
        };

    vle_desired_fluid_function_equation_fixed eq(fluid, target_inner_energy, get_inner_energy);

    // Определяем начальное приближение по температуре
    double T_initial;
    if (std::isfinite(initial_temperature)) {
        T_initial = initial_temperature;
    }
    else {
        // Используем псевдокритическую температуру как начальное приближение
        double T_critical = fluid->get_pseudocritical_temperature();
        T_initial = T_critical;
    }

    // Настройки решателя
    fixed_solver_parameters_t<1, 0, golden_section_search> solver_parameters;
    solver_parameters.line_search.function_decrement_factor = 10;
    solver_parameters.line_search.iteration_count = 30;
    solver_parameters.line_search_fail_action = line_search_fail_action_t::TreatAsFail;

    // Адаптированные настройки из estimate_temperature_for_inner_energy:
    // minimum[0] = 100 -> minimum = 100.0
    solver_parameters.constraints.minimum = get_min_antoine_bound(fluid);// 100.0;



    solver_parameters.constraints.maximum = std::numeric_limits<double>::quiet_NaN();

    // argument_increment_norm = 1e-8 (из оригинальной версии)
    solver_parameters.argument_increment_norm = 1e-8;

    // Настройки, которые не удается адаптировать к fixed_newton_raphson:
    // solver_parameters.constraints.relative_boundaries[0] = 100;  // специфично для gauss_newton_bounded
    // solver_parameters.argument_increment_relative = true;  // специфично для gauss_newton_bounded

    fixed_solver_result_t<1> solver_result;
    fixed_newton_raphson<1>::solve_dense(
        eq, T_initial, solver_parameters, &solver_result, nullptr);

    if (solver_result.result_code == numerical_result_code_t::Converged) {
        return solver_result.argument;
    }
    else {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

// НЕ ИСПОЛЬЗУЕТСЯ, НО ТУТ ЕСТЬ ЧТО-ТО ПРО МОНОКОМПОНЕНТ
//virtual double fluid_raoult_dalton_t::get_temperature_with_enthalpy_mass(double pressure, double enthalpy_mass) const override
//{
//    const auto& components = get_components();
//
//    auto last_flash_result = get_last_flash_result();
//    if (components.size() != 1) {
//
//        double initial_temperature = std::isfinite(last_flash_result.temperature)
//            ? last_flash_result.temperature
//            : KELVIN_OFFSET;
//        //double initial_temperature = KELVIN_OFFSET;
//        desired_enthalpy_fixed ph_flash(this, pressure, enthalpy_mass, initial_temperature);
//        double ret = ph_flash.solve();
//        if (!std::isfinite(ret)) {
//            crashdump_and_throw("PH-flash_gtwem" + std::to_string(clock()), ph_flash.get_stub_data());
//        }
//        return ret;
//    }
//
//    // Специальная логика для однокомпонентного потока
//    double boiling_temperature = get_dew_point_at_given_pressure(pressure);
//    double H1 = get_enthalpy_td_mass_liquid(pressure, boiling_temperature);
//    double H2 = get_enthalpy_td_mass_vapor(pressure, boiling_temperature);
//
//    if (enthalpy_mass >= H1 && enthalpy_mass < H2)
//        return boiling_temperature;
//
//    double initial_temperature = std::isfinite(last_flash_result.temperature)
//        ? last_flash_result.temperature
//        : boiling_temperature;
//    //double initial_temperature = boiling_temperature;
//    equation_solver_parameters_t parameters;
//    equation_solver_result_t solver_result;
//
//    if (enthalpy_mass < H1) {
//        // расчет для жидкости
//        auto liquid_enthalpy = [&](double temperature) {
//            return get_enthalpy_td_molar_liquid(pressure, temperature);
//            };
//
//        scalar_equation_wrapper_t eq(liquid_enthalpy, initial_temperature);
//        solve_gauss_newton_bounded(eq, eq.estimation(), parameters, &solver_result);
//        if (!solver_result.converged)
//            throw std::logic_error("Not converged");
//
//        return solver_result.argument(0);
//
//    }
//    else { // здесь обязательно будет (enthalpy_molar > H2)
//        // расчет для газа
//        auto vapor_enthalpy = [&](double temperature) {
//            return get_enthalpy_td_molar_vapor(pressure, temperature);
//            };
//
//        scalar_equation_wrapper_t eq(vapor_enthalpy, initial_temperature);
//        solve_gauss_newton_bounded(eq, eq.estimation(), parameters, &solver_result);
//        if (!solver_result.converged)
//            throw std::logic_error("Not converged");
//
//        return solver_result.argument(0);
//    }
//}

}
