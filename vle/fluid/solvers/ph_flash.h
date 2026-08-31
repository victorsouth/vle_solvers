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
struct ph_flash_bisection : public fixed_system_t<1> {
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
    virtual function_type residuals(const var_type& x) override;
    /// \brief
    /// конструктор ph_flash_bisection
    /// \param f исследуемый флюид
    /// \param de энтальпия для которой необходимо найти температуру флюида
    /// \param p заданное давление флюида
    /// \param t начальное приближение по температуре, в общем случае может быть не задано
    ph_flash_bisection(fluid_t* f, double de, double p, double t = std::numeric_limits<double>::quiet_NaN());
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
    );
    /// \brief
    /// решить PH flash методом бисекции
    /// \param numerical_result - результат
    /// \param analysis - анализ
    /// \param desired_precision - желаемая точность по невязке
    /// \return
    double solve(
        fixed_bisection_result_t<1>* numerical_result = nullptr,
        fixed_bisection_result_analysis_t<1>* analysis = nullptr
    );
    /// @brief Возвращает данные для изолированного вызова
    ph_flash_stub_data_t get_stub_data() const;
};

/// @brief Уравнение для поиска температуры по заданной энтальпии 
/// с учетом фазовых переходов
struct ph_flash_newton : fixed_system_t<1>
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
    ph_flash_newton(const fluid_t* fluid, double pressure, double target_enthalpy_mass,
        double _temperature_initial = std::numeric_limits<double>::quiet_NaN(),
        double _prec = std::numeric_limits<double>::epsilon() * 10000.);
    /// @brief Энтальпия расчитывается по начальной, подводимому теплу и вносу энтальпии расходом 
    ph_flash_newton(const fluid_t* fluid, double p,
        double initial_entalpy, double Q, double mass_flow,
        double temperature_initial);
    /// @brief Невязки уравнения
    virtual double residuals(const double& temperature) override;
    /// @brief Оптимизация расчета, учитывающая, что jacobian запускается после residuals
    virtual double jacobian_dense(const double& temperature) override;
    /// @brief Возвращает начальное приближение по температуре
    var_type estimation() const;
    /// @brief Расчет PH-flash
    /// @param numerical_result 
    /// @param analysis_result 
    /// @return Искомая температура или NaN, если расчет не сошелся
    double solve(fixed_solver_result_t<1>* numerical_result = nullptr,
        fixed_solver_result_analysis_t<1>* analysis_result = nullptr);
    /// @brief Возвращает данные для изолированного вызова
    ph_flash_stub_data_t get_stub_data() const;
};


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
    desired_liquid_vol_fraction_equation_fixed(double volume, double temperature, double liquid_volume, const fluid_t* fluid);
    /// @brief Невязки уравнения
    virtual double residuals(const double& pressure) override;
    /// @brief Костылим расчет производной.
    /// В областях, близких однофазным производная обнуляется, что не дает считать численный метод
    virtual double jacobian_dense(const double& pressure) override;
    /// @brief Решение уравнения методом Ньютона-Рафсона
    /// @param solver_result Результат работы численного метода (опционально)
    /// @return Давление, при котором достигается заданная объемная доля жидкости
    double solve(fixed_solver_result_t<1>* solver_result = nullptr);
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
    double initial_temperature = std::numeric_limits<double>::quiet_NaN());

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
//        ph_flash_newton ph_flash(this, pressure, enthalpy_mass, initial_temperature);
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
