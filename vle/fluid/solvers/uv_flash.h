#pragma once

namespace vlelib {
;




/// @brief вспомогательная структура для системы уравнений целевой плотности и целевой внутренней энергии
/// сохраняется в случае неудачного завершения расчёта системы уравнений целевой плотности и целевой внутренней энергии
/// служит для загрузки и воспроизведения упавшего расчёта
struct uv_flash_over_RR_stubdata_t {
    double specific_volume_molar;///< Целевой мольный объём
    double inner_energy_molar;///< Целевая удельная мольная внутренняя энергия
    fluid_stubdata_t fluid;///< Состав, для которого выполняется расчёт
    bool use_density;///< Использовать плотность в уравнении заполнения
    /// @brief Начальное приближение по давлению, используется в расчете в sovle
    /// Рассчитывается в конструкторе с учетом возможно заданного извне нач. приближения
    /// и обеспечения двухфазной области
    double initial_pressure = std::numeric_limits<double>::quiet_NaN();
    /// @brief Начальное приближение по давлению, используется в расчете в sovle
    /// @brief Рассчитывается в конструкторе с учетом T_min и возможно заданного извне нач. приближения
    /// Рассчитывается в конструкторе с учетом возможного известного начального приближений
    /// и обеспечения двухфазной области
    double initial_temperature = std::numeric_limits<double>::quiet_NaN();
#ifdef VLELIB_SERIALIZATION_SUPPORT
    /// @brief Для сериализации
    friend class boost::serialization::access;

    /// @brief Запись в заданный файл с помощью сериализатора нужного формата
    template <typename Serializer = boost::archive::xml_oarchive>
    void to_file(const std::string& filename = "uv_flash_failed.xml") const {
        std::ofstream ofs(filename);
        if (ofs.is_open()) {
            boost::archive::xml_oarchive oa(ofs);
            const auto& mock_data = *this;
            oa << BOOST_SERIALIZATION_NVP(mock_data);
        }
        else {
            throw std::runtime_error("Failed to write UV-flash mockdata");
        }
        ofs.flush();
        ofs.close();
    }
    /// @brief Интрузивная сериализация/десериализация
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_NVP(specific_volume_molar);
        ar& BOOST_SERIALIZATION_NVP(inner_energy_molar);
        ar& BOOST_SERIALIZATION_NVP(use_density);
        ar& BOOST_SERIALIZATION_NVP(fluid);
    }
#endif
};

/// @brief Используется для теста сходимости. 
/// Включает в себя обычную стабдату, а также эталонные значения P, T, которые известны только в тестах
struct uv_flash_over_RR_convergence_test_stubdata_t
{
    /// @brief Параметры алгоритма UV-flash
    uv_flash_over_RR_stubdata_t stubdata;
    /// @brief Эталонное даление
    double etalon_pressure;
    /// @brief Эталонная температура
    double etalon_temperature;
    /// @brief Эталонная доля отгона
    double etalon_vapor_volumetric_fraction;

#ifdef VLELIB_SERIALIZATION_SUPPORT
    /// @brief Для сериализации
    friend class boost::serialization::access;
    /// @brief Интрузивная сериализация/десериализация
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        stubdata.serialize<Archive>(ar, version);
        ar& BOOST_SERIALIZATION_NVP(etalon_pressure);
        ar& BOOST_SERIALIZATION_NVP(etalon_temperature);
        ar& BOOST_SERIALIZATION_NVP(etalon_vapor_volumetric_fraction);
    }
#endif
};


/// @brief Опции алгоритма UV-flash для попыток расчета
struct uv_flash_solution_attempt_options {
    /// @brief начальная оценка
    std::array<double, 2> initial_estimation;
    /// @brief использовать quadprog
    bool use_quadprog;
    /// \brief конструктор
    uv_flash_solution_attempt_options(const std::array<double, 2>& initial_estimation, bool use_quadprog)
        : initial_estimation(initial_estimation)
        , use_quadprog(use_quadprog)
    {
    }
};


/// @brief Система уравнений целевой плотности и целевой внутренней энергии
/// Уравнение заполнения
/// Уравнение целевой энергии
/// Решается относительно нормированных P, T, вложенный Речфорд-Райс
/// Используется нормировка давления и температуры (см. "Эффективный алгоритм UV-flash")
class uv_flash_over_RR : public fixed_system_t<2> {
private:
    /// @brief Состав, для которого выполняется расчет
    const fluid_rault_dalton_t* fluid;
    /// @brief Целевой мольный объем 
    double specific_volume_molar;
    /// @brief Целевая удельная мольная внутренняя энергия
    double inner_energy_molar;
    /// @brief Начальное приближение по давлению, используется в расчете в sovle
    /// Рассчитывается в конструкторе с учетом возможно заданного извне нач. приближения
    /// и обеспечения двухфазной области
    double initial_pressure;
    /// @brief Начальное приближение по давлению, используется в расчете в sovle
    /// @brief Рассчитывается в конструкторе с учетом T_min и возможно заданного извне нач. приближения
    /// Рассчитывается в конструкторе с учетом возможного известного начального приближений 
    /// и обеспечения двухфазной области
    double initial_temperature;
    /// @brief Псевдокритическая температура (считается в конструкторе по составу)
    double T_critical;
    /// @brief Минимальная температура поиска
    double T_min;
    /// @brief Использовать плотность в уравнении заполнения
    bool use_density;
private:
    using fixed_system_t<2>::var_type;
    using fixed_system_t<2>::function_type;
    using fixed_system_t<2>::matrix_value;
protected:
    /// @brief При необходимости корректирует начальное приближение по давлению, 
    /// обеспечивая значения из двухфазной области
    /// @return Откорректированное давление, принадлежащее двухфазной области
    static double ensure_two_phase_pressure(const fluid_rault_dalton_t* fluid, double pressure, double temperature);

public:
    /// @brief Нормировка параметров P,T в параметры \phi, \theta
    /// Аналогично нестатическому norm
    /// @return (Нормированное давление (\phi), Приведенная температура (\theta))
    static std::pair<double, double> norm(const fluid_rault_dalton_t* fluid, double P, double T);
    /// @brief Нормировка параметров P, T в параметры (см. "Эффективный алгоритм UV-flash"):
    /// - Нормированное давление (\phi), Приведенная температура (\theta)
    /// @return (Нормированное давление (\phi), Приведенная температура (\theta))
    std::pair<double, double> norm(double P, double T) const;
    /// @brief По нормированным значениям давления, температуры возвращает человеческие P, T в СИ
    /// Статический метод для внешнего использования без создания экземпляра класса
    /// Аналогично нестатическому denorm
    /// @return (Давление, Температура)
    static std::pair<double, double> denorm(const fluid_rault_dalton_t* fluid, double phi, double theta);
    /// @brief По нормированным значениям давления, температуры возвращает  P, T в СИ (см. "Эффективный алгоритм UV-flash")
    /// @param phi Нормированное давление (\phi)
    /// @param theta Нормированная температура (\theta)
    /// @return (Давление, Температура)
    std::pair<double, double> denorm(double phi, double theta) const;
public:
    /// @brief Минимальное значение приведенной температуры
    double get_min_theta() const;
public:
    /// @brief Конструктор.
    /// Если входные initial_pressure initial_temperature не заданы, 
    /// то задаются значения, соответствющие {0.5, 0.5} 
    uv_flash_over_RR(double volume, double molar_amount, double _inner_energy_molar, 
        const fluid_rault_dalton_t* _fluid, bool _use_density, 
        double _initial_pressure = std::numeric_limits<double>::quiet_NaN(), 
        double _initial_temperature = std::numeric_limits<double>::quiet_NaN());
    /// @brief Невязки по двум уравнениям (заполнение, энергия)
    /// Уравнение заполнение считается по-разному в зависимости от опции use_density
    var_type residuals(const var_type& w) override;
    /// @brief Якобиан системы уравнений, учитывается фикс 
    /// для расчета произодной по нормирвоанного давлению на газовой и жидкой границе 
    matrix_value jacobian_dense(const var_type& x) override;
    /// @brief Специфический критерий успешного завершения расчета
    virtual bool custom_success_criteria(const var_type& r, const var_type& x);
    /// @brief Тривиально выдает значения initial_pressure, initial_temperature, рассчитанные в конструкторе,
    /// лишь трививально пронормировав их
    var_type estimation() const;
private:
    /// @brief Подготовка параметров солвера Ньютона для задачи UV-flash
    fixed_solver_parameters_t<2, 0, golden_section_search> prepare_solver_parameters(
        bool perform_analysis) const;
    fixed_solver_result_t<2> run_solution_attemps(
        const std::array <uv_flash_solution_attempt_options, 4>& options,
        const fixed_solver_parameters_t<2, 0, golden_section_search>& solver_parameters,
        fixed_solver_result_analysis_t<2>* analysis_result);
public:
    /// @brief Расчет задачи UV-flash
    /// @param analysis_result Аналитика сходимости, только для отладки расчета
    /// В боевом режиме нужно этот указатель задать nullptr, тогда затратный расчет аналитики не проводится
    /// @return Результат расчета системы уравнения численным методом
    fixed_solver_result_t<2> solve(fixed_solver_result_analysis_t<2>* analysis_result = nullptr);

    /// @brief Параметры изолированного вызова (Mock)
    /// @brief Минимально необходимый набор свойств для проведения расчета uv_flash
    const uv_flash_over_RR_stubdata_t get_mock_data() const;
};


}
