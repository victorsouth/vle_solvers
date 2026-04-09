#pragma once
#include "fluid_base.h"
namespace vlelib {
;

struct rachford_rice_result_t;


/// @brief Расчеты парожидкостного равновесия
/// для идеального газа по объединенному закону Раулю-Дальтону и уравнению Антуана
class fluid_rault_dalton_t : public fluid_t {
protected:
public:
    /// \brief обнуляет last_flash_result в производных классах с кэшированием
    //virtual void invalidate_calculation()override{last_flash_result.invalidate_calculation();}

    /// @brief Копирующий конструктор с проверкой целостности мемоизации
    fluid_rault_dalton_t(const fluid_rault_dalton_t& other)
        : fluid_t(other)
    {
        /*if (other.last_flash_result.has_integrity) {
        last_flash_result = other.last_flash_result;
        }*/
    }
    /// @brief Конструктор на основе вектора компонентов
    fluid_rault_dalton_t(const std::vector<const component_properties_t*>& components)
        : fluid_t(components)
    {
    }
    /// @brief Конструктор на основе векторов компонентов (std::vector) и мольных долей (Eigen::VectorXd)
    fluid_rault_dalton_t(const std::vector<const component_properties_t*>& components, const Eigen::VectorXd& molar_fraction)
        : fluid_t(components, molar_fraction)
    {
    }
    /// @brief Конструктор на основе векторов (std::vector) компонентов и мольных долей
    fluid_rault_dalton_t(const std::vector<const component_properties_t*>& components, const std::vector<double>& molar_fraction)
        : fluid_t(components, molar_fraction)
    {
    }
    /// @brief Конструктор на основе векторов компонентов (std::vector), мольных 
    /// долей (Eigen::VectorXd) и матрицы бинарных коэффициентов (Eigen::MatrixXd)
    fluid_rault_dalton_t(const std::vector<const component_properties_t*>& components,
        const Eigen::VectorXd& molar_fraction,
        const Eigen::MatrixXd& binary_coeffs)
        : fluid_t(components, molar_fraction, binary_coeffs)
    {
    }

    /// @brief Копирует флюид вместе с данными мемоизации,
    /// вызывает копирующий конструктор, в котором проверяется целостность данных мемоизации
    virtual std::unique_ptr<fluid_t> create_copy(bool copy_memoization = true) const override
    {
        if (copy_memoization) {
            auto result = std::make_unique<fluid_rault_dalton_t>(*this);
            return std::move(result);
        }
        else {
            return create_copy(get_molar_fraction());
        }
    }
    /// @brief Копирует флюид без мемоизации, т.к. меняется состав, значит расчет некорректен
    virtual std::unique_ptr<fluid_t> create_copy(const Eigen::VectorXd& new_molar_fraction) const override
    {
        auto result = std::make_unique<fluid_rault_dalton_t>(
            get_components(), new_molar_fraction, get_binary_coeffs_ref());
        return std::move(result);
    }
    virtual std::unique_ptr<fluid_t> create_copy(const std::vector<double>& new_molar_fraction) const override
    {
        auto result = std::make_unique<fluid_rault_dalton_t>(
            get_components(),
            Eigen::VectorXd::Map(new_molar_fraction.data(), static_cast<Eigen::Index>(new_molar_fraction.size())),
            get_binary_coeffs_ref());
        return std::move(result);
    }
public:
    /// @brief Вектор плотностей газа. Идеально-газовый расчет
    virtual Eigen::VectorXd get_densities_vapor(double pressure, double temperature) const override
    {
        Eigen::VectorXd result(get_components_count());
        const auto& components = get_components();

        if (pressure != 0)
        {
            for (size_t index = 0; index < (size_t)result.size(); ++index) {
                result(index) = density_ideal_gas(pressure, temperature, components[index]->molar_mass);
            }
        }
        else //pressure == 0
        {
            for (size_t index = 0; index < (size_t)result.size(); ++index)
                result(index) = 0;
        }
        return result;
    }

    /// @brief Расчет плотности чистого вещества в жидкой фазе
    /// Можно сделать отдельной функцией, но статический метод
    /// отражает применный в fluid_rault_dalton_t способ расчета
    /// @param component Параметры чистого вещества
    static double get_density_liquid(const component_properties_t& component,
                                     double pressure, double temperature);
    /// @brief Возвращает вектор плотностей чистых компонентов в жидком состоянии
    /// Учитывается поправка на давление (коэф. сжимаемости) и температуру (ф-ла Мановяна)
    virtual Eigen::VectorXd get_densities_liquid(double pressure, double temperature) const override;
    /// @brief Возвращает объемные доли компонентов до flash-расчета
    /// Используется при заполнении объема по количествую вещества и энтальпии
    /// Описание - см. Таблицу пересчета количеств
    Eigen::VectorXd get_liquid_volumetric_fracs(double pressure, double temperature) const;
    /// @brief Возвращает плотность смеси в предположении жидкого состояния
    /// (именно смеси, не отдельной жидкой фазы!)
    /// @param pressure Учитывается через коэффициент сжимаемости
    /// @param temperature Учитывается через формулу Мановяна
    double get_density_as_liquid(double pressure, double temperature, double* molar_volume_ptr = nullptr) const;

    /// @brief Расчет частной производной плотности смеси по давлению
    /// (в предположении жидкого состояния)
    double get_density_as_liquid_pressure_derivative(double pressure, double temperature) const;

protected:
    /// @brief Получает результат расчёта флюида, если он состоит только из жидкости
    const flash_calculation_result_t flash_liquid_only(double pressure, double temperature,
                                                       flash_calculation_result_t& result
                                                       ) const;

    /// @brief Получает результат расчёта флюида, если он состоит только из газа
    const flash_calculation_result_t flash_vapor_only(double pressure, double temperature,
                                                      flash_calculation_result_t& result) const;
public:
    /// @brief Генерирует результат расчета в двухфазном случае (бывшая экспериментальная версия)
    /// @param pressure Давление
    /// @param temperature Температура
    /// @param rr_result Корректный результа расчета Речфорда-Райса
    /// @return result Результат расчета (заполненный объект структуры flash_calculation_result_t)
    void build_twophase_result2(double pressure, double temperature, const rachford_rice_result_t& rr_result,
                                flash_calculation_result_t& result) const;

    /// @brief Функция от парожидкостного равновесия. Публичная для работы two_phase_result_builder
    template <AmountType amount_type, typename GasFunction, typename LiquidFunction>
    amounts_per_phase flash_function(GasFunction gasf, LiquidFunction liquidf,const flash_calculation_result_t& flash_result) const;

public:

    virtual void flash_unsafe(
            double pressure, double temperature,
            double initial_estimation,
            flash_calculation_result_t& result
            ) const override;

    //virtual double get_density(double pressure, double temperature) const override
    //{
    //    const flash_calculation_result_t flash_result = flash(pressure, temperature);
    //    if (!std::isfinite(flash_result.density.mix))
    //        throw logic_error("flash_result.density.mix is nan");
    //    return flash_result.density.mix;
    //}

    /// @brief Точка росы при данной температуре (документ PT-flash, раздел 4.5)
    virtual double get_dew_point_at_given_temperature(double temperature) const override;
    /// @brief Точка росы (начала конденсации) при данном давлении
    virtual double get_dew_point_at_given_pressure(double pressure) const override;
    /// \brief Точка росы по воде при данном давлении
    virtual double get_water_dew_point_at_given_pressure(double pressure) const override;
    /// @brief Давление начала кипения при данной температуре (документ PT-flash, раздел 4.4)
    virtual double get_bubble_point_at_given_temperature(double temperature) const override;
    /// @brief Температура начала кипения при данном давлении
    virtual double get_bubble_point_at_given_pressure(double pressure) const override;

    /// @brief Изменить объем жидкости за счет газа, сохранив при этом составы жидкости и газа
    /// Если жидкости нет, то выдать ошибку
    /// @param pressure
    /// @param temperature
    /// @param liquid_volume_fraction
    virtual void change_liquid_volume_fraction(
            double pressure, double temperature, double liquid_volume_fraction) override;
    /// @brief Проверить жидкое накопление и если оно есть, то посчитать давление
    /// @param volume
    /// @param temperature
    /// @param molar_amount
    /// @return Давление для режима жидкого накопления, NaN если такого режима нет
    double fill_and_check_liquid_accumulation(
            double volume, double temperature, double molar_amount) const;
    /// @brief Заполнить заданный объем веществом с заданной температурой так,
    /// чтобы объем жидкости был равен заданному
    /// @param volume
    /// @param temperature
    /// @param liquid_volume
    /// @return Давление в заданном объеме
    virtual double fill_volume_with_liquid(double total_volume, double temperature, double liquid_volume) const override;
    /// @brief Заполнить заданный объем веществом с заданной температурой так,
    /// чтобы давление было равно заданному
    /// @return Количество вещества
    virtual double fill_volume_with_pressure(double pressure, double temperature,
                                             double total_volume) const override;

    /// @brief Заполнить заданный объем заданным количеством вещества данного состава с заданной температурой
    /// @return давление и температура в заданном объеме
    virtual std::pair<double, double> fill_volume_with_total_moles2(
            double volume, double inner_energy_molar, double molar_amount,
            double initial_pressure = std::numeric_limits<double>::quiet_NaN(),
            double initial_temperature = std::numeric_limits<double>::quiet_NaN()) const override;

    /// @brief Рассчитывает температуру по внутренней энергии в предположении жидкого состояния
    /// @tparam amount_type Задает тип внутренней энергии
    /// @param inner_energy Удельная внутренняя энергия (мольная или массовая)
    /// @return Температура
    template <AmountType amount_type>
    double find_liquid_temperature_with_inner_energy(double inner_energy) const;

};   // end class fluid_rault_dalton_t


/// @brief Промежуточные вычисления, требуемые для заполнения flash_calculation_result_t
/// Это в основном векторные параметры
struct intermediate_twophase_data {
    /// @brief Плотность всех компонентов в парообразном состояний
    Eigen::VectorXd densities_vapor;
    /// @brief Плотность всех компонентов в жидкофазном состояний
    Eigen::VectorXd densities_liquid;
    /// @brief Молярная масса всех компонентов 
    Eigen::VectorXd M;
    /// @brief Произведения yi*Mi по всем компонентам
    Eigen::VectorXd yiMi;
    /// @brief Произведения xi*Mi по всем компонентам
    Eigen::VectorXd xiMi;
    /// @brief Объемы чистых компонентов по газу на один моль газа
    Eigen::VectorXd w_vap;
    /// @brief Объемы чистых компонентов по жидкости на один моль газа
    Eigen::VectorXd w_liq;
    /// @brief Конструктор, внутри считает все требуемые промежуточные параметры
    intermediate_twophase_data(
        const rachford_rice_result_t& rr_result,
        const Eigen::VectorXd& densities_vapor,
        const Eigen::VectorXd& densities_liquid, const Eigen::VectorXd& M);
};

/// @brief Производные расчеты для флюида, выполняемые после расчета парожидкостного равновесия.
/// По сути - заполняет flash_calculation_result_t
/// Расчет равновесия подразумевался по Раулю-Дальтону
class two_phase_result_builder {
public:
    /// @brief Генерирует свойства смеси в двухфазном случае
    /// @param[in] pressure Давление
    /// @param[in] temperature Температура
    /// @param[in] rr_result Результа расчета Речфорда-Райса
    /// @param[in] flash_in Объект fluid_rault_dalton_t с составом смеси
    /// @param[out] result Заполненный объект структуры flash_calculation_result_t
    void build(double pressure, double temperature, const rachford_rice_result_t& rr_result
        , const fluid_rault_dalton_t& flash_in
        , flash_calculation_result_t& result) const;
private:
    /// @brief Промежуточная процедура, записывающая критические параметры
    /// @param[in] flash_in Объект fluid_rault_dalton_t с составом смеси
    /// @param[out] result Заполняемый объект структуры flash_calculation_result_t
    void set_critical_parameters(const fluid_rault_dalton_t& flash_in,
        flash_calculation_result_t& result) const;

    /// @brief Промежуточная функция, возвращающая молярные массы фаз смеси и смеси
    /// @result Объект структуры amounts_per_phase с молярными массами фаз и смеси
    /// @param[in] intermediate Объект intermediate_twophase_data с результатами промежуточных расчётов
    /// @param[in] flash_in Объект fluid_rault_dalton_t с составом смеси
    amounts_per_phase get_molar_masses(const intermediate_twophase_data& intermediate
        , const fluid_rault_dalton_t& flash_in) const;

    /// @brief Промежуточная функция, возвращающая кортеж из массовой доли отгона и объемной доли отгона
    /// @result Кортеж из массовой доли отгона и объемной доли отгона
    /// @param[in] omega Мольная доля отгона
    /// @param[in] molar_mass Молярные массы фаз и смеси
    /// @param[in] molar_volume Молярные объёмы фаз и смеси
    std::pair<double, double> get_vapor_fractions(
        double omega, const amounts_per_phase& molar_mass, const amounts_per_phase& molar_volume) const;

    /// @brief Промежуточная функция, возвращающая молярные объёмы фаз и смеси
    /// @result Объект структуры amounts_per_phase с молярными объёмами фаз и смеси
    /// @param[in] intermediate Объект intermediate_twophase_data с результатами промежуточных расчётов
    /// @param[in] flash_in Объект fluid_rault_dalton_t с составом смеси
    amounts_per_phase get_molar_volumes(double omega
        , const intermediate_twophase_data& intermediate
        , const fluid_rault_dalton_t& flash_in) const;

    /// @brief Промежуточная функция, возвращающая плотности фаз и смеси
    /// @result Объект структуры amounts_per_phase с плотностями фаз и смеси
    /// @param[in] molar_mass Молярные массы фаз и смеси
    /// @param[in] molar_volume Молярные объёмы фаз и смеси
    amounts_per_phase get_densities(const amounts_per_phase& molar_mass
        , const amounts_per_phase& molar_volume) const;

    /// @brief Промежуточная функция, возвращающая энтальпии фаз и смеси
    /// @result Объект структуры amounts_per_phase с энтальпиями фаз и смеси
    /// @param[in] flash_result Заполняемый объект структуры flash_calculation_result_t
    /// @param[in] flash_in Объект fluid_rault_dalton_t с составом смеси
    amounts_molar_and_mass get_enthalpies(const flash_calculation_result_t& flash_result
        , const fluid_rault_dalton_t& flash_in) const;

    /// @brief Промежуточная функция, возвращающая внутренние энергии фаз и смеси
    /// @result Объект структуры amounts_per_phase с внутренними энергиями фаз и смеси
    /// @param[in] flash_result Заполняемый объект структуры flash_calculation_result_t
    /// @param[in] flash_in Объект fluid_rault_dalton_t с составом смеси
    amounts_molar_and_mass get_inner_energies(const flash_calculation_result_t& flash_result
        , const fluid_rault_dalton_t& flash_in) const;

};

/// @brief Возвращает среднюю по составу минимальную температурную границу
/// области определения модели давления насыщенных паров Антуана
/// Средняя берется по коцентрациям компонентов в составе
/// @param fluid Состав
/// @return Минимальная граница T_min
double get_min_antoine_bound(const fluid_t* fluid);

/// @brief Возвращает среднюю по составу максимальную температурную границу
/// области определения модели давления насыщенных паров Антуана
/// Средняя берется по коцентрациям компонентов в составе
/// @param fluid Состав
/// @return Максимальная граница T_max
double get_max_antoine_bound(const fluid_rault_dalton_t* fluid);

template <AmountType amount_type>
extern double get_energy_gas(const component_properties_t* c, double pressure, double temperature);

template <AmountType amount_type>
double get_energy_liq(const component_properties_t* c, double pressure, double temperature);

template <AmountType amount_type>
double get_enthalpy_gas(const component_properties_t* c, double pressure, double temperature);

template <AmountType amount_type>
double get_enthalpy_liq(const component_properties_t* c, double pressure, double temperature);


template <AmountType amount_type, typename GasFunction, typename LiquidFunction>
amounts_per_phase vlelib::fluid_rault_dalton_t::flash_function(GasFunction gasf, LiquidFunction liquidf,const flash_calculation_result_t& flash_result) const
{
    //const flash_calculation_result_t& flash_result = last_flash_result;

    amounts_per_phase func_result;

    if (flash_result.fluid_vapor != nullptr) {
        Eigen::VectorXd gas_values = get_function_by_components(flash_result.pressure, flash_result.temperature, gasf);
        Eigen::VectorXd frac_vapor = amount_type == AmountType::Mass
            ? flash_result.fluid_vapor->get_mass_fraction()
                              : flash_result.fluid_vapor->get_molar_fraction();
        func_result.vapor = frac_vapor.dot(gas_values);
    }
    else {
        func_result.vapor = 0;
    }

    if (flash_result.fluid_liquid != nullptr) {
        Eigen::VectorXd liquid_values = get_function_by_components(flash_result.pressure, flash_result.temperature, liquidf);
        Eigen::VectorXd frac_liquid = amount_type == AmountType::Mass
                               ? flash_result.fluid_liquid->get_mass_fraction()
                               : flash_result.fluid_liquid->get_molar_fraction();

        func_result.liquid = frac_liquid.dot(liquid_values);
    }
    else {
        func_result.liquid = 0;
    }

    double vapor_frac = amount_type == AmountType::Mass
                        ? flash_result.vapor_mass_fraction
                        : flash_result.flash;
    func_result.mix =
            vapor_frac * func_result.vapor
            + (1 - vapor_frac) * func_result.liquid;

    return func_result;
}

}   // end namespace vlelib
