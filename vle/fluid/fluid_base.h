#pragma once

namespace vlelib {
;
/// @brief Состояние флюида и flash-расчета
struct fluid_state_t {
    /// @brief Концентрации
    std::vector<double> concentration;
    /// @brief Доля отгона последнего расчета
    double flash;
    /// @brief Давление последнего расчета
    double pressure;
    /// @brief Температура последнего расчета
    double temperature;
    /// @brief чтение данных класса из входного потока
    /// производится проверка корректности входных данных по имени
    /// @param stream входной поток
    virtual void deserialize_text(std::istream& stream){
        std::string str;
        stream>>str;
        if(str!=vle_solvers::get_class_as_string(*this))throw std::logic_error("wrong type:"+str);
        fixed_solvers::load_vector(stream,concentration);
        stream>>flash>>pressure>>temperature;
        /// вызываем исключение в случае неудачного чтения из потока
        if(stream.fail()){
            throw std::runtime_error("wrong data");
        }

    }
    /// @brief запись данных класса в выходной поток
    /// @param stream выходной поток
    virtual void serialize_text(std::ostream& stream)const{
        stream<< vle_solvers::get_class_as_string(*this)<<std::endl;
        fixed_solvers::save_vector(stream,concentration);
        stream<<flash<<' '<<pressure<<' '<<temperature<<'\n';
    }
    virtual ~fluid_state_t() = default;
};

/// @brief Переопределение оператора << для вывода данных состояния флюида в строку
/// @param[in] os Входная строка для дополнения
/// @param[in] s Объект с состоянием флюида для записи в строку
/// @return Выходная строка с данными состояния флюида
inline std::ostream& operator<<(std::ostream&os,const fluid_state_t& s){
    s.serialize_text(os);
    return os;
}

/// @brief Переопределение оператора >> для получения данных состояния флюида из строки
/// @param[in] is Входная строка
/// @param[out] s Объект с состоянием флюида, прочитанным из строки
/// @return Возвращает входную строку
inline std::istream& operator>>(std::istream&is,fluid_state_t& s){
    s.deserialize_text(is);
    return is;
}


/// @brief Данные для изолированного вызова задачи PT-flash
struct pt_flash_stub_data_t
{
    /// @brief Давление задачи PH-flash
    double pressure;
    /// @brief Температура
    double temperature;
    /// @brief Данные для изолированного вызова по смеси
    fluid_stubdata_t fluid;

#ifdef VLELIB_SERIALIZATION_SUPPORT
    /// @brief Для сериализации
    friend class boost::serialization::access;

    /// @brief Запись в заданный файл с помощью сериализатора нужного формата
    template <typename Serializer = boost::archive::xml_oarchive>
    void to_file(const std::string& filename = "pе_flash_failed.xml") const {
        std::ofstream ofs(filename);
        if (ofs.is_open()) {
            Serializer oa(ofs);
            const auto& mock_data = *this;
            oa << BOOST_SERIALIZATION_NVP(mock_data);
        }
        else {
            throw std::runtime_error("Failed to write PT-flash mockdata");
        }
        ofs.flush();
        ofs.close();
    }
    /// @brief Интрузивная сериализация/десериализация
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_NVP(pressure);
        ar& BOOST_SERIALIZATION_NVP(temperature);
        ar& BOOST_SERIALIZATION_NVP(fluid);
    }
    static pt_flash_stub_data_t from_file(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            throw std::runtime_error("Cannot open file");
        }
        boost::archive::xml_iarchive oa(ifs);
        pt_flash_stub_data_t PT_data;
        oa >> BOOST_SERIALIZATION_NVP(PT_data);
        return PT_data;
    }
#endif
};


/// @brief Функция генерирует исключение
/// @tparam StubData
/// @param crashed_algorithm Алгоритм, на котором упал расчет
/// @param stub_data Стабдата вылетевшего расчета
template <typename StubData>
inline void crashdump_and_throw(
        const std::string& crashed_algorithm,
        const StubData& stub_data)
{
#ifdef VLELIB_SERIALIZATION_SUPPORT
    try {
        stub_data.template to_file<boost::archive::xml_oarchive>(
                    crashed_algorithm + "_failed.xml");
    }
    catch (std::exception&) {
        throw std::runtime_error(
                    crashed_algorithm + " not converged, failed to write XML data");
    }
    throw std::runtime_error(crashed_algorithm + " not converged, XML data written");
#else
    throw std::runtime_error(crashed_algorithm + " not converged, XML dump disabled");
#endif
}

/// @brief Фазовое состояние вещества
enum class state_of_matter_t { Undefined, Gas, Liquid, TwoPhase, Critical };

/// @brief Результаты flash-расчета
struct flash_calculation_result_t {
    /// @brief Целостность данных
    bool has_integrity{ false };
    /// @brief Паровая фаза
    std::shared_ptr<fluid_t> fluid_vapor;
    /// @brief Жидкая фаза
    std::shared_ptr<fluid_t> fluid_liquid;
    /// @brief Давление последнего расчета (мемоизация)
    double pressure{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Температура последнего расчета (мемоизация)
    double temperature{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Молярная масса смеси, молярная масса газовой фазы (на 1 моль газа), молярная масса жидкой фазы (на один моль жидкости)
    amounts_per_phase molar_mass;
    /// @brief молярная доля отгона
    double flash{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief объемная доля отгона
    double vapor_volumetric_fraction{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief массовая доля отгона
    double vapor_mass_fraction{ std::numeric_limits<double>::quiet_NaN() };

    /// @brief Фазовое состояние вещества
    state_of_matter_t state_of_matter{ state_of_matter_t::Undefined };
    /// @brief Критическая или псевдокритическая температура
    double critical_temperature{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Критическая или псевдокритическая температура
    double critical_pressure{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Давление насыщенных паров жидкости (если жидкая фаза есть)
    double saturated_pressure_liquid{ std::numeric_limits<double>::quiet_NaN() };

    /// @brief Плотность смеси, пара, жидкости
    amounts_per_phase density;
    /// @brief Молярный объем пара, жидкости, смеси (на 1 моль пара, жидкости, смеси соответственно)
    amounts_per_phase molar_volume;
    /// @brief Фактор сжимаемости Z пара, жидкости, смеси (корни EOS PR; для идеального газа vapor = 1).
    amounts_per_phase z_factor;
    /// @brief Коэффициенты равновесия по PR: K_i = phi_l,i / phi_v,i (фугитивности смеси при том же
    ///        мольном составе, что и у flash, Z жидкости и Z пара — минимальный и максимальный корни куба).
    std::vector<double> k_value;
    /// @brief Энтальпия пара, жидкости, смеси (по массе и по молям)
    amounts_molar_and_mass enthalpy;
    /// @brief Внутренняя энергия для газа, жидкости (по массе и по молям)
    amounts_molar_and_mass inner_energy;
    /// @brief Смесевой объёмный сдвиг жидкости sum_i x_i dv_i (м^3/моль): тот же состав x,
    ///        что при расчёте molar_volume.liquid в PR flash
    double liquid_volume_shift_mix{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Возвращает тип состояния флюида
    flash_type_t get_flash_status() const;
    /// @brief Состоит ли флюид только из жидкости
    bool is_liquid_only() const;
    /// @brief Состоит ли флюид только из газа
    bool is_gas_only() const;
    /// @brief Является ли флюид двухфазным
    bool is_two_phase() const;

    /// @brief Возвращает мольную долю остатка
    double get_liquid_molar_fraction() const;
    /// @brief Возвращает объемную долю остатка
    double get_liquid_volume_fraction() const;
    /// @brief Возвращает массовую долю остатка
    double get_liquid_mass_fraction() const;

    /// @brief Проверяет несколько условий для принятия решения, что результат сформирован (является валидным)
    /// @param _pressure Давление должно совпасть с давлением флюида
    /// @param _temperature Температура должно совпасть с температурой флюида
    /// @result Также проверяется has_integrity и конечность давления и температуры флюида
    bool was_calculated(double _pressure, double _temperature) const;
    /// @brief Состояние флюида делается невалидным
    void invalidate_calculation();
};

/// @brief Класс для хранения и потокобезопасного доступа к данных флюида и результатам flash-расчета
///
/// Класс намеренно объединяет несколько зон ответственности:
/// - потокобезопасное хранение и доступ к данным флюида
/// - интерфейс для выполнения flash-расчетов
/// - кэширование результатов расчетов с потокобезопасным доступом
class fluid_fundamental_data_t {
private:
    /// @brief Параметры компонентов флюида. Ссылка на БД
    /// (никогда не меняется, запрещаем на уровне интерфейса)
    const std::vector<const component_properties_t*> components_;
    /// @brief Коэффициенты бинарного взаимодействия для компонентов данного флюида
    /// (никогда не меняется, запрещаем на уровне интерфейса)
    std::shared_ptr<const Eigen::MatrixXd> binary_coeffs_;
private:
    /// @brief Мольный состав
    Eigen::VectorXd concentration_;
    /// @brief Общий мьютекс для состава
    mutable std::recursive_mutex concentration_mutex;
private:
    /// @brief Мемоизация (кэш) последнего flash-расчета
    mutable flash_calculation_result_t last_flash_result_;
    /// @brief Мьютекс на кэш последнего flash-расчета
    mutable std::recursive_mutex flash_and_cache_mutex;
public:
    /// @brief Конструктор копирования.
    /// Копирует состав и список компонентов
    /// НЕ копирует мемоизацию и мьютексы
    explicit fluid_fundamental_data_t(const fluid_fundamental_data_t& other);
    /// @brief Инициализация по перечню компонентов.
    /// Состав инициализируется поровну между компонентами
    fluid_fundamental_data_t(const std::vector<const component_properties_t*>& components);
    /// @brief Инициализация по переченю компонентов и их концентрации (Eigen::VectorXd)
    fluid_fundamental_data_t(const std::vector<const component_properties_t*>& components,
                             const Eigen::VectorXd& components_concentration);
    /// @brief Инициализация по переченю компонентов и их концентрации (std::vector)
    fluid_fundamental_data_t(const std::vector<const component_properties_t*>& components,
                             const std::vector<double>& components_concentration);
    /// @brief Инициализация по переченю компонентов, их концентрации (Eigen::VectorXd) 
    /// и матрице бинарных коэффициентов
    fluid_fundamental_data_t(const std::vector<const component_properties_t*>& components,
        const Eigen::VectorXd& components_concentration,
        const Eigen::MatrixXd& binary_coeffs);
public:
    /// @brief (потокобезопасно) Возвращает вектор мольных концентраций
    const Eigen::VectorXd get_molar_fraction() const;
    /// @brief (потокобезопасно) Задает состав смеси , версия std::vector
    void set_molar_fraction(const std::vector<double>& fraction);
    /// @brief (потокобезопасно) Задает состав смеси (потокобезопасно), версия Eigen::VectorXd
    void set_molar_fraction(const Eigen::VectorXd& fraction);
    /// @brief (потокобезопасно) Возвращает последний результат расчета
    flash_calculation_result_t get_last_flash_result() const;
    /// @brief (потокобезопасно) Основная расчетная задача парожидкостного равновесия
    /// Допускает одновременный вызов flash из нескольких потоков - они выполняется последовательно
    /// @param pressure Рабочее давление
    /// @param temperature Рабочая температура
    /// @param initial_estimation Начальное приближение
    /// @return Результат расчета равновесия
    const flash_calculation_result_t flash(double pressure, double temperature,
                                           double initial_estimation = std::numeric_limits<double>::quiet_NaN()) const;

protected:
    /// @brief Сама эта функция имеет возможность не беспокоиться о потокобезопасности,
    /// Она вызывается после залочивания мьютекса на состав и на last_flash_result
    virtual void flash_unsafe(
            double pressure, double temperature,
            double initial_estimation,
            flash_calculation_result_t& last_flash_result
            ) const = 0;

public:
    /// @brief Возвращает перечень компонентов
    /// Не убирать ссылку, по ней инициализируется fluid_components_functions_t
    const std::vector<const component_properties_t*>& get_components() const;
    /// @brief Возвращает количество компонентов
    size_t get_components_count() const;
    /// @brief Возвращает бинарные коэффициенты компонентов
    const Eigen::MatrixXd& get_binary_coeffs_ref() const;
public:
    /// @brief (потокобезопасно)
    /// @return 
    std::wstring get_compound_name() const {
        std::wstringstream result;
        int i{0};
        for (const auto& component : components_) {
            result << component->name << std::setprecision(3)<<concentration_[i++];
        }
        return result.str();
    }

    /// @brief (потокобезопасно) Возвращает данные для изолированных вызовов
    const fluid_stubdata_t get_mock_data() const;
    /// @brief Заполняет текущее состояние
    /// Технически потокобезопасно, но содержательно можно сломать например параллельным вызовом
    /// set_molar_state
    /// @param state Указатель на структуру, в которую заполнится состояние
    void fill_state(fluid_state_t* state) const;
    /// @brief Задает ранее сохраненное состояние
    /// Технически потокобезопасно, но содержательно можно сломать например параллельным вызовом
    /// set_molar_state
    void set_state(const fluid_state_t& state);
};

/// @brief Геттеры, являющиеся функциями списка компонентов и их параметров
/// Геттеры не являются функциями концентраций, только функции компонентов
/// Нет численных методов
class fluid_components_functions_t {
    /// @brief Здесь исходим из того, что перечень компонентов меняться не будет
    const vector<const component_properties_t*>& components_;
public:
    /// @brief Запрещаем копирование, чтобы не копировать ссылку
    fluid_components_functions_t(const fluid_components_functions_t&) = delete;
    /// @brief Инциализация по списку компонентов
    fluid_components_functions_t(const std::vector<const component_properties_t*>& components);
public:
    /// @brief Возвращает вектор плотностей паровой фазы при рабочих условиях
    virtual Eigen::VectorXd get_densities_vapor(double pressure, double temperature) const = 0;
    /// @brief Возвращает вектор плотностей жидкой фазы при рабочих условиях
    virtual Eigen::VectorXd get_densities_liquid(double pressure, double temperature) const = 0;
public:
    /// @brief Не используется?
    Eigen::VectorXd get_heat_vaporization_by_component(double pressure, double temperature) const;
    /// @brief Вектор молярных масс по каждому компоненту
    Eigen::VectorXd get_molar_masses() const;
    /// @brief Удельные мольные теплоемкости Cp компонентов
    Eigen::VectorXd get_Cp_molar_vapor_by_components(double pressure, double temperature) const;
    /// @brief Удельные массовые теплоемкости Cp компонентов
    Eigen::VectorXd get_Cp_mass_vapor_by_components(double pressure, double temperature) const;
    /// @brief Удельные массовые теплоемкости Cp компонентов в жидкой фазе
    /// Реально зависят только от температуры
    Eigen::VectorXd get_Cp_molar_liquid_by_components(double pressure, double temperature) const;
    /// @brief Получает вектор жидких энтальпий
    /// @tparam amount_type Тип используемых единиц количества вещества (мольные, массовые)
    /// @param pressure Давление
    /// @param temperature Температура
    /// @return Вектор жидких энтальпий
    template <AmountType amount_type>
    Eigen::VectorXd get_enthalpy_liquid_by_components(double pressure, double temperature) const;
    /// @brief Производные по P, T жидких энтальпий по компонентам
    template <AmountType amount_type>
    std::pair<Eigen::VectorXd, Eigen::VectorXd> get_enthalpy_liquid_derivatives_by_components(
            double pressure, double temperature) const
    {
        constexpr double eps = 1e-8;

        auto enthalpy_P = [&](double P) {
            return get_enthalpy_liquid_by_components<amount_type>(P, temperature);
        };
        auto enthalpy_T = [&](double T) {
            return get_enthalpy_liquid_by_components<amount_type>(pressure, T);
        };

        Eigen::VectorXd derivative_P = two_sided_derivative(enthalpy_P, pressure, eps);
        Eigen::VectorXd derivative_T = two_sided_derivative(enthalpy_T, temperature, eps);

        return std::make_pair(std::move(derivative_P), std::move(derivative_T));
    }
    /// @brief Термодинамическая энтальпия для газа
    template <AmountType amount_type>
    Eigen::VectorXd get_enthalpy_vapor_by_components(double pressure, double temperature) const
    {
        size_t components_count = components_.size();
        const auto& components = components_;

        Eigen::VectorXd result(components_count);
        for (int index = 0; index < result.size(); ++index) {
            result(index) = components[index]->get_enthalpy_gas<amount_type>(temperature);
        }
        return result;
    }
    /// @brief Производные по P, T паровых энтальпий по компонентам
    template <AmountType amount_type>
    std::pair<Eigen::VectorXd, Eigen::VectorXd> get_enthalpy_vapor_derivatives_by_components(double pressure, double temperature) const
    {
        constexpr double eps = 1e-8;

        auto enthalpy_P = [&](double P) {
            return get_enthalpy_vapor_by_components<amount_type>(P, temperature);
        };
        auto enthalpy_T = [&](double T) {
            return get_enthalpy_vapor_by_components<amount_type>(pressure, T);
        };

        Eigen::VectorXd derivative_P = two_sided_derivative(enthalpy_P, pressure, eps);
        Eigen::VectorXd derivative_T = two_sided_derivative(enthalpy_T, temperature, eps);

        return std::make_pair(std::move(derivative_P), std::move(derivative_T));
    }
    /// @brief Запускает расчет функции f от свойств каждого компонента
    template <typename Function>
    Eigen::VectorXd get_function_by_components(double pressure, double temperature, Function f) const {
        size_t components_count = components_.size();
        const auto& components = components_;

        Eigen::VectorXd result(components_count);
        for (int index = 0; index < result.size(); ++index) {
            result(index) = f(components[index], pressure, temperature);
        }
        return result;
    }
};

/// @brief Геттеры, являющиеся функциями от состава и свойств компонентов. 
/// Нет численных методов, однопроходные расчеты
class fluid_composition_functions_t {
    /// @brief Ссылка базовые данные флюида
    const fluid_fundamental_data_t& composition;
    /// @brief Ссылка на геттеры компонентов
    const fluid_components_functions_t& components_getters;
public:
    /// @brief Удаляем конструктор копирования, чтобы не копировать поля-ссылки
    fluid_composition_functions_t(const fluid_composition_functions_t&) = delete;
    /// @brief Инициализация по базовым данным флюида и его геттерам компонентов
    explicit fluid_composition_functions_t(const fluid_fundamental_data_t& composition,
                                           const fluid_components_functions_t& components_getters
                                           );
public:
    /// @brief Давления насыщенных паров всех компонентов при заданной температуре по модели Антуана
    /// Формально не зависит от состава, но для компонентов с нулевой концентрацией возвращает нуль!
    virtual Eigen::VectorXd get_saturated_pressures(double temperature) const;
    /// @brief Возвращает вектор K-значений для всех компонентов
    /// Формально не зависит от состава, но в get_saturated_pressures
    /// есть важный костыль, учитывающий состав, поэтому K-values тоже тут
    Eigen::VectorXd get_K_values(double pressure, double temperature) const;
    /// @brief Массовые (не мольные) доли компонентов смеси
    virtual Eigen::VectorXd get_mass_fraction() const;
    /// @brief Возвращает молярную массу смеси
    double get_molar_mass() const;
    /// @brief Газовая постоянная 8.31/M
    double get_gas_constant() const;
    /// @brief Возвращает псевдокритическую температуру
    /// - среднюю критических температуру, взвешенную по концентрациям
    double get_pseudocritical_temperature() const;
    /// @brief Возвращает псевдокритическое давление
    /// - среднее критическое давление, взвешенную по концентрациям
    /// @brief Расчет (псевдо)критического давления
    /// @return Псевдо(еритическое) давление
    double get_pseudocritical_pressure() const;
    /// Возвращает мольный объем флюида, считая, что он находится в газообразном состоянии
    double get_molar_volume_vapor(double pressure, double temperature) const;
    /// Возвращает мольный объем флюида, считая, что он находится в жидком состоянии
    double get_molar_volume_liquid(double pressure, double temperature) const;
    /// @brief Термодинамическая энтальпия смеси в предположении газообразного фазового сосотояния
    /// @param pressure Игнорируется, реализация для идеального газа
    double get_enthalpy_td_mass_as_vapor(double /*pressure*/, double temperature) const;
    /// @brief Термодинамическая энтальпия смеси в предположении жидкофазного сосотояния
    /// @param pressure Игнорируется, реализация для идеального газа
    double get_enthalpy_td_mass_as_liquid(double pressure, double temperature) const;
    /// @brief Расчет удельной мольной внутренней в предположении, что вся смесь в паровом фазовом состоянии
    template <AmountType amount_type>
    double get_inner_energy_as_vapor(double pressure, double temperature) const;
    /// @brief Возвращает среднюю по составу минимальную температурную границу
    /// области определения модели давления насыщенных паров Антуана
    /// Средняя берется по коцентрациям компонентов в составе
    double get_min_antoine_bound() const;
    /// @brief Возвращает среднюю по составу максимальную температурную границу
    /// области определения модели давления насыщенных паров Антуана
    /// Средняя берется по коцентрациям компонентов в составе
    double get_max_antoine_bound() const;

};

/// @brief Критерий фазового состояния (газ, жидкость, жидкость и пар)
class fluid_phase_criteria_t {
    /// @brief Состав
    const fluid_fundamental_data_t& composition;
    /// @brief Геттеры от состава
    const fluid_composition_functions_t& getters;
public:
    /// @brief Удаляем конструктор копирования
    fluid_phase_criteria_t(const fluid_phase_criteria_t&) = delete;
    /// @brief Инициализация состава и геттеров от состава
    explicit fluid_phase_criteria_t(const fluid_fundamental_data_t& composition,
                                    const fluid_composition_functions_t& getters);
    /// @brief Величина критерия на жидкость
    double liquid_only_criteria(double pressure, double temperature) const;
    /// @brief Проверка на жидкость
    /// [01. Антуан Рауль-Дальтон\01b. Документы - актуальные\2022.03.30V03.27 Редактирование]
    bool is_liquid_only(double pressure, double temperature) const;
    /// @brief Величина критерия на газ
    double vapor_only_criteria(double pressure, double temperature) const;
    /// @brief Критерий однофазного газового состояния многокомпонентной смеси
    /// [01. Антуан Рауль-Дальтон\01b. Документы - актуальные\2022.03.30V03.27 Редактирование]
    /// @param pressure Давление
    /// @param temperature Температура
    /// @return Булево значение, которое показывает находится ли вещество в однофазном газовом состоянии
    bool is_vapor_only(double pressure, double temperature) const;
public:
    /// @brief Точка росы при данной температуре
    virtual double get_dew_point_at_given_temperature(double temperature) const = 0;
    /// @brief Точка росы при данном давлении
    virtual double get_dew_point_at_given_pressure(double pressure) const = 0;
    /// @brief Точка росы по воде при данном давлении
    virtual double get_water_dew_point_at_given_pressure(double pressure) const = 0;
    /// @brief Давления начала кипения при данной температуре
    virtual double get_bubble_point_at_given_temperature(double temperature) const = 0;
    /// @brief Давления начала кипения при данном давлении
    virtual double get_bubble_point_at_given_pressure(double pressure) const = 0;
};

/// @brief Геттеры, требующие однократного flash-расчета,
/// являющиеся функциями составов жидкой и/или газообразной фазы
//TODO: Все эти методы обязаны быть методами flash_calculation_result_t
class fluid_flash_functions_t {
private:
    /// @brief Ссылка на класс для flash-расчета
    const fluid_fundamental_data_t& fluid_fundamental;
    /// @brief Ссылка на геттеры компонентов
    const fluid_components_functions_t& components;
    /// @brief Ссылка на геттеры от состава
    const fluid_composition_functions_t& composition;
public:
    /// @brief Удаляем конструктор копирования
    fluid_flash_functions_t(const fluid_flash_functions_t&) = delete;
    /// @brief Инициализация состава и геттеров от состава
    explicit fluid_flash_functions_t(const fluid_fundamental_data_t& fluid_fundamental,
                                     const fluid_components_functions_t& components,
                                     const fluid_composition_functions_t& composition);


    /// @brief Расчет теплоты фазового перехода. Возвращает NaN, если смесь полностью газовая
    double get_heat_vaporization_mass(double pressure, double temperature) const;

    /// @brief Удельная мольная изобарная теплоемкость смеси
    virtual double get_heat_capacity_molar(double pressure, double temperature) const;

    /// @brief Удельная массовая изобарная теплоемкость смеси
    double get_heat_capacity_mass(double pressure, double temperature) const;

    /// @brief Изохорная теплоемкость смеси, мольная
    double get_heat_capacity_isochoric(double pressure, double temperature) const;

    /// @brief Показатель адиабаты
    double get_adiabatic_exponent(double pressure, double temperature) const;
};

/// @brief Методы для создания копий флюидов
class fluid_copy_functions_t {
private:
    /// @brief Создает копию потока с теми же веществами и компонентным составом
    virtual std::unique_ptr<fluid_t> create_copy(bool copy_memoization = true) const = 0;
    /// @brief Создает копию потока с теми же веществами, но другим компонентным составом
    /// @param new_molar_fraction Новый компонентный состав
    virtual std::unique_ptr<fluid_t> create_copy(const Eigen::VectorXd& new_molar_fraction) const = 0;
    /// @brief Создает копию потока с теми же веществами, но другим компонентным составом
    /// @param new_molar_fraction Новый компонентный состав
    virtual std::unique_ptr<fluid_t> create_copy(const std::vector<double>& new_molar_fraction) const = 0;
public:
    /// @brief Создает копию потока с теми же веществами и компонентным составом. Только в конструкторах!
    virtual std::unique_ptr<fluid_t> create_initial_copy(bool copy_memoization = true) const {
        return create_copy(copy_memoization);
    }
    /// @brief Создает копию потока с теми же веществами, но другим компонентным составом. Только в конструкторах!
    /// @param new_molar_fraction Новый компонентный состав
    virtual std::unique_ptr<fluid_t> create_initial_copy(const std::vector<double>& new_molar_fraction) const {
        return create_copy(new_molar_fraction);
    }
    /// @brief Создает копию потока с теми же веществами, но другим компонентным составом
    /// @param new_molar_fraction Новый компонентный состав
    virtual std::unique_ptr<fluid_t> create_initial_copy(const Eigen::VectorXd& new_molar_fraction) const {
        return create_copy(new_molar_fraction);
    }
};

/// @brief Методы для решения задач фазового равновесия при заданном объеме.
/// Могут не являться VT-flash задачами и требовать итеративного решения
class fluid_fill_functions_t {
public:
    /// @brief Заполнить заданный объем веществом с заданной температурой так, чтобы объем жидкости был равен заданному
    /// Выдать давление в заданном объеме
    virtual double fill_volume_with_liquid(double total_volume, double temperature, double liquid_volume) const = 0;

    /// @brief Заполнить заданный объем веществом с заданной температурой так, чтобы давление было равно
    /// @param pressure Давление
    /// @param temperature Температура
    /// @param total_volume Объем емкости
    /// @return Количество вещества
    virtual double fill_volume_with_pressure(double pressure, double temperature,
        double total_volume) const = 0;

    /// @brief Заполнить заданный объем заданным количеством вещества данного состава с заданной температурой
    /// @return давление и температура в заданном объеме
    virtual std::pair<double, double> fill_volume_with_total_moles2(
        double volume, double inner_energy_molar, double molar_amount,
        double initial_pressure = std::numeric_limits<double>::quiet_NaN(),
        double initial_temperature = std::numeric_limits<double>::quiet_NaN()) const = 0;

    /// @brief Изменить объем жидкости за счет газа, сохранив при этом составы жидкости и газа
    /// Если жидкости нет, то выдать ошибку
    virtual void change_liquid_volume_fraction(
        double pressure, double temperature, double liquid_volume_fraction) = 0;
};

/// @brief Абстрактный базовый класс флюида
class fluid_t 
        : public fluid_fundamental_data_t
        , public fluid_components_functions_t
        , public fluid_composition_functions_t
        , public fluid_phase_criteria_t
        , public fluid_flash_functions_t
        , public fluid_copy_functions_t
        , public fluid_fill_functions_t
{
public:
    /// @brief Конструктор копирования
    fluid_t(const fluid_t& other);
    /// @brief Инициализация по перечню компонентов.
    /// Состав инициализируется поровну между компонентами
    fluid_t(const std::vector<const component_properties_t*>& components);
    /// @brief Инициализация по переченю компонентов и их концентрации (Eigen::VectorXd)
    fluid_t(const std::vector<const component_properties_t*>& components,
        const Eigen::VectorXd& components_concentration);
    /// @brief Инициализация по переченю компонентов и их концентрации (std::vector)
    fluid_t(const std::vector<const component_properties_t*>& components,
        const std::vector<double>& components_concentration);
    /// @brief Инициализация по переченю компонентов, их концентрации и матрице бинарных коэффициентов
    fluid_t(const std::vector<const component_properties_t*>& components,
        const Eigen::VectorXd& components_concentration,
        const Eigen::MatrixXd& binary_coeffs);

    virtual ~fluid_t() = default;
};


}
