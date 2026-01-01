#include "../vle_solvers.h"

#define BOOST_EXCEPTION_DISABLE
#include "tdb_serialize_2024_09_24.h"

//#include <common/common_hydraulics.h>
//using hydraulics::celcium2kelvin;

const components_database_t components_database;

/// @brief Структура для обвязки инициализации глобавльной базы данных по компонентам
struct db_initializer_t {
    /// @brief Инициализация глобальной базы из глобальной строки
    db_initializer_t(const components_database_t& _db)
    {
        components_database_t& db = const_cast<components_database_t&>(_db);
        db = serializer_2024_09_24::deserialize_from_string<components_database_t>(std::string(thermo_db_serialized));

        override_glycole_data(db);
        override_methane_data(db);
        calc_extrapolation_coeff(db);
    }

    /// @brief Исправление свойств метана
    void override_methane_data(components_database_t& db)
    {
        auto& component = db[L"CH4_"] = db.at(L"CH4");
        auto& model = component.antoine_model;

        model.antoine_coefficients = {
            31.35, -1307.52, 0, -3.26134, 0.000029418, 2
        };

        model.formula = antoine_formula::ExtLnKPaKelvin;
        model.min_bound = vle_solvers::celcium2kelvin(-182.15);
        model.max_bound = vle_solvers::celcium2kelvin(-82.75299683);

    }

    /// @brief Исправление свойств гликоля
    void override_glycole_data(components_database_t& db)
    {
        db[L"C4H10O3_"] = db.at(L"C4H10O3");

        auto& model = db.at(L"C4H10O3_").antoine_model;

        model.formula = antoine_formula::LgMmHgCelcium;
        model.antoine_coefficients = { 7.65732, 2065.8762, 186.657 };
        model.min_bound = vle_solvers::celcium2kelvin(123.66);
        model.max_bound = vle_solvers::celcium2kelvin(274.94);


    }

    /// @brief Расчёт коэффициента экстраполяции модели Антуана
    void calc_extrapolation_coeff(components_database_t& db)
    {
        for (auto& [name, data] : db)
        {
            double my_estimation = data.estimate_antoine_extrapolation_coeff();
            data.antoine_model.extrapolation_coefficient = my_estimation;
        }

    }

};

db_initializer_t db_init(components_database);
