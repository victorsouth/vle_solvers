#include "../vle_solvers.h"

#define BOOST_EXCEPTION_DISABLE
#include "tdb_serialize_2026_02_02.h"

//#include <common/common_hydraulics.h>
//using hydraulics::celcium2kelvin;

const components_database_t components_database;

/// @brief Структура для обвязки инициализации глобавльной базы данных по компонентам
struct db_initializer_t {
    /// @brief Инициализация глобальной базы из глобальной строки
    db_initializer_t(const components_database_t& _db)
    {
        components_database_t& db = const_cast<components_database_t&>(_db);
        db = serializer_2026_02_02::deserialize_from_string<components_database_t>(std::string(thermo_db_serialized));

        calc_extrapolation_coeff(db);
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
