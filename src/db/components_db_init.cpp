#include "components_db.h"

#define BOOST_EXCEPTION_DISABLE
#include "tdb_serialize_2024_05_24.h"


namespace vle_solvers
{
;



const components_database_t components_database;

struct db_initializer_t {
    db_initializer_t(const components_database_t& _db)
    {
        extern const char* thermo_db_serialized;
        components_database_t& db = const_cast<components_database_t&>(_db);
        db = serializer_2024_05_24::deserialize_from_string<components_database_t>(string(thermo_db_serialized));

        override_glycole_data(db);
        override_methane_data(db);
        calc_extrapolation_coeff(db);
    }

    void override_methane_data(components_database_t& db)
    {
        auto& component = db[L"CH4_"] = db.at(L"CH4");
        auto& model = component.antoine_model;

        model.antoine_coefficients = {
            31.35, -1307.52, 0, -3.26134, 0.000029418, 2
        };

        model.formula = antoine_formula::ExtLnKPaKelvin;
        model.min_bound = celcium2kelvin(-182.15);
        model.max_bound = celcium2kelvin(-82.75299683);

    }

    void override_glycole_data(components_database_t& db)
    {
        db[L"C4H10O3_"] = db.at(L"C4H10O3");

        auto& model = db.at(L"C4H10O3_").antoine_model;

        model.formula = antoine_formula::LgMmHgCelcium;
        model.antoine_coefficients = { 7.65732, 2065.8762, 186.657 };
        model.min_bound = celcium2kelvin(123.66);
        model.max_bound = celcium2kelvin(274.94);


    }

    void calc_extrapolation_coeff(components_database_t& db)
    {
        for (auto& [name, data] : db)
        {
            double my_estimation = data.estimate_antoine_exptraploation_coeff();
            data.antoine_model.extrapolation_coefficient = my_estimation;
        }

    }

};

db_initializer_t db_init(components_database);


}
