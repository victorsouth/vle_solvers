#include "../vle_solvers.h"

#define BOOST_EXCEPTION_DISABLE
#include "tdb_serialize_2026_02_02.h"


const components_database_t components_database_by_formula;


/// @brief Структура для обвязки инициализации глобавльной базы данных по компонентам
struct db_initializer_t {
    /// @brief Инициализация глобальной базы из глобальной строки
    db_initializer_t(const components_database_t& _db)
    {
        components_database_t& db = const_cast<components_database_t&>(_db);
        db = serializer_2026_02_02::deserialize_from_string<components_database_t>(
            std::string(thermo_db_serialized_by_formula));

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


static thermo_db_t make_db_old()
{
    // база данных собранная по химической формуле
    db_initializer_t db_init(components_database_by_formula);
    thermo_db_t db = thermo_db_t();
    return db;
}



static thermo_db_t make_db_new()
{
    // база данных собранная по химической формуле
    //thermo_db_t db = thermo_db_t(std::string(thermo_db_serialized_by_casno));
    //return db;
    return {};
}


const thermo_db_t components_database = make_db_old();
//const components_database_t hypocomponents_database{};



thermo_db_t::thermo_db_t(const std::string& thermo_db_by_casno)
{
    throw("thermo_db_t::thermo_db_t(const std::string& thermo_db_by_casno) not realized");
    components = serializer_2026_02_02::deserialize_from_string<components_database_t>(
        std::string(thermo_db_by_casno));
    formula2cas_mapping_var1.clear();
    for (const auto& component_by_casno : components) {
        const auto& chemical_formula = component_by_casno.second.name;
        formula2cas_mapping_var1.insert({ chemical_formula, component_by_casno.first});
    }
    init_bips(bip_records_global);
    init_extrapolation_coeff();
}



thermo_db_t::thermo_db_t()
{
    // так как собираем здесь поэлементно, то очищаем.
    components.clear();
    formula2cas_mapping_var1.clear();
    // собираем из текущий базы компонентов (которая по формулам)
    for (const auto& component_by_formula : components_database_by_formula) {
        const auto& chemical_formula = component_by_formula.first;
        const auto& casno = component_by_formula.second.CASno;
        components[casno] = component_by_formula.second;
        formula2cas_mapping_var1.insert({ chemical_formula, casno });
    }
    init_bips(bip_records_global);
    init_extrapolation_coeff();
}



void thermo_db_t::init_extrapolation_coeff()
{
    for (auto& [name, data] : components)
    {
        double my_estimation = data.estimate_antoine_extrapolation_coeff();
        data.antoine_model.extrapolation_coefficient = my_estimation;
    }
}



void thermo_db_t::init_bips(const bip_records_t& bip_records)
{
    for (const auto& rec : bip_records) {
        std::set<std::wstring> key = { rec.cas1, rec.cas2 };
        bips[key] = rec.bip_value;
    }
}

const components_database_t& thermo_db_t::get_component_casno_database() const
{
    return components;
}

const component_properties_t& thermo_db_t::get_component_by_formula(const std::wstring& formula) const
{
    return get_component_by_casno(get_casno_by_formula(formula));
}

const component_properties_t& thermo_db_t::get_component_by_casno(const std::wstring& casno) const
{
    int ncomp = components.count(casno);
    if (!ncomp) {
        throw std::runtime_error("CASno not found");
    }
    return components.at(casno);
}

const std::wstring& thermo_db_t::get_casno_by_formula(const std::wstring& formula) const
{
    if (formula.empty()) {
        throw std::runtime_error("Cannot find formula");
    }
    std::size_t ncomps = formula2cas_mapping_var1.count(formula);
    if (ncomps == 0) {
        throw std::runtime_error("Cannot find formula");
    }
    if (ncomps > 1) {
        throw std::runtime_error("Fatal error! There are several components for chemical formula ");
    }
    else {
        auto iter = formula2cas_mapping_var1.find(formula);
        return iter->second;
    };
}



double thermo_db_t::get_bip_pair_formula(const std::wstring& formula1, const std::wstring& formula2) const
{
    std::wstring cas1 = get_casno_by_formula(formula1);
    std::wstring cas2 = get_casno_by_formula(formula2);
    return get_bip_pair_casno(cas1, cas2);
}

double thermo_db_t::get_bip_pair_casno(const std::wstring& cas1, const std::wstring& cas2) const
{
    const auto key = std::set<std::wstring>{ cas1, cas2 };
    auto iter = bips.find(key);
    if (iter == bips.end()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return iter->second;
}