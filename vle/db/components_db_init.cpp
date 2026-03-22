#include "../vle_solvers.h"

#define BOOST_EXCEPTION_DISABLE
#include "tdb_serialize_2026_02_02.h"

//#include <common/common_hydraulics.h>
//using hydraulics::celcium2kelvin;

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
    calc_extrapolation_coeff();
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
    calc_extrapolation_coeff();
}



void thermo_db_t::calc_extrapolation_coeff()
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
        auto key = binary_casno_t::make_key(rec.cas1_m, rec.cas2_m);
        bips[key] = rec.k_ij_m;
    }
}



const components_database_t& thermo_db_t::get_component_casno_database() const
{
    return components;
}



std::optional<std::wstring>  thermo_db_t::get_casno_by_formula(const std::wstring& formula) const
{
    if (formula.empty()) {
        return std::nullopt;
    }
    int ncomps = formula2cas_mapping_var1.count(formula);
    if (ncomps == 0) {
        return std::nullopt;
    }
    if (ncomps > 1) {
        throw std::wstring(L"Fatal error! There are several components for chemical formula ")
            + formula;
    }
    else {
        auto iter = formula2cas_mapping_var1.find(formula);
        return iter->second;
    };
}



const component_properties_t* thermo_db_t::get_component_by_formula(
    const std::wstring& formula) const
{
    // будет исключение, если не найдём
    auto casno = get_casno_by_formula(formula);
    if (!casno) {
        return nullptr;
    }
    return get_component_by_casno(casno.value());
}



const component_properties_t* thermo_db_t::get_component_by_casno(
    const std::wstring& casno) const
{
    int ncomp = components.count(casno);
    if ( !ncomp ) {
        return nullptr;
    }
    return &components.at(casno);
}



std::optional<double> thermo_db_t::get_bip_pair_formula(const binary_formula_t& pair_formula) const
{
    auto key = binary_casno_t::make_key(pair_formula);
    return get_bip_pair_casno(key);
}



std::optional<double>  thermo_db_t::get_bip_pair_casno(const binary_casno_t& pair_casno) const
{
    auto key = binary_casno_t::make_key(pair_casno);
    auto iter = bips.find(key);
    if (iter == bips.end()) {
        return std::nullopt;
    }
    return iter->second;
}



bool binary_casno_t::operator<(const binary_casno_t& other) const
{
    return std::tie(cas1_m, cas2_m) < std::tie(other.cas1_m, other.cas2_m);
}



binary_casno_t binary_casno_t::make_key(const std::wstring& cas1, const std::wstring& cas2)
{
    if (cas1 < cas2) {
        return { cas1, cas2 };
    }
    else {
        return { cas2, cas1 };
    }
}



binary_casno_t binary_casno_t::make_key(binary_casno_t pair_casno)
{
    return make_key(pair_casno.cas1_m, pair_casno.cas2_m);
}



binary_casno_t binary_casno_t::make_key(binary_formula_t pair_formula)
{
    return make_key(pair_formula.get_binary_cas());
}



binary_casno_t binary_formula_t::get_binary_cas() const
{
    try {
        auto cas1 = components_database.get_casno_by_formula(formula1);
        auto cas2 = components_database.get_casno_by_formula(formula2);
        if (!cas1 || !cas2) {
            throw std::wstring(L"Fatal error! Not found one of chemical formula: ")
                + formula1 + L", " + formula2;
        }
        return binary_casno_t::make_key(*cas1, *cas2);
    }
    catch (...) {
        // специализируем, что исключение при поиске бинарного ключа из химических формул
        throw std::wstring (L"Fatal error! CAS number not unique for either " + formula1
            + L" or for " + formula2);
    }
}



bool binary_formula_t::operator<(const binary_formula_t& other) const {
    auto pair = get_binary_cas();
    auto pai_other = other.get_binary_cas();
    return pair.operator<(pai_other);
}

