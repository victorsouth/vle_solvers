#pragma once

/// @brief Тест расчета давления насыщенных паров для гликоля
TEST(ThermoEngine, ResearchGlycole) {
    const auto& data = components_database.at(L"C4H10O3");

    double t1 = data.get_saturated_pressure(vle_solvers::celcium2kelvin(-200));
    double t2 = data.get_saturated_pressure(vle_solvers::celcium2kelvin(-170));
    double t3 = data.get_saturated_pressure(vle_solvers::celcium2kelvin(-10));

    // Проверяем, что расчеты выполнены без ошибок
    ASSERT_FALSE(std::isnan(t1));
    ASSERT_FALSE(std::isnan(t2));
    ASSERT_FALSE(std::isnan(t3));
}

/// @brief Тест сравнения теплоемкости метана
TEST(ThermoDB, Compare1) {
    const auto& c = components_database.at(L"CH4");

    double density_std = 0.7168;

    double Cp = c.get_Cp_gas_mass(KELVIN_OFFSET);

    std::vector<std::wstring> component_list{ L"CH4" };
    auto fluid = vlelib::create_fluid<vlelib::fluid_rault_dalton_t>(component_list);

    double Cpm = fluid->get_heat_capacity_mass(1e5, KELVIN_OFFSET);

    // Проверяем, что расчеты выполнены без ошибок
    ASSERT_FALSE(std::isnan(Cp));
    ASSERT_FALSE(std::isnan(Cpm));
}

TEST(ThermoDB, StdString) {
    std::map<std::string,component_properties_t> other_db;
    other_db["CH4"]=components_database.at(L"CH4");
    std::vector<std::string> component_list{ "CH4" };
    auto fluid = vlelib::create_fluid<vlelib::fluid_rault_dalton_t,std::map<std::string,component_properties_t>,std::string>(
                component_list,{},other_db,{},[](const std::string &str)->std::string{return str;});

}
struct nameCasNo{
    std::string name,CasNo;
};
bool operator<(const struct nameCasNo&a,const struct nameCasNo&b){
    return a.CasNo<b.CasNo;
}
TEST(ThermoDB, ComplexStruct) {
    const auto& c = components_database.at(L"CH4");
    nameCasNo key{fixed_solvers::wide2string(c.name),fixed_solvers::wide2string(c.CASno)};
    std::map<nameCasNo,component_properties_t> other_db;
    other_db[key]=c;
    other_db[{"CH4","89898956"}]=c;
    std::vector<nameCasNo> component_list{ key };
    {
        auto fluid = vlelib::create_fluid<vlelib::fluid_rault_dalton_t,std::map<nameCasNo,component_properties_t>,nameCasNo>(
                component_list,{},other_db,{},[](const nameCasNo &str)->std::string{return str.name+str.CasNo;});
    }
#if 0
    nameCasNo not_found{"swds","zsdfasdf"};
    std::vector<nameCasNo> component_list2{not_found};
    ASSERT_THROW({//что-то тут не собирается
        auto fluid = vlelib::create_fluid<vlelib::fluid_rault_dalton_t,std::map<nameCasNo,component_properties_t>,nameCasNo>(
                     component_list2,{},other_db,{},[](const nameCasNo &str)->std::string{return str.name+str.CasNo;});
    },std::domain_error);
#endif
}
