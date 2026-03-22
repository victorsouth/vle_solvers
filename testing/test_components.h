#pragma once

/// @brief Тест расчета давления насыщенных паров для гликоля
TEST(ThermoEngine, ResearchGlycole) {


    const auto data = components_database.get_component_by_formula(L"C4H10O3");
    ASSERT_FALSE(data==nullptr);

    double t1 = data->get_saturated_pressure(vle_solvers::celcium2kelvin(-200));
    double t2 = data->get_saturated_pressure(vle_solvers::celcium2kelvin(-170));
    double t3 = data->get_saturated_pressure(vle_solvers::celcium2kelvin(-10));

    // Проверяем, что расчеты выполнены без ошибок
    ASSERT_FALSE(std::isnan(t1));
    ASSERT_FALSE(std::isnan(t2));
    ASSERT_FALSE(std::isnan(t3));
}

/// @brief Тест сравнения теплоемкости метана
TEST(ThermoDB, Compare1) {
    const auto component_ptr = components_database.get_component_by_formula(L"CH4");
    ASSERT_FALSE(component_ptr == nullptr);
    const auto& c = *component_ptr;

    double density_std = 0.7168;

    double Cp = c.get_Cp_gas_mass(KELVIN_OFFSET);

    std::vector<std::wstring> component_list{ L"CH4" };
    auto fluid = vlelib::create_fluid<vlelib::fluid_rault_dalton_t>(component_list);

    double Cpm = fluid->get_heat_capacity_mass(1e5, KELVIN_OFFSET);

    // Проверяем, что расчеты выполнены без ошибок
    ASSERT_FALSE(std::isnan(Cp));
    ASSERT_FALSE(std::isnan(Cpm));
}
