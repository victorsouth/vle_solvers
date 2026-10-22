#pragma once

/// @brief Тест для примера того, как работать с базой данных
TEST(DataBase, ComponentPropertiesInitialization)
{
	// Инициализируем компонент
	std::wstring component_name = L"CO2";

	// Инициализируем свойства для заданного компонента
	const component_properties_t& component_properties = components_database.get_component_by_formula(component_name);
	ASSERT_FALSE(std::isnan(component_properties.molar_mass));
}

/// @brief Компонент находится одинаково по CAS-номеру и по химической формуле.
TEST(ThermoDbComponentsGetter, ResolvesCasAndFormulaToSameComponent)
{
    // Arrange: метан из глобальной БД
    const component_properties_t& methane =
        components_database.get_component_by_formula(L"CH4");

    // Act: запрос по CAS, по формуле. Смешанный запрос CAS+формула
    std::vector<const component_properties_t*> by_cas =
        components_database.get_components({ methane.CASno });
    std::vector<const component_properties_t*> by_formula =
        components_database.get_components({ L"CH4" });
    std::vector<const component_properties_t*> mixed =
        components_database.get_components({ methane.CASno, L"CH4" });

    // Assert: во всех случаях БД отдает один и тот же компонент
    ASSERT_EQ(by_cas.size(), 1u);
    ASSERT_EQ(by_formula.size(), 1u);
    EXPECT_EQ(by_cas[0], by_formula[0]);
    EXPECT_EQ(by_cas[0]->CASno, methane.CASno);
    ASSERT_EQ(mixed.size(), 2u);
    EXPECT_EQ(mixed[0], mixed[1]);
}

/// @brief Неизвестный CAS или формула с изотопами (неединственным CAS) при запросе компонентов дают ошибку.
TEST(ThermoDbComponentsGetter, ThrowsForUnknownAndAmbiguousFormula)
{
    // Arrange: локальная БД с двумя веществами на одну формулу
    component_properties_t isomer_a =
        components_database.get_component_by_formula(L"CH4");
    component_properties_t isomer_b = isomer_a;
    isomer_a.name = L"CX";
    isomer_a.CASno = L"111-11-1";
    isomer_b.name = L"CX";
    isomer_b.CASno = L"222-22-2";
    components_database_t pure_db{
        { isomer_a.CASno, isomer_a },
        { isomer_b.CASno, isomer_b }
    };
    thermo_db_t local_db(pure_db, bip_records_t{}, components_database_t{});

    // Act: запрос по неизвестному CAS
    // Assert: формируется исключение
    EXPECT_THROW(
        components_database.get_components({ L"NOT_A_REAL_COMPONENT_ID" }),
        std::runtime_error);

    // Act: запрос по неоднозначной формуле
    // Assert: формируется исключение
    EXPECT_THROW(local_db.get_components({ L"CX" }), std::runtime_error);
}

/// @brief Смесь, заданная только формулами или только CAS-номерами, даёт один и тот же набор свойств компонентов.
TEST(ThermoDbComponentsGetter, ReturnsSamePropertiesForFormulasAndCas)
{
    // Arrange: смесь CH4–C2H6–N2 в виде формул и соответствующих CAS
    std::vector<std::wstring> formulas{ L"CH4", L"C2H6", L"N2" };
    std::vector<std::wstring> casno_list =
        components_database.get_casno_by_formulas(formulas);

    // Act: получение свойств по формулам и по CAS
    std::vector<const component_properties_t*> by_formula =
        components_database.get_components(formulas);
    std::vector<const component_properties_t*> by_cas =
        components_database.get_components(casno_list);

    // Assert: одинаковый состав и порядок компонентов
    ASSERT_EQ(by_formula.size(), by_cas.size());
    ASSERT_EQ(by_formula.size(), formulas.size());
    for (std::size_t i = 0; i < by_formula.size(); ++i) {
        EXPECT_EQ(by_formula[i], by_cas[i]);
        EXPECT_EQ(by_formula[i]->CASno, casno_list[i]);
    }
}

