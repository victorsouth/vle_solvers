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

