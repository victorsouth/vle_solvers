#pragma once

using namespace vle_solvers;

/// @brief Тест для примера того, как работать с базой данных
TEST(DataBase, ComponentPropertiesInitialization)
{
	// Инициализируем компонент
	wstring component_name = L"CO2";

	// Объявляем переменную со свойствами компонентов
	component_properties_t component_properties;

	// Инициализируем свойства для заданного компонента
	component_properties = components_database.at(component_name);

	ASSERT_FALSE(std::isnan(component_properties.molar_mass));
}