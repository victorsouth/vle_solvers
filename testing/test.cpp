// Подключение бибилотеки vle_solvers.h
#include "../src/vle_solvers.h"

// Используем необходимое пространство имен
using namespace vle_solvers;


int main()
{
	// Создает переменную, которая будет содержать базу данных компонентов и инициализируем ее
	components_database_t components_database;
	db_initializer_t db(components_database);

	// Выбираем любой компонент, для которого необходимо определить свойства
	wstring component = L"CO2";

	// Возвращаем свойства компонента
	component_properties_t component_properties = components_database.at(component);

	return 1;
}