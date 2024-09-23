// Подключение бибилотеки vle_solvers.h
#include "../src/vle_solvers.h"

// Используем необходимое пространство имен
using namespace vle_solvers;


int main()
{
	// Выбираем любой компонент, для которого необходимо определить свойства
	wstring component = L"CO2";

	// Возвращаем свойства компонента
	component_properties_t component_properties = components_database.at(component);

	return 1;
}