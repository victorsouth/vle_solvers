#pragma once
#include "tdb_serialize_common.h"

#include <optional>

/* Сериализация 18.11.2020, добавление параметров:
       корреляции Ван-Вельцена B, T0, 
       газокинетический диаметр,
       равновесная энергия
*/
class serializer_2024_09_24 {
private:
    static void serialize_van_velzen(const van_velzen_viscosity_correlation& van_velzen,
        boost::property_tree::ptree& props)
    {
        props.put("van_velzen_B", van_velzen.B);
        props.put("van_velzen_T0", van_velzen.T0);
    }

    static void serialize_antoine(const antoine_model_t& antoine,
        boost::property_tree::ptree& props)
    {
        props.put("antoine_min_bound", antoine.min_bound);
        props.put("antoine_max_bound", antoine.max_bound);

        props.put("antoine_extrapolation_coefficient", antoine.extrapolation_coefficient);
        props.put("antoine_formula", static_cast<int>(antoine.formula));

        props.put("antoine_A", antoine.antoine_coefficients[0]);
        props.put("antoine_B", antoine.antoine_coefficients[1]);
        props.put("antoine_C", antoine.antoine_coefficients[2]);
    }
    static void serialize_td_functions(const thermodynamic_functions_t& functions,
        boost::property_tree::ptree& props)
    {
        const auto& ranges = functions.get_ranges();

        for (size_t index = 0; index < ranges.size(); ++index) {
            boost::property_tree::ptree tdf_diapason;
            tdf_diapason.put("range_start", ranges[index].range_start);
            tdf_diapason.put("range_end", ranges[index].range_end);
            tdf_diapason.put("enthalpy", ranges[index].coefficients.enthalpy);
            tdf_diapason.put("entropy", ranges[index].coefficients.entropy);

            boost::property_tree::ptree capacity_coeffs;
            for (double c : ranges[index].coefficients.heat_capacity) {
                boost::property_tree::ptree cell;
                cell.put_value(c);
                capacity_coeffs.push_back(std::make_pair("", cell));
            }

            tdf_diapason.add_child("heat_capacity_vapor", capacity_coeffs);

            props.add_child(
                std::string("td_function_diapason_") + fixed_solvers::int2str(index + 1),
                tdf_diapason);
        }
    }
    static void serialize_heat_capacity_liquid(
        const fixed_solvers::ranged_polynom_t<heat_capacity_coefficients_t>& heat_capacity_liquid,
        boost::property_tree::ptree& props)
    {
        const auto& ranges = heat_capacity_liquid.get_ranges();

        for (size_t index = 0; index < ranges.size(); ++index) {
            boost::property_tree::ptree diapason;
            diapason.put("range_start", ranges[index].range_start);
            diapason.put("range_end", ranges[index].range_end);

            boost::property_tree::ptree capacity_coeffs;
            for (double c : ranges[index].coefficients) {
                boost::property_tree::ptree cell;
                cell.put_value(c);
                capacity_coeffs.push_back(std::make_pair("", cell));
            }

            diapason.add_child("heat_capacity", capacity_coeffs);

            props.add_child(
                std::string("heat_capacity_liquid_") + fixed_solvers::int2str(index + 1),
                diapason);
        }
    }

    static void serialize_component(const std::wstring& name, const component_properties_t& properties,
        boost::property_tree::ptree& root)
    {
        boost::property_tree::ptree props;

        props.put("molar_mass", properties.molar_mass);
        props.put("acentric_factor", properties.acentric_factor);
        props.put("critical_pressure", properties.critical_pressure);
        props.put("critical_molarvolume", properties.critical_molarvolume);
        props.put("critical_temperature", properties.critical_temperature);
        props.put("density_liquid_20", properties.density_liquid_20);
        props.put("normal_boiling_temperature", properties.normal_boiling_temperature);
        props.put("condensation_heat_molar", properties.condensation_heat_molar);
        props.put("elastic_modulus", properties.elastic_modulus);
        props.put("gas_kinetic_diameter", properties.gas_kinetic_diameter);
        props.put("equilibrium_energy", properties.equilibrium_energy);

        serialize_antoine(properties.antoine_model, props);
        serialize_td_functions(properties.functions, props);
        serialize_heat_capacity_liquid(properties.heat_capacity_liquid, props);
        serialize_van_velzen(properties.viscosity_correlation, props);
        root.add_child(fixed_solvers::wide2string(name), props);

    }
    static std::vector<std::string> string_list_by_lines(const std::string& s) {
        std::stringstream ss(s);

        std::vector<std::string> result;
        std::string line;
        while (std::getline(ss, line, '\n')) {
            result.push_back(line);
            //cout << to << endl;
        }
        return result;
    }
    static std::string serialize_as_string(const boost::property_tree::ptree& root) {

        std::stringstream ss;
        boost::property_tree::write_json(ss, root);

        std::vector<std::string> strings = string_list_by_lines(ss.str());

        std::stringstream result;

        for (std::string& s : strings) {
            fixed_solvers::string_replace(s, "\"", "\\\""); // одинарные кавычки '"' на '\"'
            s = std::string("\"") + s + std::string("\""); // кавычки в начало и в конец

            result << s << std::endl;
        }

        return result.str();
    }


    static thermodynamic_functions_t deserialize_td_functions(const boost::property_tree::ptree& props)
    {
        vector<fixed_solvers::function_range_t<thermodynamic_functions_coefficients_t>> ranges;

        for (size_t index = 1; true; ++index) {
            //try 
            {
                std::string child_name = std::string("td_function_diapason_") + fixed_solvers::int2str(index);
                auto child = props.get_child_optional(child_name);
                if (!child)
                {
                    break;
                }

                boost::property_tree::ptree tdf_diapason = child.value();// props.get_child(child_name);

                fixed_solvers::function_range_t<thermodynamic_functions_coefficients_t> range;

                range.range_start = get_double_value(tdf_diapason.get_child("range_start"));
                range.range_end = get_double_value(tdf_diapason.get_child("range_end"));

                range.coefficients.enthalpy = tdf_diapason.get<double>("enthalpy");
                range.coefficients.entropy = tdf_diapason.get<double>("entropy");

                boost::property_tree::ptree heat_capacity = tdf_diapason.get_child("heat_capacity_vapor");

                for (const auto& heat_capacity_coeff : heat_capacity) {
                    range.coefficients.heat_capacity.push_back(heat_capacity_coeff.second.get_value<double>());
                }

                ranges.push_back(range);
            }
            /*catch (...) {
                break;
            }*/
        }

        return ranges;
    }

    static fixed_solvers::ranged_polynom_t<heat_capacity_coefficients_t> deserialize_heat_capacity_liquid(
        const boost::property_tree::ptree& props)
    {
        vector<fixed_solvers::function_range_t<heat_capacity_coefficients_t>> ranges;

        for (size_t index = 1; true; ++index) {
            try {


                std::string child_name = std::string("heat_capacity_liquid_") + fixed_solvers::int2str(index);
                auto child = props.get_child_optional(child_name);
                if (!child) {
                    break;
                }

                boost::property_tree::ptree diapason = child.value();
                    //props.get_child(string("heat_capacity_liquid_") + int2str(index));


                fixed_solvers::function_range_t<heat_capacity_coefficients_t> range;

                range.range_start = get_double_value(diapason.get_child("range_start"));
                range.range_end = get_double_value(diapason.get_child("range_end"));

                boost::property_tree::ptree heat_capacity = diapason.get_child("heat_capacity");

                for (const auto& heat_capacity_coeff : heat_capacity) {
                    range.coefficients.push_back(heat_capacity_coeff.second.get_value<double>());
                }

                ranges.push_back(range);
            }
            catch (...) {
                break;
            }
        }

        return ranges;
    }

    static antoine_model_t deserialize_antoine(const boost::property_tree::ptree& props)
    {
        antoine_model_t antoine;
        antoine.min_bound = get_double_value(props.get_child("antoine_min_bound"));
        antoine.max_bound = get_double_value(props.get_child("antoine_max_bound"));

        antoine.extrapolation_coefficient = get_double_value(props.get_child("antoine_extrapolation_coefficient"));
        antoine.formula = static_cast<antoine_formula>(props.get<int>("antoine_formula"));

        antoine.antoine_coefficients[0] = props.get<double>("antoine_A");
        antoine.antoine_coefficients[1] = props.get<double>("antoine_B");
        antoine.antoine_coefficients[2] = props.get<double>("antoine_C");

        return antoine;
    }

    static van_velzen_viscosity_correlation deserialize_van_velzen(const boost::property_tree::ptree& props)
    {
        van_velzen_viscosity_correlation van_velzen;
        van_velzen.B = get_double_value(props.get_child("van_velzen_B"));
        van_velzen.T0 = get_double_value(props.get_child("van_velzen_T0"));
        return van_velzen;
    }


    static component_properties_t deserialize_component(const std::wstring& name, const boost::property_tree::ptree& props)
    {
        component_properties_t properties;

        properties.name = name;
        properties.molar_mass = props.get<double>("molar_mass");
        properties.acentric_factor = props.get<double>("acentric_factor");
        properties.critical_pressure = props.get<double>("critical_pressure");
        properties.critical_molarvolume = props.get<double>("critical_molarvolume");
        properties.critical_temperature = props.get<double>("critical_temperature");
        properties.density_liquid_20 = props.get<double>("density_liquid_20");
        properties.condensation_heat_molar = props.get<double>("condensation_heat_molar");
        properties.elastic_modulus = props.get<double>("elastic_modulus");
        properties.normal_boiling_temperature = props.get<double>("normal_boiling_temperature");

        properties.gas_kinetic_diameter = get_double_value(props.get_child("gas_kinetic_diameter"));//props.get<double>("gas_kinetic_diameter");
        properties.equilibrium_energy = get_double_value(props.get_child("equilibrium_energy"));//props.get<double>("equilibrium_energy");

        properties.antoine_model = deserialize_antoine(props);
        properties.functions = deserialize_td_functions(props);
        properties.heat_capacity_liquid = deserialize_heat_capacity_liquid(props);
        properties.viscosity_correlation = deserialize_van_velzen(props);

        return properties;
    }
    template<typename components_database_type>
    static components_database_type deserialize(const boost::property_tree::ptree& root)
    {
        components_database_type db;

        for (const auto& child : root) {
            std::wstring name = fixed_solvers::string2wide(child.first);
            db[name] = deserialize_component(name, child.second);
        }
        return db;
    }

 public:
    /// @brief Сериализует базу данных компонентов в JSON файл
    template<typename components_database_type>
    static void serialize_as_json(const std::wstring& filename, const components_database_type& db) {
        boost::property_tree::ptree root;
        for (const auto& kvp : db) {
            serialize_component(kvp.first, kvp.second, root);
        }

        std::stringstream ss;
        boost::property_tree::write_json(fixed_solvers::wide2string(filename), root);
    }
    /// @brief Сериализует базу данных компонентов по формату, 
    /// используемому в components_db_data.cpp
    template<typename components_database_type>
    static void serialize_as_cpp(const std::wstring& path, const components_database_type& db) {
        boost::property_tree::ptree root;
        for (const auto& kvp : db) {
            serialize_component(kvp.first, kvp.second, root);
        }
        std::string str_result =
            std::string("const char* thermo_db_serialized = \n") +
            serialize_as_string(root) + ";";

        std::ofstream code_file(fixed_solvers::wide2string(path) + "components_db_data.cpp");
        code_file << str_result;
        code_file.flush();
    }
    /// @brief Десериализует базу данных компонентов из JSON строки
    template<typename components_database_type>
    static components_database_type deserialize_from_string(const std::string& serialized_str)
    {
        boost::property_tree::ptree root;

        std::stringstream serialized_stream;
        serialized_stream << serialized_str;

        boost::property_tree::read_json(serialized_stream, root);

        return deserialize<components_database_type>(root);
    }

    /// @brief Десериализует базу данных компонентов из JSON файла
    template<typename components_database_type>
    static components_database_type deserialize(const std::wstring& filename)
    {
        boost::property_tree::ptree root;
        boost::property_tree::read_json(fixed_solvers::wide2string(filename), root);

        return deserialize<components_database_type>(root);
    }

};
