#pragma once

#include <string>
#include <boost/type_index.hpp>

namespace vle_solvers {

/// @brief Получает имя класса из типа объекта
/// @param t объект для определения типа
/// @return строка с именем класса
template<typename T>
inline std::string get_class_as_string(const T& t) {
    std::string str = boost::typeindex::type_id<decltype(t)>().pretty_name();
#ifdef _MSC_VER
    auto ppos = str.find_last_of(" \t");
    if (ppos != str.npos)
        str = str.substr(ppos + 1);
#endif
    return str;
}

}
