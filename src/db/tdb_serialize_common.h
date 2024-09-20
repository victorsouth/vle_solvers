#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "components_db.h"

namespace vle_solvers
{
;

inline double get_double_value(const boost::property_tree::ptree& node) {
    auto str = node.get_value<std::string>();
    if (str == "nan") {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else if (str == "-inf") {
        return -std::numeric_limits<double>::infinity();
    }
    else if (str == "inf") {
        return std::numeric_limits<double>::infinity();
    }
    else {
        return node.get_value<double>();
    }
}

}