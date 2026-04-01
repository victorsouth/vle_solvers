#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <numeric>
#include <utility>
#include <string>
#include <mutex>
#include <shared_mutex>
#include <iomanip>
#include <fstream>
#include <memory>
#include <Eigen/Dense>

#include <fixed/fixed.h>
#include <fixed/fixed_bisection.h>

#include "physical_constants.h"

#include "helpers/physical_helpers.h"


/// включение локали для чтения nan inf -inf
#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#include <boost/type_index.hpp>
#ifdef VLELIB_SERIALIZATION_SUPPORT
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#endif



#include "components_db.h"
#include "fluid/fluid_common.h"
#include "fluid/fluid_base.h"
#include "fluid/fluid_raoult_dalton.h"
#include "fluid/fluid_equations.h"

