#pragma once


namespace vle_solvers {

inline double celcium2kelvin(double celcium_value) {
    return celcium_value + KELVIN_OFFSET;
}

inline double kelvin2celcium(double kelvin_value) {
    return kelvin_value - KELVIN_OFFSET;
}

inline double circle_area(double diameter) {
    return M_PI * diameter * diameter / 4;
}

}

