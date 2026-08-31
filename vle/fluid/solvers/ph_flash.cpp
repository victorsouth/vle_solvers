#include "../../vle_solvers.h" // очень некрасивый инклуд...

namespace vlelib {
;

ph_flash_bisection::ph_flash_bisection(fluid_t* f, double de, double p, double t)
    : fluid(f), given_enthalpy(de), given_pressure(p), initial_temperature(t)
{
}

ph_flash_bisection::function_type ph_flash_bisection::residuals(const var_type& x)
{
    return given_enthalpy - fluid->flash(given_pressure, x).td_functions.enthalpy.mass.mix;
}

double ph_flash_bisection::solve(
    double desired_precision,
    fixed_bisection_result_t<1>* numerical_result,
    fixed_bisection_result_analysis_t<1>* analysis)
{
    fixed_bisection_result_t<1> res;
    if (numerical_result == nullptr)
        numerical_result = &res;

    fixed_bisectional_parameters_t p;
    double T_critical = fluid->get_pseudocritical_temperature();
    if (fluid->is_ideal_gas())
    {
        //p.argument_limit_min = 10.0; // исходный вариант
        p.argument_limit_min = std::max(std::min(10., initial_temperature / 2.), 2.); // вариант ЮП
    }
    else {
        // PREOS / PR: T=10 K ломает Michelsen; не ниже ~0.25 Tc.
        p.argument_limit_min = std::max(10.0, 0.25 * T_critical);
    }
    p.argument_limit_max = 5000;
    // correcting limits by estimation if applicable
    if (std::isfinite(initial_temperature)) {
        double enthalpy_at_initial_temperature = fluid->flash(given_pressure, initial_temperature).td_functions.enthalpy.mass.mix;
        if (enthalpy_at_initial_temperature < given_enthalpy)p.argument_limit_min = initial_temperature;
        if (enthalpy_at_initial_temperature > given_enthalpy)p.argument_limit_max = initial_temperature;
    }
    {
        double enthalpy_at_critical_temperature = fluid->flash(given_pressure, T_critical).td_functions.enthalpy.mass.mix;
        if (enthalpy_at_critical_temperature < given_enthalpy)p.argument_limit_min = T_critical;
        if (enthalpy_at_critical_temperature > given_enthalpy)p.argument_limit_max = T_critical;
    }
    p.argument_precision = desired_precision;
    p.residual_precision = desired_precision;
    p.argument_history = true;
    p.residual_history = true;
    p.solution_type = fixed_bisectional_solution_type::Combined;
    p.secant_treshhold_iterations = 0;
    //p.secant_treshhold_max = 10;
    //p.secant_treshhold_min = 0.001;
    p.verbose = false;
    fixed_bisectional<1>::solve(p, *this, numerical_result, analysis);

    return numerical_result->argument;
}

double ph_flash_bisection::solve(
    fixed_bisection_result_t<1>* numerical_result,
    fixed_bisection_result_analysis_t<1>* analysis)
{
    return solve(std::numeric_limits<double>::epsilon() * 10000., numerical_result, analysis);
}

ph_flash_stub_data_t ph_flash_bisection::get_stub_data() const
{
    ph_flash_stub_data_t result;
    result.pressure = given_pressure;
    result.enthalpy_mass = given_enthalpy;
    result.fluid = fluid->get_mock_data();
    return result;
}

ph_flash_newton::ph_flash_newton(const fluid_t* fluid, double pressure, double target_enthalpy_mass,
    double _temperature_initial, double _prec)
    : fluid(fluid)
    , target_enthalpy(target_enthalpy_mass)
    , pressure(pressure)
    , temperature_initial(_temperature_initial)
    , desired_precission(_prec)
{
    if (!std::isfinite(temperature_initial))
    {
        auto fluid_rd = dynamic_cast<const fluid_rault_dalton_t*>(fluid);
        if (fluid_rd == nullptr)
            throw std::runtime_error("Can estimate min temperature only for Raoult-fluid Dalton");
        double T_critical = fluid->get_pseudocritical_temperature();
        double T_initial = T_critical;
        double T_min = get_min_antoine_bound(fluid_rd);
        temperature_initial = std::max(T_min, 0.5 * T_initial);

    }
}

ph_flash_newton::ph_flash_newton(const fluid_t* fluid, double p,
    double initial_entalpy, double Q, double mass_flow,
    double temperature_initial)
    : ph_flash_newton(fluid, p, initial_entalpy + Q / mass_flow, temperature_initial)
{
}

double ph_flash_newton::residuals(const double& temperature)
{
    const auto vle = fluid->flash(pressure, temperature);
    double r = vle.td_functions.enthalpy.mass.mix - target_enthalpy;
    return r;
}

double ph_flash_newton::jacobian_dense(const double& temperature)
{
    double e = epsilon * std::max(1.0, abs(temperature));
    const auto vle = fluid->flash(pressure, temperature);//= fluid->get_last_flash_result();//спорно
    double mass_enthalpy_0 = vle.td_functions.enthalpy.mass.mix;

    const auto vle2 = fluid->flash(pressure, temperature + e);
    double mass_enthalpy_eps = vle2.td_functions.enthalpy.mass.mix;

    function_type J = (mass_enthalpy_eps - mass_enthalpy_0) / e;
    return J;
}

ph_flash_newton::var_type ph_flash_newton::estimation() const
{
    return temperature_initial;
}

double ph_flash_newton::solve(fixed_solver_result_t<1>* numerical_result,
    fixed_solver_result_analysis_t<1>* analysis_result)
{
    fixed_solver_result_t<1> solver_result_carrier;
    if (numerical_result == nullptr)
        numerical_result = &solver_result_carrier;


    double T_initial = estimation();
    fixed_solver_parameters_t<1, 0, golden_section_search> solver_parameters;
    solver_parameters.line_search.function_decrement_factor = 10;
    solver_parameters.line_search.iteration_count = 30;
    solver_parameters.line_search_fail_action = line_search_fail_action_t::TreatAsFail;

    solver_parameters.constraints.minimum = temperature_minimum;
    solver_parameters.constraints.maximum = temperature_maximum;

    solver_parameters.argument_increment_norm = desired_precission;
    if (analysis_result != nullptr) {
        solver_parameters.analysis.argument_history = true;
        solver_parameters.analysis.steps = true;
        solver_parameters.analysis.line_search_explore = true;
    }

    fixed_newton_raphson<1>::solve_dense(
        *this, T_initial, solver_parameters, numerical_result, analysis_result);

    if (numerical_result->result_code == numerical_result_code_t::Converged) {
        return numerical_result->argument;
    }
    else {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

ph_flash_stub_data_t ph_flash_newton::get_stub_data() const
{
    ph_flash_stub_data_t result;
    result.pressure = pressure;
    result.enthalpy_mass = target_enthalpy;
    result.fluid = fluid->get_mock_data();
    return result;
}

}
