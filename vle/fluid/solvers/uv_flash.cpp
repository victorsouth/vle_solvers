//#include "../../vle_solvers.h" // очень некрасивый инклуд...
#include <cmath>
#include <Eigen/Dense>
#include <fixed/helpers/math_helpers.h>
#include "uv_flash.h"
#include "rachford_rice.h"

namespace vlelib {
;

double uv_flash_over_RR::ensure_two_phase_pressure(const fluid_rault_dalton_t* fluid, 
    double pressure, double temperature)
{
    double P_bubble = fluid->get_bubble_point_at_given_temperature(temperature);
    double P_dew = fluid->get_dew_point_at_given_temperature(temperature);

    double result = std::max(P_dew, std::min(P_bubble, pressure));

    //return 0.5 * (P_dew + P_bubble);
    ///return P_bubble;
    return result;
}

std::pair<double, double> uv_flash_over_RR::norm(double P, double T) const
{
    double P_bubble = fluid->get_bubble_point_at_given_temperature(T);
    double P_dew = fluid->get_dew_point_at_given_temperature(T);

    double phi = (P - P_dew) / (P_bubble - P_dew);
    double theta = T / T_critical;

    return std::make_pair(phi, theta);
}

std::pair<double, double> uv_flash_over_RR::norm(const fluid_rault_dalton_t* fluid, double P, double T)
{
    double T_critical = fluid->get_pseudocritical_temperature();
    double P_bubble = fluid->get_bubble_point_at_given_temperature(T);
    double P_dew = fluid->get_dew_point_at_given_temperature(T);

    double phi = (P - P_dew) / (P_bubble - P_dew);
    double theta = T / T_critical;

    return std::make_pair(phi, theta);
}


std::pair<double, double> uv_flash_over_RR::denorm(double phi, double theta) const
{
    double T = theta * T_critical;
    double P_bubble = fluid->get_bubble_point_at_given_temperature(T);
    double P_dew = fluid->get_dew_point_at_given_temperature(T);
    double P = (P_bubble - P_dew) * phi + P_dew;
    return std::make_pair(P, T);
}

std::pair<double, double> uv_flash_over_RR::denorm(const fluid_rault_dalton_t* fluid, double phi, double theta)
{
    double T_critical = fluid->get_pseudocritical_temperature();
    double T = theta * T_critical;
    double P_bubble = fluid->get_bubble_point_at_given_temperature(T);
    double P_dew = fluid->get_dew_point_at_given_temperature(T);
    double P = (P_bubble - P_dew) * phi + P_dew;
    return std::make_pair(P, T);
}

double uv_flash_over_RR::get_min_theta() const
{
    return T_min / T_critical;
}


uv_flash_over_RR::uv_flash_over_RR(double volume, double molar_amount, 
    double _inner_energy_molar, const fluid_rault_dalton_t* _fluid, 
    bool _use_density, double _initial_pressure /*= std::numeric_limits<double>::quiet_NaN()*/, 
    double _initial_temperature /*= std::numeric_limits<double>::quiet_NaN()*/) 
    : fluid(_fluid)
    , inner_energy_molar(_inner_energy_molar)
    , initial_pressure(_initial_pressure)
    , initial_temperature(_initial_temperature)
    , use_density(_use_density)
{
    T_critical = fluid->get_pseudocritical_temperature();
    T_min = get_min_antoine_bound(fluid);
    specific_volume_molar = volume / molar_amount;
    epsilon = 1e-6;

    if (std::isfinite(initial_pressure) && std::isfinite(initial_temperature)) {
        // сначала поправим температуру, затем для нее при необходимости загоним давление в двухфазную область
        initial_temperature = std::max(initial_temperature, T_min);
        // убеждаемся, что давление лежит в двухфазной области при заданной температуре
        initial_pressure = ensure_two_phase_pressure(
            fluid, initial_pressure, initial_temperature);
    }
    else {
        std::tie(initial_pressure, initial_temperature) = denorm(0.5, 0.5);
        initial_temperature = std::max(initial_temperature, T_min);
    }
}


/// @brief Определение доли отгона
/// @param fluid 
/// @param P 
/// @param T 
/// @return 
double calc_omega_extrapolation(const fluid_t* fluid, double P, double T)
{
    auto calc_omega_simple = [&](double P, double T) {
        rachford_rice2_t rr(fluid, P, T);
        auto rr_solver_result = rr.solve();

        double omega = rr_solver_result.argument;
        return omega;
        };

    double P_bubble = fluid->get_bubble_point_at_given_temperature(T);
    double P_dew = fluid->get_dew_point_at_given_temperature(T);

    if (P > P_bubble) {
        double eps = 1e-4;
        double dp = P_bubble * eps;
        double omega_plus = calc_omega_simple(P_bubble + dp, T);
        double omega_minus = calc_omega_simple(P_bubble - dp, T);
        double domega_dp = (omega_plus - omega_minus) / (2 * eps);
        double omega = sqrt(sqrt(P - P_bubble)) * domega_dp;
        return omega;
    }
    //if (P < P_dew) {
    //    double eps = 1e-4;
    //    double dp = P_bubble * eps;
    //    double omega_plus = calc_omega_simple(P_dew + dp, T);
    //    double omega_minus = calc_omega_simple(P_dew - dp, T);
    //    double domega_dp = (omega_plus - omega_minus) / (2 * eps);
    //    double omega = 1 + (P - P_dew) * domega_dp;
    //    return omega;
    //}

    return calc_omega_simple(P, T);
}


uv_flash_over_RR::var_type uv_flash_over_RR::residuals(const var_type& w)
{
    auto [P, T] = denorm(w[0], w[1]);

    double omega = calc_omega_extrapolation(fluid, P, T);
    //omega = std::max(0.0, std::min(1.0, omega)); // убрал 16.04.2025, без нее хуже сходимость, ц.ф. иногда не может в нуль обратиться

    rachford_rice2_t rr(fluid, P, T);
    auto rr_result = rr.build_result(omega);

    flash_calculation_result_t flash;

    fluid->build_twophase_result2(P, T, rr_result, flash);

    var_type r;

    if (use_density) {
        double density_target = fluid->get_molar_mass() / specific_volume_molar;
        r[0] = density_target - flash.density.mix;
        r[0] /= density_target;
    }
    else {
        double specific_volume_molar_calc = fluid->get_molar_mass() / flash.density.mix;
        r[0] = specific_volume_molar - specific_volume_molar_calc;
        r[0] /= specific_volume_molar;
    }

    r[1] = inner_energy_molar - flash.inner_energy.molar.mix;
    r[1] /= inner_energy_molar;

    return r;
}

uv_flash_over_RR::matrix_value uv_flash_over_RR::jacobian_dense(const var_type& x)
{
    var_type arg = x;

    matrix_value J;

    for (int component = 0; component < x.size(); ++component) {
        double e = epsilon * std::max(1.0, abs(arg[component]));

        double x_plus = x[component] + e;
        double x_minus = x[component] - e;

        if (component == 0) { // "давление" phi
            x_plus = std::min(1.0, x_plus);
            x_minus = std::max(0.0, x_minus);
        }
        // Аналогичного ограничения на температуру не делается, т.к.
        // при небольшом выходе за нижний диапазон Антуана ничего не случится
        // при выходе в закритическую температуру для Рауля тоже ничего не случится

        arg[component] = x_plus;
        function_type f_plus = residuals(arg);
        arg[component] = x_minus;
        function_type f_minus = residuals(arg);
        arg[component] = x[component];

        function_type Jcol = (f_plus - f_minus) / (x_plus - x_minus);
        for (size_t row = 0; row < static_cast<size_t>(x.size()); ++row) {
            J[row][component] = Jcol[row];
        }
    }
    return J;
}

bool uv_flash_over_RR::custom_success_criteria(const var_type& r, const var_type& x)
{
    //return false;
    constexpr double residual_border = 1e-5;
    for (double residual : r) {
        if (std::abs(residual) > residual_border) {
            return false;
        }
    }

    return true;
}

uv_flash_over_RR::var_type uv_flash_over_RR::estimation() const
{
    auto [phi, theta] = norm(initial_pressure, initial_temperature);
    var_type x0{ phi, theta };
    return x0;
}

fixed_solver_parameters_t<2, 0, golden_section_search> uv_flash_over_RR::prepare_solver_parameters(
    bool perform_analysis) const
{
    fixed_solver_parameters_t<2, 0, golden_section_search> solver_parameters;
    solver_parameters.argument_increment_norm = 1e-7;
    solver_parameters.residuals_norm = 1e-4; // исследовать невязки и раскомментить

    solver_parameters.step_constraint_as_optimization = false;
    solver_parameters.residuals_norm_allow_force_success = false;
    solver_parameters.line_search.iteration_count = 10;
    solver_parameters.line_search.function_decrement_factor = 10;
    solver_parameters.line_search_fail_action = line_search_fail_action_t::TreatAsFail;

    // Ограничения давления
    solver_parameters.constraints.minimum[0] = 0.0;
    solver_parameters.constraints.maximum[0] = 1.0;
    // Ограничения приведенной температуры
    solver_parameters.constraints.minimum[1] = get_min_theta();
    solver_parameters.constraints.maximum[1] = 2.0;


    if (perform_analysis) {
        solver_parameters.analysis.steps = true;
        solver_parameters.analysis.line_search_explore = true;
        solver_parameters.analysis.argument_history = true;
    }

    return solver_parameters;

}

fixed_solver_result_t<2> uv_flash_over_RR::run_solution_attemps(
    const std::array<uv_flash_solution_attempt_options, 4>& options, 
    const fixed_solver_parameters_t<2, 0, golden_section_search>& solver_parameters,
    fixed_solver_result_analysis_t<2>* analysis_result)
{
    fixed_solver_result_t<2> solver_result;

    // отдавать или не отдавать исключения от промежуточных попыток 
    // (или кидать )
    bool allow_intermediate_exceptions = false; 

    for (size_t index = 0; index < options.size(); ++index)
    {
        try
        {
            fixed_newton_raphson<2>::solve_dense(
                *this, options[index].initial_estimation, solver_parameters,
                &solver_result, analysis_result);
            if (solver_result.result_code == numerical_result_code_t::Converged) {
                return solver_result;
            }
        }
        catch (...)
        {
            // отдаем только последнее исключение
            bool is_final_attempt = index == options.size() - 1;
            if (is_final_attempt || allow_intermediate_exceptions)
                throw std::current_exception();
        }
    }
    return solver_result;
}


fixed_solver_result_t<2> uv_flash_over_RR::solve(
    fixed_solver_result_analysis_t<2>* analysis_result /*= nullptr*/)
{
    fixed_solver_parameters_t<2, 0, golden_section_search>
        solver_parameters = prepare_solver_parameters(analysis_result != nullptr);


    var_type x0 = estimation();
    var_type x1 = { 0.5, 1.0 };
    var_type x2 = { 0.0, std::max(0.5, get_min_theta()) };
    // на случай, если x0 взято из предшествовавшего расчета
    var_type x3 = { 0.5, std::max(0.5, get_min_theta()) }; 
    
    std::array <uv_flash_solution_attempt_options, 4> options{
        uv_flash_solution_attempt_options(x0, false),
        uv_flash_solution_attempt_options(x1, false),
        uv_flash_solution_attempt_options(x2, false),
        uv_flash_solution_attempt_options(x3, false),
    };

    return run_solution_attemps(options, solver_parameters, analysis_result);
}

const uv_flash_over_RR_stubdata_t uv_flash_over_RR::get_mock_data() const
{
    uv_flash_over_RR_stubdata_t data;
    data.fluid = fluid->get_mock_data();
    data.initial_pressure = initial_pressure;
    data.initial_temperature = initial_temperature;
    data.specific_volume_molar = specific_volume_molar;
    data.inner_energy_molar = inner_energy_molar;
    data.use_density = use_density;
    return data;
}

}
