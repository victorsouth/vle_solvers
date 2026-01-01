#include "../../vle_solvers.h" // очень некрасивый инклуд...

namespace vlelib {
;


double rr_derivative(double split, const Eigen::VectorXd& molar_fraction, const Eigen::VectorXd& K_values)
{
    double result = 0;
    for (int index = 0; index < molar_fraction.size(); ++index) {
        double Ki = K_values(index);
        double zi = molar_fraction(index);
        if (zi != 0) {
            double numerator = -zi * fixed_solvers::sqr(Ki - 1.0);

            // При подходе к разрыву, отвечающему текущей Ki 
            // функция терпит разрыв, производная имеет гигантское значение
            // ограничиваем производную, ограничив близость нуля для 
            // выражения под знаком квадрата в знаменателе 
            // (см. формулу производной по ссылке)
            double denumerator_sqrt = 1.0 + split * (Ki - 1.0);
            constexpr double eps = 1e-8;
            if (std::abs(denumerator_sqrt) < eps) {
                denumerator_sqrt = eps * fixed_solvers::pseudo_sgn(denumerator_sqrt);
            }
            
            double denumerator = fixed_solvers::sqr(denumerator_sqrt);
            

            double derivative =  numerator / denumerator;
            result += derivative;
        }
    }
    return result;
}


rachford_rice2_t::rachford_rice2_t(const fluid_t* fluid, double pressure,
    double temperature, double vapor_fraction_initial)
    : fluid(fluid)
    , K_values(fluid->get_K_values(pressure, temperature))
    , pressure(pressure)
    , temperature(temperature)
    , vapor_fraction_initial(vapor_fraction_initial)
{

}

std::pair<double, double> rachford_rice2_t::get_k_boundaries() const
{
    double minimum = std::numeric_limits<double>::max();
    double maximum = std::numeric_limits<double>::min();

    // Идем по каждому компоненту в смеси
    for (size_t i = 0; i < fluid->get_components_count(); ++i)
    {
        if (fluid->get_molar_fraction()[i] != 0.0)
        {
            // Записываем минимальное значение и обновляем, если найдем меньше
            if (K_values[i] < minimum)
            {
                minimum = K_values[i];
            }

            // Записываем максимальное значение и обновляем, если найдем больше
            if(K_values[i] > maximum)
            {
                maximum = K_values[i];
            }

        }
    }

    if (minimum == std::numeric_limits<double>::max() || maximum == std::numeric_limits<double>::min())
    {
        throw std::runtime_error("Сoncentrations are zero. Minimum or Maximum k_value not found");
    }

    return std::make_pair(minimum, maximum);
}

std::pair<double, double> rachford_rice2_t::get_omega_boundaries() const {
    auto [K_min, K_max] = get_k_boundaries();

    if (K_min < 1 && K_max > 1) {
        // Полюса Речфорда по обе стороны от нуля
        // есть ограничение и сверху, и снизу
        double omega_max = 1 / (1 - K_min) - std::numeric_limits<double>::epsilon();
        double omega_min = 1 / (1 - K_max) + std::numeric_limits<double>::epsilon();
        return std::make_pair(omega_min, omega_max);
        /*solver_parameters.constraints.maximum = omega_max;
        solver_parameters.constraints.minimum = omega_min;
        solver_parameters.constraints.ensure_constraints(initial_value);*/
    }
    else {
        return std::make_pair(
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN());
    }
    //return std::make_pair(omega_min, omega_max);
}

double rachford_rice2_t::residuals(const double& split)
{
    const auto& molar_fraction = fluid->get_molar_fraction();
    return rr_equation(split, molar_fraction, K_values);
}

double rachford_rice2_t::jacobian_dense(const double& split)
{
    const auto& molar_fraction = fluid->get_molar_fraction();
    return rr_derivative(split, molar_fraction, K_values);
}

rachford_rice_result_t rachford_rice2_t::build_result(double vapor_split) const
{
    const Eigen::VectorXd& molar_fraction = fluid->get_molar_fraction();

    rachford_rice_result_t result;
    result.vapor_split = vapor_split;

    // Чтобы при negative flash были x, y, суммарно равные 1:
    double omega_bound = std::min(1.0, std::max(vapor_split, 0.0));
    double pressure_dew = fluid->get_dew_point_at_given_temperature(temperature);
    double pressure_bubble = fluid->get_bubble_point_at_given_temperature(temperature);

    if (vapor_split >= 1) {
        Eigen::VectorXd K_values_dew_border = fluid->get_K_values(pressure_dew, temperature);
        result.y = molar_fraction;
        result.x.resize(molar_fraction.size());
        for (size_t index = 0; index < static_cast<size_t>(molar_fraction.size()); ++index) {
            // значение x[i] может повредиться, если z[i] = 0 и соответственно K[i] = 0
            if (molar_fraction(index) != 0) {
                result.x(index) = result.y(index) / K_values_dew_border(index);
            }
            else {
                result.x(index) = 0;
            }
        }
    }
    else if (vapor_split <= 0) {
        Eigen::VectorXd K_values_bubble_border = fluid->get_K_values(pressure_bubble, temperature);
        result.x = molar_fraction;
        result.y = result.x.cwiseProduct(K_values_bubble_border);
    }
    else {
        std::tie(result.x, result.y)
            = rr_calc_compositions(omega_bound, molar_fraction, K_values);
    }
    return result;
}

fixed_solver_result_t<1> rachford_rice2_t::solve(fixed_solver_result_analysis_t<1>* solver_analysis)
{
    fixed_solver_parameters_t<1, 0, golden_section_search> solver_parameters;
    solver_parameters.constraints.relative_boundary = 1.0;
    solver_parameters.argument_increment_norm = 1e-10;
    solver_parameters.line_search_fail_action = line_search_fail_action_t::TreatAsFail;
    solver_parameters.line_search.iteration_count = 100;

    double initial_value = 0.5;
    if (std::isfinite(vapor_fraction_initial))
        initial_value = vapor_fraction_initial;

    std::tie(
        solver_parameters.constraints.minimum, 
        solver_parameters.constraints.maximum
    ) = get_omega_boundaries();

    solver_parameters.constraints.ensure_constraints(initial_value);

    fixed_solver_result_t<1> solver_result;
    if (solver_analysis) {
        solver_parameters.analysis.argument_history = true;
        solver_parameters.analysis.line_search_explore = true;
        solver_parameters.analysis.steps = true;
    }
    fixed_newton_raphson<1>::solve_dense(*this, initial_value,
        solver_parameters, &solver_result, solver_analysis);

    return solver_result;
}


}
