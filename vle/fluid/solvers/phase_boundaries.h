#pragma once

namespace vlelib { 
;



/// @brief Система уравнений для расчета температуры конденсации при заданном давлении
class dew_point_at_given_pressure : public fixed_system_t<1>
{
private:
    /// @brief Давление, при котором ищется температура кипения
    double pressure;
    /// @brief Флюид, для которого ищется температура кипения
    const fluid_t* fluid;
public:
    dew_point_at_given_pressure(double pressure, const fluid_t* fluid) 
        : pressure(pressure)
        , fluid(fluid)
    {

    }

    /// @brief Функция невязок
    /// @param temperature Искомая температура
    virtual double residuals(const double& temperature) override
    {
        return fluid->vapor_only_criteria(pressure, temperature) - 1;
    }

    /// @brief Функция расчета температуры конденсации
    /// @param solver_result Результат расчета
    /// @param analysis_result Аналитика расчета
    /// @return Если НР сошелся, то температуру.
    double solve(fixed_solver_result_t<1>* solver_result, fixed_solver_result_analysis_t<1>* analysis_result = nullptr)
    {
        double initial_estimation = 273.15;
        fixed_solver_parameters_t<1, 0> solver_parameters;
        solver_parameters.constraints.maximum = 2 * fluid->get_pseudocritical_temperature();
        solver_parameters.constraints.minimum = fluid->get_min_antoine_bound();
        solver_parameters.constraints.ensure_constraints(initial_estimation);
        solver_parameters.argument_increment_norm = 1e-7;

        if (analysis_result != nullptr) {
            solver_parameters.analysis.argument_history = true;
            solver_parameters.analysis.steps = true;
        }

        try
        {
            fixed_newton_raphson<1>::solve_dense(
                *this, initial_estimation, solver_parameters, solver_result, analysis_result);
            if (solver_result->result_code == numerical_result_code_t::Converged) {
                return solver_result->argument;
            }
            else
            {
                throw std::runtime_error("Dew point calculation NOT converged");
            }
        }
        catch (...)
        {
            throw std::runtime_error("Dew point calculation NOT converged");
        }
    }

};

/// @brief Система уравнений для расчета температуры кипения при заданном давлении
class bubble_point_at_given_pressure : public fixed_system_t<1>
{
private:
    /// @brief Давление, при котором ищется температура кипения
    double pressure;
    /// @brief Флюид, для которого ищется температура кипения
    const fluid_t* fluid;
public:
    bubble_point_at_given_pressure(double pressure, const fluid_t* fluid)
        : pressure(pressure)
        , fluid(fluid)
    {
    }
    /// @brief Функция невязок
    /// @param temperature Искомая температура
    virtual double residuals(const double& temperature) override {
        return fluid->liquid_only_criteria(pressure, temperature) - 1;
    }
    /// @brief Функция расчета температуры кипения
    /// @param solver_result Результат расчета
    /// @param analysis_result Аналитика расчета
    /// @return Если НР сошелся, то температуру.
    double solve(fixed_solver_result_t<1>* solver_result,
        fixed_solver_result_analysis_t<1>* analysis_result = nullptr)
    {
        double initial_estimation = 273.15;
        fixed_solver_parameters_t<1, 0> solver_parameters;
        solver_parameters.constraints.maximum = 2 * fluid->get_pseudocritical_temperature();
        solver_parameters.constraints.minimum = fluid->get_min_antoine_bound();
        solver_parameters.constraints.ensure_constraints(initial_estimation);
        solver_parameters.argument_increment_norm = 1e-7;

        if (analysis_result != nullptr) {
            solver_parameters.analysis.argument_history = true;
            solver_parameters.analysis.steps = true;
        }

        try
        {
            fixed_newton_raphson<1>::solve_dense(
                *this, initial_estimation, solver_parameters, solver_result, analysis_result);
            if (solver_result->result_code == numerical_result_code_t::Converged) {
                return solver_result->argument;
            }
            else
            {
                throw std::runtime_error("Bubble point calculation NOT converged_1");
            }
        }
        catch (...)
        {
            throw std::runtime_error("Bubble point calculation NOT converged_1");
        }
    }
};



}
