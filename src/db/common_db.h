#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <locale>

namespace vle_solvers
{
;

#define KELVIN_OFFSET 273.15
#define M_R 8.31441
#define M_PI 3.14159265358979323846


inline double celcium2kelvin(double celcium_value) {
    return celcium_value + KELVIN_OFFSET;
}

inline double kelvin2celcium(double kelvin_value) {
    return kelvin_value - KELVIN_OFFSET;
}

template <typename Coefficients>
double polyval(const Coefficients& poly_coeffs, double x)
{
    // НЕ используется формат Matlab! 
    // индекс коэффициента является показателем степени

    const size_t n = poly_coeffs.size();

    if (n == 0)
        return 0;

    double result = 0;
    double power = 1;
    for (size_t index = 0; index < n; ++index) {
        result += power * poly_coeffs[index];
        power *= x;
    }
    return result;
};

template<class T>
inline std::string int2str(T i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

template <typename T> inline int pseudo_sgn(T val) {
    return 2 * (val >= T(0)) - 1;
}

inline void string_replace(std::string& str, const std::string& from, const std::string& to)
{
    if (from.empty())
        return;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

inline std::vector<double> solve_realpoly3_vieta(std::vector<double> poly_coeffs)
{
    using std::pow;
    // https://ru.wikipedia.org/wiki/Тригонометрическая_формула_Виета

    // 1. Делаем коэффициент при старшей степени единичным 
    // c + b*x + a*x^2 + x^3 = 0
    double higher_order_coeff = poly_coeffs.back();
    for (double& coeff : poly_coeffs) {
        coeff /= higher_order_coeff;
    }
    double a = poly_coeffs[2];
    double b = poly_coeffs[1];
    double c = poly_coeffs[0];

    // 2. Определение количества действительных корней
    double Q = (a * a - 3 * b) / 9;
    double R = (2 * pow(a, 3) - 9 * a * b + 27 * c) / 54;
    double S = pow(Q, 3) - R * R;

    // 3. Расчет корней
    if (fabs(S) < std::numeric_limits<double>::epsilon())//epsilon_S)
    {
        //Ю.П. 2022.05.26
        //if(0)std::cerr<<"two roots:"<<-2.*cbrt(R) - a/3.<<'\t'<<cbrt(R) - a/3.<<std::endl;
        return {
            -2. * cbrt(R) - a / 3.,
            cbrt(R) - a / 3.
        };
        //throw logic_error("3rd order polynom singular case (S = 0), not implemented");
    }
    else if (S > 0)
    {
        // проверено, работает 27.12.2018
        // 3 действительных корня
        double phi = 1.0 / 3.0 * acos(R / pow(Q, 1.5));
        return {
            -2 * sqrt(Q) * cos(phi) - a / 3,
            -2 * sqrt(Q) * cos(phi + 2.0 / 3.0 * M_PI) - a / 3,
            -2 * sqrt(Q) * cos(phi - 2.0 / 3.0 * M_PI) - a / 3,
        };
    }
    else // S < 0
    {
        // 1 действительный корень, 2 комплексных
        if (fabs(Q) < std::numeric_limits<double>::epsilon()) {
            //Ю.П. 2022.05.26
            return{ -cbrt(c - (a * a * a) / 27.) - a / 3. };
            //throw logic_error("3rd order polynom singular case (Q = 0), not implemented");
        }
        else if (Q > 0)
        {
            // проверено, работает 27.12.2018
            double phi = 1.0 / 3.0 * acosh(fabs(R) / pow(Q, 1.5)); // квадратный корень из Q в кубе!
            return {
                -2.0 * pseudo_sgn(R) * sqrt(Q) * cosh(phi) - a / 3
            };
        }
        else // Q < 0
        {
            // проверено, работает 6.02.2019
            double phi = 1.0 / 3.0 * asinh(fabs(R) / pow(fabs(Q), 1.5));
            return {
                -2.0 * pseudo_sgn(R) * sqrt(fabs(Q)) * sinh(phi) - a / 3
            };
        }
    }
}

inline std::string wide2string(const std::wstring& str)
{
    std::string ascii_string(str.size(), '\0');
    std::transform(str.begin(), str.end(), ascii_string.begin(),
        [](wchar_t s) { return static_cast<char>(s); });
    return ascii_string;
}

template <typename Coefficients>
std::vector<double> poly_integral_coefficients(const Coefficients& poly_coeffs)
{
    // НЕ используется формат Matlab! в отличие от функции 
    // индекс коэффициента является показателем степени

    const size_t n = poly_coeffs.size();
    std::vector<double> result(n + 1);

    //double result = 0;
    //double power = 1;

    for (size_t index = 1; index < n + 1; ++index)
    {
        result[index] = poly_coeffs[index - 1] / index;
    }
    return result;
}


inline std::wstring string2wide(const std::string& str)
{
#ifndef __GNUC__
    try {
        std::wstring wstr(str.size(), 0);
        std::use_facet<std::ctype<wchar_t> >(std::locale("rus_rus.1251")).widen
        (&str[0], &str[0] + str.size(), &wstr[0]);
        return wstr;
    }
    catch (...)
    {
        std::wstring wide_string(str.begin(), str.end());
        return wide_string;
    }
#else
    std::wstring wstr = UTF8_to_wchar(str.c_str());
    return wstr;
    //return wstring(str.begin(),str.end());
#endif
}

template <typename Coeffs>
struct function_range_t {
    double range_start;
    double range_end;
    Coeffs coefficients;
};

inline bool value_in_range(double value, double range_begin, double range_end)
{
    if (range_begin > range_end)
        return value_in_range(value, range_end, range_begin);

    return value >= range_begin && value <= range_end;
}

template <typename Coeffs>
class ranged_function_t {
protected:
    std::vector<function_range_t<Coeffs>> ranges;
public:
    ranged_function_t() = default;
    ranged_function_t& operator= (const ranged_function_t&) = default;

    ranged_function_t(const std::vector<function_range_t<Coeffs>>& _ranges)
        : ranges(_ranges)
    {
        std::sort(ranges.begin(), ranges.end(),
            [](const function_range_t<Coeffs>& r1, const function_range_t<Coeffs>& r2)
            {
                return r1.range_start < r2.range_start;
            });



    }
    size_t get_range_index(double x) const
    {
        // find_if ищет первый элемент, удовлетворяющий условию 
        // (проверил по хелпу https://en.cppreference.com/w/cpp/algorithm/find)
        auto range_coefficients = std::find_if(ranges.begin(), ranges.end(),
            [&](const function_range_t<Coeffs>& r) {
                return x < r.range_end;
            }
        );

        if (range_coefficients == ranges.end())
            throw std::logic_error("approximation range not found");

        if (x > range_coefficients->range_end)
            throw std::logic_error("approximation range not found");

        return range_coefficients - ranges.begin();
    }
    /// @brief Геттер для ranges
    const std::vector<function_range_t<Coeffs>>& get_ranges() const {
        return ranges;
    }
    /// @brief Возвращает полный диапазон - от начала младшего диапазона до конца старшего
    std::pair<double, double> get_whole_range() const {
        if (ranges.empty())
            throw std::runtime_error("cannot get whole range for empty range list");
        
        // пользуемся тем, что диапазоны отсортированы в конструкторе
        return std::make_pair(
            ranges.front().range_start,
            ranges.back().range_end
            );
    }
};


template <typename Coeffs>
class ranged_polynom_t : public ranged_function_t<Coeffs> {
public:
    using ranged_function_t<Coeffs>::ranges;
    using ranged_function_t<Coeffs>::get_range_index;
protected:
    double gain, offset;
    std::vector<std::vector<double>> polynom_integral;

    // граничные значения
    std::vector<std::pair<double, double>> boundary_values;
public:
    ranged_polynom_t() = default;
    ranged_polynom_t& operator= (const ranged_polynom_t&) = default;

    ranged_polynom_t(const std::vector<function_range_t<Coeffs>>& _ranges)
    {};

    ranged_polynom_t(double gain, double offset, const std::vector<function_range_t<Coeffs>>& _ranges,
        double dy_error = 1e-8)
        : ranged_function_t<Coeffs>(_ranges)
        , gain(gain)
        , offset(offset)
    {
        std::transform(ranges.begin(), ranges.end(), std::back_inserter(boundary_values),
            [&](const function_range_t<Coeffs>& range) {
                std::pair<double, double> result;
                result.first = get_polynom_value(range, range.range_start);
                result.second = get_polynom_value(range, range.range_end);
                return result;
            });

        // проверка отсутствия разрыва на границах
        for (size_t index = 0; index < boundary_values.size() - 1; ++index) {
            double y_prev = boundary_values[index].second;
            double y_next = boundary_values[index + 1].first;
            if (fabs(y_prev - y_next) > dy_error) {
                throw std::logic_error("dy between ranges is too large");
            }
        }
    };

    virtual double get_polynom_value(const function_range_t<Coeffs>& range, double x) const
    {
        double result = polyval(range.coefficients, x);
        result = result * gain + offset;
        return result;
    };

    virtual double get_polynom_value(double x) const
    {
        size_t range_index = get_range_index(x);
        const auto& range = ranges[range_index];
        return get_polynom_value(range, x);
    };

    virtual double get_polynom_value_integral(double x) const
    {
        if (polynom_integral.empty()) {
            auto& polynom_int = const_cast<std::vector<std::vector<double>>&>(polynom_integral);

            for (const auto& range : this->ranges) {
                auto c = poly_integral_coefficients(range.coefficients);
                polynom_int.emplace_back(c);
            }
        }

        size_t range_index = this->get_range_index(x);
        double result = polyval(polynom_integral[range_index], x);

        result = result * gain + offset * x;

        return result;
    };
    // поиск диапазон по значению функции
    size_t get_inv_range_index(double y) const
    {
        for (size_t index = 0; index < boundary_values.size(); ++index) {
            if (value_in_range(y, boundary_values[index].first, boundary_values[index].second)) {
                return index;
            }
        }
        throw std::logic_error("approximation range not found");
    };

    virtual double get_inv_polynom_value(double y) const
    {
        if (ranges.empty())
            throw std::logic_error("No polynom ranges defined");

        size_t range_index = get_inv_range_index(y);

        const auto& range = ranges[range_index];

        std::vector<double> roots;
        size_t polynom_order = range.coefficients.size() - 1;
        switch (polynom_order)
        {
        case 1:
        {
            const auto& c = range.coefficients;
            double result = (y - c[0]) / c[1];
            return result;  // поскольку корень единственный, можно дальше не фильтровать roots
        }
        case 2:
            throw std::logic_error("inv polynom order 2 not implemented");
            break;
        case 3:
        {
            auto equation = range.coefficients;
            equation[0] -= y; // слева полином с коэффициентами, а справа значение y
            roots = solve_realpoly3_vieta(equation);
            break;
        }
        default:
            throw std::logic_error("inv polynom higher order not implemented");
        }

        std::vector<double> root_selected;
        copy_if(roots.begin(), roots.end(), std::back_inserter(root_selected),
            [&](double x)
            {
                return x >= range.range_start && x <= range.range_end;
            });

        if (root_selected.size() != 1)
            throw std::logic_error("wrong root count");
        return root_selected.front();
    };

};


}