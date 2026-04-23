#include <vle/vle_solvers.h>

#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include <gtest/gtest.h>

#include <iostream>
#include <map>

// Вспомогательная функция для вывода ошибок
inline void print_errors(const std::wstring& prefix, double a, double b, double error_border) {
    double relative_error = (a - b) / a;
    std::wcout << prefix << a << '\t' << b << " err:" << relative_error << std::endl;

    ASSERT_LE(std::abs(relative_error), error_border);
}

#include "test_database.h"
#include "test_components.h"
#include "test_bips.h"
#include "test_bips_verification.h"
#include "test_fluid.h"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
#if defined(_WIN32) && !defined(__MINGW32__)
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
    int res = RUN_ALL_TESTS();
    return res;
}
