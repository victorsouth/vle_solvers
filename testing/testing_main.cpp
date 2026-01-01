#include <vle/vle_solvers.h>

#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include <gtest/gtest.h>

#include "test_database.h"

#include <iostream>


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
#if defined(_WIN32) && !defined(__MINGW32__)
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
    int res = RUN_ALL_TESTS();
    return res;
}
