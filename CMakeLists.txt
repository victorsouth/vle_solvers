cmake_minimum_required(VERSION 3.10)
if(POLICY CMP0091)
    cmake_policy(SET CMP0091 NEW)
endif()

project(vle_solvers)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
message("CXX standard: " ${CMAKE_CXX_STANDARD})

if(MSVC)
    cmake_policy(SET CMP0167 NEW)
    add_compile_options("/MP")     # включает многопроцессорную компиляцию
    add_compile_options("/bigobj")
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
if(MINGW)
    cmake_policy(SET CMP0167 NEW)
    add_compile_options("-Wa,-mbig-obj")
endif()

if(NOT TARGET fixed_solvers::fixed_solvers)
find_package(fixed_solvers REQUIRED)
endif()
find_package(Boost REQUIRED COMPONENTS serialization)
option(VLESOLVERS_WITH_PCH "Use Prcompiled Headers" ON)
# Установка путей к файлам в папках в переменные
set(HEADER
    vle/vle_solvers.h
    vle/components_db.h
    vle/physical_constants.h
)

set(HEADER_DB
    vle/db/db_bip.h
    vle/db/antoine_model.h
    vle/db/component_properties.h
    vle/db/fluid_stubdata.h
    vle/db/thermo_db.h
    vle/db/td_functions.h
    vle/db/tdb_serialize_common.h
    vle/db/tdb_serialize_2026_02_02.h
)

set(SOURCE_DB
    vle/db/antoine_model.cpp
    vle/db/component_properties.cpp
    vle/db/db_bip.cpp
    vle/db/components_db_data.cpp
	vle/db/components_db_data_BIP.cpp
    vle/db/thermo_db.cpp
    vle/db/td_functions.cpp
)

set(HEADER_FLUID
    vle/fluid/fluid_base.h
    vle/fluid/fluid_common.h
    vle/fluid/fluid_equations.h
    vle/fluid/fluid_raoult_dalton.h
)

set(SOURCE_FLUID
    vle/fluid/fluid_base.cpp
    vle/fluid/fluid_common.cpp
    vle/fluid/fluid_equations.cpp
    vle/fluid/fluid_raoult_dalton.cpp
)

set(HEADER_FLUID_SOLVERS
    vle/fluid/solvers/ph_flash.h
    vle/fluid/solvers/rachford_rice.h
    vle/fluid/solvers/uv_flash.h
    vle/fluid/solvers/phase_boundaries.h
)

set(SOURCE_FLUID_SOLVERS
    vle/fluid/solvers/ph_flash.cpp
    vle/fluid/solvers/rachford_rice.cpp
    vle/fluid/solvers/uv_flash.cpp
)

set(HEADER_HELPERS
    vle/helpers/physical_helpers.h
    vle/helpers/string_helpers.h
)

# Организация папок с исходниками в сгенерированном проекте VS
source_group("Header Files" FILES ${HEADER})
source_group("Header Files/DB" FILES ${HEADER_DB})
source_group("Header Files/Fluid" FILES ${HEADER_FLUID})
source_group("Header Files/Fluid/Solvers" FILES ${HEADER_FLUID_SOLVERS})
source_group("Header Files/Helpers" FILES ${HEADER_HELPERS})

source_group("Source Files/DB" FILES ${SOURCE_DB})
source_group("Source Files/Fluid" FILES ${SOURCE_FLUID})
source_group("Source Files/Fluid/Solvers" FILES ${SOURCE_FLUID_SOLVERS})

# Создание библиотеки
add_library(${PROJECT_NAME} 
    STATIC
    ${HEADER}
    ${HEADER_DB}
    ${HEADER_FLUID}
    ${HEADER_FLUID_SOLVERS}
    ${HEADER_HELPERS}
    ${SOURCE_DB}
    ${SOURCE_FLUID}
    ${SOURCE_FLUID_SOLVERS}
)

# Линковка с требуемыми библиотеками
target_link_libraries(${PROJECT_NAME} 
    PUBLIC 
    fixed_solvers::fixed_solvers
    Boost::boost
    Boost::serialization
)
   
# include директории
target_include_directories(${PROJECT_NAME}
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vle>
    $<INSTALL_INTERFACE:include/vle_solvers>
    $<INSTALL_INTERFACE:include/vle_solvers/vle>
)
add_library(vle_solvers::vle_solvers ALIAS vle_solvers)
if(VLESOLVERS_WITH_PCH)
    target_precompile_headers(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:vle/vle_solvers.h> )
    target_compile_definitions(${PROJECT_NAME} PUBLIC VLESOLVERS_WITH_PCH)
endif()

option(VLE_SOLVERS_INSTALL "" ON)

# установка
if(VLE_SOLVERS_INSTALL)
    include(CMakePackageConfigHelpers)
    include(GNUInstallDirs)
    install(TARGETS vle_solvers
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        EXPORT vle_solversTargets
    )
    install(EXPORT vle_solversTargets NAMESPACE vle_solvers::
    FILE vle_solversTargets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/vle_solvers
    )

    # generate the config file that is includes the exports
    configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/vle_solversConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/vle_solvers"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )
    install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/vle_solversConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/vle_solvers
    )

    # Установка заголовочных файлов (headers)
    install(FILES ${HEADER} DESTINATION include/vle_solvers/vle/)
    install(FILES ${HEADER_DB} DESTINATION include/vle_solvers/vle/db)
    install(FILES ${HEADER_FLUID} DESTINATION include/vle_solvers/vle/fluid)
    install(FILES ${HEADER_FLUID_SOLVERS} DESTINATION include/vle_solvers/vle/fluid/solvers)
    install(FILES ${HEADER_HELPERS} DESTINATION include/vle_solvers/vle/helpers)

endif()

option(VLE_SOLVERS_BUILD_TESTS "" ON)
# Сборка тестов
if(VLE_SOLVERS_BUILD_TESTS)

    find_package(Threads)
    find_package(GTest REQUIRED)

    set(TESTING_SOURCE
        testing/testing_main.cpp
    )

    set(TESTING_HEADERS
        testing/test_database.h
        testing/test_fluid.h
        testing/test_components.h
        testing/test_bips.h
        testing/test_bips_verification.h
    )

    source_group("Header Files" FILES ${TESTING_HEADERS})
    source_group("Source Files" FILES ${TESTING_SOURCE})

    add_executable(testing_vle_solvers 
        ${TESTING_SOURCE}
        ${TESTING_HEADERS}
    )
    target_link_libraries(testing_vle_solvers vle_solvers::vle_solvers GTest::gtest)
    target_include_directories(testing_vle_solvers
        PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/testing>
    )
    #target_precompile_headers(testing_vle_solvers REUSE_FROM ${PROJECT_NAME})
    install(TARGETS testing_vle_solvers DESTINATION bin)

endif()


