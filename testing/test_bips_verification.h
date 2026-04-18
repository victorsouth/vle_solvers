#pragma once

/// @brief Отображение CAS-номера в пару (формула, CAS-номер).
/// Используется как эталонная база для проверки корректности загрузки БД.
using CAS_formula_map_t = std::unordered_map<std::wstring, std::pair<std::wstring, std::wstring>>;

/// @brief Возвращает эталонную карту CAS->(формула, CAS),
/// используемую в тестах для проверки корректности базы компонентов.
const CAS_formula_map_t get_target_component_CAS_map();

/// @brief Проверяет, что количество компонентов в базе совпадает с эталонным.
TEST(BinaryCoefficients, ReturnsSameComponentCountAsReference)
{
    auto ref = get_target_component_CAS_map();
    const auto& db = components_database.get_component_casno_database();

    ASSERT_EQ(db.size(), ref.size());
}

/// @brief Проверяет, что каждый CAS-номер из эталонной базы присутствует в компонентной базе.
TEST(BinaryCoefficients, FindsAllComponentsByCAS)
{
    auto ref = get_target_component_CAS_map();
    const auto& db = components_database.get_component_casno_database();

    for (const auto& [cas, pair] : ref) {
        ASSERT_NE(db.find(cas), db.end())
            << "Component with CAS " << cas << " not found";
    }
}

/// @brief Проверяет, что каждая химическая формула соответствует ровно одному компоненту.
/// Ожидается исключение при попытке получить компонент по дублирующейся формуле.
TEST(BinaryCoefficients, ThrowsIfFormulaIsDuplicated)
{
    const auto& db = components_database.get_component_casno_database();

    for (const auto& [cas, comp] : db) {
        ASSERT_NO_THROW({
            components_database.get_component_by_formula(comp.name);
            }) << "Formula " << comp.name << " is duplicated";
    }
}

/// @brief Проверяет, что бинарная база взаимодействий содержит только известные CAS-номера.
TEST(BinaryCoefficients, BIPUsesOnlyKnownCASNumbers)
{
    auto ref = get_target_component_CAS_map();
    const auto& bip_db = components_database.get_bip_db();
    for (const auto& [pair_cas, _] : bip_db) {
        for (const auto& cas : pair_cas) {
            ASSERT_NE(ref.find(cas), ref.end())
                << "Unknown CAS " << cas;
        }
    }
}

const bip_records_t bip_records_verification_Nishiumi() {
    return {
              {L"CH4", L"CH4", 0.0000000}
            , { L"CH4", L"C2H6", 0.0052191 }
            , { L"CH4", L"C3H8", 0.0159089 }
            , { L"CH4", L"n_C4H10", 0.0260425 }
            , { L"C2H6", L"CH4", 0.0052191 }
            , { L"C2H6", L"C2H6", 0.0000000 }
            , { L"C2H6", L"C3H8", 0.0056181 }
            , { L"C2H6", L"n_C4H10", 0.0119787 }
            , { L"C3H8", L"CH4", 0.0159089 }
            , { L"C3H8", L"C2H6", 0.0056181 }
            , { L"C3H8", L"C3H8", 0.0000000 }
            , { L"C3H8", L"n_C4H10", 0.0032230 }
            , { L"n_C4H10", L"CH4", 0.0260425 }
            , { L"n_C4H10", L"C2H6", 0.0119787 }
            , { L"n_C4H10", L"C3H8", 0.0032230 }
            , { L"n_C4H10", L"n_C4H10", 0.0000000 }
    };
}

/// @brief Набор эталонных данных для верификации матрицы BIP.
struct bip_matrix_verification_t {
    /// @brief Формулы компонентов
    std::vector<std::wstring> components_formulas;
    /// @brief матрица бинарных коэффициентов
    Eigen::MatrixXd bip_matrix;

    /// @brief Возвращает эталонную матрицу BIP, сформированную по данным Нисиуми.
    static bip_matrix_verification_t get_Nishiumi() {
        const bip_records_t bip_records_verification = bip_records_verification_Nishiumi();

        std::vector<std::wstring> component_formulas;
        auto add_unique = [&](const std::wstring& component_formula) {
            if (std::find(component_formulas.begin(), component_formulas.end(), component_formula)
                == component_formulas.end()) {
                component_formulas.push_back(component_formula);
            }
            };

        for (const auto& rec : bip_records_verification) {
            add_unique(rec.cas1);
            add_unique(rec.cas2);
        }

        const size_t size = component_formulas.size();
        Eigen::MatrixXd bip_matrix = Eigen::MatrixXd::Zero(size, size);

        std::unordered_map<std::wstring, size_t> component_index;
        component_index.reserve(size);
        for (size_t index = 0; index < size; ++index) {
            component_index[component_formulas[index]] = index;
        }

        for (const auto& rec : bip_records_verification) {
            const size_t i = component_index.at(rec.cas1);
            const size_t j = component_index.at(rec.cas2);
            bip_matrix(i, j) = rec.bip_value;
        }

        return { component_formulas, bip_matrix };
    }

    /// @brief Возвращает эталонную матрицу BIP для корреляции Чуи–Праусница.
    static bip_matrix_verification_t get_ChuehPrausnitz() {
        bip_matrix_verification_t bip_matrix_verification_ChuehPrausnitz = {
            { L"CH4", L"C2H6", L"C3H8", L"n_C4H10" },
            (Eigen::MatrixXd(4, 4) <<
                0.0000000, 0.0068442, 0.0214428, 0.0367701,
                0.0068442, 0.0000000, 0.0041499, 0.0122406,
                0.0214428, 0.0041499, 0.0000000, 0.0021642,
                0.0367701, 0.0122406, 0.0021642, 0.0000000
            ).finished()
        };
        return bip_matrix_verification_ChuehPrausnitz;
    }

    /// @brief Возвращает эталонную матрицу BIP для корреляции Гао.
    static bip_matrix_verification_t get_Gao() {
        bip_matrix_verification_t bip_matrix_verification_Gao = {
            { L"CH4", L"C2H6", L"C3H8", L"n_C4H10" },
            (Eigen::MatrixXd(4, 4) <<
                0.0000000, 0.0084334, 0.0164518, 0.0238250,
                0.0084334, 0.0000000, 0.0014037, 0.0041805,
                0.0164518, 0.0014037, 0.0000000, 0.0007479,
                0.0238250, 0.0041805, 0.0007479, 0.0000000
            ).finished()
        };
        return bip_matrix_verification_Gao;
    }
};

/// @brief Проверяет, что бинарные коэффициенты по умолчанию совпадают с табличными
/// значениями Нисиуми.
TEST(ThermoDB, BIPDefaultMatchesNishiumi)
{
    const bip_records_t bip_records_verification = bip_records_verification_Nishiumi();

    for (auto [formula1, formula2, bip] : bip_records_verification) {
        if (formula1 == formula2) {
            continue;
        }
        ASSERT_NEAR(components_database.get_bip_pair_formula(formula1, formula2), bip, 1e-7);
    }
}

/// @brief Проверяет способность верифицировать данные по корреляции Чуи–Праусница на примере от АМ
TEST(BinaryCoefficients, DISABLED_DeniesRaoultDalton)
{
    //// // На будущее
    //// // fluid_rault_dalton_t должен кидать исключение
    //// //ASSERT_THROW(
    //// //    components_database.create_fluid<vlelib::fluid_rault_dalton_t>(
    //// //        components, plan),
    //// //    std::runtime_error
    //// //);
    ////
    //// // А вот peng_robinson должен работать
    //// auto fluid = components_database.create_fluid<vlelib::fluid_peng_robinson_t>(
    ////     components, plan);
}

/// @brief Проверяет способность верифицировать данные по корреляции Чуи–Праусница на примере от АМ
TEST(BinaryCoefficients, VerifiesChuehPrausnitz)
{
    // Arrange
    bip_matrix_verification_t bip_matrix_verification = bip_matrix_verification_t::get_ChuehPrausnitz();
    std::vector<std::wstring> components = bip_matrix_verification.components_formulas;
    bip_estimation_plan_t plan = generate_bip_estimation_plan_with_correlation(
        components, bip_correlation_t::ChuehPrausnitz);

    // Act
    auto fluid = components_database.create_fluid<vlelib::fluid_rault_dalton_t>(
        components, {}, plan);
    const Eigen::MatrixXd& bip_matrix_evaluated = fluid->get_binary_coeffs_ref();

    // Assert
    ASSERT_TRUE(bip_matrix_evaluated.isApprox(
        bip_matrix_verification.bip_matrix, 1e-6));
}


/// @brief Проверяет способность верифицировать данные по корреляции Гао на примере от АМ
TEST(BinaryCoefficients, VerifiesGao)
{
    // Arrange
    bip_matrix_verification_t bip_matrix_verification = bip_matrix_verification_t::get_Gao();
    std::vector<std::wstring> components = bip_matrix_verification.components_formulas;
    bip_estimation_plan_t plan = generate_bip_estimation_plan_with_correlation(
        components, bip_correlation_t::Gao);

    // Act
    auto fluid = components_database.create_fluid<vlelib::fluid_rault_dalton_t>(
        components, {}, plan);
    const Eigen::MatrixXd& bip_matrix_evaluated = fluid->get_binary_coeffs_ref();

    // Assert
    ASSERT_TRUE(bip_matrix_evaluated.isApprox(
        bip_matrix_verification.bip_matrix, 1e-6));
}


/// @brief Проверяет способность не испортить данные из корреляции Нисиуми
TEST(BinaryCoefficients, VerifiesNishiumi)
{
    // Arrange
    bip_matrix_verification_t bip_matrix_verification = bip_matrix_verification_t::get_Nishiumi();
    std::vector<std::wstring> components = bip_matrix_verification.components_formulas;
    bip_estimation_plan_t plan = generate_bip_estimation_plan_use_db_only(components);

    // Act
    auto fluid = components_database.create_fluid<vlelib::fluid_rault_dalton_t>(
        components, {}, plan);
    const Eigen::MatrixXd& bip_matrix_evaluated = fluid->get_binary_coeffs_ref();

    // Assert
    ASSERT_TRUE(bip_matrix_evaluated.isApprox(
        bip_matrix_verification.bip_matrix, 1e-6));
}

const CAS_formula_map_t get_target_component_CAS_map()
{
    std::unordered_map<std::wstring, pair< std::wstring, std::wstring> > target_component_CAS_map =
    {
        //  \"AR\", \"argon\"
        //  \"C2H2\", \"acetylene\"
        //  \"C2H4\", \"ethylene\"
        //  \"C2H4O\", \"acetaldehyde\"
        //  \"C2H5OH\", \"ethanol\"
        //  \"C2H6\", \"ethane\"
        //  \"C2H6O2\", \"ethylene glycol\"
        //  \"C3H6\", \"propylene\"
        //  \"C3H6O\", \"acetone\"
        //  \"C3H8\", \"propane\"
        {L"7440-37-1", {L"AR", L"argon"}},
        {L"74-86-2", {L"C2H2", L"acetylene"}},
        {L"74-85-1", {L"C2H4", L"ethylene"}},
        {L"75-07-0", {L"C2H4O", L"acetaldehyde"}},
        {L"64-17-5", {L"C2H5OH", L"ethanol"} },
        {L"74-84-0", {L"C2H6", L"ethane"}},
        {L"107-21-1", {L"C2H6O2", L"ethylene glycol"}},
        {L"115-07-1", {L"C3H6", L"propylene"}},
        {L"67-64-1",   {L"C3H6O",   L"acetone"} },
        {L"74-98-6",   {L"C3H8",    L"propane"}},

        //  \"C4H10O3\", \"diethylene glycol\"
        //  \"C4H6\", \"1,3-butadiene\"
        //  \"C6H6\", \"benzene\"
        //  \"C7H14\", \"methylcyclohexane\"
        //  \"C7H8\", \"toluene\"
        //  \"C9H12-1\", \"n-propylbenzene\"
        //  \"C9H12-2\", \"mezitelen\"
        //  \"C9H18\", \"n-propylcyclohexane\"
        //  \"CH2O\", \"formaldehyde\"
        //  \"CH3OH\", \"methanol\"
        {L"111-46-6",  {L"C4H10O3", L"diethylene glycol"}},
        {L"106-99-0",  {L"C4H6",     L"1,3-butadiene"}},
        {L"71-43-2",   {L"C6H6", L"benzene"}},
        {L"108-87-2",  {L"C7H14",    L"methylcyclohexane"}},
        {L"108-88-3",  {L"C7H8",     L"toluene"}},
        {L"103-65-1",  {L"C9H12-1",  L"n-propylbenzene"}},
        {L"108-67-8",  {L"C9H12-2",  L"mezitelen"}},
        {L"1678-92-8", {L"C9H18",   L"n-propylcyclohexane"}},
        {L"50-00-0",   {L"CH2O",     L"formaldehyde"}},
        {L"67-56-1",   {L"CH3OH",   L"methanol"}},


        //  \"CH4\", \"methane\"
        //  \"CH4S\", \"methanethiol\"
        //  \"CO\", \"carbon monoxide\"
        //  \"CO2\", \"carbon dioxide\"
        //  \"COS\", \"carbonyl sulfide\"
        //  \"CS2\", \"carbon disulfide\"
        //  \"H2\", \"hydrogen\"
        //  \"H2O\", \"water\"
        //  \"H2S\", \"hydrogen sulfide\"
        //  \"He\", \"helium\"
        {L"74-82-8",   {L"CH4",     L"methane"}},
        {L"74-93-1",   {L"CH4S",    L"methanethiol"}},
        {L"630-08-0",  {L"CO",      L"carbon monoxide"}},
        {L"124-38-9",  {L"CO2",     L"carbon dioxide"}},
        {L"463-58-1",  {L"COS",     L"carbonyl sulfide"}},
        {L"75-15-0",   {L"CS2",     L"carbon disulfide"}},
        {L"1333-74-0", {L"H2",      L"hydrogen"}},
        {L"7732-18-5", {L"H2O",      L"water"}},
        {L"7783-06-4", {L"H2S",     L"hydrogen sulfide"}},
        {L"7440-59-7", {L"He",       L"helium"}},

        //  \"N2\", \"nitrogen\"
        //  \"NH3\", \"ammonia\"
        //  \"NO\", \"nitric oxide\"
        //  \"NO2\", \"nitrogen dioxide\"
        //  \"O2\", \"oxygen\"
        //  \"SO2\", \"sulfur dioxide\"
        //  \"Xe\", \"xenon\"
        //  \"c_C4H8\", \"cyclobutane\"
        //  \"c_C5H10\", \"cyclopentane\"
        //  \"c_C6H12\", \"cyclohexane\"
        {L"7727-37-9", {L"N2",       L"nitrogen"}},
        {L"7664-41-7", {L"NH3",      L"ammonia"}},
        {L"10102-43-9",{L"NO",       L"nitric oxide"}},
        {L"10102-44-0",{L"NO2",     L"nitrogen dioxide"}},
        {L"7782-44-7", {L"O2",      L"oxygen"}},
        {L"7446-09-5", {L"SO2",     L"sulfur dioxide"}},
        {L"7440-63-3", {L"Xe",       L"xenon"}},
        {L"287-23-0",  {L"c_C4H8",   L"cyclobutane"}},
        {L"287-92-3",  {L"c_C5H10",  L"cyclopentane"}},
        {L"110-82-7",  {L"c_C6H12", L"cyclohexane"}},


        //  \"cis_C10H18\", \"cis-decalin\"
        //  \"i_C10H22\", \"nonane, 2-methyl\"
        //  \"i_C12H26\", \"heptane, 2,2,4,6,6-pentamethyl\"
        //  \"i_C16H34\", \"nonane 2,2,4,4,6,8,8-heptamethyl\"
        //  \"i_C4H10\", \"isobutane\"
        //  \"i_C5H12\", \"isopentane\"
        //  \"i_C8H18\", \"2,2,4-trimethylpentane\"
        //  \"n_C10H22\", \"n-decane\"
        //  \"n_C11H24\", \"n-undecane\"
        //  \"n_C12H26\", \"n-dodecane\"
        {L"493-01-6",   {L"cis_C10H18", L"cis-decalin"}},
        {L"871-83-0",   {L"i_C10H22",   L"2-methylnonane"}},
        {L"13475-82-6", {L"i_C12H26",   L"2,2,4,6,6-pentamethylheptane"}},
        {L"4390-04-9",  {L"i_C16H34",   L"2,2,4,4,6,8,8-heptamethylnonane"}},
        {L"75-28-5",    {L"i_C4H10",    L"isobutane"}},
        {L"78-78-4",    {L"i_C5H12",    L"isopentane"}},
        {L"540-84-1",   {L"i_C8H18",    L"2,2,4-trimethylpentane"}},
        {L"124-18-5",   {L"n_C10H22",   L"n-decane"}},
        {L"1120-21-4",  {L"n_C11H24",   L"n-undecane"}},
        {L"112-40-3",   {L"n_C12H26",   L"n-dodecane"}},

        //  \"n_C13H28\", \"n-tridecane\"
        //  \"n_C14H30\", \"n-tetradecane\"
        //  \"n_C15H32\", \"n-pentadecane\"
        //  \"n_C16H34\", \"n-hexadecane\"
        //  \"n_C17H36\", \"n-heptadecane\"
        //  \"n_C18H38\", \"n-octadecane\"
        //  \"n_C3H8O\", \"propanol\"
        //  \"n_C4H10\", \"n-butane\"
        //  \"n_C4H10O\", \"n-butanol\"
        //  \"n_C5H12\", \"n-pentane\"
        {L"629-50-5",  {L"n_C13H28", L"n-tridecane"}},
        {L"629-59-4",  {L"n_C14H30", L"n-tetradecane"}},
        {L"629-62-9",  {L"n_C15H32", L"n-pentadecane"}},
        {L"544-76-3",  {L"n_C16H34", L"n-hexadecane"}},
        {L"629-78-7",  {L"n_C17H36", L"n-heptadecane"}},
        {L"593-45-3",  {L"n_C18H38", L"n-octadecane"}},
        {L"71-23-8",   {L"n_C3H8O",  L"1-propanol"}},
        {L"106-97-8",  {L"n_C4H10",  L"n-butane"}},
        {L"71-36-3",   {L"n_C4H10O", L"1-butanol"}},
        {L"109-66-0",  {L"n_C5H12",  L"n-pentane"}},

        //  \"n_C6H14\", \"n-hexane\"
        //  \"n_C7H16\", \"n-heptane\"
        //  \"n_C8H18\", \"n-octane\"
        //  \"n_C9H20\", \"n-nonane\"
        //  \"p_C8H10\", \"p xylene\"
        //  \"t_C10H18\", \"trans-decalin\"
        {L"110-54-3", {L"n_C6H14",  L"n-hexane"}},
        {L"142-82-5", {L"n_C7H16",  L"n-heptane"}},
        {L"111-65-9", {L"n_C8H18",  L"n-octane"}},
        {L"111-84-2", {L"n_C9H20",  L"n-nonane"}},
        {L"106-42-3", {L"p_C8H10",  L"p-xylene"}},
        {L"493-02-7", {L"t_C10H18", L"trans-decalin"}}
    };
    return target_component_CAS_map;
};
