#pragma once
#ifndef __FLUID_STUBDATA__
#define __FLUID_STUBDATA__

/// @brief Состав потока
struct fluid_stubdata_t {
    /// @brief Список чистых компонентов из БД
    std::vector<std::wstring> component_list;
    /// @brief Мольный состав смеси
    std::vector<double> molar_fraction;

#ifdef VLELIB_SERIALIZATION_SUPPORT
    /// @brief Для сериализации
    friend class boost::serialization::access;
    /// @brief Интрузивная сериализация/десериализация
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_NVP(component_list);
        ar& BOOST_SERIALIZATION_NVP(molar_fraction);
    }
#endif
};

#endif
