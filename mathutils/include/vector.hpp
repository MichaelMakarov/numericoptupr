#pragma once
#include <array>

template <std::size_t size>
using vec = std::array<double, size>;

template <std::size_t size>
constexpr vec<size> operator+(vec<size> const &left, vec<size> const &right)
{
    vec<size> res{};
    for (std::size_t i{}; i < size; ++i)
    {
        res[i] = left[i] + right[i];
    }
    return res;
}

template <std::size_t size>
constexpr vec<size> operator-(vec<size> const &left, vec<size> const &right)
{
    vec<size> res{};
    for (std::size_t i{}; i < size; ++i)
    {
        res[i] = left[i] - right[i];
    }
    return res;
}

template <std::size_t size>
constexpr vec<size> operator*(vec<size> const &left, double right)
{
    vec<size> res{};
    for (std::size_t i{}; i < size; ++i)
    {
        res[i] = left[i] * right;
    }
    return res;
}

template <std::size_t size>
constexpr vec<size> operator*(double left, vec<size> const &right)
{
    return right * left;
}

template <std::size_t size>
constexpr vec<size> operator/(vec<size> const &left, double right)
{
    return left * (1 / right);
}