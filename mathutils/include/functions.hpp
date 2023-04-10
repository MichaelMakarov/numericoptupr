#pragma once

template <typename T>
T sign(T x)
{
    constexpr T zero{};
    return (zero < x) - (x < zero);
}

bool is_equal(double left, double right, double eps);