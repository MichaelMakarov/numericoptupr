#pragma once
#include <cstddef>

class spline
{
    class spline_impl *_impl;
    std::size_t _degree;

public:
    spline(std::size_t degree);
    spline(spline &&) noexcept;
    ~spline();
    spline &operator=(spline &&) noexcept;
    std::size_t degree() const { return _degree; }
    double operator()(double x) const;

    friend spline make_spline(double const *x, double const *y, std::size_t count);
};
