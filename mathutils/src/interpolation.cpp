#include <interpolation.hpp>
#include <gsl/gsl_spline.h>
#include <algorithm>

class spline_impl
{
    gsl_interp_accel *_acc;
    gsl_spline *_spl;

public:
    spline_impl(std::size_t degree)
    {
        _acc = gsl_interp_accel_alloc();
        _spl = gsl_spline_alloc(gsl_interp_cspline, degree + 1);
    }
    ~spline_impl()
    {
        gsl_spline_free(_spl);
        gsl_interp_accel_free(_acc);
    }
    void init(double const *x, double const *y, std::size_t count)
    {
        gsl_spline_init(_spl, x, y, count);
    }
    double eval(double x) const
    {
        return gsl_spline_eval(_spl, x, _acc);
    }
};

spline::spline(std::size_t degree) : _impl{new spline_impl{degree}}, _degree{degree}
{
}

spline::spline(spline &&other) noexcept
{
    std::swap(_impl, other._impl);
}

spline::~spline()
{
    delete _impl;
}

spline &spline::operator=(spline &&other) noexcept
{
    std::swap(_impl, other._impl);
    return *this;
}

double spline::operator()(double x) const
{
    return _impl->eval(x);
}

spline make_spline(double const *x, double const *y, std::size_t count)
{
    spline spl{count};
    spl._impl->init(x, y, count);
    return spl;
}