#include <functions.hpp>
#include <cmath>

bool is_equal(double left, double right, double eps)
{
    double diff = std::abs(left - right);
    return diff < eps || diff / std::abs(left) < eps;
}