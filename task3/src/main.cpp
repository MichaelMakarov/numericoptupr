#include <linalg.hpp>
#include <integrator.hpp>
#include <functions.hpp>

#include <vector>
#include <fstream>
#include <iostream>

std::ostream &operator<<(std::ostream &os, vector const &v);

class equation
{
    double u;

public:
    equation(double u_ = 0) : u{u_} {}
    vector operator()(vector const &x, double) const
    {
        vector dx(x.size());
        dx[0] = 1;
        dx[1] = x[2];
        dx[2] = sign(x[4]);
        dx[3] = 0;
        dx[4] = -x[3];
        return dx;
    }
};

constexpr size_t residual_dim{2};

vector residuals(vector const &x)
{
    vector dx(residual_dim);
    dx[0] = x[1];
    dx[1] = x[2];
    return dx;
}

auto integrate(vector const &x0, size_t count, std::ostream &os)
{
    std::vector<vector> arr(count);
    arr.front() = x0;
    double dt = 1. / (count - 1);
    for (std::size_t i{1}; i < count; ++i)
    {
        double tk = dt * i;
        double tn = tk - dt;
        rk4(tn, arr[i - 1], equation(), tk, arr[i]);
        os << "t = " << tk << ' ' << arr[i] << std::endl;
    }
    return arr.back();
}

void print(std::ostream &os, vector const &x0, vector const &dx);

void solve(vector x0, std::size_t count, std::ostream &os)
{
    constexpr double mult{1.01};
    double prevres{};
    std::size_t iteration{};
    while (true)
    {
        vector x01{x0}, x02{x0};
        x01[3] *= mult;
        x02[4] *= mult;
        os << "Iteration №" << ++iteration << std::endl;
        vector xk = integrate(x0, count, os);
        vector dx = residuals(xk);
        print(os, x0, dx);
        double curres = dx * dx;
        if (is_equal(curres, prevres, 0.01))
        {
            break;
        }
        prevres = curres;
        vector xk1 = integrate(x01, count, os);
        vector dx1 = residuals(xk1);
        print(os, x01, dx1);
        vector xk2 = integrate(x02, count, os);
        vector dx2 = residuals(xk2);
        print(os, x02, dx2);
        matrix mx(residual_dim, residual_dim);
        double db{x01[3] - x0[3]};
        mx[0][0] = (dx1[1] - dx[1]) / db;
        mx[1][0] = (dx1[2] - dx[2]) / db;
        db = x02[4] - x0[4];
        mx[0][1] = (dx2[1] - dx[1]) / db;
        mx[1][1] = (dx2[2] - dx[2]) / db;
        vector dx0 = solve(mx, dx);
        x0[3] += dx0[0];
        x0[4] += dx0[1];
    }
}

int main()
{
    constexpr double x0{4};
    constexpr double y0{2};
    constexpr double b1{-5};
    constexpr double b2{-3};
    try
    {
        std::ofstream fout{"output.txt"};
        if (!fout.is_open())
        {
            throw std::runtime_error("Failed to open output file.");
        }
        vector x(5);
        x[0] = 0;
        x[1] = x0;
        x[2] = y0;
        x[3] = b1;
        x[4] = b2;
        solve(x, 100, fout);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    return 0;
}

std::ostream &operator<<(std::ostream &os, vector const &v)
{
    for (size_t i{}; i < v.size(); ++i)
    {
        os << v[i] << ' ';
    }
    return os;
}

void print(std::ostream &os, vector const &x0, vector const &dx)
{
    os << "x0: " << x0 << std::endl;
    os << "dx: " << dx << std::endl;
}