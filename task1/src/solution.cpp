#include <algorithm>
#include <solution.hpp>
#include <integrator.hpp>
#include <functions.hpp>
#include <fstream>

struct coordinate_func
{
    double u;
    double operator()(double, double) const
    {
        return u;
    }
};

struct impulse_func
{
    double x;
    double operator()(double, double) const
    {
        return 2 * x;
    }
};

void integrate_coordinate(std::vector<node> &nodes)
{
    for (std::size_t i{1}; i < nodes.size(); ++i)
    {
        auto &prev = nodes[i - 1];
        auto &curr = nodes[i];
        rk4(prev.t, prev.x, coordinate_func{prev.u}, curr.t, curr.x);
    }
}

void integrate_impulse(std::vector<node> &nodes)
{
    for (std::size_t i{nodes.size() - 1}; i >= 1; --i)
    {
        auto &prev = nodes[i];
        auto &curr = nodes[i - 1];
        rk4(prev.t, prev.p, impulse_func{prev.x}, curr.t, curr.p);
    }
}

double functional_value(std::vector<node> const &nodes)
{
    double res{};
    for (std::size_t i{1}; i < nodes.size(); ++i)
    {
        auto &prev = nodes[i - 1];
        auto &curr = nodes[i];
        res += 0.5 * (prev.x * prev.x + curr.x * curr.x) * (curr.t - prev.t);
    }
    return res;
}

auto make_nodes(double x0, double t0, double tk, double dt)
{
    std::size_t count = static_cast<std::size_t>((tk - t0) / dt + 1);
    dt = (tk - t0) / (count - 1);
    std::vector<node> nodes(count);
    auto &first = nodes.front();
    first.x = x0;
    first.u = -sign(x0);
    first.t = t0;
    for (std ::size_t i{1}; i < nodes.size(); ++i)
    {
        auto &prev = nodes[i - 1];
        auto &curr = nodes[i];
        curr.t = prev.t + dt;
        curr.u = prev.u;
    }
    return nodes;
}

void print(std::ostream &os, std::vector<node> const &nodes)
{
    for (auto &n : nodes)
    {
        os << "t = " << n.t << ' '
           << "x = " << n.x << ' '
           << "p = " << n.p << ' '
           << "u = " << n.u << std::endl;
    }
}

solution solve(double x0, double t0, double tk, double dt, double eps, std::size_t maxiter, char const *filepath)
{
    solution sol;
    sol.nodes = make_nodes(x0, t0, tk, dt);
    std::ofstream fout{filepath};
    double prev{};
    for (std::size_t iter{1}; iter <= maxiter; ++iter)
    {
        integrate_coordinate(sol.nodes);
        integrate_impulse(sol.nodes);
        double curr = functional_value(sol.nodes);
        fout << "----------------------------\n";
        fout << "iteration №" << iter << ' ';
        fout << "(J = " << curr << ")\n";
        print(fout, sol.nodes);
        if (is_equal(curr, prev, eps))
        {
            sol.funcval = curr;
            break;
        }
        for (auto &n : sol.nodes)
        {
            n.u += 1e-2 * n.p;
            n.u = std::min(1., std::max(-1., n.u));
        }
        prev = curr;
    }
    return sol;
}

solution::solution(solution &&other) noexcept : nodes{std::move(other.nodes)}, funcval{other.funcval}
{
}

solution &solution::operator=(solution &&other) noexcept
{
    nodes.swap(other.nodes);
    funcval = other.funcval;
    return *this;
}