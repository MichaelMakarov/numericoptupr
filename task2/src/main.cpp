#include <iostream>
#include <vector.hpp>
#include <linalg.hpp>
#include <integrator.hpp>
#include <cmath>
#include <vector>
#include <fstream>

/**
 * @brief Вектор (x1, x2, x3, x4, dx1, dx2, dx3, dx4)
 *
 */
using vec8 = vec<8>;
/**
 * @brief Вектор (x1, x2, x3, x4)
 *
 */
using vec4 = vec<4>;

template <typename node, double node::*time>
auto make_grid(double tn, double tk, std::size_t count)
{
    double dt = (tk - tn) / (count - 1);
    std::vector<node> nodes(count);
    for (std::size_t i{}; auto &n : nodes)
    {
        n.*time = tn + i++ * dt;
    }
    return nodes;
}

struct xnode
{
    double t;
    vector v;
    matrix m;
};

struct ynode
{
    double q;
    vector y;
    vector vn;
    vector vk;
    matrix mn;
    matrix mk;
};

vector residual(double xn, double yn, double xk, double yk);

std::ostream &operator<<(std::ostream &os, vector const &v)
{
    for (std::size_t i{}; i < v.size(); ++i)
    {
        os << ' ' << v[i];
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, ynode const &n)
{
    auto res = residual(n.vn[0], n.vn[1], n.vk[0], n.vk[1]);
    return os << "q = " << n.q
              << " y =" << n.y
              << " r = " << res
              << " |r|^2 = " << (res * res);
}

void print_nodes(std::ostream &os, std::vector<ynode> const &nodes)
{
    std::size_t index{};
    os << "Computed nodes:\n";
    for (auto &n : nodes)
    {
        os << ++index << ' ' << n << std::endl;
    }
    os << "---------------------------------\n";
    auto iter = std::min_element(std::begin(nodes), std::end(nodes), [](ynode const &left, ynode const &right)
                                 {
                                    auto lres = residual(left.vn[0], left.vn[1], left.vk[0], left.vk[1]);
                                    auto rres = residual(right.vn[0], right.vn[1], right.vk[0], right.vk[1]);
                                    return lres * lres < rres * rres; });
    os << "selected q = " << iter->q << std::endl;
}

constexpr double a1{2};
constexpr double a2{};
constexpr double b1{1.0738644361};
constexpr double b2{-1.0995343576};

vector residual(double xn, double yn, double xk, double yk)
{
    vector v(4);
    v[0] = xn - a1;
    v[1] = yn - b1;
    v[2] = xk - a2;
    v[3] = yk - b2;
    return v;
}

/**
 * R - вектор краевых условий:
 * x1(0) - a1
 * x2(0) - b1
 * x1(T) - a2
 * x2(T) - b2
 *
 * f(x1, x2, x3, x4) вектор-функция системы ДУ:
 * f1 = x3
 * f2 = x4
 * f3 = -x1 / (x1 * x1 + x2 * x2)^(3/2)
 * f4 = -x2 / (x1 * x1 + x2 * x2)^(3/2)
 *
 * A - матрица линеаризованной системы ДУ:
 * строка 1: 0                                                          0                                                       1   0
 * строка 2: 0                                                          0                                                       0   1
 * строка 3: (3 * x1^2 / (x1^2 + x2^2) - 1) / (x1^2 + x2^2)^(3/2)       3 * x1 * x2 / (x1^2 + x2^2)^(5/2)                       0   0
 * строка 4: 3 * x1 * x2 / (x1^2 + x2^2)^(5/2)                          (3 * x2^2 / (x1^2 + x2^2) - 1) / (x1^2 + x2^2)^(3/2)    0   0
 */

inline double mixed_derivative(double x1, double x2, double r2, double r3_2)
{
    return 3 * x1 * x2 / (r2 * r3_2);
}

inline double major_derivative(double x, double r2, double r3_2)
{
    return (3 * x * x / r2 - 1) / r3_2;
}

vector equation(vector const &in, double)
{
    double r2 = in[0] * in[0] + in[1] * in[1];
    double r3_2 = std::sqrt(r2 * r2 * r2);
    vector out(in.size());
    out[0] = in[2];
    out[1] = in[3];
    out[2] = -in[0] / r3_2;
    out[3] = -in[1] / r3_2;
    return out;
}

class xequation_maker
{
    vector const &_v;

public:
    xequation_maker(vector const &v) : _v{v} {}
    matrix operator()(matrix const &in, double) const
    {
        double r2 = _v[0] * _v[0] + _v[1] * _v[1];
        double r3_2 = std::sqrt(r2 * r2 * r2);
        double der12 = mixed_derivative(_v[0], _v[1], r2, r3_2);
        double der1 = major_derivative(_v[0], r2, r3_2);
        double der2 = major_derivative(_v[1], r2, r3_2);
        matrix out(in.rows(), in.columns());
        out[0][2] = out[1][3] = 1;
        out[2][0] = der1;
        out[3][1] = der2;
        out[2][1] = out[3][0] = der12;
        return out;
    }
};

template <typename iterator>
void integrate(vector v0, double t0, iterator begin, iterator end)
{
    matrix m0(v0.size(), v0.size());
    for (std::size_t i{}; i < v0.size(); ++i)
        m0[i][i] = 1;
    while (begin != end)
    {
        rk4(t0, v0, &equation, begin->t, begin->v);
        rk4(t0, m0, xequation_maker(v0), begin->t, begin->m);
        t0 = begin->t;
        v0 = begin->v;
        m0 = begin->m;
        ++begin;
    }
}

auto solve_interior_equation(vector const &v0, double t0, double tn, double tk, std::size_t count)
{
    auto nodes = make_grid<xnode, &xnode::t>(tn, tk, count);
    auto begin = std::begin(nodes);
    auto end = std::end(nodes);
    auto refiter = std::upper_bound(begin, end, t0, [](double t, xnode const &n)
                                    { return t < n.t; });
    if (refiter == end)
        refiter = begin;
    using reverse_iterator_t = std::remove_reference_t<decltype(nodes)>::reverse_iterator;
    integrate(v0, t0, refiter, end);
    integrate(v0, t0, reverse_iterator_t(refiter), reverse_iterator_t(begin));
    return std::pair{nodes.front(), nodes.back()};
}

matrix get_dr_dxn()
{
    matrix m(4, 4);
    m[0][0] = m[1][1] = 1;
    return m;
}

matrix get_dr_dxk()
{
    matrix m(4, 4);
    m[2][0] = m[3][1] = 1;
    return m;
}

matrix inv_phi_derivative(matrix const &Xn, matrix const &Xk)
{
    const matrix dr_dxn = get_dr_dxn();
    const matrix dr_dxk = get_dr_dxk();
    matrix dphi = dr_dxn * Xn + dr_dxk * Xk;
    return inverse(dphi);
}

class yequation_maker
{
    vector _phi;
    matrix _dphi;

public:
    yequation_maker(vector const &phi, matrix const &dphi) : _phi{phi}, _dphi{dphi}
    {
    }
    vector operator()(vector const &, double) const
    {
        return (_dphi * _phi) * -1;
    }
};

auto solve_exterior_equation(vector y0, double t0, double tn, double tk, std::size_t count)
{
    vector phi;
    auto nodes = make_grid<ynode, &ynode::q>(0, 1, count);
    double q0{};
    for (auto &n : nodes)
    {
        auto [nn, nk] = solve_interior_equation(y0, t0, tn, tk, count);
        if (phi.size() == 0)
        {
            phi = residual(nn.v[0], nn.v[1], nk.v[0], nk.v[1]);
        }
        matrix dphi = inv_phi_derivative(nn.m, nk.m);
        rk4(q0, y0, yequation_maker(phi, dphi), n.q, n.y);
        n.vn = std::move(nn.v);
        n.mn = std::move(nn.m);
        n.vk = std::move(nk.v);
        n.mk = std::move(nk.m);
        q0 = n.q;
        y0 = n.y;
    }
    std::ofstream fout{"output.txt"};
    print_nodes(fout, nodes);
}

int main()
{
    vector y0(4);
    y0[0] = 2;
    y0[1] = 0;
    y0[2] = -0.5;
    y0[3] = 0.5;
    solve_exterior_equation(y0, 0, 0, 7, 70);
    return 0;
}