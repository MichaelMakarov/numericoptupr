#pragma once
#include <vector>

struct node
{
    double t;
    double x;
    double p;
    double u;
};

struct solution
{
    std::vector<node> nodes;
    double funcval{};

    solution() = default;
    solution(solution &&other) noexcept;
    solution &operator=(solution &&other) noexcept;
};

solution solve(double x0, double t0, double tk, double dt, double eps, std::size_t maxiter, char const *filepath);
