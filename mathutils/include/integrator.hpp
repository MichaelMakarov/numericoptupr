#pragma once

template <typename T, typename V, typename F>
void rk4(T const &t0, V const &v0, F &&func, T const &tk, V &vk)
{
    double step = tk - t0, step_2 = 0.5 * step, step_6 = step_2 / 3;
    T t = (t0 + tk) / 2;
    V k1 = func(v0, t0);
    V k2 = func(v0 + k1 * step_2, t);
    V k3 = func(v0 + k2 * step_2, t);
    V k4 = func(v0 + k3 * step, tk);
    vk = v0 + (k1 + (k2 + k3) * 2.0 + k4) * step_6;
}