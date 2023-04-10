#pragma once
#include <cstddef>

class vector;

class matrix
{
    void *_ptr;
    std::size_t _rows, _cols;

public:
    matrix() noexcept;
    matrix(size_t rows, size_t columns);
    matrix(matrix const &);
    matrix(matrix &&) noexcept;
    ~matrix();
    matrix &operator=(matrix const &);
    matrix &operator=(matrix &&) noexcept;
    std::size_t rows() const { return _rows; }
    std::size_t columns() const { return _cols; }
    double *operator[](std::size_t row);
    double const *operator[](std::size_t row) const;

    friend vector solve(matrix const &mx, vector const &vc);
};

class vector
{
    void *_ptr;
    std::size_t _size;

public:
    vector() noexcept;
    vector(std::size_t size, double elem = 0);
    vector(vector const &);
    vector(vector &&) noexcept;
    ~vector();
    vector &operator=(vector const &);
    vector &operator=(vector &&) noexcept;
    std::size_t size() const { return _size; }
    double &operator[](std::size_t i);
    double operator[](std::size_t i) const;

    friend vector solve(matrix const &mx, vector const &vc);
};