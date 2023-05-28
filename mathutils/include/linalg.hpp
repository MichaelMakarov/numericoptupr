#pragma once
#include <cstddef>

class vector;

class matrix
{
    void *_ptr;

    void clear();
    void copy(matrix const &);

public:
    matrix() noexcept;
    matrix(size_t rows, size_t columns);
    matrix(matrix const &);
    matrix(matrix &&) noexcept;
    ~matrix();
    matrix &operator=(matrix const &);
    matrix &operator=(matrix &&) noexcept;
    std::size_t rows() const;
    std::size_t columns() const;
    double *operator[](std::size_t row);
    double const *operator[](std::size_t row) const;

    friend matrix operator+(matrix const &, matrix const &);
    friend matrix operator-(matrix const &, matrix const &);
    friend matrix operator*(matrix const &, matrix const &);
    friend matrix operator*(matrix const &, double);
    friend matrix operator/(matrix const &, double);
    friend vector operator*(matrix const &, vector const &);
    friend matrix inverse(matrix const &);

    friend vector solve(matrix const &mx, vector const &vc);
};

class vector
{
    void *_ptr;

    void clear();
    void copy(vector const &);

public:
    vector() noexcept;
    vector(std::size_t size, double elem = 0);
    vector(vector const &);
    vector(vector &&) noexcept;
    ~vector();
    vector &operator=(vector const &);
    vector &operator=(vector &&) noexcept;
    std::size_t size() const;
    double &operator[](std::size_t i);
    double operator[](std::size_t i) const;

    friend vector operator+(vector const &, vector const &);
    friend vector operator-(vector const &, vector const &);
    friend vector operator*(vector const &, double);
    friend vector operator/(vector const &, double);
    friend vector operator*(matrix const &, vector const &);
    friend double operator*(vector const &, vector const &);

    friend vector solve(matrix const &mx, vector const &vc);
};