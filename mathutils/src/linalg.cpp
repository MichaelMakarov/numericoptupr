#include <linalg.hpp>
#include <stdexcept>
#include <gsl/gsl_linalg.h>

matrix::matrix() noexcept : _ptr{nullptr}, _rows{}, _cols{} {}

matrix::matrix(size_t rows, size_t columns) : _ptr{nullptr}, _rows{rows}, _cols{columns}
{
    if (rows * columns > 0)
    {
        auto ptr = gsl_matrix_alloc(rows, columns);
        if (!ptr)
        {
            throw std::bad_alloc();
        }
        gsl_matrix_set_zero(ptr);
        _ptr = ptr;
    }
}

matrix::matrix(matrix const &other) : matrix(other._rows, other._cols)
{

    auto this_ptr = reinterpret_cast<gsl_matrix *>(_ptr);
    auto other_ptr = reinterpret_cast<gsl_matrix *>(other._ptr);
    gsl_matrix_memcpy(this_ptr, other_ptr);
}

matrix::matrix(matrix &&other) noexcept : _ptr{other._ptr}, _rows{other._rows}, _cols{other._cols}
{
    other._ptr = nullptr;
}

matrix::~matrix()
{
    if (_ptr)
    {
        auto this_ptr = reinterpret_cast<gsl_matrix *>(_ptr);
        gsl_matrix_free(this_ptr);
    }
}

matrix &matrix::operator=(matrix const &other)
{
    matrix::~matrix();
    matrix::matrix(other);
    return *this;
}

matrix &matrix::operator=(matrix &&other) noexcept
{
    std::swap(_rows, other._rows);
    std::swap(_cols, other._cols);
    std::swap(_ptr, other._ptr);
    return *this;
}

double *matrix::operator[](std::size_t row)
{
    auto this_ptr = reinterpret_cast<gsl_matrix *>(_ptr);
    return this_ptr->data + row * this_ptr->tda;
}

double const *matrix::operator[](std::size_t row) const
{
    auto this_ptr = reinterpret_cast<gsl_matrix *>(_ptr);
    return this_ptr->data + row * this_ptr->tda;
}

vector::vector() noexcept : _ptr{nullptr}, _size{} {}

vector::vector(std::size_t size, double elem) : _size{size}
{
    if (size > 0)
    {
        auto ptr = gsl_vector_alloc(size);
        gsl_vector_set_all(ptr, elem);
    }
}

vector::vector(vector const &other) : vector(other._size)
{
    auto this_ptr = reinterpret_cast<gsl_vector *>(_ptr);
    auto other_ptr = reinterpret_cast<gsl_vector *>(other._ptr);
    gsl_vector_memcpy(this_ptr, other_ptr);
}

vector::vector(vector &&other) noexcept : _ptr{other._ptr}, _size{other._size}
{
    other._ptr = nullptr;
}

vector::~vector()
{
    if (_ptr)
    {
        auto this_ptr = reinterpret_cast<gsl_vector *>(_ptr);
        gsl_vector_free(this_ptr);
    }
}

vector &vector::operator=(vector const &other)
{
    vector::~vector();
    vector::vector(other);
    return *this;
}

vector &vector::operator=(vector &&other) noexcept
{
    std::swap(_ptr, other._ptr);
    std::swap(_size, other._size);
    return *this;
}

double &vector::operator[](std::size_t i)
{
    auto this_ptr = reinterpret_cast<gsl_vector *>(_ptr);
    return this_ptr->data[i * this_ptr->stride];
}

double vector::operator[](std::size_t i) const
{
    auto this_ptr = reinterpret_cast<gsl_vector *>(_ptr);
    return this_ptr->data[i * this_ptr->stride];
}

vector solve(matrix const &mx, vector const &vc)
{
    vector rs(vc._size);
    auto mx_ptr = reinterpret_cast<gsl_matrix *>(mx._ptr);
    auto vc_ptr = reinterpret_cast<gsl_vector *>(vc._ptr);
    auto rs_ptr = reinterpret_cast<gsl_vector *>(rs._ptr);
    auto pm_ptr = gsl_permutation_alloc(vc._size);
    int s;
    gsl_linalg_LU_decomp(mx_ptr, pm_ptr, &s);
    gsl_linalg_LU_solve(mx_ptr, pm_ptr, vc_ptr, rs_ptr);
    return rs;
}