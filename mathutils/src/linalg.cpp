#include <linalg.hpp>
#include <stdexcept>
#include <gsl/gsl_linalg.h>

template <typename T>
auto cast_ptr(void *ptr)
{
    return reinterpret_cast<T *>(ptr);
}

matrix::matrix() noexcept : _ptr{nullptr} {}

matrix::matrix(size_t rows, size_t columns) : matrix()
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

matrix::matrix(matrix const &other) : matrix()
{
    copy(other);
}

matrix::matrix(matrix &&other) noexcept : _ptr{other._ptr}
{
    std::swap(_ptr, other._ptr);
}

matrix::~matrix()
{
    clear();
}

void matrix::clear()
{
    if (_ptr)
    {
        gsl_matrix_free(cast_ptr<gsl_matrix>(_ptr));
    }
}

void matrix::copy(matrix const &other)
{
    if (other._ptr)
    {
        auto other_ptr = cast_ptr<gsl_matrix>(other._ptr);
        auto this_ptr = gsl_matrix_alloc(other_ptr->size1, other_ptr->size2);
        gsl_matrix_memcpy(this_ptr, other_ptr);
        _ptr = this_ptr;
    }
}

matrix &matrix::operator=(matrix const &other)
{
    clear();
    copy(other);
    return *this;
}

matrix &matrix::operator=(matrix &&other) noexcept
{
    std::swap(_ptr, other._ptr);
    return *this;
}

double *matrix::operator[](std::size_t row)
{
    auto this_ptr = cast_ptr<gsl_matrix>(_ptr);
    return this_ptr->data + row * this_ptr->tda;
}

double const *matrix::operator[](std::size_t row) const
{
    auto this_ptr = cast_ptr<gsl_matrix>(_ptr);
    return this_ptr->data + row * this_ptr->tda;
}

std::size_t matrix::rows() const
{
    return _ptr ? cast_ptr<gsl_matrix>(_ptr)->size1 : 0;
}

std::size_t matrix::columns() const
{
    return _ptr ? cast_ptr<gsl_matrix>(_ptr)->size2 : 0;
}

matrix operator+(matrix const &left, matrix const &right)
{
    matrix result(left);
    if (0 != gsl_matrix_add(cast_ptr<gsl_matrix>(result._ptr), cast_ptr<gsl_matrix>(right._ptr)))
        throw std::runtime_error("Failed to add matrices.");
    return result;
}

matrix operator-(matrix const &left, matrix const &right)
{
    matrix result(left);
    if (0 != gsl_matrix_sub(cast_ptr<gsl_matrix>(result._ptr), cast_ptr<gsl_matrix>(right._ptr)))
        throw std::runtime_error("Failed to sub matrices.");
    return result;
}

matrix operator*(matrix const &left, matrix const &right)
{
    auto left_ptr = cast_ptr<gsl_matrix>(left._ptr);
    auto right_ptr = cast_ptr<gsl_matrix>(right._ptr);
    matrix result(left_ptr->size1, right_ptr->size2);
    auto result_ptr = cast_ptr<gsl_matrix>(result._ptr);
    if (0 != gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, left_ptr, right_ptr, 0, result_ptr))
    {
        throw std::runtime_error("Failed to multiply matrices.");
    }
    return result;
}

matrix operator*(matrix const &in, double v)
{
    matrix out(in);
    gsl_matrix_scale(cast_ptr<gsl_matrix>(in._ptr), v);
    return out;
}

matrix operator/(matrix const &in, double v)
{
    return in * (1 / v);
}

class permutation
{
    gsl_permutation *_ptr;

public:
    permutation(std::size_t size)
    {
        _ptr = gsl_permutation_alloc(size);
        if (!_ptr)
        {
            throw std::bad_alloc();
        }
    }
    ~permutation()
    {
        gsl_permutation_free(_ptr);
    }
    friend matrix inverse(matrix const &);
};

matrix inverse(matrix const &mx)
{
    matrix buf(mx);
    matrix inv(mx);
    permutation perm(mx.rows());
    int s;
    auto buf_ptr = cast_ptr<gsl_matrix>(buf._ptr);
    auto inv_ptr = cast_ptr<gsl_matrix>(inv._ptr);
    auto perm_ptr = cast_ptr<gsl_permutation>(perm._ptr);
    if (0 != gsl_linalg_LU_decomp(buf_ptr, perm_ptr, &s))
    {
        throw std::runtime_error("Failed to compute LU decomposition.");
    }
    if (0 != gsl_linalg_LU_invert(buf_ptr, perm_ptr, inv_ptr))
    {
        throw std::runtime_error("Failed to compute matrix LU inversion.");
    }
    return inv;
}

vector::vector() noexcept : _ptr{nullptr} {}

vector::vector(std::size_t size, double elem)
{
    if (size > 0)
    {
        auto ptr = gsl_vector_alloc(size);
        gsl_vector_set_all(ptr, elem);
        _ptr = ptr;
    }
}

vector::vector(vector const &other) : vector()
{
    copy(other);
}

vector::vector(vector &&other) noexcept : vector()
{
    std::swap(_ptr, other._ptr);
}

vector::~vector()
{
    clear();
}

void vector::clear()
{
    if (_ptr)
    {
        gsl_vector_free(cast_ptr<gsl_vector>(_ptr));
    }
}

void vector::copy(vector const &other)
{
    if (other._ptr)
    {
        auto other_ptr = cast_ptr<gsl_vector>(other._ptr);
        auto this_ptr = gsl_vector_alloc(other_ptr->size);
        gsl_vector_memcpy(this_ptr, other_ptr);
        _ptr = this_ptr;
    }
}

vector &vector::operator=(vector const &other)
{
    clear();
    copy(other);
    return *this;
}

vector &vector::operator=(vector &&other) noexcept
{
    std::swap(_ptr, other._ptr);
    return *this;
}

std::size_t vector::size() const
{
    return _ptr ? cast_ptr<gsl_vector>(_ptr)->size : 0;
}

double &vector::operator[](std::size_t i)
{
    auto this_ptr = cast_ptr<gsl_vector>(_ptr);
    return this_ptr->data[i * this_ptr->stride];
}

double vector::operator[](std::size_t i) const
{
    auto this_ptr = cast_ptr<gsl_vector>(_ptr);
    return this_ptr->data[i * this_ptr->stride];
}

vector operator+(vector const &left, vector const &right)
{
    vector result(left);
    gsl_vector_add(cast_ptr<gsl_vector>(result._ptr), cast_ptr<gsl_vector>(right._ptr));
    return result;
}

vector operator-(vector const &left, vector const &right)
{
    vector result(left);
    gsl_vector_sub(cast_ptr<gsl_vector>(result._ptr), cast_ptr<gsl_vector>(right._ptr));
    return result;
}

vector operator*(vector const &left, double right)
{
    vector result(left);
    gsl_vector_scale(cast_ptr<gsl_vector>(result._ptr), right);
    return result;
}

vector operator/(vector const &left, double right)
{
    return left * (1 / right);
}

double operator*(vector const &left, vector const &right)
{
    double res{};
    gsl_blas_ddot(cast_ptr<gsl_vector>(left._ptr), cast_ptr<gsl_vector>(right._ptr), &res);
    return res;
}

vector operator*(matrix const &left, vector const &right)
{
    vector result(left.rows());
    if (0 != gsl_blas_dgemv(CblasNoTrans, 1,
                            cast_ptr<gsl_matrix>(left._ptr),
                            cast_ptr<gsl_vector>(right._ptr), 0,
                            cast_ptr<gsl_vector>(result._ptr)))
        throw std::runtime_error("Failed to multiply matrix and vector.");
    return result;
}

// vector solve(matrix const &mx, vector const &vc)
// {
//     vector rs(vc._size);
//     auto mx_ptr = reinterpret_cast<gsl_matrix *>(mx._ptr);
//     auto vc_ptr = reinterpret_cast<gsl_vector *>(vc._ptr);
//     auto rs_ptr = reinterpret_cast<gsl_vector *>(rs._ptr);
//     auto pm_ptr = gsl_permutation_alloc(vc._size);
//     int s;
//     gsl_linalg_LU_decomp(mx_ptr, pm_ptr, &s);
//     gsl_linalg_LU_solve(mx_ptr, pm_ptr, vc_ptr, rs_ptr);
//     return rs;
// }