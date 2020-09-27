#pragma once

#include <iostream>

template <int... n>
struct tensor;

template <int n>
struct tensor<n> {
  static constexpr int size = n;
  static constexpr int rank = 1;
  constexpr auto& operator()(size_t i) { return values[i]; }
  constexpr auto operator()(size_t i) const { return values[i]; }
  double values[n];
};

template <int n1, int n2>
struct tensor<n1, n2> {
  static constexpr int size = n1 * n2;
  static constexpr int rank = 2;
  constexpr auto row(size_t i) const { 
    tensor<n2> r{};
    for (int j = 0; j < n2; j++){ r(j) = values[i][j]; }
    return r; 
  }
  constexpr void set_row(size_t i, tensor<n2> r) { 
    for (int j = 0; j < n2; j++){ values[i][j] = r(j); }
  }
  constexpr auto col(size_t j) const { 
    tensor<n1> c{};
    for (int i = 0; i < n1; i++){ c(i) = values[i][j]; }
    return c; 
  }
  constexpr void set_col(size_t j, tensor<n1> c) { 
    for (int i = 0; i < n1; i++){ values[i][j] = c(i); }
  }
  constexpr auto& operator()(size_t i, size_t j) { return values[i][j]; }
  constexpr auto operator()(size_t i, size_t j) const { return values[i][j]; }
  double values[n1][n2];
};

template <int n1, int n2, int n3>
struct tensor<n1, n2, n3> {
  static constexpr int size = n1 * n2 * n3;
  static constexpr int rank = 3;
  constexpr auto& operator()(size_t i, size_t j, size_t k) {
    return values[i][j][k];
  }
  constexpr auto operator()(size_t i, size_t j, size_t k) const {
    return values[i][j][k];
  }
  double values[n1][n2][n3];
};

template <int n1, int n2, int n3, int n4>
struct tensor<n1, n2, n3, n4> {
  static constexpr int size = n1 * n2 * n3 * n4;
  static constexpr int rank = 4;
  constexpr auto& operator()(size_t i, size_t j, size_t k, size_t l) {
    return values[i][j][k][l];
  }
  constexpr auto operator()(size_t i, size_t j, size_t k, size_t l) const {
    return values[i][j][k][l];
  }
  double values[n1][n2][n3][n4];
};

struct TransposeTag {};

template <int m, int n>
auto transpose(const tensor<m, n>& A) {
  tensor<n, m> AT{};
  for (auto [i, j] : IndexSpace<n, m>{}) {
    AT(i, j) = A(j, i);
  }
  return AT;
}

template <int m, int n>
auto operator^(const tensor<m, n> A, const TransposeTag) {
  return transpose(A);
}

template <int... n>
auto operator*(const tensor<n...> a, const double b) {
  tensor<n...> c{};
  for (int i = 0; i < (n * ...); i++) {
    ((double*)c.values)[i] = ((double*)a.values)[i] * b;
  }
  return c;
}

template <int... n>
auto operator*(const double a, const tensor<n...>& b) {
  tensor<n...> c{};
  for (int i = 0; i < (n * ...); i++) {
    ((double*)c.values)[i] = a * ((double*)b.values)[i];
  }
  return c;
}

template <int... n>
auto operator/(const tensor<n...> a, const double b) {
  tensor<n...> c{};
  for (int i = 0; i < (n * ...); i++) {
    ((double*)c.values)[i] = ((double*)a.values)[i] / b;
  }
  return c;
}

template <int... n>
auto operator/(const double a, const tensor<n...>& b) {
  tensor<n...> c{};
  for (int i = 0; i < (n * ...); i++) {
    ((double*)c.values)[i] = a / ((double*)b.values)[i];
  }
  return c;
}

template <int... n>
auto operator+(tensor<n...>& A, const tensor<n...>& B) {
  tensor<n...> C{};
  for (int i = 0; i < (n * ...); i++) {
    ((double*)C.values)[i] = ((double*)A.values)[i] + ((double*)B.values)[i];
  }
  return C;
}

template <int... n>
auto operator-(tensor<n...>& A) {
  tensor<n...> minusA{};
  for (int i = 0; i < (n * ...); i++) {
    ((double*)minusA.values)[i] = -((double*)A.values)[i];
  }
  return minusA;
}

template <int... n>
auto operator-(tensor<n...>& A, const tensor<n...>& B) {
  tensor<n...> C{};
  for (int i = 0; i < (n * ...); i++) {
    ((double*)C.values)[i] = ((double*)A.values)[i] - ((double*)B.values)[i];
  }
  return C;
}

template <int... n>
auto operator+=(tensor<n...>& A, const tensor<n...>& B) {
  for (int i = 0; i < (n * ...); i++) {
    ((double*)A.values)[i] += ((double*)B.values)[i];
  }
}

template <int... n>
auto operator-=(tensor<n...>& A, const tensor<n...>& B) {
  for (int i = 0; i < (n * ...); i++) {
    ((double*)A.values)[i] -= ((double*)B.values)[i];
  }
}

template <int n>
auto operator*(const tensor<n>& a, const tensor<n>& b) {
  double inner_product = 0.0;
  for (int i = 0; i < n; i++) {
    inner_product += a(i) * b(i);
  }
  return inner_product;
}

template <int m, int n>
auto operator*(const tensor<m, n>& A, const tensor<n>& b) {
  tensor<m> Ab{};
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      Ab(i) += A(i, j) * b(j);
    }
  }
  return Ab;
}

template <int m, int n>
auto operator*(const tensor<m>& a, const tensor<m, n>& B) {
  tensor<m> aB{};
  for (auto [i, j] : IndexSpace<n, m>{}) {
    aB(i) += a(j) * B(j, i);
  }
  return aB;
}

template <int m, int n, int p>
auto operator*(const tensor<m, n>& A, const tensor<n, p>& B) {
  tensor<m, p> AB{};
  for (auto [i, j, k] : IndexSpace<m, p, n>{}) {
    AB(i, j) += A(i, k) * B(k, j);
  }
  return AB;
}

template <int dim>
constexpr tensor<dim, dim> Identity() {
  tensor<dim, dim> I{};
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      I(i, j) = (i == j);
    }
  }
  return I;
}

template <int dim>
constexpr tensor<dim, dim> eye = Identity<dim>();

template <int dim>
constexpr tensor<dim, dim> sym(tensor<dim, dim> A) {
  tensor<dim, dim> A_sym{};
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      A_sym(i, j) = 0.5 * (A(i, j) + A(j, i));
    }
  }
  return A_sym;
}

template <int dim>
constexpr double tr(tensor<dim, dim> A) {
  double trA = 0.0;
  for (int i = 0; i < dim; i++) {
    trA += A(i, i);
  }
  return trA;
}

constexpr double det(const tensor<2, 2> A) {
  return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
}

constexpr double det(const tensor<3, 3>& A) {
  return A(0, 0) * A(1, 1) * A(2, 2) + A(0, 1) * A(1, 2) * A(2, 0) +
         A(0, 2) * A(1, 0) * A(2, 1) - A(0, 0) * A(1, 2) * A(2, 1) -
         A(0, 1) * A(1, 0) * A(2, 2) - A(0, 2) * A(1, 1) * A(2, 0);
}

template <int n>
constexpr tensor<n> linear_solve(const tensor<n, n>& A_,
                                  const tensor<n>& b_) {
  constexpr auto abs = [](double x) { return (x < 0) ? -x : x; };

  tensor<n, n> A = A_;
  tensor<n> b = b_;
  tensor<n> x{};

  for (int i = 0; i < n; i++) {
    // Search for maximum in this column
    double max_val = abs(A(i, i));

    int max_row = i;
    for (int j = i + 1; j < n; j++) {
      if (abs(A(j, i)) > max_val) {
        max_val = abs(A(j, i));
        max_row = j;
      }
    }

    // Swap maximum row with current row
    double tmp = b(max_row);
    b(max_row) = b(i);
    b(i) = tmp;
    for (int j = 0; j < n; j++) {
      tmp = A(max_row, j);
      A(max_row, j) = A(i, j);
      A(i, j) = tmp;
    }

    // zero entries below in this column
    for (int j = i + 1; j < n; j++) {
      double c = -A(j, i) / A(i, i);

      for (int k = i + 1; k < n; k++) {
        A(j, k) += c * A(i, k);
      }
      b(j) += c * b(i);
      A(j, i) = 0;
    }
  }

  // Solve equation Ax=b for an upper triangular matrix A
  for (int i = n - 1; i >= 0; i--) {
    x(i) = b(i) / A(i, i);
    for (int j = i - 1; j >= 0; j--) {
      b(j) -= A(j, i) * x(i);
    }
  }

  return x;
}

constexpr tensor<2, 2> inv(const tensor<2, 2>& A) {
  double inv_detA(1.0 / det(A));

  tensor<2, 2> invA{};

  invA(0, 0) = A(1, 1) * inv_detA;
  invA(0, 1) = -A(0, 1) * inv_detA;
  invA(1, 0) = -A(1, 0) * inv_detA;
  invA(1, 1) = A(0, 0) * inv_detA;

  return invA;
}

constexpr tensor<3, 3> inv(const tensor<3, 3>& A) {
  double inv_detA(1.0 / det(A));

  tensor<3, 3> invA{};

  invA(0, 0) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) * inv_detA;
  invA(0, 1) = (A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2)) * inv_detA;
  invA(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * inv_detA;
  invA(1, 0) = (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2)) * inv_detA;
  invA(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) * inv_detA;
  invA(1, 2) = (A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2)) * inv_detA;
  invA(2, 0) = (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)) * inv_detA;
  invA(2, 1) = (A(0, 1) * A(2, 0) - A(0, 0) * A(2, 1)) * inv_detA;
  invA(2, 2) = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0)) * inv_detA;

  return invA;
}

template <int n>
constexpr tensor<n, n> inv(const tensor<n, n>& A) {
  tensor<n, n> invA{};
  for (int j = 0; j < n; j++) {
    auto e_j = make_tensor<n>([j](int i) { return i == j; });
    auto col = linear_solve(A, e_j);
    for (int i = 0; i < n; i++) {
      invA(i, j) = col(i);
    }
  }
  return invA;
}

template <int dim>
constexpr tensor<dim, dim> dev(tensor<dim, dim> A) {
  double delta = tr(A) / dim;
  tensor<dim, dim> A_dev{};
  for (int i = 0; i < dim; i++) {
    A_dev(i, i) -= delta;
  }
  return A_dev;
}

template <int... n>
constexpr double norm(tensor<n...> u) {
  double norm_u = 0.0;
  for (int i = 0; i < (n * ...); i++) {
    norm_u += ((double*)c.values)[i] * ((double*)c.values)[i];
  }
  return sqrt(norm_u);
}

template <typename T>
constexpr auto normalize(T u) {
  return u / norm(u);
}

// vector-vector products
template <int m>
constexpr double dot(tensor<m> u, tensor<m> v) {
  double uTv = 0.0;
  for (int i = 0; i < m; i++) {
    uTv += u(i) * v(i);
  }
  return uTv;
}

template <int m, int n>
constexpr tensor<m, n> outer(tensor<m> u, tensor<n> v) {
  tensor<m, n> uvT{};
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      uvT(i, j) = u(i) * v(j);
    }
  }
  return uvT;
}

template <int m, int n, int p>
constexpr tensor<m, n, p> outer(tensor<m> u, tensor<n, p> v) {
  tensor<m, n, p> w{};
  for (auto [i, j, k] : IndexSpace<m, n, p>{}) {
    w(i, j, k) = u(i) * v(j, k);
  }
  return w;
}

template <int m, int n, int p>
constexpr tensor<m, n, p> outer(tensor<m, n> u, tensor<p> v) {
  tensor<m, n, p> w{};
  for (auto [i, j, k] : IndexSpace<m, n, p>{}) {
    w(i, j, k) = u(i, j) * v(k);
  }
  return w;
}

template <int m, int n, int p, int q>
constexpr tensor<m, n, p, q> outer(tensor<m, n> u, tensor<p, q> v) {
  tensor<m, n, p, q> w{};
  for (auto [i, j, k, l] : IndexSpace<m, n, p, q>{}) {
    w(i, j, k, l) = u(i, j) * v(k, l);
  }
  return w;
}

template < int n >
std::ostream & operator<<(std::ostream & out, tensor<n> v) {
  out << '{' << v(0);
  for (int i = 1; i < n; i++) {
    out << ", " << v(i);
  }
  out << '}';
  return out;
}

template < int m, int n >
std::ostream & operator<<(std::ostream & out, tensor<m, n> A) {
  out << '{' << '\n';
  for (int i = 0; i < m; i++) {
    out << "  {" << A(i,0);
    for (int j = 1; j < n; j++) {
      out << ", " << A(i,j);
    }
    out << '}' << '\n';
  }
  out << '}';
  return out;
}

template <int... n>
constexpr auto indices(const tensor<n...> &) {
  return IndexSpace<n...>{};
}

// adapted from https://stackoverflow.com/a/61040973
namespace impl {
template <typename>
struct is_a_tensor : public std::false_type {};

template <int... n>
struct is_a_tensor<tensor<n...>> : public std::true_type {};
}  // namespace impl

template <typename T>
using is_a_tensor = impl::is_a_tensor<std::decay_t<T>>;