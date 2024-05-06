#define EIGEN_NO_AUTOMATIC_RESIZING 1
#define EIGEN_DONT_ALIGN 1
#define EIGEN_NO_DEBUG 1
#define EIGEN_UNROLLING_LIMIT 0
#define EIGEN_DONT_VECTORIZE 1

#include <iostream>
#include <stdexcept>

#include <Eigen/Dense>

// Enzyme voodoo
extern int enzyme_allocated
         , enzyme_const
         , enzyme_dup
         , enzyme_duponneed
         , enzyme_out
         , enzyme_tape;

template <typename Retval, typename... Args>
Retval __enzyme_autodiff(Retval (*)(Args...), auto...);

#define APPROX_EQ(LHS, RHS, TOL) {\
    if (std::abs((LHS) - (RHS)) > (TOL)) {\
        throw std::logic_error{\
            #LHS " (which is " + std::to_string(LHS) + ")"\
            " != " #RHS " (which is " + std::to_string(RHS) + ")"\
            " within " + std::to_string(TOL) + " tolerance"\
        };\
    }\
}

double square(double x) { return x * x; }

namespace eg = Eigen;

double matnorm(const eg::MatrixXd* M) {
    return M->norm();
}
double diffnorm(const eg::MatrixXd* U, const eg::MatrixXd* V) {
    auto diff = 3 * (*U) - *V;
    return diff.norm();
}

int main(int argc, char** argv) {
    const std::size_t sz = (argc < 2) ? 3 : std::stoi(argv[1]);

    // Test single-argument scalar function
    std::cout
        << "f(x)  = " << square(4)
        << '\n'
        << "f'(x) = " << __enzyme_autodiff(square, 4.0)
        << '\n';

    // Test single-argument matrix->scalar function
    eg::MatrixXd M = eg::MatrixXd::Constant(sz, sz, 1.0);
    eg::MatrixXd dfdM = eg::MatrixXd::Constant(sz, sz, 0.0);

    std::cout << "Matrix norm: " << matnorm(&M) << '\n';
    __enzyme_autodiff(matnorm, enzyme_dup, &M, &dfdM);
    std::cout << "df/dM:\n" << dfdM << '\n';

    // Test two-argument (matrix,matrix)->scalar function
    const eg::MatrixXd N = eg::MatrixXd::Constant(sz, sz, 2.0);
    eg::MatrixXd dfdN = eg::MatrixXd::Constant(sz, sz, 0.0);
    dfdM = dfdN; // IMPORTANT: Enzyme appends to buffers, so they must be zeroed

    std::cout << "Matrix difference norm: " << diffnorm(&M, &N) << '\n';
    __enzyme_autodiff(diffnorm, enzyme_dup, &M, &dfdM, enzyme_dup, &N, &dfdN);
    std::cout << "df/dM:\n" << dfdM << '\n';
    std::cout << "df/dN:\n" << dfdN << '\n';
    return 0;
}