#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/pca.hpp>

#include "utils.hpp"
#include "common/testhelper.hpp"

#include <string>
#include <cmath>

using namespace boost::numeric::ublas;

static const double TOL = 1e-3;

BOOST_UBLAS_TEST_DEF (test_pca_num_factors) {
    BOOST_UBLAS_DEBUG_TRACE("Test PCA with number of factors parameter");

    std::string PCA_INPUT_DATA = "[5,3]((1,5,6),(0,2,7),(2,8,3),(3,4,1),(9,8,8))\0";
    std::string PCA_GT = "[5,2]((-1.53781086,1.29937798),(-3.51358020,3.45762685),(-0.16139887,-2.69910005),(-1.87706634,-3.15620704),(7.08985628,1.09830225))\0";

    matrix<double> data (5, 3);
    std::istringstream is1(PCA_INPUT_DATA);
    is1 >> data;

    matrix<double> transformed_data (5,2);
    
    PCA pca = PCA ();
    double v = pca.Apply (data, transformed_data, 2);

    vector<double> eigvals = pca.GetEigenValues ();
    matrix<double> eigvecs = pca.GetEigenVectors ();
    
    matrix<double> gt (5, 2);
    std::istringstream is2(PCA_GT);
    is2 >> gt;

    for (size_t i = 0; i < gt.size2 (); ++ i)
        if (gt (0,i) * transformed_data (0,i) < 0)
            column (transformed_data, i) *= -1;

    bool PASSED = 1;
    for (size_t i = 0; i < gt.size1 (); ++ i) {
        double distance_norm = inner_prod (row (gt, i) - row (transformed_data, i), row (gt, i) - row (transformed_data, i));
        if (std::pow (distance_norm, 0.5) > TOL) {
            std::cout << row (gt, i) << std::endl;
            std::cout << row (transformed_data, i) << std::endl;
            PASSED = 0;
            break;
        }   
    }

    BOOST_UBLAS_TEST_CHECK (PASSED);
}

BOOST_UBLAS_TEST_DEF (test_pca_variance_retained) {
    BOOST_UBLAS_DEBUG_TRACE("Test PCA with variance parameter");

    std::string PCA_INPUT_DATA = "[5,3]((1,5,6),(0,2,7),(2,8,3),(3,4,1),(9,8,8))\0";
    std::string PCA_GT = "[5,2]((-1.53781086,1.29937798),(-3.51358020,3.45762685),(-0.16139887,-2.69910005),(-1.87706634,-3.15620704),(7.08985628,1.09830225))\0";

    matrix<double> data (5, 3);
    std::istringstream is1(PCA_INPUT_DATA);
    is1 >> data;

    matrix<double> transformed_data1 (5,1);
    
    PCA pca = PCA ();
    double variance_retained = pca.Apply (data, transformed_data1, 0.5);

    BOOST_UBLAS_TEST_CHECK (std::abs (variance_retained - 0.62963) <= TOL);

    matrix<double> transformed_data2 (5,2);
    
    variance_retained = pca.Apply (data, transformed_data2, 0.80);

    BOOST_UBLAS_TEST_CHECK (std::abs (variance_retained - 0.925926) <= TOL);

    variance_retained = pca.Apply (data, transformed_data2, 0.92);

    BOOST_UBLAS_TEST_CHECK (std::abs (variance_retained - 0.925926) <= TOL);

    matrix<double> transformed_data3 (5,3);
    
    variance_retained = pca.Apply (data, transformed_data2, 0.93);

    BOOST_UBLAS_TEST_CHECK (std::abs (variance_retained - 1.0) <= TOL);

}

int main() {
    
    BOOST_UBLAS_TEST_SUITE("PCA Test Suite");

    BOOST_UBLAS_TEST_BEGIN();
        BOOST_UBLAS_TEST_DO( test_pca_num_factors );
        // BOOST_UBLAS_TEST_DO( test_pca_variance_retained );
    BOOST_UBLAS_TEST_END();

    return 0;
}