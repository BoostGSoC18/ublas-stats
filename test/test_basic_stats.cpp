#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/covariance_matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "utils.hpp"
#include "common/testhelper.hpp"

using namespace boost::numeric::ublas;

static const double TOL = 1e-6;

std::string VECTOR_INPUT_SET [] = {"[3](1,4,1)\0", "[3](1.023, 2.349, 7.839)\0"};
int vector_input_set_size = 2;

std::string MATRIX_INPUT_SET [] = {"[3,3]((2,3,1),(2,3,4),(4,3,5))\0", "[2,4]((1,2,0,3),(2,1.99,0,2.0001))\0"};
int matrix_input_set_size = 2;

/*
template<class Pred>
void test_stat_method (std::string input [], int input_size, std::string output [], Pred stat_method) {

}

template<class Pred>
void test_stat_axis_method (std::string input [], int input_size, std::string output [], Pred stat_axis_method, int axis) {

}
*/

BOOST_UBLAS_TEST_DEF (test_min) {
    BOOST_UBLAS_DEBUG_TRACE("Vector Min");

    double VECTOR_MIN_RESULT_SET [] = {1, 1.023};
    for (int i = 0; i < vector_input_set_size; ++ i) {
        // std::cout << "Test vector " << i << ": ";
        vector<double> v;
        std::istringstream is(VECTOR_INPUT_SET [i]);
        is >> v;
        // assertTrue ("", min (v) == VECTOR_MIN_RESULT_SET [i]);
        BOOST_UBLAS_TEST_CHECK (min (v) == VECTOR_MIN_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Min");

    double MATRIX_MIN_RESULT_SET [] = {1, 0};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << ": ";
        matrix<double> m;
        std::istringstream is(MATRIX_INPUT_SET [i]);
        is >> m;
        BOOST_UBLAS_TEST_CHECK (min (m) == MATRIX_MIN_RESULT_SET [i]);
        // assertTrue ("", min (m) == MATRIX_MIN_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Min Axis");

    std::string MATRIX_MIN_AXIS_0_RESULT_SET [] = {"[3](2,3,1)\0", "[4](1,1.99,0,2.0001)\0"};
    std::string MATRIX_MIN_AXIS_1_RESULT_SET [] = {"[3](1,2,3)\0","[2](0,0)\0"};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << std::endl;

        int _fail;
        matrix<double> m;
        std::istringstream is1(MATRIX_INPUT_SET [i]);
        is1 >> m;

        vector<double> gt0;
        std::istringstream is2(MATRIX_MIN_AXIS_0_RESULT_SET [i]);
        is2 >> gt0;
        
        _fail = 0;
        vector<double> resv0 = min (m, 0);
        for (unsigned int j = 0; j < resv0.size (); ++ j)
            if (resv0 (j) != gt0 (j)) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 0: ", _fail == 0);
        
        _fail = 0;
        vector<double> gt1;
        std::istringstream is3(MATRIX_MIN_AXIS_1_RESULT_SET [i]);
        is3 >> gt1;
        
        vector<double> resv1 = min (m, 1);
        for (unsigned int j = 0; j < resv1.size (); ++ j)
            if (resv1 (j) != gt1 (j)) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 1: ", _fail == 0);
    }
}

BOOST_UBLAS_TEST_DEF (test_max) {
    BOOST_UBLAS_DEBUG_TRACE("Vector Max");

    double VECTOR_MAX_RESULT_SET [] = {4, 7.839};
    for (int i = 0; i < vector_input_set_size; ++ i) {
        // std::cout << "Test vector " << i << ": ";
        vector<double> v;
        std::istringstream is(VECTOR_INPUT_SET [i]);
        is >> v;
        BOOST_UBLAS_TEST_CHECK (max (v) == VECTOR_MAX_RESULT_SET [i])
        // assertTrue ("", max (v) == VECTOR_MAX_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Max");

    double MATRIX_MAX_RESULT_SET [] = {5, 3};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << ": ";
        matrix<double> m;
        std::istringstream is(MATRIX_INPUT_SET [i]);
        is >> m;
        BOOST_UBLAS_TEST_CHECK (max (m) == MATRIX_MAX_RESULT_SET [i]);
        // assertTrue ("", max (m) == MATRIX_MAX_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Max Axis");

    std::string MATRIX_MAX_AXIS_0_RESULT_SET [] = {"[3](4,3,5)\0", "[4](2,2,0,3)\0"};
    std::string MATRIX_MAX_AXIS_1_RESULT_SET [] = {"[3](3,4,5)\0","[2](3,2.0001)\0"};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << std::endl;

        int _fail;
        matrix<double> m;
        std::istringstream is1(MATRIX_INPUT_SET [i]);
        is1 >> m;

        vector<double> gt0;
        std::istringstream is2(MATRIX_MAX_AXIS_0_RESULT_SET [i]);
        is2 >> gt0;
        
        _fail = 0;
        vector<double> resv0 = max (m, 0);
        for (unsigned int j = 0; j < resv0.size (); ++ j)
            if (resv0 (j) != gt0 (j)) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 0: ", _fail == 0);
        
        _fail = 0;
        vector<double> gt1;
        std::istringstream is3(MATRIX_MAX_AXIS_1_RESULT_SET [i]);
        is3 >> gt1;
        
        vector<double> resv1 = max (m, 1);
        for (unsigned int j = 0; j < resv1.size (); ++ j)
            if (resv1 (j) != gt1 (j)) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 1: ", _fail == 0);
    }
}

BOOST_UBLAS_TEST_DEF (test_sum) {
    BOOST_UBLAS_DEBUG_TRACE("Vector Sum");

    double VECTOR_SUM_RESULT_SET [] = {6, 11.211};
    for (int i = 0; i < vector_input_set_size; ++ i) {
        // std::cout << "Test vector " << i << ": ";
        vector<double> v;
        std::istringstream is(VECTOR_INPUT_SET [i]);
        is >> v;
        BOOST_UBLAS_TEST_CHECK (std::abs(sum (v) - VECTOR_SUM_RESULT_SET [i]) <= TOL);
        // assertTrue ("", sum (v) == VECTOR_SUM_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Sum");

    double MATRIX_SUM_RESULT_SET [] = {27, 11.9901};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << ": ";
        matrix<double> m;
        std::istringstream is(MATRIX_INPUT_SET [i]);
        is >> m;
        BOOST_UBLAS_TEST_CHECK (std::abs(sum (m) - MATRIX_SUM_RESULT_SET [i]) <= TOL);
        // assertTrue ("", sum (m) == MATRIX_SUM_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Sum Axis");

    std::string MATRIX_SUM_AXIS_0_RESULT_SET [] = {"[3](8,9,10)\0", "[4](3,3.99,0,5.0001)\0"};
    std::string MATRIX_SUM_AXIS_1_RESULT_SET [] = {"[3](6,9,12)\0","[2](6, 5.9901)\0"};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << std::endl;

        int _fail;
        matrix<double> m;
        std::istringstream is1(MATRIX_INPUT_SET [i]);
        is1 >> m;

        vector<double> gt0;
        std::istringstream is2(MATRIX_SUM_AXIS_0_RESULT_SET [i]);
        is2 >> gt0;
        
        _fail = 0;
        vector<double> resv0 = sum (m, 0);
        for (unsigned int j = 0; j < resv0.size (); ++ j)
            if (std::abs(resv0 (j) - gt0 (j)) > TOL) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 0: ", _fail == 0);
        
        _fail = 0;
        vector<double> gt1;
        std::istringstream is3(MATRIX_SUM_AXIS_1_RESULT_SET [i]);
        is3 >> gt1;
        
        vector<double> resv1 = sum (m, 1);
        for (unsigned int j = 0; j < resv1.size (); ++ j)
            if (std::abs (resv1 (j) - gt1 (j)) > TOL) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 1: ", _fail == 0);
    }
}

BOOST_UBLAS_TEST_DEF (test_mean) {
    BOOST_UBLAS_DEBUG_TRACE("Vector Mean");

    double VECTOR_MEAN_RESULT_SET [] = {2, 3.737};
    for (int i = 0; i < vector_input_set_size; ++ i) {
        // std::cout << "Test vector " << i << ": ";
        vector<double> v;
        std::istringstream is(VECTOR_INPUT_SET [i]);
        is >> v;
        BOOST_UBLAS_TEST_CHECK (std::abs(mean (v) - VECTOR_MEAN_RESULT_SET [i]) <= TOL);
        // assertTrue ("", mean (v) == VECTOR_MEAN_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Mean");

    double MATRIX_MEAN_RESULT_SET [] = {3, 1.4987625};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << ": ";
        matrix<double> m;
        std::istringstream is(MATRIX_INPUT_SET [i]);
        is >> m;
        BOOST_UBLAS_TEST_CHECK (std::abs (mean (m) - MATRIX_MEAN_RESULT_SET [i]) <= TOL);
        // assertTrue ("", mean (m) == MATRIX_MEAN_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Mean Axis");

    std::string MATRIX_MEAN_AXIS_0_RESULT_SET [] = {"[3](2.66666667,3,3.3333333)\0", "[4](1.5,1.995,0,2.50005)\0"};
    std::string MATRIX_MEAN_AXIS_1_RESULT_SET [] = {"[3](2,3,4)\0","[2](1.5, 1.497525)\0"};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << std::endl;

        int _fail;
        matrix<double> m;
        std::istringstream is1 (MATRIX_INPUT_SET [i]);
        is1 >> m;

        vector<double> gt0;
        std::istringstream is2 (MATRIX_MEAN_AXIS_0_RESULT_SET [i]);
        is2 >> gt0;
        
        _fail = 0;
        vector<double> resv0 = mean (m, 0);
        for (unsigned int j = 0; j < resv0.size (); ++ j)
            if (std::abs (resv0 (j) - gt0 (j)) > TOL) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 0: ", _fail == 0);
        
        _fail = 0;
        vector<double> gt1;
        std::istringstream is3 (MATRIX_MEAN_AXIS_1_RESULT_SET [i]);
        is3 >> gt1;
        
        vector<double> resv1 = mean (m, 1);
        for (unsigned int j = 0; j < resv1.size (); ++ j)
            if (std::abs (resv1 (j) - gt1 (j)) > TOL) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 1: ", _fail == 0);
    }
}

BOOST_UBLAS_TEST_DEF (test_mode) {
    BOOST_UBLAS_DEBUG_TRACE("Vector Mode");

    double VECTOR_MODE_RESULT_SET [] = {1, 1.023};
    for (int i = 0; i < vector_input_set_size; ++ i) {
        // std::cout << "Test vector " << i << ": ";
        vector<double> v;
        std::istringstream is(VECTOR_INPUT_SET [i]);
        is >> v;
        BOOST_UBLAS_TEST_CHECK (mode (v) == VECTOR_MODE_RESULT_SET [i]);
        // assertTrue ("", mode (v) == VECTOR_MODE_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Mode");

    double MATRIX_MODE_RESULT_SET [] = {3, 0};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << ": ";
        matrix<double> m;
        std::istringstream is(MATRIX_INPUT_SET [i]);
        is >> m;
        BOOST_UBLAS_TEST_CHECK (mode (m) == MATRIX_MODE_RESULT_SET [i]);
        // assertTrue ("", mode (m) == MATRIX_MODE_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Mode Axis");

    std::string MATRIX_MODE_AXIS_0_RESULT_SET [] = {"[3](2,3,1)\0", "[4](1,1.99,0,2.0001)\0"};
    std::string MATRIX_MODE_AXIS_1_RESULT_SET [] = {"[3](1,2,3)\0","[2](0,0)\0"};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << std::endl;

        int _fail;
        matrix<double> m;
        std::istringstream is1 (MATRIX_INPUT_SET [i]);
        is1 >> m;

        vector<double> gt0;
        std::istringstream is2 (MATRIX_MODE_AXIS_0_RESULT_SET [i]);
        is2 >> gt0;
        
        _fail = 0;
        vector<double> resv0 = mode (m, 0);
        for (unsigned int j = 0; j < resv0.size (); ++ j)
            if (resv0 (j) != gt0 (j)) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 0: ", _fail == 0);
        
        _fail = 0;
        vector<double> gt1;
        std::istringstream is3 (MATRIX_MODE_AXIS_1_RESULT_SET [i]);
        is3 >> gt1;
        
        vector<double> resv1 = mode (m, 1);
        for (unsigned int j = 0; j < resv1.size (); ++ j)
            if (resv1 (j) != gt1 (j)) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 1: ", _fail == 0);
    }
}

BOOST_UBLAS_TEST_DEF (test_median) {
    BOOST_UBLAS_DEBUG_TRACE("Vector Median");

    double VECTOR_MEDIAN_RESULT_SET [] = {1, 2.349};
    for (int i = 0; i < vector_input_set_size; ++ i) {
        // std::cout << "Test vector " << i << ": ";
        vector<double> v;
        std::istringstream is(VECTOR_INPUT_SET [i]);
        is >> v;
        BOOST_UBLAS_TEST_CHECK (median (v) == VECTOR_MEDIAN_RESULT_SET [i]);
        // assertTrue ("", median (v) == VECTOR_MEDIAN_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Median");

    double MATRIX_MEDIAN_RESULT_SET [] = {3, 1.995};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << ": ";
        matrix<double> m;
        std::istringstream is(MATRIX_INPUT_SET [i]);
        is >> m;
        BOOST_UBLAS_TEST_CHECK (median (m) == MATRIX_MEDIAN_RESULT_SET [i]);
        // assertTrue ("", median (m) == MATRIX_MEDIAN_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Median Axis");

    std::string MATRIX_MEDIAN_AXIS_0_RESULT_SET [] = {"[3](2,3,4)\0", "[4](1.5,1.995,0,2.50005)\0"};
    std::string MATRIX_MEDIAN_AXIS_1_RESULT_SET [] = {"[3](2,3,4)\0","[2](1.5,1.995)\0"};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << std::endl;

        int _fail;
        matrix<double> m;
        std::istringstream is1 (MATRIX_INPUT_SET [i]);
        is1 >> m;

        vector<double> gt0;
        std::istringstream is2 (MATRIX_MEDIAN_AXIS_0_RESULT_SET [i]);
        is2 >> gt0;
        
        _fail = 0;
        vector<double> resv0 = median (m, 0);
        for (unsigned int j = 0; j < resv0.size (); ++ j)
            if (resv0 (j) != gt0 (j)) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 0: ", _fail == 0);
        
        _fail = 0;
        vector<double> gt1;
        std::istringstream is3 (MATRIX_MEDIAN_AXIS_1_RESULT_SET [i]);
        is3 >> gt1;
        
        vector<double> resv1 = median (m, 1);
        for (unsigned int j = 0; j < resv1.size (); ++ j)
            if (resv1 (j) != gt1 (j)) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 1: ", _fail == 0);
    }
}

BOOST_UBLAS_TEST_DEF (test_variance) {
    BOOST_UBLAS_DEBUG_TRACE("Vector Variance");

    double VECTOR_VARIANCE_RESULT_SET [] = {2, 8.706248000};
    for (int i = 0; i < vector_input_set_size; ++ i) {
        // std::cout << "Test vector " << i << ": ";
        vector<double> v;
        std::istringstream is(VECTOR_INPUT_SET [i]);
        is >> v;
        BOOST_UBLAS_TEST_CHECK (std::abs(variance (v) - VECTOR_VARIANCE_RESULT_SET [i]) <= TOL);
        // assertTrue ("", variance (v) == VECTOR_VARIANCE_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Variance");

    double MATRIX_VARIANCE_RESULT_SET [] = {1.33333333, 0.99877346984375};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << ": ";
        matrix<double> m;
        std::istringstream is(MATRIX_INPUT_SET [i]);
        is >> m;
        BOOST_UBLAS_TEST_CHECK (std::abs (variance (m) - MATRIX_VARIANCE_RESULT_SET [i]) <= TOL);
        // assertTrue ("", variance (m) == MATRIX_VARIANCE_RESULT_SET [i]);
    }

    BOOST_UBLAS_DEBUG_TRACE("Matrix Variance Axis");

    std::string MATRIX_VARIANCE_AXIS_0_RESULT_SET [] = {"[3](0.88888889,0,2.88888889)\0", "[4](0.25,0.000025,0,0.24995)\0"};
    std::string MATRIX_VARIANCE_AXIS_1_RESULT_SET [] = {"[3](0.66666667,0.66666667,0.66666667)\0","[2](1.25, 0.74754388)\0"};
    for (int i = 0; i < matrix_input_set_size; ++ i) {
        // std::cout << "Test matrix " << i << std::endl;

        int _fail;
        matrix<double> m;
        std::istringstream is1 (MATRIX_INPUT_SET [i]);
        is1 >> m;

        vector<double> gt0;
        std::istringstream is2 (MATRIX_VARIANCE_AXIS_0_RESULT_SET [i]);
        is2 >> gt0;
        
        _fail = 0;
        vector<double> resv0 = variance (m, 0);
        for (unsigned int j = 0; j < resv0.size (); ++ j)
            if (std::abs (resv0 (j) - gt0 (j)) > TOL) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 0: ", _fail == 0);
        
        _fail = 0;
        vector<double> gt1;
        std::istringstream is3 (MATRIX_VARIANCE_AXIS_1_RESULT_SET [i]);
        is3 >> gt1;
        
        vector<double> resv1 = variance (m, 1);
        for (unsigned int j = 0; j < resv1.size (); ++ j)
            if (std::abs (resv1 (j) - gt1 (j)) > TOL) {
                _fail = 1;
                break;
            }
        BOOST_UBLAS_TEST_CHECK (_fail == 0);
        // assertTrue ("axis 1: ", _fail == 0);
    }
}

/*
std::string VECTOR_INPUT_SET [] = {"[3](1,4,1)\0", "[3](1.023, 2.349, 7.839)\0"};
int vector_input_set_size = 2;

std::string MATRIX_INPUT_SET [] = {"[3,3]((2,3,1),(2,3,4),(4,3,5))\0", "[2,4]((1,2,0,3),(2,1.99,0,2.0001))\0"};
int matrix_input_set_size = 2;
*/
BOOST_UBLAS_TEST_DEF (test_vector_covariance) {
    BOOST_UBLAS_DEBUG_TRACE("Vector Covariance");

    vector<double> v1;
    std::istringstream is1 (VECTOR_INPUT_SET [0]);
    is1 >> v1;

    vector<double> v2;
    std::istringstream is2 (VECTOR_INPUT_SET [1]);
    is2 >> v2;
    
    double cov_gt = -1.388;
    BOOST_UBLAS_TEST_CHECK (std::abs (covariance (v1, v2) - cov_gt) <= TOL);

    // Covariance of identical vectors is just the variance.
    BOOST_UBLAS_TEST_CHECK (std::abs (covariance (v1, v1) - variance (v1)) <= TOL );
    BOOST_UBLAS_TEST_CHECK (std::abs (covariance (v2, v2) - variance (v2)) <= TOL );
}

BOOST_UBLAS_TEST_DEF (test_covariance_matrix) {
    BOOST_UBLAS_DEBUG_TRACE("Covariance Matrix");

    std::string COVARIANCE_MATRIX_ROWVAR_0_RESULT_SET [] = {"[3,3]((0.88888889,0,1.11111111),(0,0,0),(1.11111111,0,2.88888889))","[4,4]((0.25,-0.0025,0,-0.249975),(-0.0025,0.000025,0,0.00249975),(0,0,0,0),(-0.249975,0.00249975,0,0.24995))"};
    std::string COVARIANCE_MATRIX_ROWVAR_1_RESULT_SET [] = {"[3,3]((0.66666667,-0.33333333,-0.66666667),(-0.33333333,0.66666667,0.33333333),(-0.66666667,0.33333333,0.66666667))","[2,2]((1.25,0.7487875),(0.7487875,0.74754388))"};
    for (int rvar = 0; rvar < 2; ++ rvar)
        for (int i = 0; i < matrix_input_set_size; ++ i) {
            // std::cout << "Test matrix " << i << ": ";
            
            int _fail;
            matrix<double> m;
            std::istringstream is(MATRIX_INPUT_SET [i]);
            is >> m;

            matrix<double> cov = covariance_matrix (m, rvar);
            // std::cout << cov << std::endl;
            
            // Cov should be symmetric.
            BOOST_UBLAS_TEST_CHECK (cov.size1 () == cov.size2 ());
            for (unsigned int j = 0; j < cov.size1 (); ++ j)
                for (unsigned int k = 0; k <= j; ++ k)
                    BOOST_UBLAS_TEST_CHECK (cov (j,k) == cov (k,j));

            if (rvar == 0) {
                BOOST_UBLAS_TEST_CHECK (cov.size1 () == m.size2 ());
                matrix<double> gt0;
                std::istringstream is0 (COVARIANCE_MATRIX_ROWVAR_0_RESULT_SET [i]);
                is0 >> gt0;
                _fail = 0;
                for (unsigned int j = 0; j < cov.size1 (); ++ j)
                    for (unsigned int k = 0; k <= j; ++ k)
                        if (std::abs (cov (j,k) - gt0 (j,k)) > TOL) {
                            _fail = 1;
                            break;
                        }
                BOOST_UBLAS_TEST_CHECK (_fail == 0);
            }
            else if (rvar == 1) {
                BOOST_UBLAS_TEST_CHECK (cov.size1 () == m.size1 ());
                matrix<double> gt1;
                std::istringstream is1 (COVARIANCE_MATRIX_ROWVAR_1_RESULT_SET [i]);
                is1 >> gt1;
                _fail = 0;
                for (unsigned int j = 0; j < cov.size1 (); ++ j)
                    for (unsigned int k = 0; k <= j; ++ k)
                        if (std::abs (cov (j,k) - gt1 (j,k)) > TOL) {
                            _fail = 1;
                            break;
                        }
                BOOST_UBLAS_TEST_CHECK (_fail == 0);
            }
        }
}

int main() {
    
    BOOST_UBLAS_TEST_SUITE("Basic Stats Test Suite");

    BOOST_UBLAS_TEST_BEGIN();
        BOOST_UBLAS_TEST_DO( test_min );
        BOOST_UBLAS_TEST_DO( test_max );
        BOOST_UBLAS_TEST_DO( test_sum );
        BOOST_UBLAS_TEST_DO( test_mean );
        BOOST_UBLAS_TEST_DO( test_mode );
        BOOST_UBLAS_TEST_DO( test_median );
        BOOST_UBLAS_TEST_DO( test_variance );
        BOOST_UBLAS_TEST_DO( test_vector_covariance );
        BOOST_UBLAS_TEST_DO( test_covariance_matrix );
    BOOST_UBLAS_TEST_END();

    return 0;
}