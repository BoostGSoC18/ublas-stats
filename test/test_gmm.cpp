#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/gmm/gmm.hpp>

#include <boost/math/distributions/normal.hpp>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "utils.hpp"
#include "common/testhelper.hpp"

#include <string>
#include <ctime>
#include <cmath>

using namespace boost::numeric::ublas;

static const double TOL = 1e-6;

BOOST_UBLAS_TEST_DEF (test_probability) {
    BOOST_UBLAS_DEBUG_TRACE("Test the Probability () method.");  

    std::string WEIGHTS_INPUT = "[2](0.3, 0.7)\0";
    std::string MEANS_INPUT = "[2](0, 3)\0";
    std::string STDEVS_INPUT = "[2](1, 2)\0";

    vector<double> weights (2);
    std::istringstream is1 (WEIGHTS_INPUT);
    is1 >> weights;

    vector<double> means (2);
    std::istringstream is2 (MEANS_INPUT);
    is2 >> means;

    vector<double> stdevs (2);
    std::istringstream is3 (STDEVS_INPUT);
    is3 >> stdevs;

    GMM gmm (2, weights, means, stdevs);

    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (0), 0.165013842603, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (1), 0.157280970937, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (2), 0.139420154321, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (3), 0.140959352664, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (-1.3), 0.0652530921706, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (5.4), 0.067965174965, TOL);
}

BOOST_UBLAS_TEST_DEF (test_probability_component) {
    BOOST_UBLAS_DEBUG_TRACE("Test the Probability () method for each component.");  

    std::string WEIGHTS_INPUT = "[2](0.3, 0.7)\0";
    std::string MEANS_INPUT = "[2](0, 3)\0";
    std::string STDEVS_INPUT = "[2](1, 2)\0";

    vector<double> weights (2);
    std::istringstream is1 (WEIGHTS_INPUT);
    is1 >> weights;

    vector<double> means (2);
    std::istringstream is2 (MEANS_INPUT);
    is2 >> means;

    vector<double> stdevs (2);
    std::istringstream is3 (STDEVS_INPUT);
    is3 >> stdevs;

    GMM gmm (2, weights, means, stdevs);

    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (0, 0), 0.398942280401, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (0, 1), 0.0647587978329, TOL);

    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (1, 0), 0.241970724519, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (1, 1), 0.12098536226, TOL);

    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (2, 0), 0.0539909665132, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (2, 1), 0.176032663382, TOL);

    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (3, 0), 0.00443184841194, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (3, 1), 0.199471140201, TOL);

    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (-1.3, 0), 0.171368592048, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (-1.3, 1), 0.0197750207947, TOL);

    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (5.4, 0), 1.85736184456e-07, TOL);
    BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetProbability (5.4, 1), 0.0970930274916, TOL);
}

BOOST_UBLAS_TEST_DEF (test_train_single_gaussian) {
    BOOST_UBLAS_DEBUG_TRACE("Test the Train () method for a single gaussian.");  

    boost::random::mt19937 gen;
    gen.seed (static_cast<unsigned int> (std::time (0)));

    boost::random::uniform_01<double> uniform_distribution;

    GMM gmm (1);

    for (size_t iterations = 0; iterations < 4; ++ iterations) {
        double mean = uniform_distribution (gen);
        double stdev = uniform_distribution (gen);

        boost::random::normal_distribution<> normal_distribution(mean, stdev);

        vector<double> data (150 * (size_t)std::pow (10, iterations / 3.0));
        for (size_t i = 0; i < data.size (); ++ i)
            data (i) = normal_distribution (gen);

        std::cout << mean << " " << stdev << std::endl;

        gmm.Train (data);

        std::cout << gmm.GetComponentMean (0) << " " << " " << std::abs (mean - gmm.GetComponentMean (0)) / mean << std::endl;
        std::cout << gmm.GetComponentStandardDeviation (0) << " " << std::abs (stdev - gmm.GetComponentStandardDeviation (0)) / stdev << std::endl;

        std::cout << std::endl;

        BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetComponentMean (0), mean, TOL);
        BOOST_UBLAS_TEST_CHECK_CLOSE (gmm.GetComponentStandardDeviation (0), stdev, TOL);
    }
}

int main() {
    
    BOOST_UBLAS_TEST_SUITE("KMeans Test Suite");

    BOOST_UBLAS_TEST_BEGIN();
        BOOST_UBLAS_TEST_DO( test_probability );
        BOOST_UBLAS_TEST_DO( test_probability_component );
        BOOST_UBLAS_TEST_DO( test_train_single_gaussian );
        // BOOST_UBLAS_TEST_DO( test_train_multiple_gaussian );
        // BOOST_UBLAS_TEST_DO( test_sample_generator );
        // BOOST_UBLAS_TEST_DO( test_classification );
    BOOST_UBLAS_TEST_END();

    return 0;
}