#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/kmeans/kmeans.hpp>
#include <boost/numeric/ublas/kmeans/random_initialization.hpp>
#include <boost/numeric/ublas/kmeans/kmeans++.hpp>
#include <boost/numeric/ublas/kmeans/bradley_fayyad_refinement.hpp>
#include <boost/numeric/ublas/kmeans/naive_kmeans.hpp>

#include "utils.hpp"
#include "common/testhelper.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

#include <string>
#include <ctime>
#include <cmath>

using namespace boost::numeric::ublas;

static const double TOL = 1e-6;

/*
0 - 12: Class 1
13 - 19: Class 2
20 - 29: Class 3
*/
std::string KMEANS_INPUT_DATA = "[30,2]((0.0,0.0),(0.3,0.4),(0.1,0.0),(0.1,0.3),(-0.2,-0.2),(-0.1,0.3),(-0.4,0.1),(0.2,-0.1),(0.3,0.0),(-0.3,-0.3),(0.1,-0.1),(0.2,-0.3),(-0.3,0.2),(10.0,10.0),(10.1,9.9),(9.9,10.0),(10.2,9.7),(10.2,9.8),(9.7,10.3),(9.9,10.1),(-10.0,5.0),(-9.8,5.1),(-9.9,4.9),(-10.0,4.9),(-10.2,5.2),(-10.1,5.1),(-10.3,5.3),(-10.0,4.8),(-9.6,5.0),(-9.8,5.1))";

BOOST_UBLAS_TEST_DEF (test_random_initialization) {
    BOOST_UBLAS_DEBUG_TRACE("Test random initialization");

    matrix<double> data (30, 2);
    std::istringstream is(KMEANS_INPUT_DATA);
    is >> data;

    matrix<double> centroids (3, data.size2 ());

    RandomInitialization ri;
    ri.Initialize (data, 3, centroids);

    for (size_t i = 0; i < centroids.size1 (); ++ i) {
        bool found = 0;
        for (size_t j = 0; j < data.size1 (); ++ j) {
            double distance_norm = inner_prod (row (centroids, i) - row (data, j), row (centroids, i) - row (data, j));
            if (std::pow (distance_norm, 0.5) < TOL) {
                found = 1;
                break;
            }
        }
        BOOST_UBLAS_TEST_CHECK (found);
    }
    
}

BOOST_UBLAS_TEST_DEF (test_kmeans_basic) {
    BOOST_UBLAS_DEBUG_TRACE("Basic kmeans test");

    matrix<double> data (30, 2);
    std::istringstream is (KMEANS_INPUT_DATA);
    is >> data;

    KMeans<> kmeans (1000);

    vector<size_t> assignments (data.size1 ());

    kmeans.Cluster (data, 3, assignments);

    size_t first_class = assignments (0);
    for (size_t i = 0; i < 13; ++ i)
        BOOST_UBLAS_TEST_CHECK (assignments (i) == first_class);

    size_t second_class = assignments (13);
    BOOST_UBLAS_TEST_CHECK (first_class != second_class);
    for (size_t i = 13; i < 20; ++ i)
        BOOST_UBLAS_TEST_CHECK (assignments (i) == second_class);

    size_t third_class = assignments (20);
    BOOST_UBLAS_TEST_CHECK (first_class != third_class);
    BOOST_UBLAS_TEST_CHECK (second_class != third_class);
    for (size_t i = 20; i < assignments.size (); ++ i)
        BOOST_UBLAS_TEST_CHECK (assignments (i) == third_class);
}

BOOST_UBLAS_TEST_DEF (test_refined_start) {
    BOOST_UBLAS_DEBUG_TRACE("Bradley-Fayyad refined start test");

    boost::random::mt19937 gen;
    gen.seed(static_cast<unsigned int>(std::time(0)));
    boost::random::normal_distribution<> normal_dist (0, 1);

    matrix<double> data (3000, 3);
    for (size_t i = 0; i < data.size1 (); ++ i)
        for (size_t j = 0; j < data.size2 (); ++ j)
            data (i, j) = normal_dist (gen);

    std::string centroids_input = "[5,3]((0,0,0),(5,0,-2),(-2,-2,-2),(-6,8,8),(1,6,1))";
    matrix<double> centroids (5,3);
    std::istringstream is (centroids_input);
    is >> centroids;

    // First Gaussian: 10000 points, centered at (0, 0, 0).
    // Second Gaussian: 2000 points, centered at (5, 0, -2).
    // Third Gaussian: 5000 points, centered at (-2, -2, -2).
    // Fourth Gaussian: 1000 points, centered at (-6, 8, 8).
    // Fifth Gaussian: 12000 points, centered at (1, 6, 1).
    for (size_t i = 1000; i < 1200; ++i) {
    // std::cout << row (data, i) << std::endl;
    // std::cout << row (centroids, 1) << std::endl;
        row (data, i) += row (centroids, 1);
    }
    for (size_t i = 1200; i < 1700; ++i)
        row (data, i) += row (centroids, 2);
    for (size_t i = 1700; i < 1800; ++i)
        row (data, i) += row (centroids, 3);
    for (size_t i = 1800; i < 3000; ++i)
        row (data, i) += row (centroids, 4);

    size_t num_clusters = 5;


    KMeans<RefinedStart> kmeans (0, 1);
    vector<size_t> assignments (data.size1 ());
    matrix<double> resultingCentroids (num_clusters, data.size2 ());
    kmeans.Cluster (data, num_clusters, resultingCentroids, assignments);

    // Calculate sum of distances from centroid means.
    double distortion = 0;
    for (size_t i = 0; i < data.size1 (); ++ i)
        distortion += inner_prod (row (data, i) - row (resultingCentroids, assignments[i]), row (data, i) - row (resultingCentroids, assignments[i]));

    std::cout << distortion << std::endl;

    // Using the refined start, the distance for this dataset is usually around
    // 13500.  Regular k-means is between 10000 and 30000 (I think the 10000
    // figure is a corner case which actually does not give good clusters), and
    // random initial starts give distortion around 22000.  So we'll require that
    // our distortion is less than 14000.
    BOOST_UBLAS_TEST_CHECK (distortion <= 14000.0);
}

BOOST_UBLAS_TEST_DEF (test_kmeans_plusplus) {
    BOOST_UBLAS_DEBUG_TRACE("KMeans++ test");

    std::string INPUT_DATA = "[4,2]((0,0),(4,0),(4,2),(0,2))";

    matrix<double> data (4, 2);
    std::istringstream is (INPUT_DATA);
    is >> data;

    KMeans<KMeansPlusPlus> kmeans;

    vector<size_t> assignments (data.size1 ());

    size_t num_clusters = 2;

    kmeans.Cluster (data, num_clusters, assignments);

    std::cout << assignments << std::endl;

    BOOST_UBLAS_TEST_CHECK (assignments[0] == assignments[3]);
    BOOST_UBLAS_TEST_CHECK (assignments[1] == assignments[2]);
}

int main() {
    
    BOOST_UBLAS_TEST_SUITE("KMeans Test Suite");

    BOOST_UBLAS_TEST_BEGIN();
        BOOST_UBLAS_TEST_DO( test_random_initialization );
        BOOST_UBLAS_TEST_DO( test_kmeans_basic );
        BOOST_UBLAS_TEST_DO( test_refined_start );
        BOOST_UBLAS_TEST_DO( test_kmeans_plusplus );
    BOOST_UBLAS_TEST_END();

    return 0;
}