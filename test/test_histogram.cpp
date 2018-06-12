#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/histogram.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "utils.hpp"
#include "common/testhelper.hpp"

using namespace boost::numeric::ublas;

static const double TOL = 1e-6;

BOOST_UBLAS_TEST_DEF (test_histogram_fixed_nbins) {
    BOOST_UBLAS_DEBUG_TRACE("Histogram with fixed number of bins");
    
    std::string VECTOR_INPUT_SET [] = {"[6](1.0,14.5,3.7,2.3,3.699,7)\0", "[10](1.0,14.5,2.349999,5.05,6,9,3.1,14.4999999,10.4600099,10)\0"};
    int vector_input_set_size = 2;

    int _fail = 0;

    vector<double> v1;
    std::istringstream is1(VECTOR_INPUT_SET [0]);
    is1 >> v1;

    // For nbins = 1, there should only be only one bin whose count = len(v).
    vector<int> hist1 = histogram (v1, 1);
    BOOST_UBLAS_TEST_CHECK (hist1.size () == 1);
    BOOST_UBLAS_TEST_CHECK (hist1 (0) == v1.size ());

    vector<double> v2;
    std::istringstream is2(VECTOR_INPUT_SET [1]);
    is2 >> v2;

    // If nbins is not given as argument, then default number of bins should be 10.
    std::string gt2 = "[10](2,1,0,2,0,1,1,1,0,2)";
    std::istringstream is3(gt2);
    vector<int> v3;
    is3 >> v3;
    vector<int> hist2 = histogram (v2);
    BOOST_UBLAS_TEST_CHECK (hist2.size () == v3.size ());
    _fail = 0;
    for (unsigned int i = 0; i < hist2.size (); ++ i)
        if (hist2 (i) != v3 (i)) {
            _fail = 1;
            break;
        }
    BOOST_UBLAS_TEST_CHECK (_fail == 0);
}

BOOST_UBLAS_TEST_DEF (test_histogram_custom_binedges) {
    BOOST_UBLAS_DEBUG_TRACE("Histogram with custom bin edges");
    
    std::string BIN_INPUT_SET [] = {"[7](1.0,2.35,5.05,6.4,10.45,11.8,14.5)\0"};
    int bin_input_set_size = 1;

    std::string VECTOR_INPUT_SET [] = {"[3](0.0,0.5,0.9)\0", "[3](14.50001,15,20)\0", "[10](1.0,14.5,2.349999,5.05,6,9,3.1,14.4999999,10.4600099,10)\0"};
    int vector_input_set_size = 3;

    int _fail = 0;

    vector<double> bin_edges;
    std::istringstream is1(BIN_INPUT_SET [0]);
    is1 >> bin_edges;

    // Bin-count vector returned should be zero filled when source vector is out of the range of the bins.
    vector<double> v1;
    std::istringstream is2(VECTOR_INPUT_SET [0]);
    is2 >> v1;
    vector<int> hist1 = histogram (v1, bin_edges);
    BOOST_UBLAS_TEST_CHECK (hist1.size () == bin_edges.size () - 1);
    _fail = 0;
    for (unsigned int i = 0; i < hist1.size (); ++ i)
        if (hist1 (i) != 0) {
            _fail = 1;
            break;
        }
    BOOST_UBLAS_TEST_CHECK (_fail == 0);
    

    vector<double> v2;
    std::istringstream is3(VECTOR_INPUT_SET [1]);
    is3 >> v2;
    vector<int> hist2 = histogram (v2, bin_edges);
    BOOST_UBLAS_TEST_CHECK (hist2.size () == bin_edges.size () - 1);
    _fail = 0;
    for (unsigned int i = 0; i < hist2.size (); ++ i)
        if (hist2 (i) != 0) {
            _fail = 1;
            break;
        }
    BOOST_UBLAS_TEST_CHECK (_fail == 0);


    // Normal custom bin-edges case involving values at bin edge borders. Check stability of bin counts.
    std::string gt1 = "[6](2,1,2,2,1,2)";
    std::istringstream is4(gt1);
    vector<int> v3;
    is4 >> v3;
    vector<double> v4;
    std::istringstream is5(VECTOR_INPUT_SET [2]);
    is5 >> v4;
    vector<int> hist3 = histogram (v4, bin_edges);
    BOOST_UBLAS_TEST_CHECK (hist3.size () == v3.size ());
    _fail = 0;
    for (unsigned int i = 0; i < hist3.size (); ++ i)
        if (hist3 (i) != v3 (i)) {
            _fail = 1;
            break;
        }
    BOOST_UBLAS_TEST_CHECK (_fail == 0);
}

int main() {
    
    BOOST_UBLAS_TEST_SUITE("Histogram Test Suite");

    BOOST_UBLAS_TEST_BEGIN();
        BOOST_UBLAS_TEST_DO( test_histogram_fixed_nbins );
        BOOST_UBLAS_TEST_DO( test_histogram_custom_binedges );
    BOOST_UBLAS_TEST_END();

    return 0;
}