#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/kmeans/kmeans.hpp>
#include <boost/numeric/ublas/kmeans/random_initialization.hpp>
#include <boost/numeric/ublas/kmeans/kmeans++.hpp>
#include <boost/numeric/ublas/kmeans/bradley_fayyad_refinement.hpp>
#include <boost/numeric/ublas/kmeans/naive_kmeans.hpp>

#include <string.h>
#include <iostream>

int main() {
    using namespace boost::numeric::ublas;
    std::string str_data = "[7,2]((1,1),(1.5,2),(3,4),(5,7),(3.5,5),(4.5,5),(3.5,4.5))";
    matrix<double> data;
    std::istringstream is(str_data);
    is >> data;

    matrix<double> data2 (150,2);
    vector<int> gt (150);
    for (int i = 0; i < data2.size1 (); ++i) {
        for (int j = 0; j < data2.size2 (); ++j)
            std::cin >> data2 (i,j);
        std::cin >> gt (i);
    }

    // KMeans<RandomInitialization, NaiveKMeans> kmeans(100);
    // KMeans<KMeansPlusPlus, NaiveKMeans> kmeans(100);
    KMeans<> kmeans(1000);
    // std::cout << "init!" << std::endl;
    int n = 3;
    vector<size_t> assignments (data2.size1 ());
    matrix<double> centroids (n, data2.size2());
    // std::cout << "init vars!" << std::endl;
    kmeans.Cluster (data2, n, centroids, assignments);
    // std::cout << "clusterd!" << std::endl;
    for (int i = 0; i < data2.size1 (); i++) {
        // std::cout << row (data2, i) << " " << assignments (i) << " " << gt (i) << std::endl;
        for (int j = 0; j < data2.size2 (); j++)
            std::cout << data2 (i,j) << " ";
        std::cout << assignments (i) << std::endl;
        // std::cout << gt (i) << std::endl;
    }
    return 0;
}