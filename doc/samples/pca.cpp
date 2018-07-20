#include <boost/numeric/ublas/pca.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <string.h>

using namespace boost::numeric::ublas;

int main(int argc, char const *argv[]) {
    
    std::string PCA_INPUT = "[3,2]((1,2),(3,4),(5,6))\0";

    matrix<double> data (3,2);
    std::istringstream is(PCA_INPUT);
    is >> data;

    std::cout << "data " << data << std::endl;

    PCA pca = PCA ();

    matrix<double> new_data (3,2);
    vector<double> eigvals (2);
    matrix<double> eigvecs (2,2);
    pca.Apply (data, new_data, eigvals, eigvecs);
    // pca.Apply (data, new_data, 1);

    std::cout << "eigvals " << pca.GetEigenValues () << std::endl;
    std::cout << "eigvecs " << pca.GetEigenVectors () << std::endl;
    std::cout << "transformed data " << new_data << std::endl;

    return 0;
}