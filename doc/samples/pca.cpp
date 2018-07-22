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

    // std::cout << "data " << data << std::endl;

    matrix<double> data2 (150,4);
    vector<int> gt (150);
    for (int i = 0; i < data2.size1 (); ++i) {
        for (int j = 0; j < data2.size2 (); ++j)
            std::cin >> data2 (i,j);
        std::cin >> gt (i);
    }

    PCA pca = PCA ();

    matrix<double> new_data (150,2);
    // pca.Apply (data, new_data);
    pca.Apply (data2, new_data, 2);

    // std::cout << "!@#" << std::endl;
    // std::cout << new_data << std::endl;

    // std::cout << "eigvals " << pca.GetEigenValues () << std::endl;
    // std::cout << "eigvecs " << pca.GetEigenVectors () << std::endl;
    // std::cout << "transformed data " << new_data << std::endl;

    for (int i = 0; i < data2.size1 (); i++)
        std::cout << new_data (i,0) << " " << new_data (i,1) << std::endl;

    return 0;
}