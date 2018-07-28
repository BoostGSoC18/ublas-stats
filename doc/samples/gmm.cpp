#include <boost/numeric/ublas/gmm/gmm.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <string.h>

using namespace boost::numeric::ublas;

int main(int argc, char const *argv[]) {
    
    GMM s (4);
    std::cout << s.Sample () << std::endl;

    return 0;
}