#include <boost/numeric/ublas/gmm/gmm.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <string.h>

using namespace boost::numeric::ublas;

int main(int argc, char const *argv[]) {
    
    GMM s (4);
    std::cout << s.Sample () << std::endl;

    /*
    0 - 12: Class 1
    13 - 19: Class 2
    20 - 29: Class 3
    */
    std::string GMM_INPUT_DATA = "[30,2]((0.0,0.0),(0.3,0.4),(0.1,0.0),(0.1,0.3),(-0.2,-0.2),(-0.1,0.3),(-0.4,0.1),(0.2,-0.1),(0.3,0.0),(-0.3,-0.3),(0.1,-0.1),(0.2,-0.3),(-0.3,0.2),(10.0,10.0),(10.1,9.9),(9.9,10.0),(10.2,9.7),(10.2,9.8),(9.7,10.3),(9.9,10.1),(-10.0,5.0),(-9.8,5.1),(-9.9,4.9),(-10.0,4.9),(-10.2,5.2),(-10.1,5.1),(-10.3,5.3),(-10.0,4.8),(-9.6,5.0),(-9.8,5.1))";

    matrix<double> data (30, 2);
    std::istringstream is(GMM_INPUT_DATA);
    is >> data;

    s.Train (column (data,0), 100);

    for (int i = 0; i < data.size1 (); ++ i)
        std::cout << i << " " << s.GetComponentLabel (data (i,0)) << std::endl;

    return 0;
}