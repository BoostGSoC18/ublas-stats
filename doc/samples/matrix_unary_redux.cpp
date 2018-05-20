//
//  Copyright (c) 2000-2002
//  Joerg Walter, Mathias Koch
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.
//

#include <boost/numeric/ublas/matrix.hpp>

int main () {
    using namespace boost::numeric::ublas;
    matrix<int> m (3, 3);
    for (unsigned i = 0; i < m.size1 (); ++ i) {
        for (unsigned j = 0; j < m.size2 (); ++ j) {
            m (i, j) = i+j;
        }
    }

    m (0, 0) = 3;
    m (1, 2) = 1;

    // std::cout << min (m) << std::endl;
    // for (unsigned i = 0; i < 3; ++ i) {
    //     std::cout << min (m, 1) (i) << std::endl;
    // }
    // std::cout << max (m) << std::endl;
    // for (unsigned i = 0; i < 3; ++ i) {
    //     std::cout << max (m, 1) (i) << std::endl;
    // }

    std::cout << sum (m) << std::endl;
    for (unsigned i = 0; i < 3; ++ i) {
        std::cout << sum (m, 0) (i) << std::endl;
    }

    // std::cout << mean (m) << std::endl;
    // std::cout << mean_iterative (m) << std::endl;
    
    // vector<double> tmp = mean (m, 0);
    // for (unsigned i = 0; i < 3; ++ i) {
    //     std::cout << mean (m, 1) (i) << std::endl;
    // }
    
    // std::cout << median (m) << std::endl;
    // std::cout << mode (m) << std::endl;
    // std::cout << variance (m) << std::endl;
    // std::cout << variance_iterative (m) << std::endl;
    std::cout << norm_1 (m) << std::endl;
    std::cout << norm_frobenius (m) << std::endl;
    std::cout << norm_inf (m) << std::endl;
}

