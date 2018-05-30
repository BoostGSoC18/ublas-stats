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

#include <boost/numeric/ublas/io.hpp>

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

    /*
    3 1 2
    1 2 1
    2 3 4
    1 1 1 2 2 2 3 3 4
    */

    // auto p = m.find1(0, 0, 0);
    // while (p != m.end1()) {
    //     auto q = p.begin();
    //     while (q != p.end()) {
    //         std::cout << *q << " ";
    //         ++ q;
    //     }
    //     std::cout << std::endl;
    //     ++ p;
    //     break;
    // }

    // vector<int> x (3);
    // x (0) = 1;
    // x (1) = 2;
    // x (2) = 3;

    // matrix< std::complex<int> > m1 (3, 3);
    // for (unsigned i = 0; i < m.size1 (); ++ i) {
    //     for (unsigned j = 0; j < m.size2 (); ++ j) {
    //         m (i, j) = std::complex<int>(i, j);
    //     }
    // }

    // vector< std::complex<int> > v1 = mean (m1, 0);


    // std::cout << min (m) << std::endl;
    // for (unsigned i = 0; i < 3; ++ i) {
    //     std::cout << min (m, 1) (i) << std::endl;
    // }
    // std::cout << max (m) << std::endl;
    // for (unsigned i = 0; i < 3; ++ i) {
    //     std::cout << max (m, 1) (i) << std::endl;
    // }

    // std::cout << sum (m) << std::endl;
    // for (auto p : sum (m, 1)) {
    //     std::cout << p << std::endl;    
    // }
    // auto p = sum (m, 1).find (0);
    // while (p != sum (m, 1).end ()) {
    //     std::cout << *p << std::endl;
    //     ++ p;
    // }
    // std::cout << sum (m, 0) << std::endl;
    // for (unsigned i = 0; i < 3; ++ i) {
    //     std::cout << sum (m, 0) (i) << std::endl;
    // }

    // std::cout << mean (m) << std::endl;
    // std::cout << mean_iterative (m) << std::endl;
    
    // vector<double> p = mean (m, 0);
    // std::cout << p << std::endl;
    // for (unsigned i = 0; i < 3; ++ i) {
    //     std::cout << mean (m, 1) (i) << std::endl;
    // }

    // for (auto &p : mean (m, 1)) {
    //     std::cout << p << std::endl;
    // }
    
    // std::cout << variance (m) << std::endl;
    // for (unsigned i = 0; i < 3; ++ i) {
    //     std::cout << variance (m, 0) (i) << std::endl;
    // }

    std::cout << median (m) << std::endl;
    std::cout << median (m, 1) << std::endl;
    // std::cout << mode (m) << std::endl;
    // std::cout << mode (m, 0) << std::endl;
    
    // std::cout << variance_iterative (m) << std::endl;
    std::cout << norm_1 (m) << std::endl;
    std::cout << norm_frobenius (m) << std::endl;
    std::cout << norm_inf (m) << std::endl;
}

