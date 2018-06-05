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

#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/histogram.hpp>

int main () {
    using namespace boost::numeric::ublas;
    vector<int> v (11);
    for (unsigned i = 0; i < v.size (); ++ i)
        v (i) = i;
    // v (0) = 1;
    // v (9) = 1;

    // std::cout << min (v) << std::endl;
    // std::cout << max (v) << std::endl;

    // std::cout << sum (v) << std::endl;
    // std::cout << mean (v) << std::endl;
    // std::cout << mean_iterative (v) << std::endl;
    // std::cout << median (v) << std::endl;
    // std::cout << mode (v) << std::endl;
    // std::cout << variance (v) << std::endl;
    // std::cout << variance_iterative (v) << std::endl;
    // std::cout << norm_1 (v) << std::endl;
    // std::cout << norm_2 (v) << std::endl;
    // std::cout << norm_inf (v) << std::endl;
    // std::cout << index_norm_inf (v) << std::endl;

    vector<double> v1 (6);
    v1 (0) = 1; v1 (1) = 14.5;
    v1 (2) = 3.7; v1 (3) = 2.3;
    v1 (4) = 3.699; v1 (5) = 7;

    std::cout << histogram (v1, 5) << std::endl;

    vector<double> v2 (6);
    v2 (0) = 1.0; v2 (1) = 3.7;
    v2 (2) = 6.4; v2 (3) = 9.1;
    v2 (4) = 11.8; v2 (5) = 14.5;

    std::cout << histogram (v1, v2) << std::endl;

}

