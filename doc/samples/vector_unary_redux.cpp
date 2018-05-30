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

int main () {
    using namespace boost::numeric::ublas;
    vector<int> v (10);
    for (unsigned i = 0; i < v.size (); ++ i)
        v (i) = i;
    v (0) = 1;
    v (9) = 1;

    std::complex<double> s(1,2);
    std::complex<double> si(4,5);
    std::cout << type_traits< std::complex<double> >::type_abs (s) << std::endl;
    std::cout << type_traits< std::complex<int> >::type_abs (si) << std::endl;
    std::cout << type_traits< double >::type_abs (v (2)) << std::endl;

    // std::cout << min (v) << std::endl;
    // std::cout << max (v) << std::endl;

    // std::cout << sum (v) << std::endl;
    std::cout << mean (v) << std::endl;
    // std::cout << mean_iterative (v) << std::endl;
    // std::cout << median (v) << std::endl;
    std::cout << mode (v) << std::endl;
    std::cout << variance (v) << std::endl;
    std::cout << variance_iterative (v) << std::endl;
    std::cout << norm_1 (v) << std::endl;
    std::cout << norm_2 (v) << std::endl;
    std::cout << norm_inf (v) << std::endl;
    std::cout << index_norm_inf (v) << std::endl;
}

