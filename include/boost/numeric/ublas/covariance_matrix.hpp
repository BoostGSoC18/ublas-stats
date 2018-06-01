// Dattatreya Mohapatra
//
//
//  Copyright (c) 2000-2010
//  Joerg Walter, Mathias Koch, Gunter Winkler, David Bellot
//  Copyright (c) 2014, Athanasios Iliopoulos
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.

#ifndef _BOOST_UBLAS_COVARIANCE_MATRIX_
#define _BOOST_UBLAS_COVARIANCE_MATRIX_
#endif

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace boost { namespace numeric { namespace ublas {

    template <class M>
    symmetric_matrix<typename M::value_type, lower> covariance_matrix (const M &m, bool rowvar = true) {
        typedef typename M::value_type value_type;
        typedef typename M::size_type size_type;

        size_type nvars = 0;
        size_type nobvs = 0;

        if (rowvar) {
            nvars = m.size1 ();
            nobvs = m.size2 ();
        }
        else {
            nvars = m.size2 ();
            nobvs = m.size1 ();
        }

        symmetric_matrix<value_type, lower> cov_matrix (nvars);

        for (size_type i = 0; i < nvars; ++ i)
            for (size_type j = 0; j <= i; ++ j)
                if (rowvar)
                    cov_matrix (i, j) = covariance (row (m, i), row (m, j));
                else
                    cov_matrix (i, j) = covariance (column (m, i), column (m, j));

        return cov_matrix;
    }

}}}