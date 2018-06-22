//  Dattatreya Mohapatra
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

#ifndef _BOOST_UBLAS_REFINED_START_
#define _BOOST_UBLAS_REFINED_START_


#include <boost/numeric/ublas/matrix.hpp>

namespace boost { namespace numeric { namespace ublas {

    class RefinedStart {
    public:
        RefinedStart () {}

        template <class MatrixType>
        static Initialize (const MatrixType &data, const size_type num_clusters, const matrix<double> &centroids) {}
    }
}

#endif