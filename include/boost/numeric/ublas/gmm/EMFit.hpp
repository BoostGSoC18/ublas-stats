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

#ifndef _BOOST_UBLAS_EMFIT_
#define _BOOST_UBLAS_EMFIT_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/math/distributions/normal.hpp>

namespace boost { namespace numeric { namespace ublas {

    template<class VectorType>
    class EMFit {
    public:
        EMFit (const VectorType &data) : data (data) {}

        void Step (std::vector<boost::math::normal> &distributions,
                   vector<double> &weights) {

        }

    private:
        VectorType data;
    };

}}}

#endif