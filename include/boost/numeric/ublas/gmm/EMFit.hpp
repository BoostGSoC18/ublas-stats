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

            matrix<double> component_probabilities (data.size (), distributions.size ());

            for (size_t i = 0; i < data.size (); ++ i) {
                for (size_t k = 0; k < distributions.size (); ++ k)
                    component_probabilities (i, k) = pdf (distributions[k], data (i));
                row (component_probabilities, i) /= sum (row (component_probabilities, i));
            }

            weights = sum (component_probabilities, 0) / data.size ();
            vector<double> new_means = prod (data, component_probabilities);
            for (size_t k = 0; k < distributions.size (); ++ k) {
                new_means (k) /= sum (column (component_probabilities, k));
                double new_variance = 0;;
                for (size_t i = 0; i < data.size (); ++ i)
                    new_variance += component_probabilities (i, k) * std::pow (data (i) - new_means (k), 2);
                new_variance /= sum (column (component_probabilities, k));
                distributions[k] = boost::math::normal (new_means (k), std::sqrt (new_variance));
            }
        }

    private:
        const VectorType data;
    };

}}}

#endif