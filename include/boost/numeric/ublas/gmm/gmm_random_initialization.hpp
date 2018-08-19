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

#ifndef _BOOST_UBLAS_GMM_RANDOM_INITIALIZATION_
#define _BOOST_UBLAS_GMM_RANDOM_INITIALIZATION_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/container/set.hpp>

#include <ctime>

namespace boost { namespace numeric { namespace ublas {

    /*
    *   \brief RandomInitialization simply selects n_components data points at random as the means for the initial
    *   set of Gaussian distributions.
    */
    class GMM_RandomInitialization {
    public:
        GMM_RandomInitialization () {
            gen.seed(static_cast<unsigned int> (std::time (0)));
        }

        /*
        *   \brief Initializes a set of Gaussian distributions from the data through uniform random selection.
        *
        *   \tparam VectorType The type of data points (int, float, double etc)
        *   \param data Data on which GMM is to be trained.
        *   \param n_components The number of Gaussian components to be evaluated.
        *
        *   \return void
        */
        template <class VectorType>
        void Initialize (const VectorType &data, const size_t n_components, std::vector<boost::math::normal> &distributions) {
            double stdev = std::sqrt (variance (data));
            boost::random::uniform_int_distribution<> uniform_dist (0, data.size () - 1);
            for (size_t i = 0; i < distributions.size (); ++ i) {
                size_t index = uniform_dist (gen);
                double mean = data (index);
                distributions[i] = boost::math::normal (mean, stdev);
            }
        }
    private:
        boost::random::mt19937 gen;
    };
}}}

#endif