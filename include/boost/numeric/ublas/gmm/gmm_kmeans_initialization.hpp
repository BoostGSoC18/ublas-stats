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

#ifndef _BOOST_UBLAS_KMEANS_INITIALIZATION_
#define _BOOST_UBLAS_KMEANS_INITIALIZATION_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/container/set.hpp>

#include <boost/numeric/ublas/kmeans/kmeans.hpp>
#include <boost/numeric/ublas/kmeans/kmeans++.hpp>
#include <boost/numeric/ublas/kmeans/random_initialization.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <ctime>

namespace boost { namespace numeric { namespace ublas {

    /*
    *   \brief KMeansInitialization simply selects n_components data points by performing KMeans clustering on the data.
    */
    class GMM_KMeansInitialization {
    public:
        GMM_KMeansInitialization () {}

        /*
        *   \brief Initializes a set of Gaussian distributions from the data using KMeans.
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
            KMeans<RandomInitialization> kmeans;
            
            matrix<typename VectorType::value_type> data_matrix (data.size (), 1);
            for (size_t i = 0; i < data.size (); ++ i)
                data_matrix (i, 0) = data (i);

            matrix<double> means (n_components, 1);
            kmeans.Cluster (data_matrix, n_components, means);
            std::cout << data_matrix << std::endl;
            std::cout << means << std::endl;
            for (size_t i = 0; i < distributions.size (); ++ i) {
                distributions[i] = boost::math::normal (means (i, 0), stdev);
            }
        }
    };
}}}

#endif