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

#ifndef _BOOST_UBLAS_RANDOM_INITIALIZATION_
#define _BOOST_UBLAS_RANDOM_INITIALIZATION_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/container/set.hpp>

#include <ctime>

namespace boost { namespace numeric { namespace ublas {

    /*
    *   \brief RandomInitialization simply selects num_clusters data points at random as the initial
    *   set of centroids.
    */
    class RandomInitialization {
    public:
        RandomInitialization () {
            gen.seed(static_cast<unsigned int> (std::time (0)));
        }

        /*
        *   \brief Initializes a set of centroids from the data through uniform random selection.
        *
        *   \tparam MatrixType The type of data points (int, float, double etc)
        *   \param data Data on which clustering is to be performed.
        *   \param num_clusters The number of clusters to be evaluated.
        *   \param cluster_centroids Container to store the generated centroids.
        *
        *   \return void
        */
        template <class MatrixType>
        void Initialize (const MatrixType &data, const size_t num_clusters, matrix<double> &centroids) {
            /*
            Since we are initializing the centroids using random data points,
            what happens if same data point becomes a centroid twice?
            This will lead to empty clusters. Should be handled by EmptyClusterPolicy
            or here only? eg set<size_t> selected_indices;
            */
            for (size_t i = 0; i < num_clusters; ++ i) {
                boost::random::uniform_int_distribution<> dist (0, data.size1 () - 1);
                size_t index = dist (gen);
                row (centroids, i) = row (data, index);
            }
        }
    private:
        boost::random::mt19937 gen;
    };
}}}

#endif