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

    class RandomInitialization {
    public:
        RandomInitialization () {}

        template <class MatrixType>
        void Initialize (const MatrixType &data, const size_t num_clusters, matrix<double> &centroids) {
            /*
            Since we are initializing the centroids using random data points,
            what happens if same data point becomes a centroid twice?
            This will lead to empty clusters. Should be handled by EmptyClusterPolicy
            or here only? eg set<size_t> selected_indices;
            */
            boost::random::mt19937 gen;
            gen.seed(static_cast<unsigned int> (std::time (0)));;
            
            for (size_t i = 0; i < num_clusters; ++ i) {
                boost::random::uniform_int_distribution<> dist (0, data.size1 () - 1);
                size_t index = dist (gen);
                row (centroids, i) = row (data, index);
                // for (size_t j = 0; j < centroids.size2 (); ++j)
                //     centroids (i, j) = data (index, j);
            }
        }
    };
}}}

#endif