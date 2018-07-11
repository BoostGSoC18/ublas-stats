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

#ifndef _BOOST_UBLAS_KMEANSPLUSPLUS_
#define _BOOST_UBLAS_KMEANSPLUSPLUS_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/container/set.hpp>

#include <ctime>

namespace boost { namespace numeric { namespace ublas {

    class KMeansPlusPlus {
    public:
        KMeansPlusPlus () {
            gen.seed(static_cast<unsigned int>(std::time(0)));
        }

        template <class MatrixType>
        void Initialize (const MatrixType &data, const size_t num_clusters, matrix<double> &centroids) {
            boost::random::uniform_int_distribution<> uniform_dist (0, data.size1 () - 1);
            
            size_t first_centroid_index = uniform_dist (gen);
            row (centroids, 0) = row (data, first_centroid_index);

            boost::container::set<size_t> indices_selected_as_centroid;
            indices_selected_as_centroid.insert (first_centroid_index);
            boost::container::set<size_t>::iterator set_it = indices_selected_as_centroid.begin ();

            vector<double> closest_centroid_distance (data.size1 ());
            for (size_t i = 0; i < data.size1 (); ++ i)
                closest_centroid_distance (i) = inner_prod (row (data, i) - row (data, *set_it), row (data, i) - row (data, *set_it));

            vector<double> weights (data.size1 ());
            for (size_t i = 1; i < num_clusters; ++ i) {
                for (size_t data_counter = 0; data_counter < data.size1 (); ++ data_counter)
                    if (indices_selected_as_centroid.find (data_counter) == indices_selected_as_centroid.end ())
                        weights (data_counter) = 0;
                    else
                        weights (data_counter) = closest_centroid_distance (data_counter);

                boost::random::discrete_distribution<> discrete_dist (weights);
                size_t new_centroid_index = discrete_dist (gen);
                row (centroids, i) = row (data, new_centroid_index);
                indices_selected_as_centroid.insert (new_centroid_index);

                for (size_t data_counter = 0; data_counter < data.size1 (); ++ data_counter) {
                    double new_centroid_distance = inner_prod (row (data, data_counter) - row (data, new_centroid_index), row (data, data_counter) - row (data, new_centroid_index));
                    if (new_centroid_distance < closest_centroid_distance (data_counter))
                        closest_centroid_distance (data_counter) = new_centroid_distance;
                }
            }            
        }
    private:
        boost::random::mt19937 gen;
    };
}}}

#endif