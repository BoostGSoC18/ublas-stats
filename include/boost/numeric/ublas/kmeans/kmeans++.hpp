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

    /*
    *   \brief KMeans++ initialization provides a more optimal set of initial centroids than random
    *   initialization. Each new centroid is selected from the data using a weighted probability
    *   distribution of the distance of each point from the closest existing centroid.
    */
    class KMeansPlusPlus {
    public:
        KMeansPlusPlus () {
            gen.seed(static_cast<unsigned int>(std::time(0)));
        }

        /*
        *   \brief Initializes a set of centroids from the data by first generating selecting a
        *   random data point as a centroid, and then selecting the one furthest from the current
        *   set of centroids as the next centroid.
        *   \tparam MatrixType The type of data points (int, float, double etc)
        *   \param data Data on which clustering is to be performed.
        *   \param num_clusters The number of clusters to be evaluated.
        *   \param cluster_centroids Container to store the generated centroids.
        */
        template <class MatrixType>
        void Initialize (const MatrixType &data, const size_t num_clusters, matrix<double> &centroids) {
            boost::random::uniform_int_distribution<> uniform_dist (0, data.size1 () - 1);
            
            // Select the first centroid.
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
                
                // To choose the next centroid, we first calculate the weight of each
                // data point by finding the distance from closest centroid.
                for (size_t data_counter = 0; data_counter < data.size1 (); ++ data_counter)
                    if (indices_selected_as_centroid.find (data_counter) == indices_selected_as_centroid.end ())
                        weights (data_counter) = 0;
                    else
                        weights (data_counter) = closest_centroid_distance (data_counter);

                // Now, we generate the next centroid from this weighted probability
                // distribution.
                boost::random::discrete_distribution<> discrete_dist (weights);
                size_t new_centroid_index = discrete_dist (gen);
                row (centroids, i) = row (data, new_centroid_index);
                indices_selected_as_centroid.insert (new_centroid_index);

                // We update the closest centroid distance for each data point by checking
                // its distance from the latest centroid.
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