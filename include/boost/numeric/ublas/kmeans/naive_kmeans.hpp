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

#ifndef _BOOST_UBLAS_NAIVE_KMEANS_
#define _BOOST_UBLAS_NAIVE_KMEANS_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/math/special_functions/pow.hpp>

#include <cmath>

namespace boost { namespace numeric { namespace ublas {

    /*
    *   \brief Implements the Naive kmeans algorith step. Macro parameters like the distance metric
    *   to be used etc are provided as template parameters.
    *   \tparam MetricType The type of distance metric to be used for evaluating
    *   node-cluster distances. Default = EuclideanDistanceMetric
    *   \tparam MatrixType The type of data points (int, float, double etc).
    *   \param data_set Data to be clustered.
    */
    template </* class MetricType = EuclideanDistanceMetric,*/
              class MatrixType = matrix<double> >
    class NaiveKMeans {
    public:
        NaiveKMeans (MatrixType &data_set) : 
                    data (data_set) {}

        /*
        *   \brief Performs a single iteration of naive kmeans algorithm. New centroids are evaluated
        *   from the cluster assignment of data points according to the old centroids. A new
        *   centroid is simply the mean of all data points belonging to the old centroid.
        *   \param centroids Old cluster centroids.
        *   \param new_centroids Container to store and return the new set of cluster centroids.
        *   \return The inertia of the current set of centroids, i.e the sum of squared distance of
        *   each node from its closest centroid.
        */
        double Iterate (const matrix<double> &centroids, matrix<double> &new_centroids) {
            vector<size_t> data_points_per_centroid  = zero_vector<size_t>(centroids.size1 ());
            new_centroids *= 0;

            vector<size_t> cluster_assignments (data.size1 ());

            double inertia = 0;

            // First we find the closest centroid for each data point and set the cluster
            // assignment accordingly.
            for (size_t i = 0; i < data.size1 (); ++ i) {
                matrix_row<MatrixType> data_row (data, i);
                size_t assigned_cluster = 0;
                
                double min_centroid_distance = inner_prod (data_row - row (centroids, 0), data_row - row (centroids, 0));
                
                for (size_t j = 0; j < centroids.size1 (); ++ j) {
                    
                    double distance = inner_prod (data_row - row (centroids, j), data_row - row (centroids, j));
                    
                    if (distance < min_centroid_distance) {
                        min_centroid_distance = distance;
                        assigned_cluster = j;
                    }
                }

                inertia += min_centroid_distance;

                cluster_assignments (i) = assigned_cluster;
                
                row (new_centroids, assigned_cluster) += data_row;
                
                data_points_per_centroid (assigned_cluster) += 1;
            }

            // Then, the new centroids are evaluated by taking the mean of data points
            // belonging to each centroid.
            for (size_t i = 0; i < new_centroids.size1 (); ++ i)
                if (data_points_per_centroid (i))
                    row (new_centroids, i) /= data_points_per_centroid (i);
            
            return inertia;
        }

    private:
        MatrixType data;
    };
}}}

#endif