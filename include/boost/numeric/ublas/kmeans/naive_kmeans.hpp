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

    template </* class MetricType,*/
              class MatrixType = matrix<double> >
    class NaiveKMeans {
    public:
        NaiveKMeans (MatrixType &data_set/*, MetricType distance_metric*/) : 
                    /*distance_metric (distance_metric),*/
                    data (data_set) {

            // for (size_t i = 0; i < data_set.size1 (); ++ i)
            //     for (size_t j = 0; j < data_set.size2 (); ++ i)
            //         data (i, j) = data_set (i, j);
        }

        double Iterate (matrix<double> &centroids, matrix<double>& new_centroids) {
            vector<int> data_points_per_centroid  = zero_vector<int>(centroids.size1 ());
            new_centroids *= 0;

            vector<int> cluster_assignments (data.size1 ());

            // std::cout << "init centroids!" << centroids << std::endl;
            // std::cout << "new centroids!" << new_centroids << std::endl;

            for (size_t i = 0; i < data.size1 (); ++ i) {
                matrix_row<MatrixType> data_row (data, i);
                size_t assigned_cluster = 0;
                /*
                Use distance_metric.Apply () here, when implemented.
                double min_centroid_distance = distance_metric.Apply (data_row, row (centroids, 0));
                */
                double min_centroid_distance = inner_prod (data_row - row (centroids, 0), data_row - row (centroids, 0));
                for (size_t j = 0; j < centroids.size1 (); ++ j) {
                    /*
                    Use distance_metric.Apply () here, when implemented.
                    double min_centroid_distance = distance_metric.Apply (data_row, row (centroids, 0));
                    */
                    double distance = inner_prod (data_row - row (centroids, j), data_row - row (centroids, j));
                    // std::cout << "distance from " << row (centroids, j) << " " << j << ": " << distance << std::endl;
                    if (distance < min_centroid_distance) {
                        min_centroid_distance = distance;
                        assigned_cluster = j;
                    }
                }
                // std::cout << "assigned cluster to " << i << ", " << data_row << ": " << assigned_cluster << std::endl;
                cluster_assignments (i) = assigned_cluster;
                
                row (new_centroids, assigned_cluster) += data_row;
                // for (size_t i = 0; i < new_centroids.size2 (); ++i)
                //     new_centroids (assigned_cluster, i) += data_row (i);

                
                data_points_per_centroid (assigned_cluster) += 1;
            }

            // std::cout << "data points per centroid: " << data_points_per_centroid << std::endl;

            for (size_t i = 0; i < new_centroids.size1 (); ++ i)
                if (data_points_per_centroid (i))
                    row (new_centroids, i) /= data_points_per_centroid (i);
                    // for (size_t j = 0; j < new_centroids.size2 (); ++j)
                    //     new_centroids (i, j) /= data_points_per_centroid (i);

            /*
            Use distance_metric.Apply () here, when implemented.
            return distance_metric.Apply (centroids, new_centroids);
            */
            double distance_norm = 0;
            for (size_t i = 0; i < centroids.size2 (); ++ i) {
                distance_norm += std::pow (inner_prod (column (centroids, i) - column (new_centroids, i), column (centroids, i) - column (new_centroids, i)), 2.0);
            }
            return std::pow (distance_norm, 0.5);

        }

    private:
        MatrixType data;
        /*MetricType distance_metric;*/
    };
}}}

#endif