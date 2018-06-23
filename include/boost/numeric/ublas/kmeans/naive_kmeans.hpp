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


#include <boost/numeric/ublas/matrix.hpp>

namespace boost { namespace numeric { namespace ublas {

    class NaiveKMeans {
    public:
        template <class MatrixType/*, class MetricType*/>
        NaiveKMeans (MatrixType &data/*, MetricType distance_metric*/) :
        /*distance_metric (distance_metric),*/
        data (data) {}

        double Iterate (matrix<double> &centroids, matrix<double>& new_centroids) {
            vector<int>  data_points_per_centroid (centroids.size1 ());

            vector<int> cluster_assignments (data.size1 ());

            for (size_t i = 0; i < data.size1 (); ++ i) {
                matrix_row<MatrixType> data_row (data, i);
                size_t assigned_cluster = 0;
                /*
                Use distance_metric.Apply () here, when implemented.
                double min_centroid_distance = distance_metric.Apply (data_row, row (cluster_centroids, 0));
                */
                double min_centroid_distance = inner_prod (data_row - row (cluster_centroids, 0));
                for (size_t j = 1; j < cluster_centroids.size1 (); ++ j) {
                    /*
                    Use distance_metric.Apply () here, when implemented.
                    double min_centroid_distance = distance_metric.Apply (data_row, row (cluster_centroids, 0));
                    */
                    double distance = inner_prod (data_row - row (cluster_centroids, j));
                    if (distance < min_centroid_distance) {
                        min_centroid_distance = distance;
                        assigned_cluster = j;
                    }
                }
                cluster_assignments (i) = assigned_cluster;
                new_centroids (assigned_cluster) += data_row;
                data_points_per_centroid (assigned_cluster) += 1;
            }

            for (size_t i = 0; i < new_centroids.size1 (); ++ i)
                new_centroids (i) /= data_points_per_centroid (i);

        }

    private:
        MatrixType data;
        /*MetricType distance_metric;*/
    }
}

#endif