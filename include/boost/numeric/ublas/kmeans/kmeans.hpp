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

#ifndef _BOOST_UBLAS_KMEANS_
#define _BOOST_UBLAS_KMEANS_

#include <boost/numeric/ublas/matrix.hpp>

namespace boost { namespace numeric { namespace ublas {

    template </*class MetricType = EuclideanDistanceMetric,*/
              class InitialPartitionPolicy = RandomInitialization,
              class EvaluateStepType = NaiveKMeans,
              class MatrixType = matrix<double>>
    class KMeans {
    public:
        KMeans (const size_t max_iterations = 1000,
                /*const MetricType distance_metric = MetricType (),*/
                const InitialPartitionPolicy partition_policy = InitialPartitionPolicy (),
                const EvaluateStepType evaluation_step = EvaluateStepType ()) :
                max_iterations (max_iterations),
                /*distance_metric (distance_metric),*/
                partition_policy (partition_policy) {}

        void Cluster (const MatrixType &data,
                      const size_t num_clusters,
                      vector<int> &cluster_assignments) {
            
            matrix<double> cluster_centroids (num_clusters, data.size2 ());
            Cluster (data, num_clusters, cluster_centroids, cluster_assignments);
        }

        void Cluster (const MatrixType &data,
              const size_t num_clusters,
              matrix<double> &cluster_centroids) {
            
            GetInitialCentroids (data, num_clusters, centroids);

            double norm = 0;
            size_t num_iterations = 0;
            bool convergence = false, iterations_finished = false;
            matrix<double> new_centroids (centroids.size1 (), centroids.size2 ());
            EvaluateStepType evaluation_step_type = EvaluateStepType (data/*, distance_metric*/);
            do {
                /*
                We don't have to repeatedly copy new_centroids to centroids
                because of this alternating step.
                */
                if (num_iterations % 2)
                    norm = evaluation_step.Iterate (new_centroids, centroids);
                else
                    norm = evaluation_step.Iterate (centroids, new_centroids);
                ++ num_iterations;

                if (num_iterations == max_iterations)
                    iterations_finished = true;
                if (norm <= 1e-6)
                    convergence = true;

            } while (!convergence && !iterations_finished);

            centroids = new_centroids;
        }

        void Cluster (const MatrixType &data,
              const size_t num_clusters,
              matrix<double> &cluster_centroids,
              vector<int> &cluster_assignments) {

            Cluster (data, num_clusters, cluster_centroids);

            vector<double> cluster_centroid_distance (num_clusters);
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
                cluster_assignments = assigned_cluster;
            }
        }

        size_t MaxIterations () const {
            return max_iterations;
        }

        size_t& MaxIterations () const {
            return max_iterations;
        }

        /*MetricType DistanceMetric () {
            return distance_metric;
        }

        MetricType& DistanceMetric () {
            return distance_metric;
        }*/

        InitialPartitionPolicy PartitionPolicy () {
            return partition_policy;
        }

        InitialPartitionPolicy& PartitionPolicy () {
            return partition_policy;
        }

    private:
        /*
        Modify this method to support partiotion policies which give
        initial cluster assignments instead of cluster centroids.
        */
        template <MatrixType = matrix<double> >
        void GetInitialCentroids (const MatrixType &data, const size_t num_clusters, const matrix<double> &centroids) {
            partition_policy.Initialize (data, num_clusters, centroids);
        }

    private:
        size_t max_iterations;
        /*MetricType distance_metric;*/
        InitialPartitionPolicy partition_policy;
    };
}

#endif