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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/ublas/kmeans/random_initialization.hpp>
#include <boost/numeric/ublas/kmeans/naive_kmeans.hpp>

namespace boost { namespace numeric { namespace ublas {

    template </*class MetricType = EuclideanDistanceMetric,*/
              class InitialPartitionPolicy = RandomInitialization,
              template<class> class EvaluateStepType = NaiveKMeans>
    class KMeans {
    public:
        KMeans (const size_t max_iterations = 100,
                const size_t n_init = 10,
                /*const MetricType distance_metric = MetricType (),*/
                const InitialPartitionPolicy partition_policy = InitialPartitionPolicy ()) :
                max_iterations (max_iterations),
                n_init (n_init),
                /*distance_metric (distance_metric),*/
                partition_policy (partition_policy) {}

        template <class MatrixType>
        void Cluster (const MatrixType &data,
                      const size_t num_clusters,
                      vector<size_t> &cluster_assignments) {
            
            matrix<double> cluster_centroids (num_clusters, data.size2 ());
            Cluster (data, num_clusters, cluster_centroids, cluster_assignments);
        }

        template <class MatrixType>
        void Cluster (const MatrixType &data,
              const size_t num_clusters,
              matrix<double> &cluster_centroids) {
            
            matrix<double> current_centroids = zero_matrix<double> (num_clusters, data.size2 ());

            double min_inertia;

            EvaluateStepType<const MatrixType> evaluation_step_type = EvaluateStepType<const MatrixType> (data/*, distance_metric*/);

            size_t num_inits = 0;
            while (num_inits++ < n_init) {

                GetInitialCentroids (data, num_clusters, current_centroids);

                // std::cout << "got init centroids!" << std::endl;

                if (max_iterations == 0) {
                    cluster_centroids = current_centroids;
                    continue;
                }

                double inertia;

                size_t num_iterations = 0;
                bool convergence = false, iterations_finished = false;
                matrix<double> new_cluster_centroids (current_centroids.size1 (), current_centroids.size2 ());
                // std::cout << "evaluate step init!" << std::endl;
                do {
                    /*
                    We don't have to repeatedly copy new_cluster_centroids to current_centroids
                    because of this alternating step.
                    */
                    // std::cout << "iterate " << num_iterations << "!" << std::endl;
                    if (num_iterations % 2)
                        inertia = evaluation_step_type.Iterate (new_cluster_centroids, current_centroids);
                    else
                        inertia = evaluation_step_type.Iterate (current_centroids, new_cluster_centroids);
                    // std::cout << "iterate " << num_iterations << "done!" << std::endl;

                    double norm = 0;
                    for (size_t i = 0; i < current_centroids.size2 (); ++ i)
                        norm += std::pow (inner_prod (column (current_centroids, i) - column (new_cluster_centroids, i), column (current_centroids, i) - column (new_cluster_centroids, i)), 2.0);

                    ++ num_iterations;

                    if (num_iterations == max_iterations)
                        iterations_finished = true;
                    if (norm <= 1e-6)
                        convergence = true;

                } while (!convergence && !iterations_finished);

                // std::cout << "evaluate iterations done!" << std::endl;

                // std::cout << new_cluster_centroids.size1() << " " << new_cluster_centroids.size2 () << std::endl;
                // std::cout << current_centroids.size1() << " " << current_centroids.size2 () << std::endl;
                
                if (num_inits == 1 || inertia < min_inertia) {
                    min_inertia = inertia;
                    cluster_centroids = new_cluster_centroids;
                }

                // for (size_t i = 0; i < cluster_centroids.size1 (); ++ i)
                //     for (size_t j = 0; j < cluster_centroids.size2 (); ++ j)
                //         cluster_centroids (i, j) = new_cluster_centroids (i, j);

                // std::cout << "new centroids assigned!" << std::endl;
            }

            // cluster_centroids /= n_init;
        }

        template <class MatrixType>
        void Cluster (const MatrixType &data,
              const size_t num_clusters,
              matrix<double> &cluster_centroids,
              vector<size_t> &cluster_assignments) {

            Cluster (data, num_clusters, cluster_centroids);

            // std::cout << "1st clusterd!" << std::endl;

            for (size_t i = 0; i < data.size1 (); ++ i) {
                matrix_row<const MatrixType> data_row (data, i);
                size_t assigned_cluster = 0;
                /*
                Use distance_metric.Apply () here, when implemented.
                double min_centroid_distance = distance_metric.Apply (data_row, row (cluster_centroids, 0));
                */
                double min_centroid_distance = inner_prod (data_row - row (cluster_centroids, 0), data_row - row (cluster_centroids, 0));
                for (size_t j = 1; j < cluster_centroids.size1 (); ++ j) {
                    /*
                    Use distance_metric.Apply () here, when implemented.
                    double min_centroid_distance = distance_metric.Apply (data_row, row (cluster_centroids, 0));
                    */
                    double distance = inner_prod (data_row - row (cluster_centroids, j), data_row - row (cluster_centroids, j));
                    if (distance < min_centroid_distance) {
                        min_centroid_distance = distance;
                        assigned_cluster = j;
                    }
                }
                cluster_assignments (i) = assigned_cluster;
            }
        }

        // size_t MaxIterations () const {
        //     return max_iterations;
        // }

        // size_t& MaxIterations () const {
        //     return max_iterations;
        // }

        /*MetricType DistanceMetric () {
            return distance_metric;
        }

        MetricType& DistanceMetric () {
            return distance_metric;
        }*/

        // InitialPartitionPolicy PartitionPolicy () {
        //     return partition_policy;
        // }

        // InitialPartitionPolicy& PartitionPolicy () {
        //     return partition_policy;
        // }

    private:
        /*
        Modify this method to support partiotion policies which give
        initial cluster assignments instead of cluster centroids.
        */
        template <class MatrixType>
        void GetInitialCentroids (const MatrixType &data, const size_t num_clusters, matrix<double> &centroids) {
            partition_policy.Initialize (data, num_clusters, centroids);
            // std::cout << centroids << std::endl;
        }

    private:
        size_t max_iterations;
        size_t n_init;
        /*MetricType distance_metric;*/
        InitialPartitionPolicy partition_policy;
    };
}}}

#endif