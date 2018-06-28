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
        KMeans (const size_t max_iterations = 10,
                /*const MetricType distance_metric = MetricType (),*/
                const InitialPartitionPolicy partition_policy = InitialPartitionPolicy ()) :
                max_iterations (max_iterations),
                /*distance_metric (distance_metric),*/
                partition_policy (partition_policy) {}

        template <class MatrixType>
        void Cluster (const MatrixType &data,
                      const size_t num_clusters,
                      vector<int> &cluster_assignments) {
            
            matrix<double> cluster_centroids (num_clusters, data.size2 ());
            Cluster (data, num_clusters, cluster_centroids, cluster_assignments);
        }

        template <class MatrixType>
        void Cluster (const MatrixType &data,
              const size_t num_clusters,
              matrix<double> &cluster_centroids) {
            
            GetInitialCentroids (data, num_clusters, cluster_centroids);

            // std::cout << "got init centroids!" << std::endl;

            double norm = 0;
            size_t num_iterations = 0;
            bool convergence = false, iterations_finished = false;
            matrix<double> new_cluster_centroids (cluster_centroids.size1 (), cluster_centroids.size2 ());
            EvaluateStepType<const MatrixType> evaluation_step_type = EvaluateStepType<const MatrixType> (data/*, distance_metric*/);
            // std::cout << "evaluate step init!" << std::endl;
            do {
                /*
                We don't have to repeatedly copy new_cluster_centroids to cluster_centroids
                because of this alternating step.
                */
                // std::cout << "iterate " << num_iterations << "!" << std::endl;
                if (num_iterations % 2)
                    norm = evaluation_step_type.Iterate (new_cluster_centroids, cluster_centroids);
                else
                    norm = evaluation_step_type.Iterate (cluster_centroids, new_cluster_centroids);
                // std::cout << "iterate " << num_iterations << "done!" << std::endl;
                ++ num_iterations;

                if (num_iterations == max_iterations)
                    iterations_finished = true;
                if (norm <= 1e-6)
                    convergence = true;

            } while (!convergence && !iterations_finished);

            // std::cout << "evaluate iterations done!" << std::endl;

            // std::cout << new_cluster_centroids.size1() << " " << new_cluster_centroids.size2 () << std::endl;
            // std::cout << cluster_centroids.size1() << " " << cluster_centroids.size2 () << std::endl;
            
            cluster_centroids = new_cluster_centroids;
            // for (size_t i = 0; i < cluster_centroids.size1 (); ++ i)
            //     for (size_t j = 0; j < cluster_centroids.size2 (); ++ j)
            //         cluster_centroids (i, j) = new_cluster_centroids (i, j);

            // std::cout << "new centroids assigned!" << std::endl;
        }

        template <class MatrixType>
        void Cluster (const MatrixType &data,
              const size_t num_clusters,
              matrix<double> &cluster_centroids,
              vector<int> &cluster_assignments) {

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
        }

    private:
        size_t max_iterations;
        /*MetricType distance_metric;*/
        InitialPartitionPolicy partition_policy;
    };
}}}

#endif