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

    /*
    *   \brief Implements the KMeans class. Macro parameters like the distance metric
    *   to be used, the initial partition method etc are provided as template parameters.
    *
    *   \tparam MetricType The type of distance metric to be used for evaluating
    *   node-cluster distances. Default = EuclideanDistanceMetric
    *
    *   \tparam InitialPartitionPolicy The policy/method to use to initialise centroids.
    *   Default = RandomInitialization.
    *
    *   \tparam EvaluateStepType The algorithm to use for a single iteration of evaluating
    *   new centroids. Default = NaiveKMeans.
    */
    template </*class MetricType = EuclideanDistanceMetric,*/
              class InitialPartitionPolicy = RandomInitialization,
              template<class> class EvaluateStepType = NaiveKMeans>
    class KMeans {
    public:
        KMeans (const size_t max_iterations = 100,
                const size_t n_init = 10) :
                max_iterations (max_iterations),
                n_init (n_init) {}

        /*
        *   \brief Performs clustering on data specified. Returns only the cluster
        *   assignments.
        *
        *   \tparam MatrixType The type of data points (int, float, double etc)
        *   \param data Data on which clustering is to be performed.
        *   \param num_clusters The number of clusters to be evaluated.
        *   \param cluster_assignments Container to store the evaluated assignments.
        *
        *   \return void
        */
        template <class MatrixType>
        void Cluster (const MatrixType &data,
                      const size_t num_clusters,
                      vector<size_t> &cluster_assignments) {
            
            matrix<double> cluster_centroids (num_clusters, data.size2 ());
            Cluster (data, num_clusters, cluster_centroids, cluster_assignments);
        }

        /*
        *   \brief Performs clustering on data specified. Returns only the cluster
        *   centroids.
        *
        *   \tparam MatrixType The type of data points (int, float, double etc)
        *   \param data Data on which clustering is to be performed.
        *   \param num_clusters The number of clusters to be evaluated.
        *   \param cluster_centroids Container to store the evaluated centroids.
        *
        *   \return void
        */
        template <class MatrixType>
        void Cluster (const MatrixType &data,
              const size_t num_clusters,
              matrix<double> &cluster_centroids) {
            
            matrix<double> current_centroids = zero_matrix<double> (num_clusters, data.size2 ());

            double min_inertia;

            InitialPartitionPolicy partition_policy = InitialPartitionPolicy ();
            EvaluateStepType<const MatrixType> evaluation_step_type = EvaluateStepType<const MatrixType> (data/*, distance_metric*/);

            // We will run kmeans on the data n_init times a different set of initial centroids
            // each time. The best set of centroids/assignments will be those which have the
            // least inertia.
            size_t num_inits = 0;
            while (num_inits++ < n_init) {

                GetInitialCentroids (data, partition_policy, num_clusters, current_centroids);

                if (max_iterations == 0) {
                    cluster_centroids = current_centroids;
                    continue;
                }

                double inertia;

                size_t num_iterations = 0;
                bool convergence = false, iterations_finished = false;
                matrix<double> new_cluster_centroids (current_centroids.size1 (), current_centroids.size2 ());
                do {
                    /*
                    We don't have to repeatedly copy new_cluster_centroids to current_centroids
                    because of this alternating step.
                    */
                    if (num_iterations % 2)
                        inertia = evaluation_step_type.Iterate (new_cluster_centroids, current_centroids);
                    else
                        inertia = evaluation_step_type.Iterate (current_centroids, new_cluster_centroids);

                    double norm = 0;
                    for (size_t i = 0; i < current_centroids.size2 (); ++ i)
                        norm += std::pow (inner_prod (column (current_centroids, i) - column (new_cluster_centroids, i), column (current_centroids, i) - column (new_cluster_centroids, i)), 2.0);

                    ++ num_iterations;

                    if (num_iterations == max_iterations)
                        iterations_finished = true;
                    if (norm <= 1e-6)
                        convergence = true;

                } while (!convergence && !iterations_finished);

                if (num_inits == 1 || inertia < min_inertia) {
                    min_inertia = inertia;
                    cluster_centroids = new_cluster_centroids;
                }
            }
        }

        /*
        *   \brief Performs clustering on data specified. Returns the cluster
        *   assignments and cluster centroids both. Calls the second Cluster method
        *   and evaluates the assignments by finding the closest cluster centroid.
        *
        *   \tparam MatrixType The type of data points (int, float, double etc)
        *   \param data Data on which clustering is to be performed.
        *   \param num_clusters The number of clusters to be evaluated.
        *   \param cluster_centroids Container to store the evaluated centroids.
        *   \param cluster_assignments Container to store the evaluated assignments.
        *
        *   \return void
        */
        template <class MatrixType>
        void Cluster (const MatrixType &data,
              const size_t num_clusters,
              matrix<double> &cluster_centroids,
              vector<size_t> &cluster_assignments) {

            Cluster (data, num_clusters, cluster_centroids);

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

    private:
        /*
        *   \brief Uses the partition policy specified to generate an initial set of
        *   centroids from the data. Modify this method to support partition policies
        *   which give initial cluster assignments instead of cluster centroids.
        *
        *   \param data Data on which clustering is to be performed.
        *   \param partition_policy The paritioning method to be used to generate
        *   centroids.
        *   \param num_clusters The number of centroids to be generated.
        *   \param centroids Container to store the generated centroids.
        *
        *   \return void
        */
        template <class MatrixType>
        void GetInitialCentroids (const MatrixType &data, InitialPartitionPolicy &partition_policy, const size_t num_clusters, matrix<double> &centroids) {
            partition_policy.Initialize (data, num_clusters, centroids);
        }

    private:
        size_t max_iterations;
        size_t n_init;
    };
}}}

#endif