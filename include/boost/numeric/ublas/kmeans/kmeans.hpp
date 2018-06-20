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
        KMeans (const size_type max_iterations = 1000,
                /*MetricType distance_metric = MetricType (),*/
                InitialPartitionPolicy partition_policy = InitialPartitionPolicy (),
                EvaluateStepType evaluation_step = EvaluateStepType ()) {}

        void Cluster (const MatrixType &data,
                      const size_type num_clusters,
                      vector<int> &cluster_assignments) {}

        void Cluster (const MatrixType &data,
              const size_type num_clusters,
              matrix<double> &cluster_centroids) {}

        void Cluster (const MatrixType &data,
              const size_type num_clusters,
              matrix<double> &cluster_centroids,
              vector<int> &cluster_assignments) {}

        size_type MaxIterations () const {
            return max_iterations;
        }

        size_type& MaxIterations () const {
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

        EvaluateStepType EvaluationStep () {
            return evaluation_step;
        }

        EvaluateStepType& EvaluationStep () {
            return evaluation_step;
        }

    private:
        size_type max_iterations;
        /*MetricType distance_metric;*/
        InitialPartitionPolicy partition_policy;
        EvaluationStep evaluation_step;
    };
}

#endif