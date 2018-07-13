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

#ifndef _BOOST_UBLAS_REFINED_START_
#define _BOOST_UBLAS_REFINED_START_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <ctime>

namespace boost { namespace numeric { namespace ublas {

    /*
    *   \brief The initialization method proposed by Bradley and Fayyad provides a more refined set
    *   of initial centroids by clustering over multiple subsamples.
    */
    class RefinedStart {
    public:

        /*
        *   \param samplings: The number of times clustering is to be performed over sub-sampled data.
        *   \param subsample: The sub-sampling factor for the dataset.
        */
        RefinedStart (size_t samplings = 100,
                      double subsample = 0.05) :
                      num_samplings(samplings),
                      subsample_factor (subsample) {

            gen.seed(static_cast<unsigned int>(std::time(0)));
        }

        /*
        *   \brief Initializes a set of centroids from the data by first generating centroids from
        *   subsampled data, and then performing kmeans on those centroids as data points.
        *
        *   \tparam MatrixType The type of data points (int, float, double etc)
        *   \param data Data on which clustering is to be performed.
        *   \param num_clusters The number of clusters to be evaluated.
        *   \param cluster_centroids Container to store the generated centroids.
        *
        *   \return void
        */
        template <class MatrixType>
        void Initialize (const MatrixType &data, const size_t num_clusters, matrix<double> &centroids) {

            uniform_dist = boost::random::uniform_int_distribution<> (0, data.size1 () - 1);

            // Calculate the number of data points to be sampled due to sub-sampling.
            size_t num_sampled_points = subsample_factor * data.size1 ();
            MatrixType sampled_data (num_sampled_points, data.size2 ());
            bool is_point_sampled [data.size1 ()] = {0};

            matrix<double> sampled_centroids (num_samplings * num_clusters, data.size2 ());

            for (size_t i = 0; i < num_samplings; ++ i) {
                // Sample data points.
                size_t cur_sampled = 0;
                while (cur_sampled < num_sampled_points) {
                    size_t sampled_index = uniform_dist (gen);
                    if (!is_point_sampled [sampled_index]) {
                        row (sampled_data, cur_sampled) = row (data, sampled_index);
                        ++ cur_sampled;
                    }
                }
                // Perform KMeans on the sampled data. The centroids obtained here are
                // stored in a separate array for another round of selection in the future.
                KMeans<> kmeans;
                kmeans.Cluster (sampled_data, num_clusters, centroids);
                subrange (sampled_centroids, i * num_clusters, (i + 1) * num_clusters, 0, data.size2 ()) = centroids;
            }
            
            // Now, from the array of set sampled centroids obtained for each subsampling process,
            // extract the final set of centroids for the entire data by performing KMeans on them.
            KMeans<> kmeans;
            kmeans.Cluster (sampled_centroids, num_clusters, centroids);
        }

    private:
        size_t num_samplings;
        double subsample_factor;

        boost::random::mt19937 gen;
        boost::random::uniform_int_distribution<> uniform_dist;
    };
}}}

#endif