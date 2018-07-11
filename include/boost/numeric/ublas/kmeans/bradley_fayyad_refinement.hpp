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

    class RefinedStart {
    public:
        RefinedStart (size_t samplings = 100,
                      double subsample = 0.05) :
                      num_samplings(samplings),
                      subsample_factor (subsample) {

            gen.seed(static_cast<unsigned int>(std::time(0)));
        }

        template <class MatrixType>
        void Initialize (const MatrixType &data, const size_t num_clusters, matrix<double> &centroids) {

            uniform_dist = boost::random::uniform_int_distribution<> (0, data.size1 () - 1);

            size_t num_sampled_points = subsample_factor * data.size1 ();
            MatrixType sampled_data (num_sampled_points, data.size2 ());
            bool is_point_sampled [data.size1 ()] = {0};

            matrix<double> sampled_centroids (num_samplings * num_clusters, data.size2 ());

            for (size_t i = 0; i < num_samplings; ++ i) {
                size_t cur_sampled = 0;
                while (cur_sampled < num_sampled_points) {
                    size_t sampled_index = uniform_dist (gen);
                    if (!is_point_sampled [sampled_index]) {
                        row (sampled_data, cur_sampled) = row (data, sampled_index);
                        ++ cur_sampled;
                    }
                }
                KMeans<> kmeans;
                kmeans.Cluster (sampled_data, num_clusters, centroids);
                subrange (sampled_centroids, i * num_clusters, (i + 1) * num_clusters, 0, data.size2 ()) = centroids;
            }
            
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