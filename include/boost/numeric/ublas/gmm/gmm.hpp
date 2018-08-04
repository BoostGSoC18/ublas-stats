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

#ifndef _BOOST_UBLAS_GMM_
#define _BOOST_UBLAS_GMM_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "EMFit.hpp"

#include <boost/math/distributions/normal.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

#include <ctime>
#include <cmath>
#include <vector>

namespace boost { namespace numeric { namespace ublas {

    class GMM {
    public:
        GMM (const size_t n_components) :
            n_components (n_components),
            distributions (n_components, boost::math::normal ()) {

            random_generator.seed (static_cast<unsigned int> (std::time (0)));

            weights = scalar_vector<double> (n_components, 1.0 / n_components);
        }

        GMM (const size_t n_components,
             vector<double> &weights,
             vector<double> &means,
             vector<double> &stdevs) :
            n_components (n_components),
            weights (weights) {


            random_generator.seed (static_cast<unsigned int> (std::time (0)));

            for (size_t i = 0; i < n_components; ++ i)
                distributions.push_back (boost::math::normal (means (i), stdevs (i)));
        }

        double Sample () {
            double gen_sample = 0;
            for (size_t i = 0; i < distributions.size (); ++ i) {
                boost::random::normal_distribution<> normal_distribution (distributions[i].mean (), distributions[i].standard_deviation ());
                gen_sample += weights (i) * normal_distribution (random_generator);
            }
        }

        template <class ValueType>
        double GetProbability (const ValueType x) {
            double weighted_pdf_sum = 0;
            for (size_t i = 0; i < distributions.size (); ++ i)
                weighted_pdf_sum += weights (i) * pdf (distributions[i], x);
            return weighted_pdf_sum;
        }

        template <class ValueType>
        double GetProbability (const ValueType x, size_t component_label) {
            return pdf (distributions[component_label], x);
        }

        template <class ValueType>
        size_t GetComponentLabel (const ValueType x) {
            size_t closest_component = 0;
            double highest_pdf = pdf (distributions[0], x);
            for (size_t i = 1; i < distributions.size (); ++ i) {
                double current_pdf = pdf (distributions[i], x);
                if (current_pdf > highest_pdf ) {
                    highest_pdf = current_pdf;
                    closest_component = i;
                }
            }
            return closest_component;
        }

        template <class VectorType,
                  class FittingPolicy = EMFit<VectorType> >
        void Train (const VectorType &data,
                    const size_t max_iterations = 1000,
                    const size_t n_trials = 10,
                    const double TOL = 1e-9) {

            size_t i_trials = 0;
            
            do {

                FittingPolicy fitter = FittingPolicy (data);

                size_t num_iterations = 0;
                bool convergence = false, iterations_finished = false;

                double current_loglikelihood = LogLikelihood (data, weights, distributions);

                InitializeComponents (data);

                do {
                    fitter.Step (distributions, weights);
                    
                    double new_loglikelihood = LogLikelihood (data, weights, distributions);

                    ++ num_iterations;
                    
                    if (num_iterations == max_iterations)
                        iterations_finished = true;
                    if (std::abs (new_loglikelihood - current_loglikelihood) < TOL)
                        convergence = true;

                } while (!convergence && !iterations_finished);

                // std::cout << distributions[0].mean () << " " << distributions[0].standard_deviation () << std::endl;

                ++ i_trials;

            } while (i_trials < 10);
        }

        double GetComponentMean (size_t component_label) {
            return distributions[component_label].mean ();
        }

        double GetComponentStandardDeviation (size_t component_label) {
            return distributions[component_label].standard_deviation ();
        }

        double GetComponentWeight (size_t component_label) {
            return weights[component_label];
        }

    private:

        template <class VectorType>
        void InitializeComponents (const VectorType &data) {
            double stdev = std::sqrt (variance (data));
            boost::random::uniform_int_distribution<> uniform_dist (0, data.size () - 1);
            for (size_t i = 0; i < distributions.size (); ++ i) {
                size_t index = uniform_dist (random_generator);
                double mean = data (index);
                distributions[i] = boost::math::normal (mean, stdev);
                weights (i) = 1.0 / n_components;
            }
        }

        template <class VectorType>
        double LogLikelihood (const VectorType &data,
                              const vector<double> &weights,
                              const std::vector<boost::math::normal> &distributions) {

            double log_likelihood = 0;
            for (size_t i = 0; i < data.size (); ++ i) {
                double likelihood = 0;
                for (size_t j = 0; j < distributions.size (); ++ j)
                    likelihood += weights (j) * pdf (distributions[j], data (i));
                if (likelihood != 0)
                    log_likelihood += std::log (likelihood);
            }

            return log_likelihood;
        }

    private:
        size_t n_components;
        vector<double> weights;
        std::vector<boost::math::normal> distributions;
        boost::random::mt19937 random_generator;
    };

}}}

#endif