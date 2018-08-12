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

    /*
    *   \brief Implements the univariate GMM class.
    *   \param n_components Number of gaussian distributions in the mixture model.
    *   The gaussian components are initialized with zero mean and unit standard
    *   deviation. Each component has an equal weight.
    */
    class GMM {
    public:
        GMM (const size_t n_components) :
            n_components (n_components),
            distributions (n_components, boost::math::normal ()) {

            random_generator.seed (static_cast<unsigned int> (std::time (0)));

            weights = scalar_vector<double> (n_components, 1.0 / n_components);
        }

        /*
        *   \brief Overloaded constructor for the GMM class.
        *   \param n_components n_components Number of gaussian distributions in the mixture model.
        *   \param weights Pre-calculated weights for each gaussian component
        *   \param means Mean of each gaussian component.
        *   \param stdevs Standard deviation of each gaussian component.
        */
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

        /*
        *   \brief Generate a random sample from this GMM. A random sample from a GMM is
        *   actually a weighted sum of random samples from each of the component Gaussians.
        *   \returns A single random double precision value, since only univariate GMMs are
        *   supported currently.
        */
        double Sample () {
            double gen_sample = 0;
            for (size_t i = 0; i < distributions.size (); ++ i) {
                boost::random::normal_distribution<> normal_distribution (distributions[i].mean (), distributions[i].standard_deviation ());
                gen_sample += weights (i) * normal_distribution (random_generator);
            }
        }

        /*
        *   \brief Generate the probability of a sample belonging to this mixture model.
        *   \tparam ValueType The type of sample (int, float, double etc).
        *   \param x Sample whose probability is to be calculated.
        *   \returns Weighted sum of PDF of all component distributions for the given sample.
        */
        template <class ValueType>
        double GetProbability (const ValueType x) {
            double weighted_pdf_sum = 0;
            for (size_t i = 0; i < distributions.size (); ++ i)
                weighted_pdf_sum += weights (i) * pdf (distributions[i], x);
            return weighted_pdf_sum;
        }

        /*
        *   \brief Generate the probability of a sample belonging to the given component distribution.
        *   \tparam ValueType The type of sample (int, float, double etc).
        *   \param x Sample whose probability is to be calculated.
        *   \param component_label The 0-indexed label of the component distribution to be evaluated.
        *   \returns PDF of the given component distribution for the given sample.
        */
        template <class ValueType>
        double GetProbability (const ValueType x, size_t component_label) {
            return pdf (distributions[component_label], x);
        }

        /*
        *   \brief Classify the given sample as belonging to one of the component distributions.
        *   \tparam ValueType The type of sample (int, float, double etc).
        *   \param x Sample to be classified.
        *   \returns The 0-indexed label of the closest component distribution, by checking with
        *   the PDF of each component distribution.
        */
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

        /*
        *   \brief Train the GMM on the given data.
        *   \tparam VectorType The type of data points (int, float, double etc).
        *   \tparam FittingPolicy The algorithm used for a single step of evaluating new weights,
        *   means and standard deviations for the component distributions. Default = EMFit.
        *   \param data Data on which the mixture model is to be fitted.
        *   \param max_iterations Maximum number of training steps (FittingPolicy) steps to be performed.
        *   \param TOL The tolerance value for distance between two parameters for convergence.
        */
        template <class VectorType,
                  class FittingPolicy = EMFit<VectorType> >
        void Train (const VectorType &data,
                    const size_t max_iterations = 1000,
                    const double TOL = 1e-9) {

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

                current_loglikelihood = new_loglikelihood;

            } while (!convergence && !iterations_finished);

        }

        /*
        *   \brief Get the mean of the given component distribution.
        *   \param component_label The 0-indexed label of the component distribution.
        *   \returns Mean of the gaussian distribution.
        */
        double GetComponentMean (size_t component_label) {
            return distributions[component_label].mean ();
        }

        /*
        *   \brief Get the standard deviation of the given component distribution.
        *   \param component_label The 0-indexed label of the component distribution.
        *   \returns Standard deviation of the gaussian distribution.
        */
        double GetComponentStandardDeviation (size_t component_label) {
            return distributions[component_label].standard_deviation ();
        }

        /*
        *   \brief Get the weight of the given component distribution.
        *   \param component_label The 0-indexed label of the component distribution.
        *   \returns Weight of the gaussian distribution in the mixture model.
        */
        double GetComponentWeight (size_t component_label) {
            return weights[component_label];
        }

    private:

        /*
        *   \brief Randomly assign a set of data points as initial means of the component
        *   distributions. Modify this method to support other partition policies like KMeans in future.
        *   \param data Data on which GMM is to be trained.
        */
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

        /*
        *   \brief Evaluate the log likelihood of the current state of mixture model trained
        *   on the given data.
        *   \param data Data on which the mixture model is being trained.
        *   \param weights Weight of each component gaussian distribution.
        *   \param distributions Container holding the each of the component distributions.
        *   \returns The log likelihood of the GMM, calculated as the weighted sum of likelihoods of individual
        *   components.
        */
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