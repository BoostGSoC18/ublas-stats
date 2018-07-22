// Dattatreya Mohapatra
// 
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

/// \file pca.hpp Contains the class and methods for performing PCA on matrices.

#ifndef _BOOST_UBLAS_PCA_
#define _BOOST_UBLAS_PCA_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/eigen_solver.hpp>

#include <boost/range/algorithm.hpp>
#include <boost/tuple/tuple.hpp>

namespace boost { namespace numeric { namespace ublas {

    /*
    *   \brief Implements the PCA class. Currently, principal factors are calculated by
    *   performing EVD on the empirical covariance matrix. Other decomposition policies
    *   can be plugged as a template parameter.
    */
    class PCA {
    public:
        PCA () {}

        /*
        *   \brief Overload for Apply (const MatrixType, matrix<double>, const size_t)
        *   to handle integer type n_components.
        */
        template <class MatrixType>
        double Apply (const MatrixType &data,
                      matrix<double> &transformed_data,
                      const int n_components) {

            return Apply (data, transformed_data, (size_t)n_components);
        }

        /*
        *   \brief Performs PCA and returns the principal components representing maximum
        *   variance from data
        *   \tparam MatrixType The type of data points (int, float, double etc).
        *   \param data Data on which PCA is to be performed.
        *   \param transformed_data Container to store and return the transformed data
        *   (principal factors as columns).
        *   \param n_components The number of principal components to return.
        */
        template <class MatrixType>
        double Apply (const MatrixType &data,
                      matrix<double> &transformed_data,
                      const size_t n_components) {

            assert (n_components > 0);
            assert (n_components <= data.size2 ());

            matrix<double> transformed_data_temp (data.size1 (), data.size2 ());
            
            Apply (data, transformed_data_temp);
            
            transformed_data = subrange (transformed_data_temp, 0, data.size1 (), 0, n_components);
            
            double variance_retained = sum (subslice (EIGEN_VALUES, 0, 1, n_components)) / sum (EIGEN_VALUES);

            return variance_retained;
        }

        /*
        *   \brief Performs PCA and returns the minimal set of principal components
        *   whose cumulative represented variance is at least variance_retained.
        *   \tparam MatrixType The type of data points (int, float, double etc).
        *   \param data Data on which PCA is to be performed.
        *   \param transformed_data Container to store and return the transformed data
        *   (principal factors as columns).
        *   \param variance_retained The minimum amount of variance to retain from data.
        */
        template <class MatrixType>
        double Apply (const MatrixType &data,
                    matrix<double> &transformed_data,
                    const double variance_retained) {

            matrix<double> transformed_data_temp (data.size1 (), data.size2 ());
            
            Apply (data, transformed_data_temp);
            
            size_t eig_counter = 0;
            double cur_eig_val_sum = 0;
            double total_eig_val_sum = sum (EIGEN_VALUES);
            while (eig_counter < EIGEN_VALUES.size ()) {
                cur_eig_val_sum += EIGEN_VALUES (eig_counter);
                double cumulative_variance = cur_eig_val_sum / total_eig_val_sum;
                if (cumulative_variance >= variance_retained)
                    break;
                ++ eig_counter;
            }

            transformed_data = subrange (transformed_data_temp, 0, data.size1 (), 0, eig_counter + 1);

            return cur_eig_val_sum / total_eig_val_sum;
        }

        /*
        *   \brief The actual method called by other Apply methods to perform PCA on
        *   data. Principal components are calculated by first calculating eigenvalues
        *   and eigenvectors of the empirical covariance matrix, then transforming data
        *   using eigenvector matrix.
        *   \tparam MatrixType The type of data points (int, float, double etc).
        *   \param data Data on which PCA is to be performed.
        *   \param transformed_data Container to store and return the transformed data
        *   (principal factors as columns).
        */
        template <class MatrixType>
        void Apply (const MatrixType &data,
                    matrix<double> &transformed_data) {

            // Initialise eigenvalue and eigenvector containers.
            EIGEN_VALUES = vector<double> (data.size2 ());
            EIGEN_VECTORS = matrix<double> (data.size2 (), data.size2 ());
            
            matrix<double> centered_data (data.size1 (), data.size2 ());
            
            // Center each column by subtracting from its mean. C = D - mean (D, axis=0)
            CenterData (data, centered_data);

            size_t n = data.size1 ();
            
            // Construct the empirical covariance matrix. B = (C.T * C) / (n - 1)
            MatrixType empirical_covariance_matrix = prod (trans (centered_data), centered_data) / (n - 1);
            
            // Get eigenvalues and eigenvectors using the eigensolver implemented in uBLAS.
            eigen_solver< matrix<double> > evd_solver (empirical_covariance_matrix, EIGVEC);
            matrix<double> eigval_matrix = evd_solver.get_real_eigenvalues ();
            
            // The eigenvectors returned by the eigensolver are not sorted. They need to be
            // sorted in order to return the set of components with maximum variance.
            vector< boost::tuple<size_t, double> > eigval_index_tuples (EIGEN_VALUES.size ());
            for (size_t i = 0; i < EIGEN_VALUES.size (); ++ i)
                eigval_index_tuples (i) = boost::make_tuple (eigval_matrix (i, i), i);
            boost::sort (eigval_index_tuples, CompareTuples);

            for (size_t i = 0; i < eigval_index_tuples.size (); ++ i)
                EIGEN_VALUES (i) = eigval_index_tuples (i).get<0> ();

            matrix<double> eigvecs_temp = evd_solver.get_real_eigenvectors ();
            
            // Order the eigenvectors according to the sorted order of their corresponding
            // eigenvalues.
            for (size_t i = 0; i < eigval_index_tuples.size (); ++ i)
                column (EIGEN_VECTORS, i) = column (eigvecs_temp, eigval_index_tuples (i).get<1> ());

            transformed_data = prod (centered_data, EIGEN_VECTORS);
        }

    public:

        /*
        *   \brief Getter method for eigenvalues of the covariance matrix.
        */
        vector<double> GetEigenValues () const {
            return EIGEN_VALUES;
        }

        /*
        *   \brief Getter method for eigenvectors of the covariance matrix.
        */
        matrix<double> GetEigenVectors () const {
            return EIGEN_VECTORS;
        }

    private:

        /*
        *   \brief Centers each column of the data matrix around its mean.
        *   \tparam MatrixType The type of data points (int, float, double etc).
        *   \param data Data to be centered.
        *   \param centered_data Container to store and return the centered data.
        */
        template <class MatrixType>
        void CenterData (const MatrixType &data, matrix<double> &centered_data) {
            for (size_t i = 0; i < data.size2 (); ++ i)
                column (centered_data, i) = column (data, i) - scalar_vector<double> (data.size1 (), mean (column (data, i)));
        }

        static bool CompareTuples (const boost::tuple<double, size_t> &i, const boost::tuple<double, size_t> &j) {
            return (i.get<0> () > j.get<0> ());
        }

    private:

        vector<double> EIGEN_VALUES;
        matrix<double> EIGEN_VECTORS;

    };

}}}

#endif