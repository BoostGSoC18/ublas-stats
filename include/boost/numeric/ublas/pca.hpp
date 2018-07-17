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

#include <boost/numeric/ublas/eigen_solver.hpp>

namespace boost { namespace numeric { namespace ublas {

    class PCA {
    public:
        PCA () {}

        template <class MatrixType>
        double Apply (const MatrixType &data,
                      matrix<double> &transformed_data,
                      const int n_components) {
            
            return Apply (data, transformed_data, (size_t)n_components);
        }

        template <class MatrixType>
        double Apply (const MatrixType &data,
                      matrix<double> &transformed_data,
                      const size_t n_components) {

            assert (n_components > 0);
            assert (n_components <= data.size2 ());

            matrix<double> transformed_data_temp (data.size1 (), data.size2 ());
            vector<double> eigen_values (data.size2 ());
            matrix<double> eigen_vectors (data.size2 (), data.size2 ());
            
            Apply (data, transformed_data_temp, eigen_values, eigen_vectors);
            
            transformed_data = subrange (transformed_data_temp, 0, data.size1 (), 0, n_components);
        }

        template <class MatrixType>
        void Apply (const MatrixType &data,
                    matrix<double> &transformed_data,
                    const double variance_retained) {

            matrix<double> transformed_data_temp (data.size1 (), data.size2 ());
            vector<double> eigen_values (data.size2 ());
            matrix<double> eigen_vectors (data.size2 (), data.size2 ());
            
            Apply (data, transformed_data_temp, eigen_values, eigen_vectors);
            
            size_t eig_counter = 0;
            double cur_eig_val_sum = 0;
            double total_eig_val_sum = sum (eigen_values);
            while (eig_counter < eigen_values.size ()) {
                cur_eig_val_sum += eigen_values (eig_counter);
                double cumulative_variance = cur_eig_val_sum / total_eig_val_sum;
                if (cumulative_variance >= variance_retained)
                    break;
                ++ eig_counter;
            }

            transformed_data = subrange (transformed_data_temp, 0, data.size1 (), 0, eig_counter + 1);
        }

        template <class MatrixType>
        void Apply (const MatrixType &data,
                    matrix<double> &transformed_data,
                    vector<double> &eigen_values,
                    matrix<double> &eigen_vectors) {

            matrix<double> centered_data (data.size1 (), data.size2 ());
            CenterData (data, centered_data);

            size_t n = data.size1 ();
            MatrixType empirical_covariance_matrix = (1 / (n - 1)) * (trans (data) * data);

            // Code to get eigen values and eigen vectors.
            eigen_solver< matrix<double> > evd_solver (empirical_covariance_matrix, EIGVEC);
            matrix<double> eigval_matrix = evd_solver.get_real_eigenvalues ();
            for (size_t i = 0; i < eigval_matrix.size1 (); ++ i)
                eigen_values (i) = eigval_matrix (i, i);
            eigen_vectors = evd_solver.get_real_eigenvectors ();

            transformed_data = centered_data * eigen_vectors;
        }

    private:

        template <class MatrixType>
        void CenterData (const MatrixType &data, matrix<double> &centered_data) {
            centered_data = trans (scalar_vector<int> (1)) * mean (data, 0);
        }

    };

}}}