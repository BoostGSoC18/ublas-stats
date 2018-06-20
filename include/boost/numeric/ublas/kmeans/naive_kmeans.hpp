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

#ifndef _BOOST_UBLAS_NAIVE_KMEANS_
#define _BOOST_UBLAS_NAIVE_KMEANS_


#include <boost/numeric/ublas/matrix.hpp>

namespace boost { namespace numeric { namespace ublas {

    class NaiveKMeans {
    public:
        template </*class MetricType,*/ class MatrixType>
        NaiveKMeans (MatrixType &data/*, MetricType distance_metric*/) {}

        double Iterate (matrix<double> &centroids, matrix<double>& new_centroids) {}

    private:
        MatrixType data;
        /*MetricType distance_metric;*/
        
    }
}

#endif