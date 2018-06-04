//  Dattatreya Mohapatra
//
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

#ifndef _BOOST_UBLAS_HISTOGRAM_
#define _BOOST_UBLAS_HISTOGRAM_
#endif

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace boost { namespace numeric { namespace ublas {

    BOOST_UBLAS_INLINE
    template<class V>
    vector<int> histogram (const V &v, int nbins = 10) {
        typedef typename V::value_type value_type;
        typedef typename V::size_type size_type;

        size_type vsize (v.size ());

        assert (vsize > 0 && "Vector is empty.");

        assert (nbins > 0 && "Number of bins should be positive.");

        double bin_size = (max (v) - min (v)) / (double)(nbins);

        assert ((bin_size > 0 || (bin_size == 0 && nbins == 1)) && "Bin size cannot be zero.");

        vector<int> bin_counts (nbins);
        for (int i = 0; i < nbins; ++ i)
            bin_counts (i) = 0;

        if (nbins == 1)
            bin_counts (0) = vsize;
        else {    
            value_type minVal = min (v);
            for (size_type i = 0; i < vsize; ++ i) {
                int bin_index = (v (i) - minVal) / bin_size;
                bin_index = (bin_index >= nbins) ? (nbins - 1) : bin_index;
                bin_counts (bin_index) += 1;
            }
        }

        return bin_counts;

    }


    BOOST_UBLAS_INLINE
    template<class V1, class V2>
    vector<int> histogram (const V1 &v, const V2 &bin_edges) {
        typedef typename V1::value_type value_type1;
        typedef typename V1::size_type size_type1;
        typedef typename V2::value_type value_type2;
        typedef typename V2::size_type size_type2;

        size_type1 size1 (v.size ());
        size_type2 size2 (bin_edges.size ());

        vector<int> bin_counts (size2);

        return bin_counts;
    }


}}}