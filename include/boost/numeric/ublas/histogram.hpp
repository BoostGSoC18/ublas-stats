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

static const double TOLER = 1e-12;

namespace boost { namespace numeric { namespace ublas {

    BOOST_UBLAS_INLINE
    template<class V>
    vector<int> histogram (const V &v, int nbins = 10) {
        typedef typename V::value_type value_type;
        typedef typename V::size_type size_type;

        size_type vsize (v.size ());

        assert (vsize > 0 && "Vector is empty.");

        assert (nbins > 0 && "Number of bins should be positive.");

        value_type bin_size = (max (v) - min (v)) / (nbins);

        assert ((bin_size > 0 || (bin_size == 0 && nbins == 1)) && "Bin size cannot be zero.");

        vector<int> bin_counts (nbins);
        for (int i = 0; i < nbins; ++ i)
            bin_counts (i) = 0;

        if (nbins == 1)
            bin_counts (0) = vsize;
        else {    
            value_type min_val = min (v);
            for (size_type i = 0; i < vsize; ++ i) {
                int bin_index = (v (i) - min_val) / bin_size;
                if ((double)(bin_index + 1) - ((v (i) - min_val) / bin_size) <= TOLER)
                    bin_index += 1;
                bin_index = (bin_index >= nbins) ? (nbins - 1) : bin_index;
                bin_counts (bin_index) += 1;
            }
        }

        return bin_counts;
    }

    BOOST_UBLAS_INLINE
    template<class V>
    vector<int> histogram (V v, const V &bin_edges) {
        typedef typename V::value_type value_type;
        typedef typename V::size_type size_type;

        size_type size1 (v.size ());

        assert (size1 > 0 && "Vector is empty.");
        
        size_type size2 (bin_edges.size ());

        assert (size2 > 1 && "Number of bins should be positive.");

        for (size_type i = 0; i < size2 - 1; ++ i)
            assert ((bin_edges (i) < bin_edges (i + 1)) &&
                "Bin edges must be monotonically increasing.");

        vector<int> bin_counts (size2 - 1);
        for (int i = 0; i < bin_counts.size (); ++ i)
            bin_counts (i) = 0;

        boost::sort (v);

        if (v (0) > bin_edges (size2 - 1) || v (size1 - 1) < bin_edges (0))
            return bin_counts;

        size_type vcounter = size_type (0);
        size_type bincounter = size_type (1);
        while (vcounter < size1 && bincounter < size2) {
            if (bin_edges (bincounter) - v (vcounter) >= TOLER) {
                bin_counts (bincounter - 1) += 1;
                ++ vcounter;
            }
            else
                ++ bincounter;
        }

        while (vcounter < size1 && v (vcounter) <= bin_edges (size2 - 1)) {
            bin_counts (size2 - 2) += 1;
            ++ vcounter;
        }

        return bin_counts;
    }

}}}