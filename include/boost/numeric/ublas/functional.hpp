//
//  Copyright (c) 2000-2009
//  Joerg Walter, Mathias Koch, Gunter Winkler
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.
//

#ifndef _BOOST_UBLAS_FUNCTIONAL_
#define _BOOST_UBLAS_FUNCTIONAL_

#include <functional>

#include <boost/core/ignore_unused.hpp>
#include <boost/unordered_map.hpp>
#include <boost/range/algorithm.hpp>

#include <boost/numeric/ublas/traits.hpp>
#ifdef BOOST_UBLAS_USE_DUFF_DEVICE
#include <boost/numeric/ublas/detail/duff.hpp>
#endif
#ifdef BOOST_UBLAS_USE_SIMD
#include <boost/numeric/ublas/detail/raw.hpp>
#else
namespace boost { namespace numeric { namespace ublas { namespace raw {
}}}}
#endif
#ifdef BOOST_UBLAS_HAVE_BINDINGS
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/atlas/cblas.hpp>
#endif

#include <boost/numeric/ublas/detail/definitions.hpp>

#include <boost/numeric/ublas/vector.hpp>



namespace boost { namespace numeric { namespace ublas {

    // Scalar functors

    // Unary
    template<class T>
    struct scalar_unary_functor {
        typedef T value_type;
        typedef typename type_traits<T>::const_reference argument_type;
        typedef typename type_traits<T>::value_type result_type;
    };

    template<class T>
    struct scalar_identity:
        public scalar_unary_functor<T> {
        typedef typename scalar_unary_functor<T>::argument_type argument_type;
        typedef typename scalar_unary_functor<T>::result_type result_type;

        static BOOST_UBLAS_INLINE
        result_type apply (argument_type t) {
            return t;
        }
    };
    template<class T>
    struct scalar_negate:
        public scalar_unary_functor<T> {
        typedef typename scalar_unary_functor<T>::argument_type argument_type;
        typedef typename scalar_unary_functor<T>::result_type result_type;

        static BOOST_UBLAS_INLINE
        result_type apply (argument_type t) {
            return - t;
        }
    };
    template<class T>
    struct scalar_conj:
        public scalar_unary_functor<T> {
        typedef typename scalar_unary_functor<T>::value_type value_type;
        typedef typename scalar_unary_functor<T>::argument_type argument_type;
        typedef typename scalar_unary_functor<T>::result_type result_type;

        static BOOST_UBLAS_INLINE
        result_type apply (argument_type t) {
            return type_traits<value_type>::conj (t);
        }
    };

    // Unary returning real
    template<class T>
    struct scalar_real_unary_functor {
        typedef T value_type;
        typedef typename type_traits<T>::const_reference argument_type;
        typedef typename type_traits<T>::real_type result_type;
    };

    template<class T>
    struct scalar_real:
        public scalar_real_unary_functor<T> {
        typedef typename scalar_real_unary_functor<T>::value_type value_type;
        typedef typename scalar_real_unary_functor<T>::argument_type argument_type;
        typedef typename scalar_real_unary_functor<T>::result_type result_type;

        static BOOST_UBLAS_INLINE
        result_type apply (argument_type t) {
            return type_traits<value_type>::real (t);
        }
    };
    template<class T>
    struct scalar_imag:
        public scalar_real_unary_functor<T> {
        typedef typename scalar_real_unary_functor<T>::value_type value_type;
        typedef typename scalar_real_unary_functor<T>::argument_type argument_type;
        typedef typename scalar_real_unary_functor<T>::result_type result_type;

        static BOOST_UBLAS_INLINE
        result_type apply (argument_type t) {
            return type_traits<value_type>::imag (t);
        }
    };

    // Binary
    template<class T1, class T2>
    struct scalar_binary_functor {
        typedef typename type_traits<T1>::const_reference argument1_type;
        typedef typename type_traits<T2>::const_reference argument2_type;
        typedef typename promote_traits<T1, T2>::promote_type result_type;
    };

    template<class T1, class T2>
    struct scalar_plus:
        public scalar_binary_functor<T1, T2> {
        typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
        typedef typename scalar_binary_functor<T1, T2>::result_type result_type;

        static BOOST_UBLAS_INLINE
        result_type apply (argument1_type t1, argument2_type t2) {
            return t1 + t2;
        }
    };
    template<class T1, class T2>
    struct scalar_minus:
        public scalar_binary_functor<T1, T2> {
        typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
        typedef typename scalar_binary_functor<T1, T2>::result_type result_type;

        static BOOST_UBLAS_INLINE
        result_type apply (argument1_type t1, argument2_type t2) {
            return t1 - t2;
        }
    };
    template<class T1, class T2>
    struct scalar_multiplies:
        public scalar_binary_functor<T1, T2> {
        typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
        typedef typename scalar_binary_functor<T1, T2>::result_type result_type;

        static BOOST_UBLAS_INLINE
        result_type apply (argument1_type t1, argument2_type t2) {
            return t1 * t2;
        }
    };
    template<class T1, class T2>
    struct scalar_divides:
        public scalar_binary_functor<T1, T2> {
        typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
        typedef typename scalar_binary_functor<T1, T2>::result_type result_type;

        static BOOST_UBLAS_INLINE
        result_type apply (argument1_type t1, argument2_type t2) {
            return t1 / t2;
        }
    };

    template<class T1, class T2>
    struct scalar_binary_assign_functor {
        // ISSUE Remove reference to avoid reference to reference problems
        typedef typename type_traits<typename boost::remove_reference<T1>::type>::reference argument1_type;
        typedef typename type_traits<T2>::const_reference argument2_type;
    };

    struct assign_tag {};
    struct computed_assign_tag {};

    template<class T1, class T2>
    struct scalar_assign:
        public scalar_binary_assign_functor<T1, T2> {
        typedef typename scalar_binary_assign_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_assign_functor<T1, T2>::argument2_type argument2_type;
#if BOOST_WORKAROUND( __IBMCPP__, <=600 )
        static const bool computed ;
#else
        static const bool computed = false ;
#endif

        static BOOST_UBLAS_INLINE
        void apply (argument1_type t1, argument2_type t2) {
            t1 = t2;
        }

        template<class U1, class U2>
        struct rebind {
            typedef scalar_assign<U1, U2> other;
        };
    };

#if BOOST_WORKAROUND( __IBMCPP__, <=600 )
    template<class T1, class T2>
    const bool scalar_assign<T1,T2>::computed = false;
#endif

    template<class T1, class T2>
    struct scalar_plus_assign:
        public scalar_binary_assign_functor<T1, T2> {
        typedef typename scalar_binary_assign_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_assign_functor<T1, T2>::argument2_type argument2_type;
#if BOOST_WORKAROUND( __IBMCPP__, <=600 )
        static const bool computed ;
#else
        static const bool computed = true ;
#endif

        static BOOST_UBLAS_INLINE
        void apply (argument1_type t1, argument2_type t2) {
            t1 += t2;
        }

        template<class U1, class U2>
        struct rebind {
            typedef scalar_plus_assign<U1, U2> other;
        };
    };

#if BOOST_WORKAROUND( __IBMCPP__, <=600 )
    template<class T1, class T2>
    const bool scalar_plus_assign<T1,T2>::computed = true;
#endif

    template<class T1, class T2>
    struct scalar_minus_assign:
        public scalar_binary_assign_functor<T1, T2> {
        typedef typename scalar_binary_assign_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_assign_functor<T1, T2>::argument2_type argument2_type;
#if BOOST_WORKAROUND( __IBMCPP__, <=600 )
        static const bool computed ;
#else
        static const bool computed = true ;
#endif

        static BOOST_UBLAS_INLINE
        void apply (argument1_type t1, argument2_type t2) {
            t1 -= t2;
        }

        template<class U1, class U2>
        struct rebind {
            typedef scalar_minus_assign<U1, U2> other;
        };
    };

#if BOOST_WORKAROUND( __IBMCPP__, <=600 )
    template<class T1, class T2>
    const bool scalar_minus_assign<T1,T2>::computed = true;
#endif

    template<class T1, class T2>
    struct scalar_multiplies_assign:
        public scalar_binary_assign_functor<T1, T2> {
        typedef typename scalar_binary_assign_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_assign_functor<T1, T2>::argument2_type argument2_type;
        static const bool computed = true;

        static BOOST_UBLAS_INLINE
        void apply (argument1_type t1, argument2_type t2) {
            t1 *= t2;
        }

        template<class U1, class U2>
        struct rebind {
            typedef scalar_multiplies_assign<U1, U2> other;
        };
    };
    template<class T1, class T2>
    struct scalar_divides_assign:
        public scalar_binary_assign_functor<T1, T2> {
        typedef typename scalar_binary_assign_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_assign_functor<T1, T2>::argument2_type argument2_type;
        static const bool computed ;

        static BOOST_UBLAS_INLINE
        void apply (argument1_type t1, argument2_type t2) {
            t1 /= t2;
        }

        template<class U1, class U2>
        struct rebind {
            typedef scalar_divides_assign<U1, U2> other;
        };
    };
    template<class T1, class T2>
    const bool scalar_divides_assign<T1,T2>::computed = true;

    template<class T1, class T2>
    struct scalar_binary_swap_functor {
        typedef typename type_traits<typename boost::remove_reference<T1>::type>::reference argument1_type;
        typedef typename type_traits<typename boost::remove_reference<T2>::type>::reference argument2_type;
    };

    template<class T1, class T2>
    struct scalar_swap:
        public scalar_binary_swap_functor<T1, T2> {
        typedef typename scalar_binary_swap_functor<T1, T2>::argument1_type argument1_type;
        typedef typename scalar_binary_swap_functor<T1, T2>::argument2_type argument2_type;

        static BOOST_UBLAS_INLINE
        void apply (argument1_type t1, argument2_type t2) {
            std::swap (t1, t2);
        }

        template<class U1, class U2>
        struct rebind {
            typedef scalar_swap<U1, U2> other;
        };
    };

    // Vector functors

    // Unary returning scalar
    template<class V>
    struct vector_scalar_unary_functor {
        typedef typename V::value_type value_type;
        typedef typename V::value_type result_type;
    };

    template<class V>
    struct vector_min: 
        public vector_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) { 
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            result_type min_val = result_type (e () (vector_size_type (0)));
            for (vector_size_type i = 0; i < size; ++ i)
                if ( type_traits<value_type>::type_abs(min_val) > type_traits<value_type>::type_abs(e () (i)) )
                    min_val = e () (i);
            return min_val;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type min_val = *it;
            while (-- size >= 0) {
                if ( type_traits<value_type>::type_abs(min_val) > type_traits<value_type>::type_abs(*it) )
                    min_val = *it;
                ++ it;
            }
            return min_val;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type min_val = *it;
            typedef typename I::difference_type vector_difference_type;
            vector_difference_type size (0);
            while (it != it_end) {
                if ( type_traits<value_type>::type_abs(min_val) > type_traits<value_type>::type_abs(*it) )
                    min_val = *it;
                ++ it;
                ++ size;
            }
            return min_val;
        }
    };

    template<class V>
    struct vector_max:
        public vector_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) { 
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            result_type max_val = result_type (e () (vector_size_type (0)));
            for (vector_size_type i = 0; i < size; ++ i)
                if ( type_traits<value_type>::type_abs(max_val) < type_traits<value_type>::type_abs(e () (i)) )
                    max_val = e () (i);
            return max_val;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type max_val = *it;
            while (-- size >= 0) {
                if ( type_traits<value_type>::type_abs(max_val) < type_traits<value_type>::type_abs(*it) )
                    max_val = *it;
                ++ it;
            }
            return max_val;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type max_val = *it;
            typedef typename I::difference_type vector_difference_type;
            vector_difference_type size (0);
            while (it != it_end) {
                if ( type_traits<value_type>::type_abs(max_val) < type_traits<value_type>::type_abs(*it) )
                    max_val = *it;
                ++ it;
                ++ size;
            }
            return max_val;
        }
    };

    template<class V>
    struct vector_sum:
        public vector_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) { 
            result_type t = result_type (0);
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            for (vector_size_type i = 0; i < size; ++ i)
                t += e () (i);
            return t;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type t = result_type (0);
            while (-- size >= 0) {
                t += *it;
                ++ it;
            }
            return t; 
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type t = result_type (0);
            while (it != it_end) {
                t += *it;
                ++ it;
            }
            return t; 
        }
    };

    template<class V>
    struct vector_mean: 
        public vector_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;
        // typedef double result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) { 
            result_type t = result_type (0);
            // std::cout << t << std::endl;
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            for (vector_size_type i = 0; i < size; ++ i)
                t += e () (i);
            return t / size;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type t = result_type (0);
            while (-- size >= 0) {
                t += *it;
                ++ it;
            }
            return t / size; 
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type t = result_type (0);
            typedef typename I::difference_type vector_difference_type;
            vector_difference_type size (0);
            while (it != it_end) {
                t += *it;
                ++ it;
                ++ size;
            }
            return t / size;
        }
    };

    template<class V>
    struct vector_mean_iterative: 
        public vector_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;
        // typedef double result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) { 
            result_type t = result_type (0);
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            for (vector_size_type i = 0; i < size; ++ i)
                t += (e () (i) - t) / (i + 1);
            return t;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type t = result_type (0);
            D i (0);
            while (++i <= size) {
                t += (*it - t) / i;
                ++ it;
            }
            return t;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type t = result_type (0);
            typedef typename I::difference_type vector_difference_type;
            vector_difference_type size (1);
            while (it != it_end) {
                t += (*it - t) / size;
                ++ it;
                ++ size;
            }
            return t;
        }
    };

    template<class V>
    struct vector_variance: 
        public vector_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;
        // typedef double result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) { 
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            result_type sumsq = result_type (0);
            result_type sum = result_type (0);
            for (vector_size_type i = 0; i < size; ++ i) {
                sumsq += e () (i) * e () (i);
                sum += e () (i);
            }
            return (sumsq - (sum * sum) / size) / size;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type sumsq = result_type (0);
            result_type sum = result_type (0);
            D i (0);
            while (++i <= size) {
                sumsq += *it * *it;
                sum += *it;
                ++ it;
            }
            return (sumsq - (sum * sum) / size) / size;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type sumsq = result_type (0);
            result_type sum = result_type (0);
            typedef typename I::difference_type vector_difference_type;
            vector_difference_type size (0);
            while (it != it_end) {
                sumsq += *it * *it;
                sum += *it;
                ++ it;
                ++ size;
            }
            return (sumsq - (sum * sum) / size) / size;
        }
    };

    template<class V>
    struct vector_variance_iterative: 
        public vector_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;
        // typedef double result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) { 
            result_type mean = result_type (0);
            result_type var = result_type (0);
            result_type del;
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            for (vector_size_type i = 0; i < size; ++ i) {
                del = e () (i) - mean;
                mean += del / (i + 1);
                var += del * (e () (i) - mean);
            }
            return var / size;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type mean = result_type (0);
            result_type n = result_type (0);
            result_type var = result_type (0);
            result_type del;
            while (++ n <= size){
                del = *it - mean;
                mean += del / n;
                var += del * (*it - mean);
                ++ it;
            }
            return var / size; 
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            typedef typename I::difference_type vector_difference_type;
            result_type mean = result_type (0);
            vector_difference_type n (0);
            result_type var = result_type (0);
            result_type del;
            while (it != it_end)  {
                ++ n;
                del = *it - mean;
                mean += del / n;
                var += del * (*it - mean);
                ++ it;
            }
            return var / n; 
        }
    };

    template<class V>
    struct vector_mode: 
        public vector_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;
        
        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) {
            typedef typename E::size_type vector_size_type;
            boost::unordered_map<result_type, vector_size_type> count_map (0);
            vector_size_type size (e ().size ());
            typename boost::unordered_map<result_type, vector_size_type>::iterator p;
            // auto p = count_map.iterator;
            for (vector_size_type i = 0; i < size; ++ i) {
                p = count_map.find (e () (i));
                if (p != count_map.end ())
                    p->second += vector_size_type (1);
                else
                    count_map.emplace (e () (i), vector_size_type (1));
            }
            
            result_type mode = result_type (0);
            vector_size_type mode_val = vector_size_type (0);
            p = count_map.begin ();
            while (p != count_map.end ()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs (p->first) < type_traits<value_type>::type_abs (mode))
                    mode = p->first;
                ++ p;
            }
            
            return mode;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            boost::unordered_map<result_type, D> count_map;
            
            typename boost::unordered_map<result_type, D>::iterator p;
            // auto p = count_map.iterator;
            while (-- size >= 0) {
                p = count_map.find (*it);
                if ( p != count_map.end ())
                    p->second += D (1);
                else
                    count_map.emplace (*it, D (1));
                ++ it;
            }
            
            result_type mode = result_type (0);
            D mode_val = D (0);
            p = count_map.begin ();
            while (p != count_map.end ()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs (p->first) < type_traits<value_type>::type_abs (mode))
                    mode = p->first;
                ++ p;
            }
            
            return mode;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            typedef typename I::difference_type vector_difference_type;

            boost::unordered_map<result_type, vector_difference_type> count_map;
            
            typename boost::unordered_map<result_type, vector_difference_type>::iterator p;
            // auto p = count_map.iterator;

            while (it != it_end) {
                p = count_map.find (*it);
                if ( p != count_map.end ())
                    p->second += vector_difference_type (1);
                else
                    count_map.emplace (*it, vector_difference_type (1));
                ++ it;
            }
            
            result_type mode = result_type (0);
            vector_difference_type mode_val = vector_difference_type (0);
            p = count_map.begin();
            while (p != count_map.end()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs (p->first) < type_traits<value_type>::type_abs (mode))
                    mode = p->first;
                ++ p;
            }
            
            return mode;
        }
    };

    template<class V>
    struct vector_median:
        public vector_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;
        
        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) {
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            vector<value_type> v (size);
            for (int i = 0; i < size; ++ i)
                v (i) = e () (i);
            boost::sort (v, compareElement);
            if (size % 2)
                return v (size / 2);
            else
                return (v (size / 2) + v ((size / 2) - 1)) / 2;
        }

        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            vector<value_type> v (size);
            D size_temp = size;
            D i = D (0);
            while (-- size_temp >= 0) {
                v (i) = *it;
                ++ it;
                ++  i;
            }
            boost::sort (v, compareElement);
            if (size % 2)
                return v (size / 2);
            else
                return (v (size / 2) + v ((size / 2) - 1)) / 2;
        }

        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            typedef typename I::difference_type vector_difference_type;
            vector_difference_type size = vector_difference_type (0);
            I it1 = it;
            while (it1 != it_end) {
                ++ it1;
                ++ size;
            }
            vector<value_type> v (size);
            vector_difference_type i = vector_difference_type (0);
            while (it != it_end) {
                v (i) = *it;
                ++ it;
                ++ i;
            }
            
            boost::sort (v, compareElement);
            if (size % 2)
                return v (size / 2);
            else
                return (v (size / 2) + v ((size / 2) - 1)) / 2;
        }

        static BOOST_UBLAS_INLINE
        bool compareElement (const value_type &A, const value_type &B) {
            return (type_traits<value_type>::type_abs (A) < type_traits<value_type>::type_abs (B));
        }
    };

    // Unary returning real scalar 
    template<class V>
    struct vector_scalar_real_unary_functor {
        typedef typename V::value_type value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef real_type result_type;
    };

    template<class V>
    struct vector_norm_1:
        public vector_scalar_real_unary_functor<V> {
        typedef typename vector_scalar_real_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_real_unary_functor<V>::real_type real_type;
        typedef typename vector_scalar_real_unary_functor<V>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) {
            real_type t = real_type ();
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            for (vector_size_type i = 0; i < size; ++ i) {
                real_type u (type_traits<value_type>::type_abs (e () (i)));
                t += u;
            }
            return t;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            real_type t = real_type ();
            while (-- size >= 0) {
                real_type u (type_traits<value_type>::norm_1 (*it));
                t += u;
                ++ it;
            }
            return t;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            real_type t = real_type ();
            while (it != it_end) {
                real_type u (type_traits<value_type>::norm_1 (*it));
                t += u;
                ++ it;
            }
            return t;
        }
    };
    template<class V>
    struct vector_norm_2:
        public vector_scalar_real_unary_functor<V> {
        typedef typename vector_scalar_real_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_real_unary_functor<V>::real_type real_type;
        typedef typename vector_scalar_real_unary_functor<V>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) {
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
#ifndef BOOST_UBLAS_SCALED_NORM
            real_type t = real_type ();
            for (vector_size_type i = 0; i < size; ++ i) {
                real_type u (type_traits<value_type>::norm_2 (e () (i)));
                t +=  u * u;
            }
            return static_cast<result_type>(type_traits<real_type>::type_sqrt (t));
#else
            real_type scale = real_type ();
            real_type sum_squares (1);
            for (vector_size_type i = 0; i < size; ++ i) {
                real_type u (type_traits<value_type>::norm_2 (e () (i)));
                if ( real_type () /* zero */ == u ) continue;
                if (scale < u) {
                    real_type v (scale / u);
                    sum_squares = sum_squares * v * v + real_type (1);
                    scale = u;
                } else {
                    real_type v (u / scale);
                    sum_squares += v * v;
                }
            }
            return static_cast<result_type>(scale * type_traits<real_type>::type_sqrt (sum_squares));
#endif
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
#ifndef BOOST_UBLAS_SCALED_NORM
            real_type t = real_type ();
            while (-- size >= 0) {
                real_type u (type_traits<value_type>::norm_2 (*it));
                t +=  u * u;
                ++ it;
            }
            return static_cast<result_type>(type_traits<real_type>::type_sqrt (t));
#else
            real_type scale = real_type ();
            real_type sum_squares (1);
            while (-- size >= 0) {
                real_type u (type_traits<value_type>::norm_2 (*it));
                if (scale < u) {
                    real_type v (scale / u);
                    sum_squares = sum_squares * v * v + real_type (1);
                    scale = u;
                } else {
                    real_type v (u / scale);
                    sum_squares += v * v;
                }
                ++ it;
            }
            return static_cast<result_type>(scale * type_traits<real_type>::type_sqrt (sum_squares));
#endif
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
#ifndef BOOST_UBLAS_SCALED_NORM
            real_type t = real_type ();
            while (it != it_end) {
                real_type u (type_traits<value_type>::norm_2 (*it));
                t +=  u * u;
                ++ it;
            }
            return static_cast<result_type>(type_traits<real_type>::type_sqrt (t));
#else
            real_type scale = real_type ();
            real_type sum_squares (1);
            while (it != it_end) {
                real_type u (type_traits<value_type>::norm_2 (*it));
                if (scale < u) {
                    real_type v (scale / u);
                    sum_squares = sum_squares * v * v + real_type (1);
                    scale = u;
                } else {
                    real_type v (u / scale);
                    sum_squares += v * v;
                }
                ++ it;
            }
            return static_cast<result_type>(scale * type_traits<real_type>::type_sqrt (sum_squares));
#endif
        }
    };

    template<class V>
    struct vector_norm_2_square :
        public vector_scalar_real_unary_functor<V> {
        typedef typename vector_scalar_real_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_real_unary_functor<V>::real_type real_type;
        typedef typename vector_scalar_real_unary_functor<V>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) {
            real_type t = real_type ();
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            for (vector_size_type i = 0; i < size; ++ i) {
                real_type u (type_traits<value_type>::norm_2 (e () (i)));
                t +=  u * u;
            }
            return t;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            real_type t = real_type ();
            while (-- size >= 0) {
                real_type u (type_traits<value_type>::norm_2 (*it));
                t +=  u * u;
                ++ it;
            }
            return t;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            real_type t = real_type ();
            while (it != it_end) {
                real_type u (type_traits<value_type>::norm_2 (*it));
                t +=  u * u;
                ++ it;
            }
            return t;
        }
    };

    template<class V>
    struct vector_norm_inf:
        public vector_scalar_real_unary_functor<V> {
        typedef typename vector_scalar_real_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_real_unary_functor<V>::real_type real_type;
        typedef typename vector_scalar_real_unary_functor<V>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) {
            real_type t = real_type ();
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            for (vector_size_type i = 0; i < size; ++ i) {
                real_type u (type_traits<value_type>::norm_inf (e () (i)));
                if (u > t)
                    t = u;
            }
            return t;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            real_type t = real_type ();
            while (-- size >= 0) {
                real_type u (type_traits<value_type>::norm_inf (*it));
                if (u > t)
                    t = u;
                ++ it;
            }
            return t;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) { 
            real_type t = real_type ();
            while (it != it_end) {
                real_type u (type_traits<value_type>::norm_inf (*it));
                if (u > t) 
                    t = u;
                ++ it;
            }
            return t; 
        }
    };

    // Unary returning index
    template<class V>
    struct vector_scalar_index_unary_functor {
        typedef typename V::value_type value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename V::size_type result_type;
    };

    template<class V>
    struct vector_index_norm_inf:
        public vector_scalar_index_unary_functor<V> {
        typedef typename vector_scalar_index_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_index_unary_functor<V>::real_type real_type;
        typedef typename vector_scalar_index_unary_functor<V>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E> &e) {
            // ISSUE For CBLAS compatibility return 0 index in empty case
            result_type i_norm_inf (0);
            real_type t = real_type ();
            typedef typename E::size_type vector_size_type;
            vector_size_type size (e ().size ());
            for (vector_size_type i = 0; i < size; ++ i) {
                real_type u (type_traits<value_type>::norm_inf (e () (i)));
                if (u > t) {
                    i_norm_inf = i;
                    t = u;
                }
            }
            return i_norm_inf;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            // ISSUE For CBLAS compatibility return 0 index in empty case
            result_type i_norm_inf (0);
            real_type t = real_type ();
            while (-- size >= 0) {
                real_type u (type_traits<value_type>::norm_inf (*it));
                if (u > t) {
                    i_norm_inf = it.index ();
                    t = u;
                }
                ++ it;
            }
            return i_norm_inf;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            // ISSUE For CBLAS compatibility return 0 index in empty case
            result_type i_norm_inf (0);
            real_type t = real_type ();
            while (it != it_end) {
                real_type u (type_traits<value_type>::norm_inf (*it));
                if (u > t) {
                    i_norm_inf = it.index ();
                    t = u;
                }
                ++ it;
            }
            return i_norm_inf;
        }
    };

    // Binary returning scalar
    template<class V1, class V2, class TV>
    struct vector_scalar_binary_functor {
        typedef TV value_type;
        typedef TV result_type;
    };

    template<class V1, class V2, class TV>
    struct vector_inner_prod:
        public vector_scalar_binary_functor<V1, V2, TV> {
        typedef typename vector_scalar_binary_functor<V1, V2, TV>::value_type value_type;
        typedef typename vector_scalar_binary_functor<V1, V2, TV>::result_type result_type;

        template<class C1, class C2>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_container<C1> &c1,
                           const vector_container<C2> &c2) {
#ifdef BOOST_UBLAS_USE_SIMD
            using namespace raw;
            typedef typename C1::size_type vector_size_type;
            vector_size_type size (BOOST_UBLAS_SAME (c1 ().size (), c2 ().size ()));
            const typename V1::value_type *data1 = data_const (c1 ());
            const typename V1::value_type *data2 = data_const (c2 ());
            vector_size_type s1 = stride (c1 ());
            vector_size_type s2 = stride (c2 ());
            result_type t = result_type (0);
            if (s1 == 1 && s2 == 1) {
                for (vector_size_type i = 0; i < size; ++ i)
                    t += data1 [i] * data2 [i];
            } else if (s2 == 1) {
                for (vector_size_type i = 0, i1 = 0; i < size; ++ i, i1 += s1)
                    t += data1 [i1] * data2 [i];
            } else if (s1 == 1) {
                for (vector_size_type i = 0, i2 = 0; i < size; ++ i, i2 += s2)
                    t += data1 [i] * data2 [i2];
            } else {
                for (vector_size_type i = 0, i1 = 0, i2 = 0; i < size; ++ i, i1 += s1, i2 += s2)
                    t += data1 [i1] * data2 [i2];
            }
            return t;
#elif defined(BOOST_UBLAS_HAVE_BINDINGS)
            return boost::numeric::bindings::atlas::dot (c1 (), c2 ());
#else
            return apply (static_cast<const vector_expression<C1> > (c1), static_cast<const vector_expression<C2> > (c2));
#endif
        }
        template<class E1, class E2>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E1> &e1,
                           const vector_expression<E2> &e2) {
            typedef typename E1::size_type vector_size_type;
            vector_size_type size (BOOST_UBLAS_SAME (e1 ().size (), e2 ().size ()));
            result_type t = result_type (0);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            for (vector_size_type i = 0; i < size; ++ i)
                t += e1 () (i) * e2 () (i);
#else
            vector_size_type i (0);
            DD (size, 4, r, (t += e1 () (i) * e2 () (i), ++ i));
#endif
            return t;
        }
        // Dense case
        template<class D, class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I1 it1, I2 it2) {
            result_type t = result_type (0);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            while (-- size >= 0)
                t += *it1 * *it2, ++ it1, ++ it2;
#else
            DD (size, 4, r, (t += *it1 * *it2, ++ it1, ++ it2));
#endif
            return t;
        }
        // Packed case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end) {
            result_type t = result_type (0);
            typedef typename I1::difference_type vector_difference_type;
            vector_difference_type it1_size (it1_end - it1);
            vector_difference_type it2_size (it2_end - it2);
            vector_difference_type diff (0);
            if (it1_size > 0 && it2_size > 0)
                diff = it2.index () - it1.index ();
            if (diff != 0) {
                vector_difference_type size = (std::min) (diff, it1_size);
                if (size > 0) {
                    it1 += size;
                    it1_size -= size;
                    diff -= size;
                }
                size = (std::min) (- diff, it2_size);
                if (size > 0) {
                    it2 += size;
                    it2_size -= size;
                    diff += size;
                }
            }
            vector_difference_type size ((std::min) (it1_size, it2_size));
            while (-- size >= 0)
                t += *it1 * *it2, ++ it1, ++ it2;
            return t;
        }
        // Sparse case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end, sparse_bidirectional_iterator_tag) {
            result_type t = result_type (0);
            if (it1 != it1_end && it2 != it2_end) {
                for (;;) {
                    if (it1.index () == it2.index ()) {
                        t += *it1 * *it2, ++ it1, ++ it2;
                        if (it1 == it1_end || it2 == it2_end)
                            break;
                    } else if (it1.index () < it2.index ()) {
                        increment (it1, it1_end, it2.index () - it1.index ());
                        if (it1 == it1_end)
                            break;
                    } else if (it1.index () > it2.index ()) {
                        increment (it2, it2_end, it1.index () - it2.index ());
                        if (it2 == it2_end)
                            break;
                    }
                }
            }
            return t;
        }
    };

    template<class V1, class V2, class TV>
    struct vector_covariance:
        public vector_scalar_binary_functor<V1, V2, TV> {
        typedef typename vector_scalar_binary_functor<V1, V2, TV>::value_type value_type;
        typedef typename vector_scalar_binary_functor<V1, V2, TV>::result_type result_type;

        template<class C1, class C2>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_container<C1> &c1, const vector_container<C2> &c2) {
#ifdef BOOST_UBLAS_USE_SIMD
            using namespace raw;
            typedef typename C1::size_type vector_size_type;
            vector_size_type size (BOOST_UBLAS_SAME (c1 ().size (), c2 ().size ()));
            const typename V1::value_type *data1 = data_const (c1 ());
            const typename V1::value_type *data2 = data_const (c2 ());
            vector_size_type s1 = stride (c1 ());
            vector_size_type s2 = stride (c2 ());
            result_type t = result_type (0);
            result_type mean1 = result_type (0);
            result_type mean2 = result_type (0);
            result_type del1, del2;
            if (s1 == 1 && s2 == 1) {
                for (vector_size_type i = 0; i < size; ++ i) {
                    del1 = (data1[i] - mean1) / (i + 1);
                    mean1 += del1;
                    del2 = (data2[i] - mean2) / (i + 1);
                    mean2 += del2;
                    t += i * del1 * del2 - t / (i + 1);
                }
            } else if (s2 == 1) {
                for (vector_size_type i = 0, i1 = 0; i < size; ++ i, i1 += s1) {
                    del1 = (data1[i1] - mean1) / (i + 1);
                    mean1 += del1;
                    del2 = (data2[i] - mean2) / (i + 1);
                    mean2 += del2;
                    t += i * del1 * del2 - t / (i + 1);
                }
            } else if (s1 == 1) {
                for (vector_size_type i = 0, i2 = 0; i < size; ++ i, i2 += s2) {
                    del1 = (data1[i] - mean1) / (i + 1);
                    mean1 += del1;
                    del2 = (data2[i2] - mean2) / (i + 1);
                    mean2 += del2;
                    t += i * del1 * del2 - t / (i + 1);
                }
            } else {
                for (vector_size_type i = 0, i1 = 0, i2 = 0; i < size; ++ i, i1 += s1, i2 += s2) {
                    del1 = (data1[i1] - mean1) / (i + 1);
                    mean1 += del1;
                    del2 = (data2[i2] - mean2) / (i + 1);
                    mean2 += del2;
                    t += i * del1 * del2 - t / (i + 1);
                }
            }
            return size / (size - 1) * t;
#else
            return apply (static_cast<const vector_expression<C1> > (c1), static_cast<const vector_expression<C2> > (c2));
#endif
        }

        template<class E1, class E2>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
            typedef typename E1::size_type vector_size_type;
            vector_size_type size (BOOST_UBLAS_SAME (e1 ().size (), e2 ().size ()));
            result_type t = result_type (0);
            result_type mean1 = result_type (0);
            result_type mean2 = result_type (0);
            result_type del1, del2;
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            for (vector_size_type i = 0; i < size; ++ i) {
                del1 = (e1 () (i) - mean1) / (i + 1);
                mean1 += del1;
                del2 = (e2 () (i) - mean2) / (i + 1);
                mean2 += del2;
                t += i * del1 * del2 - t / (i + 1);
            }
#else
            vector_size_type i (0);
            DD (size, 4, r, del1 = (e1 () (i) - mean1) / (i + 1),
                mean1 += del1,
                del2 = (e2 () (i) - mean2) / (i + 1),
                mean2 += del2,
                t += i * del1 * del2 - t / (i + 1), ++ i);
#endif
            // return size / (size - 1) * t;
            return t;
        }

        // Dense case
        template<class D, class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I1 it1, I2 it2) {
            result_type t = result_type (0);
            result_type i = result_type (0);
            result_type mean1 = result_type (0);
            result_type mean2 = result_type (0);
            result_type del1, del2;
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            while (++ i <= size) {
                del1 = (*it1 - mean1) / i;
                mean1 += del1;
                del2 = (*it2 - mean2) / i;
                mean2 += del2;
                t += (i - 1) * del1 * del2 - t / i;
                ++ it1;
                ++ it2;
            }
#else
            DD (size, 4, r, del1 = (*it1 - mean1) / (i + 1),
                mean1 += del1,
                del2 = (*it2 - mean2) / (i + 1),
                mean2 += del2,
                t += i * del1 * del2 - t / (i + 1),
                ++ it1, ++ it2, ++ i);
#endif
            return size / (size - 1) * t;
        }

        // Packed case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end) {
            result_type t = result_type (0);
            result_type i = result_type (0);
            result_type mean1 = result_type (0);
            result_type mean2 = result_type (0);
            result_type del1, del2;
            typedef typename I1::difference_type vector_difference_type;
            vector_difference_type it1_size (it1_end - it1);
            vector_difference_type it2_size (it2_end - it2);
            vector_difference_type diff (0);
            if (it1_size > 0 && it2_size > 0)
                diff = it2.index () - it1.index ();
            if (diff != 0) {
                vector_difference_type size = (std::min) (diff, it1_size);
                if (size > 0) {
                    it1 += size;
                    it1_size -= size;
                    diff -= size;
                }
                size = (std::min) (- diff, it2_size);
                if (size > 0) {
                    it2 += size;
                    it2_size -= size;
                    diff += size;
                }
            }
            vector_difference_type size ((std::min) (it1_size, it2_size));
            while (++ i <= size) {
                del1 = (*it1 - mean1) / i;
                mean1 += del1;
                del2 = (*it2 - mean2) / i;
                mean2 += del2;
                t += (i - 1) * del1 * del2 - t / i;
                ++ it1;
                ++ it2;
            }
            return size / (size - 1) * t;
        }

        // Sparse case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end, sparse_bidirectional_iterator_tag) {
            typedef typename I1::difference_type vector_difference_type;
            result_type t = result_type (0);
            vector_difference_type i (1);
            result_type mean1 = result_type (0);
            result_type mean2 = result_type (0);
            result_type del1, del2;
            if (it1 != it1_end && it2 != it2_end) {
                for (;;) {
                    if (it1.index () == it2.index ()) {
                        del1 = (*it1 - mean1) / i;
                        mean1 += del1;
                        del2 = (*it2 - mean2) / i;
                        mean2 += del2;
                        t += (i - 1) * del1 * del2 - t / i;
                        ++ i;
                        ++ it1;
                        ++ it2;
                        if (it1 == it1_end || it2 == it2_end)
                            break;
                    } else if (it1.index () < it2.index ()) {
                        increment (it1, it1_end, it2.index () - it1.index ());
                        if (it1 == it1_end)
                            break;
                    } else if (it1.index () > it2.index ()) {
                        increment (it2, it2_end, it1.index () - it2.index ());
                        if (it2 == it2_end)
                            break;
                    }
                }
            }
            return i / (i - 1) * t;
        }
    };

    // Matrix functors

    // Unary returning scalar
    template<class M>
    struct matrix_scalar_unary_functor {
        typedef typename M::value_type value_type;
        typedef typename M::value_type result_type;
    };

    template<class M>
    struct matrix_min:
        public matrix_scalar_unary_functor<M> {
        typedef typename matrix_scalar_unary_functor<M>::value_type value_type;
        typedef typename matrix_scalar_unary_functor<M>::result_type result_type;
        
        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) { 
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());
            result_type min_val = result_type ( e () (matrix_size_type (0), matrix_size_type (0)));
            for (matrix_size_type i = 0; i < size1; ++ i)
                for (matrix_size_type j = 0; j < size2; ++ j)
                    if ( type_traits<value_type>::type_abs(min_val) > type_traits<value_type>::type_abs(e () (i, j)) )
                        min_val = e () (i, j);
            return min_val;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type min_val = *it;
            while (-- size >= 0) {
                if ( type_traits<value_type>::type_abs(min_val) > type_traits<value_type>::type_abs(*it) )
                    min_val = *it;
                ++ it;
            }
            return min_val;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type min_val = *it;
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size (0);
            while (it != it_end) {
                if ( type_traits<value_type>::type_abs(min_val) > type_traits<value_type>::type_abs(*it) )
                    min_val = *it;
                ++ it;
                ++ size;
            }
            return min_val;
        }
    };

    template<class M>
    struct matrix_max:
        public matrix_scalar_unary_functor<M> {
        typedef typename matrix_scalar_unary_functor<M>::value_type value_type;
        typedef typename matrix_scalar_unary_functor<M>::result_type result_type;
        
        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) { 
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());
            result_type max_val = result_type ( e () (matrix_size_type (0), matrix_size_type (0)));
            for (matrix_size_type i = 0; i < size1; ++ i)
                for (matrix_size_type j = 0; j < size2; ++ j)
                    if ( type_traits<value_type>::type_abs(max_val) < type_traits<value_type>::type_abs(e () (i, j)) )
                        max_val = e () (i, j);
            return max_val;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type max_val = *it;
            while (-- size >= 0) {
                if ( type_traits<value_type>::type_abs(max_val) < type_traits<value_type>::type_abs(*it) )
                    max_val = *it;
                ++ it;
            }
            return max_val;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type max_val = *it;
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size (0);
            while (it != it_end) {
                if ( type_traits<value_type>::type_abs(max_val) < type_traits<value_type>::type_abs(*it) )
                    max_val = *it;
                ++ it;
                ++ size;
            }
            return max_val;
        }
    };

    template<class M>
    struct matrix_sum:
        public matrix_scalar_unary_functor<M> {
        typedef typename matrix_scalar_unary_functor<M>::value_type value_type;
        typedef typename matrix_scalar_unary_functor<M>::result_type result_type;
        // typedef double result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) { 
            result_type t = result_type (0);
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());
            for (matrix_size_type i = 0; i < size1; ++ i) {
                for (matrix_size_type j = 0; j < size2; ++ j) {
                    t += e () (i, j);
                }
            }
            return t;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type t = result_type (0);
            while (-- size >= 0) {
                t += *it;
                ++ it;
            }
            return t; 
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type t = result_type (0);
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size (0);
            while (it != it_end) {
                t += *it;
                ++ it;
                ++ size;
            }
            return t;
        }
    };

    template<class M>
    struct matrix_mean:
        public matrix_scalar_unary_functor<M> {
        typedef typename matrix_scalar_unary_functor<M>::value_type value_type;
        typedef typename matrix_scalar_unary_functor<M>::result_type result_type;
        // typedef double result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) { 
            result_type t = result_type (0);
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());
            for (matrix_size_type i = 0; i < size1; ++ i) {
                for (matrix_size_type j = 0; j < size2; ++ j) {
                    t += e () (i, j);
                }
            }
            return t / (size1 * size2);
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type t = result_type (0);
            while (-- size >= 0) {
                t += *it;
                ++ it;
            }
            return t / size; 
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type t = result_type (0);
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size (0);
            while (it != it_end) {
                t += *it;
                ++ it;
                ++ size;
            }
            return t / size;
        }
    };

    template<class M>
    struct matrix_variance:
        public matrix_scalar_unary_functor<M> {
        typedef typename matrix_scalar_unary_functor<M>::value_type value_type;
        typedef typename matrix_scalar_unary_functor<M>::result_type result_type;
        // typedef double result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) { 
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());
            matrix_size_type num_elements = size1 * size2;
            result_type sumsq = result_type (0);
            result_type sum = result_type (0);
            for (matrix_size_type i = 0; i < size1; ++ i) {
                for (matrix_size_type j = 0; j < size2; ++ j) {
                    sumsq += e () (i, j) * e () (i, j);
                    sum += e () (i, j);
                }
            }
            return (sumsq - (sum * sum) / num_elements) / (num_elements);
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            result_type sumsq = result_type (0);
            result_type sum = result_type (0);
            D i (0);
            while (++i <= size) {
                sumsq += *it * *it;
                sum += *it;
                ++ it;
            }
            return (sumsq - (sum * sum) / size) / size;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type sumsq = result_type (0);
            result_type sum = result_type (0);
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size (0);
            while (it != it_end) {
                sumsq += *it * *it;
                sum += *it;
                ++ it;
                ++ size;
            }
            return (sumsq - (sum * sum) / size) / size;
        }
    };

    template<class V>
    struct matrix_mode:
        public matrix_scalar_unary_functor<V> {
        typedef typename matrix_scalar_unary_functor<V>::value_type value_type;
        typedef typename matrix_scalar_unary_functor<V>::result_type result_type;
        
        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) {
            typedef typename E::size_type matrix_size_type;
            boost::unordered_map<result_type, matrix_size_type> count_map;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());
            typename boost::unordered_map<result_type, matrix_size_type>::iterator p;
            // auto p = count_map.iterator;
            for (matrix_size_type i = 0; i < size1; ++ i)
                for (matrix_size_type j = 0; j < size2; ++ j) {
                    p = count_map.find (e () (i, j));
                    if (p != count_map.end ())
                        p->second += matrix_size_type(1);
                    else
                        count_map.emplace (e () (i, j), matrix_size_type (1));
                }
            
            result_type mode = result_type (0);
            matrix_size_type mode_val = matrix_size_type(0);
            p = count_map.begin();
            while (p != count_map.end()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs(p->first) < type_traits<value_type>::type_abs(mode))
                    mode = p->first;
                ++ p;
            }
            
            return mode;
        }
        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) { 
            boost::unordered_map<result_type, D> count_map;
            
            typename boost::unordered_map<result_type, D>::iterator p;
            // auto p = count_map.iterator;
            while (-- size >= 0) {
                p = count_map.find (*it);
                if ( p != count_map.end ())
                    p->second += D(1);
                else
                    count_map.emplace (*it, D (1));
                ++ it;
            }
            
            result_type mode = result_type (0);
            D mode_val = D(0);
            p = count_map.begin();
            while (p != count_map.end()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs(p->first) < type_traits<value_type>::type_abs(mode))
                    mode = p->first;
                ++ p;
            }
            
            return mode;
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            typedef typename I::difference_type matrix_difference_type;

            boost::unordered_map<result_type, matrix_difference_type> count_map;
            
            typename boost::unordered_map<result_type, matrix_difference_type>::iterator p;
            // auto p = count_map.iterator;

            while (it != it_end) {
                p = count_map.find (*it);
                if ( p != count_map.end ())
                    p->second += matrix_difference_type(1);
                else
                    count_map.emplace (*it, matrix_difference_type (1));
                ++ it;
            }
            
            result_type mode = result_type (0);
            matrix_difference_type mode_val = matrix_difference_type(0);
            p = count_map.begin();
            while (p != count_map.end()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs(p->first) < type_traits<value_type>::type_abs(mode))
                    mode = p->first;
                ++ p;
            }
            
            return mode;
        }
    };

    template<class V>
    struct matrix_median:
        public matrix_scalar_unary_functor<V> {
        typedef typename vector_scalar_unary_functor<V>::value_type value_type;
        typedef typename vector_scalar_unary_functor<V>::result_type result_type;
        
        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) {
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());
            matrix_size_type num_elements = size1 * size2;
            vector<value_type> v (num_elements);
            matrix_size_type k (0);
            for (matrix_size_type i = 0; i < size1; ++ i)
                for (matrix_size_type j = 0; j < size2; ++ j)
                    v (k ++) = e () (i, j);
            boost::sort (v, compareElement);
            if (num_elements % 2)
                return v (num_elements / 2);
            else
                return (v (num_elements / 2) + v ((num_elements / 2) - 1)) / 2;
        }

        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            vector<value_type> v (size);
            D size_temp = size;
            D i (0);
            while (-- size_temp >= 0) {
                v (i) = *it;
                ++ it;
                ++  i;
            }
            boost::sort (v, compareElement);
            if (size % 2)
                return v (size / 2);
            else
                return (v (size / 2) + v ((size / 2) - 1)) / 2;
        }

        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size = matrix_difference_type (0);
            I it1 = it;
            while (it1 != it_end) {
                ++ it1;
                ++ size;
            }
            vector<value_type> v (size);
            matrix_difference_type i (0);
            while (it != it_end) {
                v (i) = *it;
                ++ it;
                ++ i;
            }
            
            boost::sort (v, compareElement);
            if (size % 2)
                return v (size / 2);
            else
                return (v (size / 2) + v ((size / 2) - 1)) / 2;
        }

        static BOOST_UBLAS_INLINE
        bool compareElement (const value_type &A, const value_type &B) {
            return (type_traits<value_type>::type_abs (A) < type_traits<value_type>::type_abs (B));
        }
    };

    // Unary returning vector of value_type TV

    template<class M, class TV>
    struct matrix_vector_unary_functor {
        typedef typename M::value_type value_type;
        typedef TV result_type;
    };

    template<class M, class TV>
    struct matrix_min_axis: 
        public matrix_vector_unary_functor<M, TV> {
        typedef typename matrix_vector_unary_functor<M, TV>::value_type value_type;
        typedef typename matrix_vector_unary_functor<M, TV>::result_type result_type;

        // template<class E>
        // static BOOST_UBLAS_INLINE
        // result_type apply (const matrix_expression<E> &e, typename E::size_type axis) { 
        //     typedef typename E::size_type matrix_size_type;
        //     matrix_size_type size1 (e ().size1 ());
        //     matrix_size_type size2 (e ().size2 ());

        //     matrix_size_type result_vector_size = matrix_size_type(0);
        //     if (axis == 0)
        //         result_vector_size = size2;
        //     else if (axis == 1)
        //         result_vector_size = size1;

        //     result_type t = result_type (result_vector_size);
        //     for (matrix_size_type i = 0; i < result_vector_size; ++ i) {
        //         t (i) = apply (e, axis, i);
        //     }
        //     return t;
        // }

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e, typename E::size_type axis, typename E::size_type index) { 
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());

            result_type vmin = result_type(0);
            
            if (axis == 0){
                vmin = result_type(e () (0, index));
                for (matrix_size_type i = 0; i < size1; ++ i)
                    if ( type_traits<value_type>::type_abs(vmin) > type_traits<value_type>::type_abs(e () (i, index)) )
                        vmin = e () (i, index);
                return vmin;
            }
            else if (axis == 1){
                vmin = result_type(e () (index, 0));
                for (matrix_size_type i = 0; i < size2; ++ i)
                    if ( type_traits<value_type>::type_abs(vmin) > type_traits<value_type>::type_abs(e () (index, i)) )
                        vmin = e () (index, i);
                return vmin;
            }
        }

        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            result_type vmin = *it;
            while (-- size >= 0) {
                if ( type_traits<value_type>::type_abs(vmin) > type_traits<value_type>::type_abs(*it) )
                    vmin = *it;
                ++ it;
            }
            return vmin;
        }

        // Packed case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type vmin = *it;
            while (it != it_end){
                if ( type_traits<value_type>::type_abs(vmin) > type_traits<value_type>::type_abs(*it) )
                    vmin = *it;
                ++ it;
            }
            return vmin;
        }

        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end, sparse_bidirectional_iterator_tag) {
            result_type vmin = *it;
            while (it != it_end){
                if ( type_traits<value_type>::type_abs(vmin) > type_traits<value_type>::type_abs(*it) )
                    vmin = *it;
                ++ it;
            }
            return vmin;
        //     result_type t = result_type (0);
        //     if (it1 != it1_end && it2 != it2_end) {
        //         size_type it1_index = it1.index2 (), it2_index = it2.index ();
        //         for (;;) {
        //             difference_type compare = it1_index - it2_index;
        //             if (compare == 0) {
        //                 t += *it1 * *it2, ++ it1, ++ it2;
        //                 if (it1 != it1_end && it2 != it2_end) {
        //                     it1_index = it1.index2 ();
        //                     it2_index = it2.index ();
        //                 } else
        //                     break;
        //             } else if (compare < 0) {
        //                 increment (it1, it1_end, - compare);
        //                 if (it1 != it1_end)
        //                     it1_index = it1.index2 ();
        //                 else
        //                     break;
        //             } else if (compare > 0) {
        //                 increment (it2, it2_end, compare);
        //                 if (it2 != it2_end)
        //                     it2_index = it2.index ();
        //                 else
        //                     break;
        //             }
        //         }
        //     }
        //     return t;
        }
    };

    template<class M, class TV>
    struct matrix_max_axis: 
        public matrix_vector_unary_functor<M, TV> {
        typedef typename matrix_vector_unary_functor<M, TV>::value_type value_type;
        typedef typename matrix_vector_unary_functor<M, TV>::result_type result_type;

        // template<class E>
        // static BOOST_UBLAS_INLINE
        // result_type apply (const matrix_expression<E> &e, typename E::size_type axis) { 
        //     typedef typename E::size_type matrix_size_type;
        //     matrix_size_type size1 (e ().size1 ());
        //     matrix_size_type size2 (e ().size2 ());

        //     matrix_size_type result_vector_size = matrix_size_type(0);
        //     if (axis == 0)
        //         result_vector_size = size2;
        //     else if (axis == 1)
        //         result_vector_size = size1;

        //     result_type t = result_type (result_vector_size);
        //     for (matrix_size_type i = 0; i < result_vector_size; ++ i) {
        //         t (i) = apply (e, axis, i);
        //     }
        //     return t;
        // }

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e, typename E::size_type axis, typename E::size_type index) { 
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());

            result_type vmax = result_type(0);
            
            if (axis == 0){
                vmax = result_type(e () (0, index));
                for (matrix_size_type i = 0; i < size1; ++ i)
                    if ( type_traits<value_type>::type_abs(vmax) < type_traits<value_type>::type_abs(e () (i, index)) )
                        vmax = e () (i, index);
                return vmax;
            }
            else if (axis == 1){
                vmax = result_type(e () (index, 0));
                for (matrix_size_type i = 0; i < size2; ++ i)
                    if ( type_traits<value_type>::type_abs(vmax) < type_traits<value_type>::type_abs(e () (index, i)) )
                        vmax = e () (index, i);
                return vmax;
            }
        }

        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            result_type vmax = *it;
            while (-- size >= 0) {
                if ( type_traits<value_type>::type_abs(vmax) < type_traits<value_type>::type_abs(*it) )
                    vmax = *it;
                ++ it;
            }
            return vmax;
        }

        // Packed case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type vmax = *it;
            while (it != it_end){
                if ( type_traits<value_type>::type_abs(vmax) < type_traits<value_type>::type_abs(*it) )
                    vmax = *it;
                ++ it;
            }
            return vmax;
        }

        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end, sparse_bidirectional_iterator_tag) {
            result_type vmax = *it;
            while (it != it_end){
                if ( type_traits<value_type>::type_abs(vmax) < type_traits<value_type>::type_abs(*it) )
                    vmax = *it;
                ++ it;
            }
            return vmax;
        //     result_type t = result_type (0);
        //     if (it1 != it1_end && it2 != it2_end) {
        //         size_type it1_index = it1.index2 (), it2_index = it2.index ();
        //         for (;;) {
        //             difference_type compare = it1_index - it2_index;
        //             if (compare == 0) {
        //                 t += *it1 * *it2, ++ it1, ++ it2;
        //                 if (it1 != it1_end && it2 != it2_end) {
        //                     it1_index = it1.index2 ();
        //                     it2_index = it2.index ();
        //                 } else
        //                     break;
        //             } else if (compare < 0) {
        //                 increment (it1, it1_end, - compare);
        //                 if (it1 != it1_end)
        //                     it1_index = it1.index2 ();
        //                 else
        //                     break;
        //             } else if (compare > 0) {
        //                 increment (it2, it2_end, compare);
        //                 if (it2 != it2_end)
        //                     it2_index = it2.index ();
        //                 else
        //                     break;
        //             }
        //         }
        //     }
        //     return t;
        }
    };

    template<class M, class TV>
    struct matrix_sum_axis: 
        public matrix_vector_unary_functor<M, TV> {
        typedef typename matrix_vector_unary_functor<M, TV>::value_type value_type;
        typedef typename matrix_vector_unary_functor<M, TV>::result_type result_type;

        // template<class E>
        // static BOOST_UBLAS_INLINE
        // result_type apply (const matrix_expression<E> &e, typename E::size_type axis) {
        //     typedef typename E::size_type matrix_size_type;
        //     matrix_size_type size1 (e ().size1 ());
        //     matrix_size_type size2 (e ().size2 ());

        //     matrix_size_type result_vector_size = matrix_size_type(0);
        //     if (axis == 0)
        //         result_vector_size = size2;
        //     else if (axis == 1)
        //         result_vector_size = size1;

        //     result_type t = result_type (result_vector_size);
        //     for (matrix_size_type i = 0; i < result_vector_size; ++ i) {
        //         t (i) = apply (e, axis, i);
        //     }
        //     return t;
        // }

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e, typename E::size_type axis, typename E::size_type index) { 
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());

            result_type vsum = result_type(0);
            
            if (axis == 0){
                for (matrix_size_type i = 0; i < size1; ++ i) {
                    vsum += e () (i, index);
                }
                return vsum;
            }
            else if (axis == 1){
                for (matrix_size_type i = 0; i < size2; ++ i) {
                    vsum += e () (index, i);
                }
                return vsum;
            }

        }

        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            result_type t = result_type (0);
            while (-- size >= 0) {
                t += *it;
                ++ it;
            }
            return t;
        }

        // Packed case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            result_type t = result_type (0);
            while (it != it_end){
                t += *it;
                ++ it;
            }
            return t;
        }

        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end, sparse_bidirectional_iterator_tag) {
            result_type t = result_type (0);
            while (it != it_end){
                t += *it;
                ++ it;
            }
            return t;
        //     result_type t = result_type (0);
        //     if (it1 != it1_end && it2 != it2_end) {
        //         size_type it1_index = it1.index2 (), it2_index = it2.index ();
        //         for (;;) {
        //             difference_type compare = it1_index - it2_index;
        //             if (compare == 0) {
        //                 t += *it1 * *it2, ++ it1, ++ it2;
        //                 if (it1 != it1_end && it2 != it2_end) {
        //                     it1_index = it1.index2 ();
        //                     it2_index = it2.index ();
        //                 } else
        //                     break;
        //             } else if (compare < 0) {
        //                 increment (it1, it1_end, - compare);
        //                 if (it1 != it1_end)
        //                     it1_index = it1.index2 ();
        //                 else
        //                     break;
        //             } else if (compare > 0) {
        //                 increment (it2, it2_end, compare);
        //                 if (it2 != it2_end)
        //                     it2_index = it2.index ();
        //                 else
        //                     break;
        //             }
        //         }
        //     }
        //     return t;
        }

    };

    template<class M, class TV>
    struct matrix_mean_axis: 
        public matrix_vector_unary_functor<M, TV> {
        typedef typename matrix_vector_unary_functor<M, TV>::value_type value_type;
        typedef typename matrix_vector_unary_functor<M, TV>::result_type result_type;

        // template<class E>
        // static BOOST_UBLAS_INLINE
        // result_type apply (const matrix_expression<E> &e, typename E::size_type axis) {
        //     typedef typename E::size_type matrix_size_type;
        //     matrix_size_type size1 (e ().size1 ());
        //     matrix_size_type size2 (e ().size2 ());

        //     matrix_size_type result_vector_size = matrix_size_type(0);
        //     if (axis == 0)
        //         result_vector_size = size2;
        //     else if (axis == 1)
        //         result_vector_size = size1;

        //     result_type t = result_type (result_vector_size);
        //     for (matrix_size_type i = 0; i < result_vector_size; ++ i) {
        //         t (i) = apply (e, axis, i);
        //     }
        //     return t;
        // }

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e, typename E::size_type axis, typename E::size_type index) { 
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());

            result_type vsum = result_type(0);
            
            if (axis == 0){
                for (matrix_size_type i = 0; i < size1; ++ i) {
                    vsum += e () (i, index);
                }
                return vsum / size1;
            }
            else if (axis == 1){
                for (matrix_size_type i = 0; i < size2; ++ i) {
                    vsum += e () (index, i);
                }
                return vsum / size2;
            }

        }

        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            result_type t = result_type (0);
            while (-- size >= 0) {
                t += *it;
                ++ it;
            }
            return t / size;
        }

        // Packed case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size (0);
            result_type t = result_type (0);
            while (it != it_end){
                t += *it;
                ++ it;
                ++ size;
            }
            return t / size;
        }

        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end, sparse_bidirectional_iterator_tag) {
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size (0);
            result_type t = result_type (0);
            while (it != it_end){
                t += *it;
                ++ it;
                ++ size;
            }
            return t / size;
        //     result_type t = result_type (0);
        //     if (it1 != it1_end && it2 != it2_end) {
        //         size_type it1_index = it1.index2 (), it2_index = it2.index ();
        //         for (;;) {
        //             difference_type compare = it1_index - it2_index;
        //             if (compare == 0) {
        //                 t += *it1 * *it2, ++ it1, ++ it2;
        //                 if (it1 != it1_end && it2 != it2_end) {
        //                     it1_index = it1.index2 ();
        //                     it2_index = it2.index ();
        //                 } else
        //                     break;
        //             } else if (compare < 0) {
        //                 increment (it1, it1_end, - compare);
        //                 if (it1 != it1_end)
        //                     it1_index = it1.index2 ();
        //                 else
        //                     break;
        //             } else if (compare > 0) {
        //                 increment (it2, it2_end, compare);
        //                 if (it2 != it2_end)
        //                     it2_index = it2.index ();
        //                 else
        //                     break;
        //             }
        //         }
        //     }
        //     return t;
        }
    };

    template<class M, class TV>
    struct matrix_variance_axis: 
        public matrix_vector_unary_functor<M, TV> {
        typedef typename matrix_vector_unary_functor<M, TV>::value_type value_type;
        typedef typename matrix_vector_unary_functor<M, TV>::result_type result_type;

        // template<class E>
        // static BOOST_UBLAS_INLINE
        // result_type apply (const matrix_expression<E> &e, typename E::size_type axis) {
        //     typedef typename E::size_type matrix_size_type;
        //     matrix_size_type size1 (e ().size1 ());
        //     matrix_size_type size2 (e ().size2 ());

        //     matrix_size_type result_vector_size = matrix_size_type(0);
        //     if (axis == 0)
        //         result_vector_size = size2;
        //     else if (axis == 1)
        //         result_vector_size = size1;

        //     result_type t = result_type (result_vector_size);
        //     for (matrix_size_type i = 0; i < result_vector_size; ++ i) {
        //         t (i) = apply (e, axis, i);
        //     }
        //     return t;
        // }

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e, typename E::size_type axis, typename E::size_type index) { 
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());

            result_type sumsq = result_type (0);
            result_type sum = result_type (0);

            if (axis == 0){
                for (matrix_size_type i = 0; i < size1; ++ i) {
                    sumsq += e () (i, index) * e () (i, index);
                    sum += e () (i, index);
                }
                return (sumsq - (sum * sum) / size1) / size1;
            }
            else if (axis == 1){
                for (matrix_size_type i = 0; i < size2; ++ i) {
                    sumsq += e () (index, i) * e () (index, i);
                    sum += e () (index, i);
                }
                return (sumsq - (sum * sum) / size2) / size2;
            }

        }

        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            result_type sumsq = result_type (0);
            result_type sum = result_type (0);
            while (-- size >= 0) {
                sumsq += *it * *it;
                sum += *it;
                ++ it;
            }
            return (sumsq - (sum * sum) / size) / size;
        }

        // Packed case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size (0);
            result_type sumsq = result_type (0);
            result_type sum = result_type (0);
            while (it != it_end){
                sumsq += *it * *it;
                sum += *it;
                ++ it;
                ++ size;
            }
            return (sumsq - (sum * sum) / size) / size;
        }

        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end, sparse_bidirectional_iterator_tag) {
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size (0);
            result_type sumsq = result_type (0);
            result_type sum = result_type (0);
            while (it != it_end){
                sumsq += *it * *it;
                sum += *it;
                ++ it;
                ++ size;
            }
            return (sumsq - (sum * sum) / size) / size;
        //     result_type t = result_type (0);
        //     if (it1 != it1_end && it2 != it2_end) {
        //         size_type it1_index = it1.index2 (), it2_index = it2.index ();
        //         for (;;) {
        //             difference_type compare = it1_index - it2_index;
        //             if (compare == 0) {
        //                 t += *it1 * *it2, ++ it1, ++ it2;
        //                 if (it1 != it1_end && it2 != it2_end) {
        //                     it1_index = it1.index2 ();
        //                     it2_index = it2.index ();
        //                 } else
        //                     break;
        //             } else if (compare < 0) {
        //                 increment (it1, it1_end, - compare);
        //                 if (it1 != it1_end)
        //                     it1_index = it1.index2 ();
        //                 else
        //                     break;
        //             } else if (compare > 0) {
        //                 increment (it2, it2_end, compare);
        //                 if (it2 != it2_end)
        //                     it2_index = it2.index ();
        //                 else
        //                     break;
        //             }
        //         }
        //     }
        //     return t;
        }
    };

    template<class M, class TV>
    struct matrix_mode_axis: 
        public matrix_vector_unary_functor<M, TV> {
        typedef typename matrix_vector_unary_functor<M, TV>::value_type value_type;
        typedef typename matrix_vector_unary_functor<M, TV>::result_type result_type;

        // template<class E>
        // static BOOST_UBLAS_INLINE
        // result_type apply (const matrix_expression<E> &e, typename E::size_type axis) {
        //     typedef typename E::size_type matrix_size_type;
        //     matrix_size_type size1 (e ().size1 ());
        //     matrix_size_type size2 (e ().size2 ());

        //     matrix_size_type result_vector_size = matrix_size_type(0);
        //     if (axis == 0)
        //         result_vector_size = size2;
        //     else if (axis == 1)
        //         result_vector_size = size1;

        //     result_type t = result_type (result_vector_size);
        //     for (matrix_size_type i = 0; i < result_vector_size; ++ i) {
        //         t (i) = apply (e, axis, i);
        //     }
        //     return t;
        // }

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e, typename E::size_type axis, typename E::size_type index) {
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());

            boost::unordered_map<result_type, matrix_size_type> count_map;
            typename boost::unordered_map<result_type, matrix_size_type>::iterator p;
            
            if (axis == 0){
                for (matrix_size_type i = 0; i < size1; ++ i) {
                    p = count_map.find (e () (i, index));
                    if (p != count_map.end ())
                        p->second += matrix_size_type (1);
                    else
                        count_map.emplace (e () (i, index), matrix_size_type (1));
                }

            }
            else if (axis == 1){
                for (matrix_size_type i = 0; i < size2; ++ i) {
                    p = count_map.find (e () (index, i));
                    if (p != count_map.end ())
                        p->second += matrix_size_type (1);
                    else
                        count_map.emplace (e () (index, i), matrix_size_type (1));
                }
            }

            result_type mode = result_type (0);
            matrix_size_type mode_val = matrix_size_type (0);
            p = count_map.begin();
            while (p != count_map.end ()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs (p->first) < type_traits<value_type>::type_abs (mode))
                    mode = p->first;
                ++ p;
            }
            
            return mode;
        }

        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            boost::unordered_map<result_type, D> count_map;
            typename boost::unordered_map<result_type, D>::iterator p;
            
            while (-- size >= 0) {
                p = count_map.find (*it);
                if (p != count_map.end ())
                    p->second += D (1);
                else
                    count_map.emplace (*it, D (1));
                ++ it;
            }

            result_type mode = result_type (0);
            D mode_val = D (0);
            p = count_map.begin();
            while (p != count_map.end ()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs (p->first) < type_traits<value_type>::type_abs (mode))
                    mode = p->first;
                ++ p;
            }
            
            return mode;
        }

        // Packed case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            typedef typename I::difference_type matrix_difference_type;
            
            boost::unordered_map<result_type, matrix_difference_type> count_map;
            typename boost::unordered_map<result_type, matrix_difference_type>::iterator p;
            
            while (it != it_end){
                p = count_map.find (*it);
                if (p != count_map.end ())
                    p->second += matrix_difference_type (1);
                else
                    count_map.emplace (*it, matrix_difference_type (1));
                ++ it;
            }

            result_type mode = result_type (0);
            matrix_difference_type mode_val = matrix_difference_type (0);
            p = count_map.begin();
            while (p != count_map.end ()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs (p->first) < type_traits<value_type>::type_abs (mode))
                    mode = p->first;
                ++ p;
            }

            return mode;
        }

        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end, sparse_bidirectional_iterator_tag) {
            typedef typename I::difference_type matrix_difference_type;
            
            boost::unordered_map<result_type, matrix_difference_type> count_map;
            typename boost::unordered_map<result_type, matrix_difference_type>::iterator p;
            
            while (it != it_end){
                p = count_map.find (*it);
                if (p != count_map.end ())
                    p->second += matrix_difference_type (1);
                else
                    count_map.emplace (*it, matrix_difference_type (1));
                ++ it;
            }

            result_type mode = result_type (0);
            matrix_difference_type mode_val = matrix_difference_type (0);
            p = count_map.begin();
            while (p != count_map.end ()) {
                if (p->second > mode_val) {
                    mode = p->first;
                    mode_val = p->second;
                }
                else if (p->second == mode_val && type_traits<value_type>::type_abs (p->first) < type_traits<value_type>::type_abs (mode))
                    mode = p->first;
                ++ p;
            }

            return mode;
        //     result_type t = result_type (0);
        //     if (it1 != it1_end && it2 != it2_end) {
        //         size_type it1_index = it1.index2 (), it2_index = it2.index ();
        //         for (;;) {
        //             difference_type compare = it1_index - it2_index;
        //             if (compare == 0) {
        //                 t += *it1 * *it2, ++ it1, ++ it2;
        //                 if (it1 != it1_end && it2 != it2_end) {
        //                     it1_index = it1.index2 ();
        //                     it2_index = it2.index ();
        //                 } else
        //                     break;
        //             } else if (compare < 0) {
        //                 increment (it1, it1_end, - compare);
        //                 if (it1 != it1_end)
        //                     it1_index = it1.index2 ();
        //                 else
        //                     break;
        //             } else if (compare > 0) {
        //                 increment (it2, it2_end, compare);
        //                 if (it2 != it2_end)
        //                     it2_index = it2.index ();
        //                 else
        //                     break;
        //             }
        //         }
        //     }
        //     return t;
        }
    };

    template<class M, class TV>
    struct matrix_median_axis: 
        public matrix_vector_unary_functor<M, TV> {
        typedef typename matrix_vector_unary_functor<M, TV>::value_type value_type;
        typedef typename matrix_vector_unary_functor<M, TV>::result_type result_type;

        // template<class E>
        // static BOOST_UBLAS_INLINE
        // result_type apply (const matrix_expression<E> &e, typename E::size_type axis) {
        //     typedef typename E::size_type matrix_size_type;
        //     matrix_size_type size1 (e ().size1 ());
        //     matrix_size_type size2 (e ().size2 ());

        //     matrix_size_type result_vector_size = matrix_size_type(0);
        //     if (axis == 0)
        //         result_vector_size = size2;
        //     else if (axis == 1)
        //         result_vector_size = size1;

        //     result_type t = result_type (result_vector_size);
        //     for (matrix_size_type i = 0; i < result_vector_size; ++ i) {
        //         t (i) = apply (e, axis, i);
        //     }
        //     return t;
        // }

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e, typename E::size_type axis, typename E::size_type index) {
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            matrix_size_type size2 (e ().size2 ());

            if (axis == 0){
                vector<value_type> v (size1);
                for (matrix_size_type i = 0; i < size1; ++ i)
                    v (i) = e () (i, index);
                boost::sort (v, compareElement);
                if (size1 % 2)
                    return v (size1 / 2);
                else
                    return (v (size1 / 2) + v ((size1 / 2) - 1)) / 2;

            }
            else if (axis == 1){
                vector<value_type> v (size2);
                for (matrix_size_type i = 0; i < size2; ++ i)
                    v (i) = e () (index, i);
                boost::sort (v, compareElement);
                if (size2 % 2)
                    return v (size2 / 2);
                else
                    return (v (size2 / 2) + v ((size2 / 2) - 1)) / 2;
            }
        }

        // Dense case
        template<class D, class I>
        static BOOST_UBLAS_INLINE
        result_type apply (D size, I it) {
            vector<value_type> v (size);
            D size_temp = size;
            D i (0);
            while (-- size_temp >= 0) {
                v (i) = *it;
                ++ it;
                ++  i;
            }
            boost::sort (v, compareElement);
            if (size % 2)
                return v (size / 2);
            else
                return (v (size / 2) + v ((size / 2) - 1)) / 2;
        }

        // Packed case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size = matrix_difference_type (0);
            I it1 = it;
            while (it1 != it_end) {
                ++ it1;
                ++ size;
            }
            vector<value_type> v (size);
            matrix_difference_type i (0);
            while (it != it_end) {
                v (i) = *it;
                ++ it;
                ++ i;
            }
            
            boost::sort (v, compareElement);
            if (size % 2)
                return v (size / 2);
            else
                return (v (size / 2) + v ((size / 2) - 1)) / 2;
        }

        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end, sparse_bidirectional_iterator_tag) {
            typedef typename I::difference_type matrix_difference_type;
            matrix_difference_type size = matrix_difference_type (0);
            I it1 = it;
            while (it1 != it_end) {
                ++ it1;
                ++ size;
            }
            vector<value_type> v (size);
            matrix_difference_type i (0);
            while (it != it_end) {
                v (i) = *it;
                ++ it;
                ++ i;
            }
            
            boost::sort (v, compareElement);
            if (size % 2)
                return v (size / 2);
            else
                return (v (size / 2) + v ((size / 2) - 1)) / 2;
        //     result_type t = result_type (0);
        //     if (it1 != it1_end && it2 != it2_end) {
        //         size_type it1_index = it1.index2 (), it2_index = it2.index ();
        //         for (;;) {
        //             difference_type compare = it1_index - it2_index;
        //             if (compare == 0) {
        //                 t += *it1 * *it2, ++ it1, ++ it2;
        //                 if (it1 != it1_end && it2 != it2_end) {
        //                     it1_index = it1.index2 ();
        //                     it2_index = it2.index ();
        //                 } else
        //                     break;
        //             } else if (compare < 0) {
        //                 increment (it1, it1_end, - compare);
        //                 if (it1 != it1_end)
        //                     it1_index = it1.index2 ();
        //                 else
        //                     break;
        //             } else if (compare > 0) {
        //                 increment (it2, it2_end, compare);
        //                 if (it2 != it2_end)
        //                     it2_index = it2.index ();
        //                 else
        //                     break;
        //             }
        //         }
        //     }
        //     return t;
        }

        static BOOST_UBLAS_INLINE
        bool compareElement (const value_type &A, const value_type &B) {
            return (type_traits<value_type>::type_abs (A) < type_traits<value_type>::type_abs (B));
        }
    };

    // Binary returning vector
    template<class M1, class M2, class TV>
    struct matrix_vector_binary_functor {
        typedef typename M1::size_type size_type;
        typedef typename M1::difference_type difference_type;
        typedef TV value_type;
        typedef TV result_type;
    };

    template<class M1, class M2, class TV>
    struct matrix_vector_prod1:
        public matrix_vector_binary_functor<M1, M2, TV> {
        typedef typename matrix_vector_binary_functor<M1, M2, TV>::size_type size_type;
        typedef typename matrix_vector_binary_functor<M1, M2, TV>::difference_type difference_type;
        typedef typename matrix_vector_binary_functor<M1, M2, TV>::value_type value_type;
        typedef typename matrix_vector_binary_functor<M1, M2, TV>::result_type result_type;

        template<class C1, class C2>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_container<C1> &c1,
                           const vector_container<C2> &c2,
                           size_type i) {
#ifdef BOOST_UBLAS_USE_SIMD
            using namespace raw;
            size_type size = BOOST_UBLAS_SAME (c1 ().size2 (), c2 ().size ());
            const typename M1::value_type *data1 = data_const (c1 ()) + i * stride1 (c1 ());
            const typename M2::value_type *data2 = data_const (c2 ());
            size_type s1 = stride2 (c1 ());
            size_type s2 = stride (c2 ());
            result_type t = result_type (0);
            if (s1 == 1 && s2 == 1) {
                for (size_type j = 0; j < size; ++ j)
                    t += data1 [j] * data2 [j];
            } else if (s2 == 1) {
                for (size_type j = 0, j1 = 0; j < size; ++ j, j1 += s1)
                    t += data1 [j1] * data2 [j];
            } else if (s1 == 1) {
                for (size_type j = 0, j2 = 0; j < size; ++ j, j2 += s2)
                    t += data1 [j] * data2 [j2];
            } else {
                for (size_type j = 0, j1 = 0, j2 = 0; j < size; ++ j, j1 += s1, j2 += s2)
                    t += data1 [j1] * data2 [j2];
            }
            return t;
#elif defined(BOOST_UBLAS_HAVE_BINDINGS)
            return boost::numeric::bindings::atlas::dot (c1 ().row (i), c2 ());
#else
            return apply (static_cast<const matrix_expression<C1> > (c1), static_cast<const vector_expression<C2> > (c2, i));
#endif
        }
        template<class E1, class E2>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E1> &e1,
                           const vector_expression<E2> &e2,
                           size_type i) {
            size_type size = BOOST_UBLAS_SAME (e1 ().size2 (), e2 ().size ());
            result_type t = result_type (0);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            for (size_type j = 0; j < size; ++ j)
                t += e1 () (i, j) * e2 () (j);
#else
            size_type j (0);
            DD (size, 4, r, (t += e1 () (i, j) * e2 () (j), ++ j));
#endif
            return t;
        }
        // Dense case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (difference_type size, I1 it1, I2 it2) {
            result_type t = result_type (0);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            while (-- size >= 0)
                t += *it1 * *it2, ++ it1, ++ it2;
#else
            DD (size, 4, r, (t += *it1 * *it2, ++ it1, ++ it2));
#endif
            return t;
        }
        // Packed case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end) {
            result_type t = result_type (0);
            difference_type it1_size (it1_end - it1);
            difference_type it2_size (it2_end - it2);
            difference_type diff (0);
            if (it1_size > 0 && it2_size > 0)
                diff = it2.index () - it1.index2 ();
            if (diff != 0) {
                difference_type size = (std::min) (diff, it1_size);
                if (size > 0) {
                    it1 += size;
                    it1_size -= size;
                    diff -= size;
                }
                size = (std::min) (- diff, it2_size);
                if (size > 0) {
                    it2 += size;
                    it2_size -= size;
                    diff += size;
                }
            }
            difference_type size ((std::min) (it1_size, it2_size));
            while (-- size >= 0)
                t += *it1 * *it2, ++ it1, ++ it2;
            return t;
        }
        // Sparse case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end,
                           sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag) {
            result_type t = result_type (0);
            if (it1 != it1_end && it2 != it2_end) {
                size_type it1_index = it1.index2 (), it2_index = it2.index ();
                for (;;) {
                    difference_type compare = it1_index - it2_index;
                    if (compare == 0) {
                        t += *it1 * *it2, ++ it1, ++ it2;
                        if (it1 != it1_end && it2 != it2_end) {
                            it1_index = it1.index2 ();
                            it2_index = it2.index ();
                        } else
                            break;
                    } else if (compare < 0) {
                        increment (it1, it1_end, - compare);
                        if (it1 != it1_end)
                            it1_index = it1.index2 ();
                        else
                            break;
                    } else if (compare > 0) {
                        increment (it2, it2_end, compare);
                        if (it2 != it2_end)
                            it2_index = it2.index ();
                        else
                            break;
                    }
                }
            }
            return t;
        }
        // Sparse packed case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &/* it2_end */,
                           sparse_bidirectional_iterator_tag, packed_random_access_iterator_tag) {
            result_type t = result_type (0);
            while (it1 != it1_end) {
                t += *it1 * it2 () (it1.index2 ());
                ++ it1;
            }
            return t;
        }
        // Packed sparse case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &/* it1_end */, I2 it2, const I2 &it2_end,
                           packed_random_access_iterator_tag, sparse_bidirectional_iterator_tag) {
            result_type t = result_type (0);
            while (it2 != it2_end) {
                t += it1 () (it1.index1 (), it2.index ()) * *it2;
                ++ it2;
            }
            return t;
        }
        // Another dispatcher
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end,
                           sparse_bidirectional_iterator_tag) {
            typedef typename I1::iterator_category iterator1_category;
            typedef typename I2::iterator_category iterator2_category;
            return apply (it1, it1_end, it2, it2_end, iterator1_category (), iterator2_category ());
        }
    };

    template<class M1, class M2, class TV>
    struct matrix_vector_prod2:
        public matrix_vector_binary_functor<M1, M2, TV> {
        typedef typename matrix_vector_binary_functor<M1, M2, TV>::size_type size_type;
        typedef typename matrix_vector_binary_functor<M1, M2, TV>::difference_type difference_type;
        typedef typename matrix_vector_binary_functor<M1, M2, TV>::value_type value_type;
        typedef typename matrix_vector_binary_functor<M1, M2, TV>::result_type result_type;

        template<class C1, class C2>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_container<C1> &c1,
                           const matrix_container<C2> &c2,
                           size_type i) {
#ifdef BOOST_UBLAS_USE_SIMD
            using namespace raw;
            size_type size = BOOST_UBLAS_SAME (c1 ().size (), c2 ().size1 ());
            const typename M1::value_type *data1 = data_const (c1 ());
            const typename M2::value_type *data2 = data_const (c2 ()) + i * stride2 (c2 ());
            size_type s1 = stride (c1 ());
            size_type s2 = stride1 (c2 ());
            result_type t = result_type (0);
            if (s1 == 1 && s2 == 1) {
                for (size_type j = 0; j < size; ++ j)
                    t += data1 [j] * data2 [j];
            } else if (s2 == 1) {
                for (size_type j = 0, j1 = 0; j < size; ++ j, j1 += s1)
                    t += data1 [j1] * data2 [j];
            } else if (s1 == 1) {
                for (size_type j = 0, j2 = 0; j < size; ++ j, j2 += s2)
                    t += data1 [j] * data2 [j2];
            } else {
                for (size_type j = 0, j1 = 0, j2 = 0; j < size; ++ j, j1 += s1, j2 += s2)
                    t += data1 [j1] * data2 [j2];
            }
            return t;
#elif defined(BOOST_UBLAS_HAVE_BINDINGS)
            return boost::numeric::bindings::atlas::dot (c1 (), c2 ().column (i));
#else
            return apply (static_cast<const vector_expression<C1> > (c1), static_cast<const matrix_expression<C2> > (c2, i));
#endif
        }
        template<class E1, class E2>
        static BOOST_UBLAS_INLINE
        result_type apply (const vector_expression<E1> &e1,
                           const matrix_expression<E2> &e2,
                           size_type i) {
            size_type size = BOOST_UBLAS_SAME (e1 ().size (), e2 ().size1 ());
            result_type t = result_type (0);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            for (size_type j = 0; j < size; ++ j)
                t += e1 () (j) * e2 () (j, i);
#else
            size_type j (0);
            DD (size, 4, r, (t += e1 () (j) * e2 () (j, i), ++ j));
#endif
            return t;
        }
        // Dense case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (difference_type size, I1 it1, I2 it2) {
            result_type t = result_type (0);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            while (-- size >= 0)
                t += *it1 * *it2, ++ it1, ++ it2;
#else
            DD (size, 4, r, (t += *it1 * *it2, ++ it1, ++ it2));
#endif
            return t;
        }
        // Packed case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end) {
            result_type t = result_type (0);
            difference_type it1_size (it1_end - it1);
            difference_type it2_size (it2_end - it2);
            difference_type diff (0);
            if (it1_size > 0 && it2_size > 0)
                diff = it2.index1 () - it1.index ();
            if (diff != 0) {
                difference_type size = (std::min) (diff, it1_size);
                if (size > 0) {
                    it1 += size;
                    it1_size -= size;
                    diff -= size;
                }
                size = (std::min) (- diff, it2_size);
                if (size > 0) {
                    it2 += size;
                    it2_size -= size;
                    diff += size;
                }
            }
            difference_type size ((std::min) (it1_size, it2_size));
            while (-- size >= 0)
                t += *it1 * *it2, ++ it1, ++ it2;
            return t;
        }
        // Sparse case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end,
                           sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag) {
            result_type t = result_type (0);
            if (it1 != it1_end && it2 != it2_end) {
                size_type it1_index = it1.index (), it2_index = it2.index1 ();
                for (;;) {
                    difference_type compare = it1_index - it2_index;
                    if (compare == 0) {
                        t += *it1 * *it2, ++ it1, ++ it2;
                        if (it1 != it1_end && it2 != it2_end) {
                            it1_index = it1.index ();
                            it2_index = it2.index1 ();
                        } else
                            break;
                    } else if (compare < 0) {
                        increment (it1, it1_end, - compare);
                        if (it1 != it1_end)
                            it1_index = it1.index ();
                        else
                            break;
                    } else if (compare > 0) {
                        increment (it2, it2_end, compare);
                        if (it2 != it2_end)
                            it2_index = it2.index1 ();
                        else
                            break;
                    }
                }
            }
            return t;
        }
        // Packed sparse case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &/* it1_end */, I2 it2, const I2 &it2_end,
                           packed_random_access_iterator_tag, sparse_bidirectional_iterator_tag) {
            result_type t = result_type (0);
            while (it2 != it2_end) {
                t += it1 () (it2.index1 ()) * *it2;
                ++ it2;
            }
            return t;
        }
        // Sparse packed case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &/* it2_end */,
                           sparse_bidirectional_iterator_tag, packed_random_access_iterator_tag) {
            result_type t = result_type (0);
            while (it1 != it1_end) {
                t += *it1 * it2 () (it1.index (), it2.index2 ());
                ++ it1;
            }
            return t;
        }
        // Another dispatcher
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end,
                           sparse_bidirectional_iterator_tag) {
            typedef typename I1::iterator_category iterator1_category;
            typedef typename I2::iterator_category iterator2_category;
            return apply (it1, it1_end, it2, it2_end, iterator1_category (), iterator2_category ());
        }
    };

    // Binary returning matrix
    template<class M1, class M2, class TV>
    struct matrix_matrix_binary_functor {
        typedef typename M1::size_type size_type;
        typedef typename M1::difference_type difference_type;
        typedef TV value_type;
        typedef TV result_type;
    };

    template<class M1, class M2, class TV>
    struct matrix_matrix_prod:
        public matrix_matrix_binary_functor<M1, M2, TV> {
        typedef typename matrix_matrix_binary_functor<M1, M2, TV>::size_type size_type;
        typedef typename matrix_matrix_binary_functor<M1, M2, TV>::difference_type difference_type;
        typedef typename matrix_matrix_binary_functor<M1, M2, TV>::value_type value_type;
        typedef typename matrix_matrix_binary_functor<M1, M2, TV>::result_type result_type;

        template<class C1, class C2>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_container<C1> &c1,
                           const matrix_container<C2> &c2,
                           size_type i, size_type j) {
#ifdef BOOST_UBLAS_USE_SIMD
            using namespace raw;
            size_type size = BOOST_UBLAS_SAME (c1 ().size2 (), c2 ().sizc1 ());
            const typename M1::value_type *data1 = data_const (c1 ()) + i * stride1 (c1 ());
            const typename M2::value_type *data2 = data_const (c2 ()) + j * stride2 (c2 ());
            size_type s1 = stride2 (c1 ());
            size_type s2 = stride1 (c2 ());
            result_type t = result_type (0);
            if (s1 == 1 && s2 == 1) {
                for (size_type k = 0; k < size; ++ k)
                    t += data1 [k] * data2 [k];
            } else if (s2 == 1) {
                for (size_type k = 0, k1 = 0; k < size; ++ k, k1 += s1)
                    t += data1 [k1] * data2 [k];
            } else if (s1 == 1) {
                for (size_type k = 0, k2 = 0; k < size; ++ k, k2 += s2)
                    t += data1 [k] * data2 [k2];
            } else {
                for (size_type k = 0, k1 = 0, k2 = 0; k < size; ++ k, k1 += s1, k2 += s2)
                    t += data1 [k1] * data2 [k2];
            }
            return t;
#elif defined(BOOST_UBLAS_HAVE_BINDINGS)
            return boost::numeric::bindings::atlas::dot (c1 ().row (i), c2 ().column (j));
#else
            boost::ignore_unused(j);
            return apply (static_cast<const matrix_expression<C1> > (c1), static_cast<const matrix_expression<C2> > (c2, i));
#endif
        }
        template<class E1, class E2>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E1> &e1,
                           const matrix_expression<E2> &e2,
                           size_type i, size_type j) {
            size_type size = BOOST_UBLAS_SAME (e1 ().size2 (), e2 ().size1 ());
            result_type t = result_type (0);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            for (size_type k = 0; k < size; ++ k)
                t += e1 () (i, k) * e2 () (k, j);
#else
            size_type k (0);
            DD (size, 4, r, (t += e1 () (i, k) * e2 () (k, j), ++ k));
#endif
            return t;
        }
        // Dense case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (difference_type size, I1 it1, I2 it2) {
            result_type t = result_type (0);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            while (-- size >= 0)
                t += *it1 * *it2, ++ it1, ++ it2;
#else
            DD (size, 4, r, (t += *it1 * *it2, ++ it1, ++ it2));
#endif
            return t;
        }
        // Packed case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end, packed_random_access_iterator_tag) {
            result_type t = result_type (0);
            difference_type it1_size (it1_end - it1);
            difference_type it2_size (it2_end - it2);
            difference_type diff (0);
            if (it1_size > 0 && it2_size > 0)
                diff = it2.index1 () - it1.index2 ();
            if (diff != 0) {
                difference_type size = (std::min) (diff, it1_size);
                if (size > 0) {
                    it1 += size;
                    it1_size -= size;
                    diff -= size;
                }
                size = (std::min) (- diff, it2_size);
                if (size > 0) {
                    it2 += size;
                    it2_size -= size;
                    diff += size;
                }
            }
            difference_type size ((std::min) (it1_size, it2_size));
            while (-- size >= 0)
                t += *it1 * *it2, ++ it1, ++ it2;
            return t;
        }
        // Sparse case
        template<class I1, class I2>
        static BOOST_UBLAS_INLINE
        result_type apply (I1 it1, const I1 &it1_end, I2 it2, const I2 &it2_end, sparse_bidirectional_iterator_tag) {
            result_type t = result_type (0);
            if (it1 != it1_end && it2 != it2_end) {
                size_type it1_index = it1.index2 (), it2_index = it2.index1 ();
                for (;;) {
                    difference_type compare = difference_type (it1_index - it2_index);
                    if (compare == 0) {
                        t += *it1 * *it2, ++ it1, ++ it2;
                        if (it1 != it1_end && it2 != it2_end) {
                            it1_index = it1.index2 ();
                            it2_index = it2.index1 ();
                        } else
                            break;
                    } else if (compare < 0) {
                        increment (it1, it1_end, - compare);
                        if (it1 != it1_end)
                            it1_index = it1.index2 ();
                        else
                            break;
                    } else if (compare > 0) {
                        increment (it2, it2_end, compare);
                        if (it2 != it2_end)
                            it2_index = it2.index1 ();
                        else
                            break;
                    }
                }
            }
            return t;
        }
    };

    // Unary returning scalar norm
    template<class M>
    struct matrix_scalar_real_unary_functor {
        typedef typename M::value_type value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef real_type result_type;
    };

    template<class M>
    struct matrix_norm_1:
        public matrix_scalar_real_unary_functor<M> {
        typedef typename matrix_scalar_real_unary_functor<M>::value_type value_type;
        typedef typename matrix_scalar_real_unary_functor<M>::real_type real_type;
        typedef typename matrix_scalar_real_unary_functor<M>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) {
            real_type t = real_type ();
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size2 (e ().size2 ());
            for (matrix_size_type j = 0; j < size2; ++ j) {
                real_type u = real_type ();
                matrix_size_type size1 (e ().size1 ());
                for (matrix_size_type i = 0; i < size1; ++ i) {
                    real_type v (type_traits<value_type>::norm_1 (e () (i, j)));
                    u += v;
                }
                if (u > t)
                    t = u;
            }
            return t; 
        }
    };

    template<class M>
    struct matrix_norm_frobenius:
        public matrix_scalar_real_unary_functor<M> {
        typedef typename matrix_scalar_real_unary_functor<M>::value_type value_type;
        typedef typename matrix_scalar_real_unary_functor<M>::real_type real_type;
        typedef typename matrix_scalar_real_unary_functor<M>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) { 
            real_type t = real_type ();
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            for (matrix_size_type i = 0; i < size1; ++ i) {
                matrix_size_type size2 (e ().size2 ());
                for (matrix_size_type j = 0; j < size2; ++ j) {
                    real_type u (type_traits<value_type>::norm_2 (e () (i, j)));
                    t +=  u * u;
                }
            }
            return type_traits<real_type>::type_sqrt (t); 
        }
    };

    template<class M>
    struct matrix_norm_inf: 
        public matrix_scalar_real_unary_functor<M> {
        typedef typename matrix_scalar_real_unary_functor<M>::value_type value_type;
        typedef typename matrix_scalar_real_unary_functor<M>::real_type real_type;
        typedef typename matrix_scalar_real_unary_functor<M>::result_type result_type;

        template<class E>
        static BOOST_UBLAS_INLINE
        result_type apply (const matrix_expression<E> &e) {
            real_type t = real_type ();
            typedef typename E::size_type matrix_size_type;
            matrix_size_type size1 (e ().size1 ());
            for (matrix_size_type i = 0; i < size1; ++ i) {
                real_type u = real_type ();
                matrix_size_type size2 (e ().size2 ());
                for (matrix_size_type j = 0; j < size2; ++ j) {
                    real_type v (type_traits<value_type>::norm_inf (e () (i, j)));
                    u += v;
                }
                if (u > t) 
                    t = u;  
            }
            return t; 
        }
    };

    // forward declaration
    template <class Z, class D> struct basic_column_major;

    // This functor defines storage layout and it's properties
    // matrix (i,j) -> storage [i * size_i + j]
    template <class Z, class D>
    struct basic_row_major {
        typedef Z size_type;
        typedef D difference_type;
        typedef row_major_tag orientation_category;
        typedef basic_column_major<Z,D> transposed_layout;

        static
        BOOST_UBLAS_INLINE
        size_type storage_size (size_type size_i, size_type size_j) {
            // Guard against size_type overflow
            BOOST_UBLAS_CHECK (size_j == 0 || size_i <= (std::numeric_limits<size_type>::max) () / size_j, bad_size ());
            return size_i * size_j;
        }

        // Indexing conversion to storage element
        static
        BOOST_UBLAS_INLINE
        size_type element (size_type i, size_type size_i, size_type j, size_type size_j) {
            BOOST_UBLAS_CHECK (i < size_i, bad_index ());
            BOOST_UBLAS_CHECK (j < size_j, bad_index ());
            detail::ignore_unused_variable_warning(size_i);
            // Guard against size_type overflow
            BOOST_UBLAS_CHECK (i <= ((std::numeric_limits<size_type>::max) () - j) / size_j, bad_index ());
            return i * size_j + j;
        }
        static
        BOOST_UBLAS_INLINE
        size_type address (size_type i, size_type size_i, size_type j, size_type size_j) {
            BOOST_UBLAS_CHECK (i <= size_i, bad_index ());
            BOOST_UBLAS_CHECK (j <= size_j, bad_index ());
            // Guard against size_type overflow - address may be size_j past end of storage
            BOOST_UBLAS_CHECK (size_j == 0 || i <= ((std::numeric_limits<size_type>::max) () - j) / size_j, bad_index ());
            detail::ignore_unused_variable_warning(size_i);
            return i * size_j + j;
        }

        // Storage element to index conversion
        static
        BOOST_UBLAS_INLINE
        difference_type distance_i (difference_type k, size_type /* size_i */, size_type size_j) {
            return size_j != 0 ? k / size_j : 0;
        }
        static
        BOOST_UBLAS_INLINE
        difference_type distance_j (difference_type k, size_type /* size_i */, size_type /* size_j */) {
            return k;
        }
        static
        BOOST_UBLAS_INLINE
        size_type index_i (difference_type k, size_type /* size_i */, size_type size_j) {
            return size_j != 0 ? k / size_j : 0;
        }
        static
        BOOST_UBLAS_INLINE
        size_type index_j (difference_type k, size_type /* size_i */, size_type size_j) {
            return size_j != 0 ? k % size_j : 0;
        }
        static
        BOOST_UBLAS_INLINE
        bool fast_i () {
            return false;
        }
        static
        BOOST_UBLAS_INLINE
        bool fast_j () {
            return true;
        }

        // Iterating storage elements
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void increment_i (I &it, size_type /* size_i */, size_type size_j) {
            it += size_j;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void increment_i (I &it, difference_type n, size_type /* size_i */, size_type size_j) {
            it += n * size_j;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void decrement_i (I &it, size_type /* size_i */, size_type size_j) {
            it -= size_j;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void decrement_i (I &it, difference_type n, size_type /* size_i */, size_type size_j) {
            it -= n * size_j;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void increment_j (I &it, size_type /* size_i */, size_type /* size_j */) {
            ++ it;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void increment_j (I &it, difference_type n, size_type /* size_i */, size_type /* size_j */) {
            it += n;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void decrement_j (I &it, size_type /* size_i */, size_type /* size_j */) {
            -- it;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void decrement_j (I &it, difference_type n, size_type /* size_i */, size_type /* size_j */) {
            it -= n;
        }

        // Triangular access
        static
        BOOST_UBLAS_INLINE
        size_type triangular_size (size_type size_i, size_type size_j) {
            size_type size = (std::max) (size_i, size_j);
            // Guard against size_type overflow - simplified
            BOOST_UBLAS_CHECK (size == 0 || size / 2 < (std::numeric_limits<size_type>::max) () / size /* +1/2 */, bad_size ());
            return ((size + 1) * size) / 2;
        }
        static
        BOOST_UBLAS_INLINE
        size_type lower_element (size_type i, size_type size_i, size_type j, size_type size_j) {
            BOOST_UBLAS_CHECK (i < size_i, bad_index ());
            BOOST_UBLAS_CHECK (j < size_j, bad_index ());
            BOOST_UBLAS_CHECK (i >= j, bad_index ());
            detail::ignore_unused_variable_warning(size_i);
            detail::ignore_unused_variable_warning(size_j);
            // FIXME size_type overflow
            // sigma_i (i + 1) = (i + 1) * i / 2
            // i = 0 1 2 3, sigma = 0 1 3 6
            return ((i + 1) * i) / 2 + j;
        }
        static
        BOOST_UBLAS_INLINE
        size_type upper_element (size_type i, size_type size_i, size_type j, size_type size_j) {
            BOOST_UBLAS_CHECK (i < size_i, bad_index ());
            BOOST_UBLAS_CHECK (j < size_j, bad_index ());
            BOOST_UBLAS_CHECK (i <= j, bad_index ());
            // FIXME size_type overflow
            // sigma_i (size - i) = size * i - i * (i - 1) / 2
            // i = 0 1 2 3, sigma = 0 4 7 9
            return (i * (2 * (std::max) (size_i, size_j) - i + 1)) / 2 + j - i;
        }

        // Major and minor indices
        static
        BOOST_UBLAS_INLINE
        size_type index_M (size_type index1, size_type /* index2 */) {
            return index1;
        }
        static
        BOOST_UBLAS_INLINE
        size_type index_m (size_type /* index1 */, size_type index2) {
            return index2;
        }
        static
        BOOST_UBLAS_INLINE
        size_type size_M (size_type size_i, size_type /* size_j */) {
            return size_i;
        }
        static
        BOOST_UBLAS_INLINE
        size_type size_m (size_type /* size_i */, size_type size_j) {
            return size_j;
        }
    };

    // This functor defines storage layout and it's properties
    // matrix (i,j) -> storage [i + j * size_i]
    template <class Z, class D>
    struct basic_column_major {
        typedef Z size_type;
        typedef D difference_type;
        typedef column_major_tag orientation_category;
        typedef basic_row_major<Z,D> transposed_layout;

        static
        BOOST_UBLAS_INLINE
        size_type storage_size (size_type size_i, size_type size_j) {
            // Guard against size_type overflow
            BOOST_UBLAS_CHECK (size_i == 0 || size_j <= (std::numeric_limits<size_type>::max) () / size_i, bad_size ());
            return size_i * size_j;
        }

        // Indexing conversion to storage element
        static
        BOOST_UBLAS_INLINE
        size_type element (size_type i, size_type size_i, size_type j, size_type size_j) {
            BOOST_UBLAS_CHECK (i < size_i, bad_index ());
            BOOST_UBLAS_CHECK (j < size_j, bad_index ());
            detail::ignore_unused_variable_warning(size_j);
            // Guard against size_type overflow
            BOOST_UBLAS_CHECK (j <= ((std::numeric_limits<size_type>::max) () - i) / size_i, bad_index ());
            return i + j * size_i;
        }
        static
        BOOST_UBLAS_INLINE
        size_type address (size_type i, size_type size_i, size_type j, size_type size_j) {
            BOOST_UBLAS_CHECK (i <= size_i, bad_index ());
            BOOST_UBLAS_CHECK (j <= size_j, bad_index ());
            detail::ignore_unused_variable_warning(size_j);
            // Guard against size_type overflow - address may be size_i past end of storage
            BOOST_UBLAS_CHECK (size_i == 0 || j <= ((std::numeric_limits<size_type>::max) () - i) / size_i, bad_index ());
            return i + j * size_i;
        }

        // Storage element to index conversion
        static
        BOOST_UBLAS_INLINE
        difference_type distance_i (difference_type k, size_type /* size_i */, size_type /* size_j */) {
            return k;
        }
        static
        BOOST_UBLAS_INLINE
        difference_type distance_j (difference_type k, size_type size_i, size_type /* size_j */) {
            return size_i != 0 ? k / size_i : 0;
        }
        static
        BOOST_UBLAS_INLINE
        size_type index_i (difference_type k, size_type size_i, size_type /* size_j */) {
            return size_i != 0 ? k % size_i : 0;
        }
        static
        BOOST_UBLAS_INLINE
        size_type index_j (difference_type k, size_type size_i, size_type /* size_j */) {
            return size_i != 0 ? k / size_i : 0;
        }
        static
        BOOST_UBLAS_INLINE
        bool fast_i () {
            return true;
        }
        static
        BOOST_UBLAS_INLINE
        bool fast_j () {
            return false;
        }

        // Iterating
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void increment_i (I &it, size_type /* size_i */, size_type /* size_j */) {
            ++ it;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void increment_i (I &it, difference_type n, size_type /* size_i */, size_type /* size_j */) {
            it += n;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void decrement_i (I &it, size_type /* size_i */, size_type /* size_j */) {
            -- it;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void decrement_i (I &it, difference_type n, size_type /* size_i */, size_type /* size_j */) {
            it -= n;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void increment_j (I &it, size_type size_i, size_type /* size_j */) {
            it += size_i;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void increment_j (I &it, difference_type n, size_type size_i, size_type /* size_j */) {
            it += n * size_i;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void decrement_j (I &it, size_type size_i, size_type /* size_j */) {
            it -= size_i;
        }
        template<class I>
        static
        BOOST_UBLAS_INLINE
        void decrement_j (I &it, difference_type n, size_type size_i, size_type /* size_j */) {
            it -= n* size_i;
        }

        // Triangular access
        static
        BOOST_UBLAS_INLINE
        size_type triangular_size (size_type size_i, size_type size_j) {
            size_type size = (std::max) (size_i, size_j);
            // Guard against size_type overflow - simplified
            BOOST_UBLAS_CHECK (size == 0 || size / 2 < (std::numeric_limits<size_type>::max) () / size /* +1/2 */, bad_size ());
            return ((size + 1) * size) / 2;
        }
        static
        BOOST_UBLAS_INLINE
        size_type lower_element (size_type i, size_type size_i, size_type j, size_type size_j) {
            BOOST_UBLAS_CHECK (i < size_i, bad_index ());
            BOOST_UBLAS_CHECK (j < size_j, bad_index ());
            BOOST_UBLAS_CHECK (i >= j, bad_index ());
            // FIXME size_type overflow
            // sigma_j (size - j) = size * j - j * (j - 1) / 2
            // j = 0 1 2 3, sigma = 0 4 7 9
            return i - j + (j * (2 * (std::max) (size_i, size_j) - j + 1)) / 2;
        }
        static
        BOOST_UBLAS_INLINE
        size_type upper_element (size_type i, size_type size_i, size_type j, size_type size_j) {
            BOOST_UBLAS_CHECK (i < size_i, bad_index ());
            BOOST_UBLAS_CHECK (j < size_j, bad_index ());
            BOOST_UBLAS_CHECK (i <= j, bad_index ());
            // FIXME size_type overflow
            // sigma_j (j + 1) = (j + 1) * j / 2
            // j = 0 1 2 3, sigma = 0 1 3 6
            return i + ((j + 1) * j) / 2;
        }

        // Major and minor indices
        static
        BOOST_UBLAS_INLINE
        size_type index_M (size_type /* index1 */, size_type index2) {
            return index2;
        }
        static
        BOOST_UBLAS_INLINE
        size_type index_m (size_type index1, size_type /* index2 */) {
            return index1;
        }
        static
        BOOST_UBLAS_INLINE
        size_type size_M (size_type /* size_i */, size_type size_j) {
            return size_j;
        }
        static
        BOOST_UBLAS_INLINE
        size_type size_m (size_type size_i, size_type /* size_j */) {
            return size_i;
        }
    };


    template <class Z>
    struct basic_full {
        typedef Z size_type;

        template<class L>
        static
        BOOST_UBLAS_INLINE
        size_type packed_size (L, size_type size_i, size_type size_j) {
            return L::storage_size (size_i, size_j);
        }

        static
        BOOST_UBLAS_INLINE
        bool zero (size_type /* i */, size_type /* j */) {
            return false;
        }
        static
        BOOST_UBLAS_INLINE
        bool one (size_type /* i */, size_type /* j */) {
            return false;
        }
        static
        BOOST_UBLAS_INLINE
        bool other (size_type /* i */, size_type /* j */) {
            return true;
        }
        // FIXME: this should not be used at all
        static
        BOOST_UBLAS_INLINE
        size_type restrict1 (size_type i, size_type /* j */) {
            return i;
        }
        static
        BOOST_UBLAS_INLINE
        size_type restrict2 (size_type /* i */, size_type j) {
            return j;
        }
        static
        BOOST_UBLAS_INLINE
        size_type mutable_restrict1 (size_type i, size_type /* j */) {
            return i;
        }
        static
        BOOST_UBLAS_INLINE
        size_type mutable_restrict2 (size_type /* i */, size_type j) {
            return j;
        }
    };

    namespace detail {
        template < class L >
        struct transposed_structure {
            typedef typename L::size_type size_type;

            template<class LAYOUT>
            static
            BOOST_UBLAS_INLINE
            size_type packed_size (LAYOUT l, size_type size_i, size_type size_j) {
                return L::packed_size(l, size_j, size_i);
            }

            static
            BOOST_UBLAS_INLINE
            bool zero (size_type i, size_type j) {
                return L::zero(j, i);
            }
            static
            BOOST_UBLAS_INLINE
            bool one (size_type i, size_type j) {
                return L::one(j, i);
            }
            static
            BOOST_UBLAS_INLINE
            bool other (size_type i, size_type j) {
                return L::other(j, i);
            }
            template<class LAYOUT>
            static
            BOOST_UBLAS_INLINE
            size_type element (LAYOUT /* l */, size_type i, size_type size_i, size_type j, size_type size_j) {
                return L::element(typename LAYOUT::transposed_layout(), j, size_j, i, size_i);
            }

            static
            BOOST_UBLAS_INLINE
            size_type restrict1 (size_type i, size_type j, size_type size1, size_type size2) {
                return L::restrict2(j, i, size2, size1);
            }
            static
            BOOST_UBLAS_INLINE
            size_type restrict2 (size_type i, size_type j, size_type size1, size_type size2) {
                return L::restrict1(j, i, size2, size1);
            }
            static
            BOOST_UBLAS_INLINE
            size_type mutable_restrict1 (size_type i, size_type j, size_type size1, size_type size2) {
                return L::mutable_restrict2(j, i, size2, size1);
            }
            static
            BOOST_UBLAS_INLINE
            size_type mutable_restrict2 (size_type i, size_type j, size_type size1, size_type size2) {
                return L::mutable_restrict1(j, i, size2, size1);
            }

            static
            BOOST_UBLAS_INLINE
            size_type global_restrict1 (size_type index1, size_type size1, size_type index2, size_type size2) {
                return L::global_restrict2(index2, size2, index1, size1);
            }
            static
            BOOST_UBLAS_INLINE
            size_type global_restrict2 (size_type index1, size_type size1, size_type index2, size_type size2) {
                return L::global_restrict1(index2, size2, index1, size1);
            }
            static
            BOOST_UBLAS_INLINE
            size_type global_mutable_restrict1 (size_type index1, size_type size1, size_type index2, size_type size2) {
                return L::global_mutable_restrict2(index2, size2, index1, size1);
            }
            static
            BOOST_UBLAS_INLINE
            size_type global_mutable_restrict2 (size_type index1, size_type size1, size_type index2, size_type size2) {
                return L::global_mutable_restrict1(index2, size2, index1, size1);
            }
        };
    }

    template <class Z>
    struct basic_lower {
        typedef Z size_type;
        typedef lower_tag triangular_type;

        template<class L>
        static
        BOOST_UBLAS_INLINE
        size_type packed_size (L, size_type size_i, size_type size_j) {
            return L::triangular_size (size_i, size_j);
        }

        static
        BOOST_UBLAS_INLINE
        bool zero (size_type i, size_type j) {
            return j > i;
        }
        static
        BOOST_UBLAS_INLINE
        bool one (size_type /* i */, size_type /* j */) {
            return false;
        }
        static
        BOOST_UBLAS_INLINE
        bool other (size_type i, size_type j) {
            return j <= i;
        }
        template<class L>
        static
        BOOST_UBLAS_INLINE
        size_type element (L, size_type i, size_type size_i, size_type j, size_type size_j) {
            return L::lower_element (i, size_i, j, size_j);
        }

        // return nearest valid index in column j
        static
        BOOST_UBLAS_INLINE
        size_type restrict1 (size_type i, size_type j, size_type size1, size_type /* size2 */) {
            return (std::max)(j, (std::min) (size1, i));
        }
        // return nearest valid index in row i
        static
        BOOST_UBLAS_INLINE
        size_type restrict2 (size_type i, size_type j, size_type /* size1 */, size_type /* size2 */) {
            return (std::max)(size_type(0), (std::min) (i+1, j));
        }
        // return nearest valid mutable index in column j
        static
        BOOST_UBLAS_INLINE
        size_type mutable_restrict1 (size_type i, size_type j, size_type size1, size_type /* size2 */) {
            return (std::max)(j, (std::min) (size1, i));
        }
        // return nearest valid mutable index in row i
        static
        BOOST_UBLAS_INLINE
        size_type mutable_restrict2 (size_type i, size_type j, size_type /* size1 */, size_type /* size2 */) {
            return (std::max)(size_type(0), (std::min) (i+1, j));
        }

        // return an index between the first and (1+last) filled row
        static
        BOOST_UBLAS_INLINE
        size_type global_restrict1 (size_type index1, size_type size1, size_type /* index2 */, size_type /* size2 */) {
            return (std::max)(size_type(0), (std::min)(size1, index1) );
        }
        // return an index between the first and (1+last) filled column
        static
        BOOST_UBLAS_INLINE
        size_type global_restrict2 (size_type /* index1 */, size_type /* size1 */, size_type index2, size_type size2) {
            return (std::max)(size_type(0), (std::min)(size2, index2) );
        }

        // return an index between the first and (1+last) filled mutable row
        static
        BOOST_UBLAS_INLINE
        size_type global_mutable_restrict1 (size_type index1, size_type size1, size_type /* index2 */, size_type /* size2 */) {
            return (std::max)(size_type(0), (std::min)(size1, index1) );
        }
        // return an index between the first and (1+last) filled mutable column
        static
        BOOST_UBLAS_INLINE
        size_type global_mutable_restrict2 (size_type /* index1 */, size_type /* size1 */, size_type index2, size_type size2) {
            return (std::max)(size_type(0), (std::min)(size2, index2) );
        }
    };

    // the first row only contains a single 1. Thus it is not stored.
    template <class Z>
    struct basic_unit_lower : public basic_lower<Z> {
        typedef Z size_type;
        typedef unit_lower_tag triangular_type;

        template<class L>
        static
        BOOST_UBLAS_INLINE
        size_type packed_size (L, size_type size_i, size_type size_j) {
            // Zero size strict triangles are bad at this point
            BOOST_UBLAS_CHECK (size_i != 0 && size_j != 0, bad_index ());
            return L::triangular_size (size_i - 1, size_j - 1);
        }

        static
        BOOST_UBLAS_INLINE
        bool one (size_type i, size_type j) {
            return j == i;
        }
        static
        BOOST_UBLAS_INLINE
        bool other (size_type i, size_type j) {
            return j < i;
        }
        template<class L>
        static
        BOOST_UBLAS_INLINE
        size_type element (L, size_type i, size_type size_i, size_type j, size_type size_j) {
            // Zero size strict triangles are bad at this point
            BOOST_UBLAS_CHECK (size_i != 0 && size_j != 0 && i != 0, bad_index ());
            return L::lower_element (i-1, size_i - 1, j, size_j - 1);
        }

        static
        BOOST_UBLAS_INLINE
        size_type mutable_restrict1 (size_type i, size_type j, size_type size1, size_type /* size2 */) {
            return (std::max)(j+1, (std::min) (size1, i));
        }
        static
        BOOST_UBLAS_INLINE
        size_type mutable_restrict2 (size_type i, size_type j, size_type /* size1 */, size_type /* size2 */) {
            return (std::max)(size_type(0), (std::min) (i, j));
        }

        // return an index between the first and (1+last) filled mutable row
        static
        BOOST_UBLAS_INLINE
        size_type global_mutable_restrict1 (size_type index1, size_type size1, size_type /* index2 */, size_type /* size2 */) {
            return (std::max)(size_type(1), (std::min)(size1, index1) );
        }
        // return an index between the first and (1+last) filled mutable column
        static
        BOOST_UBLAS_INLINE
        size_type global_mutable_restrict2 (size_type /* index1 */, size_type /* size1 */, size_type index2, size_type size2) {
            BOOST_UBLAS_CHECK( size2 >= 1 , external_logic() );
            return (std::max)(size_type(0), (std::min)(size2-1, index2) );
        }
    };

    // the first row only contains no element. Thus it is not stored.
    template <class Z>
    struct basic_strict_lower : public basic_unit_lower<Z> {
        typedef Z size_type;
        typedef strict_lower_tag triangular_type;

        template<class L>
        static
        BOOST_UBLAS_INLINE
        size_type packed_size (L, size_type size_i, size_type size_j) {
            // Zero size strict triangles are bad at this point
            BOOST_UBLAS_CHECK (size_i != 0 && size_j != 0, bad_index ());
            return L::triangular_size (size_i - 1, size_j - 1);
        }

        static
        BOOST_UBLAS_INLINE
        bool zero (size_type i, size_type j) {
            return j >= i;
        }
        static
        BOOST_UBLAS_INLINE
        bool one (size_type /*i*/, size_type /*j*/) {
            return false;
        }
        static
        BOOST_UBLAS_INLINE
        bool other (size_type i, size_type j) {
            return j < i;
        }
        template<class L>
        static
        BOOST_UBLAS_INLINE
        size_type element (L, size_type i, size_type size_i, size_type j, size_type size_j) {
            // Zero size strict triangles are bad at this point
            BOOST_UBLAS_CHECK (size_i != 0 && size_j != 0 && i != 0, bad_index ());
            return L::lower_element (i-1, size_i - 1, j, size_j - 1);
        }

        static
        BOOST_UBLAS_INLINE
        size_type restrict1 (size_type i, size_type j, size_type size1, size_type size2) {
            return basic_unit_lower<Z>::mutable_restrict1(i, j, size1, size2);
        }
        static
        BOOST_UBLAS_INLINE
        size_type restrict2 (size_type i, size_type j, size_type size1, size_type size2) {
            return basic_unit_lower<Z>::mutable_restrict2(i, j, size1, size2);
        }

        // return an index between the first and (1+last) filled row
        static
        BOOST_UBLAS_INLINE
        size_type global_restrict1 (size_type index1, size_type size1, size_type index2, size_type size2) {
            return basic_unit_lower<Z>::global_mutable_restrict1(index1, size1, index2, size2);
        }
        // return an index between the first and (1+last) filled column
        static
        BOOST_UBLAS_INLINE
        size_type global_restrict2 (size_type index1, size_type size1, size_type index2, size_type size2) {
            return basic_unit_lower<Z>::global_mutable_restrict2(index1, size1, index2, size2);
        }
    };


    template <class Z>
    struct basic_upper : public detail::transposed_structure<basic_lower<Z> >
    { 
        typedef upper_tag triangular_type;
    };

    template <class Z>
    struct basic_unit_upper : public detail::transposed_structure<basic_unit_lower<Z> >
    { 
        typedef unit_upper_tag triangular_type;
    };

    template <class Z>
    struct basic_strict_upper : public detail::transposed_structure<basic_strict_lower<Z> >
    { 
        typedef strict_upper_tag triangular_type;
    };


}}}

#endif
