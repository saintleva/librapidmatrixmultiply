#ifndef __LIBRAPIDMATRIXMULTIPLY_SIMPLEDEDUCING_H__
#define __LIBRAPIDMATRIXMULTIPLY_SIMPLEDEDUCING_H__

/*
 * Copyright (C) Anton Liaukevich 2009 <leva.dev@gmail.com>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <ublasaux/TypeReplacer.h>

namespace boost { namespace numeric { namespace ublas { namespace Rapid {


template<class MatrixLike1, class MatrixLike2>
class SimpleDeducing: public TypeReplacer {
private:

    /* Extract matrix type from operations (derived from matrix_expression) */

    template<class Object> struct UncoverOperation {
        typedef Object Answer;
    };

    template<class Matrix, class Functor>
    struct UncoverOperation< matrix_unary1<Matrix,Functor> > {
        typedef Matrix Answer;
    };

    template<class Matrix, class Functor>
    struct UncoverOperation< matrix_unary2<Matrix,Functor> > {
        typedef Matrix Answer;
    };

    template<class Matrix1, class Matrix2, class Functor>
    struct UncoverOperation< matrix_binary<Matrix1,Matrix2,Functor> > {
        typedef typename SimpleDeducing<Matrix1,Matrix2>::Answer Answer; // recursive!

    };

    template<class Scalar, class Matrix, class Functor>
    struct UncoverOperation< matrix_binary_scalar1<Scalar,Matrix,Functor> > {

        typedef typename promote_traits<Scalar, typename Matrix::value_type>::promote_type Item;
        typedef typename Replace<Matrix,Item>::Answer Answer;

    };

    template<class Matrix, class Scalar, class Functor>
    struct UncoverOperation< matrix_binary_scalar2<Matrix,Scalar,Functor> > {

        typedef typename promote_traits<typename Matrix::value_type, Scalar>::promote_type Item;
        typedef typename Replace<Matrix,Item>::Answer Answer;

    };

    template<class Matrix1, class Matrix2, class Functor>
    struct UncoverOperation< matrix_matrix_binary<Matrix1,Matrix2,Functor> > {
        typedef typename SimpleDeducing<Matrix1,Matrix2>::Answer Answer; // recursive!

    };

    /* Extract matrix type from proxy backends */

    template<class Object> struct UncoverProxy {
        typedef Object Answer;
    };

    template<class Matrix>
    struct UncoverProxy< matrix_range<Matrix> > {
        typedef Matrix Answer;
    };

    template<class Matrix>
    struct UncoverProxy< matrix_slice<Matrix> > {
        typedef Matrix Answer;
    };

public:

    typedef typename UncoverProxy<typename UncoverOperation<MatrixLike1>::Answer>::Answer Matrix1;
    typedef typename UncoverProxy<typename UncoverOperation<MatrixLike2>::Answer>::Answer Matrix2;

    typedef typename promote_traits<typename Matrix1::value_type,
                                    typename Matrix2::value_type>::promote_type Item;
    #if 0
    typedef typename Replace<Matrix2,Item>::Answer Answer;
    #else
    typedef Matrix2 Answer;
    #endif

}; //template class SimpleDeducing


}}}} //namespace boost::numeric::ublas::Rapid

#endif //__LIBRAPIDMATRIXMULTIPLY_SIMPLEDEDUCING_H__
