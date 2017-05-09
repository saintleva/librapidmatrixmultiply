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

#include <boost/type_traits/remove_cv.hpp> //TODO: ?
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <ublasaux/TypeReplacer.h>
#include <boost/static_assert.hpp> //TODO: find error and remove it

namespace boost { namespace numeric { namespace ublas { namespace Rapid {


template<class MatrixLike1, class MatrixLike2>
class SimpleDeducing: public TypeReplacer {
private:

    template<class Matrix1, class Matrix2>
    struct MatrixBinary {
        typedef typename promote_traits<typename Matrix1::value_type,
                                        typename Matrix2::value_type>::promote_type Item;
        typedef typename Replace<Matrix1,Item>::Answer Result;
    };

    /* Extract matrix type from operations (derived from matrix_expression) */

    template<class Object> struct Uncover {
        typedef Object Answer;
    };

    #if 0
    /* For removing const & volatile qualifiers at every stage */

    template<class Object>
    struct Uncover< const Object > {
        typedef typename Uncover<Object>::Answer Answer;
    };

    template<class Object>
    struct Uncover< volatile Object > {
        typedef typename Uncover<Object>::Answer Answer;
    };
    #endif

    /* //TODO */

    #if 0
    template<class Matrix, class Functor>
    struct Uncover< matrix_unary1<Matrix,Functor> > {
        typedef Matrix Current;
        typedef typename Uncover<Current>::Answer Answer;
    };

    template<class Matrix, class Functor>
    struct Uncover< matrix_unary2<Matrix,Functor> > {
        typedef Matrix Current;
        typedef typename Uncover<Current>::Answer Answer;
    };
    #endif

    #if 1
    template<class Matrix1, class Matrix2, class Functor>
    struct Uncover< matrix_binary<Matrix1,Matrix2,Functor> > {

        typedef typename Uncover<Matrix1>::Answer ExpandedMatrix1;
        typedef typename Uncover<Matrix2>::Answer ExpandedMatrix2;
        #if 1
        typedef typename MatrixBinary<ExpandedMatrix1,ExpandedMatrix2>::Result Answer;
        #else
        typedef ExpandedMatrix1 Answer;
        #endif

    };
    #endif

    #if 0
    template<class Scalar, class Matrix, class Functor>
    struct Uncover< matrix_binary_scalar1<Scalar,Matrix,Functor> > {

        typedef typename promote_traits<Scalar, typename Matrix::value_type>::promote_type Item;
        typedef typename Replace<Matrix,Item>::Answer Current;
        typedef typename Uncover<Current>::Answer Answer;

    };

    template<class Matrix, class Scalar, class Functor>
    struct Uncover< matrix_binary_scalar2<Matrix,Scalar,Functor> > {

        typedef typename promote_traits<typename Matrix::value_type, Scalar>::promote_type Item;
        typedef typename Replace<Matrix,Item>::Answer Current;
        typedef typename Uncover<Current>::Answer Answer;

    };

    template<class Matrix1, class Matrix2, class Functor>
    struct Uncover< matrix_matrix_binary<Matrix1,Matrix2,Functor> > {
        typedef typename SimpleDeducing<Matrix1,Matrix2>::Result Current; // recursive!
        typedef typename Uncover<Current>::Answer Answer;
    };
    #endif

    template<class Matrix>
    struct Uncover< matrix_range<Matrix> > {
//        typedef Matrix Current;
        typedef typename Uncover<typename remove_cv<Matrix>::type>::Answer Answer;

    };

    template<class Matrix>
    struct Uncover< matrix_slice<Matrix> > {
//        typedef Matrix Current;
        typedef typename Uncover<Matrix>::Answer Answer;
    };

public:

    #if 0
    typedef typename Uncover<MatrixLike1>::Answer Matrix1;
    typedef typename Uncover<MatrixLike2>::Answer Matrix2;
    #else
    typedef typename Uncover<typename remove_cv<MatrixLike1>::type>::Answer Matrix1;
    typedef typename Uncover<typename remove_cv<MatrixLike2>::type>::Answer Matrix2;
    #endif

    #if 1
    typedef typename MatrixBinary<Matrix1,Matrix2>::Result Result;
    #else
    typedef Matrix2 Result;
    #endif

}; //template class SimpleDeducing


}}}} //namespace boost::numeric::ublas::Rapid

#endif //__LIBRAPIDMATRIXMULTIPLY_SIMPLEDEDUCING_H__
