#ifndef __LIBRAPIDMATRIXMULTIPLY_STRASSENMULTIPLIER_H__
#define __LIBRAPIDMATRIXMULTIPLY_STRASSENMULTIPLIER_H__

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

#include "exceptions.h"
#include "SimpleDeducing.h"
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace boost { namespace numeric { namespace ublas { namespace Rapid {


template<class Matrix1, class Matrix2>
struct DefinitionProd {

    static inline
    typename SimpleDeducing<Matrix1,Matrix2>::Result
    calculate(const matrix_expression<Matrix1>& exp1, const matrix_expression<Matrix2>& exp2)
    {
        return prod(exp1, exp2);
    }

};

template<
         template<class> class IntermediateMatrixPolicy,
         template<class> class SizingPolicy,
         template<class,class> class DeducingPolicy = SimpleDeducing,
         template<class,class> class SmallProdBackend = DefinitionProd
        >
class StrassenMultiplier {
public:

    /* Nested types */

    template<class Size_>
    class SizeChecker: public SizingPolicy<Size_> {
    public:
        /* Types */

        typedef Size_ Size;
        typedef SizingPolicy<Size> Policy;

        static void checkMNxNM2t(Size m1, Size n1, Size m2, Size n2)
        {
            if (m1 != n2 || m2 != n1)
                throw MatricesAreNotSizeInverted(m1, n1, m2, n2);

            Size m = m1,
                 n = n1;
            while (!isNeedToSwitch(m, n))
            {
                if (isOdd(m) || isOdd(n))
                    throw MatrixSizeIsInsufficientlyRound(m, n);
                m /= 2;
                n /= 2;
            }
        }

    private:
        /* Auxiliary functions */

        static inline bool isOdd(Size number)
        {
            return number % 2 == 1;
        }

    }; //class SizeChecker

    /* User-directed multiply functions */

    template<class Matrix1, class Matrix2>
    static inline
    typename DeducingPolicy<Matrix1,Matrix2>::Result
    prodMNxNM2t(const matrix_expression<Matrix1>& exp1, const matrix_expression<Matrix2>& exp2)
    {
        typedef typename DeducingPolicy<Matrix1,Matrix2>::Result Result;
        typedef typename Result::size_type Size;

        SizeChecker<Size>::checkMNxNM2t(exp1().size1(), exp1().size2(),
                                        exp2().size1(), exp2().size2());

        return doProdMNxNM2t(exp1, exp2);
    }

private:

    template<class Matrix1, class Matrix2>
    static
    typename DeducingPolicy<Matrix1,Matrix2>::Result
    doProdMNxNM2t(const matrix_expression<Matrix1>& exp1, const matrix_expression<Matrix2>& exp2)
    {
        typedef typename DeducingPolicy<Matrix1,Matrix2>::Result Result;
        typedef typename Result::size_type Size;
        typedef typename Result::value_type Item;

        const Size m = exp1().size1(),
                   n = exp1().size2();

        if (SizingPolicy<Size>::isNeedToSwitch(m, n))
            return SmallProdBackend<Matrix1,Matrix2>::calculate(exp1, exp2);

        Matrix1 matrix1 = exp1();
        Matrix2 matrix2 = exp2();

        const range lowerMRange(0, m/2),
                    upperMRange(m/2, m),
                    lowerNRange(0, n/2),
                    upperNRange(n/2, n);

        /* Build submatrices of the operands */
        const matrix_range<const Matrix1> a1(matrix1, lowerMRange, lowerNRange),
                                          b1(matrix1, lowerMRange, upperNRange),
                                          c1(matrix1, upperMRange, lowerNRange),
                                          d1(matrix1, upperMRange, upperNRange);
        const matrix_range<const Matrix2> a2(matrix2, lowerNRange, lowerMRange),
                                          b2(matrix2, lowerNRange, upperMRange),
                                          c2(matrix2, upperNRange, lowerMRange),
                                          d2(matrix2, upperNRange, upperMRange);

        typedef typename IntermediateMatrixPolicy<Item>::Matrix Intermediate;

        /* Recursively calculate blocks of the resulting matrix */
        Intermediate blocks[7];
        /**
        * @todo: Optimize it by canceling converting of range to matrix
        */
        #if 0
        blocks[0] = doProdMNxNM2t(Intermediate(a1),      Intermediate(b2 - d2));
        blocks[1] = doProdMNxNM2t(Intermediate(a1 + b1), Intermediate(d2)     );
        blocks[2] = doProdMNxNM2t(Intermediate(c1 + d1), Intermediate(a2)     );
        blocks[3] = doProdMNxNM2t(Intermediate(d1),      Intermediate(c2 - a2));
        blocks[4] = doProdMNxNM2t(Intermediate(a1 + d1), Intermediate(a2 + d2));
        blocks[5] = doProdMNxNM2t(Intermediate(b1 - d1), Intermediate(c2 + d2));
        blocks[6] = doProdMNxNM2t(Intermediate(a1 - c1), Intermediate(a2 + b2));
        #else
        blocks[0] = doProdMNxNM2t(a1,
                                  Intermediate(b2 - d2));
        blocks[1] = doProdMNxNM2t(Intermediate(a1 + b1), d2     );
        blocks[2] = doProdMNxNM2t(Intermediate(c1 + d1), a2     );
        blocks[3] = doProdMNxNM2t(d1,                    Intermediate(c2 - a2));
        blocks[4] = doProdMNxNM2t(Intermediate(a1 + d1), a2 + d2);
        blocks[5] = doProdMNxNM2t(b1 - d1, c2 + d2);
        blocks[6] = doProdMNxNM2t(a1 - c1, a2 + b2);
        #endif

        /* Create the resulting matrix and its proxies (ranges) */
        Result result(m, m);
        matrix_range<Result> a3(result, lowerMRange, lowerMRange),
                             b3(result, lowerMRange, upperMRange),
                             c3(result, upperMRange, lowerMRange),
                             d3(result, upperMRange, upperMRange);

        /* Calculate submatrices of the result matrix */
        a3 = blocks[3] + blocks[4] - blocks[1] + blocks[5];
        b3 = blocks[0] + blocks[1];
        c3 = blocks[2] + blocks[3];
        d3 = blocks[4] + blocks[0] - blocks[2] - blocks[6];

        return result;
    }

}; //template class StrassenMultiplier


}}}} //namespace boost::numeric::ublas::Rapid

#endif //__LIBRAPIDMATRIXMULTIPLY_STRASSENMULTIPLIER_H__
