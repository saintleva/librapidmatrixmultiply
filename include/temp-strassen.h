#ifndef __LIBRAPIDMATRIXMULTIPLY_TEMPSTRASSEN_H__
#define __LIBRAPIDMATRIXMULTIPLY_TEMPSTRASSEN_H__

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

#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace RapidMatrixMultiply {


template<class Matrix>
Matrix strassenMultiplyK2n_(const Matrix& matrix1, const Matrix& matrix2,
                            typename Matrix::size_type gaussLimit)
{
    namespace Ublas = boost::numeric::ublas;

    typedef typename Matrix::size_type Size;
    Size size = matrix1.size1();

    if (size <= gaussLimit)
        //TODO: have I to explicitly indicate the namespace?
        return Ublas::prod(matrix1, matrix2);

    Size halfSize = size/2;
    Ublas::range lowerRange(0, halfSize),
                 upperRange(halfSize, size);

    /* Build submatrices of operands */
    Ublas::matrix_range<const Matrix> a1(matrix1, lowerRange, lowerRange),
                                      b1(matrix1, lowerRange, upperRange),
                                      c1(matrix1, upperRange, lowerRange),
                                      d1(matrix1, upperRange, upperRange);
    Ublas::matrix_range<const Matrix> a2(matrix2, lowerRange, lowerRange),
                                      b2(matrix2, lowerRange, upperRange),
                                      c2(matrix2, upperRange, lowerRange),
                                      d2(matrix2, upperRange, upperRange);

    /* Recursively calculate blocks of the resulting matrix */
    Matrix blocks[7];
    /**
     * @todo: Optimize it by canceling converting of range to matrix
     */
    blocks[0] = strassenMultiplyK2n_(Matrix(a1),      Matrix(b2 - d2), gaussLimit);
    blocks[1] = strassenMultiplyK2n_(Matrix(a1 + b1), Matrix(d2),      gaussLimit);
    blocks[2] = strassenMultiplyK2n_(Matrix(c1 + d1), Matrix(a2),      gaussLimit);
    blocks[3] = strassenMultiplyK2n_(Matrix(d1),      Matrix(c2 - a2), gaussLimit);
    blocks[4] = strassenMultiplyK2n_(Matrix(a1 + d1), Matrix(a2 + d2), gaussLimit);
    blocks[5] = strassenMultiplyK2n_(Matrix(b1 - d1), Matrix(c2 + d2), gaussLimit);
    blocks[6] = strassenMultiplyK2n_(Matrix(a1 - c1), Matrix(a2 + b2), gaussLimit);

    /* Create the resulting matrix and its proxies (ranges) */
    Matrix result(size, size);
    Ublas::matrix_range<Matrix> a3(result, lowerRange, lowerRange),
                                b3(result, lowerRange, upperRange),
                                c3(result, upperRange, lowerRange),
                                d3(result, upperRange, upperRange);

    /* Calculate submatrices of the result matrix */
    a3 = blocks[3] + blocks[4] - blocks[1] + blocks[5];
    b3 = blocks[0] + blocks[1];
    c3 = blocks[2] + blocks[3];
    d3 = blocks[4] + blocks[0] - blocks[2] - blocks[6];

    return result;
}

template<class Matrix>
Matrix strassenMultiplyNx2n_(const Matrix& matrix1, const Matrix& matrix2,
                             typename Matrix::size_type gaussLimit)
{
    namespace Ublas = boost::numeric::ublas;

    typedef typename Matrix::size_type Size;
    Size n = matrix1.size1();
    BOOST_UBLAS_CHECK(matrix1.size2() == 2*n, Ublas::bad_size("First matrix is not Nx2N"));
    BOOST_UBLAS_CHECK(matrix2.size1() == 2*n, Ublas::bad_size("Second matrix is not Nx2N"));
    BOOST_UBLAS_CHECK(matrix2.size2() == n, Ublas::bad_size("Second matrix is not Nx2N"));

    if (n <= gaussLimit)
        return Ublas::prod(matrix1, matrix2);

    Ublas::range lowerRange(0, n),
                 upperRange(n, 2*n);

    /* Build submatrices of operands */
    Ublas::matrix_range<const Matrix> a1(matrix1, lowerRange, lowerRange),
                                      b1(matrix1, lowerRange, upperRange);
    Ublas::matrix_range<const Matrix> a2(matrix2, lowerRange, lowerRange),
                                      b2(matrix2, upperRange, lowerRange);

    Matrix tempA1 = a1,
           tempB1 = b1;
    Matrix tempA2 = a2,
           tempB2 = b2;

    return strassenMultiplyK2n_(tempA1, tempA2, gaussLimit) +
           strassenMultiplyK2n_(tempB1, tempB2, gaussLimit);
}

template<class Matrix>
Matrix strassenMultiplySameRects_(const Matrix& matrix1, const Matrix& matrix2,
                                  typename Matrix::size_type gaussLimit)
{
    namespace Ublas = boost::numeric::ublas;

    typedef typename Matrix::size_type Size;
    Size m = matrix1.size1(),
         n = matrix1.size2();

    BOOST_UBLAS_CHECK(matrix2.size1() == n, Ublas::bad_size("Matrices are different by size"));
    BOOST_UBLAS_CHECK(matrix2.size2() == m, Ublas::bad_size("Matrices are different by size"));

    if (m <= gaussLimit || n <= gaussLimit)
        return Ublas::prod(matrix1, matrix2);

    Ublas::range lowerMRange(0, m/2),
                 upperMRange(m/2, m),
                 lowerNRange(0, n/2),
                 upperNRange(n/2, n);

    /* Build submatrices of operands */
    Ublas::matrix_range<const Matrix> a1(matrix1, lowerMRange, lowerNRange),
                                      b1(matrix1, lowerMRange, upperNRange),
                                      c1(matrix1, upperMRange, lowerNRange),
                                      d1(matrix1, upperMRange, upperNRange);
    Ublas::matrix_range<const Matrix> a2(matrix2, lowerNRange, lowerMRange),
                                      b2(matrix2, lowerNRange, upperMRange),
                                      c2(matrix2, upperNRange, lowerMRange),
                                      d2(matrix2, upperNRange, upperMRange);

    /* Calculate left and right multipliers of P[i] */
    Matrix left[7]  = { a1, a1 + b1, c1 + d1, d1, a1 + d1, b1 - d1, a1 - c1 },
           right[7] = { b2 - d2, d2, a2, c2 - a2, a2 + d2, c2 + d2, a2 + b2 };

    /* Calculate P[i] */
    Matrix prod[7];
    for (std::size_t i = 0; i < 7; ++i)
        prod[i] = strassenMultiplySameRects_(left[i], right[i], gaussLimit);

    Matrix result(m, m);
    Ublas::matrix_range<Matrix> a3(result, lowerMRange, lowerMRange),
                                b3(result, lowerMRange, upperMRange),
                                c3(result, upperMRange, lowerMRange),
                                d3(result, upperMRange, upperMRange);

    /* Calculate submatrices of the result matrix */
    a3 = prod[3] + prod[4] - prod[1] + prod[5];
    b3 = prod[0] + prod[1];
    c3 = prod[2] + prod[3];
    d3 = prod[4] + prod[0] - prod[2] - prod[6];

    return result;
}


} //namespace RapidMatrixMultiply

#endif //__LIBRAPIDMATRIXMULTIPLY_TEMPSTRASSEN_H__
