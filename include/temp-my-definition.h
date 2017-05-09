#ifndef __LIBRAPIDMATRIXMULTIPLY_TEMPMYDEFINITION_H__
#define __LIBRAPIDMATRIXMULTIPLY_TEMPMYDEFINITION_H__

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

namespace RapidMatrixMultiply {


template<class Matrix>
Matrix myDefinitionMultiply(const Matrix& a, const Matrix& b)
{
    typedef typename Matrix::size_type Size;
    typedef typename Matrix::value_type Item;
    Size size = a.size1();

    Matrix result(size, size);
    for (Size i = Size(); i < size; ++i)
        for (Size j = Size(); j < size; ++j)
        {
            Item current = Item();
            for (Size k = Size(); k < size; ++k)
                current += a(i, k) * b(k, j);
            result(i, j) = current;
        }

    return result;
}


} //namespace RapidMatrixMultiply

#endif //__LIBRAPIDMATRIXMULTIPLY_TEMPMYDEFINITION_H__
