#ifndef __LIBRAPIDMATRIXMULTIPLY_EXCEPTIONS_H__
#define __LIBRAPIDMATRIXMULTIPLY_EXCEPTIONS_H__

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

#include <cstddef>
#include <string>

namespace boost { namespace numeric { namespace ublas { namespace Rapid {


struct Error {
public:
    typedef std::string String;
    virtual String what() const throw() = 0;

protected:
    virtual ~Error() throw();
};

class SizeError: public Error {
public:
    typedef std::size_t Size;

protected:
    virtual ~SizeError() throw();
};

class MatricesHaveUncompartibleSizes: public SizeError {
public:
    inline MatricesHaveUncompartibleSizes(Size m1, Size n1, Size m2, Size n2):
        m1_(m1), n1_(n1), m2_(m2), n2_(n2) {}

    virtual ~MatricesHaveUncompartibleSizes() throw();

    inline Size getM1() const
    {
        return m1_;
    }

    inline Size getN1() const
    {
        return n1_;
    }

    inline Size getM2() const
    {
        return m2_;
    }

    inline Size getN2() const
    {
        return n2_;
    }

private:
    /* Fields */

    Size m1_, n1_, m2_, n2_;
}; //class MatricesHaveUncompartibleSizes

struct MatricesAreNotCoordinated: public MatricesHaveUncompartibleSizes {
    inline MatricesAreNotCoordinated(Size m1, Size n1, Size m2, Size n2):
        MatricesHaveUncompartibleSizes(m1, n1, m2, n2) {}
    virtual ~MatricesAreNotCoordinated() throw();
    virtual String what() const throw();
};

struct MatricesAreNotSizeInverted: public MatricesHaveUncompartibleSizes {
    inline MatricesAreNotSizeInverted(Size m1, Size n1, Size m2, Size n2):
        MatricesHaveUncompartibleSizes(m1, n1, m2, n2) {}
    virtual ~MatricesAreNotSizeInverted() throw();
    virtual String what() const throw();
};

class MatrixSizeIsInsufficientlyRound: public SizeError {
public:
    inline MatrixSizeIsInsufficientlyRound(Size m, Size n):
        m_(m), n_(n) {}

    virtual ~MatrixSizeIsInsufficientlyRound() throw();

    inline Size getM() const
    {
        return m_;
    }

    inline Size getN() const
    {
        return n_;
    }

    virtual String what() const throw();

private:
    /* Fields */

    Size m_, n_;
}; //class MatrixSizeIsInsufficientlyRound


}}}} //namespace boost::numeric::ublas::Rapid

#endif //__LIBRAPIDMATRIXMULTIPLY_EXCEPTIONS_H__
