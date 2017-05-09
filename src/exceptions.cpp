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
#include <boost/format.hpp>
#include <iostream> //TODO: debug the program and remove it

using namespace boost::numeric::ublas::Rapid;


/* class Error */

Error::~Error() throw() {}


/* class SizeError */

SizeError::~SizeError() throw() {}


/* class MatricesHaveUncompartibleSizes */

MatricesHaveUncompartibleSizes::~MatricesHaveUncompartibleSizes() throw() {}


/* class MatricesAreNotCoordinated */

MatricesAreNotCoordinated::~MatricesAreNotCoordinated() throw() {}

Error::String MatricesAreNotCoordinated::what() const throw()
{
    return (boost::format("Matrices %1%x%2% and %3%x%4% are not coordinated")
            % getM1() % getN1() % getM2() % getN2()).str();
}


/* class MatricesAreNotSizeInverted */

MatricesAreNotSizeInverted::~MatricesAreNotSizeInverted() throw() {}

Error::String MatricesAreNotSizeInverted::what() const throw()
{
    return (boost::format("Matrices %1%x%2% and %3%x%4% are not 90-degrees-rotated (by size)")
            % getM1() % getN1() % getM2() % getN2()).str();
}


/* class MatrixSizeIsInsufficientlyRound */

MatrixSizeIsInsufficientlyRound::~MatrixSizeIsInsufficientlyRound() throw() {}

Error::String MatrixSizeIsInsufficientlyRound::what() const throw()
{
    return (boost::format("Matrix size %1%x%2% is insufficiently round") % getM() % getN()).str();
}
