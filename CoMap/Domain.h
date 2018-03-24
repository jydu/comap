//
// File: Domain.h
// Created by: Julien Dutheil
// Created on: Wed Feb  4 19:12:44 2004
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

using namespace bpp;

/**
 * @brief Defines a set of intervals.
 *
 * Adapted from the Java PSOL library.
 * This class is no longer supported by the Bio++ libraries and is now included here.
 */ 
class Domain
{
  private:
    Vdouble bounds_;
    Vdouble midPoints_;
  
  public:
    /**
     * @brief Build a new Domain object given the lower and upper bound and a number of classes.
     *
     * If a > b, a will be taken as the upper bound and b as the lower.
     * Bounds will be computed to be equally distributed between a and b.
     * Midpoints will be computed to be in the midle of each interval.
     *
     * @param a One bound.
     * @param b Another bound.
     * @param n The number of classes.
     * @throw Exception If the number of classes is 0.
     */
    Domain(double a = 0, double b = 0, size_t n = 1);
    
    /**
     * @brief Build a new Domain object given the classes bounds.
     * 
     * Midpoints will be computed to be in the midle of each interval.
     *
     * @param bounds The classes bounds.
     * @throw Exception If bound are not in increasing order.
     */
    Domain(const Vdouble& bounds);
    
    /**
     * @brief Build a new Domain object given the classes bounds and midpoints.
     *
     * @param bounds    The classes bounds.
     * @param midPoints The classes midpoints.
     * @throw Exception If bound are not in increasing order,
     * or midpoints are not between corresponding bounds,
     * or number of midpoints is not equal to number of bounds - 1.
     */
    Domain(const Vdouble& bounds, const Vdouble& midPoints);
  
    virtual ~Domain() {}
  
  public:
    /**
     * @return The lower bound of the domain.
     */
    virtual double getLowerBound() const { return *bounds_.begin(); }

    /**
     * @return The upper bound of the domain.
     */
    virtual double getUpperBound() const { return *bounds_.rbegin(); }
    
    /**
     * @return The lower midpoint of the domain.
     */
    virtual double getLowerValue() const { return *midPoints_.begin(); }

    /**
     * @return The upper midpoint of the domain.
     */
    virtual double getUpperValue() const { return *midPoints_.rbegin(); }

    /**
     * @return The number of classes = number of intervals = number of midpoints.
     */
    virtual size_t getSize() const { return midPoints_.size(); }

    /**
     * @return All bounds.
     */
    virtual const Vdouble& getBounds() const { return bounds_; }
  
    /**
     * @return One bound.
     * @param i The index of the bound \f$\in [0,n]\f$ where \f$n\f$ is the number of classes.
     */
    virtual double getBound(size_t i) const { return bounds_[i]; }
    
    /**
     * @return One midpoint.
     * @param i The index of the midpoint \f$\in [0,n[\f$ where \f$n\f$ is the number of classes.
     */
    virtual double getValue(size_t i) const { return midPoints_[i]; }

    /**
     * @return The nearest midpoint of a given value.
     * @param x The value to check.
     * @throw OutOfRangeException if the value is not in the domain interval.
     */
    virtual double getNearestValue(double x) const;
    
    /**
     * @return The index of the interval/class a given value belong to
     * @param x The value to check.
     * @throw OutOfRangeException if the value is not in the domain interval.
     */
    virtual size_t getIndex(double x) const;
    
};

#endif  //_DOMAIN_H_

