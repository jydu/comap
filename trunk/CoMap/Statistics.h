//
// File: Statistics.h
// Created by: Julien Dutheil
// Created on: Mon Feb  2 14:50:40 2004
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004, 2005, 2006)

This software is a computer program whose purpose is to map substitutions
on a tree and to detect co-evolving positions in a dataset.

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


#ifndef _STATISTICS_H_
#define _STATISTICS_H_

// From NumCalc:
#include <NumCalc/VectorTools.h>

using namespace bpp;

/**
 * @brief Coevolution measure.
 *
 * A coevolution statistic can be computed for a pair of sites, or for a group of higher size.
 * It is computed from substitution vectors, passed as input.
 * A vector of weights (one weight per branch) can be specified, depending on the measure and its implementation.
 */
class Statistic
{
	public:
		Statistic() {}
		virtual ~Statistic() {}

	public:
    /**
     * @brief Get the value of the statistic for a pair of vectors.
     *
     * @param v1 The first vector.
     * @param v2 The second vector.
     * @return The value of the statistic for this pair.
     * @throw DimensionException If the two vectors do not have the same dimension or do not match the dimension of the weights if set up. 
     */
		virtual double getValueForPair(const Vdouble & v1, const Vdouble & v2) const throw (DimensionException) = 0;
    /**
     * @brief Get the value of the statistic for a set of vectors.
     *
     * @param v A set of vectors.
     * @return The value of the statistic for this group.
     * @throw DimensionException If the two vectors do not have the same dimension or do not match the dimension of the weights if set up. 
     */
    virtual double getValueForGroup(const vector<const Vdouble *> & v) const throw (DimensionException) = 0;

    /**
     * @brief Set weights for this statistic.
     *
     * @param w A vector of weights.
     */
    virtual void setWeights(const Vdouble & w) = 0;
    /**
     * @brief Delete the weights associated to this statistic, if any.
     */
    virtual void deleteWeights() = 0;
    /**
     * @brief Tell if weights are associated to this statistic.
     *
     * @return y/n.
     */
    virtual bool hasWeights() const = 0;
    /**
     * @brief Get the weights associated to this statistic, if any.
     *
     * @return A pointer toward weights, NULL if no weights are attached to this statistic.
     */
    virtual const Vdouble * getWeights() const = 0;
};

class AbstractMinimumStatistic: public Statistic
{
  protected:
    Vdouble * _weights;

  public:
    AbstractMinimumStatistic(): _weights(NULL) {}
    virtual ~AbstractMinimumStatistic() { if(_weights) delete _weights; }

  public:
    virtual double getValueForGroup(const vector<const Vdouble *> & v) const throw (DimensionException)
    {
      double mini = -log(0), val = 0.;
      for(unsigned int i = 1; i < v.size(); i++)
      {
        for(unsigned int j = 0; j < i; j++)
        {
          val = getValueForPair(*v[i], *v[j]);
          if(val < mini) mini = val;
        }
      }
      return mini;
    }
    void setWeights(const Vdouble & w)
    { 
      if(_weights) delete _weights;
      _weights = new Vdouble(w);
      *_weights /= VectorTools::sum(w);
    }
    void deleteWeights()
    {
      if(_weights) delete _weights;
      _weights = NULL;
    }
    bool hasWeights() const
    {
      return _weights != NULL;
    }
    const Vdouble * getWeights() const { return _weights; }

};

class CorrelationStatistic: public AbstractMinimumStatistic
{
	public:
		double getValueForPair(const Vdouble & v1, const Vdouble & v2) const throw (DimensionException)
    {
      if(_weights)
        return VectorTools::cor<double, double>(v1, v2, *_weights, false);
      else 
			  return VectorTools::cor<double, double>(v1, v2);
		}
};

class CovarianceStatistic: public AbstractMinimumStatistic
{
	public:
		double getValueForPair(const Vdouble & v1, const Vdouble & v2) const throw (DimensionException)
    {
      if(_weights)
        return VectorTools::cov<double, double>(v1, v2, *_weights, false, false);
      else
			  return VectorTools::cov<double, double>(v1, v2);
		}
};

class CosinusStatistic: public AbstractMinimumStatistic
{
	public:
		double getValueForPair(const Vdouble & v1, const Vdouble & v2) const throw (DimensionException)
    {
      if(_weights)
        return VectorTools::cos<double, double>(v1, v2, *_weights);
      else
			  return VectorTools::cos<double, double>(v1, v2);
		}
};

class CosubstitutionNumberStatistic: public AbstractMinimumStatistic
{
	public: 
		double getValueForPair(const Vdouble & v1, const Vdouble & v2) const throw (DimensionException)
    {
			if(v1.size() != v2.size()) throw new Exception("CosubstitutionNumberStatistic::getValueForPair: vectors must have the same size.");
			unsigned int n = v1.size();
			double c = 0;
			for(unsigned int i = 0; i < n; i++)
      {
				if(v1[i] >= 1. && v2[i] >= 1.) c++;
			}
			return c;
		}
};

class CompensationStatistic: public AbstractMinimumStatistic
{
  public:
    double getValueForPair(const Vdouble & v1, const Vdouble & v2) const throw (DimensionException)
    {
 		  if(v1.size() != v2.size()) throw DimensionException("CompensationStatistic::getValueForPair.", v2.size(), v1.size());
      double sumsq1 = 0., sumsq2 = 0., sumsq3 = 0., w = 1.;
      for(unsigned int i = 0; i < v1.size(); i++)
      {
        if(_weights) w = (*_weights)[i];
        sumsq1 += pow(v1[i], 2) * w;
        sumsq2 += pow(v2[i], 2) * w;
        sumsq3 += pow(v1[i] + v2[i], 2) * w;
      }
      return 1. - sqrt(sumsq3) / (sqrt(sumsq1) + sqrt(sumsq2));
    }
    double getValueForGroup(const vector<const Vdouble *> & v) const throw (DimensionException)
    {
      for(unsigned int j = 1; j < v.size(); j++)
      {
 		    if(v[0]->size() != v[j]->size()) throw DimensionException("CompensationStatistic::getValueForGroup.", v[0]->size(), v[j]->size());
      }
      vector<double> sumsq1(v.size(), 0);
      double s = 0., sumsq2 = 0., w = 1.;
      for(unsigned int i = 0; i < v[0]->size(); i++)
      {
        if(_weights) w = (*_weights)[i];
        s = 0.;
        for(unsigned int j = 0; j < v.size(); j++)
        {
          sumsq1[j] += pow((*v[j])[i], 2) * w;
          s += (*v[j])[i];
        }
        sumsq2 += pow(s, 2) * w;
      }
      double sumnorms = 0.;
      for(unsigned int j = 0; j < v.size(); j++)
      {
        sumnorms += sqrt(sumsq1[j]);
      }
      return 1. - sqrt(sumsq2) / (sumnorms);
    }
};

#endif	//_STATISTICS_H_

