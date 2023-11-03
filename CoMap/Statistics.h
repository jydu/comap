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

#include "Domain.h"
#include <Bpp/Numeric/VectorTools.h>

using namespace bpp;
using namespace std;

/**
 * @brief Coevolution measure.
 *
 * A coevolution statistic can be computed for a pair of sites, or for a group of higher size.
 * It is computed from substitution vectors, passed as input. The vectors have a second dimension allowing detailed counts if > 1.
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
     * @param v1 The first vector (with putative detailed mutation counts).
     * @param v2 The second vector (with putative detailed mutation counts).
     * @return The value of the statistic for this pair.
     * @throw DimensionException If the two vectors do not have the same dimension or do not match the dimension of the weights if set up. 
     */
    virtual double getValueForPair(const VVdouble& v1, const VVdouble& v2) const = 0;
    /**
     * @brief Get the value of the statistic for a set of vectors.
     *
     * @param v A set of vectors.
     * @return The value of the statistic for this group.
     * @throw DimensionException If the two vectors do not have the same dimension or do not match the dimension of the weights if set up. 
     */
    virtual double getValueForGroup(const vector<const VVdouble*>& v) const = 0;

    /**
     * @brief Set weights for this statistic.
     *
     * @param w A vector of weights.
     */
    virtual void setWeights(const Vdouble& w) = 0;
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
    virtual const Vdouble* getWeights() const = 0;
};

class AbstractMinimumStatistic: public Statistic
{
  protected:
    Vdouble* weight_s;

  public:
    AbstractMinimumStatistic(): weight_s(0) {}
    AbstractMinimumStatistic(const AbstractMinimumStatistic& stat) : weight_s(new Vdouble(*stat.weight_s)) {}
    AbstractMinimumStatistic& operator=(const AbstractMinimumStatistic& stat) {
      weight_s = new Vdouble(*stat.weight_s);
      return *this;
    }
    virtual ~AbstractMinimumStatistic() { if (weight_s) delete weight_s; }

  public:
    virtual double getValueForGroup(const vector<const VVdouble*>& v) const override
    {
      double mini = -log(0), val = 0.;
      for (size_t i = 1; i < v.size(); ++i)
      {
        for (size_t j = 0; j < i; ++j)
        {
          val = getValueForPair(*v[i], *v[j]);
          if (val < mini) mini = val;
        }
      }
      return mini;
    }

    void setWeights(const Vdouble& w) override
    { 
      if(weight_s) delete weight_s;
      weight_s = new Vdouble(w);
      *weight_s /= VectorTools::sum(w);
    }

    void deleteWeights() override
    {
      if (weight_s) delete weight_s;
      weight_s = 0;
    }

    bool hasWeights() const override
    {
      return weight_s != 0;
    }
    const Vdouble* getWeights() const override { return weight_s; }

    static Vdouble getUD(const VVdouble& v, unsigned int type) {
      Vdouble u(v.size());
      for (size_t i = 0; i < v.size(); ++i) {
        u[i] = v[i][type];
      }
      return u;
    }

};

class CorrelationStatistic: public AbstractMinimumStatistic
{
  public:
    double getValueForPair(const VVdouble& v1, const VVdouble& v2) const
    {
      if (weight_s)
        return VectorTools::cor<double, double>(getUD(v1, 0), getUD(v2, 0), *weight_s, false);
      else 
        return VectorTools::cor<double, double>(getUD(v1, 0), getUD(v2, 0));
    }
};

class CorrectedCorrelationStatistic: public AbstractMinimumStatistic
{
  private:
    Vdouble meanVector1_;
    Vdouble meanVector2_;

  public:
    CorrectedCorrelationStatistic(const Vdouble& meanVector): meanVector1_(meanVector), meanVector2_(meanVector) {}
    CorrectedCorrelationStatistic(const Vdouble& meanVector1, const Vdouble& meanVector2): meanVector1_(meanVector1), meanVector2_(meanVector2) {}
    CorrectedCorrelationStatistic(): meanVector1_(), meanVector2_() {}

  public:
    double getValueForPair(const VVdouble& v1, const VVdouble& v2) const
    {
      if (weight_s)
        return VectorTools::cor<double, double>(getUD(v1, 0) - meanVector1_, getUD(v2, 0) - meanVector2_, *weight_s, false);
      else 
        return VectorTools::cor<double, double>(getUD(v1, 0) - meanVector1_, getUD(v2, 0) - meanVector2_);
    }

    void setMeanVector(const Vdouble& meanVector) {
      meanVector1_ = meanVector;
      meanVector2_ = meanVector;
    }
    void setMeanVectors(const Vdouble& meanVector1, const Vdouble& meanVector2) {
      meanVector1_ = meanVector1;
      meanVector2_ = meanVector2;
    }
};

class CovarianceStatistic: public AbstractMinimumStatistic
{
  public:
    double getValueForPair(const VVdouble& v1, const VVdouble& v2) const
    {
      if (weight_s)
        return VectorTools::cov<double, double>(getUD(v1, 0), getUD(v2, 0), *weight_s, false, false);
      else
        return VectorTools::cov<double, double>(getUD(v1, 0), getUD(v2, 0));
    }
};

class CosinusStatistic: public AbstractMinimumStatistic
{
  public:
    double getValueForPair(const VVdouble& v1, const VVdouble& v2) const
    {
      if(weight_s)
        return VectorTools::cos<double, double>(getUD(v1, 0), getUD(v2, 0), *weight_s);
      else
        return VectorTools::cos<double, double>(getUD(v1, 0), getUD(v2, 0));
    }
};

class CosubstitutionNumberStatistic: public AbstractMinimumStatistic
{
  public: 
    double getValueForPair(const VVdouble& v1, const VVdouble& v2) const
    {
      if (v1.size() != v2.size())
        throw new Exception("CosubstitutionNumberStatistic::getValueForPair: vectors must have the same size.");
      size_t n = v1.size();
      double c = 0;
      for (size_t i = 0; i < n; i++)
      {
        if (VectorTools::sum(v1[i]) >= 1. && VectorTools::sum(v2[i]) >= 1.) c++;
      }
      return c;
    }
};

class CompensationStatistic: public AbstractMinimumStatistic
{
  public:
    double getValueForPair(const VVdouble& v1, const VVdouble& v2) const override
    {
       if (v1.size() != v2.size())
        throw DimensionException("CompensationStatistic::getValueForPair.", v2.size(), v1.size());
      double sumsq1 = 0., sumsq2 = 0., sumsq3 = 0., w = 1.;
      for (size_t i = 0; i < v1.size(); i++)
      {
        double sv1 = VectorTools::sum(v1[i]);
        double sv2 = VectorTools::sum(v2[i]);
        if (weight_s) w = (*weight_s)[i];
        sumsq1 += pow(sv1, 2) * w;
        sumsq2 += pow(sv2, 2) * w;
        sumsq3 += pow(sv1 + sv2, 2) * w;
      }
      return 1. - sqrt(sumsq3) / (sqrt(sumsq1) + sqrt(sumsq2));
    }

    double getValueForGroup(const vector<const VVdouble*>& v) const override
    {
      for (size_t j = 1; j < v.size(); j++)
      {
         if (v[0]->size() != v[j]->size())
          throw DimensionException("CompensationStatistic::getValueForGroup.", v[0]->size(), v[j]->size());
      }
      vector<double> sumsq1(v.size(), 0);
      double s = 0., sumsq2 = 0., w = 1.;
      for (size_t i = 0; i < v[0]->size(); i++)
      {
        if (weight_s) w = (*weight_s)[i];
        s = 0.;
        for (size_t j = 0; j < v.size(); j++)
        {
          double sv = VectorTools::sum((*v[j])[i]);
          sumsq1[j] += pow(sv, 2) * w;
          s += sv;
        }
        sumsq2 += pow(s, 2) * w;
      }
      double sumnorms = 0.;
      for (size_t j = 0; j < v.size(); j++)
      {
        sumnorms += sqrt(sumsq1[j]);
      }
      return 1. - sqrt(sumsq2) / (sumnorms);
    }
};

class MutualInformationStatistic: public AbstractMinimumStatistic
{
  public:
    double getValueForPair(const VVdouble& v1, const VVdouble& v2) const
    {
      //Forn now weights are ignored
      return VectorTools::miContinuous<double, double>(getUD(v1, 0), getUD(v2, 0));
    }
};

class DiscreteMutualInformationStatistic: public AbstractMinimumStatistic
{
  private:
    Domain domain_;

  public:
    DiscreteMutualInformationStatistic(const Vdouble& bounds) :
      AbstractMinimumStatistic(), domain_(bounds) {}

  public:
    double getValueForPair(const VVdouble& v1, const VVdouble& v2) const
    {
      vector<size_t> c1(v1.size());
      vector<size_t> c2(v2.size());
      for (size_t i = 0; i < c1.size(); i++)
      {
        c1[i] = domain_.getIndex(VectorTools::sum(v1[i]));
        c2[i] = domain_.getIndex(VectorTools::sum(v2[i]));
      }
      return VectorTools::miDiscrete<size_t, double>(c1, c2); //Weights ignored for now.
    }
    //NB: also define a group statistic if it works...
};


#endif  //_STATISTICS_H_

