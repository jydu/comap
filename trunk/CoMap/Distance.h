//
// File: Distance.h
// Created by: Julien Dutheil
// Created on: Tue Aug 30 15:58 2005
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

#ifndef _DISTANCE_H_
#define _DISTANCE_H_

#include "Statistics.h"

// From Utils:
#include <Utils/Number.h>

// From NumCalc:
#include <NumCalc/VectorExceptions.h>
#include <NumCalc/VectorTools.h>
using namespace VectorOperators;

// From PhylLib:
#include <Phyl/Node.h>
#include <Phyl/ProbabilisticSubstitutionMapping.h>

// From the STL:
#include <cmath>
using namespace std;

/**
 * @brief Compute a distance from two vectors.
 */
class Distance
{
	public:
		Distance() {}
		virtual ~Distance() {}

	public:
		virtual double getDistanceForPair(const Vdouble & v1, const Vdouble & v2) const throw (DimensionException) = 0;
		virtual double getDistanceForGroup(const vector<const Vdouble *> & v) const throw (DimensionException) = 0;
    virtual void setWeights(const Vdouble & weights) = 0;
    virtual void deleteWeights() = 0;
    virtual const Vdouble * getWeights() const = 0;
    virtual bool hasWeights() const = 0;
    virtual void setStatisticAsProperty(Node & node, const ProbabilisticSubstitutionMapping & mapping) const = 0;
};

class AbstractDistance : public Distance 
{
  protected:
    vector<double> * _weights;
    
	public:
		AbstractDistance(): _weights(NULL) {}
		virtual ~AbstractDistance() { if(_weights != NULL) delete _weights; }

	public:
    void setWeights(const Vdouble & weights)
    { 
      _weights = new Vdouble(weights);
      double s = VectorTools::sum(*_weights);
      *_weights /= s;
    }
    void deleteWeights() { delete _weights; _weights = NULL; }
    const Vdouble * getWeights() const { return _weights; }
    bool hasWeights() const { return _weights != NULL; }
    void setStatisticAsProperty(Node & node, const ProbabilisticSubstitutionMapping & mapping) const
    {
      _setStatisticAsProperty(node, mapping);
    }

  protected:
    //Use statistic from the tree:
    double _setStatisticAsProperty(Node & node, const ProbabilisticSubstitutionMapping & mapping) const
    {
      double height = 0;
      if(!node.isLeaf())
      {
        for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
        {
          Node * son = node.getSon(i);
          height = _setStatisticAsProperty(*son, mapping) + son->getDistanceToFather();
        }
        node.setNodeProperty("Stat", Number<double>(2*height));
      }
      return height;
    }
};

class AbstractMaximumDistance: public AbstractDistance
{
  public:
		double getDistanceForGroup(const vector<const Vdouble *> & v) const throw (DimensionException)
    {
      double maxi = log(0), val = 0;
      for(unsigned int i = 1; i < v.size(); i++)
      {
        for(unsigned int j = 0; j < i; j++)
        {
          val = getDistanceForPair(*v[i], *v[j]);
          if(val > maxi) maxi = val;
        }
      }
      return maxi;
    }
};

class EuclidianDistance: public AbstractMaximumDistance
{
	public:
		EuclidianDistance() {}
		virtual ~EuclidianDistance() {}

	public:
		double getDistanceForPair(const Vdouble & v1, const Vdouble & v2) const
			throw (DimensionException)
		{
			if(v1.size() != v2.size()) throw DimensionException("EuclidianDistance::d(...).", v2.size(), v1.size());
      if(_weights != NULL && _weights->size() != v1.size()) throw DimensionException("EuclidianDistance::d(...).", v1.size(), _weights->size());
			double d = 0;
			for(unsigned int i = 0; i < v1.size(); i++)
      {
				d += (_weights == NULL) ? std::pow(v2[i] - v1[i], 2) : (*_weights)[i] * std::pow(v2[i] - v1[i], 2) ;
			}
			return sqrt(d);
		}

};

/*
class CorrelationDistance: public AbstractDistance
{

	public:
		CorrelationDistance() {}
		virtual ~CorrelationDistance() {}

	public:
		double d(const vector<double> & v1, const vector<double> & v2) const
			throw (DimensionException)
		{
      return (_weights == NULL) ? 1. - getCorrelation(v1, v2) : 1. - getWeightedCorrelation(v1, v2);
		}

  protected:
    double getCorrelation(const vector<double> & v1, const vector<double> & v2) const
			throw (DimensionException)
    {
			if(v1.size() != v2.size()) throw DimensionException("CorrelationDistance::getCorrelation(...).", v2.size(), v1.size());
			double sumx = 0., sumy = 0., sumx2 = 0., sumy2 = 0., sumxy = 0.;
      for(unsigned int i = 0; i < v1.size(); i++)
      {
			  sumx += v1[i];
			  sumy += v2[i];
			  sumxy += v1[i]*v2[i];
			  sumx2 += std::pow(v1[i], 2.);
			  sumy2 += std::pow(v2[i], 2.);
			}
			double n = (double)v1.size();
			double mx = sumx/n;
			double my = sumy/n;
			return (sumxy/n - mx * my)
			  / (
			     sqrt(sumx2/n - std::pow(mx, 2.))
				  *sqrt(sumy2/n - std::pow(my, 2.))
			  );
    }
    
    double getWeightedCorrelation(const vector<double> & v1, const vector<double> & v2) const
			throw (DimensionException)
    {
			if(v1.size() != v2.size()) throw DimensionException("CorrelationDistance::getWeightedCorrelation(...).", v2.size(), v1.size());
			double wsumx = 0., wsumy = 0., wsumx2 = 0., wsumy2 = 0., wsumxy = 0.;
      for(unsigned int i = 0; i < v1.size(); i++)
      {
			  wsumx += (*_weights)[i] * v1[i];
			  wsumy += (*_weights)[i] * v2[i];
			  wsumxy += (*_weights)[i] * v1[i] * v2[i];
			  wsumx2 += (*_weights)[i] * std::pow(v1[i], 2.);
			  wsumy2 += (*_weights)[i] * std::pow(v2[i], 2.);
			}
			return (wsumxy - wsumx * wsumy)
			  / (
			     sqrt(wsumx2 - std::pow(wsumx, 2.))
				  *sqrt(wsumy2 - std::pow(wsumy, 2.))
			  );
    }    

};

class SquareCorrelationDistance: public CorrelationDistance
{

	public:
		SquareCorrelationDistance() {}
		virtual ~SquareCorrelationDistance() {}

	public:
		double d(const vector<double> & v1, const vector<double> & v2) const
			throw (DimensionException)
		{
      return (_weights == NULL) ?
        1. - std::pow(getCorrelation(v1, v2), 2.) :
        1. - std::pow(getWeightedCorrelation(v1, v2), 2.);
    }

};

class CompensationDistance: public AbstractDistance
{

	public:
		CompensationDistance() {}
		virtual ~CompensationDistance() {}

	public:
		double d(const vector<double> & v1, const vector<double> & v2) const
			throw (DimensionException)
		{
      return (_weights == NULL) ? getCompensation(v1, v2) : getWeightedCompensation(v1, v2);
		}

  protected:
    double getCompensation(const vector<double> & v1, const vector<double> & v2) const
			throw (DimensionException)
    {
			if(v1.size() != v2.size()) throw DimensionException("CompensationDistance::getCompensation(...).", v2.size(), v1.size());
      double sumsq1 = 0., sumsq2 = 0., sumsq3 = 0.;
      for(unsigned int i = 0; i < v1.size(); i++)
      {
        sumsq1 += pow(v1[i],2);
        sumsq2 += pow(v2[i],2);
        sumsq3 += pow(v1[i]+v2[i],2);
      }
      return sqrt(sumsq3) / (sqrt(sumsq1) + sqrt(sumsq2));
      
			//vector<double> v3 = v1 + v2;
      //double sumsq = 0.;
      //for(unsigned int i = 0; i < v1.size(); i++)
      //{
			//  sumsq += pow(v1[i] + v2[i], 2.);
			//}
			//double n = (double)v1.size();
			//return sumsq / n;
    }
    
    double getWeightedCompensation(const vector<double> & v1, const vector<double> & v2) const
			throw (DimensionException)
    {
			if(v1.size() != v2.size()) throw DimensionException("CompensationDistance::getWeightedCompensation(...).", v2.size(), v1.size());
      double wsumsq1 = 0., wsumsq2 = 0., wsumsq3 = 0.;
      for(unsigned int i = 0; i < v1.size(); i++)
      {
        wsumsq1 += pow(v1[i],2)*(*_weights)[i];
        wsumsq2 += pow(v2[i],2)*(*_weights)[i];
        wsumsq3 += pow(v1[i]+v2[i],2)*(*_weights)[i];
      }
      return sqrt(wsumsq3) / (sqrt(wsumsq1) + sqrt(wsumsq2));

      //vector<double> v3 = v1 + v2;
      //double wsumsq = 0.;
      //for(unsigned int i = 0; i < v1.size(); i++)
      //{
			//  wsumsq += (*_weights)[i] * pow(v1[i] + v2[i], 2.);
			//}
			//return wsumsq;
    }

};*/

class StatisticBasedDistance: public Distance
{
  protected:
    Statistic * _stat;
    double _comp;

  public:
    StatisticBasedDistance(Statistic * stat, double comp = 0.): _stat(stat), _comp(comp) {}
    virtual ~StatisticBasedDistance() {}

  public:
		double getDistanceForPair(const Vdouble & v1, const Vdouble & v2) const throw (DimensionException)
    {
      return _comp - _stat->getValueForPair(v1, v2);
    }
		double getDistanceForGroup(const vector<const Vdouble *> & v) const throw (DimensionException)
    {
      return _comp - _stat->getValueForGroup(v);
    }
    void setWeights(const vector<double> & weights) { _stat->setWeights(weights); }
    void deleteWeights() { _stat->deleteWeights(); }
    const vector<double> * getWeights() const { return _stat->getWeights(); }
    bool hasWeights() const { return _stat->hasWeights(); }
    void setStatisticAsProperty(Node & node, const ProbabilisticSubstitutionMapping & mapping) const
    {
      _setStatisticAsProperty(node, mapping);
    }

    Statistic * getStatistic() { return _stat; }
    const Statistic * getStatistic() const { return _stat; }
    
  protected:
    double _setStatisticAsProperty(Node & node, const ProbabilisticSubstitutionMapping & mapping) const
    {
      double height = 0;
      if(!node.isLeaf())
      {
        for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
        {
          Node * son = node.getSon(i);
          height = _setStatisticAsProperty(*son, mapping) + son->getDistanceToFather();
        }
        node.setNodeProperty("Stat", Number<double>(_comp - 2*height));
      }
      return height;
    }

};

class CompensationDistance: public Distance
{
  protected:
    CompensationStatistic _stat;

  public:
    CompensationDistance() {}
    virtual ~CompensationDistance() {}

  public:
		double getDistanceForPair(const Vdouble & v1, const Vdouble & v2) const throw (DimensionException)
    {
      return 1. - _stat.getValueForPair(v1, v2);
    }
		double getDistanceForGroup(const vector<const Vdouble *> & v) const throw (DimensionException)
    {
      return 1. - _stat.getValueForGroup(v);
    }
    void setStatisticAsProperty(Node & node, const ProbabilisticSubstitutionMapping & mapping) const
    {
      Vdouble sigma(mapping.getNumberOfBranches(), 0.);
      double sumNorms = 0;
      _setStatisticAsProperty(node, mapping, sigma, sumNorms);
    }

    void setWeights(const vector<double> & weights) { _stat.setWeights(weights); }
    void deleteWeights() { _stat.deleteWeights(); }
    const vector<double> * getWeights() const { return _stat.getWeights(); }
    virtual bool hasWeights() const { return _stat.hasWeights(); }

  protected:
    void _setStatisticAsProperty(Node & node, const ProbabilisticSubstitutionMapping & mapping, Vdouble & sigma, double & sumNorms) const
    {
      if(node.isLeaf())
      {
        sigma = mapping[TextTools::to<unsigned int>(node.getName())];
        sumNorms = VectorTools::norm<double,double>(mapping[TextTools::to<unsigned int>(node.getName())]);
      }
      else
      {
        for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
        {
          Vdouble sigmaGroup(mapping.getNumberOfBranches(), 0.);
          double sumNormsGroup = 0;
          _setStatisticAsProperty(*node.getSon(i), mapping, sigmaGroup, sumNormsGroup);
          sigma += sigmaGroup;
          sumNorms += sumNormsGroup;
        }
        node.setNodeProperty("Stat", Number<double>(1. - VectorTools::norm<double,double>(sigma) / sumNorms));
      }
    }

};

#endif //_DISTANCE_H_

