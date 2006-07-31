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
using namespace VectorStatTools;
using namespace VectorFunctions;

class Statistic
{
	public:
		Statistic() {}
		virtual ~Statistic() {}

	public:
		virtual double getValue(const Vdouble & v1, const Vdouble & v2) const = 0;
};

class CorrelationStatistic: public Statistic
{
	public:
		double getValue(const Vdouble & v1, const Vdouble & v2) const {
			return cor(v1, v2);
		}
};

class CovarianceStatistic: public Statistic
{
	public:
		double getValue(const Vdouble & v1, const Vdouble & v2) const {
			return cov(v1, v2);
		}
};

class CosinusStatistic: public Statistic
{
	public:
		double getValue(const Vdouble & v1, const Vdouble & v2) const {
			return cos(v1, v2);
		}
};

class CosubstitutionNumberStatistic: public Statistic
{
	public: 
		double getValue(const Vdouble & v1, const Vdouble & v2) const {
			if(v1.size() != v2.size()) throw new Exception("CosubstitutionNumberStatistic::getValue: vectors must have the same size.");
			unsigned int n = v1.size();
			double c = 0;
			for(unsigned int i = 0; i < n; i++) {
				if(v1[i] >= 1. && v2[i] >= 1.) c++;
			}
			return c;
		}
};

#endif	//_STATISTICS_H_

