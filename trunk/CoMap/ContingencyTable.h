//
// File: ContingencyTable.h
// Created by: Julien Dutheil
// Created on: Fri Mar 17 2006
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

#ifndef _CONTINGENCYTABLE_H_
#define _CONTINGENCYTABLE_H_

// From the STL:
#include <vector>
using namespace std;

// From Utils:
#include <Utils/exceptions>

class ContingencyTable
{
  protected:
    unsigned int _groupSizeMax;
    double _normUnit;
    double _dmaxUnit;
    unsigned int _normClassNumber;
    unsigned int _dmaxClassNumber;
    vector< vector< vector<unsigned int> > > _table;
    unsigned int _count;
    
  public:
    ContingencyTable(
        unsigned int groupSizeMax,
        double normUnit,
        unsigned int normClassNumber,
        double dmaxUnit,
        unsigned int dmaxClassNumber)
    {
      _groupSizeMax = groupSizeMax;
      _normUnit = normUnit;
      _normClassNumber = normClassNumber;
      _dmaxUnit = dmaxUnit;
      _dmaxClassNumber = dmaxClassNumber;
      _table.resize(_groupSizeMax);
      for(unsigned int i = 0; i < _groupSizeMax; i++)
      {
        _table[i].resize(_normClassNumber);
        for(unsigned int j = 0; j < _normClassNumber; j++)
        {
          _table[i][j].resize(_dmaxClassNumber);
        }
      }
    }

    virtual ~ContingencyTable() {}

  public:
    void addToCount(unsigned int groupSize, double norm, double dmax) throw (Exception);
    unsigned int getCount(unsigned int groupSize, double norm, double dmax) const throw (Exception);
    double getDensity(unsigned int groupSize, double norm, double dmax) const throw (Exception);
    double getPValue(unsigned int groupSize, double norm, double dmax) const throw (Exception);
    unsigned int getTotalCount() const { return _count; }
    unsigned int getGroupSizeIndex(unsigned int groupSize) const throw (Exception);
    unsigned int getNormIndex(double norm) const throw (Exception);
    unsigned int getDMaxIndex(double dmax) const throw (Exception);
};

#endif // _CONTINGENCYTABLE_H_
    
