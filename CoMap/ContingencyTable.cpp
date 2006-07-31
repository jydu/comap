//
// File: ContingencyTable.cpp
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

#include "ContingencyTable.h"

// From the STL:
#include <cmath>

// From Utils:
#include <Utils/TextTools.h>

void ContingencyTable::addToCount(unsigned int groupSize, double norm, double dmax) throw (Exception)
{
  try {
    unsigned int i = getGroupSizeIndex(groupSize);
    unsigned int j = getNormIndex(norm);
    unsigned int k = getDMaxIndex(dmax);
    _table[i][j][k]++;
  } catch (Exception & ex) {
    _count++;
    throw ex;
  }
  _count++;
}

unsigned int ContingencyTable::getCount(unsigned int groupSize, double norm, double dmax) const throw (Exception)
{
  unsigned int i = getGroupSizeIndex(groupSize);
  unsigned int j = getNormIndex(norm);
  unsigned int k = getDMaxIndex(dmax);
  return _table[i][j][k];
}

double ContingencyTable::getDensity(unsigned int groupSize, double norm, double dmax) const throw (Exception)
{
  unsigned int i = getGroupSizeIndex(groupSize);
  unsigned int j = getNormIndex(norm);
  unsigned int k = getDMaxIndex(dmax);
  return (double)_table[i][j][k] / (double)_count;
}

double ContingencyTable::getPValue(unsigned int groupSize, double norm, double dmax) const throw (Exception)
{
  unsigned int i = getGroupSizeIndex(groupSize);
  unsigned int j = getNormIndex(norm);
  unsigned int k = getDMaxIndex(dmax);
  return std::min(1., (double)(_table[i][j][k]+1) / (double)_count);
}

unsigned int ContingencyTable::getGroupSizeIndex(unsigned int groupSize) const throw (Exception)
{
  if(groupSize > _groupSizeMax) throw Exception("Group size " + TextTools::toString(groupSize) + " is out of table range.");
  return groupSize - 2; // Smaller group has size = 2 and is at index 0.
}

unsigned int ContingencyTable::getNormIndex(double norm) const throw (Exception)
{
  unsigned int n = (unsigned int)floor(norm / _normUnit);
  if(n >= _normClassNumber) throw Exception("Norm " + TextTools::toString(norm) + " is out of table range.");
  return n; // n ranges from 0 to _normClassNumber - 1.
}
    
unsigned int ContingencyTable::getDMaxIndex(double dmax) const throw (Exception)
{
  unsigned int n = (unsigned int)floor(dmax / _dmaxUnit);
  if(n >= _dmaxClassNumber) throw Exception("Distance max " + TextTools::toString(dmax) + " is out of table range.");
  return n; // n ranges from 0 to _dmaxClassNumber - 1.
}

