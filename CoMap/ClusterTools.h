//
// File: TestTools.h
// Created by: Julien Dutheil
// Created on: Thu Mar 2 2006
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

#ifndef _CLUSTERTOOLS_H_
#define _CLUSTERTOOLS_H_

#include "Distance.h"

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>

// From Phylib:
#include <Phyl/TreeTemplate.h>
#include <Phyl/HomogeneousSequenceSimulator.h>
#include <Phyl/AgglomerativeDistanceMethod.h>
#include <Phyl/SubstitutionCount.h>
#include <Phyl/ProbabilisticSubstitutionMapping.h>

class Group: public vector<string>
{
  protected:
    double _height;
    vector< vector<unsigned int> > _subgroups;
    vector<double> _subgroupsHeights;
    map<string, const Clonable *> _properties;

  public:
    Group() { _height = 0; }
    virtual ~Group() {}

  public:
    void setHeight(double height) { _height = height; }
    const double getHeight() const { return  _height; }
    void addSubgroup(const vector<unsigned int> & subgroup, double height)
    {
      _subgroups.push_back(subgroup);
      _subgroupsHeights.push_back(height);
    }
    const vector< vector<unsigned int> > & getSubgroups() const
    {
      return _subgroups;
    }
    const vector<unsigned int> & getSubgroup(unsigned int i) const
    {
      return _subgroups[i];
    }
    const vector<string> getSubgroupNames(unsigned int i) const
    {
      vector<string> names;
      for(unsigned int j = 0; j < _subgroups[i].size(); j++)
      {
        names.push_back((*this)[_subgroups[i][j]]);
      }
      return names;
    }
    double getSubgroupHeight(unsigned int i) const { return _subgroupsHeights[i]; }
    string toString() const
    {
      if(size() == 0) return "[]";
      string text = "[" + at(0);
      for(unsigned int i = 1; i < size(); i++)
      {
        text += ";" + at(i);
      }
      text += "]";
      return text;
    }
    // In case we are using indices instead of site names in tree:
    string toString(const vector<string> & tln) const
    {
      if(size() == 0) return "[]";
      string text = "[" + tln[TextTools::to<unsigned int>(at(0))];
      for(unsigned int i = 1; i < size(); i++)
      {
        text += ";" + tln[TextTools::to<unsigned int>(at(i))];
      }
      text += "]";
      return text;
    }
    string toString(unsigned int i) const
    {
      if(_subgroups[i].size() == 0) return "[]";
      string text = "[" + at(_subgroups[i][0]);
      for(unsigned int j = 1; j < _subgroups[i].size(); j++)
      {
        text += ";" + at(_subgroups[i][j]);
      }
      text += "]";
      return text;
    }

    const Clonable * getProperty(const string & name) { return _properties[name]; }

    void setProperty(const string & name, const Clonable * value) { _properties[name] = value; }

    vector<string> getPropertyNames() const { return MapTools::getKeys(_properties); }
};

class ClusterTools
{
  public:
    ClusterTools() {}
    virtual ~ClusterTools() {}

  public:

    static vector<Group> getGroups(const TreeTemplate<Node> * tree);

    static vector<const Node *> getSubtreesWithSize(const TreeTemplate<Node> & tree, unsigned int size);

    static vector<Group> getGroupsWithSize(const TreeTemplate<Node> & tree, unsigned int size);
    
    static void computeGlobalDistanceDistribution(
        const HomogeneousSequenceSimulator & simulator,
	      const SubstitutionCount & nijt,
        const Distance & distance,
        AgglomerativeDistanceMethod & clustering,
        const vector<double> & scales,
        unsigned int sizeOfDataSet,
        unsigned int nrep,
        unsigned int maxGroupSize,
        ofstream * out = NULL);

    static void translate(TreeTemplate<Node> & tree, const vector<string> & tln)
    {
      translate(tree.getRootNode(), tln);
    }

    //Only min norm for now.
    static void computeNormProperties(TreeTemplate<Node> & tree, const ProbabilisticSubstitutionMapping & mapping);
    
  private:
    static Group getGroups(const Node & subtree, vector<Group> & groups);
    static Group getGroup(const Node & subtree);
    static unsigned int getSubtreesWithSize(const Node * subtree, unsigned int size, vector<const Node *> & subtrees);
    static void translate(Node * node, const vector<string> & tln)
    {
      if(node->isLeaf()) node->setName(tln[TextTools::to<unsigned int>(node->getName())]);
      for(unsigned int i = 0; i < node->getNumberOfSons(); i++) translate(node->getSon(i), tln);
    }
    static void _computeNormProperties(Node * node, const ProbabilisticSubstitutionMapping & mapping, double & minNorm);

};

#endif // _CLUSTERTOOLS_H_

