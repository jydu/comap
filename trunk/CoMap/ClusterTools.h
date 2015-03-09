//
// File: ClusterTools.h
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

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

// From Phylib:
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Likelihood/DRTreeLikelihood.h>
#include <Bpp/Phyl/Simulation/SequenceSimulator.h>
#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include <Bpp/Phyl/Mapping/SubstitutionCount.h>
#include <Bpp/Phyl/Mapping/ProbabilisticSubstitutionMapping.h>

using namespace bpp;

class Group
{
  protected:
    vector<string> names_;
    double height_;
    vector< vector<size_t> > subgroups_;
    vector<double> subgroupsHeights_;
    map<string, const Clonable*> properties_;

  public:
    Group() : names_(), height_(0), subgroups_(), subgroupsHeights_(), properties_() {}
    virtual ~Group() {}

  public:
    void setHeight(double height) { height_ = height; }
    const double getHeight() const { return  height_; }
    void addSubgroup(const vector<size_t>& subgroup, double height)
    {
      subgroups_.push_back(subgroup);
      subgroupsHeights_.push_back(height);
    }
    const vector< vector<size_t> > & getSubgroups() const
    {
      return subgroups_;
    }
    const vector<size_t> & getSubgroup(size_t i) const
    {
      return subgroups_[i];
    }
    const vector<string> getSubgroupNames(size_t i) const
    {
      vector<string> names;
      for(size_t j = 0; j < subgroups_[i].size(); j++)
      {
        names.push_back(names_[subgroups_[i][j]]);
      }
      return names;
    }
    double getSubgroupHeight(size_t i) const { return subgroupsHeights_[i]; }
    string toString() const
    {
      if (names_.size() == 0) return "[]";
      string text = "[" + names_[0];
      for (size_t i = 1; i < names_.size(); i++)
      {
        text += ";" + names_[i];
      }
      text += "]";
      return text;
    }
    // In case we are using indices instead of site names in tree:
    string toString(const vector<string>& tln) const
    {
      if (names_.size() == 0) return "[]";
      string text = "[" + tln[TextTools::to<size_t>(names_[0])];
      for(size_t i = 1; i < names_.size(); i++)
      {
        text += ";" + tln[TextTools::to<size_t>(names_[i])];
      }
      text += "]";
      return text;
    }
    string toString(size_t i) const
    {
      if (subgroups_[i].size() == 0) return "[]";
      string text = "[" + names_[subgroups_[i][0]];
      for(size_t j = 1; j < subgroups_[i].size(); j++)
      {
        text += ";" + names_[subgroups_[i][j]];
      }
      text += "]";
      return text;
    }

    size_t size() const { return names_.size(); }

    const string& operator[](size_t i) const { return names_[i]; }
    string& operator[](size_t i) { return names_[i]; }

    void add(const string& name) { names_.push_back(name); }

    const Clonable* getProperty(const string& name) { return properties_[name]; }

    void setProperty(const string& name, const Clonable* value) { properties_[name] = value; }

    vector<string> getPropertyNames() const { return MapTools::getKeys(properties_); }
};

class ClusterTools
{
  public:
    ClusterTools() {}
    virtual ~ClusterTools() {}

  public:

    static vector<Group> getGroups(const TreeTemplate<Node> * tree);

    static vector<const Node *> getSubtreesWithSize(const TreeTemplate<Node> & tree, size_t size);

    static vector<Group> getGroupsWithSize(const TreeTemplate<Node> & tree, size_t size);
    
    static void computeGlobalDistanceDistribution(
        DRTreeLikelihood & drtl,
        const SequenceSimulator & simulator,
	      SubstitutionCount & nijt,
        const Distance & distance,
        AgglomerativeDistanceMethod & clustering,
        size_t sizeOfDataSet,
        size_t nrep,
        size_t maxGroupSize,
        ofstream* out = 0);

    static void translate(TreeTemplate<Node> & tree, const vector<string> & tln)
    {
      translate(tree.getRootNode(), tln);
    }

    //Only min norm for now.
    static void computeNormProperties(TreeTemplate<Node> & tree, const ProbabilisticSubstitutionMapping & mapping);
    
  private:
    static Group getGroups(const Node & subtree, vector<Group> & groups);
    static Group getGroup(const Node & subtree);
    static size_t getSubtreesWithSize(const Node * subtree, size_t size, vector<const Node *> & subtrees);
    static void translate(Node * node, const vector<string> & tln)
    {
      if(node->isLeaf()) node->setName(tln[TextTools::to<size_t>(node->getName())]);
      for(size_t i = 0; i < node->getNumberOfSons(); i++) translate(node->getSon(i), tln);
    }
    static void computeNormProperties_(Node * node, const ProbabilisticSubstitutionMapping & mapping, double & minNorm);

};

#endif // CLUSTERTOOLS_H__

