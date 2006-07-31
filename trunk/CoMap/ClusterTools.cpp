//
// File: TestTools.cpp
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

#include "ClusterTools.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>
using namespace VectorOperators;
using namespace VectorFunctions;

// From SeqLib:
#include <Seq/VectorSiteContainer.h>

// From PhylLib:
#include <Phyl/DRHomogeneousTreeLikelihood.h>
#include <Phyl/SubstitutionMappingTools.h>

#include <iostream>
#include <fstream>
using namespace std;


vector<Group> ClusterTools::getGroups(const TreeTemplate<Node> * tree)
{
  vector<Group> groups;
  getGroups(*tree->getRootNode(), groups);
  return groups;
}
  
Group ClusterTools::getGroups(const Node & subtree, vector<Group> & groups)
{
  Group group;
  if(subtree.isLeaf()) {
    group.push_back(subtree.getName());
  } else {
    unsigned int subcount = 0;
    vector<unsigned int> thisgroup;
    for(unsigned int i = 0; i < subtree.getNumberOfSons(); i++)
    {
      const Node * son = subtree.getSon(i);
      // Get group for each son node:
      Group sonGroup = getGroups(*son, groups); 
      vector<unsigned int> subgroup(sonGroup.size());
      unsigned int index = subcount;
      for(unsigned int j = 0; j < sonGroup.size(); j++)
      {
        group.push_back(sonGroup[j]);
        subgroup[j] = subcount;
        thisgroup.push_back(subcount);
        subcount++;
      }
      if(subgroup.size() > 1) {
        //group.addSubgroup(subgroup, sonGroup.getHeight());
        // Recursively add subgroups:
        for(unsigned int j = 0; j < sonGroup.getSubgroups().size(); j++)
        {
          group.addSubgroup(sonGroup.getSubgroup(j) + index, sonGroup.getSubgroupHeight(j));
        }
      }
      group.setHeight(sonGroup.getHeight() + son->getDistanceToFather());
    }
    group.addSubgroup(thisgroup, group.getHeight());
    
    // Add this group:
    groups.push_back(group);
  }
  return group;
}

Group ClusterTools::getGroup(const Node & subtree)
{
  Group group;
  if(subtree.isLeaf()) {
    group.push_back(subtree.getName());
  } else {
    unsigned int subcount = 0;
    vector<unsigned int> thisgroup;
    for(unsigned int i = 0; i < subtree.getNumberOfSons(); i++)
    {
      const Node * son = subtree.getSon(i);
      // Get group for each son node:
      Group sonGroup = getGroup(*son); 
      vector<unsigned int> subgroup(sonGroup.size());
      unsigned int index = subcount;
      for(unsigned int j = 0; j < sonGroup.size(); j++)
      {
        group.push_back(sonGroup[j]);
        subgroup[j] = subcount;
        thisgroup.push_back(subcount);
        subcount++;
      }
      if(subgroup.size() > 1) {
        //group.addSubgroup(subgroup, sonGroup.getHeight());
        // Recursively add subgroups:
        for(unsigned int j = 0; j < sonGroup.getSubgroups().size(); j++)
        {
          group.addSubgroup(sonGroup.getSubgroup(j) + index, sonGroup.getSubgroupHeight(j));
        }
      }
      group.setHeight(sonGroup.getHeight() + son->getDistanceToFather());
    }
    group.addSubgroup(thisgroup, group.getHeight());
  }
  return group;
}

unsigned int ClusterTools::getSubtreesWithSize(const Node * subtree, unsigned int size, vector<const Node *> & subtrees)
{
  if(subtree->isLeaf()) return 1;
  unsigned int sizeOfThisNode = 0;
  vector<unsigned int> sonSizes(subtree->getNumberOfSons());
  for(unsigned int i = 0; i < subtree->getNumberOfSons(); i++)
  {
    sonSizes[i] = getSubtreesWithSize(subtree->getSon(i), size, subtrees);
    sizeOfThisNode += sonSizes[i];
  }
  if(sizeOfThisNode == size) subtrees.push_back(subtree);
  if(sizeOfThisNode > size) {
    for(unsigned int i = 0; i < sonSizes.size(); i++)
    {
      if(sonSizes[i] < size && sonSizes[i] > 1) subtrees.push_back(subtree->getSon(i));
    }
  }
  return sizeOfThisNode;
}

vector<const Node *> ClusterTools::getSubtreesWithSize(const TreeTemplate<Node> & tree, unsigned int size)
{
  vector<const Node *> subtrees;
  getSubtreesWithSize(tree.getRootNode(), size, subtrees);
  return subtrees;
}

vector<Group> ClusterTools::getGroupsWithSize(const TreeTemplate<Node> & tree, unsigned int size)
{
  vector<const Node *> nodes = getSubtreesWithSize(tree, size);
  vector<Group> groups(nodes.size());
  for(unsigned int i = 0; i < nodes.size(); i++)
  {
    groups[i] = getGroup(* nodes[i]);
  }
  return groups;
}

void ClusterTools::computeGlobalDistanceDistribution(
  const HomogeneousSequenceSimulator & simulator,
	const SubstitutionCount & nijt,
  const Distance & distance,
  AgglomerativeDistanceMethod & clustering,
  const vector<double> & scales,
  unsigned int sizeOfDataSet,
  unsigned int nrep,
  ContingencyTable & table,
  ofstream * out)
{
  DRHomogeneousTreeLikelihood drhtl(
	  const_cast<TreeTemplate<Node> *>(simulator.getTree()),
		const_cast<SubstitutionModel *>(simulator.getSubstitutionModel()),
		const_cast<DiscreteDistribution *>(simulator.getRateDistribution()),
		true, false);
  vector<string> siteNames(sizeOfDataSet);
  for(unsigned int i = 0; i < sizeOfDataSet; i++)
  {
    siteNames[i] = TextTools::toString(i);
  }
  
  if(out != NULL) *out << "Size\tNmin\tDmax\tDelta" << endl;

  for(unsigned int k = 0; k < nrep; k++)
  {
    SiteContainer * sites = simulator.simulate(sizeOfDataSet);
    drhtl.setData(*sites);
    drhtl.computeTreeLikelihood();
    ProbabilisticSubstitutionMapping * mapping = SubstitutionMappingTools::computeSubstitutionVectors(drhtl, nijt, true);
    unsigned int nbBranches = mapping->getNumberOfBranches();
    
    //Mean vector:
    vector<double> meanVector(nbBranches);
    for(unsigned int j = 0; j < nbBranches; j++)
    {
		  double sum = 0;
      for(unsigned int i = 0; i < sizeOfDataSet; i++)
      {
        sum += (*mapping)(j, i);
      }
      meanVector[j] = sum / sizeOfDataSet;
    }

    //Scale vectors:
		for(unsigned int j = 0; j < nbBranches; j++)
    {
      double scale = scales[j];
		  for(unsigned int i = 0; i < sizeOfDataSet; i++)
      {
			  (*mapping)(j, i) *= scale;
			}
		}
    
    //Compute norms:
    vector<double> norms(sizeOfDataSet);
    for(unsigned int i = 0; i < sizeOfDataSet; i++)
    {
      norms[i] = norm((*mapping)[i]);
    }
    
    //Compute distance matrix:
	  DistanceMatrix * mat = new DistanceMatrix(siteNames);
	  for(unsigned int i = 0; i < sizeOfDataSet; i++)
    {
		  (*mat)(i,i) = 0.;
      Vdouble * vec = &((*mapping)[i]);
		  for(unsigned int j = 0; j < i; j++)
      {
			  (*mat)(i,j) = (*mat)(j,i) = distance.d(*vec, (*mapping)[j]);
		  }
	  }

    //Perform clustering:
    clustering.setDistanceMatrix(*mat);
    clustering.computeTree(true);
    TreeTemplate<Node> * tree = dynamic_cast<TreeTemplate<Node> *>(clustering.getTree());
    
    //Now parse each group:
    vector<Group> groups = getGroups(tree);
    for(unsigned int i = 0; i < groups.size(); i++)
    {
      Group * group = &groups[i];
      
      // Compute minimal norm:
      double normMin = -log(0.);
      for(unsigned int j = 0; j < group->size(); j++)
      {
        double norm = norms[TextTools::to<unsigned int>(group->at(j))];
        if(norm < normMin) normMin = norm;
      }
      
      // Compute distance from mean vector:
      vector<double> groupMeanVector(mapping->getNumberOfBranches(), 0.);
      for(unsigned int j = 0; j < group->size(); j++)
      {
        groupMeanVector += (*mapping)[TextTools::to<unsigned int>((*group)[j])]; 
      }
      double distFromMeanVector = distance.d(groupMeanVector/group->size(), meanVector);
      
      //try
      //{
      //  table.addToCount(group->size(), normMin, max(0., group->getHeight() * 2.));
      //}
      //catch(Exception & ex)
      //{
      //  cout << "WARNING: This group is out of bound: " << group->size() << "," << normMin << "," << (group->getHeight() * 2.) << "." << endl; 
      //}
      if(out != NULL) *out << group->size() << "\t" << normMin << "\t" << (group->getHeight() * 2.) <<"\t" << distFromMeanVector << endl;
    }

    //Housekeeping:
    delete sites;
    delete mapping;
    delete mat;
    delete tree;
  }
}

