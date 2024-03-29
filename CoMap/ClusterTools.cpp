//
// File: ClusterTools.cpp
// Created by: Julien Dutheil
// Created on: Thu Mar 2 2006
//

/*
Copyright or � or Copr. CNRS, (November 16, 2004, 2005, 2006)

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
#include "Cluster.h"

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// From bpp-phyl:
#include <Bpp/Phyl/Legacy/Likelihood/DRHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Mapping/SubstitutionMappingTools.h>

#include <iostream>
#include <fstream>
using namespace std;


vector<Group> ClusterTools::getGroups(const TreeTemplate<Node>* tree)
{
  vector<Group> groups;
  getGroups(*tree->getRootNode(), groups);
  return groups;
}
  
Group ClusterTools::getGroups(const Node& subtree, vector<Group>& groups)
{
  Group group;
  if (subtree.isLeaf())
  {
    group.add(subtree.getName());
  }
  else
  {
    size_t subcount = 0;
    vector<size_t> thisgroup;
    for (size_t i = 0; i < subtree.getNumberOfSons(); i++)
    {
      const Node * son = subtree.getSon(i);
      // Get group for each son node:
      Group sonGroup = getGroups(*son, groups); 
      vector<size_t> subgroup(sonGroup.size());
      size_t index = subcount;
      for (size_t j = 0; j < sonGroup.size(); j++)
      {
        group.add(sonGroup[j]);
        subgroup[j] = subcount;
        thisgroup.push_back(subcount);
        subcount++;
      }
      if(subgroup.size() > 1)
      {
        //group.addSubgroup(subgroup, sonGroup.getHeight());
        // Recursively add subgroups:
        for(size_t j = 0; j < sonGroup.getSubgroups().size(); j++)
        {
          group.addSubgroup(sonGroup.getSubgroup(j) + index, sonGroup.getSubgroupHeight(j));
        }
      }
      group.setHeight(sonGroup.getHeight() + son->getDistanceToFather());
    }
    vector<string> propNames = subtree.getNodePropertyNames();
    for(size_t i = 0; i < propNames.size(); i++)
    {
      group.setProperty(propNames[i], subtree.getNodeProperty(propNames[i]));
    }
    group.addSubgroup(thisgroup, group.getHeight());
    
    // Add this group:
    groups.push_back(group);
  }
  return group;
}

//Useful?
Group ClusterTools::getGroup(const Node& subtree)
{
  Group group;
  if (subtree.isLeaf())
  {
    group.add(subtree.getName());
  }
  else
  {
    size_t subcount = 0;
    vector<size_t> thisgroup;
    for (size_t i = 0; i < subtree.getNumberOfSons(); i++)
    {
      const Node * son = subtree.getSon(i);
      // Get group for each son node:
      Group sonGroup = getGroup(*son); 
      vector<size_t> subgroup(sonGroup.size());
      size_t index = subcount;
      for (size_t j = 0; j < sonGroup.size(); j++)
      {
        group.add(sonGroup[j]);
        subgroup[j] = subcount;
        thisgroup.push_back(subcount);
        subcount++;
      }
      if (subgroup.size() > 1)
      {
        //group.addSubgroup(subgroup, sonGroup.getHeight());
        // Recursively add subgroups:
        for(size_t j = 0; j < sonGroup.getSubgroups().size(); j++)
        {
          group.addSubgroup(sonGroup.getSubgroup(j) + index, sonGroup.getSubgroupHeight(j));
        }
      }
      group.setHeight(sonGroup.getHeight() + son->getDistanceToFather());
    }
    group.addSubgroup(thisgroup, group.getHeight());
  }
  vector<string> propNames = subtree.getNodePropertyNames();
  for(size_t i = 0; i < propNames.size(); i++)
  {
    group.setProperty(propNames[i], subtree.getNodeProperty(propNames[i]));
  }
  return group;
}

size_t ClusterTools::getSubtreesWithSize(const Node* subtree, size_t size, vector<const Node*>& subtrees)
{
  if(subtree->isLeaf()) return 1;
  size_t sizeOfThisNode = 0;
  vector<size_t> sonSizes(subtree->getNumberOfSons());
  for(size_t i = 0; i < subtree->getNumberOfSons(); i++)
  {
    sonSizes[i] = getSubtreesWithSize(subtree->getSon(i), size, subtrees);
    sizeOfThisNode += sonSizes[i];
  }
  if(sizeOfThisNode == size) subtrees.push_back(subtree);
  if(sizeOfThisNode > size) {
    for(size_t i = 0; i < sonSizes.size(); i++)
    {
      if(sonSizes[i] < size && sonSizes[i] > 1) subtrees.push_back(subtree->getSon(i));
    }
  }
  return sizeOfThisNode;
}

vector<const Node *> ClusterTools::getSubtreesWithSize(const TreeTemplate<Node>& tree, size_t size)
{
  vector<const Node *> subtrees;
  getSubtreesWithSize(tree.getRootNode(), size, subtrees);
  return subtrees;
}

vector<Group> ClusterTools::getGroupsWithSize(const TreeTemplate<Node>& tree, size_t size)
{
  vector<const Node *> nodes = getSubtreesWithSize(tree, size);
  vector<Group> groups(nodes.size());
  for(size_t i = 0; i < nodes.size(); i++)
  {
    groups[i] = getGroup(* nodes[i]);
  }
  return groups;
}

void ClusterTools::computeGlobalDistanceDistribution(
  shared_ptr<DRTreeLikelihoodInterface> drtl,
  const SequenceSimulatorInterface& simulator,
  shared_ptr<SubstitutionCountInterface> nijt,
  const Distance& distance,
  AgglomerativeDistanceMethodInterface& clustering,
  size_t sizeOfDataSet,
  size_t nrep,
  size_t maxGroupSize,
  ofstream* out)
{
  size_t nbTypes = nijt->getNumberOfSubstitutionTypes();
  vector<string> siteNames(sizeOfDataSet);
  for (size_t i = 0; i < sizeOfDataSet; i++)
  {
    siteNames[i] = TextTools::toString(i);
  }
  
  if (out)
    *out << "Rep\tGroup\tSize\tDmax\tStat\tNmin" << endl;

  for (size_t k = 0; k < nrep; k++)
  {
    ApplicationTools::displayGauge(k, nrep-1, '>');
    auto sites = simulator.simulate(sizeOfDataSet);
    drtl->setData(*sites);
    drtl->initialize();
    auto mapping = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, nijt, false);
    size_t nbBranches = mapping->getNumberOfBranches();
    
    //Mean vector:
    vector<double> meanVector(nbBranches);
    for (size_t j = 0; j < nbBranches; j++)
    {
      double sum = 0;
      for (size_t i = 0; i < sizeOfDataSet; i++)
        for (size_t t = 0; t < nbTypes; t++)
          sum += (*mapping)(j, i, t);
      meanVector[j] = sum / static_cast<double>(sizeOfDataSet);
    }

    //Compute distance matrix:
	  unique_ptr<DistanceMatrix> mat(new DistanceMatrix(siteNames));
	  for (size_t i = 0; i < sizeOfDataSet; ++i)
    {
		  (*mat)(i, i) = 0.;
      VVdouble* vec = &((*mapping)[i]);
		  for (size_t j = 0; j < i; ++j)
      {
			  (*mat)(i,j) = (*mat)(j,i) = distance.getDistanceForPair(*vec, (*mapping)[j]);
		  }
	  }

    //Perform clustering:
    try
    {
      dynamic_cast<SumClustering&>(clustering).setMapping(*mapping);
    }
    catch(exception& e) {}
    
    clustering.setDistanceMatrix(*mat);
    clustering.computeTree();
    TreeTemplate<Node> clusteringTree(clustering.tree());

    //Add information to tree:
    computeNormProperties(clusteringTree, *mapping);
    distance.setStatisticAsProperty(*clusteringTree.getRootNode(), *mapping);
    
    //Now parse each group:
    vector<Group> groups = getGroups(&clusteringTree);
    for (size_t i = 0; i < groups.size(); i++)
    {
      Group* group = &groups[i];
      if (group->size() > maxGroupSize) continue;
      
      //// Compute distance from mean vector:
      //vector<double> groupMeanVector(mapping->getNumberOfBranches(), 0.);
      //for(size_t j = 0; j < group->size(); j++)
      //{
      //  groupMeanVector += (*mapping)[TextTools::to<size_t>((*group)[j])]; 
      //}
      //double distFromMeanVector = distance.getDistanceForPair(groupMeanVector/group->size(), meanVector);
      
      if (out) *out << k
        << "\t" << group->toString()
        << "\t" << group->size()
        << "\t" << (group->getHeight() * 2.)
        << "\t" << (dynamic_cast<const Number<double> *>(group->getProperty("Stat"))->getValue())
        << "\t" << (dynamic_cast<const Number<double> *>(group->getProperty("Nmin"))->getValue())
        //<< "\t" << distFromMeanVector
        << endl;
    }
  }
  ApplicationTools::message->endLine();
}

void ClusterTools::computeNormProperties(TreeTemplate<Node>& tree, const LegacyProbabilisticSubstitutionMapping& mapping)
{
  double min;
  computeNormProperties_(tree.getRootNode(), mapping, min);
}
    
void ClusterTools::computeNormProperties_(Node* node, const LegacyProbabilisticSubstitutionMapping& mapping, double & minNorm)
{
  minNorm = -log(0.);
  if (node->isLeaf())
  {
    minNorm = LegacySubstitutionMappingTools::computeNormForSite(mapping, TextTools::to<size_t>(node->getName()));
  }
  else
  {
    for (size_t i = 0; i < node->getNumberOfSons(); ++i)
    {
      double minNormSon;
      computeNormProperties_(node->getSon(i), mapping, minNormSon);
      if (minNormSon < minNorm) minNorm = minNormSon;
    }
    node->setNodeProperty("Nmin", Number<double>(minNorm));
  }
}    

