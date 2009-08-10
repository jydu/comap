//
// File: Cluster.cpp
// Created by: Julien Dutheil
// Created on: Tue Aug 31 09:29 2005
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

#include "Cluster.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>

// From PhylLib:
#include <Phyl/NodeTemplate.h>
#include <Phyl/Tree.h>
#include <Phyl/TreeTemplate.h>

// From the STL:
#include <cmath>
#include <iostream>

using namespace std;
	
const string SimpleClustering::COMPLETE = "Complete"; 
const string SimpleClustering::SINGLE   = "Single"; 
const string SimpleClustering::AVERAGE  = "Average"; 
const string SimpleClustering::MEDIAN   = "Median"; 
const string SimpleClustering::WARD     = "Ward"; 
const string SimpleClustering::CENTROID = "Centroid"; 

TreeTemplate<Node> * SimpleClustering::getTree() const
{
	Node* root = TreeTemplateTools::cloneSubtree<Node>(*dynamic_cast<TreeTemplate<NodeTemplate<ClusterInfos> > *>(_tree) -> getRootNode());
	return new TreeTemplate<Node>(root);
}

vector<unsigned int> SimpleClustering::getBestPair() throw (Exception)
{
	vector<unsigned int> bestPair(2);
	double distMin = -std::log(0.);
	for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++)
  {
		unsigned int id = i->first;
		map<unsigned int, Node *>::iterator j = i;
		j++;
		for(; j != _currentNodes.end(); j++)
    {
			unsigned int jd = j -> first;
			double dist = _matrix(id, jd);
			if(dist < distMin)
      {
				distMin = dist;
				bestPair[0] = id;
				bestPair[1] = jd;
			}
		}
	}
	if(distMin == -std::log(0.))
  {
    cout << "---------------------------------------------------------------------------------" << endl;
	  for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++)
    {
		  unsigned int id = i -> first;
		  map<unsigned int, Node *>::iterator j = i;
		  j++;
		  for(; j != _currentNodes.end(); j++) {
			  unsigned int jd = j -> first;
			  double dist = _matrix(id, jd);
        cout << dist << "\t";
		  }
      cout << endl;
	  }
    cout << "---------------------------------------------------------------------------------" << endl;
    
    throw Exception("Unexpected error: no minimum found in the distance matrix.");
	}

	return bestPair;	
}
vector<double> SimpleClustering::computeBranchLengthsForPair(const vector<unsigned int> & pair)
{
	vector<double> d(2);
	double dist = _matrix(pair[0], pair[1]) / 2.;
	d[0] = dist - dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[0]])->getInfos().length; 
	d[1] = dist - dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[1]])->getInfos().length; 
	return d;
}

double SimpleClustering::computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos)
{
	double w1, w2, w3, w4;
	if(_method == "Single")
  {
		w1 = .5;
		w2 = .5;
		w3 = 0.;
		w4 = -.5;
	} else if(_method == "Complete") {
		w1 = .5;
		w2 = .5;
		w3 = 0.;
		w4 = .5;
	} else if(_method == "Median") {
		w1 = .5;
		w2 = .5;
		w3 = -0.25;
		w4 = 0.;
	} else if(_method == "Average") {
		double n1 = dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[0]]) -> getInfos().numberOfLeaves;
		double n2 = dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[1]]) -> getInfos().numberOfLeaves;
		w1 = n1/(n1+n2);
		w2 = n2/(n1+n2);
		w3 = 0.;
		w4 = 0.;
	} else if(_method == "Ward") {
		double n1 = dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[0]]) -> getInfos().numberOfLeaves;
		double n2 = dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[1]]) -> getInfos().numberOfLeaves;
		double n3 = dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pos])     -> getInfos().numberOfLeaves;
		w1 = (n1+n3)/(n1+n2+n3);
		w2 = (n2+n3)/(n1+n2+n3);
		w3 = -n3/(n1+n2+n3);
		w4 = 0.;
	} else if(_method == "Centroid") {
		double n1 = dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[0]]) -> getInfos().numberOfLeaves;
		double n2 = dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[1]]) -> getInfos().numberOfLeaves;
		w1 = n1/(n1+n2);
		w2 = n2/(n1+n2);
		w3 = -n1*n2/pow(n1+n2, 2.);
		w4 = 0.;
	} else throw Exception("SimpleClustering::computeBranchLengthsForPair. unknown method '" + _method + "'.");
	double d1 = _matrix(pair[0], pos);
	double d2 = _matrix(pair[1], pos);
	double d3 = _matrix(pair[0], pair[1]);
	return w1*d1 + w2*d2 + w3*d3 + w4*std::abs(d1-d2);
}

void SimpleClustering::finalStep(int idRoot)
{
	NodeTemplate<ClusterInfos>* root = new NodeTemplate<ClusterInfos>(idRoot);
	map<unsigned int, Node*>::iterator it = _currentNodes.begin();
	unsigned int i1 = it->first;
	Node* n1        = it->second;
	it++;
	unsigned int i2 = it->first;
	Node* n2        = it->second;
	double d = _matrix(i1, i2) / 2;
	root->addSon(n1);
	root->addSon(n2);
	n1->setDistanceToFather(d - dynamic_cast<NodeTemplate<ClusterInfos>*>(n1)->getInfos().length); 
	n2->setDistanceToFather(d - dynamic_cast<NodeTemplate<ClusterInfos>*>(n2)->getInfos().length); 
	_tree = new TreeTemplate<NodeTemplate<ClusterInfos> >(root);
}

Node * SimpleClustering::getLeafNode(int id, const string& name)
{
	ClusterInfos infos;
	infos.numberOfLeaves = 1;
	infos.length = 0.;
	NodeTemplate<ClusterInfos> * leaf = new NodeTemplate<ClusterInfos>(id, name);
	leaf->setInfos(infos);
	return leaf;
}

Node * SimpleClustering::getParentNode(int id, Node* son1, Node* son2)
{
	ClusterInfos infos;
	infos.numberOfLeaves = 
		dynamic_cast<NodeTemplate<ClusterInfos> *>(son1)->getInfos().numberOfLeaves
	+ dynamic_cast<NodeTemplate<ClusterInfos> *>(son2)->getInfos().numberOfLeaves;
	infos.length = dynamic_cast<NodeTemplate<ClusterInfos> *>(son1)->getInfos().length + son1 -> getDistanceToFather();
	Node* parent = new NodeTemplate<ClusterInfos>(id);
	dynamic_cast<NodeTemplate<ClusterInfos> *>(parent)->setInfos(infos);
	parent->addSon(son1);
	parent->addSon(son2);
	return parent;
}

//---------------------------------------------------------------------------------------------

TreeTemplate<Node>* SumClustering::getTree() const
{
	Node * root = TreeTemplateTools::cloneSubtree<Node>(* dynamic_cast<TreeTemplate<NodeTemplate<ClusterInfos> > *>(_tree) -> getRootNode());
	return new TreeTemplate<Node>(root);
}

vector<unsigned int> SumClustering::getBestPair() throw (Exception)
{
	vector<unsigned int> bestPair(2);
	double distMin = -std::log(0.);
	for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++)
  {
		unsigned int id = i -> first;
		map<unsigned int, Node *>::iterator j = i;
		j++;
		for(; j != _currentNodes.end(); j++)
    {
			unsigned int jd = j->first;
			double dist = _matrix(id, jd);
			if(dist < distMin)
      {
				distMin = dist;
				bestPair[0] = id;
				bestPair[1] = jd;
			}
		}
	}
	// actualize vectors:
	_mapping[bestPair[0]] += _mapping[bestPair[1]];
	return bestPair;	
}
vector<double> SumClustering::computeBranchLengthsForPair(const vector<unsigned int>& pair)
{
	vector<double> d(2);
	double dist = _matrix(pair[0], pair[1]) / 2.;
	d[0] = dist - dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[0]])->getInfos().length; 
	d[1] = dist - dynamic_cast<NodeTemplate<ClusterInfos> *>(_currentNodes[pair[1]])->getInfos().length; 
	//d[0] = dist; 
	//d[1] = dist; 
	return d;
}

double SumClustering::computeDistancesFromPair(const vector<unsigned int>& pair, const vector<double>& branchLengths, unsigned int pos)
{
	return _distance->getDistanceForPair(_mapping[pos], _mapping[pair[0]]);
}

void SumClustering::finalStep(int idRoot)
{
	NodeTemplate<ClusterInfos>* root = new NodeTemplate<ClusterInfos>(idRoot);
	map<unsigned int, Node*>::iterator it = _currentNodes.begin();
	unsigned int i1 = it->first;
	Node* n1        = it->second;
	it++;
	unsigned int i2 = it->first;
	Node* n2        = it->second;
	double d = _matrix(i1, i2) / 2;
	root->addSon(n1);
	root->addSon(n2);
	n1->setDistanceToFather(d - dynamic_cast<NodeTemplate<ClusterInfos> *>(n1)->getInfos().length); 
	n2->setDistanceToFather(d - dynamic_cast<NodeTemplate<ClusterInfos> *>(n2)->getInfos().length); 
	//n1->setDistanceToFather(d); 
	//n2->setDistanceToFather(d); 
	_tree = new TreeTemplate< NodeTemplate<ClusterInfos> >(root);
}

Node * SumClustering::getLeafNode(int id, const string & name)
{
	ClusterInfos infos;
	infos.numberOfLeaves = 1;
	infos.length = 0.;
	NodeTemplate<ClusterInfos> * leaf = new NodeTemplate<ClusterInfos>(id, name);
	leaf->setInfos(infos);
	return leaf;
}

Node * SumClustering::getParentNode(int id, Node * son1, Node * son2)
{
	ClusterInfos infos;
	infos.numberOfLeaves = 
		dynamic_cast<NodeTemplate<ClusterInfos> *>(son1)->getInfos().numberOfLeaves
	+ dynamic_cast<NodeTemplate<ClusterInfos> *>(son2)->getInfos().numberOfLeaves;
	infos.length = dynamic_cast<NodeTemplate<ClusterInfos> *>(son1)->getInfos().length + son1->getDistanceToFather();
	Node * parent = new NodeTemplate<ClusterInfos>(id);
	dynamic_cast<NodeTemplate<ClusterInfos> *>(parent)->setInfos(infos);
	parent->addSon(son1);
	parent->addSon(son2);
	return parent;
}

//---------------------------------------------------------------------------------------------


