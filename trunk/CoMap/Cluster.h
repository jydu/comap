//
// File: Cluster.h
// Created by: Julien Dutheil
// Created on: Tue Aug 30 17:19 2005
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

#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include "Distance.h"

// From PhyLib:
#include <Phyl/AbstractAgglomerativeDistanceMethod.h>
#include <Phyl/ProbabilisticSubstitutionMapping.h>

struct ClusterInfos
{
	unsigned int numberOfLeaves;
	double length;
};

class SimpleClustering:
  public AbstractAgglomerativeDistanceMethod
{
	public:
		static const string COMPLETE; 
		static const string SINGLE; 
		static const string AVERAGE; 
		static const string MEDIAN; 
		static const string WARD; 
		static const string CENTROID; 
	
	protected:
		string method_;
	
	public:
		SimpleClustering(const string & method, bool verbose = false): AbstractAgglomerativeDistanceMethod(verbose), method_(method) {}
		SimpleClustering(const string & method, const DistanceMatrix & matrix, bool verbose = false) throw (Exception): AbstractAgglomerativeDistanceMethod(matrix, verbose), method_(method) 
		{
			computeTree(true);
		}
		~SimpleClustering() {}

	public:
		TreeTemplate<Node> * getTree() const;

	protected:
		vector<unsigned int> getBestPair() throw (Exception);
		vector<double> computeBranchLengthsForPair(const vector<unsigned int> & pair);
		double computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos);
		void finalStep(int idRoot);	
		virtual Node * getLeafNode(int id, const string & name);
		virtual Node * getParentNode(int id, Node * son1, Node * son2);
};

/**
 * @brief Specific clustering.
 *
 * Compute the distance between two groups as the pairwise distance from the
 * two sums of vectors for the two groups.
 */
class SumClustering : public AbstractAgglomerativeDistanceMethod
{
	protected:
		ProbabilisticSubstitutionMapping mapping_;
		const Distance* distance_;
	
	public:
		SumClustering(const ProbabilisticSubstitutionMapping& mapping, const Distance* distance, const DistanceMatrix& matrix) throw (Exception) :
      AbstractAgglomerativeDistanceMethod(matrix), mapping_(mapping), distance_(distance)
		{
			computeTree(true);
		}
    SumClustering(const SumClustering& sc) :
      AbstractAgglomerativeDistanceMethod(sc), mapping_(sc.mapping_), distance_(sc.distance_) 
    {}
    SumClustering& operator=(const SumClustering& sc)
    {
      AbstractAgglomerativeDistanceMethod::operator=(sc);
      mapping_ = sc.mapping_;
      distance_ = sc.distance_;
      return *this;
    }
		virtual ~SumClustering() {}

	public:
		TreeTemplate<Node>* getTree() const;
		void setMapping(const ProbabilisticSubstitutionMapping & mapping)
		{
			mapping_ = mapping;
		}

	protected:
		/**
		 * @brief Returns the pair with minimum distance and actualizes the vectors.
		 *
		 * The vector at position bestPair[0] is now the sum of vectors bestPair[0] and bestPair[1].
		 * It is then used for computation of distances.
		 */
		vector<unsigned int> getBestPair() throw (Exception);
		vector<double> computeBranchLengthsForPair(const vector<unsigned int> & pair);
		double computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos);
		void finalStep(int idRoot);	
		virtual Node * getLeafNode(int id, const string & name);
		virtual Node * getParentNode(int id, Node * son1, Node * son2);
};

#endif //_CLUSTER_H_

