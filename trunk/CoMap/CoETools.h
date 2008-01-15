//
// File: CoETools.h
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

#ifndef _COETOOLS_H_
#define _COETOOLS_H_

#include "Statistics.h"

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/VectorTools.h>

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/SubstitutionModel.h>
#include <Phyl/DRTreeLikelihood.h>
#include <Phyl/MutationProcess.h>
#include <Phyl/DRTreeLikelihood.h>
#include <Phyl/DiscreteRatesAcrossSitesTreeLikelihood.h>
#include <Phyl/SubstitutionCount.h>
#include <Phyl/ProbabilisticSubstitutionMapping.h>
#include <Phyl/SequenceSimulator.h>

using namespace bpp;

// From the STL:
#include <map>
#include <vector>
#include <deque>

using namespace std;

class CandidateSite
{
  protected:
    unsigned int _index;
    double _normMin, _normMax;

  public:
    CandidateSite(unsigned int index): _index(index), _normMin(0), _normMax(0) {}

  public:
    void setNormRange(double min, double max)
    {
      _normMin = min;
      _normMax = max;
    }
    bool checkNorm(double norm) const
    {
      return norm >= _normMin && norm <= _normMax;
    }
    unsigned int getIndex() const { return _index; }
};

class CandidateGroup:
  public vector<CandidateSite>
{
  protected:
    double _statistic;

  public:
    CandidateGroup(): _statistic(0) {}

  public:
    double getStatisticValue() const { return _statistic; }
    void setStatisticValue(double value) { _statistic = value; }
    void computeStatisticValue(const Statistic & stat, const ProbabilisticSubstitutionMapping & mapping) throw (Exception)
    {
      if(size() == 0) throw Exception("CandidateGroup::computeStatisticValue. Group is empty!");
      vector<const Vdouble *> group;
      for(unsigned int i = 0; i < size(); i++)
      {
        unsigned int index = (*this)[i].getIndex();
        group.push_back(&mapping[index]);
      }
      _statistic = stat.getValueForGroup(group);
    }
    void computeNormRanges(double omega, const ProbabilisticSubstitutionMapping & mapping)
    {
      for(unsigned int i = 0; i < size(); i++)
      {
        double norm = VectorTools::norm<double, double>(mapping[(*this)[i].getIndex()]);
        (*this)[i].setNormRange(norm - omega, norm + omega);
      }
    }

};

class CandidateGroupSet
{
  protected:
    const Statistic * _statistic;
    vector<CandidateGroup> _candidates;

    //Simulated site(s) for each group, each site...
    //_simulations[i][j][k] is the kth simulated site for site j in group i.
    vector< vector< deque<const Vdouble *> > > _simulations;
    //Contain the number of simulations for each group.
    vector<unsigned int> _n1, _n2;
    //The minimum number of simulation for each group:
    unsigned int _minSim;

    unsigned int _nbCompleted;

    unsigned int _verbose;

    //Iterator:
    mutable unsigned int _groupPos, _sitePos;

  public:
    CandidateGroupSet(const Statistic * statistic, unsigned int minSim, unsigned int verbose = 1): _statistic(statistic), _minSim(minSim), _nbCompleted(0), _verbose(verbose), _groupPos(0), _sitePos(0) {}

  public:
    /**
     * @brief Analyse a set of simulated substitution maps.
     *
     * All simulated sites will be checked against the set.
     * @param simMap Simulated sites.
     * @eturn true if the minimum number of simulations is reached for every group.
     */
    bool analyseSimulations(const ProbabilisticSubstitutionMapping & simMap);

    void addCandidate(const CandidateGroup & group)
    {
      _candidates.push_back(group);
      _simulations.push_back(vector< deque<const Vdouble *> >(group.size()));
      _n1.push_back(0);
      _n2.push_back(0);
    }

    double getPValueForGroup(unsigned int groupIndex) const
    {
      return ((double)_n1[groupIndex] + 1.) / ((double)_n2[groupIndex] + 1.);
    }

    double getN1ForGroup(unsigned int groupIndex) const { return _n1[groupIndex]; }
    double getN2ForGroup(unsigned int groupIndex) const { return _n2[groupIndex]; }

    unsigned int size() const { return _candidates.size(); }

    const CandidateGroup & operator[](unsigned int index) const { return _candidates[index]; }

    void resize(unsigned int size)
    {
      _candidates.resize(size);
      _simulations.resize(size);
      _n1.resize(size);
      _n2.resize(size);
    }

    unsigned int getVerbose() const { return _verbose; }

 protected:

    /**
     * @brief Candidate site iterator.
     *
     * This iterator loop over all sites and over all groups for which _n2 is lovaer than _minSim.
     *
     * @return A 2-vector, the first dimenson giving the group index, and the second one the index of the site in the group.
     * @throw Exception If there is no more candidate site!
     */
    vector<unsigned int> nextCandidateSite() const throw (Exception);
    
    vector<unsigned int> currentCandidateSite() const throw (Exception);

    /**
     * @brief Add a new simulated vector to the set.
     *
     * This will check each site in the group, and actualize _n1 and _n2.
     *
     * @param groupIndex Position of the group.
     * @param siteIndex  Position of the site in the group.
     * @param v          A pointer toward the substitution vector to add.
     */
    void addSimulatedSite(unsigned int groupIndex, unsigned int siteIndex, const Vdouble * v) throw (IndexOutOfBoundsException);

    /**
     * @brief Remove all pointers in _simulations, but do not modify _n1 and _n2.
     */
    void resetSimulations()
    {
      for(unsigned int i = 0; i < _simulations.size(); i++)
        for(unsigned int j = 0; j < _simulations[i].size(); j++)
          _simulations[i][j].clear();
    }

};



class CoETools
{
	public:
		CoETools() {};
		virtual ~CoETools() {};
	
	public:
		
		static void readData(
			TreeTemplate<Node> *   tree,
			Alphabet *             &alphabet,
			VectorSiteContainer *  &allSites,
			VectorSiteContainer *  &sites,
			SubstitutionModel *    &model,
			SubstitutionModelSet * &modelSet,
			DiscreteDistribution * &rDist,
			DRTreeLikelihood *     &tl,
			map<string, string>    &params,
			const string           &suffix = "");

		static ProbabilisticSubstitutionMapping * getVectors(
			const DRTreeLikelihood & drtl, 
      SubstitutionCount      & substitutionCount,
      const SiteContainer    & completeSites,
			map<string, string>    & params,
			const string           & suffix = "");
			
	  /**
 	   * This write a file with some information for each site in the selected alignment
     * (ie selected sites without gaps).
     * Information includes:
     * - rate classe
     * - constant y/n
     * This may be helpful in order to know what are the sites realy used in the analysis.
     * The output is in Co-E' output format, ie a site set followed by properties:
     * Here the site set contains only one site and the properties are
     * - rate class [integer]
     * - constant [boolean (0 = false)]
     */
		static void writeInfos(
			const SiteContainer & completeSites,
			const DiscreteRatesAcrossSitesTreeLikelihood & ras,
			map<string, string> & params,
			const string & suffix = "");

		static int getMinRateClass(map<string, string> & params, string suffix = "");
		
		static double  getMinRate(map<string, string> & params, string suffix = "");
		
		static int getMaxRateClassDiff(map<string, string> & params);
		
		static double getMaxRateDiff(map<string, string> & params);

		static double getStatisticMin(map<string, string> & params);
		
		static bool haveToPerformIndependantComparisons(map<string, string> & params);
		
		static const Statistic * getStatistic(map<string, string> & params) throw (Exception);
		
		static SubstitutionCount * getSubstitutionCount(
			const Alphabet* alphabet,
			const SubstitutionModel* model,
			const DiscreteDistribution* rDist,
			map<string, string> & params,
			string suffix = "");

    static void computeIntraStats(
			const DiscreteRatesAcrossSitesTreeLikelihood & tl,
			const SiteContainer & completeSites,
			ProbabilisticSubstitutionMapping & mapping,
			const Statistic & statistic,
			map<string, string> & params);

		static void computeInterStats(
			const DiscreteRatesAcrossSitesTreeLikelihood & tl1,
			const DiscreteRatesAcrossSitesTreeLikelihood & tl2,
			const SiteContainer & completeSites1,
			const SiteContainer & completeSites2,
			ProbabilisticSubstitutionMapping & mapping1,
			ProbabilisticSubstitutionMapping & mapping2,
			const Statistic & statistic,
			map<string, string> & params);
		
		static void computeIntraNullDistribution(
        DRTreeLikelihood & drtl,
  			const SequenceSimulator& seqSim,
	  		SubstitutionCount & nijt,
		  	const Statistic & statistic,
			  map<string, string> & params);
	
		static void computeInterNullDistribution(
        DRTreeLikelihood & drtl1,
        DRTreeLikelihood & drtl2,
			  const SequenceSimulator& seqSim1,
			  const SequenceSimulator& seqSim2,
			  SubstitutionCount & nijt1,
			  SubstitutionCount & nijt2,
			  const Statistic & statistic,
			  map<string, string> & params);

    static void computePValuesForCandidateGroups(
        CandidateGroupSet & candidates,
        DRTreeLikelihood & drtl,
			  const SequenceSimulator& seqSim,
        SubstitutionCount & nijt,
			  map<string, string> & params);
	
};

#endif	//_COETOOLS_H_

