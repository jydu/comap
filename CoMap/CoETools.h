//
// File: CoETools.h
// Created by: Julien Dutheil
// Created on: Mon Feb  2 14:50:40 2004
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

#ifndef _COETOOLS_H_
#define _COETOOLS_H_

#include "Statistics.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Legacy/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Legacy/Likelihood/DRTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Mapping/ProbabilisticSubstitutionMapping.h>
#include <Bpp/Phyl/Legacy/Mapping/SubstitutionMappingTools.h>
#include <Bpp/Phyl/Simulation/SequenceSimulator.h>

using namespace bpp;

// From the STL:
#include <map>
#include <vector>
#include <deque>
#include <memory>

using namespace std;

class CandidateSite
{
  protected:
    size_t index_;
    double normMin_, normMax_;

  public:
    CandidateSite(size_t index): index_(index), normMin_(0), normMax_(0) {}

  public:
    void setNormRange(double min, double max)
    {
      normMin_ = min;
      normMax_ = max;
    }
    bool checkNorm(double norm) const
    {
      return (norm >= normMin_ && norm <= normMax_);
    }
    size_t getIndex() const { return index_; }
};

class CandidateGroup
{
  protected:
    vector<CandidateSite> sites_;
    double statistic_;
    bool analysable_;

  public:
    CandidateGroup(): sites_(), statistic_(0), analysable_(true) {}

  public:
    double getStatisticValue() const { return statistic_; }
    void setStatisticValue(double value) { statistic_ = value; }
    void computeStatisticValue(const Statistic& stat, const LegacyProbabilisticSubstitutionMapping& mapping)
    {
      if (!analysable_) throw Exception("CandidateGroup::computeStatisticValue. Group is not analysable.");
      if (sites_.size() == 0) throw Exception("CandidateGroup::computeStatisticValue. Group is empty!");
      vector<const VVdouble*> group(sites_.size());
      for (size_t i = 0; i < sites_.size(); i++)
      {
        size_t index = sites_[i].getIndex();
	group[i] = &mapping[index];
      }
      statistic_ = stat.getValueForGroup(group);
    }
    void computeNormRanges(double omega, const LegacyProbabilisticSubstitutionMapping& mapping)
    {
      if (!analysable_) throw Exception("CandidateGroup::computeNormRanges. Group is not analyzable.");
      if (sites_.size() == 0) throw Exception("CandidateGroup::computeNormRanges. Group is empty!");
      for (size_t i = 0; i < sites_.size(); i++)
      {
        //double norm = VectorTools::norm<double, double>(mapping[sites_[i].getIndex()]);
        double norm = LegacySubstitutionMappingTools::computeNormForSite(mapping, sites_[i].getIndex());
        sites_[i].setNormRange(norm - omega, norm + omega);
      }
    }
    void setAnalysable(bool yn) { analysable_ = yn; }
    bool isAnalysable() const { return analysable_; }
    size_t size() const { return sites_.size(); }

    const CandidateSite& operator[](unsigned int i) const { return sites_[i]; }
    CandidateSite& operator[](unsigned int i) { return sites_[i]; }

    void addSite(const CandidateSite& cs) { sites_.push_back(cs); }
};

class CandidateGroupSet
{
  protected:
    const Statistic* statistic_;
    vector<CandidateGroup> candidates_;

    //Simulated site(s) for each group, each site...
    //_simulations[i][j][k] is the kth simulated site for site j in group i.
    vector< vector< deque<const VVdouble*> > > simulations_;
    //Contain the number of simulations for each group.
    vector<unsigned int> n1_, n2_;
    //The minimum number of simulation for each group:
    unsigned int minSim_;

    unsigned int nbCompleted_;

    unsigned int nbAnalysable_;
    
    unsigned int nbTrials_;
    
    unsigned int verbose_;

    //Iterator:
    mutable unsigned int groupPos_, sitePos_;

  public:
    CandidateGroupSet(const Statistic* statistic, unsigned int minSim, unsigned int verbose = 1):
      statistic_(statistic),
      candidates_(),
      simulations_(),
      n1_(), n2_(),
      minSim_(minSim),
      nbCompleted_(0),
      nbAnalysable_(0),
      nbTrials_(0),
      verbose_(verbose),
      groupPos_(0),
      sitePos_(0)
    {}

    CandidateGroupSet(const CandidateGroupSet& cgs) :
      statistic_(cgs.statistic_),
      candidates_(cgs.candidates_),
      simulations_(cgs.simulations_),
      n1_(cgs.n1_), n2_(cgs.n2_),
      minSim_(cgs.minSim_),
      nbCompleted_(cgs.nbCompleted_),
      nbAnalysable_(cgs.nbAnalysable_),
      nbTrials_(cgs.nbTrials_),
      verbose_(cgs.verbose_),
      groupPos_(cgs.groupPos_),
      sitePos_(cgs.sitePos_)
    {}

    CandidateGroupSet& operator=(const CandidateGroupSet& cgs)
    {
      statistic_ = cgs.statistic_;
      candidates_ = cgs.candidates_;
      simulations_ = cgs.simulations_;
      n1_ = cgs.n1_;
      n2_ = cgs.n2_;
      minSim_ = cgs.minSim_;
      nbCompleted_ = cgs.nbCompleted_;
      nbAnalysable_ = cgs.nbAnalysable_;
      nbTrials_ = cgs.nbTrials_;
      verbose_ = cgs.verbose_;
      groupPos_ = cgs.groupPos_;
      sitePos_ = cgs.sitePos_;
      return *this;
    }

  public:
    /**
     * @brief Analyse a set of simulated substitution maps.
     *
     * All simulated sites will be checked against the set.
     * @param simMap Simulated sites.
     * @eturn true if the minimum number of simulations is reached for every group.
     */
    bool analyseSimulations(const LegacyProbabilisticSubstitutionMapping& simMap);

    void addCandidate(const CandidateGroup& group)
    {
      candidates_.push_back(group);
      simulations_.push_back(vector< deque<const VVdouble*> >(group.size()));
      n1_.push_back(0);
      n2_.push_back(0);
      if (group.isAnalysable()) nbAnalysable_++;
    }

    double getPValueForGroup(unsigned int groupIndex) const
    {
      return (static_cast<double>(n1_[groupIndex]) + 1.) / (static_cast<double>(n2_[groupIndex]) + 1.);
    }

    double getN1ForGroup(unsigned int groupIndex) const { return n1_[groupIndex]; }
    double getN2ForGroup(unsigned int groupIndex) const { return n2_[groupIndex]; }

    size_t size() const { return candidates_.size(); }

    const CandidateGroup& operator[](unsigned int index) const { return candidates_[index]; }

    void resize(unsigned int newSize)
    {
      candidates_.resize(newSize);
      simulations_.resize(newSize);
      n1_.resize(newSize);
      n2_.resize(newSize);
      nbAnalysable_ = 0;
      for (unsigned int i = 0; i < newSize; i++)
        if (candidates_[i].isAnalysable()) nbAnalysable_++;
    }

    unsigned int getVerbose() const { return verbose_; }
    
    unsigned int getNumberOfTrials() const { return nbTrials_; }

 protected:

    /**
     * @brief Candidate site iterator.
     *
     * This iterator loop over all sites and over all groups for which n2_ is lovaer than minSim_.
     *
     * @return A 2-vector, the first dimenson giving the group index, and the second one the index of the site in the group.
     * @throw Exception If there is no more candidate site!
     */
    vector<unsigned int> nextCandidateSite() const;
    
    vector<unsigned int> currentCandidateSite() const;

    /**
     * @brief Add a new simulated vector to the set.
     *
     * This will check each site in the group, and actualize n1_ and n2_.
     *
     * @param groupIndex Position of the group.
     * @param siteIndex  Position of the site in the group.
     * @param v          A pointer toward the substitution vector to add.
     * @return 'true' if a group was completed.
     */
    bool addSimulatedSite(unsigned int groupIndex, unsigned int siteIndex, const VVdouble* v);

    /**
     * @brief Remove all pointers in simulations_, but do not modify n1_ and n2_.
     */
    void resetSimulations()
    {
      for (size_t i = 0; i < simulations_.size(); ++i)
        for (size_t j = 0; j < simulations_[i].size(); ++j)
          simulations_[i][j].clear();
    }

};



class CoETools
{
  public:
    CoETools() {};
    virtual ~CoETools() {};
  
  public:
    
    static void readData(
      shared_ptr<TreeTemplate<Node>>            & tree,
      shared_ptr<Alphabet>                      & alphabet,
      shared_ptr<GeneticCode>                   & geneticCode,
      shared_ptr<VectorSiteContainer>           & allSites,
      shared_ptr<VectorSiteContainer>           & sites,
      shared_ptr<SubstitutionModelInterface>    & model,
      shared_ptr<SubstitutionModelSet>          & modelSet,
      shared_ptr<DiscreteDistributionInterface> & rDist,
      shared_ptr<DRTreeLikelihoodInterface>     & tl,
      map<string, string>                       & params,
      const string                              & suffix = "");

    static std::unique_ptr<LegacyProbabilisticSubstitutionMapping> getVectors(
      std::shared_ptr<const DRTreeLikelihoodInterface> drtl, 
      std::shared_ptr<SubstitutionCountInterface> substitutionCount,
      const SiteContainerInterface    & completeSites,
      map<string, string>             & params,
      const string                    & suffix = "");
      
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
      const SiteContainerInterface& completeSites,
      const DiscreteRatesAcrossSitesTreeLikelihoodInterface& ras,
      const std::vector<double>& norms,
      map<string, string>& params,
      const string& suffix = "");

    static int getMinRateClass(map<string, string>& params, string suffix = "");
    
    static double  getMinRate(map<string, string>& params, string suffix = "");
    
    static int getMaxRateClassDiff(map<string, string>& params);
    
    static double getMaxRateDiff(map<string, string>& params);

    static double getStatisticMin(map<string, string>& params);
    
    static bool haveToPerformIndependantComparisons(map<string, string>& params);
    
    static unique_ptr<Statistic> getStatistic(
	map<string, string>& params,
       	shared_ptr<const Alphabet> alphabet,
       	shared_ptr<const SubstitutionCountInterface> nijt);
    
    /**
     * This subroutine will call computeIntraNullDistribution on request.
     */
    static void computeIntraStats(
      const DRTreeLikelihoodInterface& tl,
      const SequenceSimulatorInterface& seqSim,
      const SiteContainerInterface& completeSites,
      LegacyProbabilisticSubstitutionMapping& mapping,
      shared_ptr<SubstitutionCountInterface> nijt,
      const Statistic& statistic,
      bool computeNull,
      map<string, string>& params);

    static void computeInterStats(
      const DiscreteRatesAcrossSitesTreeLikelihoodInterface& tl1,
      const DiscreteRatesAcrossSitesTreeLikelihoodInterface& tl2,
      const SiteContainerInterface& completeSites1,
      const SiteContainerInterface& completeSites2,
      LegacyProbabilisticSubstitutionMapping& mapping1,
      LegacyProbabilisticSubstitutionMapping& mapping2,
      const Statistic& statistic,
      map<string, string>& params);
    
    static vector< vector<double> >* computeIntraNullDistribution(
        shared_ptr<DRTreeLikelihoodInterface> drtl,
        const Domain* rateDomain,
        const SequenceSimulatorInterface& seqSim,
        shared_ptr<SubstitutionCountInterface> nijt,
        const Statistic& statistic,
        map<string, string>& params);
  
    static void computeInterNullDistribution(
        shared_ptr<DRTreeLikelihoodInterface> drtl1,
        shared_ptr<DRTreeLikelihoodInterface> drtl2,
        const SequenceSimulatorInterface& seqSim1,
        const SequenceSimulatorInterface& seqSim2,
        shared_ptr<SubstitutionCountInterface> nijt1,
        shared_ptr<SubstitutionCountInterface> nijt2,
        const Statistic& statistic,
        map<string, string>& params);

    static void computePValuesForCandidateGroups(
        CandidateGroupSet& candidates,
        shared_ptr<DRTreeLikelihoodInterface> drtl,
        const SequenceSimulatorInterface& seqSim,
        shared_ptr<SubstitutionCountInterface> nijt,
        map<string, string>& params,
        unsigned int maxTrials);
  
};

#endif  //_COETOOLS_H_

