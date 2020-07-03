//
// File: Mica.cpp
// Created by: Julien Dutheil
// Created on: Sat Feb 28 07:42 2009
//

/*
Copyright or © or Copr. CNRS, (November 16, 2009)

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

// From the STL:

#include <iostream>
#include <cstdlib>
#include <stdexcept>

using namespace std;
#include "Domain.h"

// From bpp-core:
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/NumCalcApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Simulation/MutationProcess.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Mapping/ProbabilisticSubstitutionMapping.h>
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>
#include <Bpp/Phyl/Mapping/UniformizationSubstitutionCount.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine(); 
  (*ApplicationTools::message << "mica param=optionfile").endLine();
  (*ApplicationTools::message << "Check the CoMap manual for a full syntax specification.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

/******************************************************************************/

void miTest(const Site& site1, const Site& site2, size_t maxNbPermutations, double& mi, double& pvalue, size_t& nbPermutations)
{
  mi = SiteTools::mutualInformation(site1, site2, true);
  Site copyOfSite1 = site1;
  Site copyOfSite2 = site2;
  double rep;
  double count = 0;
  if (SiteTools::isConstant(site1, true) || SiteTools::isConstant(site2, true))
  {
    pvalue = 1.;
    nbPermutations = 0;
  }
  else
  {
    size_t i;
    for (i = 0; count < 5 && i < maxNbPermutations; i++)
    {
      copyOfSite1.shuffle();
      copyOfSite2.shuffle();
      rep = SiteTools::mutualInformation(copyOfSite1, copyOfSite2, true);
      if (rep >= mi) count++;
    }
    pvalue = static_cast<double>(count + 1) / (static_cast<double>(i + 1));
    nbPermutations = i;
  }
}

/******************************************************************************/


vector<double> computeNorms(const ProbabilisticSubstitutionMapping& mapping)
{
  size_t nbVectors = mapping.getNumberOfSites();
  vector<double> vect(nbVectors);
  for (size_t i = 0; i < nbVectors; i++)
    vect[i] = SubstitutionMappingTools::computeNormForSite(mapping, i);
  return vect;
}

/******************************************************************************/

class CoevolutionStatistic
{
  public:
    CoevolutionStatistic() {}
    virtual ~CoevolutionStatistic() {}

    virtual shared_ptr<CoevolutionStatistic> clone(shared_ptr<const SiteContainer> newSites) const = 0;

  public:
    virtual size_t getNumberOfStatistics() const = 0;
    virtual vector<string> getStatisticsNames() const = 0;
    virtual string getMainStatisticsName() const = 0; //Several statistics can be computed, say which one is the main one (for computing pvalues, etc.).
    virtual double getValue(size_t site1, size_t site2, const string& statName) const = 0;
    virtual double getValue(size_t site1, size_t site2) const = 0; //call main statistic by default
    virtual double getRate(size_t site) const = 0;
    virtual const vector<double>& getRates() const = 0;
    virtual double getMinRate(size_t site1, size_t site2) const = 0;
};


class AbstractCoevolutionStatistic:
    public virtual CoevolutionStatistic
{
  protected:
    shared_ptr<const SiteContainer> sites_;
    size_t nbSites_;
    mutable map<string, size_t> statNames_;
    bool storeResults_;
    mutable map<size_t, vector<double> > statistics_;
    vector<double> rates_;

  public:
    AbstractCoevolutionStatistic(shared_ptr<const SiteContainer> sites, const vector<string>& statNames, bool storeResults = true):
         sites_(sites), nbSites_(sites->getNumberOfSites()), statNames_(), storeResults_(storeResults), statistics_(), rates_(sites->getNumberOfSites())
    {
      for (size_t i = 0; i < statNames.size(); ++i) {
	statNames_[statNames[i]] = i;
      }
    }

    virtual ~AbstractCoevolutionStatistic() {}

  public:
    size_t getNumberOfStatistics() const { return statNames_.size(); }

    vector<string> getStatisticsNames() const { return MapTools::getKeys(statNames_); }

    double getValue(size_t site1, size_t site2, const string& statName) const { 
      auto iter = statNames_.find(statName);
      if (iter != statNames_.end()) {
        size_t key = site1 + (site2 * nbSites_);
        auto iter2 = statistics_.find(key);
	if (iter2 != statistics_.end()) {
	  return iter2->second[iter->second];
	} else {
	  if (!storeResults_) {
            //We remove any previous calculation for other site pairs to save memory
	    statistics_.clear();
	  }
          computeStatistics_(site1, site2);
	  return statistics_[key][iter->second];
	}
      } else {
        throw Exception("CoevolutionStatistic. Statistic not found: " + statName);
      }
    }
    double getValue(size_t site1, size_t site2) const {
      return getValue(site1, site2, getMainStatisticsName());
    }
    
    double getRate(size_t site) const { return rates_[site]; }
    const vector<double>& getRates() const { return rates_; }
    double getMinRate(size_t site1, size_t site2) const { return min(rates_[site1], rates_[site2]); }

  protected:
    virtual void computeStatistics_(size_t site1, size_t site2) const = 0;
};

class MiCoevolutionStatistic:
   public virtual AbstractCoevolutionStatistic
{
  public:
    MiCoevolutionStatistic(shared_ptr<const SiteContainer> sites):
        AbstractCoevolutionStatistic(sites, {"MI", "Hjoint", "Hmin"})
    {
      computeEntropies_();
    }
    
    ~MiCoevolutionStatistic() {}
    
    shared_ptr<CoevolutionStatistic> clone(shared_ptr<const SiteContainer> newSites) const {
      return(shared_ptr<MiCoevolutionStatistic>(new MiCoevolutionStatistic(newSites)));
    }

  public:
    string getMainStatisticsName() const { return "MI"; }
    
    void computeStatistics_(size_t site1, size_t site2) const {
      size_t key = site1 + (site2 * nbSites_);
      statistics_[key].resize(3);
      statistics_[key][statNames_["MI"]]     = SiteTools::mutualInformation(sites_->getSite(site1), sites_->getSite(site2), true);
      statistics_[key][statNames_["Hjoint"]] = SiteTools::jointEntropy(sites_->getSite(site1), sites_->getSite(site2), true);
      statistics_[key][statNames_["Hmin"]]   = std::min(rates_[site1], rates_[site2]);
    }

  private:
    //We use the entropy of a site as a measure of its rate:
    void computeEntropies_() {
      for (size_t i = 0; i < nbSites_; ++i) {
        rates_[i] = SiteTools::entropy(sites_->getSite(i), true);
      }
    }
};



class IndexCorrelationCoevolutionStatistic:
   public virtual AbstractCoevolutionStatistic
{
  protected:
    shared_ptr<const AlphabetIndex1> index_;
    VVdouble recodedSites_; //Convert the alignment into a matrix of index values
    char test_; //0: positive, 1: negative, 2 (or anything else): both

  public:
    IndexCorrelationCoevolutionStatistic(shared_ptr<const SiteContainer> sites, shared_ptr<const AlphabetIndex1> index, char test = 0):
        AbstractCoevolutionStatistic(sites, {"Cor", "MinVar"}), index_(index), test_(test)
    {
      if (test != 0 && test != 1) {
        statNames_["Sign"] = 2;
      }
      computeRecodedSites_();
      computeVariances_();
    }
    
    ~IndexCorrelationCoevolutionStatistic() {}
    
    shared_ptr<CoevolutionStatistic> clone(shared_ptr<const SiteContainer> newSites) const {
      return(shared_ptr<IndexCorrelationCoevolutionStatistic>(new IndexCorrelationCoevolutionStatistic(newSites, index_)));
    }

  public:
    string getMainStatisticsName() const { return "Cor"; }
    
    void computeStatistics_(size_t site1, size_t site2) const {
      size_t key = site1 + (site2 * nbSites_);
      double val = VectorTools::cor<double, double>(recodedSites_[site1], recodedSites_[site2]);
      if (isnan(val)) val = 1; //This happens when variances equal 0. Then the proterties are conserved and the sites are perfectly correlated.
      if (test_ == 0) {
        statistics_[key].resize(2);
        statistics_[key][statNames_["Cor"]] = val;
      } else if (test_ == 1) {
        statistics_[key].resize(2);
        statistics_[key][statNames_["Cor"]] = -val;
      } else {
        statistics_[key].resize(3);
        statistics_[key][statNames_["Cor"]] = abs(val);
	int s = 0;
	if (val < 0) s = -1;
	else if (val > 0) s = 1;
        statistics_[key][statNames_["Sign"]] = s;
      }
      statistics_[key][statNames_["MinVar"]] = min(rates_[site1],rates_[site2]);
    }

  private:
    void computeRecodedSites_() {
      size_t nbSeqs = sites_->getNumberOfSequences();
      recodedSites_.resize(nbSites_);
      ApplicationTools::displayTask("Converting Aln. according to ppty");
      ApplicationTools::message->endLine();
      int state;
      vector<int> states;
      double v;
      for (size_t i = 0; i < nbSites_; ++i) {
        ApplicationTools::displayGauge(i + 1, nbSites_);
	recodedSites_[i].resize(nbSeqs);
	for (size_t j = 0; j < nbSeqs; ++j) {
	  state = ((*sites_)(j, i));
	  states = sites_->getAlphabet()->getAlias(state);
          v = 0;
	  for (auto x : states) {
	    v += index_->getIndex(x);
	  }
          recodedSites_[i][j] = v / static_cast<double>(states.size()); //Average over all compatible states
	}
      }
      ApplicationTools::displayTaskDone();
    }

    //We use the variance of a site as a measure of its rate:
    void computeVariances_() {
      for (size_t i = 0; i < nbSites_; ++i) {
        rates_[i] = VectorTools::var<double, double>(recodedSites_[i]);
      }
    }
};


vector<size_t> getRandomIndex(size_t sampleSize, size_t sizeMax) {
  vector<size_t> sample(sampleSize);
  for (size_t i = 0; i < sampleSize; ++i) {
    sample[i] = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(sizeMax);
  }
  return sample;
}


/******************************************************************************/

/**
 * @brief The main function.
 */
int main(int argc, char *argv[])
{
  cout << endl;
  cout << endl;
  cout << "***********************************************************" << endl;
  cout << "* This is Mica+        version 1.5.5b      date: 02/07/20 *" << endl;
  cout << "*                                                         *" << endl;
  cout << "*         Mutual Information Coevolution Analysis         *" << endl;
  cout << "*                 (+ other statistics)                    *" << endl;
  cout << "*                                                         *" << endl;
  cout << "* Author: Julien Y. Dutheil                               *" << endl;
  cout << "***********************************************************" << endl;
  cout << endl;
  
  try
  {
  
  BppApplication mica(argc, argv, "Mica");
  mica.startTimer();

  // **************************
  // * Retrieving parameters: *
  // **************************
  
  if (argc == 1)
  { 
    // No argument, show display some help and leave.
    help();
    exit(-1);
  }

  shared_ptr<Alphabet> alphabet(SequenceApplicationTools::getAlphabet(mica.getParams(), "", false));

  shared_ptr<VectorSiteContainer> allSites(SequenceApplicationTools::getSiteContainer(alphabet.get(), mica.getParams()));
  
  shared_ptr<VectorSiteContainer> sites(SequenceApplicationTools::getSitesToAnalyse(*allSites, mica.getParams()));
  allSites.reset();

  bool removeConst = ApplicationTools::getBooleanParameter("input.remove_const", mica.getParams(), true);
  if (removeConst) {
    size_t n = sites->getNumberOfSites();
    for (size_t i = n; i > 0; --i) {
      if (SiteTools::isConstant(sites->getSite(i - 1), true))
        sites->deleteSite(i - 1);
    }
    ApplicationTools::displayResult("Number of conserved sites ignored", n - sites->getNumberOfSites());
  }

  size_t nbSites = sites->getNumberOfSites();
  size_t nbSeqs  = sites->getNumberOfSequences();
 
  ApplicationTools::displayResult("Number of sequences", nbSeqs);
  ApplicationTools::displayResult("Number of sites", nbSites);
  

  //Shall we use a model?
  shared_ptr<TreeTemplate<Node>>   tree;
  shared_ptr<DRTreeLikelihood>     tl;
  shared_ptr<SubstitutionModel>    model;
  shared_ptr<SubstitutionModelSet> modelSet;
  shared_ptr<DiscreteDistribution> rDist;
  shared_ptr<SubstitutionCount>    simple;
  shared_ptr<vector<double>>       norms;
  bool withModel = ApplicationTools::getBooleanParameter("use_model", mica.getParams(), false, "", true, false);
  
  ApplicationTools::displayBooleanResult("Model of sequence evolution", withModel);
  if (withModel)
  {
    // Get the initial tree
    shared_ptr<Tree> treeTmp(PhylogeneticsApplicationTools::getTree(mica.getParams()));
    tree.reset(new TreeTemplate<Node>(*treeTmp));
    treeTmp.reset();
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));

    string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", mica.getParams(), "no", "", true, false);
    ApplicationTools::displayResult("Heterogeneous model", nhOpt);

    if (nhOpt == "no")
    {  
      model.reset(PhylogeneticsApplicationTools::getSubstitutionModel(alphabet.get(), 0, sites.get(), mica.getParams()));
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist.reset(new ConstantDistribution(1.));
      }
      else
      {
        rDist.reset(PhylogeneticsApplicationTools::getRateDistribution(mica.getParams()));
      }
      tl.reset(new DRHomogeneousTreeLikelihood(*tree, *sites, model.get(), rDist.get(), true, false));
    }
    else if (nhOpt == "one_per_branch")
    {
      model.reset(PhylogeneticsApplicationTools::getSubstitutionModel(alphabet.get(), 0, sites.get(), mica.getParams()));
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist.reset(new ConstantDistribution(1.));
      }
      else
      {
        rDist.reset(PhylogeneticsApplicationTools::getRateDistribution(mica.getParams()));
      }
      vector<double> rateFreqs;
      if (model->getNumberOfStates() != alphabet->getSize())
      {
        //Markov-Modulated Markov Model...
        size_t n = (size_t)(model->getNumberOfStates() / alphabet->getSize());
        rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                     // we should assume a rate distribution for the root also!!!  
      }
      map<string, string> sharedParams;
      shared_ptr<FrequencySet> rootFreqs(PhylogeneticsApplicationTools::getRootFrequencySet(alphabet.get(), 0, sites.get(), mica.getParams(), sharedParams, rateFreqs));
      
      string descGlobal = ApplicationTools::getStringParameter("nonhomogeneous_one_per_branch.shared_parameters", mica.getParams(), "", "", true, 1);

      NestedStringTokenizer nst(descGlobal, "[", "]", ",");
      const deque<string>& descGlobalParameters = nst.getTokens();

      map<string, vector<Vint> > globalParameters;
      for (const auto& desc:descGlobalParameters)
      {
        size_t post = desc.rfind("_");
        if (post == std::string::npos || post == desc.size() - 1 || desc[post + 1] != '[')
          globalParameters[desc]={};
        else
        {
          string key = desc.substr(0,post);
          Vint sint = NumCalcApplicationTools::seqFromString(desc.substr(post + 2, desc.size() - post - 3));
          if (globalParameters.find(key) == globalParameters.end())
            globalParameters[key] = vector<Vint>(1, sint);
          else
            globalParameters[key].push_back(sint);
        }
      }

      for (const auto& globpar:globalParameters)
      {
        ApplicationTools::displayResult("Global parameter", globpar.first);
        if (globpar.second.size()==0)
        {
          string all="All nodes";
          ApplicationTools::displayResult(" shared between nodes", all);
        }
        else
          for (const auto& vint:globpar.second)
            ApplicationTools::displayResult(" shared between nodes", VectorTools::paste(vint,","));
      }

      modelSet.reset(SubstitutionModelSetTools::createNonHomogeneousModelSet(model.get(), rootFreqs.get(), tree.get(), sharedParams, globalParameters)); 
      model.reset(modelSet->getSubstitutionModel(0)->clone());
      tl.reset(new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet.get(), rDist.get(), false));
    }
    else if (nhOpt == "general")
    {
      modelSet.reset(PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet.get(), 0, sites.get(), mica.getParams()));
      if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist.reset(new ConstantDistribution(1.));
      }
      else
      {
        rDist.reset(PhylogeneticsApplicationTools::getRateDistribution(mica.getParams()));
      }
      model.reset(modelSet->getSubstitutionModel(0)->clone());
      tl.reset(new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet.get(), rDist.get(), false));
    }
    else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
    tl->initialize();
 
    double logL = tl->getValue();
    if (std::isinf(logL))
    {
      // This may be due to null branch lengths, leading to null likelihood!
      ApplicationTools::displayWarning("!!! Warning!!! Likelihood is zero.");
      ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
      ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
      ParameterList pl = tl->getBranchLengthsParameters();
      for (size_t i = 0; i < pl.size(); i++)
      {
        if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
      }
      tl->matchParametersValues(pl);
      logL = tl->getValue();
    }
    if (std::isinf(logL))
    {
      ApplicationTools::displayError("!!! Unexpected likelihood == 0.");
      ApplicationTools::displayError("!!! Looking at each site:");
      for (size_t i = 0; i < sites->getNumberOfSites(); i++)
      {
        (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
      }
      ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
      exit(-1);
    }
    tree.reset(new TreeTemplate<Node>(tl->getTree()));

    //Get the substitution mapping in order to compute the rates:
    //TODO: we could use a weighted vector here, eventually...

    TotalSubstitutionRegister* reg = new TotalSubstitutionRegister(model.get());
    simple.reset(new UniformizationSubstitutionCount(model.get(), reg));
    shared_ptr<ProbabilisticSubstitutionMapping> mapping(SubstitutionMappingTools::computeSubstitutionVectors(*tl, *simple, true));
    norms.reset(new vector<double>(computeNorms(*mapping)));
  }
  

  //Get the statistic to use:
  string coevStatDesc = ApplicationTools::getStringParameter("coevolution_statistic", mica.getParams(), "MI", "", true, 0);
  string coevStatName;
  map<string, string> coevStatArgs;
  KeyvalTools::parseProcedure(coevStatDesc, coevStatName, coevStatArgs);
  unique_ptr<CoevolutionStatistic> coevStat;
  ApplicationTools::displayResult("Coevolution statistic:", coevStatName);
  if (coevStatName == "MI") {
    coevStat.reset(new MiCoevolutionStatistic(sites));
  } else if (coevStatName == "Correlation") {
    // This is Neher's approach
    string indexDesc = ApplicationTools::getStringParameter("index", coevStatArgs, "GranthamVolume");
    shared_ptr<const AlphabetIndex1> aaIndex(SequenceApplicationTools::getAlphabetIndex1(alphabet.get(), indexDesc, "Biochemical property:", true));
    
    string altDesc = ApplicationTools::getStringParameter("alternative", coevStatArgs, "both");
    char test;
    if (altDesc == "positive") test = 0;
    else if (altDesc == "negative") test = 1;
    else if (altDesc == "both") test = 2;
    else throw Exception("ERROR: 'alternative' should be either 'positive', 'negative' or 'both'.");

    coevStat.reset(new IndexCorrelationCoevolutionStatistic(sites, aaIndex, test));
  } else {
    throw Exception("Unknown statistic: " + coevStatName);
  }


 
  //Compute all pairwise statistics values:
  
  string path = ApplicationTools::getAFilePath("output.file", mica.getParams(), true, false);
  ApplicationTools::displayResult("Output file", path);

  vector<double> rates = coevStat->getRates();
  vector<double> averageStat;
  ApplicationTools::displayTask("Computing average statistics for every site", true);
  for (size_t i = 0; i < nbSites; ++i) {
    ApplicationTools::displayGauge(i, nbSites - 1);
    double sum = 0;
    for (size_t j = 0; j < nbSites; ++j) {
      if (j != i) {
        sum += coevStat->getValue(i, j);
      }
    }
    averageStat.push_back(sum / static_cast<double>(nbSites - 1));
  }
  ApplicationTools::displayTaskDone();
  double fullAverageStat = VectorTools::mean<double, double>(averageStat);

  //Some variables we'll need:
  const Site *site1, *site2;
  double stat, minRate, apc, rcw, nmin = 0, perm;
  size_t maxNbPermutations = 0;

  //Get the null distribution of statistics:
  string nullMethod = ApplicationTools::getStringParameter("null.method", mica.getParams(), "none", "", true, false);
  ApplicationTools::displayResult("Null distribution", nullMethod);

  bool computePValues = false;
  //These variables will only be used if p-values are computed:
  size_t nbRateClasses = 0;
  Domain* rateDomain = 0;
  vector< vector<double> >* simValues = 0;

  if (nullMethod != "none")
  {
    //Shall we compute built-in p-values?
    if (nullMethod == "z-score")
      computePValues = true;
    else if (nullMethod == "permutations")
      computePValues = false;
    else
      computePValues = ApplicationTools::getBooleanParameter("null.compute_pvalues", mica.getParams(), true);
    if (computePValues)
    {
      nbRateClasses = ApplicationTools::getParameter<size_t>("null.nb_rate_classes", mica.getParams(), 10);
      ApplicationTools::displayResult("Number of sub-distributions", nbRateClasses);
      if (withModel)
        rateDomain = new Domain(0, VectorTools::max(*norms), nbRateClasses);
      else
        rateDomain = new Domain(0, VectorTools::max(rates), nbRateClasses);
      simValues = new vector< vector<double> >(nbRateClasses);
    }

    if (nullMethod == "nonparametric-bootstrap")
    {
      string simpath = ApplicationTools::getAFilePath("null.output.file", mica.getParams(), false, false);
      bool outputToFile = (simpath != "none") && (!TextTools::isEmpty(simpath));
      ApplicationTools::displayTask("Computing null distribution", true);
   
      unique_ptr<ofstream> simout;
      if (outputToFile) {
        ApplicationTools::displayResult("Null output file", simpath);
        simout.reset(new ofstream(simpath.c_str(), ios::out));
	vector<string> statNames = coevStat->getStatisticsNames();
	for (size_t i = 0; i < statNames.size(); ++i) {
          *simout << statNames[i];
	  if (i < statNames.size() - 1)
            *simout << "\t";
	}
        if (withModel)
          *simout << "\tNmin";
        *simout << endl;
      }
    
      size_t nbRepCPU = ApplicationTools::getParameter<size_t>("null.nb_rep_CPU", mica.getParams(), 10);
      size_t nbRepRAM = ApplicationTools::getParameter<size_t>("null.nb_rep_RAM", mica.getParams(), 100);
  
      shared_ptr<OutputStream> os = ApplicationTools::warning;
      ApplicationTools::warning.reset();
      for (size_t i = 0; i < nbRepCPU; i++)
      {
        //Generate data set:
        vector<size_t> index1 = getRandomIndex(nbRepRAM, nbSites);
        vector<size_t> index2 = getRandomIndex(nbRepRAM, nbSites);
  
        for (size_t j = 0; j < nbRepRAM; j++)
        {
          ApplicationTools::displayGauge(i * nbRepRAM + j, nbRepCPU * nbRepRAM - 1, '>');
          if (outputToFile) {
	    for (size_t k = 0; k < coevStat->getNumberOfStatistics(); ++k) {
              string statName = coevStat->getStatisticsNames()[k];
              stat = coevStat->getValue(index1[j], index2[j], statName);
              *simout << stat;
	      if (k != coevStat->getNumberOfStatistics())
	        *simout << "\t";
	    }
	  }
	  stat = coevStat->getValue(index1[j], index2[j]);
          if (withModel) {
            nmin = min((*norms)[index1[j]], (*norms)[index2[j]]);
            if (outputToFile)
              *simout << "\t" << nmin;
            if (computePValues) {
              try {
                size_t cat = rateDomain->getIndex(nmin);
                (*simValues)[cat].push_back(stat);
              } catch (OutOfRangeException& oore) {}
            }
          } else {
            if (computePValues) {
              try {
		minRate = coevStat->getMinRate(index1[j], index2[j]);
                size_t cat = rateDomain->getIndex(minRate);
                (*simValues)[cat].push_back(stat);
              } catch (OutOfRangeException& oore) {}
            }
          }
          if (outputToFile)
            *simout << endl;
        }
      }
      ApplicationTools::warning = os;
  
      if (outputToFile) {
        simout->close();
      }
      ApplicationTools::displayTaskDone();  
    }
    else if (nullMethod == "parametric-bootstrap")
    {
      if (!withModel)
        throw Exception("You need to specify a model of sequence evolution in order to use a parametric bootstrap approach!");
      bool continuousSim = ApplicationTools::getBooleanParameter("simulations.continuous", mica.getParams(), false, "", true, false);
      ApplicationTools::displayResult("Rate distribution for simulations", (continuousSim ? "continuous" : "discrete"));
      shared_ptr<SequenceSimulator> simulator;
      if (modelSet)
      {
        simulator.reset(new NonHomogeneousSequenceSimulator(modelSet.get(), rDist.get(), tree.get()));
        dynamic_cast<NonHomogeneousSequenceSimulator *>(simulator.get())->enableContinuousRates(continuousSim);
      }
      else
      {
        simulator.reset(new HomogeneousSequenceSimulator(model.get(), rDist.get(), tree.get()));
        dynamic_cast<HomogeneousSequenceSimulator *>(simulator.get())->enableContinuousRates(continuousSim);
      }
   
      string simpath = ApplicationTools::getAFilePath("null.output.file", mica.getParams(), false, false);
      bool outputToFile = (simpath != "none") && (!TextTools::isEmpty(simpath));
      ApplicationTools::displayTask("Computing null distribution", true);
 
      unique_ptr<ofstream> simout;
      if (outputToFile) {
        ApplicationTools::displayResult("Null output file", simpath);
        simout.reset(new ofstream(simpath.c_str(), ios::out));
	vector<string> statNames = coevStat->getStatisticsNames();
	for (size_t i = 0; i < statNames.size(); ++i) {
          *simout << statNames[i];
	  if (i < statNames.size() - 1)
            *simout << "\t";
	}
        *simout << "\tNmin" << endl;
      }

      size_t nbRepCPU = ApplicationTools::getParameter<size_t>("null.nb_rep_CPU", mica.getParams(), 10);
      size_t nbRepRAM = ApplicationTools::getParameter<size_t>("null.nb_rep_RAM", mica.getParams(), 100);
  
      shared_ptr<OutputStream> os = ApplicationTools::warning;
      ApplicationTools::warning.reset();
      for (size_t i = 0; i < nbRepCPU; i++)
      {
        //Generate data set. We simulate two times nbRepRAM and study all pairs (j, nbRepRAM +j), with 0 < j < nbRepRAM.
        shared_ptr<SiteContainer> simSites(simulator->simulate(2 * nbRepRAM));
        tl->setData(*simSites);
        tl->initialize();
        unique_ptr<ProbabilisticSubstitutionMapping> simMapping(SubstitutionMappingTools::computeSubstitutionVectors(*tl, *simple, false));
        vector<double> simNorms = computeNorms(*simMapping);
        shared_ptr<CoevolutionStatistic> simStat = coevStat->clone(simSites);
        for (size_t j = 0; j < nbRepRAM; j++)
        {
          ApplicationTools::displayGauge(i * nbRepRAM + j, nbRepCPU * nbRepRAM - 1, '>');

          if (outputToFile) {
	    for (size_t k = 0; k < simStat->getNumberOfStatistics(); ++k) {
              string statName = simStat->getStatisticsNames()[k];
              stat = simStat->getValue(j, nbRepRAM + j, statName);
              *simout << stat;
	      if (k != simStat->getNumberOfStatistics())
	        *simout << "\t";
	    }
	  }
	  stat = simStat->getValue(j, nbRepRAM + j);
          nmin = min(simNorms[j], simNorms[nbRepRAM + j]);
          if (outputToFile)
            *simout << "\t" << nmin << endl;
  
          if (computePValues) {
            try {
              size_t cat = rateDomain->getIndex(nmin);
              (*simValues)[cat].push_back(stat);
            } catch(OutOfRangeException& oore) {}
          }

        }
  
      }
      ApplicationTools::warning = os;
  
      if (outputToFile)
        simout->close();
      ApplicationTools::displayTaskDone(); 
    }
    else if (nullMethod == "z-score")
    {
      string zScoreStat = ApplicationTools::getStringParameter("null.method_zscore.stat", mica.getParams(), "MIp", "", true, false);
      ApplicationTools::displayResult("Compute p-value for", zScoreStat);
      short correctStat = 0;
      if (zScoreStat == "MIp" || zScoreStat == "APC")
        correctStat = 1;
      else if (zScoreStat == "MIc" || zScoreStat == "RWC")
        correctStat = 2;
      else if (zScoreStat != "MI" || zScoreStat == "Raw")
        throw Exception("Unkown statistic, should be 'MI'/'Raw', 'MIp'/'APC' or 'MIc'/'RWC'.");
      ApplicationTools::displayTask("Computing total distribution", true);
      shared_ptr<OutputStream> os = ApplicationTools::warning;
      ApplicationTools::warning.reset();
      size_t c = 0;
      for (size_t i = 0; i < nbSites - 1; i++)
      {
        for (size_t j = i + 1; j < nbSites; j++)
        {
          ApplicationTools::displayGauge(c++, nbSites * (nbSites - 1) / 2 - 1, '>');
          stat = coevStat->getValue(i, j);
  
          if (withModel) {
            minRate  = min((*norms)[i], (*norms)[j]);
	  } else {
	    minRate = coevStat->getMinRate(i, j);
	  }
          try {
            size_t cat = rateDomain->getIndex(minRate);
            if (correctStat == 1) {
              apc = averageStat[i] * averageStat[j] / fullAverageStat;
              (*simValues)[cat].push_back(stat - apc);
            } else if (correctStat == 2) {
              rcw = averageStat[i] * averageStat[j] / 2.;
              (*simValues)[cat].push_back(stat / rcw);
            } else {
              (*simValues)[cat].push_back(stat);
            }
          } catch (OutOfRangeException& oore) {}
        }
      }
      ApplicationTools::warning = os;
      ApplicationTools::displayTaskDone();  
    }
    else if (nullMethod == "permutations")
    {
      maxNbPermutations = ApplicationTools::getParameter<size_t>("null.max_number_of_permutations", mica.getParams(), 1000);
      if (maxNbPermutations == 0)
      {
        throw Exception("Permutation number should be greater than 0!");
      }
      else
      {
        ApplicationTools::displayResult("Maximum number of permutations", maxNbPermutations);
      }
    } 
    else
      throw Exception("Unvalid null distribution method specified: " + nullMethod);
  }

  //We need to sort observations for a better efficiency:
  if (computePValues) {
    for (size_t i = 0; i < nbRateClasses; ++i) {
      sort((*simValues)[i].begin(), (*simValues)[i].end());
    }
  }

  //here comes the real stuff:
  ApplicationTools::displayTask("Computing all statistics", true);

  ofstream out(path.c_str(), ios::out);
  out << "Group";

  for (auto statName: coevStat->getStatisticsNames()) {
    out << "\t" << statName;
  }
  out << "\tAPC\tRCW";
  if (withModel)
    out << "\tNmin";
  if (maxNbPermutations > 0)
    out << "\tPerm.p.value\tPerm.nb";
  if (computePValues)
    out << "\tBs.p.value\tBs.nb";
  out << endl;

  size_t nbPerm;
  size_t c = 0;
  for (size_t i = 0; i < nbSites - 1; ++i)
  {
    site1 = &sites->getSite(i);
    for (size_t j = i + 1; j < nbSites; ++j)
    {
      ApplicationTools::displayGauge(c++, (nbSites - 1) * (nbSites) / 2 - 1);
      site2 = &sites->getSite(j);
      miTest(*site1, *site2, maxNbPermutations, stat, perm, nbPerm);

      out << "[" << site1->getPosition() << ";" << site2->getPosition() << "]";
      for (auto statName: coevStat->getStatisticsNames()) {
        out << "\t" << coevStat->getValue(i, j, statName);
      }
      apc = averageStat[i] * averageStat[j] / fullAverageStat;
      rcw = averageStat[i] * averageStat[j] / 2.;
      out << "\t" << apc << "\t" << rcw;
      
      if (withModel) {
        nmin = min((*norms)[i], (*norms)[j]);
        out << "\t" << nmin;
	minRate = nmin;
      } else {
        minRate = coevStat->getMinRate(i, j);
      }
      //Permutations p-value:
      if (maxNbPermutations > 0)
        out << "\t" << perm << "\t" << nbPerm;
      
      //P-values:
      if (computePValues) {
        //Bootstrap:
        try {
          size_t cat = rateDomain->getIndex(minRate);
          size_t nsim = (*simValues)[cat].size();
          size_t count;
          for (count = 0; count < nsim && (*simValues)[cat][count] < stat; ++count) {}
          double pvalue = static_cast<double>(nsim - count + 1) / static_cast<double>(nsim + 1);
          out << "\t" << pvalue << "\t" << nsim;
        } catch (OutOfRangeException& oore) {
          out << "\tNA\t0";
        }
        //MIp p-value:
        //TODO
      }
      
      out << endl;
    }
  }
  out.close();
  ApplicationTools::displayTaskDone();

  if (computePValues) {
    delete simValues;
    delete rateDomain;
  }

  mica.done();
  }
  catch (exception & e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}
