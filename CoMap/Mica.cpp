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
#include <Bpp/Utils/AttributesTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
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
#include <Bpp/Phyl/Legacy/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Legacy/Likelihood/NonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/DRHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/RASTools.h>
#include <Bpp/Phyl/Legacy/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Simulation/MutationProcess.h>
#include <Bpp/Phyl/Legacy/Simulation/NonHomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Mapping/UniformizationSubstitutionCount.h>
#include <Bpp/Phyl/Legacy/Mapping/ProbabilisticSubstitutionMapping.h>
#include <Bpp/Phyl/Legacy/Mapping/SubstitutionMappingTools.h>

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


vector<double> computeNorms(const LegacyProbabilisticSubstitutionMapping& mapping)
{
  size_t nbVectors = mapping.getNumberOfSites();
  vector<double> vect(nbVectors);
  for (size_t i = 0; i < nbVectors; i++)
    vect[i] = LegacySubstitutionMappingTools::computeNormForSite(mapping, i);
  return vect;
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
  cout << "* This is Mica         version 1.6.0a      date: 01/11/23 *" << endl;
  cout << "*         Mutual Information Coevolution Analysis         *" << endl;
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

  auto allSites = SequenceApplicationTools::getSiteContainer(alphabet, mica.getParams());
  
  shared_ptr<SiteContainerInterface> sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, mica.getParams());
  allSites.reset();

  bool removeConst = ApplicationTools::getBooleanParameter("input.remove_const", mica.getParams(), true);
  if (removeConst) {
    size_t n = sites->getNumberOfSites();
    for (size_t i = n; i > 0; --i) {
      if (SiteTools::isConstant(sites->site(i - 1), true))
        sites->deleteSite(i - 1);
    }
    ApplicationTools::displayResult("Number of conserved sites ignored", n - sites->getNumberOfSites());
  }

  size_t nbSites = sites->getNumberOfSites();
  size_t nbSeqs  = sites->getNumberOfSequences();
 
  ApplicationTools::displayResult("Number of sequences", nbSeqs);
  ApplicationTools::displayResult("Number of sites", nbSites);
  

  //Shall we use a model?
  shared_ptr<TreeTemplate<Node>>            tree;
  shared_ptr<DRTreeLikelihoodInterface>     tl;
  shared_ptr<SubstitutionModelInterface>    model;
  shared_ptr<SubstitutionModelSet>          modelSet;
  shared_ptr<DiscreteDistributionInterface> rDist;
  shared_ptr<SubstitutionCountInterface>    simple;
  shared_ptr<vector<double>>                norms;
  bool withModel = ApplicationTools::getBooleanParameter("use_model", mica.getParams(), false, "", true, false);
  
  ApplicationTools::displayBooleanResult("Model of sequence evolution", withModel);
  if (withModel)
  {
    // Get the initial tree
    shared_ptr<Tree> treeTmp = PhylogeneticsApplicationTools::getTree(mica.getParams());
    tree.reset(new TreeTemplate<Node>(*treeTmp));
    treeTmp.reset();
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));

    string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", mica.getParams(), "no", "", true, false);
    ApplicationTools::displayResult("Heterogeneous model", nhOpt);

    if (nhOpt == "no")
    {
      map<string, string> sharedParams;
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, nullptr, sites, mica.getParams(), sharedParams);
      if (model->getNumberOfStates() > model->alphabet().getSize())
      {
        //Markov-modulated Markov model!
        rDist.reset(new ConstantDistribution(1.));
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mica.getParams());
      }
      tl.reset(new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, false));
    }
    else if (nhOpt == "one_per_branch")
    {
      map<string, string> sharedParams;
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, nullptr, sites, mica.getParams(), sharedParams);
      if (model->getNumberOfStates() > model->alphabet().getSize())
      {
        //Markov-modulated Markov model!
        rDist.reset(new ConstantDistribution(1.));
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mica.getParams());
      }
      vector<double> rateFreqs;
      if (model->getNumberOfStates() != alphabet->getSize())
      {
        //Markov-Modulated Markov Model...
        size_t n =(size_t)(model->getNumberOfStates() / alphabet->getSize());
        rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                     // we should assume a rate distribution for the root also!!!  
      }
      shared_ptr<FrequencySetInterface> rootFreqs = PhylogeneticsApplicationTools::getRootFrequencySet(alphabet, nullptr, *sites, mica.getParams(), sharedParams, rateFreqs);
      
      string descGlobal = ApplicationTools::getStringParameter("nonhomogeneous_one_per_branch.shared_parameters", mica.getParams(), "", "", true, 1);

      NestedStringTokenizer nst(descGlobal, "[", "]", ",");
      const deque<string>& descGlobalParameters = nst.getTokens();

      map<string, vector<Vint>> globalParameters;
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

      modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, *tree, sharedParams, globalParameters); 
      model.reset(modelSet->substitutionModel(0).clone());
      tl.reset(new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false));
    }
    else if (nhOpt == "general")
    {
      modelSet = LegacyPhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, nullptr, sites, mica.getParams());
      if (modelSet->getNumberOfStates() > modelSet->alphabet().getSize())
      {
        //Markov-modulated Markov model!
        rDist.reset(new ConstantDistribution(1.));
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mica.getParams());
      }
      model.reset(modelSet->substitutionModel(0).clone());
      tl.reset(new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false));
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
        (*ApplicationTools::error << "Site " << sites->site(i).getCoordinate() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
      }
      ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
      exit(-1);
    }
    tree.reset(new TreeTemplate<Node>(tl->tree()));

    //Get the substitution mapping in order to compute the rates:

    auto reg = make_shared<TotalSubstitutionRegister>(model->getStateMap());
    simple.reset(new UniformizationSubstitutionCount(model, reg));
    auto mapping = LegacySubstitutionMappingTools::computeSubstitutionVectors(tl, simple, true);
    norms.reset(new vector<double>(computeNorms(*mapping)));
  }
  
  //Compute all pairwise MI values:
  
  string path = ApplicationTools::getAFilePath("output.file", mica.getParams(), true, false);
  ApplicationTools::displayResult("Output file", path);

  vector<double> averageMI;
  vector<double> entropy;
  ApplicationTools::displayTask("Computing average MIs", true);
  for (size_t i = 0; i < nbSites; ++i) {
    auto& site1 =  sites->site(i);
    ApplicationTools::displayGauge(i, nbSites - 1);
    double sum = 0;
    for (size_t j = 0; j < nbSites; ++j) {
      if (j != i) {
        auto& site2 =  sites->site(j);
        sum += SiteTools::mutualInformation(site1, site2, true); 
      }
    }
    entropy.push_back(SiteTools::entropy(site1, true));
    averageMI.push_back(sum / static_cast<double>(nbSites - 1));
  }
  ApplicationTools::displayTaskDone();
  double fullAverageMI = VectorTools::mean<double, double>(averageMI);

  //Some variables we'll need:
  double stat, apc, rcw, nmin = 0, hj, hm, perm;
  size_t maxNbPermutations = 0;

  //Get the null distribution of MI values:
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
        rateDomain = new Domain(0, VectorTools::max(entropy), nbRateClasses);
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
        *simout << "MI\tHjoint\tHmin";
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
        auto index1 = make_shared<vector<size_t>>();
         auto sites1 = SiteContainerTools::sampleSites(*sites, nbRepRAM, index1);
  
        auto index2 = make_shared<vector<size_t>>();
        auto sites2 = SiteContainerTools::sampleSites(*sites, nbRepRAM, index2);
  
        for (size_t j = 0; j < nbRepRAM; j++)
        {
          ApplicationTools::displayGauge(i * nbRepRAM + j, nbRepCPU * nbRepRAM - 1, '>');
          auto& site1 = sites1->site(j);
          auto& site2 = sites2->site(j);
          stat  = SiteTools::mutualInformation(site1, site2, true);
          hj    = SiteTools::jointEntropy(site1, site2, true);
          hm    = std::min(entropy[(*index1)[j]], entropy[(*index2)[j]]);
  
          if (outputToFile)
            *simout << stat << "\t" << hj << "\t" << hm;
          if (withModel) {
            nmin  = min((*norms)[(*index1)[j]], (*norms)[(*index2)[j]]);
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
                size_t cat = rateDomain->getIndex(hm);
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
      shared_ptr<NonHomogeneousSequenceSimulator> simulator;
      if (modelSet)
      {
        simulator.reset(new NonHomogeneousSequenceSimulator(modelSet, rDist, tree));
        simulator->enableContinuousRates(continuousSim);
      }
      else
      {
        simulator.reset(new NonHomogeneousSequenceSimulator(model, rDist, tree));
        simulator->enableContinuousRates(continuousSim);
      }
   
      string simpath = ApplicationTools::getAFilePath("null.output.file", mica.getParams(), false, false);
      bool outputToFile = (simpath != "none") && (!TextTools::isEmpty(simpath));
      ApplicationTools::displayTask("Computing null distribution", true);
 
      unique_ptr<ofstream> simout;
      if (outputToFile) {
        ApplicationTools::displayResult("Null output file", simpath);
        simout.reset(new ofstream(simpath.c_str(), ios::out));
    
        *simout << "MI\tHjoint\tHmin\tNmin" << endl;
      }

      size_t nbRepCPU = ApplicationTools::getParameter<size_t>("null.nb_rep_CPU", mica.getParams(), 10);
      size_t nbRepRAM = ApplicationTools::getParameter<size_t>("null.nb_rep_RAM", mica.getParams(), 100);
  
      shared_ptr<OutputStream> os = ApplicationTools::warning;
      ApplicationTools::warning.reset();
      for (size_t i = 0; i < nbRepCPU; i++)
      {
        //Generate data set:
        auto sites1 = simulator->simulate(nbRepRAM);
        tl->setData(*sites1);
        tl->initialize();
        vector<double> norms1;
        auto mapping1 = LegacySubstitutionMappingTools::computeSubstitutionVectors(tl, simple, false);
        norms1 = computeNorms(*mapping1);
  
        auto sites2 = simulator->simulate(nbRepRAM);
        tl->setData(*sites2);
        tl->initialize();
        vector<double> norms2;
        auto mapping2 = LegacySubstitutionMappingTools::computeSubstitutionVectors(tl, simple, false);
        norms2 = computeNorms(*mapping2);
  
        for (size_t j = 0; j < nbRepRAM; j++)
        {
          ApplicationTools::displayGauge(i * nbRepRAM + j, nbRepCPU * nbRepRAM - 1, '>');
          auto& site1 = sites1->site(j);
          auto& site2 = sites2->site(j);
          stat  = SiteTools::mutualInformation(site1, site2, true);
          hj    = SiteTools::jointEntropy(site1, site2, true);
          hm    = std::min(entropy[i], entropy[j]);
          nmin  = min(norms1[j], norms2[j]);
  
          if (computePValues) {
            try {
              size_t cat = rateDomain->getIndex(nmin);
              (*simValues)[cat].push_back(stat);
            } catch(OutOfRangeException& oore) {}
          }

          if (outputToFile)
            *simout << stat << "\t" << hj << "\t" << hm << "\t" << nmin << endl;
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
      if (zScoreStat == "MIp")
        correctStat = 1;
      else if (zScoreStat == "MIc")
        correctStat = 2;
      else if (zScoreStat != "MI")
        throw Exception("Unkown statistic, should be 'MI', 'MIp' or 'MIc'.");
      ApplicationTools::displayTask("Computing total distribution", true);
      shared_ptr<OutputStream> os = ApplicationTools::warning;
      ApplicationTools::warning.reset();
      size_t c = 0;
      for (size_t i = 0; i < nbSites - 1; i++)
      {
        auto& site1 = sites->site(i);
        for (size_t j = i + 1; j < nbSites; j++)
        {
          auto& site2 = sites->site(j);
          ApplicationTools::displayGauge(c++, nbSites * (nbSites - 1) / 2 - 1, '>');
          stat  = SiteTools::mutualInformation(site1, site2, true);
          hj    = SiteTools::jointEntropy(site1, site2, true);
          hm    = std::min(entropy[i], entropy[j]);
  
          if (withModel) {
            nmin  = min((*norms)[i], (*norms)[j]);
            try {
              size_t cat = rateDomain->getIndex(nmin);
              if (correctStat == 1) {
                apc = averageMI[i] * averageMI[j] / fullAverageMI;
                (*simValues)[cat].push_back(stat - apc);
              } else if (correctStat == 2) {
                rcw = averageMI[i] * averageMI[j] / 2.;
                (*simValues)[cat].push_back(stat / rcw);
              } else {
                (*simValues)[cat].push_back(stat);
              }
            } catch (OutOfRangeException& oore) {}
          } else {
            try {
              size_t cat = rateDomain->getIndex(hm);
              if (correctStat == 1) {
                apc = averageMI[i] * averageMI[j] / fullAverageMI;
                (*simValues)[cat].push_back(stat - apc);
              } else if (correctStat == 2) {
                rcw = averageMI[i] * averageMI[j] / 2.;
                (*simValues)[cat].push_back(stat / rcw);
              } else {
                (*simValues)[cat].push_back(stat);
              }
            } catch (OutOfRangeException& oore) {}
          }
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
  ApplicationTools::displayTask("Computing all MI scores", true);

  ofstream out(path.c_str(), ios::out);
  out << "Group\tMI\tAPC\tRCW\tHjoint\tHmin";
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
    auto& site1 = sites->site(i);
    for (size_t j = i + 1; j < nbSites; ++j)
    {
      ApplicationTools::displayGauge(c++, (nbSites - 1) * (nbSites) / 2 - 1);
      auto& site2 = sites->site(j);
      miTest(site1, site2, maxNbPermutations, stat, perm, nbPerm);

      out << "[" << site1.getCoordinate() << ";" << site2.getCoordinate() << "]\t" << stat << "\t";
      apc = averageMI[i] * averageMI[j] / fullAverageMI;
      rcw = averageMI[i] * averageMI[j] / 2.;
      hj  = SiteTools::jointEntropy(site1, site2, true);
      hm  = std::min(entropy[i], entropy[j]);
      out << apc << "\t" << rcw << "\t" << hj << "\t" << hm;
      
      if (withModel) {
        nmin = min((*norms)[i], (*norms)[j]);
        out << "\t" << nmin;
      }
      //Permutations p-value:
      if (maxNbPermutations > 0)
        out << "\t" << perm << "\t" << nbPerm;
      
      //P-values:
      if (computePValues) {
        //Bootstrap:
        try {
          size_t cat = withModel ? rateDomain->getIndex(nmin) : rateDomain->getIndex(hm);
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
