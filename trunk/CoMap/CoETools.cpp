//
// File: CoETools.cpp
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

#include "CoETools.h"
#include "AnalysisTools.h"
#include "Domain.h"
#include "IntervalData.h"

#include <Bpp/Text.all>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/BppString.h>
#include <Bpp/Numeric/ParameterExceptions.h>

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/AlphabetIndex.all>

// From PhylLib:
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Mapping.all>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>

using namespace bpp;

// From the STL:
#include <fstream>
#include <iomanip>
using namespace std;

/******************************************************************************/

void CoETools::readData(
  TreeTemplate<Node>*          &tree,
  Alphabet*                    &alphabet,
  VectorSiteContainer*         &allSites,
  VectorSiteContainer*         &sites,
  SubstitutionModel*           &model,
  SubstitutionModelSet*        &modelSet,
  DiscreteDistribution*        &rDist,
  DRTreeLikelihood*            &tl,
  map<string, string>          &params,
  const string                 &suffix)
{
  alphabet = SequenceApplicationTools::getAlphabet(params, suffix, true);
  allSites = SequenceApplicationTools::getSiteContainer(alphabet, params, suffix, false);
  sites    = SequenceApplicationTools::getSitesToAnalyse(*allSites, params, suffix, true, true, true);

  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", params, "no", "", true, 1);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  model = 0;
  modelSet = 0;

  if (nhOpt == "no")
  {  
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, 0, sites, params);
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      //rDist = new ConstantDistribution(1.);
      throw Exception("Covarion models not supported for now :(");
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    }
    tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, false);
  }
  else if (nhOpt == "one_per_branch")
  {
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, 0, sites, params);
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      //rDist = new ConstantDistribution(1.);
      throw Exception("Covarion models not supported for now :(");
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    }
    vector<double> rateFreqs;
    if (model->getNumberOfStates() != alphabet->getSize())
    {
      //Markov-Modulated Markov Model...
      size_t n =static_cast<size_t>(model->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1./static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                    // we should assume a rate distribution for the root also!!!  
    }
    FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, 0, sites, params, rateFreqs);
    vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", params, ',', "", "", true, 1);
    modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
    model = modelSet->getModel(0)->clone();
    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false);
  }
  else if (nhOpt == "general")
  {
    modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, 0, 0, params);
    if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      //rDist = new ConstantDistribution(1.);
      throw Exception("Covarion models not supported for now :(");
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    }
    model = modelSet->getModel(0)->clone();
    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false); 
  }
  else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
  
  tl->initialize();

  ApplicationTools::displayTask("Tree likelihood");
  double ll = tl->getValue();
  if (isinf(ll))
  {
    ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
    ApplicationTools::displayError("!!! You should consider reestimating all branch lengths parameters.");
    ApplicationTools::displayError("!!! Site-specific likelihood have been written in file DEBUG_likelihoods.txt .");
    ofstream debug ("DEBUG_likelihoods.txt", ios::out);
    for (size_t i = 0; i < sites->getNumberOfSites(); i++)
    {
      debug << "Position " << sites->getSite(i).getPosition() << " = " << tl->getLogLikelihoodForASite(i) << endl; 
    }
    debug.close();
    exit(-1);
  }
  (ApplicationTools::message->setPrecision(20) << ll).endLine();
  
  string optimization = ApplicationTools::getStringParameter("optimization", params, "None", suffix, true, false);
  if (optimization != "None")
  {
    ApplicationTools::displayResult("Optimization", string("yes"));
    PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), params, suffix, true, true);
    TreeTemplate<Node> *treeOpt = new TreeTemplate<Node>(tl->getTree());
    PhylogeneticsApplicationTools::writeTree(*treeOpt, params, "output.", suffix);
  
    // Actualize tree.
    // Substitution model and rate distribution are actualized automatically by
    // the likelihood function.
    delete tree;
    tree = treeOpt;
    
    // Print parameters:
    ApplicationTools::displayResult("Final likelihood, -lnL =", TextTools::toString(tl -> getLogLikelihood(), 20));
  }
  string tags = ApplicationTools::getAFilePath("output.tags.file", params, false, false, suffix, false, "none", 2);
  if (tags != "none")
  {
    TreeTemplate<Node> treeCopy(*tree);
    vector<Node *> nodes = treeCopy.getInnerNodes();
    for (unsigned int i = 0; i < nodes.size(); i++) {
      nodes[i]->setNodeProperty("name", BppString(TextTools::toString(nodes[i]->getId())));
    }
    nodes = treeCopy.getLeaves();

    // Writing leaf names translation:
    string tlnPath = ApplicationTools::getAFilePath("output.tags.translation", params, false, false, suffix, false);
    if (tlnPath != "none") {
      ofstream tln(tlnPath.c_str(), ios::out);
      tln << "Name\tId" << endl;
      for (unsigned int i = 0; i < nodes.size(); i++)
      {
        tln << nodes[i]->getName() << "\t" << nodes[i]->getId() << endl;
      }
      tln.close();
    }

    // Translate names:
    for (unsigned int i = 0; i < nodes.size(); ++i) {
      nodes[i]->setName(TextTools::toString(nodes[i]->getId()));
    }
    Newick newick;
    newick.enableExtendedBootstrapProperty("name");
    newick.write(treeCopy, tags, true);
  }

  bool removeConst = ApplicationTools::getBooleanParameter("input.remove_const", params, true, "", false, 1);
  if (removeConst) {
    ApplicationTools::displayTask("Remove conserved positions", true);
    size_t n = sites->getNumberOfSites();
    for (size_t i = n; i > 0; --i) {
      ApplicationTools::displayGauge(n - i, n - 1, '=');
      if (SiteTools::isConstant(sites->getSite(i - 1), true))
        sites->deleteSite(i - 1);
    }
    ApplicationTools::message->endLine();
    ApplicationTools::displayResult("Number of conserved sites ignored", n - sites->getNumberOfSites());
    tl->setData(*sites);
    tl->initialize();
  }

}
  
/******************************************************************************/
  
ProbabilisticSubstitutionMapping* CoETools::getVectors(
  const DRTreeLikelihood& drtl,
  SubstitutionCount     & substitutionCount,
  const SiteContainer   & completeSites,
  map<string, string>   & params,
  const string          & suffix)
{
  ProbabilisticSubstitutionMapping* substitutions = 0;
  string inputVectorsFilePath = ApplicationTools::getAFilePath("input.vectors.file", params, false, true, suffix, false, "none", 1);

  if (inputVectorsFilePath != "none")
  {
    //We try to load the substitution vector directly from file:
    size_t nbSites = drtl.getNumberOfSites();
    ApplicationTools::displayResult("Substitution mapping in file:", inputVectorsFilePath);
    ifstream sc(inputVectorsFilePath.c_str(), ios::in);
    substitutions = new ProbabilisticSubstitutionMapping(drtl.getTree(), &substitutionCount, nbSites);
    SubstitutionMappingTools::readFromStream(sc, *substitutions, 0);
  }
  else
  {
    //We compute the substitutions vector:

    string outputVectorsFilePath = ApplicationTools::getAFilePath("output.vectors.file", params, false, false, suffix, false, "none", 1);
    ApplicationTools::displayResult("Output mapping to file" + suffix, outputVectorsFilePath);

    bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true, "", true, 4); //These two options are realy for benchmarking only!
    bool joint   = ApplicationTools::getBooleanParameter("nijt.joint"  , params, true, "", true, 4);
    if (average) {
      if (joint) {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectors(drtl, substitutionCount, 0);
      } else {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl, substitutionCount, 0);
      }
    } else {
      if (joint) {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl, substitutionCount, 0);
      } else {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl, substitutionCount, 0);
      }
    }
    if (outputVectorsFilePath != "none" && outputVectorsFilePath != "None") {
      ofstream outputVectors(outputVectorsFilePath.c_str(), ios::out);
      SubstitutionMappingTools::writeToStream(*substitutions, completeSites, 0, outputVectors);
      outputVectors.close();
    }

  }
  return substitutions;
}

/******************************************************************************/

int CoETools::getMinRateClass(map<string, string>& params, string suffix)
{
  int minRateClass = ApplicationTools::getIntParameter("statistic.min_rate_class", params, 0, suffix, true, 2);
  if (minRateClass > 0)
    ApplicationTools::displayMessage(
        "Only sites with posterior rate class >= " +
        TextTools::toString(minRateClass) +
        " will be compared.");
  return minRateClass;  
}

/******************************************************************************/

double CoETools::getMinRate(map<string, string>& params, string suffix)
{
  double minRate = ApplicationTools::getDoubleParameter("statistic.min_rate", params, 0., suffix, true, 2);
  if (minRate > 0.)
    ApplicationTools::displayMessage(
        "Only sites with posterior rate > = " +
        TextTools::toString(minRate) +
        " will be compared.");
  return minRate;
}

/******************************************************************************/

int CoETools::getMaxRateClassDiff(map<string, string>& params)
{
  int maxRateClassDiff = ApplicationTools::getIntParameter("statistic.max_rate_class_diff", params, -1, "", false, 2);
  if (maxRateClassDiff >= 0) 
    ApplicationTools::displayMessage(
      "Only pairs of sites with difference in posterior rate class <= " +
      TextTools::toString(maxRateClassDiff) +
      " will be compared.");
  return maxRateClassDiff;  
}

/******************************************************************************/

double CoETools::getMaxRateDiff(map<string, string>& params)
{
  double maxRateDiff = ApplicationTools::getDoubleParameter("statistic.max_rate_diff", params, -1., "", false, 2);
  if (maxRateDiff >= 0.)
    ApplicationTools::displayMessage(
        "Only pairs of sites with difference in posterior rate <= " +
        TextTools::toString(maxRateDiff) +
        " will be compared.");
  return maxRateDiff;
}

/******************************************************************************/

double CoETools::getStatisticMin(map<string, string>& params)
{
  double minStatistic = ApplicationTools::getDoubleParameter("statistic.min", params, 0, "", false, 2);
  if (minStatistic > 0)
    ApplicationTools::displayMessage(
      "Only pairs of sites with abs(statistic) >= " +
      TextTools::toString(minStatistic) +
      " will be written.");
  return minStatistic;
}

/******************************************************************************/

bool CoETools::haveToPerformIndependantComparisons(map<string, string>& params)
{
  bool indepComp = ApplicationTools::getBooleanParameter("independant_comparisons", params, false, "", false, false);
  if (indepComp)
    ApplicationTools::displayMessage(
      "Only independant comparisons will be performed.");
  return indepComp;
}

/******************************************************************************/

void CoETools::writeInfos(
  const SiteContainer & completeSites,
  const DiscreteRatesAcrossSitesTreeLikelihood & ras,
  map<string, string> & params,
  const string & suffix)
{
  string outFile = ApplicationTools::getAFilePath("output.infos", params, false, false, suffix, true, "none", 1);
  if (outFile == "none") return;

  // Get the rate class with maximum posterior probability:
  vector<size_t> classes = ras.getRateClassWithMaxPostProbOfEachSite();
  // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
  Vdouble rates = ras.getPosteriorRateOfEachSite();
  Vdouble logLn = ras.getLogLikelihoodForEachSite();

  ApplicationTools::displayResult("Alignment information logfile", outFile);

  ofstream out(outFile.c_str(), ios::out);
  out << "Group\tIsComplete\tIsConstant\tRC\tPR\tlogLn" << endl;

  for (unsigned int i = 0; i < completeSites.getNumberOfSites(); i++)
  {
    const Site * currentSite = &completeSites.getSite(i);
    int currentSitePosition = currentSite->getPosition();
    int isCompl = (SiteTools::isComplete(* currentSite) ? 1 : 0);
    int isConst = (SiteTools::isConstant(* currentSite, true) ? 1 : 0);
    out << "[" << currentSitePosition << "]\t";
    out << isCompl << "\t";
    out << isConst << "\t";
    out << classes[i] << "\t";
    out << rates[i] << "\t";
    out << logLn[i] << endl;
  }
}

/******************************************************************************/

Statistic* CoETools::getStatistic(map<string, string>& params, const Alphabet* alphabet) throw (Exception)
{
  string statistic = ApplicationTools::getStringParameter("statistic", params, "Correlation");
  string name;
  map<string,string> args;
  KeyvalTools::parseProcedure(statistic, name, args);
  if (name == "Cosinus")
  {
    return new CosinusStatistic();
  }
  else if (name == "Correlation")
  {
    return new CorrelationStatistic();
  }
  else if (name == "CorrectedCorrelation")
  {
    return new CorrectedCorrelationStatistic();
  }
  else if (name == "Covariance")
  {
    return new CovarianceStatistic();
  }
  else if (name == "Cosubstitution")
  {
    return new CosubstitutionNumberStatistic();
  }
  else if (name == "Compensation")
  {
    string nijtOption = ApplicationTools::getStringParameter("nijt", params, "simule", "", true);
    bool sym = ApplicationTools::getBooleanParameter("nijt_aadist.sym", params, true, "", true); 
    if (nijtOption != "aadist" || sym)
    {
      throw Exception("Compensation distance must be used with 'nijt=aadist' and 'nijt_aadist.sym=no' options.");
    }
    else
    {
		  return new CompensationStatistic();
    }
  }
  else if (name == "MI")
  {
    string nijtOption = ApplicationTools::getStringParameter("nijt", params, "simple", "", true, false);
    if (nijtOption == "label") { //Warning, does not work with a second data set with a different alphabet :s
      bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true, "", true, 4); //For benchmarking only!
      if (average) {
        throw Exception("MI distance with 'nijt=label' can't be used with 'nijt.average=yes'.");
      }
      unsigned int n = alphabet->getSize() * (alphabet->getSize() - 1);
      vector<double> b(n + 2);
      b[0] = -0.5;
      for (unsigned int i = 0; i < n + 1; ++i)
        b[i + 1] = b [i] + 1;
      return new DiscreteMutualInformationStatistic(b);  
    } else if (nijtOption == "simple" || nijtOption == "laplace" || nijtOption == "prob_one_jump") {
      double threshold = ApplicationTools::getDoubleParameter("threshold", args, 0.99, "", true);
      vector<double> b(3);
      b[0] = 0.; b[1] = threshold; b[2] = 10000.;
		  return new DiscreteMutualInformationStatistic(b);
    } else {
      throw Exception("MI distance can only be used with 'nijt=label', 'nijt=simple', 'nijt=laplace', 'nijt=prob_one_jump' options for now.");
    }
  }
  else
  {
    throw Exception("Unknown statistic used: " + statistic);
  }
}

/******************************************************************************/

void CoETools::computeIntraStats(
  const DRTreeLikelihood& tl,
  const SequenceSimulator& seqSim,
  const SiteContainer& completeSites,
  ProbabilisticSubstitutionMapping& mapping,
	SubstitutionCount& nijt,
  const Statistic& statistic,
  bool computeNull,
  map<string, string>& params)
{
  //Vdouble branchLengths = vectors[0];
  //vectors.erase(vectors.begin());//remove branch lengths.
  
  string statFilePath = ApplicationTools::getAFilePath("statistic.output.file", params, true, false, "", false, "statistics.txt", 1);
  ofstream statOut(statFilePath.c_str(), ios::out);

  //Getting parameters:
    
  int minRateClass     = getMinRateClass(params);
  int maxRateClassDiff = getMaxRateClassDiff(params);
  double minRate       = getMinRate(params);
  double maxRateDiff   = getMaxRateDiff(params);
  double minStatistic  = getStatisticMin(params);

  size_t nbSites = mapping.getNumberOfSites();
  vector<double> norms = AnalysisTools::computeNorms(mapping);

  //Compute null distribution?
  //tl would be modified after the simulations, so we copy it before.
  vector< vector<double> >* simValues = 0; 
  size_t nbRateClasses = 1;
  auto_ptr<Domain> rateDomain;
  if (computeNull)
  { 
    nbRateClasses = ApplicationTools::getParameter<unsigned int>("statistic.null.nb_rate_classes", params, 10, "", false, 1);
    ApplicationTools::displayResult("Number of sub-distributions", nbRateClasses);
    rateDomain.reset(new Domain(0, VectorTools::max(norms), nbRateClasses));
    auto_ptr<DRTreeLikelihood> tlCopy(tl.clone());
    simValues = computeIntraNullDistribution(
        *tlCopy,
        rateDomain.get(),
        seqSim,
        nijt,
        statistic,
        params);
    //We need to sort observations for a better efficiency:
    for (size_t i = 0; i < simValues->size(); ++i) {
      sort((*simValues)[i].begin(), (*simValues)[i].end());
    }
  }

  ApplicationTools::displayMessage("\n\n-*- Compute statistics -*- \n");
  ApplicationTools::displayMessage(
    TextTools::toString(nbSites) + 
    " sites => " +
    TextTools::toString(nbSites * (nbSites + 1) / 2) +
    " pairs to compute!");

  statOut << "Group\tStat\tRCmin\tPRmin\tNmin";
  if (computeNull)
    statOut << "\tPValue\tNsim";
  statOut << endl;

  ApplicationTools::displayTask("Analyse each site pair", true);
  
  vector<size_t> classes = tl.getRateClassWithMaxPostProbOfEachSite();
  Vdouble rates = tl.getPosteriorRateOfEachSite();

  for (size_t i = 0; i < nbSites; i++)
  {
    int    iClass = static_cast<int>(classes[i]);
    double iRate  = rates[i];
    if (iClass < minRateClass) continue;
    if (iRate  < minRate     ) continue;
    double iNorm  = norms[i];
    ApplicationTools::displayGauge(i, nbSites - 1);
    for (size_t j = i + 1; j < nbSites; j++)
    {
      int    jClass = static_cast<int>(classes[j]);
      double jRate  = rates[j];
      if (jClass < minRateClass) continue;
      if (jRate  < minRate     ) continue;
      double jNorm  = norms[j];
    
      //Sites which are in too different rate classes are not compared:
      if (maxRateClassDiff >= 0  && NumTools::abs(jClass - iClass) > maxRateClassDiff) continue;
      if (maxRateDiff      >= 0. && NumTools::abs(jRate  - iRate ) > maxRateDiff)      continue;

      double stat = statistic.getValueForPair(mapping[i], mapping[j]);
      if (NumTools::abs(stat) < minStatistic) continue;

      double minNorm = min(iNorm, jNorm);

      //Then print to file:
      statOut << "[";
      statOut << completeSites.getSite(i).getPosition();
      statOut << ";";
      statOut << completeSites.getSite(j).getPosition();
      statOut << "]\t";
      statOut << stat;
      statOut << "\t";
      statOut << min(iClass, jClass);
      statOut << "\t";
      statOut << min(iRate, jRate);
      statOut << "\t";
      statOut << minNorm;
      if (computeNull) {
        try {
          size_t cat = rateDomain->getIndex(minNorm);
          size_t nsim = (*simValues)[cat].size();
          unsigned int count;
          for (count = 0; count < nsim && (*simValues)[cat][count] < stat; ++count) {}
          double pvalue = static_cast<double>(nsim - count + 1) / static_cast<double>(nsim + 1);
          statOut << "\t" << pvalue << "\t" << nsim;
        } catch (OutOfRangeException& oore) {
          statOut << "\tNA\t0";
        }
      }
      statOut << endl;
    }
  }

  ApplicationTools::displayTaskDone();
  statOut.close();
}

/******************************************************************************/

void CoETools::computeInterStats(
  const DiscreteRatesAcrossSitesTreeLikelihood& tl1,
  const DiscreteRatesAcrossSitesTreeLikelihood& tl2,
  const SiteContainer& completeSites1,
  const SiteContainer& completeSites2,
  ProbabilisticSubstitutionMapping& mapping1,
  ProbabilisticSubstitutionMapping& mapping2,
  const Statistic& statistic,
  map<string, string>& params)
{
  // Compute statistics from data:

  bool indepComp = haveToPerformIndependantComparisons(params);
  if (indepComp && mapping1.getNumberOfSites() != mapping2.getNumberOfSites())
  {
    ApplicationTools::displayError("When performing independant comparisons, the two datasets must have the same length.");
    exit(-1);
  }

  string statFilePath = ApplicationTools::getAFilePath("statistic.output.file", params, true, false, "", false, "statistics.txt", 1);
  ofstream statOut(statFilePath.c_str(), ios::out);

  //Getting parameters:
  int minRateClass1    = getMinRateClass(params);
  int minRateClass2    = getMinRateClass(params, "2");
  int maxRateClassDiff = getMaxRateClassDiff(params);
  double minRate1      = getMinRate(params);
  double minRate2      = getMinRate(params, "2");
  double maxRateDiff   = getMaxRateDiff(params);
  double minStatistic  = getStatisticMin(params);

  size_t nbSites1 = mapping1.getNumberOfSites();
  size_t nbSites2 = mapping2.getNumberOfSites();
  vector<double> norms1 = AnalysisTools::computeNorms(mapping1);
  vector<double> norms2 = AnalysisTools::computeNorms(mapping2);


  ApplicationTools::displayMessage(
    TextTools::toString(nbSites1) +
    " sites * " +
    TextTools::toString(nbSites2) +
    " = " +
    TextTools::toString(indepComp ? nbSites1 : nbSites1 * nbSites2) +
    " pairs to compute!");

  statOut << "Group\tStat\tRCmin\tPRmin\tNmin" << endl;

  ApplicationTools::displayTask("Analyse each site pair", true);
    
  vector<size_t> classes1 = tl1.getRateClassWithMaxPostProbOfEachSite();
  vector<size_t> classes2 = tl2.getRateClassWithMaxPostProbOfEachSite();
  Vdouble rates1   = tl1.getPosteriorRateOfEachSite();
  Vdouble rates2   = tl2.getPosteriorRateOfEachSite();

  for (size_t i = 0; i < nbSites1; i++)
  {
    int    iClass = static_cast<int>(classes1[i]);
    double iRate  = rates1[i];
    if (iClass < minRateClass1) continue;
    if (iRate  < minRate1     ) continue;
    double iNorm  = norms1[i];
    ApplicationTools::displayGauge(i, nbSites1 - 1);
      
    size_t begin = indepComp ? i : 0;
    size_t end   = indepComp ? i + 1 : nbSites2;
    for (size_t j = begin; j < end; j++)
    {
      int    jClass = static_cast<int>(classes2[j]);
      double jRate  = rates2[j];
      if (jClass < minRateClass2) continue;
      if (jRate  < minRate2     ) continue;
      double jNorm  = norms2[i];
    
      //Sites which are in too different rate classes are not compared:
      if (maxRateClassDiff >= 0  && NumTools::abs(jClass - iClass) > maxRateClassDiff) continue;
      if (maxRateDiff      >= 0. && NumTools::abs(jRate  - iRate ) > maxRateDiff     ) continue;

      //double stat = table[i + 1][j + 1];
      double stat = statistic.getValueForPair(mapping1[i], mapping2[j]);
      if (NumTools::abs(stat) < minStatistic) continue;

      //Then print to file:
      statOut << "[";
      statOut << completeSites1.getSite(i).getPosition();
      statOut << ";";
      statOut << completeSites2.getSite(j).getPosition();
      statOut << "]\t";
      statOut << stat;
      statOut << "\t";
      statOut << min(iClass, jClass);
      statOut << "\t";
      statOut << min(iRate, jRate);
      statOut << "\t";
      statOut << min(iNorm, jNorm);
      statOut << endl;
    }
  }

  ApplicationTools::displayTaskDone();
  statOut.close();
}

/******************************************************************************/

vector< vector<double> >* CoETools::computeIntraNullDistribution(
    DRTreeLikelihood& drtl,
    const Domain* rateDomain,
    const SequenceSimulator& seqSim,
    SubstitutionCount& nijt,
    const Statistic& statistic,
    map<string, string>& params)
{
  string path = ApplicationTools::getAFilePath("statistic.null.output.file", params, false, false, "", false, "none", 1);
  auto_ptr<ofstream> outFile;
  if (path != "none") {
    outFile.reset(new ofstream(path.c_str(), ios::out));
    ApplicationTools::displayResult("Write simulation results to", path);
  }
 
  vector< vector<double> >* simValues = 0;
   
  unsigned int nbRepCPU = ApplicationTools::getParameter<unsigned int>("statistic.null.nb_rep_CPU", params, 100, "", false, 1);
  unsigned int nbRepRAM = ApplicationTools::getParameter<unsigned int>("statistic.null.nb_rep_RAM", params, 1000, "", false, 1);
  bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true, "", true, 3); //This is for testing purpose only
  bool joint   = ApplicationTools::getBooleanParameter("nijt.joint", params, true, "", true, 3);

  bool computePValues = ApplicationTools::getBooleanParameter("statistic.null.compute_pvalue", params, true, "", false, 2);
 
  if (computePValues)
  {
    simValues = new vector< vector<double> >(rateDomain->getSize());
  }

  ApplicationTools::displayResult("Nb of simulations to perform", nbRepCPU * nbRepRAM);
  AnalysisTools::getNullDistributionIntraDR(drtl, seqSim, nijt, statistic, outFile.get(), simValues, rateDomain, nbRepCPU, nbRepRAM, average, joint, true);
  ApplicationTools::message->endLine();

  if (outFile.get())
    outFile->close();
  return simValues;
}

/******************************************************************************/

void CoETools::computeInterNullDistribution(
    DRTreeLikelihood& drtl1,
    DRTreeLikelihood& drtl2,
    const SequenceSimulator& seqSim1,
    const SequenceSimulator& seqSim2,
    SubstitutionCount& nijt1,
    SubstitutionCount& nijt2,
    const Statistic& statistic,
    map<string, string>& params)
{
  string path = ApplicationTools::getAFilePath("statistic.null.output.file", params, true, false, "", false, "statistics.null.txt", 1);
  ofstream outFile(path.c_str(), ios::out);
  
  int nbRepCPU = ApplicationTools::getIntParameter("statistic.null.nb_rep_CPU", params, 10, "", false, 1);
  int nbRepRAM = ApplicationTools::getIntParameter("statistic.null.nb_rep_RAM", params, 1000, "", false, 1);
  bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true, "", true, 3); //For testing purpose only
  bool joint   = ApplicationTools::getBooleanParameter("nijt.joint", params, true, "", true, 3);

  ApplicationTools::displayMessage("Compute statistic under null hypothesis...");
  AnalysisTools::getNullDistributionInterDR(drtl1, drtl2, seqSim1, seqSim2, nijt1, nijt2, statistic, outFile, nbRepCPU, nbRepRAM, average, joint, true);
  outFile.close();
}

/******************************************************************************/

vector<unsigned int> CandidateGroupSet::nextCandidateSite() const throw (Exception)
{
  if (nbCompleted_ == size()) throw Exception("CandidateGroupSet::nextCandidateSite. enough simulations!!");
  //Site increment:
  if (n2_[groupPos_] < minSim_)
  {
    sitePos_++;
    if (sitePos_ >= (*this)[groupPos_].size())
    {
      groupPos_++;
      if (groupPos_ >= size()) groupPos_ = 0;
      sitePos_ = 0;
    }
  }
  unsigned int startSearch = groupPos_;
  if (n2_[groupPos_] >= minSim_ || !candidates_[groupPos_].isAnalysable())
  {
    while (n2_[groupPos_] >= minSim_ || !candidates_[groupPos_].isAnalysable())
    {
      groupPos_++;

      if (groupPos_ >= size()) groupPos_ = 0;
      if (groupPos_ == startSearch)
      {
        //No more site to complete!
        throw Exception("DEBUG: something wrong happened, this message should never appear!");
      }
    }
    sitePos_ = 0;
  }
  vector<unsigned int> pos(2);
  pos[0] = groupPos_;
  pos[1] = sitePos_;
  return pos;
}

/******************************************************************************/

vector<unsigned int> CandidateGroupSet::currentCandidateSite() const throw (Exception)
{
  if (nbCompleted_ == size()) throw Exception("CandidateGroupSet::nextCandidateSite. enough simulations!!");
  vector<unsigned int> pos(2);
  pos[0] = groupPos_;
  pos[1] = sitePos_;
  return pos;
}

/******************************************************************************/

bool CandidateGroupSet::analyseSimulations(const ProbabilisticSubstitutionMapping& mapping)
{
  Vdouble norms = AnalysisTools::computeNorms(mapping);
  //Analyse each site in the set:
  bool test = true, testNorm, first, testFree = true;
  vector<unsigned int> pos, start;
  for (unsigned int i = 0; test && i < mapping.getNumberOfSites(); i++)
  {
    first = true;
    testNorm = false;
    while (test && !testNorm)
    {
      pos = nextCandidateSite();
      if (first)
      {
        start = pos;
        first = false;
      }
      else
      {
        if (currentCandidateSite() == start) //We looped over all set, drop this simulated site:
          break;
      }
      testNorm = candidates_[pos[0]][pos[1]].checkNorm(norms[i]);
      if (testNorm)
      {
        //cout << pos[0] << "\t" << pos[1] << "\t" << simulations_[pos[0]][pos[1]].size() << endl;
        if (addSimulatedSite(pos[0], pos[1], &mapping[i]))
        {
          testFree = false; //At least one groupe was completed!
        }
        if (nbCompleted_ == nbAnalysable_) test = false;
      }
    }
  }
  if (testFree)
  {
    if (verbose_)
      ApplicationTools::displayMessage("This simulation set did not provide sites with appropriate norm :(");
    nbTrials_++;
  }

  //Reset pointers: (we need to do that since we do not recopy the vectors!)
  resetSimulations();
  return test;
}

/******************************************************************************/

bool CandidateGroupSet::addSimulatedSite(unsigned int groupIndex, unsigned int siteIndex, const VVdouble* v) throw (IndexOutOfBoundsException)
{
  if (groupIndex >= simulations_.size())
    throw IndexOutOfBoundsException("CandidateGroupSet::addSimulatedSite. Bad group index.", groupIndex, 0, simulations_.size());
  if (siteIndex >= simulations_[groupIndex].size())
    throw IndexOutOfBoundsException("CandidateGroupSet::addSimulatedSite. Bad site index.", siteIndex, 0, simulations_[groupIndex].size());
  vector< deque<const VVdouble*> >* group = &simulations_[groupIndex];
  (*group)[siteIndex].push_back(v);
  //Test if the group is complete:
  bool test = true;
  for (unsigned int i = 0; test && i < group->size(); ++i)
  {
    if ((*group)[i].size() == 0) test = false;
  }
  if (test)
  {
    vector<const VVdouble*> groupVectors;
    for (unsigned int i = 0; i < group->size(); i++)
    {
      groupVectors.push_back((*group)[i][0]);
      (*group)[i].pop_front();
    }
    n2_[groupIndex]++;
    double stat = statistic_->getValueForGroup(groupVectors);
    if (stat >= (*this)[groupIndex].getStatisticValue()) n1_[groupIndex]++;
    if (n2_[groupIndex] == minSim_)
    {
      nbCompleted_++;
      if (verbose_ == 1)
      {
        ApplicationTools::displayGauge(nbCompleted_, size(), '=');
      }
      if (verbose_ > 1)
      {
        ApplicationTools::displayResult("Group completed", TextTools::toString(groupIndex) + " (" + TextTools::toString(nbCompleted_) + "/" + TextTools::toString(nbAnalysable_) + ")");
      }
    }
  }
  return test;
}

/******************************************************************************/

void CoETools::computePValuesForCandidateGroups(
    CandidateGroupSet& candidates,
    DRTreeLikelihood& drtl,
    const SequenceSimulator& seqSim,
    SubstitutionCount& nijt,
	  map<string, string>& params,
    unsigned int maxTrials)
{
  unsigned int repRAM = ApplicationTools::getParameter<unsigned int>("candidates.null.nb_rep_RAM", params, 1000, "", true, 1);
  bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true, "", true, 3); //For testing purpose only
  bool joint   = ApplicationTools::getBooleanParameter("nijt.joint", params, true, "", true, 3);

  bool test = true;
  while(test)
  {
    if (candidates.getVerbose() >= 2)
      ApplicationTools::displayResult("Simulate ", TextTools::toString(repRAM) + " sites.");
    SiteContainer* sites = seqSim.simulate(repRAM);
    drtl.setData(*sites);
    drtl.initialize();
    ProbabilisticSubstitutionMapping* mapping;
		if (average)
    {
			if (joint)
      {
				mapping = SubstitutionMappingTools::computeSubstitutionVectors(drtl, nijt, false);
			}
      else
      {
				mapping = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl, nijt, false);
			}
		}
    else
    {
			if (joint)
      {
				mapping = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl, nijt, false);
			}
      else
      {
				mapping = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl, nijt, false);
			}
		}
		test = candidates.analyseSimulations(*mapping) && (candidates.getNumberOfTrials() < maxTrials);

    delete sites;
    delete mapping;
  }
}

/******************************************************************************/

