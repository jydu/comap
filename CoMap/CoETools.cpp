//
// File: CoeTools.cpp
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

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/KeyvalTools.h>
#include <Utils/StringTokenizer.h>
#include <Utils/ApplicationTools.h>
#include <Utils/BppString.h>

// From NumCalc:
#include <NumCalc/Domain.h>
#include <NumCalc/IntervalData.h>
#include <NumCalc/ParameterExceptions.h>

// From SeqLib:
#include <Seq/SiteTools.h>
#include <Seq/SequenceApplicationTools.h>
#include <Seq/GranthamAAChemicalDistance.h>
#include <Seq/MiyataAAChemicalDistance.h>
#include <Seq/KleinAANetChargeIndex.h>
#include <Seq/AAChargeIndex.h>
#include <Seq/GranthamAAPolarityIndex.h>
#include <Seq/GranthamAAVolumeIndex.h>
#include <Seq/SimpleIndexDistance.h>
#include <Seq/AAIndex1Entry.h>
#include <Seq/AAIndex2Entry.h>

// From PhylLib:
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/PatternTools.h>
#include <Phyl/OptimizationTools.h>
#include <Phyl/SubstitutionCount.h>
#include <Phyl/AnalyticalSubstitutionCount.h>
#include <Phyl/SimpleSubstitutionCount.h>
#include <Phyl/IndexToCount.h>
#include <Phyl/SubstitutionMappingTools.h>
#include <Phyl/Newick.h>
#include <Phyl/SubstitutionModelSet.h>
#include <Phyl/SubstitutionModelSetTools.h>
#include <Phyl/DRHomogeneousTreeLikelihood.h>
#include <Phyl/DRNonHomogeneousTreeLikelihood.h>

using namespace bpp;

// From the STL:
#include <fstream>
#include <iomanip>
using namespace std;

/******************************************************************************/

void CoETools::readData(
  TreeTemplate<Node> *          tree,
  Alphabet *                    &alphabet,
  VectorSiteContainer *         &allSites,
  VectorSiteContainer *         &sites,
  SubstitutionModel *           &model,
  SubstitutionModelSet *        &modelSet,
  DiscreteDistribution *        &rDist,
  DRTreeLikelihood *            &tl,
  map<string, string>           &params,
  const string                  &suffix)
{
  alphabet = SequenceApplicationTools::getAlphabet(params, suffix, true);
  allSites = SequenceApplicationTools::getSiteContainer(alphabet, params, suffix, false);
  sites    = SequenceApplicationTools::getSitesToAnalyse(*allSites, params, suffix, true, true, true);

  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", params, "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  model = 0;
  modelSet = 0;

  if (nhOpt == "no")
  {  
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, params);
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
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, params);
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
      unsigned int n =(unsigned int)(model->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                    // we should assume a rate distribution for the root also!!!  
    }
    FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, sites, params, rateFreqs);
    vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", params, ',', "");
    modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
    model = modelSet->getModel(0)->clone();
    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false);
  }
  else if (nhOpt == "general")
  {
    modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet,0, params);
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
    for (unsigned int i = 0; i < allSites->getNumberOfSites(); i++)
    {
      debug << "Position " << i+1 << " = " << tl->getLogLikelihoodForASite(i) << endl; 
    }
    debug.close();
    exit(-1);
  }
  (ApplicationTools::message->setPrecision(20) << ll).endLine();
  
  bool optimize = ApplicationTools::getBooleanParameter("optimization", params, true, suffix, true, false);
  if (optimize)
  {
    ApplicationTools::displayResult("Optimization", (optimize ? "yes" : "no"));
    PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), params, suffix, true, true);
    TreeTemplate<Node> *treeOpt = new TreeTemplate<Node>(tl->getTree());
    PhylogeneticsApplicationTools::writeTree(*treeOpt, params, "output.", suffix);
  
    // Actualize tree.
    // Substitution model and rate distribution are actualized automatically by
    // the likelihood function.
    tree  = treeOpt;
    
    // Print parameters:
    ApplicationTools::displayResult("Final likelihood, -lnL =", TextTools::toString(tl -> getLogLikelihood(), 20));
  }
  string tags = ApplicationTools::getAFilePath("output.tags.file", params, false, false, suffix, false);
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

  bool removeConst = ApplicationTools::getBooleanParameter("input.remove_const", params, true);
  if (removeConst) {
    ApplicationTools::displayTask("Remove conserved positions", true);
    unsigned int n = sites->getNumberOfSites();
    for (unsigned int i = n; i > 0; --i) {
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
  ProbabilisticSubstitutionMapping * substitutions = 0;
  string inputVectorsFilePath = ApplicationTools::getAFilePath("input.vectors.file", params, false, false, suffix, false);

  if (inputVectorsFilePath != "none")
  {
    //We try to load the substitution vector directly from file:
    int nbSites = drtl.getNumberOfSites();
    ApplicationTools::displayResult("Substitution mapping in file:", inputVectorsFilePath);
    ifstream sc(inputVectorsFilePath.c_str(), ios::in);
    substitutions = new ProbabilisticSubstitutionMapping(drtl.getTree(), nbSites);
    SubstitutionMappingTools::readFromStream(sc, *substitutions);
  }
  else
  {
    //We compute the substitutions vector:

    string outputVectorsFilePath = ApplicationTools::getAFilePath("output.vectors.file", params, true, false, suffix, false);
    ApplicationTools::displayResult("Output mapping to file" + suffix, outputVectorsFilePath);

    bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true);
    bool joint   = ApplicationTools::getBooleanParameter("nijt.joint"  , params, true);
    if (average)
    {
      if (joint)
      {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectors(drtl, substitutionCount);
      }
      else
      {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl, substitutionCount);
      }
    }
    else
    {
      if (joint)
      {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl, substitutionCount);
      }
      else
      {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl, substitutionCount);
      }
    }
    if (outputVectorsFilePath != "none")
    {
      ofstream outputVectors(outputVectorsFilePath.c_str(), ios::out);
      SubstitutionMappingTools::writeToStream(*substitutions, completeSites, outputVectors);
      outputVectors.close();
    }

  }
  return substitutions;
}

/******************************************************************************/

int CoETools::getMinRateClass(map<string, string>& params, string suffix)
{
  int minRateClass = ApplicationTools::getIntParameter("statistic.min_rate_class", params, 0, suffix, true);
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
  double minRate = ApplicationTools::getDoubleParameter("statistic.min_rate", params, 0., suffix, true);
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
  int maxRateClassDiff = ApplicationTools::getIntParameter("statistic.max_rate_class_diff", params, -1);
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
  double maxRateDiff = ApplicationTools::getDoubleParameter("statistic.max_rate_diff", params, -1.);
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
  double minStatistic = ApplicationTools::getDoubleParameter("statistic.min", params, 0);
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
  bool indepComp = ApplicationTools::getBooleanParameter("independant_comparisons", params, false);
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
  string outFile = ApplicationTools::getAFilePath("output.infos", params, false, false, suffix, true);
  if (outFile == "none") return;

  // Get the rate class with maximum posterior probability:
  vector<unsigned int> classes = ras.getRateClassWithMaxPostProbOfEachSite();
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

const Statistic* CoETools::getStatistic(map<string, string>& params) throw (Exception)
{
  string statistic = ApplicationTools::getStringParameter("statistic", params, "none");
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
    string nijtOption = ApplicationTools::getStringParameter("nijt", params, "simule", "", true);
    if (nijtOption == "simple" || nijtOption == "laplace" || nijtOption == "prob_one_jump")
    {
      double threshold = ApplicationTools::getDoubleParameter("threshold", args, 0.99, "", true);
      vector<double> b(3);
      b[0] = 0.; b[1] = threshold; b[2] = 10000.;
		  return new DiscreteMutualInformationStatistic(b);
    }
    else
    {
      throw Exception("MI distance can only be used with 'nijt=simple', 'nijt=laplace', 'nijt=prob_one_jump' options for now.");
    }
  }
  else
  {
    throw Exception("Unknown statistic used: " + statistic);
  }
}

/******************************************************************************/

SubstitutionCount* CoETools::getSubstitutionCount(
  const Alphabet* alphabet,
  const SubstitutionModel* model,
  const DiscreteDistribution* rDist,
  map<string, string>& params,
  string suffix)
{
  SubstitutionCount* substitutionCount = 0;
  string nijtOption = ApplicationTools::getStringParameter("nijt", params, "simule", suffix, true);

  if (nijtOption == "laplace")
  {
    int trunc = ApplicationTools::getIntParameter("nijt_laplace.trunc", params, 10, suffix, true);
    substitutionCount = new AnalyticalSubstitutionCount(model, trunc);
  }
  else if (nijtOption == "simple")
  {
    substitutionCount = new SimpleSubstitutionCount(alphabet);
  }
  else if (nijtOption == "aadist")
  {
    if (!AlphabetTools::isProteicAlphabet(alphabet))
    {
      ApplicationTools::displayError("Chemical distance can only be used with protein data.");
      exit(-1);
    }
    string dist = ApplicationTools::getStringParameter("nijt_aadist.type", params, "grantham", suffix, true);
    bool sym = ApplicationTools::getBooleanParameter("nijt_aadist.sym", params, true, suffix, true);
    if (dist == "grantham")
    {
      GranthamAAChemicalDistance* M = new GranthamAAChemicalDistance();
      M->setSymmetric(sym);
      if (!sym) M->setPC1Sign(true);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if (dist == "miyata")
    {
      MiyataAAChemicalDistance* M = new MiyataAAChemicalDistance();
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if (dist == "grantham.polarity")
    {
      GranthamAAPolarityIndex I;
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if (dist == "grantham.volume")
    {
      GranthamAAVolumeIndex I;
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if (dist == "klein.charge")
    {
      KleinAANetChargeIndex I;
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if (dist == "charge")
    {
      AAChargeIndex I;
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if (dist == "user1")
    {
      string aax1FilePath = ApplicationTools::getAFilePath("nijt_aadist.type_user1.file", params, true, true, suffix, false);
      ifstream aax1File(aax1FilePath.c_str(), ios::in);
      AAIndex1Entry I(aax1File);
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      aax1File.close();
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if (dist == "user2")
    {
      string aax2FilePath = ApplicationTools::getAFilePath("nijt_aadist.type_user2.file", params, true, true, suffix, false);
      ifstream aax2File(aax2FilePath.c_str(), ios::in);
      AAIndex2Entry * M = new AAIndex2Entry(aax2File, sym);
      aax2File.close();
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else
    {
      ApplicationTools::displayError("Invalid distance '" + dist + ", in 'nijt_aadist' parameter.");
      exit(-1);
    }
  }
  else if (nijtOption == "prob_one_jump")
  {
    substitutionCount = new OneJumpSubstitutionCount(model);
  }
  else
  {
    ApplicationTools::displayError("Invalid option '" + nijtOption + ", in 'nijt' parameter.");
    exit(-1);
  }

  // Send results:
  return substitutionCount;
}

/******************************************************************************/

void CoETools::computeIntraStats(
  const DiscreteRatesAcrossSitesTreeLikelihood& tl,
  const SiteContainer& completeSites,
  ProbabilisticSubstitutionMapping& mapping,
  const Statistic& statistic,
  map<string, string>& params)
{
  //Vdouble branchLengths = vectors[0];
  //vectors.erase(vectors.begin());//remove branch lengths.
  
  string statFilePath = ApplicationTools::getAFilePath("statistic.output.file", params, true, false);
  ofstream statOut(statFilePath.c_str(), ios::out);

  //Getting parameters:
    
  int minRateClass     = getMinRateClass(params);
  int maxRateClassDiff = getMaxRateClassDiff(params);
  double minRate       = getMinRate(params);
  double maxRateDiff   = getMaxRateDiff(params);
  double minStatistic  = getStatisticMin(params);

  unsigned int nbSites = mapping.getNumberOfSites();
  vector<double> norms = AnalysisTools::computeNorms(mapping);

  ApplicationTools::displayMessage(
    TextTools::toString(nbSites) + 
    " sites => " +
    TextTools::toString(nbSites * (nbSites + 1) / 2) +
    " pairs to compute!");

  statOut << "Group\tStat\tRCmin\tPRmin\tNmin" << endl;

  ApplicationTools::displayTask("Analyse each site pair", true);
  
  vector<unsigned int> classes = tl.getRateClassWithMaxPostProbOfEachSite();
  Vdouble rates   = tl.getPosteriorRateOfEachSite();

  for (unsigned int i = 0; i < nbSites; i++)
  {
    int    iClass = classes[i];
    double iRate  = rates[i];
    if (iClass < minRateClass) continue;
    if (iRate  < minRate     ) continue;
    double iNorm  = norms[i];
    ApplicationTools::displayGauge(i, nbSites - 1);
    for (unsigned int j = i + 1; j < nbSites; j++)
    {
      int    jClass = classes[j];
      double jRate  = rates[j];
      if (jClass < minRateClass) continue;
      if (jRate  < minRate     ) continue;
      double jNorm  = norms[j];
    
      //Sites which are in too different rate classes are not compared:
      if (maxRateClassDiff >= 0  && NumTools::abs(jClass - iClass) > maxRateClassDiff) continue;
      if (maxRateDiff      >= 0. && NumTools::abs(jRate  - iRate ) > maxRateDiff)      continue;

      double stat = statistic.getValueForPair(mapping[i], mapping[j]);
      if (NumTools::abs(stat) < minStatistic) continue;

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
      statOut << min(iNorm, jNorm);
      statOut << endl;
    }
  }

  ApplicationTools::displayTaskDone();
  statOut.close();
}

/******************************************************************************/

void CoETools::computeInterStats(
  const DiscreteRatesAcrossSitesTreeLikelihood & tl1,
  const DiscreteRatesAcrossSitesTreeLikelihood & tl2,
  const SiteContainer & completeSites1,
  const SiteContainer & completeSites2,
  ProbabilisticSubstitutionMapping & mapping1,
  ProbabilisticSubstitutionMapping & mapping2,
  const Statistic & statistic,
  map<string, string> & params)
{
  // Compute statistics from data:

  bool indepComp = haveToPerformIndependantComparisons(params);
  if (indepComp && mapping1.getNumberOfSites() != mapping2.getNumberOfSites())
  {
    ApplicationTools::displayError("When performing independant comparisons, the two datasets must have the same length.");
    exit(-1);
  }

  string statFilePath = ApplicationTools::getAFilePath("statistic.output.file", params, true, false);
  ofstream statOut(statFilePath.c_str(), ios::out);

  //Getting parameters:
  int minRateClass1    = getMinRateClass(params);
  int minRateClass2    = getMinRateClass(params, "2");
  int maxRateClassDiff = getMaxRateClassDiff(params);
  double minRate1      = getMinRate(params);
  double minRate2      = getMinRate(params, "2");
  double maxRateDiff   = getMaxRateDiff(params);
  double minStatistic  = getStatisticMin(params);

  unsigned int nbSites1 = mapping1.getNumberOfSites();
  unsigned int nbSites2 = mapping2.getNumberOfSites();
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
    
  vector<unsigned int> classes1 = tl1.getRateClassWithMaxPostProbOfEachSite();
  vector<unsigned int> classes2 = tl2.getRateClassWithMaxPostProbOfEachSite();
  Vdouble rates1   = tl1.getPosteriorRateOfEachSite();
  Vdouble rates2   = tl2.getPosteriorRateOfEachSite();

  for (unsigned int i = 0; i < nbSites1; i++)
  {
    int    iClass = classes1[i];
    double iRate  = rates1[i];
    if (iClass < minRateClass1) continue;
    if (iRate  < minRate1     ) continue;
    double iNorm  = norms1[i];
    ApplicationTools::displayGauge(i, nbSites1 - 1);
      
    unsigned int begin = indepComp ? i : 0;
    unsigned int end   = indepComp ? i + 1 : nbSites2;
    for (unsigned int j = begin; j < end; j++)
    {
      int    jClass = classes2[j];
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

void CoETools::computeIntraNullDistribution(
    DRTreeLikelihood & drtl,
    const SequenceSimulator& seqSim,
    SubstitutionCount & nijt,
    const Statistic & statistic,
    map<string, string> & params)
{
  const DiscreteDistribution* rDist = drtl.getRateDistribution();

  string path = ApplicationTools::getAFilePath("statistic.null.output.file", params, true, false);
  ofstream outFile(path.c_str(), ios::out);
  
  ApplicationTools::displayMessage("Compute statistic under null hypothesis...");
  
  unsigned int nbRepCPU = ApplicationTools::getParameter<unsigned int>("statistic.null.nb_rep_CPU", params, 10);
  unsigned int nbRepRAM = ApplicationTools::getParameter<unsigned int>("statistic.null.nb_rep_RAM", params, 100);
  bool cumul   = ApplicationTools::getBooleanParameter("statistic.null.cumul", params, false);
  bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true);
  bool joint   = ApplicationTools::getBooleanParameter("nijt.joint", params, true);

  if (cumul)
  {
    // Building domain:
    double lowerSB = ApplicationTools::getDoubleParameter("statistic.null.lower", params, -1.);
    double upperSB = ApplicationTools::getDoubleParameter("statistic.null.upper", params, 1.);
    int     nbSInt = ApplicationTools::getIntParameter   ("statistic.null.nb_int", params, 20);
    Domain statDomain(lowerSB, upperSB, nbSInt);
    Domain rateDomain = rDist->getDomain();
  
    // Simulate:
    vector<IntervalData *> id; 
    id = AnalysisTools::getNullDistributionIntraDR(drtl, seqSim, nijt, statistic, statDomain, rateDomain, nbRepCPU, nbRepRAM, average, joint, true);
  
    // Print to file:
    for (unsigned int i = 0; i < rateDomain.getSize(); i++)
    {
      outFile << "# Distribution with minimum rate " << rateDomain.getValue(i) << endl;
      id[i]->print(outFile);
      outFile << endl;
    }

    // Free memory:
    for (unsigned int i = 0; i < id.size(); i++) delete id[i];
  }
  else
  {
    AnalysisTools::getNullDistributionIntraDR(drtl, seqSim, nijt, statistic, outFile, nbRepCPU, nbRepRAM, average, joint, true);
  }
  outFile.close();
}

/******************************************************************************/

void CoETools::computeInterNullDistribution(
    DRTreeLikelihood & drtl1,
    DRTreeLikelihood & drtl2,
    const SequenceSimulator& seqSim1,
    const SequenceSimulator& seqSim2,
    SubstitutionCount & nijt1,
    SubstitutionCount & nijt2,
    const Statistic & statistic,
    map<string, string> & params)
{
  string path = ApplicationTools::getAFilePath("statistic.null.output.file", params, true, false);
  ofstream outFile(path.c_str(), ios::out);
  
  int nbRepCPU = ApplicationTools::getIntParameter("statistic.null.nb_rep_CPU", params, 10);
  int nbRepRAM = ApplicationTools::getIntParameter("statistic.null.nb_rep_RAM", params, 100);
  bool cumul   = ApplicationTools::getBooleanParameter("statistic.null.cumul", params, false);
  bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true);
  bool joint   = ApplicationTools::getBooleanParameter("nijt.joint", params, true);

  ApplicationTools::displayMessage("Compute statistic under null hypothesis...");
  if (cumul)
  {
    // Building domain:
    double lowerSB = ApplicationTools::getDoubleParameter("statistic.null.lower", params, -1.);
    double upperSB = ApplicationTools::getDoubleParameter("statistic.null.upper", params, 1.);
    int     nbSInt = ApplicationTools::getIntParameter   ("statistic.null.nb_int", params, 20);
    Domain statDomain(lowerSB, upperSB, nbSInt);
    Vdouble bounds = ApplicationTools::getVectorParameter<double>("statistic.null.rate.bounds", params, ',', "");
    Domain rateDomain = Domain(bounds);
    ApplicationTools::displayMessage("The following rate domain will be used:");
    for (unsigned int i = 0; i < rateDomain.getSize(); i++)
    {
      *ApplicationTools::message
        << "["
        << rateDomain.getBound(i)
        << ", "
        << rateDomain.getBound(i+1)
        << "[, "
        << rateDomain.getValue(i);
      ApplicationTools::message->endLine();
    }
  
    // Simulate:
    vector<IntervalData *> id;
    id = AnalysisTools::getNullDistributionInterDR(
      drtl1, drtl2,
      seqSim1, seqSim2,
      nijt1, nijt2,
      statistic,
      statDomain, rateDomain,
      nbRepCPU, nbRepRAM,
      average, joint,
      true);
  
    // Print to file:
    for (unsigned int i = 0; i < rateDomain.getSize(); i++)
    {
      outFile << "# Distribution with minimum rate " << rateDomain.getValue(i) << endl;
      id[i]->print(outFile);
      outFile << endl;
    }
  
    // Free memory:
    for (unsigned int i = 0; i < id.size(); i++) delete id[i];
  }
  else
  {
    AnalysisTools::getNullDistributionInterDR(drtl1, drtl2, seqSim1, seqSim2, nijt1, nijt2, statistic, outFile, nbRepCPU, nbRepRAM, average, joint, true);
  }
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
  if (n2_[groupPos_] >= minSim_)
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

bool CandidateGroupSet::analyseSimulations(const ProbabilisticSubstitutionMapping & mapping)
{
  Vdouble norms = AnalysisTools::computeNorms(mapping);
  //Analyse each site in the set:
  bool test = true, testNorm, first, testFree = true;
  vector<unsigned int> pos, start;
  for (unsigned int i = 0; test && i < mapping.getNumberOfSites(); i++)
  {
    first = true;
    testNorm = false;
    while(test && !testNorm)
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
        if(addSimulatedSite(pos[0], pos[1], &mapping[i]))
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

bool CandidateGroupSet::addSimulatedSite(unsigned int groupIndex, unsigned int siteIndex, const Vdouble * v) throw (IndexOutOfBoundsException)
{
  if (groupIndex >= simulations_.size()) throw IndexOutOfBoundsException("CandidateGroupSet::addSimulatedSite. Bad group index.", groupIndex, 0, simulations_.size());
  if (siteIndex >= simulations_[groupIndex].size()) throw IndexOutOfBoundsException("CandidateGroupSet::addSimulatedSite. Bad site index.", siteIndex, 0, simulations_[groupIndex].size());
  vector< deque<const Vdouble*> >* group = &simulations_[groupIndex];
  (*group)[siteIndex].push_back(v);
  //Test if the group is complete:
  bool test = true;
  for (unsigned int i = 0; test && i < group->size(); i++)
  {
    if ((*group)[i].size() == 0) test = false;
  }
  if (test)
  {
    vector<const Vdouble *> groupVectors;
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
    CandidateGroupSet & candidates,
    DRTreeLikelihood & drtl,
    const SequenceSimulator& seqSim,
    SubstitutionCount & nijt,
	  map<string, string> & params,
    unsigned int maxTrials)
{
  unsigned int repRAM = ApplicationTools::getParameter<unsigned int>("candidates.null.nb_rep_RAM", params, 1000, "", true, true);
  bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true);
  bool joint   = ApplicationTools::getBooleanParameter("nijt.joint", params, true);

  bool test = true;
  while(test)
  {
    if (candidates.getVerbose() >= 2)
      ApplicationTools::displayResult("Simulate ", TextTools::toString(repRAM) + " sites.");
    SiteContainer * sites = seqSim.simulate(repRAM);
    drtl.setData(*sites);
    drtl.initialize();
    ProbabilisticSubstitutionMapping * mapping;
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

