//
// File: CoeTools.cpp
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

#include "CoETools.h"
#include "AnalysisTools.h"

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/StringTokenizer.h>
#include <Utils/ApplicationTools.h>
#include <Utils/String.h>

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
#include <Phyl/HomogeneousSequenceSimulator.h>
#include <Phyl/Newick.h>

// From the STL:
#include <fstream>
#include <iomanip>
using namespace std;

/******************************************************************************/

void CoETools::readData(
  TreeTemplate<Node> *          &tree,
  Alphabet *                    &alphabet,
  VectorSiteContainer *         &allSites,
  VectorSiteContainer *         &sites,
  SubstitutionModel *           &model,
  DiscreteDistribution *        &rDist,
  DRHomogeneousTreeLikelihood * &tl,
  map<string, string>           &params,
  const string                  &suffix)
{
  alphabet = SequenceApplicationTools::getAlphabet(params, suffix, true);
  allSites = SequenceApplicationTools::getSiteContainer(alphabet, params, suffix, false);
  sites    = SequenceApplicationTools::getSitesToAnalyse(*allSites, params, suffix, true, true, true);
  model    = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, params, suffix, true);
  rDist    = PhylogeneticsApplicationTools::getRateDistribution(params, suffix, true);
  
  tl = new DRHomogeneousTreeLikelihood(
    *tree,
    *sites,
    model,
    rDist, true, true);
  tl->initialize();

  ApplicationTools::displayTask("Tree likelihood");
  double ll = tl->getValue();
  if(isinf(ll))
  {
    ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
    ApplicationTools::displayError("!!! You should consider reestimating all branch lengths parameters.");
    exit(-1);
  }
  *ApplicationTools::message << setprecision(20) << ll << endl;
  
  bool optimize = ApplicationTools::getBooleanParameter("optimization", params, true, suffix, true, false);
  if(optimize)
  {
    ApplicationTools::displayResult("Optimization", (optimize ? "yes" : "no"));
    PhylogeneticsApplicationTools::optimizeParameters(tl, params, suffix, true, true);
    TreeTemplate<Node> *treeOpt = dynamic_cast<TreeTemplate<Node> *>(tl->getTree()->clone());
    PhylogeneticsApplicationTools::writeTree(*treeOpt, params, suffix);
  
    // Actualize tree.
    // Substitution model and rate distribution are actualized automatically by
    // the likelihood function.
    tree  = treeOpt;
    
    // Print parameters:
    ApplicationTools::displayResult("Final likelihood, -lnL =", TextTools::toString(tl -> getLogLikelihood(), 20));
    model->getParameters().printParameters(*ApplicationTools::message);
    rDist->getParameters().printParameters(*ApplicationTools::message);
  }
  string tags = ApplicationTools::getAFilePath("output.tags.file", params, false, false, suffix, false);
  if(tags != "none")
  {
    TreeTemplate<Node> treeCopy(*tree);
    vector<Node *> nodes = treeCopy.getInnerNodes();
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      nodes[i]->setNodeProperty("name", String(TextTools::toString(nodes[i]->getId())));
    }
    nodes = treeCopy.getLeaves();

    // Writing leaf names translation:
    string tlnPath = ApplicationTools::getAFilePath("output.tags.translation", params, false, false, suffix, false);
    if(tlnPath != "none")
    {
      ofstream tln(tlnPath.c_str(), ios::out);
      tln << "Name\tId" << endl;
      for(unsigned int i = 0; i < nodes.size(); i++)
      {
        tln << nodes[i]->getName() << "\t" << nodes[i]->getId() << endl;
      }
      tln.close();
    }

    // Translate names:
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      nodes[i]->setName(TextTools::toString(nodes[i]->getId()));
    }
    Newick newick;
    newick.enableExtendedBootstrapProperty("name");
    newick.write(treeCopy, tags, true);
  }
}
  
/******************************************************************************/
  
ProbabilisticSubstitutionMapping * CoETools::getVectors(
  const Alphabet * alphabet,
        TreeTemplate<Node> & tree,
  const SiteContainer & completeSites,
  const SiteContainer & sites,
        SubstitutionModel & model,
        DiscreteDistribution & rDist,
  const SubstitutionCount & substitutionCount,
  map<string, string> & params,
  const string & suffix)
{
  ProbabilisticSubstitutionMapping * substitutions = NULL;
  string inputVectorsFilePath = ApplicationTools::getAFilePath("input.vectors.file", params, false, false, suffix, false);

  if(inputVectorsFilePath != "none")
  {
    //We try to load the substitution vector directly from file:
    int nbSites = completeSites.getNumberOfSites();
    ApplicationTools::displayResult("Substitution mapping in file:", inputVectorsFilePath);
    ifstream sc(inputVectorsFilePath.c_str(), ios::in);
    substitutions = new ProbabilisticSubstitutionMapping(tree, nbSites);
    SubstitutionMappingTools::readFromStream(sc, *substitutions);
  }
  else
  {
    //We compute the substitutions vector:

    string outputVectorsFilePath = ApplicationTools::getAFilePath("output.vectors.file", params, true, false, suffix, false);
    ApplicationTools::displayResult("Output mapping to file" + suffix, outputVectorsFilePath);

    ApplicationTools::displayMessage("Compute all substitution numbers for each site.");
    
    bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true);
    bool joint   = ApplicationTools::getBooleanParameter("nijt.joint"  , params, true);
    DRHomogeneousTreeLikelihood drhtl(tree, completeSites, &model, &rDist, true);
    drhtl.initialize();
    if(average)
    {
      if(joint)
      {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectors(drhtl, substitutionCount);
      }
      else
      {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drhtl, substitutionCount);
      }
    }
    else
    {
      if(joint)
      {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drhtl, substitutionCount);
      }
      else
      {
        substitutions = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drhtl, substitutionCount);
      }
    }
    if(outputVectorsFilePath != "none")
    {
      ofstream outputVectors(outputVectorsFilePath.c_str(), ios::out);
      SubstitutionMappingTools::writeToStream(*substitutions, completeSites, outputVectors);
      outputVectors.close();
    }

  }
  return substitutions;
}

/******************************************************************************/

int CoETools::getMinRateClass(map<string, string> & params, string suffix)
{
  int minRateClass = ApplicationTools::getIntParameter("statistic.min_rate_class", params, 0, suffix, true);
  if(minRateClass > 0)
    ApplicationTools::displayMessage(
        "Only sites with posterior rate class >= " +
        TextTools::toString(minRateClass) +
        " will be compared.");
  return minRateClass;  
}

/******************************************************************************/

double CoETools::getMinRate(map<string, string> & params, string suffix)
{
  double minRate = ApplicationTools::getDoubleParameter("statistic.min_rate", params, 0., suffix, true);
  if(minRate > 0.)
    ApplicationTools::displayMessage(
        "Only sites with posterior rate > = " +
        TextTools::toString(minRate) +
        " will be compared.");
  return minRate;
}

/******************************************************************************/

int CoETools::getMaxRateClassDiff(map<string, string> & params)
{
  int maxRateClassDiff = ApplicationTools::getIntParameter("statistic.max_rate_class_diff", params, -1);
  if(maxRateClassDiff >= 0) 
    ApplicationTools::displayMessage(
      "Only pairs of sites with difference in posterior rate class <= " +
      TextTools::toString(maxRateClassDiff) +
      " will be compared.");
  return maxRateClassDiff;  
}

/******************************************************************************/

double CoETools::getMaxRateDiff(map<string, string> & params)
{
  double maxRateDiff = ApplicationTools::getDoubleParameter("statistic.max_rate_diff", params, -1.);
  if(maxRateDiff >= 0.)
    ApplicationTools::displayMessage(
        "Only pairs of sites with difference in posterior rate <= " +
        TextTools::toString(maxRateDiff) +
        " will be compared.");
  return maxRateDiff;
}

/******************************************************************************/

double CoETools::getStatisticMin(map<string, string> & params)
{
  double minStatistic = ApplicationTools::getDoubleParameter("statistic.min", params, 0);
  if(minStatistic > 0)
    ApplicationTools::displayMessage(
      "Only pairs of sites with abs(statistic) >= " +
      TextTools::toString(minStatistic) +
      " will be written.");
  return minStatistic;
}

/******************************************************************************/

bool CoETools::haveToPerformIndependantComparisons(map<string, string> & params) {
  bool indepComp = ApplicationTools::getBooleanParameter("independant_comparisons", params, false);
  if(indepComp)
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
  if(outFile == "none") return;

  // Get the rate class with maximum posterior probability:
  vector<unsigned int> classes = ras.getRateClassWithMaxPostProbOfEachSite();
  // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
  Vdouble rates = ras.getPosteriorRateOfEachSite();
  Vdouble logLn = ras.getLogLikelihoodForEachSite();

  ApplicationTools::displayResult("Alignment information logfile", outFile);

  ofstream out(outFile.c_str(), ios::out);
  out << "is.complete\tis.constant\trc\tpr\tlogLn" << endl;

  for(unsigned int i = 0; i < completeSites.getNumberOfSites(); i++) {
    const Site * currentSite = completeSites.getSite(i);
    int currentSitePosition = currentSite -> getPosition();
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

const Statistic * CoETools::getStatistic(map<string, string> & params)
{
  string statistic = ApplicationTools::getStringParameter("statistic", params, "none");
  if(statistic == "cosinus") {
    return new CosinusStatistic();
  } else if(statistic == "correlation") {
    return new CorrelationStatistic();
  } else if(statistic == "covariance") {
    return new CovarianceStatistic();
  } else if(statistic == "cosubstitution") {
    return new CosubstitutionNumberStatistic();
  } else {
    return NULL;
  }
}

/******************************************************************************/

SubstitutionCount * CoETools::getSubstitutionCount(
  const Alphabet * alphabet,
  const TreeTemplate<Node> & tree,
  const MutationProcess & process,
  const DiscreteDistribution & rDist,
  map<string, string> & params,
  string suffix)
{
  SubstitutionCount * substitutionCount = NULL;
  string nijtOption = ApplicationTools::getStringParameter("nijt", params, "simule", suffix, true);

   if(nijtOption == "laplace")
  {
    int trunc = ApplicationTools::getIntParameter("nijt_laplace.trunc", params, 10, suffix, true);
    substitutionCount = new AnalyticalSubstitutionCount(process.getSubstitutionModel(), trunc);
  }
  else if(nijtOption == "simple")
  {
    substitutionCount = new SimpleSubstitutionCount(alphabet);
  }
  else if(nijtOption == "aadist")
  {
    if(!AlphabetTools::isProteicAlphabet(alphabet))
    {
      ApplicationTools::displayError("Chemical distance can only be used with protein data.");
      exit(-1);
    }
    string dist = ApplicationTools::getStringParameter("nijt_aadist.type", params, "grantham", suffix, true);
    bool sym = ApplicationTools::getBooleanParameter("nijt_aadist.sym", params, true, suffix, true);
    if(dist == "grantham")
    {
      GranthamAAChemicalDistance * M = new GranthamAAChemicalDistance();
      M -> setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if(dist == "miyata")
    {
      MiyataAAChemicalDistance * M = new MiyataAAChemicalDistance();
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if(dist == "grantham.polarity")
    {
      GranthamAAPolarityIndex I;
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if(dist == "grantham.volume")
    {
      GranthamAAVolumeIndex I;
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if(dist == "klein.charge")
    {
      KleinAANetChargeIndex I;
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if(dist == "charge")
    {
      AAChargeIndex I;
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if(dist == "user1")
    {
      string aax1FilePath = ApplicationTools::getAFilePath("nijt_aadist.type_user1.file", params, true, true, suffix, false);
      ifstream aax1File(aax1FilePath.c_str(), ios::in);
      AAIndex1Entry I(aax1File);
      SimpleIndexDistance<double> * M = new SimpleIndexDistance<double>(I);
      aax1File.close();
      M->setSymmetric(sym);
      substitutionCount = new IndexToCount(M, true); // M will be deleted when this substitutionsCount will be deleted.
    }
    else if(dist == "user2")
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
  const DiscreteRatesAcrossSitesTreeLikelihood & tl,
  const SiteContainer & completeSites,
  ProbabilisticSubstitutionMapping & mapping,
  const Statistic & statistic,
  map<string, string> & params)
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

  statOut << "Group\tstatistic\tmin.rc\tmin.pr\tNmin" << endl;

  ApplicationTools::displayTask("Analyse each site pair");
  
  vector<unsigned int> classes = tl.getRateClassWithMaxPostProbOfEachSite();
  Vdouble rates   = tl.getPosteriorRateOfEachSite();

  for(unsigned int i = 0; i < nbSites; i++)
  {
    int    iClass = classes[i];
    double iRate  = rates[i];
    if(iClass < minRateClass) continue;
    if(iRate  < minRate     ) continue;
    double iNorm  = norms[i];
    *ApplicationTools::message << ".";
    ApplicationTools::message->flush();
    for(unsigned int j = i + 1; j < nbSites; j++)
    {
      int    jClass = classes[j];
      double jRate  = rates[j];
      if(jClass < minRateClass) continue;
      if(jRate  < minRate     ) continue;
      double jNorm  = norms[j];
    
      //Sites which are in too different rate classes are not compared:
      if(maxRateClassDiff >= 0  && NumTools::abs(jClass - iClass) > maxRateClassDiff) continue;
      if(maxRateDiff      >= 0. && NumTools::abs(jRate  - iRate ) > maxRateDiff)      continue;

      double stat = statistic.getValueForPair(mapping[i], mapping[j]);
      if(NumTools::abs(stat) < minStatistic) continue;

      //Then print to file:
      statOut << "[";
      statOut << completeSites.getSite(i)->getPosition();
      statOut << ",";
      statOut << completeSites.getSite(j)->getPosition();
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
  if(indepComp && mapping1.getNumberOfSites() != mapping2.getNumberOfSites())
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

  statOut << "Group\tstatistic\tmin.rc\tmin.pr\tNmin" << endl;

  ApplicationTools::displayTask("Analyse each site pair");
    
  vector<unsigned int> classes1 = tl1.getRateClassWithMaxPostProbOfEachSite();
  vector<unsigned int> classes2 = tl2.getRateClassWithMaxPostProbOfEachSite();
  Vdouble rates1   = tl1.getPosteriorRateOfEachSite();
  Vdouble rates2   = tl2.getPosteriorRateOfEachSite();

  for(unsigned int i = 0; i < nbSites1; i++)
  {
    int    iClass = classes1[i];
    double iRate  = rates1[i];
    if(iClass < minRateClass1) continue;
    if(iRate  < minRate1     ) continue;
    double iNorm  = norms1[i];
    *ApplicationTools::message << ".";
    ApplicationTools::message->flush();
      
    unsigned int begin = indepComp ? i : 0;
    unsigned int end   = indepComp ? i + 1 : nbSites2;
    for(unsigned int j = begin; j < end; j++)
    {
      int    jClass = classes2[j];
      double jRate  = rates2[j];
      if(jClass < minRateClass2) continue;
      if(jRate  < minRate2     ) continue;
      double jNorm  = norms2[i];
    
      //Sites which are in too different rate classes are not compared:
      if(maxRateClassDiff >= 0  && NumTools::abs(jClass - iClass) > maxRateClassDiff) continue;
      if(maxRateDiff      >= 0. && NumTools::abs(jRate  - iRate ) > maxRateDiff     ) continue;

      //double stat = table[i + 1][j + 1];
      double stat = statistic.getValueForPair(mapping1[i], mapping2[j]);
      if(NumTools::abs(stat) < minStatistic) continue;

      //Then print to file:
      statOut << "[";
      statOut << completeSites1.getSite(i)->getPosition();
      statOut << ",";
      statOut << completeSites2.getSite(j)->getPosition();
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
  const MutationProcess & process,
  const DiscreteDistribution & rDist,
  const TreeTemplate<Node> & tree,
  const SubstitutionCount & nijt,
  const Statistic & statistic,
  map<string, string> & params,
  bool useContinuousRates)
{
  HomogeneousSequenceSimulator seqSim(&process, &rDist, &tree);
  seqSim.enableContinuousRates(useContinuousRates);
  string path = ApplicationTools::getAFilePath("statistic.null.output.file", params, true, false);
  ofstream outFile(path.c_str(), ios::out);
  
  ApplicationTools::displayMessage("Compute statistic under null hypothesis...");
  
  unsigned int nbRepCPU = ApplicationTools::getParameter<unsigned int>("statistic.null.nb_rep_CPU", params, 10);
  
  // Drop this for now:
  //bool reestimate = ApplicationTools::getBooleanParameter("statistic.null.reestimate", params, true);
  //if(!reestimate) {
  //  AnalysisTools::getNullDistributionIntraWithoutReestimatingCounts(seqSim, statistic, outFile, nbRepCPU, true);
  //  outFile.close();
  //  ApplicationTools::displayTaskDone();
  //  return;
  //}
  
  unsigned int nbRepRAM = ApplicationTools::getParameter<unsigned int>("statistic.null.nb_rep_RAM", params, 100);
  bool cumul   = ApplicationTools::getBooleanParameter("statistic.null.cumul", params, false);
  bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true);
  bool joint   = ApplicationTools::getBooleanParameter("nijt.joint", params, true);

  if(cumul)
  {
    // Building domain:
    double lowerSB = ApplicationTools::getDoubleParameter("statistic.null.lower", params, -1.);
    double upperSB = ApplicationTools::getDoubleParameter("statistic.null.upper", params, 1.);
    int     nbSInt = ApplicationTools::getIntParameter   ("statistic.null.nb_int", params, 20);
    Domain statDomain(lowerSB, upperSB, nbSInt);
    Domain rateDomain = rDist.getDomain();
  
    // Simulate:
    vector<IntervalData *> id; 
    id = AnalysisTools::getNullDistributionIntraDR(seqSim, nijt, statistic, statDomain, rateDomain, nbRepCPU, nbRepRAM, average, joint, true);
  
    // Print to file:
    for(unsigned int i = 0; i < rateDomain.getSize(); i++)
    {
      outFile << "# Distribution with minimum rate " << rateDomain.getValue(i) << endl;
      id[i]->print(outFile);
      outFile << endl;
    }

    // Free memory:
    for(unsigned int i = 0; i < id.size(); i++) delete id[i];
  }
  else
  {
    AnalysisTools::getNullDistributionIntraDR(seqSim, nijt, statistic, outFile, nbRepCPU, nbRepRAM, average, joint, true);
  }
  outFile.close();
}

/******************************************************************************/

void CoETools::computeInterNullDistribution(
  const MutationProcess & process1,
  const MutationProcess & process2,
  const DiscreteDistribution & rDist1,
  const DiscreteDistribution & rDist2,
  const TreeTemplate<Node> & tree1,
  const TreeTemplate<Node> & tree2,
  const SubstitutionCount & nijt1,
  const SubstitutionCount & nijt2,
  const Statistic & statistic,
  map<string, string> & params,
  bool useContinuousRates)
{
  HomogeneousSequenceSimulator seqSim1(&process1, &rDist1, &tree1);
  seqSim1.enableContinuousRates(useContinuousRates);
  HomogeneousSequenceSimulator seqSim2(&process2, &rDist2, &tree2);
  seqSim2.enableContinuousRates(useContinuousRates);
  string path = ApplicationTools::getAFilePath("statistic.null.output.file", params, true, false);
  ofstream outFile(path.c_str(), ios::out);
  
  int nbRepCPU = ApplicationTools::getIntParameter("statistic.null.nb_rep_CPU", params, 10);
  int nbRepRAM = ApplicationTools::getIntParameter("statistic.null.nb_rep_RAM", params, 100);
  bool cumul   = ApplicationTools::getBooleanParameter("statistic.null.cumul", params, false);
  bool average = ApplicationTools::getBooleanParameter("nijt.average", params, true);
  bool joint   = ApplicationTools::getBooleanParameter("nijt.joint", params, true);

  ApplicationTools::displayMessage("Compute statistic under null hypothesis...");
  if(cumul)
  {
    // Building domain:
    double lowerSB = ApplicationTools::getDoubleParameter("statistic.null.lower", params, -1.);
    double upperSB = ApplicationTools::getDoubleParameter("statistic.null.upper", params, 1.);
    int     nbSInt = ApplicationTools::getIntParameter   ("statistic.null.nb_int", params, 20);
    Domain statDomain(lowerSB, upperSB, nbSInt);
    Vdouble bounds = ApplicationTools::getVectorParameter<double>("statistic.null.rate.bounds", params, ',', "");
    Domain rateDomain = Domain(bounds);
    ApplicationTools::displayMessage("The following rate domain will be used:");
    for(unsigned int i = 0; i < rateDomain.getSize(); i++)
    {
      *ApplicationTools::message
        << "["
        << rateDomain.getBound(i)
        << ", "
        << rateDomain.getBound(i+1)
        << "[, "
        << rateDomain.getValue(i)
        << endl;
    }
  
    // Simulate:
    vector<IntervalData *> id;
    id = AnalysisTools::getNullDistributionInterDR(
      seqSim1, seqSim2,
      nijt1, nijt2,
      statistic,
      statDomain, rateDomain,
      nbRepCPU, nbRepRAM,
      average, joint,
      true);
  
    // Print to file:
    for(unsigned int i = 0; i < rateDomain.getSize(); i++)
    {
      outFile << "# Distribution with minimum rate " << rateDomain.getValue(i) << endl;
      id[i]->print(outFile);
      outFile << endl;
    }
  
    // Free memory:
    for(unsigned int i = 0; i < id.size(); i++) delete id[i];
  }
  else
  {
    AnalysisTools::getNullDistributionInterDR(seqSim1, seqSim2, nijt1, nijt2, statistic, outFile, nbRepCPU, nbRepRAM, average, joint, true);
  }
  outFile.close();
}

/******************************************************************************/
