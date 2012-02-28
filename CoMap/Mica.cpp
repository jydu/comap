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

#include <Bpp/Utils/AttributesTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/Prob.all>
#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
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
#include <Bpp/Phyl/Mapping.all>

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

void miTest(const Site& site1, const Site& site2, unsigned int maxNbPermutations, double& mi, double& pvalue, unsigned int& nbPermutations)
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
    unsigned int i;
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
  unsigned int nbVectors = mapping.getNumberOfSites();
  vector<double> vect(nbVectors);
  for (unsigned int i = 0; i < nbVectors; i++)
    vect[i] = SubstitutionMappingTools::computeNormForSite(mapping, i);
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
	cout << "* This is Mica         version 1.0.0       date: 28/03/11 *" << endl;
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

  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(mica.getParams(), "", false);

  VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, mica.getParams());
  
  VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, mica.getParams());
  delete allSites;

  bool removeConst = ApplicationTools::getBooleanParameter("input.remove_const", mica.getParams(), true);
  if (removeConst) {
    unsigned int n = sites->getNumberOfSites();
    for (unsigned int i = n; i > 0; --i) {
      if (SiteTools::isConstant(sites->getSite(i - 1), true))
        sites->deleteSite(i - 1);
    }
    ApplicationTools::displayResult("Number of conserved sites ignored", n - sites->getNumberOfSites());
  }

  unsigned int nbSites = sites->getNumberOfSites();
  unsigned int nbSeqs  = sites->getNumberOfSequences();
 
  ApplicationTools::displayResult("Number of sequences", nbSeqs);
  ApplicationTools::displayResult("Number of sites", nbSites);
  

  //Shall we use a model?
  TreeTemplate<Node>  * tree     = 0;
  DRTreeLikelihood    * tl       = 0;
  SubstitutionModel   * model    = 0;
  SubstitutionModelSet* modelSet = 0;
  DiscreteDistribution* rDist    = 0;
  SubstitutionCount   * simple   = 0;
  vector<double>      * norms    = 0;
  bool withModel = ApplicationTools::getBooleanParameter("use_model", mica.getParams(), false, "", true, false);
  
  ApplicationTools::displayBooleanResult("Model of sequence evolution", withModel);
  if (withModel)
  {
    // Get the initial tree
    Tree* tree_tmp = PhylogeneticsApplicationTools::getTree(mica.getParams());
    tree = new TreeTemplate<Node>(*tree_tmp);
    delete tree_tmp;
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));

    string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", mica.getParams(), "no", "", true, false);
    ApplicationTools::displayResult("Heterogeneous model", nhOpt);

    if (nhOpt == "no")
    {  
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, mica.getParams());
      if(model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mica.getParams());
      }
      tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, false);
    }
    else if(nhOpt == "one_per_branch")
    {
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, mica.getParams());
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mica.getParams());
      }
      vector<double> rateFreqs;
      if (model->getNumberOfStates() != alphabet->getSize())
      {
        //Markov-Modulated Markov Model...
        unsigned int n =(unsigned int)(model->getNumberOfStates() / alphabet->getSize());
        rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                     // we should assume a rate distribution for the root also!!!  
      }
      FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, sites, mica.getParams(), rateFreqs);
      vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", mica.getParams(), ',', "");
      modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
      model = modelSet->getModel(0)->clone();
      tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false);
    }
    else if (nhOpt == "general")
    {
      modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, sites, mica.getParams());
      if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mica.getParams());
      }
      model = modelSet->getModel(0)->clone();
      tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false);
    }
    else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
    tl->initialize();
 
    delete tree;
    
    double logL = tl->getValue();
    if (isinf(logL))
    {
      // This may be due to null branch lengths, leading to null likelihood!
      ApplicationTools::displayWarning("!!! Warning!!! Likelihood is zero.");
      ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
      ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
      ParameterList pl = tl->getBranchLengthsParameters();
      for (unsigned int i = 0; i < pl.size(); i++)
      {
        if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
      }
      tl->matchParametersValues(pl);
      logL = tl->getValue();
    }
    if (isinf(logL))
    {
      ApplicationTools::displayError("!!! Unexpected likelihood == 0.");
      ApplicationTools::displayError("!!! Looking at each site:");
      for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
      {
        (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
      }
      ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
      exit(-1);
    }
    tree = new TreeTemplate<Node>(tl->getTree());

    //Get the substitution mapping in order to compute the rates:

    TotalSubstitutionRegister* reg = new TotalSubstitutionRegister(model->getAlphabet());
    simple = new SimpleSubstitutionCount(reg); 
    ProbabilisticSubstitutionMapping* mapping = 
      SubstitutionMappingTools::computeSubstitutionVectors(*tl, *simple, true);
    norms = new vector<double>(computeNorms(*mapping));
    delete mapping;
  }
  
  //Compute all pairwise MI values:
  
  string path = ApplicationTools::getAFilePath("output.file", mica.getParams(), true, false);
  ApplicationTools::displayResult("Output file", path);

  vector<double> averageMI;
  vector<double> entropy;
  ApplicationTools::displayTask("Computing average MIs", true);
  for (unsigned int i = 0; i < nbSites; ++i) {
    const Site* site1 =  &sites->getSite(i);
    ApplicationTools::displayGauge(i, nbSites - 1);
    double sum = 0;
    for (unsigned int j = 0; j < nbSites; ++j) {
      if (j != i) {
        const Site* site2 =  &sites->getSite(j);
        sum += SiteTools::mutualInformation(*site1, *site2, true); 
      }
    }
    entropy.push_back(SiteTools::entropy(*site1, true));
    averageMI.push_back(sum / static_cast<double>(nbSites - 1));
  }
  ApplicationTools::displayTaskDone();
  double fullAverageMI = VectorTools::mean<double, double>(averageMI);

  //Some variables we'll need:
  const Site *site1, *site2;
  double stat, apc, rcw, nmin = 0, hj, hm, perm;
  unsigned int maxNbPermutations = 0;

  //Get the null distribution of MI values:
  string nullMethod = ApplicationTools::getStringParameter("null.method", mica.getParams(), "none", "", true, false);
  ApplicationTools::displayResult("Null distribution", nullMethod);

  bool computePValues = false;
  //These variables will only be used if p-values are computed:
  unsigned int nbRateClasses = 0;
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
      nbRateClasses = ApplicationTools::getParameter<unsigned int>("null.nb_rate_classes", mica.getParams(), 10);
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
   
      auto_ptr<ofstream> simout;
      if (outputToFile) {
        ApplicationTools::displayResult("Null output file", simpath);
        simout.reset(new ofstream(simpath.c_str(), ios::out));
        *simout << "MI\tHjoint\tHmin";
        if (withModel)
          *simout << "\tNmin";
        *simout << endl;
      }
    
      unsigned int nbRepCPU = ApplicationTools::getParameter<unsigned int>("null.nb_rep_CPU", mica.getParams(), 10);
      unsigned int nbRepRAM = ApplicationTools::getParameter<unsigned int>("null.nb_rep_RAM", mica.getParams(), 100);
  
      OutputStream* os = ApplicationTools::warning;
      ApplicationTools::warning = 0;
      for (unsigned int i = 0; i < nbRepCPU; i++)
      {
        //Generate data set:
        vector<unsigned int> index1;
        SiteContainer* sites1 = SiteContainerTools::sampleSites(*sites, nbRepRAM, &index1);
  
        vector<unsigned int> index2;
        SiteContainer* sites2 = SiteContainerTools::sampleSites(*sites, nbRepRAM, &index2);
  
        for (unsigned int j = 0; j < nbRepRAM; j++)
        {
          ApplicationTools::displayGauge(i * nbRepRAM + j, nbRepCPU * nbRepRAM - 1, '>');
          site1 = &sites1->getSite(j);
          site2 = &sites2->getSite(j);
          stat  = SiteTools::mutualInformation(*site1, *site2, true);
          hj    = SiteTools::jointEntropy(*site1, *site2, true);
          hm    = std::min(entropy[index1[j]], entropy[index2[j]]);
  
          if (outputToFile)
            *simout << stat << "\t" << hj << "\t" << hm;
          if (withModel) {
            nmin  = min((*norms)[index1[j]], (*norms)[index2[j]]);
            if (outputToFile)
              *simout << "\t" << nmin;
            if (computePValues) {
              try {
                unsigned int cat = rateDomain->getIndex(nmin);
                (*simValues)[cat].push_back(stat);
              } catch (OutOfRangeException& oore) {}
            }
          } else {
            if (computePValues) {
              try {
                unsigned int cat = rateDomain->getIndex(hm);
                (*simValues)[cat].push_back(stat);
              } catch (OutOfRangeException& oore) {}
            }
          }
          if (outputToFile)
            *simout << endl;
        }
  
        delete sites1;
        delete sites2;
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
      SequenceSimulator* simulator;
      if (modelSet)
      {
        simulator = new NonHomogeneousSequenceSimulator(modelSet, rDist, tree);
        dynamic_cast<NonHomogeneousSequenceSimulator *>(simulator)->enableContinuousRates(continuousSim);
      }
      else
      {
        simulator = new HomogeneousSequenceSimulator(model, rDist, tree);
        dynamic_cast<HomogeneousSequenceSimulator *>(simulator)->enableContinuousRates(continuousSim);
      }
   
      string simpath = ApplicationTools::getAFilePath("null.output.file", mica.getParams(), false, false);
      bool outputToFile = (simpath != "none") && (!TextTools::isEmpty(simpath));
      ApplicationTools::displayTask("Computing null distribution", true);
 
      auto_ptr<ofstream> simout;
      if (outputToFile) {
        ApplicationTools::displayResult("Null output file", simpath);
        simout.reset(new ofstream(simpath.c_str(), ios::out));
    
        *simout << "MI\tHjoint\tHmin\tNmin" << endl;
      }

      unsigned int nbRepCPU = ApplicationTools::getParameter<unsigned int>("null.nb_rep_CPU", mica.getParams(), 10);
      unsigned int nbRepRAM = ApplicationTools::getParameter<unsigned int>("null.nb_rep_RAM", mica.getParams(), 100);
  
      OutputStream* os = ApplicationTools::warning;
      ApplicationTools::warning = 0;
      for (unsigned int i = 0; i < nbRepCPU; i++)
      {
        //Generate data set:
        SiteContainer* sites1 = simulator->simulate(nbRepRAM);
        tl->setData(*sites1);
        tl->initialize();
        vector<double> norms1;
        ProbabilisticSubstitutionMapping* mapping1 = 
           SubstitutionMappingTools::computeSubstitutionVectors(*tl, *simple, false);
        norms1 = computeNorms(*mapping1);
        delete mapping1;
  
        SiteContainer* sites2 = simulator->simulate(nbRepRAM);
        tl->setData(*sites2);
        tl->initialize();
        vector<double> norms2;
        ProbabilisticSubstitutionMapping* mapping2 = 
           SubstitutionMappingTools::computeSubstitutionVectors(*tl, *simple, false);
        norms2 = computeNorms(*mapping2);
        delete mapping2;
  
        for (unsigned int j = 0; j < nbRepRAM; j++)
        {
          ApplicationTools::displayGauge(i * nbRepRAM + j, nbRepCPU * nbRepRAM - 1, '>');
          site1 = &sites1->getSite(j);
          site2 = &sites2->getSite(j);
          stat  = SiteTools::mutualInformation(*site1, *site2, true);
          hj    = SiteTools::jointEntropy(*site1, *site2, true);
          hm    = std::min(entropy[i], entropy[j]);
          nmin  = min(norms1[j], norms2[j]);
  
          if (computePValues) {
            try {
              unsigned int cat = rateDomain->getIndex(nmin);
              (*simValues)[cat].push_back(stat);
            } catch(OutOfRangeException& oore) {}
          }

          if (outputToFile)
            *simout << stat << "\t" << hj << "\t" << hm << "\t" << nmin << endl;
        }
  
        delete sites1;
        delete sites2;
      }
      ApplicationTools::warning = os;
  
      if (outputToFile)
        simout->close();
      ApplicationTools::displayTaskDone(); 
      delete simulator;
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
      OutputStream* os = ApplicationTools::warning;
      ApplicationTools::warning = 0;
      unsigned int c = 0;
      for (unsigned int i = 0; i < nbSites - 1; i++)
      {
        site1 = &sites->getSite(i);
        for (unsigned int j = i + 1; j < nbSites; j++)
        {
          site2 = &sites->getSite(j);
          ApplicationTools::displayGauge(c++, nbSites * (nbSites - 1) / 2 - 1, '>');
          stat  = SiteTools::mutualInformation(*site1, *site2, true);
          hj    = SiteTools::jointEntropy(*site1, *site2, true);
          hm    = std::min(entropy[i], entropy[j]);
  
          if (withModel) {
            nmin  = min((*norms)[i], (*norms)[j]);
            try {
              unsigned int cat = rateDomain->getIndex(nmin);
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
              unsigned int cat = rateDomain->getIndex(hm);
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
      maxNbPermutations = ApplicationTools::getParameter<unsigned int>("null.max_number_of_permutations", mica.getParams(), 1000);
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
    for (unsigned int i = 0; i < nbRateClasses; ++i) {
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

  unsigned int nbPerm;
  unsigned int c = 0;
  for (unsigned int i = 0; i < nbSites - 1; ++i)
  {
    site1 = &sites->getSite(i);
    for (unsigned int j = i + 1; j < nbSites; ++j)
    {
      ApplicationTools::displayGauge(c++, (nbSites - 1) * (nbSites) / 2 - 1);
      site2 = &sites->getSite(j);
      miTest(*site1, *site2, maxNbPermutations, stat, perm, nbPerm);

      out << "[" << site1->getPosition() << ";" << site2->getPosition() << "]\t" << stat << "\t";
      apc = averageMI[i] * averageMI[j] / fullAverageMI;
      rcw = averageMI[i] * averageMI[j] / 2.;
      hj  = SiteTools::jointEntropy(*site1, *site2, true);
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
          unsigned int cat = withModel ? rateDomain->getIndex(nmin) : rateDomain->getIndex(hm);
          unsigned int nsim = (*simValues)[cat].size();
          unsigned int count;
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

  delete alphabet;
  delete sites;
  delete tl;
  if (model)    delete model;
  if (modelSet) delete modelSet;
  delete rDist;
  delete tree;

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
