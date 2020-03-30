//
// File: CoMap.cpp
// Created by: Julien Dutheil
// Created on: Thu Apr 6 11:06 2006
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004-2018)

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

// From bpp-core:
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/App/BppApplication.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/Likelihood/MarginalAncestralStateReconstruction.h>
#include <Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Mapping/WeightedSubstitutionCount.h>
#include <Bpp/Phyl/AncestralStateReconstruction.h>
#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include <Bpp/Phyl/Io/PhylipDistanceMatrixFormat.h>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

// CoMap's include files:
#include "CoETools.h"
#include "Distance.h"
#include "Cluster.h"
#include "ClusterTools.h"

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "comap parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "      Refer to the CoMap Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
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
  cout << "* This is CoMap        version 1.5.5       date: 25/04/18 *" << endl;
  cout << "*     A C++ shell program to detect co-evolving sites.    *" << endl;
  cout << "***********************************************************" << endl;
  cout << endl;
  
  try
  {
  
  // **************************
  // * Retrieving parameters: *
  // **************************
  
  if(argc == 1)
  { 
    // No argument, show display some help and leave.
    help();
    return 0;
  }

  BppApplication comap(argc, argv, "CoMap");
  comap.startTimer();

  ApplicationTools::displayMessage("\n\n-*- Retrieve data and model -*-\n");

  unique_ptr<Tree> tmpTree(PhylogeneticsApplicationTools::getTree(comap.getParams()));
  shared_ptr<TreeTemplate<Node>> tree1(new TreeTemplate<Node>(*tmpTree));
  tmpTree.reset();
  ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree1->getNumberOfLeaves()));
  ApplicationTools::displayResult("Number of sons at root", TextTools::toString(tree1->getRootNode()->getNumberOfSons()));
  
  // Get data 1:
  //Tree * tree1 = new Tree(* tree); // Copy tree.
  shared_ptr<Alphabet> alphabet1;
  shared_ptr<GeneticCode> geneticCode1;
  shared_ptr<VectorSiteContainer> allSites1, sites1;
  shared_ptr<SubstitutionModel> model1;
  shared_ptr<SubstitutionModelSet> modelSet1;
  shared_ptr<DiscreteDistribution> rDist1;
  shared_ptr<DRTreeLikelihood> tl1;
  CoETools::readData(tree1, alphabet1, geneticCode1, allSites1, sites1, model1, modelSet1, rDist1, tl1, comap.getParams(), "");
 
  ApplicationTools::displayResult("Number of sites in file", allSites1->getNumberOfSites());
  ApplicationTools::displayResult("Number of sites to analyse", sites1->getNumberOfSites());
  ApplicationTools::displayResult("Number of site patterns", tl1->getLikelihoodData()->getNumberOfDistinctSites());
  
  bool continuousSim = ApplicationTools::getBooleanParameter("simulations.continuous", comap.getParams(), false, "", true, 1);
  ApplicationTools::displayResult("Rate distribution for simulations", (continuousSim ? "continuous" : "discrete"));

  ApplicationTools::displayMessage("\n\n-*- Get substitution vectors -*-\n");

  // Getting the substitutions count function:
  shared_ptr<SubstitutionCount> substitutionCount1(PhylogeneticsApplicationTools::getSubstitutionCount(alphabet1.get(), model1.get(), comap.getParams(), "", true, 1));
    
  // Getting the substitution vectors:
  shared_ptr<ProbabilisticSubstitutionMapping> mapping1(CoETools::getVectors(*tl1, *substitutionCount1, *sites1, comap.getParams()));

  // Compute norms:
  size_t nbSites1 = mapping1->getNumberOfSites();
  vector<double> norms1(nbSites1);
  for (size_t i = 0; i < nbSites1; i++)
  {
    norms1[i] = SubstitutionMappingTools::computeNormForSite(*mapping1, i);
  }
  CoETools::writeInfos(*sites1, *tl1, norms1, comap.getParams());

  string analysis = ApplicationTools::getStringParameter("analysis", comap.getParams(), "pairwise", "", false, 0);
  ApplicationTools::displayResult("Analysis type", analysis);

  // Reconstruct ancestral sequences (not used in the analysis):
  string reconstruction = ApplicationTools::getStringParameter("asr.method", comap.getParams(), "none", "", true, 1);
  ApplicationTools::displayResult("Ancestral state reconstruction method", reconstruction);

  unique_ptr<AncestralStateReconstruction> asr;
  unique_ptr<DRTreeLikelihood> tl1copy;
  if (reconstruction == "none") {
    //do nothing
  } else if (reconstruction == "marginal") {
    //We reconstruct ancestral states for selected sites only.
    //Corresponding positions in the original alignment can be obtained from the "counts" output file.
    tl1copy.reset(tl1->clone());
    tl1copy->setData(*sites1);
    tl1copy->initialize();
    asr.reset(new MarginalAncestralStateReconstruction(tl1copy.get()));
  } else
    throw Exception("Unknown ancestral state reconstruction method: " + reconstruction);

  if (asr) {
    unique_ptr<SiteContainer> asSites1(asr->getAncestralSequences());
  
    //Add existing sequence to output:
    SequenceContainerTools::append(*asSites1, *sites1);

    //Write output:
    if (ApplicationTools::getStringParameter("output.sequence.file", comap.getParams(), "none", "", true, 1) != "none") {
      SequenceApplicationTools::writeAlignmentFile(*asSites1, comap.getParams());
    }
  }



  // Coevolution analysis:
  if (analysis == "none")
  {
    //No coevolution analysis, only perform mapping!
  }
  else
  {
    // Building a simulator object:
    shared_ptr<SequenceSimulator> seqSim1;
    if (modelSet1)
    {
      seqSim1.reset(new NonHomogeneousSequenceSimulator(modelSet1.get(), rDist1.get(), tree1.get()));
      dynamic_cast<NonHomogeneousSequenceSimulator *>(seqSim1.get())->enableContinuousRates(continuousSim);
    }
    else
    {
      seqSim1.reset(new HomogeneousSequenceSimulator(model1.get(), rDist1.get(), tree1.get()));
      dynamic_cast<HomogeneousSequenceSimulator *>(seqSim1.get())->enableContinuousRates(continuousSim);
    }

    // *********************
    // * Pairwise analysis *
    // *********************

    if (analysis == "pairwise")
    {
      shared_ptr<Statistic> statistic(CoETools::getStatistic(comap.getParams(), alphabet1.get(), substitutionCount1.get()));

      bool computeNullHyp = false;
      computeNullHyp = ApplicationTools::getBooleanParameter("statistic.null", comap.getParams(), true, "", false, 1);
  
  
      // *******************************************
      // * The same for a putative second dataset: *
      // *******************************************
      if (comap.getParams().find("input.sequence.file2") != comap.getParams().end() && comap.getParams()["input.sequence.file2"] != "none")
      {
        shared_ptr<TreeTemplate<Node>> tree2;
        if (comap.getParams().find("input.tree.file2") != comap.getParams().end() && comap.getParams()["input.tree.file2"] != "none")
        {
          unique_ptr<Tree> tmpTree2(PhylogeneticsApplicationTools::getTree(comap.getParams(), "input.", "2", true));
          tree2.reset(new TreeTemplate<Node>(*tmpTree2));
          if (!TreeTools::haveSameTopology(*tree1, *tree2))
            throw Exception("The second tree must have the same topology as the first tree.");
          ApplicationTools::displayResult("# number of leaves", TextTools::toString(tree2->getNumberOfLeaves()));
          ApplicationTools::displayResult("# number of sons at root", TextTools::toString(tree2->getRootNode()->getNumberOfSons()));
          //We haverage weights:
        }
        else
        {
          tree2.reset(new TreeTemplate<Node>(*tree1)); // Copy tree.
        }

        shared_ptr<Alphabet> alphabet2;
        shared_ptr<GeneticCode> geneticCode2;
        shared_ptr<VectorSiteContainer> allSites2, sites2;
        shared_ptr<SubstitutionModel> model2;
        shared_ptr<SubstitutionModelSet> modelSet2;
        shared_ptr<DiscreteDistribution> rDist2;
        shared_ptr<DRTreeLikelihood> tl2;
        ApplicationTools::displayMessage("\nLoading second dataset...\n");
        CoETools::readData(tree2, alphabet2, geneticCode2, allSites2, sites2, model2, modelSet2, rDist2, tl2, comap.getParams(), "2");
        ApplicationTools::displayResult("Number of sites in file", allSites2->getNumberOfSites());
        ApplicationTools::displayResult("Number of sites to analyse", sites2->getNumberOfSites());
        ApplicationTools::displayResult("Number of site patterns", tl2->getLikelihoodData()->getNumberOfDistinctSites());

        ApplicationTools::displayMessage("\n... and get its substitution vectors.\n");
    
        // Building a simulator object:
        shared_ptr<SequenceSimulator> seqSim2;
        if (modelSet2)
        {
          seqSim2.reset(new NonHomogeneousSequenceSimulator(modelSet2.get(), rDist2.get(), tree2.get()));
          dynamic_cast<NonHomogeneousSequenceSimulator *>(seqSim2.get())->enableContinuousRates(continuousSim);
        }
        else
        {
          seqSim2.reset(new HomogeneousSequenceSimulator(model2.get(), rDist2.get(), tree2.get()));
          dynamic_cast<HomogeneousSequenceSimulator *>(seqSim2.get())->enableContinuousRates(continuousSim);
        }

        // Getting the substitutions count function:
        SubstitutionCount* substitutionCount2 = PhylogeneticsApplicationTools::getSubstitutionCount(alphabet2.get(), model2.get(), comap.getParams());
    
        // Getting the substitution vectors:
        ProbabilisticSubstitutionMapping* mapping2 = CoETools::getVectors(*tl2, *substitutionCount2, *sites2, comap.getParams(), "2");
 
        // Compute norms:
        size_t nbSites2 = mapping1->getNumberOfSites();
        vector<double> norms2(nbSites2);
        for (size_t i = 0; i < nbSites2; i++)
        {
          norms2[i] = SubstitutionMappingTools::computeNormForSite(*mapping2, i);
        }

        CorrectedCorrelationStatistic* cstat = dynamic_cast<CorrectedCorrelationStatistic*>(statistic.get());
        if (cstat) {
          Vdouble mv1(mapping1->getNumberOfBranches());
          //Compute mean vector:
          for (size_t i = 0; i < mapping1->getNumberOfSites(); ++i) {
            mv1 += SubstitutionMappingTools::computeTotalSubstitutionVectorForSitePerBranch(*mapping1, i);
          }
          mv1 /= static_cast<double>(mapping1->getNumberOfSites());
          Vdouble mv2(mapping2->getNumberOfBranches());
          //Compute mean vector:
          for (size_t i = 0; i < mapping2->getNumberOfSites(); ++i) {
            mv2 += SubstitutionMappingTools::computeTotalSubstitutionVectorForSitePerBranch(*mapping2, i);
          }
          mv2 /= static_cast<double>(mapping2->getNumberOfSites());
          cstat->setMeanVectors(mv1, mv2);
        }
   
        ApplicationTools::displayMessage("\n\n-*- Compute statistics -*-\n");
    
        // *************************
        // * Now the stat stuff... *
        // *************************
    
        ApplicationTools::displayMessage("Compares data set 1 with data set 2.");
        CoETools::computeInterStats(
          *tl1,
          *tl2,
          *sites1,
          *sites2,
          *mapping1,
          *mapping2,
          *statistic,
          comap.getParams());
      
        CoETools::writeInfos(*sites2, *tl2, norms2, comap.getParams(), "2");
      
        //tl2 will be modified after the simulations!
        if (computeNullHyp)
        {
          CoETools::computeInterNullDistribution(
            *tl1,
            *tl2,
            *seqSim1,
            *seqSim2,
            *substitutionCount1,
            *substitutionCount2,
            *statistic,
            comap.getParams());
        }

        ApplicationTools::displayTaskDone();
      }
      else
      {
        CorrectedCorrelationStatistic* cstat = dynamic_cast<CorrectedCorrelationStatistic*>(statistic.get());
        if (cstat) {
          Vdouble mv1(mapping1->getNumberOfBranches());
          //Compute mean vector:
          for (size_t i = 0; i < mapping1->getNumberOfSites(); ++i) {
            mv1 += SubstitutionMappingTools::computeTotalSubstitutionVectorForSitePerBranch(*mapping1, i);
          }
          mv1 /= static_cast<double>(mapping1->getNumberOfSites());
          cstat->setMeanVector(mv1);
        }

        ApplicationTools::displayMessage("\n\n-*- Perform pairwise analysis -*-\n");
        
        CoETools::computeIntraStats(
          *tl1,
          *seqSim1,
          *sites1,
          *mapping1,
          *substitutionCount1,
          *statistic,
          computeNullHyp,
          comap.getParams());
      }
    }






    // ***********************
    // * Clustering analysis *
    // ***********************
  
    else if (analysis == "clustering")
    {
      ApplicationTools::displayMessage("\n\n-*- Perform clustering analysis -*-\n");
  
      string clusteringMethod = ApplicationTools::getStringParameter("clustering.method", comap.getParams(), "complete", "", false, 0);
      if (clusteringMethod != "none")
      { 
        // Get site positions:
        vector<string> siteNames(nbSites1);
        for (size_t i = 0; i < nbSites1; i++)
        {
          string siteName = TextTools::toString(sites1->getSite(i).getPosition());
          siteNames[i] = siteName; 
        }
  
        // Get the distance to use
    
        string distanceMethod = ApplicationTools::getStringParameter("clustering.distance", comap.getParams(), "cor", "", true, 1);
        unique_ptr<Distance> dist;
        if (distanceMethod == "Euclidian" || distanceMethod == "euclidian")
        {
          dist.reset(new EuclidianDistance());
        }
        else if (distanceMethod == "Correlation" || distanceMethod == "cor")
        {
          Statistic* cor = new CorrelationStatistic();
          dist.reset(new StatisticBasedDistance(cor, 1.));
        }
        else if (distanceMethod == "Compensation" || distanceMethod == "comp")
        {
          WeightedSubstitutionCount* wsc = dynamic_cast<WeightedSubstitutionCount*>(substitutionCount1.get());
          if (!wsc)
            throw Exception("Compensation distance must be used with a mapping procedure allowing weights, e.g. 'nijt=Uniformization(weight=Diff(index1=Volume, symmetrical=no))'.");
          if (!wsc->hasWeights())
            throw Exception("Compensation distance must be used with a mapping procedure with weights, e.g. 'nijt=Uniformization(weight=Diff(index1=Volume, symmetrical=no))'.");
          if (wsc->getWeights()->isSymmetric())
            throw Exception("Compensation distance must be used with a mapping procedure allowing non-symmetric weights, e.g. 'nijt=Uniformization(weight=Diff(index1=Volume, symmetrical=no))'.");
          dist.reset(new CompensationDistance());
        }
        else
        {
          ApplicationTools::displayError("Unknown distance method.");
          exit(-1);
        }
        //dist->setWeights(weights);
        ApplicationTools::displayResult("Distance to use", distanceMethod);

        // Compute the distance matrix.
        unique_ptr<DistanceMatrix> mat(new DistanceMatrix(siteNames));
        for (size_t i = 0; i < nbSites1; ++i)
        {
          (*mat)(i,i) = 0.;
          for (size_t j = 0; j < i; ++j)
          {
            (*mat)(i,j) = (*mat)(j,i) = dist->getDistanceForPair((*mapping1)[i], (*mapping1)[j]);
          }
        }

        string matrixFile = ApplicationTools::getAFilePath("clustering.output.matrix.file", comap.getParams(), false, false, "", false, "none", 1);
        if (matrixFile != "none")
        {
          //Write out matrix to a file, in Phylip format:
          PhylipDistanceMatrixFormat phylip;
          phylip.writeDistanceMatrix(*mat, matrixFile, true);
          ApplicationTools::displayResult("Wrote matrix to file", matrixFile);
        }
  
        // Clustering
  
        // The leaves in the clustering tree are the position in the mapping.
        // We will translate them before writing it to file.
        vector<string> matNames(nbSites1);
        for (size_t i = 0; i < nbSites1; ++i)
          matNames[i] = TextTools::toString(i);
        mat->setNames(matNames);
   
        unique_ptr<AgglomerativeDistanceMethod> clustering;
        if (clusteringMethod == "complete")
        {
          clustering.reset(new HierarchicalClustering(HierarchicalClustering::COMPLETE, *mat, false));
        }
        else if (clusteringMethod == "single")
        {
          clustering.reset(new HierarchicalClustering(HierarchicalClustering::SINGLE, *mat, false));
        }
        else if (clusteringMethod == "average")
        {
          clustering.reset(new HierarchicalClustering(HierarchicalClustering::AVERAGE, *mat, false));
        }
        //else if(clusteringMethod == "sum")
        //{
        //  clustering = new SumClustering(*mapping1, *dist, *mat);
        //}
        else
        {
          ApplicationTools::displayError("Unknown clustering method.");
          exit(-1);
        }
        ApplicationTools::displayResult("Clustering method", clusteringMethod);
  
        // Build tree:
        TreeTemplate<Node> clusteringTree(* clustering->getTree());

        // Add information to tree:
        ClusterTools::computeNormProperties(clusteringTree, *mapping1);
        dist->setStatisticAsProperty(*clusteringTree.getRootNode(), *mapping1);
  
        // Dumping groups to file, with more or less information, depending on the method used.
        string groupsPath = ApplicationTools::getAFilePath("clustering.output.groups.file", comap.getParams(), true, false, "", false, "groups_output_stats.txt", 0);
        ApplicationTools::displayResult("Site clusters output file", groupsPath);
        vector<string> colNames(6);
        colNames[0] = "Group";
        colNames[1] = "Size";
        colNames[2] = "IsConstant";
        colNames[3] = "Dmax";
        colNames[4] = "Stat";
        colNames[5] = "Nmin";
        //colNames[6] = "Delta"; //Distance From Mean Vector
        DataTable groupsData(colNames);

        // A few infos we will need:
        vector<double> rates = tl1->getPosteriorRateOfEachSite();
        vector<bool> isConst(nbSites1);
        for (size_t i = 0; i < nbSites1; ++i)
          isConst[i] = SiteTools::isConstant(sites1->getSite(i), true);
  
        vector<Group> groups = ClusterTools::getGroups(&clusteringTree);
        //vector<double> minNorms(groups.size());
  
        unsigned int maxGroupSize = ApplicationTools::getParameter<unsigned int>("clustering.maximum_group_size", comap.getParams(), 10, "", true, 2);
        for (size_t i = 0; i < groups.size(); ++i)
        {
          Group* group = &groups[i];
          if (group->size() > maxGroupSize) continue; 

          // Does the group contain a constant site?
          bool test = false;
          for (size_t j = 0; j < group->size() && !test; ++j)
          {
            if (isConst[TextTools::to<unsigned int>((*group)[j])]) test = true;
          }

          //// Compute distance from mean vector:
          //vector<double> groupMeanVector(nbBranches, 0.);
          //for(unsigned int j = 0; j < group->size(); j++)
          //{
          //  groupMeanVector += (*mapping1)[TextTools::to<unsigned int>((*group)[j])]; 
          //}
          //double distFromMeanVector = dist->getDistanceForPair(groupMeanVector/group->size(), meanVector);   

          //?minNorms[i] = minNorm;
          // Store results:
          vector<string> row(6);
          row[0] = group->toString(siteNames);
          row[1] = TextTools::toString(group->size());
          row[2] = test ? "yes" : "no";
          row[3] = TextTools::toString(group->getHeight() * 2.);
          row[4] = TextTools::toString(dynamic_cast<const Number<double> *>(group->getProperty("Stat"))->getValue());
          row[5] = TextTools::toString(dynamic_cast<const Number<double> *>(group->getProperty("Nmin"))->getValue());
          //row[6] = TextTools::toString(distFromMeanVector);
          groupsData.addRow(row);
        }

        //Write detected groups to file:
        ofstream groupsFile(groupsPath.c_str(), ios::out);
        DataTable::write(groupsData, groupsFile);
        groupsFile.close();

        //Write clustering tree to file:
        string treeFile = ApplicationTools::getAFilePath("clustering.output.tree.file", comap.getParams(), false, false, "", false, "none", 1);
        if (treeFile != "none")
        {
          // First we retranslate leaves names:
          ClusterTools::translate(clusteringTree, siteNames);
          Newick newick;
          newick.writeTree(clusteringTree, treeFile, true);
          ApplicationTools::displayResult("Wrote tree to file", treeFile);
        }
    
        // Now test each group:
        ApplicationTools::displayMessage("\n\n-*- Compute null distribution of clusters -*-\n");
  
        bool testGroupsGlobal = ApplicationTools::getBooleanParameter("clustering.null", comap.getParams(), false, "", true, 0);
        ApplicationTools::displayResult("Maximum group size to test", TextTools::toString(maxGroupSize));
        if (testGroupsGlobal)
        {       
          string simPath = ApplicationTools::getAFilePath("clustering.null.output.file", comap.getParams(), true, false, "", false, "groups_output_null.txt", 0);
          unsigned int nrep = ApplicationTools::getParameter<unsigned int>("clustering.null.number", comap.getParams(), 1, "", true, 1);
          ApplicationTools::displayResult("Number of simulations", TextTools::toString(nrep));
          ApplicationTools::displayResult("Simulations output file", simPath);
          ofstream* out = 0;
          if (simPath != "none") out = new ofstream(simPath.c_str(), ios::out);
          ApplicationTools::displayTask("Simulating groups", true);
          ClusterTools::computeGlobalDistanceDistribution(*tl1, *seqSim1, *substitutionCount1, *dist, *clustering, nbSites1, nrep, maxGroupSize, out);
          ApplicationTools::displayTaskDone();
          if (out) out->close();
        }
      }
    }

  



    // ****************************
    // * Candidate sites analysis *
    // ****************************
  
    else if (analysis == "candidates")
    {
      // *****************************
      // * Candidate groups analysis *
      // *****************************
      // We only deal with the one data set case for now.
    
      shared_ptr<const Statistic> statistic(CoETools::getStatistic(comap.getParams(), alphabet1.get(), substitutionCount1.get()));
  
      string groupsPath = ApplicationTools::getAFilePath("candidates.input.file", comap.getParams(), false, true);
      if (groupsPath != "none")
      {
        ApplicationTools::displayResult("Candidate groups are in file", groupsPath);

        ifstream groupsFile(groupsPath.c_str(), ios::in);

        //Store all positions and norm ranges in site container for each site in each candidate group:
        unsigned int minSim = ApplicationTools::getParameter<unsigned int>("candidates.null.min", comap.getParams(), 1000, "", true, true);
        ApplicationTools::displayResult("Minimum number of simulations", TextTools::toString(minSim));
        unsigned int verbose = ApplicationTools::getParameter<unsigned int>("candidates.null.verbose", comap.getParams(), 1, "", true, true);
        ApplicationTools::displayResult("Verbose level", TextTools::toString(verbose));
        CandidateGroupSet candidates(statistic.get(), minSim, verbose);
    
        //Create positions index first:
        map<int, size_t> posIndex;
        vector<int> positions = sites1->getSitePositions();
        for (size_t i = 0; i < positions.size(); i++)
          posIndex[positions[i]] = i;

        //Parse file:
        string groupsColumnSep = ApplicationTools::getStringParameter("candidates.input.column_sep", comap.getParams(), "\t", "", true, true);
        shared_ptr<DataTable> table(DataTable::read(groupsFile, groupsColumnSep, true, -1));
        groupsFile.close();
        string groupsColumnName = ApplicationTools::getStringParameter("candidates.input.column_name", comap.getParams(), "Group", "", true, true);
        vector<string> groups = table->getColumn(groupsColumnName);

        //Get range value:
        double omega = ApplicationTools::getDoubleParameter("candidates.omega", comap.getParams(), 0.25, "", true, true);
        ApplicationTools::displayResult("Norm interval", TextTools::toString(omega));
        if (omega < 0)
        {
          ApplicationTools::displayWarning("Norm range parameter 'omega' must be positive... |omega| was used instead.");
          omega = -omega;
        }

        //Now parse each group, compute index and corresponding norm ranges and statistic value:
        ApplicationTools::displayTask("Compute statistic for each group", true);
        string group, tmp, strok = "0123456789;,";
        for (size_t i = 0; i < groups.size(); i++)
        {
          if (groups.size() > 100) ApplicationTools::displayGauge(i, groups.size() - 1, '=');
          tmp = TextTools::removeWhiteSpaces(groups[i]);
          group = "";
          //Clean group description:
          for (unsigned int j = 0; j < tmp.length(); j++)
          {
            if (strok.find(tmp[j]) != string::npos) group += tmp[j];
          }
          //Now parse group for each site:
          StringTokenizer st(group, ";,", false);
          int pos;
          if (st.numberOfRemainingTokens() <= 1)
          {
            throw Exception("Error, group " + TextTools::toString(i) + " has " + TextTools::toString(st.numberOfRemainingTokens()) + "sites.");
          }
          CandidateGroup candidate;
          while (st.hasMoreToken())
          {
            pos = TextTools::toInt(st.nextToken());
            if (posIndex.find(pos) == posIndex.end())
            {
              candidate.setAnalysable(false);
              ApplicationTools::displayWarning("Position " + TextTools::toString(pos) + " is not included in the selected sites. The group will be ignored (line " + TextTools::toString(i+1) + " in input file).");
              break;
            }
            else
            {
              candidate.addSite(CandidateSite(posIndex[pos]));
            }
          }
          if (candidate.isAnalysable())
          {
            candidate.computeStatisticValue(*statistic, *mapping1);
            candidate.computeNormRanges(omega, *mapping1);
          }
          candidates.addCandidate(candidate);
        }
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Number of groups to test", candidates.size());

        if (candidates.size() == 0)
          throw Exception("ERROR!!! No group can be tested!");

        //Now compute p-values:
        unsigned int nbMaxTrials = ApplicationTools::getParameter<unsigned int>("candidates.nb_max_trials", comap.getParams(), 10, "", true, true);
        CoETools::computePValuesForCandidateGroups(candidates, *tl1, *seqSim1, *substitutionCount1, comap.getParams(), nbMaxTrials);

        //Now fill data table:
        vector<string> stats(table->getNumberOfRows()), pvalues(table->getNumberOfRows());
        for (unsigned int i = 0; i < table->getNumberOfRows(); i++)
        {
          stats[i] = candidates[i].isAnalysable()
            ? TextTools::toString(candidates[i].getStatisticValue(), 6)
            : "NA";
          pvalues[i] = candidates[i].isAnalysable()
            ? TextTools::toString(candidates.getPValueForGroup(i), 6)
            : "NA";
        }
        table->addColumn("Stat", stats);
        table->addColumn("p-value", pvalues);

        //Write result to file:
        string outputPath = ApplicationTools::getAFilePath("candidates.output.file", comap.getParams(), true, false);
        groupsColumnSep = ApplicationTools::getStringParameter("candidates.output.column_sep", comap.getParams(), groupsColumnSep, "", true, true);
        ofstream outputFile(outputPath.c_str(), ios::out);
        DataTable::write(*table, outputFile, groupsColumnSep);
        outputFile.close();
        ApplicationTools::displayResult("Wrote results in file", outputPath);
      }
    }






    else throw Exception("Unknown analysis type: " + analysis);
  }

  // ***********
  // * The End *
  // ***********
  
  comap.done();
  
  ApplicationTools::displayMessage("Bye bye ;-)");  

  }
  catch(Exception e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

