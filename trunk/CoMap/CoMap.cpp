//
// File: CoMap.cpp
// Created by: Julien Dutheil
// Created on: Thu Apr 6 11:06 2006
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

// From the STL:

#include <iostream>
#include <cstdlib>
#include <stdexcept>

using namespace std;

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

// From NumCalc:
#include <NumCalc/DataTable.h>

// From SeqLib:
#include <Seq/SiteTools.h>

// From PhylLib:
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/HomogeneousSequenceSimulator.h>
#include <Phyl/NonHomogeneousSequenceSimulator.h>
#include <Phyl/MutationProcess.h>
#include <Phyl/RASTools.h>
#include <Phyl/PhylipDistanceMatrixFormat.h>
#include <Phyl/AgglomerativeDistanceMethod.h>
#include <Phyl/Newick.h>

using namespace bpp;

// CoMap's include files:
#include "CoETools.h"
#include "Distance.h"
#include "Cluster.h"
#include "ClusterTools.h"

/******************************************************************************/

void help()
{
  cout << "Parameters reminder. For a complete description, see the help files." << endl;
 	cout << "____________________________________________________________________" << endl;
	cout << "param           | parameters file." << endl;
	cout << "alphabet        | the datatype (RNA, DNA or Protein)." << endl;
  cout << "tree_file       | the tree file." << endl;
  cout << "sequence_file   | the sequence file." << endl;
  cout << "output.vectors  | the vectors output file." << endl;
  cout << "output.tree     | the tree with estimated branches lengths output file." << endl;
	cout << endl;
	PhylogeneticsApplicationTools::printOptimizationHelp(false, false);
	cout << endl;
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
	cout << "* This is CoMap        version 1.4.0       date: 14/02/09 *" << endl;
	cout << "*     A C++ shell program to detect co-evolving sites.    *" << endl;
	cout << "***********************************************************" << endl;
	cout << endl;
  
	try
  {
  
  ApplicationTools::startTimer();

	// **************************
	// * Retrieving parameters: *
	// **************************
	
	if(argc == 1)
  { 
    // No argument, show display some help and leave.
		help();
		exit(-1);
	}

  cout << "Parsing options:" << endl;
  
  map<string, string> params = AttributesTools::parseOptions(argc, argv);

	ApplicationTools::displayMessage("\n\n-*- Retrieve data and model -*-\n");

	TreeTemplate<Node> * tree1 = PhylogeneticsApplicationTools::getTree(params);
	ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree1->getNumberOfLeaves()));
	ApplicationTools::displayResult("Number of sons at root", TextTools::toString(tree1->getRootNode()->getNumberOfSons()));
	
	// Get data 1:
	//Tree * tree1 = new Tree(* tree); // Copy tree.
	Alphabet *alphabet1;
	VectorSiteContainer *allSites1, *sites1;
	SubstitutionModel *model1;
	SubstitutionModelSet *modelSet1;
	DiscreteDistribution *rDist1;
	DRTreeLikelihood *tl1;
	CoETools::readData(tree1, alphabet1, allSites1, sites1, model1, modelSet1, rDist1, tl1, params, "");
	ApplicationTools::displayResult("Number of sites in file", allSites1->getNumberOfSites());
	ApplicationTools::displayResult("Number of sites to analyse", sites1->getNumberOfSites());
	ApplicationTools::displayResult("Number of site patterns", tl1->getLikelihoodData()->getNumberOfDistinctSites());
  
  bool continuousSim = ApplicationTools::getBooleanParameter("simulations.continuous", params, false, "", true, false);
	ApplicationTools::displayResult("Rate distribution for simulations", (continuousSim ? "continuous" : "discrete"));

	ApplicationTools::displayMessage("\n\n-*- Get substitution vectors -*-\n");

	// Building a simulator object:
  SequenceSimulator* seqSim1 = NULL;
  if(modelSet1)
  {
    seqSim1 = new NonHomogeneousSequenceSimulator(modelSet1, rDist1, tree1);
    dynamic_cast<NonHomogeneousSequenceSimulator *>(seqSim1)->enableContinuousRates(continuousSim);
  }
  else
  {
    seqSim1 = new HomogeneousSequenceSimulator(model1, rDist1, tree1);
    dynamic_cast<HomogeneousSequenceSimulator *>(seqSim1)->enableContinuousRates(continuousSim);
  }

	// Getting the substitutions count function:
	SubstitutionCount* substitutionCount1 = CoETools::getSubstitutionCount(alphabet1, model1, rDist1, params);
		
	// Getting the substitution vectors:
	ProbabilisticSubstitutionMapping* mapping1 = CoETools::getVectors(*tl1, *substitutionCount1, *sites1, params);

  string analysis = ApplicationTools::getStringParameter("analysis", params, "pairwise");
  ApplicationTools::displayResult("Analysis type", analysis);

  if(analysis == "pairwise")
  {
	  const Statistic * statistic = CoETools::getStatistic(params);

    bool computeNullHyp = false;
	  computeNullHyp = ApplicationTools::getBooleanParameter("statistic.null", params, false);
	
  
    // *******************************************
	  // * The same for a putative second dataset: *
	  // *******************************************
	  if(params.find("sequence.file2") != params.end() && params["sequence.file2"] != "none")
    {
		  TreeTemplate<Node> * tree2;
		  if(params.find("tree.file2") != params.end() && params["tree.file2"] != "none")
      {
		    ApplicationTools::displayMessage("WARNING!!! Second tree file specified.\n Tree 1 and Tree 2 must differ only by their branch lengths, otherwise results may be uninterpretable.\n");			
  	    tree2 = PhylogeneticsApplicationTools::getTree(params, "2", true);
	      ApplicationTools::displayResult("# number of leaves", TextTools::toString(tree2 -> getNumberOfLeaves()));
    	  ApplicationTools::displayResult("# number of sons at root", TextTools::toString(tree2 -> getRootNode() -> getNumberOfSons()));
		  }
      else
      {
			  tree2 = new TreeTemplate<Node>(* tree1); // Copy tree.
		  }

		  Alphabet *alphabet2;
		  VectorSiteContainer *allSites2, *sites2;
		  SubstitutionModel *model2;
		  SubstitutionModelSet *modelSet2;
		  DiscreteDistribution *rDist2;
		  DRTreeLikelihood *tl2;
		  ApplicationTools::displayMessage("\nLoading second dataset...\n");
		  CoETools::readData(tree2, alphabet2, allSites2, sites2, model2, modelSet2, rDist2, tl2, params, "2");
	    ApplicationTools::displayResult("Number of sites in file", allSites2->getNumberOfSites());
	    ApplicationTools::displayResult("Number of sites to analyse", sites2->getNumberOfSites());
	    ApplicationTools::displayResult("Number of site patterns", tl2->getLikelihoodData()->getNumberOfDistinctSites());

		  ApplicationTools::displayMessage("\n... and get its substitution vectors.\n");
		
		  // Building a simulator object:
      SequenceSimulator* seqSim2 = NULL;
      if(modelSet2)
      {
        seqSim2 = new NonHomogeneousSequenceSimulator(modelSet2, rDist2, tree2);
        dynamic_cast<NonHomogeneousSequenceSimulator *>(seqSim2)->enableContinuousRates(continuousSim);
      }
      else
      {
        seqSim2 = new HomogeneousSequenceSimulator(model2, rDist2, tree2);
        dynamic_cast<HomogeneousSequenceSimulator *>(seqSim2)->enableContinuousRates(continuousSim);
      }

		  // Getting the substitutions count function:
		  SubstitutionCount* substitutionCount2 = CoETools::getSubstitutionCount(alphabet2, model2, rDist2, params);
		
		  // Getting the substitution vectors:
		  ProbabilisticSubstitutionMapping * mapping2 = CoETools::getVectors(*tl2, *substitutionCount2, *sites2, params, "2");
		
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
			  params);
			
		  CoETools::writeInfos(*sites2, *tl2, params, "2");
		  
      //tl2 will be modified after the simulations!
      if(computeNullHyp)
      {
        CoETools::computeInterNullDistribution(
          *tl1,
          *tl2,
			    *seqSim1,
			    *seqSim2,
			    *substitutionCount1,
          *substitutionCount2,
			    *statistic,
			    params);
      }

      delete tree2;
		  delete allSites2;
		  delete sites2;
		  if(model2) delete model2;
		  if(modelSet2) delete modelSet2;
		  delete rDist2;
		  delete tl2;
		  delete substitutionCount2;
      delete seqSim2;
		  ApplicationTools::displayTaskDone();
	  }
    else
    {
			ApplicationTools::displayMessage("\n\n-*- Compute statistics -*- \n");
		
			ApplicationTools::displayMessage("Analyse dataset.");
			CoETools::computeIntraStats(
				*tl1,
				*sites1,
				*mapping1,
				*statistic,
				params);

		  CoETools::writeInfos(*sites1, *tl1, params);
      
      //tl2 will be modified after the simulations!
      if(computeNullHyp)
      {
			  CoETools::computeIntraNullDistribution(
          *tl1,
				  *seqSim1,
				  *substitutionCount1,
				  *statistic,
				  params);
		  }
	  }

    delete statistic;
  }






  // ***********************
	// * Clustering analysis *
	// ***********************
	
  else if(analysis == "clustering")
  {
		CoETools::writeInfos(*sites1, *tl1, params);

    ApplicationTools::displayMessage("\n\n-*- Perform clustering analysis -*-\n");
	
    string clusteringMethod = ApplicationTools::getStringParameter("clustering.method", params, "none", "", false, false);
	  if(clusteringMethod != "none")
    {
      unsigned int nbBranches = mapping1->getNumberOfBranches();
      unsigned int nbSites    = mapping1->getNumberOfSites();
      bool scale = ApplicationTools::getBooleanParameter("clustering.scale", params, false, "", true, true);
      vector<double> scales(nbBranches, 1.);
      vector<double> weights(nbBranches, 1.);
	    vector<double> brLens = mapping1->getBranchLengths();
  	  double minLen = ApplicationTools::getDoubleParameter("clustering.min_length", params, false, "", 0.000001, false);
      vector<double> meanVector(nbBranches);
      for(unsigned int j = 0; j < nbBranches; j++)
      {
	  	  double sum = 0;
        for(unsigned int i = 0; i < nbSites; i++)
        {
          sum += (*mapping1)(j, i);
        }
        meanVector[j] = sum / nbSites;
      }
  
      // Scale if required (i.e., devide each vector by the mean vector):
	    for(unsigned int j = 0; j < nbBranches; j++)
      {
        double len = brLens[j];
        if(len <= minLen) weights[j] = 0.; // Branch ignored, considered as a multifurcation.
        if(scale)
        {
          double scale = 0;
          if(len > minLen)
          {
            double m = meanVector[j];
            if(m > 0) scale = 1./m; // if mean==0, hence the susbstitution number is 0 too, the scale does not matter... 0/0=>0.
            //else scale = 0.;
          }// else: weight[j] == 0, so the position will not be used in distance calculation whatever...
          for(unsigned int i = 0; i < nbSites; i++)
          {
            (*mapping1)(j, i) *= scale ;
          }
          scales[j] = scale;
        }
	    }  
	    ApplicationTools::displayResult("Scale by row", scale ? "yes" : "no");
  
      // Compute norms:
      vector<double> norms(nbSites);
      vector<string> siteNames(nbSites);
      for(unsigned int i = 0; i < nbSites; i++)
      {
        string siteName = TextTools::toString(sites1->getSite(i)->getPosition());
        siteNames[i] = siteName; 
        norms[i] = VectorTools::norm<double, double>((*mapping1)[i]);
      }
  
	    // Get the distance to use
	
	    string distanceMethod = ApplicationTools::getStringParameter("clustering.distance", params, "cor", "", true, true);
	    Distance * dist = NULL;
	    if(distanceMethod == "euclidian")
      {
		    dist = new EuclidianDistance();
	    }
      else if(distanceMethod == "cor")
      {
        Statistic * cor = new CorrelationStatistic();
		    dist = new StatisticBasedDistance(cor, 1.);
      }
      else if(distanceMethod == "comp")
      {
        string nijtOption = ApplicationTools::getStringParameter("nijt", params, "simule", "", true);
        bool sym = ApplicationTools::getBooleanParameter("nijt_aadist.sym", params, true, "", true); 
        if(nijtOption != "aadist" || sym)
        {
          throw Exception("Compensation distance must be used with 'nijt=aadist' and 'nijt_aadist.sym=no' options.");
        }
        else
        {
		      dist = new CompensationDistance();
        }
      }
      else
      {
		    ApplicationTools::displayError("Unknown distance method.");
		    exit(-1);
	    }
      dist->setWeights(weights);
	    ApplicationTools::displayResult("Distance to use", distanceMethod);

	    // Compute the distance matrix.
	    DistanceMatrix * mat = new DistanceMatrix(siteNames);
	    for(unsigned int i = 0; i < nbSites; i++)
      {
		    (*mat)(i,i) = 0.;
		    for(unsigned int j = 0; j < i; j++)
        {
		  	  (*mat)(i,j) = (*mat)(j,i) = dist->getDistanceForPair((*mapping1)[i], (*mapping1)[j]);
		    }
	    }

	    string matrixFile = ApplicationTools::getAFilePath("clustering.output.matrix.file", params, false, false, "", false);
	    if(matrixFile != "none")
      {
		    //Write out matrix to a file, in Phylip format:
		    PhylipDistanceMatrixFormat phylip;
		    phylip.write(*mat, matrixFile, true);
		    ApplicationTools::displayResult("Wrote matrix to file", matrixFile);
	    }
  
  	  // Clustering
	
      // The leaves in the clustering tree are the position in the mapping.
      // We will translate them before writing it to file.
      vector<string> matNames(nbSites);
      for(unsigned int i = 0; i < nbSites; i++) matNames[i] = TextTools::toString(i);
      mat->setNames(matNames);
	 
	    AgglomerativeDistanceMethod * clustering = NULL;
	    if(clusteringMethod == "complete")
      {
		    clustering = new SimpleClustering(SimpleClustering::COMPLETE, *mat, false);
	    }
      else if(clusteringMethod == "single")
      {
		    clustering = new SimpleClustering(SimpleClustering::SINGLE, *mat, false);
	    }
      else if(clusteringMethod == "average")
      {
		    clustering = new SimpleClustering(SimpleClustering::AVERAGE, *mat, false);
	    }
      //else if(clusteringMethod == "sum")
      //{
  	  //	clustering = new SumClustering(*mapping1, *dist, *mat);
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
	    string groupsPath = ApplicationTools::getAFilePath("clustering.output.groups.file", params, true, false, "", false);
      ApplicationTools::displayResult("Site clusters output file", groupsPath);
      vector<string> colNames(6);
      colNames[0] = "Group";
      colNames[1] = "Size";
      colNames[2] = "Const";
      colNames[3] = "Dmax";
      colNames[4] = "Stat";
      colNames[5] = "Nmin";
      //colNames[6] = "Delta"; //Distance From Mean Vector
      DataTable groupsData(colNames);

      // A few infos we will need:
	    vector<double> rates = tl1->getPosteriorRateOfEachSite();
      vector<bool> isConst(nbSites);
      for(unsigned int i = 0; i < nbSites; i++)
        isConst[i] = SiteTools::isConstant(*(sites1->getSite(i)), true);
  
      vector<Group> groups = ClusterTools::getGroups(&clusteringTree);
      //vector<double> minNorms(groups.size());
  
      for(unsigned int i = 0; i < groups.size(); i++)
      {
        Group * group = &groups[i];
      
        // Does the group contain a constant site?
        bool test = false;
        for(unsigned int j = 0; j < group->size() && !test; j++)
        {
          if(isConst[TextTools::to<unsigned int>((*group)[j])]) test = true;
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
	    string treeFile = ApplicationTools::getAFilePath("clustering.output.tree.file", params, false, false, "", false);
	    if(treeFile != "none")
      {
	      // First we retranslate leaves names:
        ClusterTools::translate(clusteringTree, siteNames);
        Newick newick;
	      newick.write(clusteringTree, treeFile, true);
	      ApplicationTools::displayResult("Wrote tree to file", treeFile);
      }
    
      // Now test each group:
      ApplicationTools::displayMessage("\n\n-*- Compute null distribution of clusters -*-\n");
  
      bool testGroupsGlobal = ApplicationTools::getBooleanParameter("clustering.null", params, false, "", true, false);
      unsigned int maxGroupSize = ApplicationTools::getParameter<unsigned int>("clustering.maximum_group_size", params, 10, "", true, false);
      ApplicationTools::displayResult("Maximum group size to test", TextTools::toString(maxGroupSize));
      if(testGroupsGlobal)
      {       
	      string simPath = ApplicationTools::getAFilePath("clustering.null.output.file", params, true, false, "", false);
        unsigned int nrep = ApplicationTools::getParameter<unsigned int>("clustering.null.number", params, 1, "", true);
        ApplicationTools::displayResult("Number of simulations", TextTools::toString(nrep));
        ApplicationTools::displayResult("Simulations output file", simPath);
        ofstream * out = NULL;
        if(simPath != "none") out = new ofstream(simPath.c_str(), ios::out);
        ApplicationTools::displayTask("Simulating groups", true);
        ClusterTools::computeGlobalDistanceDistribution(*tl1, *seqSim1, *substitutionCount1, *dist, *clustering, scales, nbSites, nrep, maxGroupSize, out);
        ApplicationTools::displayTaskDone();
        if(out != NULL) out->close();
      }

      // Clustering memory freeing:
	    delete dist;
	    delete mat;
	    delete clustering;
    }
  }

  



  // ****************************
	// * Candidate sites analysis *
	// ****************************
	
  else if(analysis == "candidates")
  {
    // *****************************
  	// * Candidate groups analysis *
  	// *****************************
    // We only deal with the one data set case for now.
	  
		CoETools::writeInfos(*sites1, *tl1, params);
    
    const Statistic * statistic = CoETools::getStatistic(params);
	
    string groupsPath = ApplicationTools::getAFilePath("candidates.input.file", params, false, true);
    if(groupsPath != "none")
    {
      ApplicationTools::displayResult("Candidate groups are in file", groupsPath);

      ifstream groupsFile(groupsPath.c_str(), ios::in);

      //Store all positions and norm ranges in site container for each site in each candidate group:
      unsigned int minSim = ApplicationTools::getParameter<unsigned int>("candidates.null.min", params, 1000, "", true, true);
      ApplicationTools::displayResult("Minimum number of simulations", TextTools::toString(minSim));
      unsigned int verbose = ApplicationTools::getParameter<unsigned int>("candidates.null.verbose", params, 1, "", true, true);
      ApplicationTools::displayResult("Verbose level", TextTools::toString(verbose));
      CandidateGroupSet candidates(statistic, minSim, verbose);
    
      //Create positions index first:
      map<int, unsigned int> posIndex;
      vector<int> positions = sites1->getSitePositions();
      for(unsigned int i = 0; i < positions.size(); i++)
        posIndex[positions[i]] = i;

      //Parse file:
      string groupsColumnSep = ApplicationTools::getStringParameter("candidates.input.column_sep", params, "\t", "", true, true);
      DataTable* table = DataTable::read(groupsFile, groupsColumnSep, true, -1);
      groupsFile.close();
      string groupsColumnName = ApplicationTools::getStringParameter("candidates.input.column_name", params, "Group", "", true, true);
      vector<string> groups = table->getColumn(groupsColumnName);

      //Get range value:
      double omega = ApplicationTools::getDoubleParameter("candidates.omega", params, 0.25, "", true, true);
      ApplicationTools::displayResult("Norm interval", TextTools::toString(omega));
      if(omega < 0)
      {
        ApplicationTools::displayWarning("Norm range parameter 'omega' must be positive... |omega| was used instead.");
        omega = -omega;
      }

      //Now parse each group, compute index and corresponding norm ranges and statistic value:
      ApplicationTools::displayTask("Compute statistic for each group", true);
      string group, tmp, strok = "0123456789;,";
      for(unsigned int i = 0; i < groups.size(); i++)
      {
        if(groups.size() >100) ApplicationTools::displayGauge(i, groups.size() - 1, '=');
        tmp = TextTools::removeWhiteSpaces(groups[i]);
        group = "";
        //Clean group description:
        for(unsigned int j = 0; j < tmp.length(); j++)
        {
          if(strok.find(tmp[j]) != string::npos) group += tmp[j];
        }
        //Now parse group for each site:
        StringTokenizer st(group, ";,", false);
        int pos;
        if(st.numberOfRemainingTokens() <= 1)
        {
          throw Exception("Error, group " + TextTools::toString(i) + " has " + TextTools::toString(st.numberOfRemainingTokens()) + "sites.");
        }
        CandidateGroup candidate;
        while(st.hasMoreToken())
        {
          pos = TextTools::toInt(st.nextToken());
          if(posIndex.find(pos) == posIndex.end())
          {
            candidate.setAnalysable(false);
            ApplicationTools::displayWarning("Position " + TextTools::toString(pos) + " is not included in the selected sites. The group will be ignored (line " + TextTools::toString(i+1) + " in input file).");
            break;
          }
          candidate.push_back(CandidateSite(posIndex[pos]));
        }
        if(candidate.isAnalysable())
        {
          candidate.computeStatisticValue(*statistic, *mapping1);
          candidate.computeNormRanges(omega, *mapping1);
        }
        candidates.addCandidate(candidate);
      }
      ApplicationTools::displayTaskDone();

      //Now compute p-values:
      CoETools::computePValuesForCandidateGroups(candidates, *tl1, *seqSim1, *substitutionCount1, params);

      //Now fill data table:
      vector<string> stats(table->getNumberOfRows()), pvalues(table->getNumberOfRows());
      for(unsigned int i = 0; i < table->getNumberOfRows(); i++)
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
      string outputPath = ApplicationTools::getAFilePath("candidates.output.file", params, true, false);
      groupsColumnSep = ApplicationTools::getStringParameter("candidates.output.column_sep", params, groupsColumnSep, "", true, true);
      ofstream outputFile(outputPath.c_str(), ios::out);
      DataTable::write(*table, outputFile, groupsColumnSep);
      outputFile.close();
      ApplicationTools::displayResult("Wrote results in file", outputPath);

      //Free memory.
      delete table;
    }

    delete statistic;
  }






  else throw Exception("Unknown analysis type: " + analysis);


  // ***********
	// * The End *
	// ***********
	
	ApplicationTools::displayTask("CoMap's done. Freeing memory");
  
  // Dataset 1 memory freeing:
  delete tree1;
	delete allSites1;
	delete sites1;
	if(model1) delete model1;
	if(modelSet1) delete modelSet1;
	delete rDist1;
	delete tl1;
	delete substitutionCount1;
  delete mapping1;
  delete seqSim1;
	
  ApplicationTools::displayTaskDone();

  ApplicationTools::displayTime("Total execution time:");
	
  ApplicationTools::displayMessage("Bye bye ;-)");	

	}
  catch(Exception e)
  {
		cout << e.what() << endl;
		exit(-1);
	}

	return (0);
}

