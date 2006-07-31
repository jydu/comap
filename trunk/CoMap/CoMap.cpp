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
using namespace std;
#include <cstdlib>
#include <stdexcept>

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
#include <Phyl/HomogeneousTreeLikelihood.h>
#include <Phyl/MutationProcess.h>
#include <Phyl/RASTools.h>
#include <Phyl/PhylipDistanceMatrixFormat.h>
#include <Phyl/AgglomerativeDistanceMethod.h>
#include <Phyl/Newick.h>

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
	PhylogeneticsApplicationTools::printSubstitutionModelHelp();
	cout << endl;
	PhylogeneticsApplicationTools::printRateDistributionHelp();
	cout << endl;
	PhylogeneticsApplicationTools::printOptimizationHelp();
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
	cout << "* This is CoMap       version 1.0.0      date: 25/07/2006 *" << endl;
	cout << "*     A C++ shell program to detect co-evolving sites.    *" << endl;
	cout << "***********************************************************" << endl;
	cout << endl;

	try
  {

	// **************************
	// * Retrieving parameters: *
	// **************************
	
	if(argc == 1)
  { // No argument, show display some help and leave.
		help();
		exit(-1);
	}

	// Get the parameters from command line:
	map<string, string> cmdParams = AttributesTools::getAttributesMap(
	AttributesTools::getVector(argc, argv), "=");

	// Look for a specified file with parameters:
	map<string, string> params;
	if(cmdParams.find("param") != cmdParams.end())
  {
		string file = cmdParams["param"];
		if(!FileTools::fileExists(file))
    {
			ApplicationTools::displayError("Parameter file not found.");
			exit(-1);
		}
    else
    {
			params = AttributesTools::getAttributesMapFromFile(file, "=");
			// Actualize attributes with ones passed to command line:
			AttributesTools::actualizeAttributesMap(params, cmdParams);
		}
	}
  else
  {
		params = cmdParams;
	}

	ApplicationTools::displayMessage("\n\nI) Retrieve data and model\n");

	TreeTemplate<Node> * tree1 = PhylogeneticsApplicationTools::getTree(params);
	ApplicationTools::displayResult("# number of leaves", TextTools::toString(tree1 -> getNumberOfLeaves()));
	ApplicationTools::displayResult("# number of sons at root", TextTools::toString(tree1 -> getRootNode() -> getNumberOfSons()));
	
	// Get data 1:
	//Tree * tree1 = new Tree(* tree); // Copy tree.
	Alphabet *alphabet1;
	VectorSiteContainer *allSites1, *sites1;
	SubstitutionModel *model1;
	DiscreteDistribution *rDist1;
	HomogeneousTreeLikelihood *tl1;
	CoETools::readData(tree1, alphabet1, allSites1, sites1, model1, rDist1, tl1, params, "");

	ApplicationTools::displayMessage("\n\nII) Get substitution vectors\n");

	// Building a MutationProcess object:
	MutationProcess * process1 = new SimpleMutationProcess(model1);

	// Getting the substitutions count function:
	SubstitutionCount * substitutionCount1 = CoETools::getSubstitutionCount(
		alphabet1,
		* tree1,
		* process1,
		* rDist1,
		params);
		
	// Getting the substitution vectors:
	ProbabilisticSubstitutionMapping * mapping1 = CoETools::getVectors(
		alphabet1,
		* tree1,
		* sites1,
		* allSites1,
		* model1,
		* rDist1,
		* substitutionCount1,
		params);

	// Getting posterior rate class distribution:
	DiscreteDistribution * prDist1 = RASTools::getPosteriorRateDistribution(* tl1);

	ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
	prDist1-> print(ApplicationTools::message);
	ApplicationTools::displayMessage("\n");

	// *******************************************
	// * The same for a putative second dataset: *
	// *******************************************
	
	const Statistic * statistic = CoETools::getStatistic(params);
	bool computeNullHyp = ApplicationTools::getBooleanParameter("statistic.null", params, false);
	
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
		DiscreteDistribution *rDist2;
		HomogeneousTreeLikelihood *tl2;
		ApplicationTools::displayMessage("\nLoading second dataset...\n");
		CoETools::readData(tree2, alphabet2, allSites2, sites2, model2, rDist2, tl2, params, "2");

		ApplicationTools::displayMessage("\n... and get its substitution vectors.\n");
		
		// Building a MutationProcess object:
		MutationProcess * process2 = new SimpleMutationProcess(model2);

		// Getting the substitutions count function:
		SubstitutionCount * substitutionCount2 = CoETools::getSubstitutionCount(
			alphabet2,
			* tree2,
			* process2,
			* rDist2,
			params);
		
		// Getting the substitution vectors:
		ProbabilisticSubstitutionMapping * mapping2 = CoETools::getVectors(
			alphabet2,
			* tree2,
			* sites2,
			* allSites2,
			* model2,
			* rDist2,
			* substitutionCount2,
			params,
			"2");
		
		// Getting posterior rate class distribution:
		DiscreteDistribution * prDist2 = RASTools::getPosteriorRateDistribution(* tl2);

		ApplicationTools::displayMessage("\nPosterior rate distribution for 2nd dataset:\n");
		prDist2 -> print(ApplicationTools::message);
		ApplicationTools::displayMessage("\n");
		
		ApplicationTools::displayMessage("\n\nIII) Compute statistics\n");
		
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
			
		if(computeNullHyp) CoETools::computeInterNullDistribution(
			*process1, *process2,
			*rDist1, *rDist2,
			*tree1, *tree2,
			*substitutionCount1,
      *substitutionCount2,
			*statistic,
			params);

		CoETools::writeInfos(*sites2, *tl2, params, "2");

		delete tree2;
		delete allSites2;
		delete sites2;
		delete model2;
		delete rDist2;
		delete prDist2;
		delete tl2;
		delete substitutionCount2;
		delete process2;
    delete mapping2;
		ApplicationTools::displayTaskDone();
	}
  else
  {
		if(statistic != NULL)
    {
			ApplicationTools::displayMessage("\n\nIII) Compute statistics\n");
		
			ApplicationTools::displayMessage("Analyse dataset.");
			CoETools::computeIntraStats(
				*tl1,
				*sites1,
				*mapping1,
				*statistic,
				params);
		}

		if(computeNullHyp)
    {
			CoETools::computeIntraNullDistribution(
				*process1,
				*rDist1,
				*tree1,
				*substitutionCount1,
				*statistic,
				params);
		}
 
		CoETools::writeInfos(*sites1, *tl1, params);
	}
	ApplicationTools::displayTaskDone();








  
  // ***********************
	// * Clustering analysis *
	// ***********************

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
        if(len > minLen) {
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
	  ApplicationTools::displayResult("Scale by row:", scale ? "yes" : "no");
  
    // Compute norms:
    vector<double> norms(nbSites);
    vector<string> siteNames(nbSites);
    for(unsigned int i = 0; i < nbSites; i++)
    {
      string siteName = "Site" + TextTools::toString(sites1->getSite(i)->getPosition());
      siteNames[i] = siteName; 
      norms[i] = norm((*mapping1)[i]);
    }
  
	  // Get the distance to use
	
	  string distanceMethod = ApplicationTools::getStringParameter("clustering.distance", params, "cor", "", true, true);
	  Distance * dist = NULL;
	  if(distanceMethod == "euclidian")
    {
		  dist = new EuclidianDistance();
	  }
    else if(distanceMethod == "sqcor")
    {
		  dist = new SquareCorrelationDistance();
	  }
    else if(distanceMethod == "cor")
    {
		  dist = new CorrelationDistance();
    }
    else
    {
		  ApplicationTools::displayError("Unknown distance method.");
		  exit(-1);
	  }
    dist->setWeights(weights);
	  ApplicationTools::displayResult("Distance to use:", distanceMethod);

	  // Compute the distance matrix.
	  DistanceMatrix * mat = new DistanceMatrix(siteNames);
	  for(unsigned int i = 0; i < nbSites; i++)
    {
		  (*mat)(i,i) = 0.;
		  for(unsigned int j = 0; j < i; j++)
      {
		  	(*mat)(i,j) = (*mat)(j,i) = dist->d((*mapping1)[i], (*mapping1)[j]);
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
		  clustering = new SimpleClustering(SimpleClustering::COMPLETE, *mat);
	  }
    else if(clusteringMethod == "single")
    {
		  clustering = new SimpleClustering(SimpleClustering::SINGLE, *mat);
	  }
    else if(clusteringMethod == "average")
    {
		  clustering = new SimpleClustering(SimpleClustering::AVERAGE, *mat);
	  }
    else if(clusteringMethod == "sum")
    {
	  	clustering = new SumClustering(*mapping1, *dist, *mat);
	  }
    else
    {
		  ApplicationTools::displayError("Unknown clustering method.");
		  exit(-1);
	  }
	  ApplicationTools::displayResult("Clustering method:", clusteringMethod);
	
	  // Build tree:
	  Tree * clusteringTree = clustering->getTree();
	
    // Dumping groups to file, with more or less information, depending on the method used.
	  string groupsPath = ApplicationTools::getAFilePath("clustering.output.groups.file", params, true, false, "", false);
    vector<string> colNames(6);
    colNames[0] = "Group";
    colNames[1] = "Size";
    colNames[2] = "Const";
    colNames[3] = "Dmax";
    colNames[4] = "Nmin";
    colNames[5] = "Delta"; //Distance From Mean Vector
    DataTable groupsData(colNames);

    // A few infos we will need:
	  vector<double> rates = tl1->getPosteriorRateOfEachSite();
    vector<bool> isConst(nbSites);
    for(unsigned int i = 0; i < nbSites; i++)
      isConst[i] = SiteTools::isConstant(*(sites1->getSite(i)), true);
  
    vector<Group> groups = ClusterTools::getGroups(dynamic_cast<TreeTemplate<Node> *>(clusteringTree));
    vector<double> minNorms(groups.size());
  
    for(unsigned int i = 0; i < groups.size(); i++)
    {
      Group * group = &groups[i];
      // Compute minimal norms:
      double minNorm = -log(0.);
      for(unsigned int j = 0; j < group->size(); j++)
      {
        double norm = norms[TextTools::to<unsigned int>((*group)[j])];
        if(norm < minNorm) minNorm = norm; 
      }
      // Does the group contain a constant site?
      bool test = false;
      for(unsigned int j = 0; j < group->size() && !test; j++)
      {
        if(isConst[TextTools::to<unsigned int>((*group)[j])]) test = true;
      }

      // Compute distance from mean vector:
      vector<double> groupMeanVector(nbBranches, 0.);
      for(unsigned int j = 0; j < group->size(); j++)
      {
        groupMeanVector += (*mapping1)[TextTools::to<unsigned int>((*group)[j])]; 
      }
      double distFromMeanVector = dist->d(groupMeanVector/group->size(), meanVector);   

      minNorms[i] = minNorm;
      // Store results:
      vector<string> row(6);
      row[0] = group->toString(siteNames);
      row[1] = TextTools::toString(group->size());
      row[2] = test ? "yes" : "no";
      row[3] = TextTools::toString(group->getHeight() * 2.);
      row[4] = TextTools::toString(minNorm);
      row[5] = TextTools::toString(distFromMeanVector);
      groupsData.addRow(row);
    }

    // Now test each group:
  
    bool testGroupsGlobal = ApplicationTools::getBooleanParameter("clustering.null", params, false, "", true, false);
    if(testGroupsGlobal)
    {
      HomogeneousSequenceSimulator simulator(process1, rDist1, tree1, true);

	    string simPath = ApplicationTools::getAFilePath("clustering.null.output.file", params, true, false, "", false);
      unsigned int nrep = ApplicationTools::getParameter<unsigned int>("clustering.null.number", params, 1, "", true);
      ofstream * out = NULL;
      if(simPath != "none") out = new ofstream(simPath.c_str(), ios::out);
      ContingencyTable table(50, 0.5, 20, 0.1, 15);
      ClusterTools::computeGlobalDistanceDistribution(simulator, *substitutionCount1, *dist, *clustering, scales, nbSites, nrep, table, out);
      if(out != NULL) out->close();
      //// Now write p-values:
      //vector<string> pvalues(groups.size());
      //for(unsigned int i = 0; i < groups.size(); i++)
      //{
      //  Group * group = &groups[i];
      //  string pvalue = "NA";
      //  try {
      //    pvalue = TextTools::toString(table.getPValue(group->size(), minNorms[i], max(0., group->getHeight() * 2.)));
      //  } catch(Exception & ex) {
      //    cout << "WARNING: Could not compute p-value for group " << group->toString() << " (value out of table range)." << endl;
      //  }
      //  pvalues[i] = pvalue;
      //}
      //groupsData.addColumn("p-values", pvalues);
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
      ClusterTools::translate(*dynamic_cast<TreeTemplate<Node> *>(clusteringTree), siteNames);
      Newick newick;
	    newick.write(*clusteringTree, treeFile, true);
	    ApplicationTools::displayResult("Wrote tree to file", treeFile);
    }
    
    // Clustering memory freeing:
	  delete dist;
	  delete mat;
	  delete clustering;
	  delete clusteringTree;
  }





  
	ApplicationTools::displayTask("CoMap's done. Freeing memory");
  
  // Dataset 1 memory freeing:
	delete tree1;
	delete allSites1;
	delete sites1;
	delete model1;
	delete rDist1;
	delete prDist1;
	delete tl1;
	delete substitutionCount1;
	delete process1;
  delete mapping1;

	ApplicationTools::displayMessage("Bye bye ;-)");	

	}
  catch(Exception e)
  {
		cout << e.what() << endl;
		exit(-1);
	}

	return (0);
}

