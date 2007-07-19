//
// File: CoETools.h
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

#ifndef _COETOOLS_H_
#define _COETOOLS_H_

#include "Statistics.h"

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/VectorTools.h>

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/SubstitutionModel.h>
#include <Phyl/HomogeneousTreeLikelihood.h>
#include <Phyl/DRHomogeneousTreeLikelihood.h>
#include <Phyl/MutationProcess.h>
#include <Phyl/DiscreteRatesAcrossSitesTreeLikelihood.h>
#include <Phyl/SubstitutionCount.h>
#include <Phyl/ProbabilisticSubstitutionMapping.h>

// From the STL:
#include <map>
using namespace std;

class CoETools
{
	public:
		CoETools() {};
		virtual ~CoETools() {};
	
	public:
		
		static void readData(
			TreeTemplate<Node> *          &tree,
			Alphabet *                    &alphabet,
			VectorSiteContainer *         &allSites,
			VectorSiteContainer *         &sites,
			SubstitutionModel *           &model,
			DiscreteDistribution *        &rDist,
			DRHomogeneousTreeLikelihood * &tl,
			map<string, string>           &params,
			const string                  &suffix = "");

		static ProbabilisticSubstitutionMapping * getVectors(
			const Alphabet * alphabet,
				    TreeTemplate<Node> & tree,
			const SiteContainer      & completeSites,
			const SiteContainer      & sites,
				  SubstitutionModel    & model,
				  DiscreteDistribution & rDist,
			const SubstitutionCount  & substitutionCount,
			map<string, string>      & params,
			const string             & suffix = "");
			
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
			const SiteContainer & completeSites,
			const DiscreteRatesAcrossSitesTreeLikelihood & ras,
			map<string, string> & params,
			const string & suffix = "");

		static int getMinRateClass(map<string, string> & params, string suffix = "");
		
		static double  getMinRate(map<string, string> & params, string suffix = "");
		
		static int getMaxRateClassDiff(map<string, string> & params);
		
		static double getMaxRateDiff(map<string, string> & params);

		static double getStatisticMin(map<string, string> & params);
		
		static bool haveToPerformIndependantComparisons(map<string, string> & params);
		
		static const Statistic * getStatistic(map<string, string> & params);
		
		static SubstitutionCount * getSubstitutionCount(
			const Alphabet * chars,
			const TreeTemplate<Node> & tree,
			const MutationProcess & process,
			const DiscreteDistribution & rDist,
			map<string, string> & params,
			string suffix = "");

    static void computeIntraStats(
			const DiscreteRatesAcrossSitesTreeLikelihood & tl,
			const SiteContainer & completeSites,
			ProbabilisticSubstitutionMapping & mapping,
			const Statistic & statistic,
			map<string, string> & params);

		static void computeInterStats(
			const DiscreteRatesAcrossSitesTreeLikelihood & tl1,
			const DiscreteRatesAcrossSitesTreeLikelihood & tl2,
			const SiteContainer & completeSites1,
			const SiteContainer & completeSites2,
			ProbabilisticSubstitutionMapping & mapping1,
			ProbabilisticSubstitutionMapping & mapping2,
			const Statistic & statistic,
			map<string, string> & params);
		
		static void computeIntraNullDistribution(
			const MutationProcess & process,
			const DiscreteDistribution & rDist,
			const TreeTemplate<Node> & tree,
			const SubstitutionCount & nijt,
			const Statistic & statistic,
			map<string, string> & params,
      bool useContinuousRates = false);
	
		static void computeInterNullDistribution(
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
      bool useContinuousRates = false);
	
};

#endif	//_COETOOLS_H_

