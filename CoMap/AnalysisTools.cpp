//
// File: AnalysisTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Nov 28 16:33:00 2003
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

#include "AnalysisTools.h"

// From the NumCalc library:
#include <NumCalc/VectorTools.h>

// From Utils
#include <Utils/ApplicationTools.h>

// From PhylLib:
#include <Phyl/SubstitutionMappingTools.h>

/******************************************************************************/

AnalysisTools::AnalysisTools() {}

AnalysisTools::~AnalysisTools() {}
	
/******************************************************************************/

VVdouble AnalysisTools::getFromStream(istream & in)
{
	string token;
	VVdouble analysis;
	in >> token;
	if(token != "Site") {
		cerr << "Error while reading result file. Bad header line." << endl;
		exit(-1);
	}
	//This is the first line of file.
	//Initialise the result array:
	while(in && token.find("Branch") == string::npos) {
		analysis.push_back(Vdouble());
	}
		
	int size = analysis.size();
	if(size == 0) {
		cerr << "Error while reading result file. Bad header line or no data." << endl;
		exit(-1);
	}

	double value;
	while(in) {
		cin >> token;
		if(token.find("Branch") == string::npos) {
			cerr << "Error while reading result file. Bad branch data line." << endl;
			exit(-1);
		}
		//begin to parse a new line:
		for(int i = 0; i < size + 1; i++) {
			if(!in) {
				cerr << "Error while reading result file. Incomplete branch data line." << endl;
				exit(-1);
			}
			cin >> value;
			analysis[i].push_back(value);
		}
	}
	return analysis;
}

/******************************************************************************/

VVdouble AnalysisTools::computeScalarProductMatrix(const VVdouble & vectors)
{
	unsigned int nbVectors = vectors.size();
	VVdouble matrix = VVdouble(nbVectors);
	for(unsigned int i = 0; i < nbVectors; i++)
  {
		//Must initialize all vector first, since we access
		//both matrix[i, j] and matrix[j, i] in the same time.
		matrix[i] = Vdouble(nbVectors, 0);
	}
	for(unsigned int i = 0; i < nbVectors; i++)
  {
		matrix[i][i] = VectorTools::scalar(vectors[i], vectors[i]);
		for(unsigned int j = i + 1; j < nbVectors; j++)
    {
			matrix[i][j] = matrix[j][i] = VectorTools::scalar(vectors[i], vectors[j]);
		}
	}
	return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeScalarProductMatrix (
	const VVdouble & vectors1,
	const VVdouble & vectors2,
	bool independantComparisons
) throw (DimensionException)
{
	unsigned int nbVectors1 = vectors1.size();
	unsigned int nbVectors2 = vectors2.size();
	if(independantComparisons && nbVectors1 != nbVectors2)
  {
		throw DimensionException(
			string("AnalysisTools::computeScalarProductMatrix.\n") +
			string("When performing independant comparisons, ") +
			string("the two datasets must have the same length."),
			nbVectors1,
			nbVectors2
		);
	}
	VVdouble matrix = VVdouble(nbVectors1);
	for(unsigned int i = 0; i < nbVectors1; i++)
  {
		//Must initialize all vector first, since we access
		matrix[i] = Vdouble(nbVectors2, 0);
	}
	for(unsigned int i = 0; i < nbVectors1; i++)
  {
		unsigned int begin = independantComparisons ? i : 0;
		unsigned int end   = independantComparisons ? i + 1 : nbVectors2;
		for(unsigned int j = begin; j < end; j++)
    {
			matrix[i][j] = VectorTools::scalar(vectors1[i], vectors2[j]);
		}
	}
	return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCosinusMatrix(const VVdouble & vectors)
{
	unsigned int nbVectors = vectors.size();
	VVdouble matrix = VVdouble(nbVectors);
	for(unsigned int i = 0; i < nbVectors; i++)
  {
		//Must initialize all vector first, since we access
		//both matrix[i, j] and matrix[j, i] in the same time.
		matrix[i] = Vdouble(nbVectors, 0);
	}
	for(unsigned int i = 0; i < nbVectors; i++)
  {
		matrix[i][i] = 1;
		for(unsigned int j = i + 1; j < nbVectors; j++)
    {
			matrix[i][j] = matrix[j][i] = VectorTools::cos(vectors[i], vectors[j]);
		}
	}
	return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCosinusMatrix(
	const VVdouble & vectors1,
	const VVdouble & vectors2,
	bool independantComparisons
) throw (DimensionException)
{
	unsigned int nbVectors1 = vectors1.size();
	unsigned int nbVectors2 = vectors2.size();
	if(independantComparisons && nbVectors1 != nbVectors2)
  {
		throw DimensionException(
			string("AnalysisTools::computeCosinusMatrix.\n") +
			string("When performing independant comparisons, ") +
			string("the two datasets must have the same length."),
			nbVectors1,
			nbVectors2
		);
	}
	VVdouble matrix = VVdouble(nbVectors1);
	for(unsigned int i = 0; i < nbVectors1; i++)
  {
		//Must initialize all vector first, since we access
		matrix[i] = Vdouble(nbVectors2, 0);
	}
	for(unsigned int i = 0; i < nbVectors1; i++)
  {
		unsigned int begin = independantComparisons ? i : 0;
		unsigned int end   = independantComparisons ? i + 1 : nbVectors2;
		for(unsigned int j = begin; j < end; j++)
    {
			matrix[i][j] = VectorTools::cos(vectors1[i], vectors2[j]);
		}
	}
	return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCorrelationMatrix(const VVdouble & vectors)
{
	unsigned int nbVectors = vectors.size();
	VVdouble matrix = VVdouble(nbVectors);
	for(unsigned int i = 0; i < nbVectors; i++)
  {
		//Must initialize all vector first, since we access
		//both matrix[i, j] and matrix[j, i] in the same time.
		matrix[i] = Vdouble(nbVectors, 0);
	}
	for(unsigned int i = 0; i < nbVectors; i++)
  {
		matrix[i][i] = 1;
		for(unsigned int j = i + 1; j < nbVectors; j++)
    {
			matrix[i][j] = matrix[j][i] = VectorTools::cor(vectors[i], vectors[j]);
		}
	}
	return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCorrelationMatrix(
	const VVdouble & vectors1,
	const VVdouble & vectors2,
	bool independantComparisons
) throw (DimensionException)
{
	unsigned int nbVectors1 = vectors1.size();
	unsigned int nbVectors2 = vectors2.size();
	if(independantComparisons && nbVectors1 != nbVectors2)
  {
		throw DimensionException(
			string("AnalysisTools::computeCorrelationMatrix.\n") +
			string("When performing independant comparisons, ") +
			string("the two datasets must have the same length."),
			nbVectors1,
			nbVectors2
		);
	}
	VVdouble matrix = VVdouble(nbVectors1);
	for(unsigned int i = 0; i < nbVectors1; i++)
  {
		//Must initialize all vector first, since we access
		matrix[i] = Vdouble(nbVectors2, 0);
	}
	for(unsigned int i = 0; i < nbVectors1; i++)
  {
		unsigned int begin = independantComparisons ? i : 0;
		unsigned int end   = independantComparisons ? i + 1 : nbVectors2;
		for(unsigned int j = begin; j < end; j++) 
    {
			matrix[i][j] = VectorTools::cor(vectors1[i], vectors2[j]);
		}
	}
	return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCovarianceMatrix(const VVdouble & vectors)
{
	unsigned int nbVectors = vectors.size();
	VVdouble matrix = VVdouble(nbVectors);
	for(unsigned int i = 0; i < nbVectors; i++)
  {
		//Must initialize all vector first, since we access
		//both matrix[i, j] and matrix[j, i] in the same time.
		matrix[i] = Vdouble(nbVectors, 0);
	}
	for(unsigned int i = 0; i < nbVectors; i++)
  {
		matrix[i][i] = VectorTools::var(vectors[i]);
		for(unsigned int j = i + 1; j < nbVectors; j++)
    {
			matrix[i][j] = matrix[j][i] = VectorTools::cov(vectors[i], vectors[j]);
		}
	}
	return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCovarianceMatrix(
	const VVdouble & vectors1,
	const VVdouble & vectors2,
	bool independantComparisons
) throw (DimensionException)
{
	unsigned int nbVectors1 = vectors1.size();
	unsigned int nbVectors2 = vectors2.size();
	if(independantComparisons && nbVectors1 != nbVectors2)
  {
		throw DimensionException(
			string("AnalysisTools::computeCovarianceMatrix.\n") +
			string("When performing independant comparisons, ") +
			string("the two datasets must have the same length."),
			nbVectors1,
			nbVectors2
		);
	}
	VVdouble matrix = VVdouble(nbVectors1);
	for(unsigned int i = 0; i < nbVectors1; i++)
  {
		//Must initialize all vector first, since we access
		matrix[i] = Vdouble(nbVectors2, 0);
	}
	for(unsigned int i = 0; i < nbVectors1; i++)
  {
		unsigned int begin = independantComparisons ? i : 0;
		unsigned int end   = independantComparisons ? i + 1 : nbVectors2;
		for(unsigned int j = begin; j < end; j++)
    {
			matrix[i][j] = VectorTools::cov(vectors1[i], vectors2[j]);
		}
	}
	return matrix;
}

/******************************************************************************/

Vdouble AnalysisTools::computeNorms(const ProbabilisticSubstitutionMapping & mapping)
{
	unsigned int nbVectors = mapping.getNumberOfSites();
	Vdouble vect(nbVectors);
	for(unsigned int i = 0; i < nbVectors; i++) vect[i] = VectorTools::norm(mapping[i]);
	return vect;
}

/******************************************************************************/

void AnalysisTools::writeMatrix(
	const VVdouble & matrix,
	const SiteContainer & sites,
	ostream & out)
{
	out << "\tMean";
	for(unsigned int i = 0; i < matrix.size() - 1; i++) {
		out << "\tSite" << sites.getSite(i) -> getPosition();
	}
	out << endl;
	for(unsigned int j = 0; j < matrix[0].size(); j++) {
		if(j == 0) out << "Mean";
		else out << "Site" << sites.getSite(j - 1) -> getPosition();
		for(unsigned int i = 0; i < matrix.size(); i++) {
			out << "\t" << matrix[i][j];
		}
		out << endl;
	}
}

/******************************************************************************/

void AnalysisTools::writeMatrix(
	const VVdouble & matrix,
	const SiteContainer & sites1,
	const SiteContainer & sites2,
	ostream & out)
{
	out << "\tMean";
	for(unsigned int i = 0; i < matrix.size() - 1; i++) {
		out << "\tSite" << sites1.getSite(i) -> getPosition();
	}
	out << endl;
	for(unsigned int j = 0; j < matrix[0].size(); j++) {
		if(j == 0) out << "Mean";
		else out << "Site" << sites2.getSite(j - 1) -> getPosition();
		for(unsigned int i = 0; i < matrix.size(); i++) {
			out << "\t" << matrix[i][j];
		}
		out << endl;
	}
}

/******************************************************************************/

vector<IntervalData *> AnalysisTools::getNullDistributionIntraDR(
	const HomogeneousSequenceSimulator & seqSim,
	const SubstitutionCount & nijt,
	const Statistic & statistic,
	const Domain & statDomain,
	const Domain & rateDomain,
	unsigned int repCPU,
	unsigned int repRAM,
	bool average,
	bool joint,
	bool verbose)
{
	unsigned int nbClasses = rateDomain.getSize();
	vector<IntervalData *> id(nbClasses);
	for(unsigned int i = 0; i < nbClasses; i++)
  {
		id[i] = new IntervalData(statDomain, "Null distribution (classe > " + TextTools::toString(rateDomain.getValue(i)) + ")");
	}
	for(unsigned int i = 0; i < repCPU; i++)
  {
		if(verbose)
    {
			*ApplicationTools::message << ".";
			ApplicationTools::message->flush();
		}
		SiteContainer * sites1 = seqSim.simulate(repRAM);
		SiteContainer * sites2 = seqSim.simulate(repRAM);
		DRHomogeneousTreeLikelihood * drhtl1 =
			new DRHomogeneousTreeLikelihood(
				*seqSim.getTree(),
				*sites1,
				const_cast<SubstitutionModel *>(seqSim.getSubstitutionModel()),
				const_cast<DiscreteDistribution *>(seqSim.getRateDistribution()),
				false, false);
    drhtl1->initialize();
		DRHomogeneousTreeLikelihood * drhtl2 =
			new DRHomogeneousTreeLikelihood(
				*seqSim.getTree(),
				*sites2,
				const_cast<SubstitutionModel *>(seqSim.getSubstitutionModel()),
				const_cast<DiscreteDistribution *>(seqSim.getRateDistribution()),
				false, false);
    drhtl2->initialize();
		//drhtl1->computeTreeLikelihood();
		//drhtl2->computeTreeLikelihood();
		Vdouble r1 = drhtl1->getRateWithMaxPostProbOfEachSite();
		Vdouble r2 = drhtl2->getRateWithMaxPostProbOfEachSite();
    ProbabilisticSubstitutionMapping * mapping1;
		ProbabilisticSubstitutionMapping * mapping2;
		if(average)
    {
			if(joint)
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectors(* drhtl1, nijt, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectors(* drhtl2, nijt, false);
			}
      else
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(* drhtl1, nijt, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(* drhtl2, nijt, false);
			}
		}
    else
    {
			if(joint)
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(* drhtl1, nijt, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(* drhtl2, nijt, false);
			}
      else
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(* drhtl1, nijt, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(* drhtl2, nijt, false);
			}
		}
   	for(unsigned int j = 0; j < repRAM; j++)
    {
			double stat = statistic.getValue((*mapping1)[j], (*mapping2)[j]);
			double minR = min(r1[j], r2[j]);
			//int minRc = min (rc1[j], rc2[j]);
			//cout << r1[j] << "\t" << r2[j] << "\t" << minR << "\t" << rateDomain.getIndex(minR) << endl;
			//cout << rc1[j] << "\t" << rc2[j] << "\t" << minRc << "\t" << endl;
			id[rateDomain.getIndex(minR)]->addValue(stat);
			//id[minRc] -> addValue(stat);
		}
		delete sites1;
		delete sites2;
		delete drhtl1;
		delete drhtl2;
    delete mapping1;
    delete mapping2;
	}
	return id;
}

/******************************************************************************/

vector<IntervalData *> AnalysisTools::getNullDistributionInterDR(
	const HomogeneousSequenceSimulator & seqSim1,
	const HomogeneousSequenceSimulator & seqSim2,
	const SubstitutionCount & nijt1,
	const SubstitutionCount & nijt2,
	const Statistic & statistic,
	const Domain & statDomain,
  const Domain & rateDomain,
	unsigned int repCPU,
  unsigned int repRAM,
	bool average,
	bool joint,
	bool verbose)
{
  unsigned int nbClasses = rateDomain.getSize();
  vector<IntervalData *> id(nbClasses);
  for(unsigned int i = 0; i < nbClasses; i++)
  {
    id[i] = new IntervalData(statDomain, "Null distribution (classe > " + TextTools::toString(rateDomain.getValue(i)) + ")");
  }
	for(unsigned int i = 0; i < repCPU; i++)
  {
		if(verbose)
    {
			*ApplicationTools::message << ".";
			ApplicationTools::message->flush();
		}
		SiteContainer * sites1 = seqSim1.simulate(repRAM);
		SiteContainer * sites2 = seqSim2.simulate(repRAM);
		DRHomogeneousTreeLikelihood * drhtl1 = 
			new DRHomogeneousTreeLikelihood(
				*seqSim1.getTree(),
				*sites1,
				const_cast<SubstitutionModel *>(seqSim1.getSubstitutionModel()),
				const_cast<DiscreteDistribution *>(seqSim1.getRateDistribution()),
        false, false
			);
    drhtl1->initialize();
		DRHomogeneousTreeLikelihood * drhtl2 = 
			new DRHomogeneousTreeLikelihood(
				*seqSim2.getTree(),
				*sites2,
				const_cast<SubstitutionModel *>(seqSim2.getSubstitutionModel()),
				const_cast<DiscreteDistribution *>(seqSim2.getRateDistribution()),
				false, false
			);
    drhtl2->initialize();
    //drhtl1->computeTreeLikelihood();
    //drhtl2->computeTreeLikelihood();
    Vdouble r1 = drhtl1->getRateWithMaxPostProbOfEachSite();
    Vdouble r2 = drhtl2->getRateWithMaxPostProbOfEachSite();

    ProbabilisticSubstitutionMapping * mapping1;
		ProbabilisticSubstitutionMapping * mapping2;
		if(average)
    {
			if(joint)
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectors(* drhtl1, nijt1, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectors(* drhtl2, nijt2, false);
			}
      else
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(* drhtl1, nijt1, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(* drhtl2, nijt2, false);
			}
		}
    else
    {
			if(joint)
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(* drhtl1, nijt1, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(* drhtl2, nijt2, false);
			}
      else
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(* drhtl1, nijt1, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(* drhtl2, nijt2, false);
			}
		}
    
		for(unsigned int j = 0; j < repRAM; j++)
    {
      double stat = statistic.getValue((*mapping1)[j], (*mapping2)[j]);
      double minR = min(r1[j], r2[j]);
      id[rateDomain.getIndex(minR)] -> addValue(stat);
    }
    delete sites1;
    delete sites2;
		delete drhtl1;
		delete drhtl2;
    delete mapping1;
    delete mapping2;
	}
	return id;
}


/******************************************************************************/

void AnalysisTools::getNullDistributionIntraDR(
	const HomogeneousSequenceSimulator & seqSim,
	const SubstitutionCount & nijt,
	const Statistic & statistic,
	ostream & out,
	unsigned int repCPU,
	unsigned int repRAM,
	bool average,
	bool joint,
	bool verbose)
{
	// Write header line:
	out << "statistic\tmin.rc\tmin.pr\tNmin" << endl;
	for(unsigned int i = 0; i < repCPU; i++)
  {
		if(verbose)
    {
			*ApplicationTools::message << ".";
			ApplicationTools::message->flush();
		}
		//cout << "Creating datasets..." << endl;
		SiteContainer * sites1 = seqSim.simulate(repRAM);
		SiteContainer * sites2 = seqSim.simulate(repRAM);
		//cout << "Computing likelihoods..." << endl;
		DRHomogeneousTreeLikelihood * drhtl1 =
			new DRHomogeneousTreeLikelihood(
				*seqSim.getTree(),
				*sites1,
				const_cast<SubstitutionModel *>(seqSim.getSubstitutionModel()),
				const_cast<DiscreteDistribution *>(seqSim.getRateDistribution()),
				false, false
			);
    drhtl1->initialize();
		DRHomogeneousTreeLikelihood * drhtl2 =
			new DRHomogeneousTreeLikelihood(
				*seqSim.getTree(),
				*sites2,
				const_cast<SubstitutionModel *>(seqSim.getSubstitutionModel()),
				const_cast<DiscreteDistribution *>(seqSim.getRateDistribution()),
				false, false
			);
    drhtl2->initialize();
		//drhtl1->computeTreeLikelihood();
		//drhtl2->computeTreeLikelihood();
		//cout << "Computing posterior rates..." << endl;
		vector<unsigned int> rc1 = drhtl1->getRateClassWithMaxPostProbOfEachSite();
		vector<unsigned int> rc2 = drhtl2->getRateClassWithMaxPostProbOfEachSite();
		Vdouble pr1 = drhtl1->getPosteriorRateOfEachSite();
		Vdouble pr2 = drhtl2->getPosteriorRateOfEachSite();
		
		//cout << "Computing vectors..." << endl;
    ProbabilisticSubstitutionMapping * mapping1;
		ProbabilisticSubstitutionMapping * mapping2;
		if(average)
    {
			if(joint)
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectors(* drhtl1, nijt, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectors(* drhtl2, nijt, false);
			}
      else
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(* drhtl1, nijt, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(* drhtl2, nijt, false);
			}
		}
    else
    {
			if(joint)
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(* drhtl1, nijt, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(* drhtl2, nijt, false);
			}
      else
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(* drhtl1, nijt, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(* drhtl2, nijt, false);
			}
		}
    Vdouble norm1 = computeNorms(*mapping1);
    Vdouble norm2 = computeNorms(*mapping2);
   	//cout << "Computing statistics..." << endl;
		for(unsigned int j = 0; j < repRAM; j++)
    {
			double stat = statistic.getValue((*mapping1)[j], (*mapping2)[j]);
			out << stat << "\t" << std::min(rc1[j], rc2[j]) << "\t" << std::min(pr1[j], pr2[j]) << "\t" << std::min(norm1[j], norm2[j]) << endl;
		}
		//cout << "Freeing memory." << endl;
		delete sites1;
		delete sites2;
		delete drhtl1;
		delete drhtl2;
    delete mapping1;
    delete mapping2;
	}
}

/******************************************************************************/

void AnalysisTools::getNullDistributionInterDR(
	const HomogeneousSequenceSimulator & seqSim1,
	const HomogeneousSequenceSimulator & seqSim2,
	const SubstitutionCount & nijt1,
	const SubstitutionCount & nijt2,
	const Statistic & statistic,
	ostream & out,
	unsigned int repCPU,
	unsigned int repRAM,
	bool average,
	bool joint,
	bool verbose)
{
	// Write header line:
	out << "statistic\tmin.rc\tmin.pr\tNmin" << endl;
	for(unsigned int i = 0; i < repCPU; i++)
  {
		if(verbose)
    {
			*ApplicationTools::message << ".";
			ApplicationTools::message->flush();
		}
		SiteContainer * sites1 = seqSim1.simulate(repRAM);
		SiteContainer * sites2 = seqSim2.simulate(repRAM);
		DRHomogeneousTreeLikelihood * drhtl1 =
			new DRHomogeneousTreeLikelihood(
				*seqSim1.getTree(),
				*sites1,
				const_cast<SubstitutionModel *>(seqSim1.getSubstitutionModel()),
				const_cast<DiscreteDistribution *>(seqSim1.getRateDistribution()),
				false, false
			);
    drhtl1->initialize();
		DRHomogeneousTreeLikelihood * drhtl2 =
			new DRHomogeneousTreeLikelihood(
				*seqSim2.getTree(),
				*sites2,
				const_cast<SubstitutionModel *>(seqSim2.getSubstitutionModel()),
				const_cast<DiscreteDistribution *>(seqSim2.getRateDistribution()),
				false, false
			);
    drhtl2->initialize();
		//drhtl1->computeTreeLikelihood();
		//drhtl2->computeTreeLikelihood();
		vector<unsigned int> rc1 = drhtl1->getRateClassWithMaxPostProbOfEachSite();
		vector<unsigned int> rc2 = drhtl2->getRateClassWithMaxPostProbOfEachSite();
		Vdouble pr1 = drhtl1->getPosteriorRateOfEachSite();
		Vdouble pr2 = drhtl2->getPosteriorRateOfEachSite();
		
    ProbabilisticSubstitutionMapping * mapping1;
		ProbabilisticSubstitutionMapping * mapping2;
		if(average)
    {
			if(joint)
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectors(* drhtl1, nijt1, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectors(* drhtl2, nijt2, false);
			}
      else
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(* drhtl1, nijt1, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(* drhtl2, nijt2, false);
			}
		}
    else
    {
			if(joint)
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(* drhtl1, nijt1, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(* drhtl2, nijt2, false);
			}
      else
      {
				mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(* drhtl1, nijt1, false);
				mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(* drhtl2, nijt2, false);
			}
		}
    Vdouble norm1 = computeNorms(*mapping1);
    Vdouble norm2 = computeNorms(*mapping2);
   	for(unsigned int j = 0; j < repRAM; j++)
    {
			double stat = statistic.getValue((*mapping1)[j], (*mapping2)[j]);
			out << stat << "\t" << std::min(rc1[j], rc2[j]) << "\t" << std::min(pr1[j], pr2[j]) << "\t" << std::min(norm1[j], norm2[j]) << endl;
		}

		delete sites1;
		delete sites2;
		delete drhtl1;
		delete drhtl2;
    delete mapping1;
    delete mapping2;
	}
}

/******************************************************************************/

void AnalysisTools::getNullDistributionIntraWithoutReestimatingCounts(
	const HomogeneousSequenceSimulator & seqSim,
	const Statistic & statistic,
	ostream & out,
	unsigned int rep,
	bool verbose)
{
	// Write header line:
	out << "statistic\tr1\tr2" << endl;
	unsigned int c = 0;
	for(unsigned int i = 0; i < rep; i++)
  {
		if(c == 1000) { cout << "."; cout.flush(); c = 0; }
		c++;
		HomogeneousSiteSimulationResult * hssr1 = 
			dynamic_cast<HomogeneousSiteSimulationResult *>(seqSim.dSimulate());
		vector<double> vector1 = hssr1->getSubstitutionVector();
		double         rate1   = hssr1->getRate();
		delete hssr1;
		HomogeneousSiteSimulationResult * hssr2 = 
			dynamic_cast<HomogeneousSiteSimulationResult *>(seqSim.dSimulate());
		vector<double> vector2 = hssr2->getSubstitutionVector();
		double         rate2   = hssr2->getRate(); 
		delete hssr2;

		double stat = statistic.getValue(vector1, vector2);
		out << stat << "\t" << rate1 << "\t" << rate2 << endl;
	}
}

/******************************************************************************/

