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
#include "IntervalData.h"

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>

// From PhylLib:
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>

using namespace bpp;

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
    
  size_t size = analysis.size();
  if (size == 0) {
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
    for(size_t i = 0; i < size + 1; i++) {
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
  size_t nbVectors = vectors.size();
  VVdouble matrix = VVdouble(nbVectors);
  for(size_t i = 0; i < nbVectors; i++)
  {
    //Must initialize all vector first, since we access
    //both matrix[i, j] and matrix[j, i] in the same time.
    matrix[i] = Vdouble(nbVectors, 0);
  }
  for(size_t i = 0; i < nbVectors; i++)
  {
    matrix[i][i] = VectorTools::scalar<double, double>(vectors[i], vectors[i]);
    for(size_t j = i + 1; j < nbVectors; j++)
    {
      matrix[i][j] = matrix[j][i] = VectorTools::scalar<double, double>(vectors[i], vectors[j]);
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
  size_t nbVectors1 = vectors1.size();
  size_t nbVectors2 = vectors2.size();
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
  for(size_t i = 0; i < nbVectors1; i++)
  {
    //Must initialize all vector first, since we access
    matrix[i] = Vdouble(nbVectors2, 0);
  }
  for(size_t i = 0; i < nbVectors1; i++)
  {
    size_t begin = independantComparisons ? i : 0;
    size_t end   = independantComparisons ? i + 1 : nbVectors2;
    for(size_t j = begin; j < end; j++)
    {
      matrix[i][j] = VectorTools::scalar<double, double>(vectors1[i], vectors2[j]);
    }
  }
  return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCosinusMatrix(const VVdouble & vectors)
{
  size_t nbVectors = vectors.size();
  VVdouble matrix = VVdouble(nbVectors);
  for(size_t i = 0; i < nbVectors; i++)
  {
    //Must initialize all vector first, since we access
    //both matrix[i, j] and matrix[j, i] in the same time.
    matrix[i] = Vdouble(nbVectors, 0);
  }
  for(size_t i = 0; i < nbVectors; i++)
  {
    matrix[i][i] = 1;
    for(size_t j = i + 1; j < nbVectors; j++)
    {
      matrix[i][j] = matrix[j][i] = VectorTools::cos<double, double>(vectors[i], vectors[j]);
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
  size_t nbVectors1 = vectors1.size();
  size_t nbVectors2 = vectors2.size();
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
  for(size_t i = 0; i < nbVectors1; i++)
  {
    //Must initialize all vector first, since we access
    matrix[i] = Vdouble(nbVectors2, 0);
  }
  for(size_t i = 0; i < nbVectors1; i++)
  {
    size_t begin = independantComparisons ? i : 0;
    size_t end   = independantComparisons ? i + 1 : nbVectors2;
    for(size_t j = begin; j < end; j++)
    {
      matrix[i][j] = VectorTools::cos<double, double>(vectors1[i], vectors2[j]);
    }
  }
  return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCorrelationMatrix(const VVdouble & vectors)
{
  size_t nbVectors = vectors.size();
  VVdouble matrix = VVdouble(nbVectors);
  for(size_t i = 0; i < nbVectors; i++)
  {
    //Must initialize all vector first, since we access
    //both matrix[i, j] and matrix[j, i] in the same time.
    matrix[i] = Vdouble(nbVectors, 0);
  }
  for(size_t i = 0; i < nbVectors; i++)
  {
    matrix[i][i] = 1;
    for(size_t j = i + 1; j < nbVectors; j++)
    {
      matrix[i][j] = matrix[j][i] = VectorTools::cor<double, double>(vectors[i], vectors[j]);
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
  size_t nbVectors1 = vectors1.size();
  size_t nbVectors2 = vectors2.size();
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
  for(size_t i = 0; i < nbVectors1; i++)
  {
    //Must initialize all vector first, since we access
    matrix[i] = Vdouble(nbVectors2, 0);
  }
  for(size_t i = 0; i < nbVectors1; i++)
  {
    size_t begin = independantComparisons ? i : 0;
    size_t end   = independantComparisons ? i + 1 : nbVectors2;
    for(size_t j = begin; j < end; j++) 
    {
      matrix[i][j] = VectorTools::cor<double, double>(vectors1[i], vectors2[j]);
    }
  }
  return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCovarianceMatrix(const VVdouble & vectors)
{
  size_t nbVectors = vectors.size();
  VVdouble matrix = VVdouble(nbVectors);
  for(size_t i = 0; i < nbVectors; i++)
  {
    //Must initialize all vector first, since we access
    //both matrix[i, j] and matrix[j, i] in the same time.
    matrix[i] = Vdouble(nbVectors, 0);
  }
  for(size_t i = 0; i < nbVectors; i++)
  {
    matrix[i][i] = VectorTools::var<double, double>(vectors[i]);
    for(size_t j = i + 1; j < nbVectors; j++)
    {
      matrix[i][j] = matrix[j][i] = VectorTools::cov<double, double>(vectors[i], vectors[j]);
    }
  }
  return matrix;
}

/******************************************************************************/

VVdouble AnalysisTools::computeCovarianceMatrix(
  const VVdouble& vectors1,
  const VVdouble& vectors2,
  bool independantComparisons
) throw (DimensionException)
{
  size_t nbVectors1 = vectors1.size();
  size_t nbVectors2 = vectors2.size();
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
  for(size_t i = 0; i < nbVectors1; i++)
  {
    //Must initialize all vector first, since we access
    matrix[i] = Vdouble(nbVectors2, 0);
  }
  for(size_t i = 0; i < nbVectors1; i++)
  {
    size_t begin = independantComparisons ? i : 0;
    size_t end   = independantComparisons ? i + 1 : nbVectors2;
    for(size_t j = begin; j < end; j++)
    {
      matrix[i][j] = VectorTools::cov<double, double>(vectors1[i], vectors2[j]);
    }
  }
  return matrix;
}

/******************************************************************************/

Vdouble AnalysisTools::computeNorms(const ProbabilisticSubstitutionMapping& mapping)
{
  size_t nbVectors = mapping.getNumberOfSites();
  Vdouble vect(nbVectors);
  for(size_t i = 0; i < nbVectors; i++)
    vect[i] = SubstitutionMappingTools::computeNormForSite(mapping, i);
  return vect;
}

/******************************************************************************/

void AnalysisTools::writeMatrix(
  const VVdouble& matrix,
  const SiteContainer& sites,
  ostream& out)
{
  out << "\tMean";
  for(size_t i = 0; i < matrix.size() - 1; i++)
  {
    out << "\tSite" << sites.getSite(i).getPosition();
  }
  out << endl;
  for(size_t j = 0; j < matrix[0].size(); j++)
  {
    if(j == 0) out << "Mean";
    else out << "Site" << sites.getSite(j - 1).getPosition();
    for(size_t i = 0; i < matrix.size(); i++)
    {
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
  for(size_t i = 0; i < matrix.size() - 1; i++)
  {
    out << "\tSite" << sites1.getSite(i).getPosition();
  }
  out << endl;
  for(size_t j = 0; j < matrix[0].size(); j++)
  {
    if(j == 0) out << "Mean";
    else out << "Site" << sites2.getSite(j - 1).getPosition();
    for(size_t i = 0; i < matrix.size(); i++)
    {
      out << "\t" << matrix[i][j];
    }
    out << endl;
  }
}

/******************************************************************************/

vector<IntervalData*> AnalysisTools::getNullDistributionIntraDR(
  DRTreeLikelihood& drtl,
  const SequenceSimulator& seqSim,
  SubstitutionCount& nijt,
  const Statistic& statistic,
  const Domain& statDomain,
  const Domain& rateDomain,
  size_t repCPU,
  size_t repRAM,
  bool average,
  bool joint,
  bool verbose)
{
  size_t nbClasses = rateDomain.getSize();
  vector<IntervalData *> id(nbClasses);
  for(size_t i = 0; i < nbClasses; i++)
  {
    id[i] = new IntervalData(statDomain, "Null distribution (classe > " + TextTools::toString(rateDomain.getValue(i)) + ")");
  }
  for(size_t i = 0; i < repCPU; i++)
  {
    if(verbose) ApplicationTools::displayGauge(i, repCPU - 1);

    SiteContainer * sites1 = seqSim.simulate(repRAM);
    drtl.setData(*sites1);
    drtl.initialize();
    Vdouble pr1 = drtl.getPosteriorRateOfEachSite();
    ProbabilisticSubstitutionMapping * mapping1;
    if(average)
    {
      if(joint)
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectors(drtl, nijt, false);
      else
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl, nijt, false);
    }
    else
    {
      if(joint)
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl, nijt, false);
      else
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl, nijt, false);
    }
    delete sites1;
    SiteContainer * sites2 = seqSim.simulate(repRAM);
    drtl.setData(*sites2);
    drtl.initialize();
    Vdouble pr2 = drtl.getPosteriorRateOfEachSite();
    ProbabilisticSubstitutionMapping * mapping2;
    if(average)
    {
      if(joint)
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectors(drtl, nijt, false);
      else
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl, nijt, false);
    }
    else
    {
      if(joint)
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl, nijt, false);
      else
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl, nijt, false);
    }
    delete sites2;

    for(size_t j = 0; j < repRAM; j++)
    {
      double stat = statistic.getValueForPair((*mapping1)[j], (*mapping2)[j]);
      double minR = min(pr1[j], pr2[j]);
      //int minRc = min (rc1[j], rc2[j]);
      //cout << r1[j] << "\t" << r2[j] << "\t" << minR << "\t" << rateDomain.getIndex(minR) << endl;
      //cout << rc1[j] << "\t" << rc2[j] << "\t" << minRc << "\t" << endl;
      id[rateDomain.getIndex(minR)]->addValue(stat);
      //id[minRc] -> addValue(stat);
    }
    delete mapping1;
    delete mapping2;
  }
  if(verbose) ApplicationTools::message->endLine();
  return id;
}

/******************************************************************************/

vector<IntervalData*> AnalysisTools::getNullDistributionInterDR(
  DRTreeLikelihood& drtl1,
  DRTreeLikelihood& drtl2,
  const SequenceSimulator& seqSim1,
  const SequenceSimulator& seqSim2,
  SubstitutionCount& nijt1,
  SubstitutionCount& nijt2,
  const Statistic& statistic,
  const Domain& statDomain,
  const Domain& rateDomain,
  size_t repCPU,
  size_t repRAM,
  bool average,
  bool joint,
  bool verbose)
{
  size_t nbClasses = rateDomain.getSize();
  vector<IntervalData*> id(nbClasses);
  for (size_t i = 0; i < nbClasses; i++)
  {
    id[i] = new IntervalData(statDomain, "Null distribution (classe > " + TextTools::toString(rateDomain.getValue(i)) + ")");
  }
  for (size_t i = 0; i < repCPU; i++)
  {
    if(verbose) ApplicationTools::displayGauge(i, repCPU - 1);

    SiteContainer * sites1 = seqSim1.simulate(repRAM);
    drtl1.setData(*sites1);
    drtl1.initialize();
    Vdouble pr1 = drtl1.getPosteriorRateOfEachSite();
    ProbabilisticSubstitutionMapping * mapping1;
    if(average)
    {
      if(joint)
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectors(drtl1, nijt1, false);
      else
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl1, nijt1, false);
    }
    else
    {
      if(joint)
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl1, nijt1, false);
      else
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl1, nijt1, false);
    }
    delete sites1;

    SiteContainer* sites2 = seqSim2.simulate(repRAM);
    drtl2.setData(*sites2);
    drtl2.initialize();
    Vdouble pr2 = drtl2.getPosteriorRateOfEachSite();
    ProbabilisticSubstitutionMapping* mapping2;
    if(average)
    {
      if(joint)
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectors(drtl2, nijt2, false);
      else
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl2, nijt2, false);
    }
    else
    {
      if(joint)
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl2, nijt2, false);
      else
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl2, nijt2, false);
    }
    delete sites2;
    
    for(size_t j = 0; j < repRAM; j++)
    {
      double stat = statistic.getValueForPair((*mapping1)[j], (*mapping2)[j]);
      double minR = min(pr1[j], pr2[j]);
      id[rateDomain.getIndex(minR)]->addValue(stat);
    }
    delete mapping1;
    delete mapping2;
  }
  if (verbose) ApplicationTools::message->endLine();
  return id;
}


/******************************************************************************/

void AnalysisTools::getNullDistributionIntraDR(
  DRTreeLikelihood& drtl,
  const SequenceSimulator& seqSim,
  SubstitutionCount& nijt,
  const Statistic& statistic,
  ostream* out,
  vector< vector<double> >* simstats,
  const Domain* rateDomain,
  size_t repCPU,
  size_t repRAM,
  bool average,
  bool joint,
  bool verbose)
{
  // Write header line:
  if (out)
    *out << "Stat\tRCmin\tPRmin\tNmin" << endl;
  if (simstats && rateDomain)
    if (rateDomain->getSize() != simstats->size())
      throw Exception("AnalysisTools::getNullDistributionIntraDR. Input vector should be of same size as rate domain.");
  if (simstats && !rateDomain)
    if (simstats->size() != 1)
      throw Exception("AnalysisTools::getNullDistributionIntraDR. Input vector should be of same size 1 as no rate domain was specified.");
  for (size_t i = 0; i < repCPU; i++)
  {
    if (verbose) ApplicationTools::displayGauge(i, repCPU - 1);
    
    SiteContainer* sites1 = seqSim.simulate(repRAM);
    drtl.setData(*sites1);
    drtl.initialize();
    vector<size_t> rc1 = drtl.getRateClassWithMaxPostProbOfEachSite();
    Vdouble pr1 = drtl.getPosteriorRateOfEachSite();
    
    ProbabilisticSubstitutionMapping* mapping1;
    if (average)
    {
      if (joint)
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectors(drtl, nijt, false);
      else
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl, nijt, false);
    }
    else
    {
      if (joint)
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl, nijt, false);
      else
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl, nijt, false);
    }
    delete sites1;
    Vdouble norm1 = computeNorms(*mapping1);
     
    SiteContainer* sites2 = seqSim.simulate(repRAM);
    drtl.setData(*sites2);
    drtl.initialize();
    vector<size_t> rc2 = drtl.getRateClassWithMaxPostProbOfEachSite();
    Vdouble pr2 = drtl.getPosteriorRateOfEachSite();
    
    ProbabilisticSubstitutionMapping* mapping2;
    if (average)
    {
      if (joint)
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectors(drtl, nijt, false);
      else
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl, nijt, false);
    }
    else
    {
      if (joint)
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl, nijt, false);
      else
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl, nijt, false);
    }
    delete sites2;
    Vdouble norm2 = computeNorms(*mapping2);
      
    for (size_t j = 0; j < repRAM; j++)
    {
      double stat = statistic.getValueForPair((*mapping1)[j], (*mapping2)[j]);
      double nmin = std::min(norm1[j], norm2[j]);
      if (out)
        *out << stat << "\t" << std::min(rc1[j], rc2[j]) << "\t" << std::min(pr1[j], pr2[j]) << "\t" << nmin << endl;
      if (simstats) {
        if (rateDomain) {
          try {
            size_t cat = rateDomain->getIndex(nmin);
            (*simstats)[cat].push_back(stat);
          } catch (OutOfRangeException& oore) {}
        } else {
          (*simstats)[0].push_back(stat);
        }
      }
    }
    //cout << "Freeing memory." << endl;
    delete mapping1;
    delete mapping2;
  }
  if (verbose)
    ApplicationTools::message->endLine();
}

/******************************************************************************/

void AnalysisTools::getNullDistributionInterDR(
  DRTreeLikelihood & drtl1,
  DRTreeLikelihood & drtl2,
  const SequenceSimulator & seqSim1,
  const SequenceSimulator & seqSim2,
  SubstitutionCount & nijt1,
  SubstitutionCount & nijt2,
  const Statistic & statistic,
  ostream & out,
  size_t repCPU,
  size_t repRAM,
  bool average,
  bool joint,
  bool verbose)
{
  // Write header line:
  out << "Stat\tRCmin\tPRmin\tNmin" << endl;
  for(size_t i = 0; i < repCPU; i++)
  {
    if(verbose) ApplicationTools::displayGauge(i, repCPU - 1);

    SiteContainer * sites1 = seqSim1.simulate(repRAM);
    drtl1.setData(*sites1);
    drtl1.initialize();
    vector<size_t> rc1 = drtl1.getRateClassWithMaxPostProbOfEachSite();
    Vdouble pr1 = drtl1.getPosteriorRateOfEachSite();
 
    ProbabilisticSubstitutionMapping * mapping1;
    if(average)
    {
      if(joint)
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectors(drtl1, nijt1, false);
      else
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl1, nijt1, false);
    }
    else
    {
      if(joint)
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl1, nijt1, false);
      else
        mapping1 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl1, nijt1, false);
    }
    delete sites1;
    Vdouble norm1 = computeNorms(*mapping1);

    SiteContainer * sites2 = seqSim2.simulate(repRAM);
    drtl2.setData(*sites2);
    drtl2.initialize();
    vector<size_t> rc2 = drtl2.getRateClassWithMaxPostProbOfEachSite();
    Vdouble pr2 = drtl2.getPosteriorRateOfEachSite();
    ProbabilisticSubstitutionMapping * mapping2;
    if(average)
    {
      if(joint)
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectors(drtl2, nijt2, false);
      else
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsMarginal(drtl2, nijt2, false);
    }
    else
    {
      if(joint)
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(drtl2, nijt2, false);
      else
        mapping2 = SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(drtl2, nijt2, false);
    }
    delete sites2;
    Vdouble norm2 = computeNorms(*mapping2);

    for(size_t j = 0; j < repRAM; j++)
    {
      double stat = statistic.getValueForPair((*mapping1)[j], (*mapping2)[j]);
      out << stat << "\t" << std::min(rc1[j], rc2[j]) << "\t" << std::min(pr1[j], pr2[j]) << "\t" << std::min(norm1[j], norm2[j]) << endl;
    }

    delete mapping1;
    delete mapping2;
  }
  if(verbose) ApplicationTools::message->endLine();
}

/******************************************************************************/

void AnalysisTools::getNullDistributionIntraWithoutReestimatingCounts(
  const NonHomogeneousSequenceSimulator& seqSim,
  const Statistic& statistic,
  ostream& out,
  size_t rep,
  bool verbose)
{
  // Write header line:
  out << "Stat\tr1\tr2" << endl;
  TotalSubstitutionRegister reg(seqSim.getSubstitutionModelSet()->getModel(0));
  for (size_t i = 0; i < rep; ++i)
  {
    ApplicationTools::displayGauge(i, rep - 1);
    RASiteSimulationResult* hssr1 = 
      dynamic_cast<RASiteSimulationResult*>(seqSim.dSimulateSite());
    VVdouble vector1 = hssr1->getSubstitutionVector(reg);
    double     rate1 = hssr1->getRate();
    delete hssr1;
    RASiteSimulationResult* hssr2 = 
      dynamic_cast<RASiteSimulationResult*>(seqSim.dSimulateSite());
    VVdouble vector2 = hssr2->getSubstitutionVector(reg);
    double     rate2 = hssr2->getRate(); 
    delete hssr2;

    double stat = statistic.getValueForPair(vector1, vector2);
    out << stat << "\t" << rate1 << "\t" << rate2 << endl;
  }
  ApplicationTools::message->endLine();
}

/******************************************************************************/

