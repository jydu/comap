//
// File: AnalysisTools.h
// Created by: Julien Dutheil
// Created on: Fri Nov 28 16:33:00 2003
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

#ifndef _ANALYSISTOOLS_H_
#define _ANALYSISTOOLS_H_

#include "Statistics.h"
#include "Domain.h"
#include "IntervalData.h"

#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>

// From bpp-phyl:
#include <Bpp/Phyl/Legacy/Simulation/NonHomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Legacy/Mapping/ProbabilisticSubstitutionMapping.h>
#include <Bpp/Phyl/Legacy/Likelihood/DRTreeLikelihood.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <vector>

using namespace std;

/**
 * @brief Some tools for the coevolution analysis.
 */
class AnalysisTools
{
  public:
    AnalysisTools();
    virtual ~AnalysisTools();
  
  public:
    
    /**
     * @brief Read a vector of dimension 2 from an input stream.
     *
     * @param in The input stream.
     * @return A VVdouble object.
     */
    static VVdouble getFromStream(istream& in);

    /**
     * @brief Compute the matrix of scalar products.
     *
     * The function takes as input a vector of vectors (VVdouble) and computes
     * the scalar product of each pair of vectors.
     *
     * @param vectors A vector of vectors.
     * @return The matrix of scalar products.
     */
    static VVdouble computeScalarProductMatrix(const VVdouble& vectors);

    /**
     * @brief Compute the matrix of scalar products.
     *
     * The function takes as input two vectors of vectors (VVdouble) and computes
     * the scalar product of each pair of vectors.
     * If independantComparisons is set to true, then the two vectors of vectors are
     * expected to be of the same length and vectors are compared 2 by 2.
     *
     * @param vectors1, vectors2 Two vectors of vectors.
     * @return The matrix of scalar products.
     */
    static VVdouble computeScalarProductMatrix(
      const VVdouble& vectors1,
      const VVdouble& vectors2,
      bool independantComparisons);

    /**
     * @brief Compute the matrix of cosinus.
     *
     * The function takes as input a vector of vectors (VVdouble) and computes
     * the cosinus of each pair of vectors.
     *
     * @param vectors A vector of vectors.
     * @return The matrix of cosinus.
     */
    static VVdouble computeCosinusMatrix(const VVdouble& vectors);
    
    /**
     * @brief Compute the matrix of cosinus.
     *
     * The function takes as input two vectors of vectors (VVdouble) and computes
     * the cosinus of each pair of vectors.
     * If independantComparisons is set to true, then the two vectors of vectors are
     * expected to be of the same length and vectors are compared 2 by 2.
     *
     * @param vectors1, vectors2 Two vectors of vectors.
     * @return The cosinus of scalar products.
     */
    static VVdouble computeCosinusMatrix(
      const VVdouble& vectors1,
      const VVdouble& vectors2,
      bool independantComparisons);
  
    /**
     * @brief Compute the matrix of correlations.
     *
     * The function takes as input a vector of vectors (VVdouble) and computes
     * the correlation of each pair of vectors.
     *
     * @param vectors A vector of vectors.
     * @return The matrix of correlations.
     */
    static VVdouble computeCorrelationMatrix(const VVdouble& vectors);

    /**
     * @brief Compute the matrix of correlations.
     *
     * The function takes as input two vectors of vectors (VVdouble) and computes
     * the correlation of each pair of vectors.
     * If independantComparisons is set to true, then the two vectors of vectors are
     * expected to be of the same length and vectors are compared 2 by 2.
     *
     * @param vectors1, vectors2 Two vectors of vectors.
     * @return The matrix of correlations.
     */
    static VVdouble computeCorrelationMatrix(
      const VVdouble& vectors1,
      const VVdouble& vectors2,
      bool independantComparisons);
  
    /**
     * @brief Compute the matrix of covariance.
     *
     * The function takes as input a vector of vectors (VVdouble) and computes
     * the covariance of each pair of vectors.
     *
     * @param vectors A vector of vectors.
     * @return The matrix of covariance.
     */
    static VVdouble computeCovarianceMatrix(const VVdouble& vectors);

    /**
     * @brief Compute the matrix of covariances.
     *
     * The function takes as input two vectors of vectors (VVdouble) and computes
     * the covariance product of each pair of vectors.
     * If independantComparisons is set to true, then the two vectors of vectors are
     * expected to be of the same length and vectors are compared 2 by 2.
     *
     * @param vectors1, vectors2 Two vectors of vectors.
     * @return The matrix of covariances.
     */
    static VVdouble computeCovarianceMatrix(
      const VVdouble& vectors1,
      const VVdouble& vectors2,
      bool independantComparisons);
  
    /**
     * @brief Get the norms of each substitution vector.
     *
     * @param mapping A substitution mapping.
     * @return A vector containing the norms of each vector in the list.
     */
    static Vdouble computeNorms(const LegacyProbabilisticSubstitutionMapping& mapping);
  
    static void writeMatrix(
      const VVdouble& matrix,
      const SiteContainerInterface& sites,
      ostream& out);
      
    static void writeMatrix(
      const VVdouble& matrix,
      const SiteContainerInterface& sites1,
      const SiteContainerInterface& sites2,
      ostream& out);
      
    /**********************************************************************/
    
    static vector<IntervalData*> getNullDistributionIntraDR(
      shared_ptr<DRTreeLikelihoodInterface> drtl,
      const SequenceSimulatorInterface& seqSim,
      shared_ptr<SubstitutionCountInterface> nijt,
      const Statistic& statistic,
      const Domain& statDomain,
      const Domain& rateDomain,
      size_t repCPU,
      size_t repRAM,
      bool average,
      bool joint,
      bool verbose = true);

    static vector<IntervalData*> getNullDistributionInterDR(
      shared_ptr<DRTreeLikelihoodInterface> drtl1,
      shared_ptr<DRTreeLikelihoodInterface> drtl2,
      const SequenceSimulatorInterface& seqSim1,
      const SequenceSimulatorInterface& seqSim2,
      shared_ptr<SubstitutionCountInterface> nijt1,
      shared_ptr<SubstitutionCountInterface> nijt2,
      const Statistic& statistic,
      const Domain& statDomain,
      const Domain& rateDomain,
      size_t repCPU,
      size_t repRAM,
      bool average,
      bool joint,
      bool verbose = true);

    /**
     * @param out If not 0, where to output simulation results.
     * @param simstats If not 0, where to store simulation results (statistic only).
     * @param rateDomain A domain object describing rate categories in order to perform conditional p-value computations.
     *                   If 0 or with one rate class, no rate conditionning will be done.
     */
    static void getNullDistributionIntraDR(
      shared_ptr<DRTreeLikelihoodInterface> drtl,
      const SequenceSimulatorInterface& seqSim,
      shared_ptr<SubstitutionCountInterface> nijt,
      const Statistic& statistic,
      ostream* out,
      VVdouble* simstats,
      const Domain* rateDomain,
      size_t repCPU,
      size_t repRAM,
      bool average,
      bool joint,
      bool verbose = true);
      
    static void getNullDistributionInterDR(
      shared_ptr<DRTreeLikelihoodInterface> drtl1,
      shared_ptr<DRTreeLikelihoodInterface> drtl2,
      const SequenceSimulatorInterface& seqSim1,
      const SequenceSimulatorInterface& seqSim2,
      shared_ptr<SubstitutionCountInterface> nijt1,
      shared_ptr<SubstitutionCountInterface> nijt2,
      const Statistic& statistic,
      ostream& out,
      size_t repCPU,
      size_t repRAM,
      bool average,
      bool joint,
      bool verbose = true);

    static void getNullDistributionIntraWithoutReestimatingCounts(
      const NonHomogeneousSequenceSimulator& seqSim,
      const Statistic& statistic,
      ostream& out,
      size_t rep,
      bool verbose = true);
};


#endif  //_ANALYSISTOOLS_H_

