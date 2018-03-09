#include "qcf.hpp"
#include "SeqData.hpp"
#include "Quartet.hpp"
#include "Bootstrap.hpp"

Bootstrap::Bootstrap(const int nReps){
  reps = nReps;
  r = new MbRandom;
}

std::vector<double> Bootstrap::operator()(Quartet& qrt){
  std::vector<double> weights(3, 0.0);
  std::vector<double> weightSums(3,0.0);
  double overallSum;
  std::vector<int> ix(qrt.seqPtr->nSites, 0);
  std::vector<double> bootResults(3, 0.0);
  for(int i = 0; i < reps; ++i){
    randomVector_(0, qrt.seqPtr->nSites - 1, ix);
    bootResults = qrt.eval2(ix);
    weightSums[0] += bootResults[0];
    weightSums[1] += bootResults[1];
    weightSums[2] += bootResults[2];
  }
  overallSum = weightSums[0] + weightSums[1] + weightSums[2];
  weights[0] = weightSums[0] / overallSum;
  weights[1] = weightSums[1] / overallSum;
  weights[2] = weightSums[2] / overallSum;
  return weights;
}

void Bootstrap::randomVector_(const int low, const int high, std::vector<int>& vec){
  for(size_t i = 0; i < vec.size(); ++i){
    vec[i] = r->sampleInteger(low,high);
  }
}
