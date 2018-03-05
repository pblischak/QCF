#include "qcf.hpp"
#include "SeqData.hpp"
#include "Quartet.hpp"
#include "MbRandom.hpp"
#include "Bootstrap.hpp"


Bootstrap::Bootstrap(int &nreps){
  reps = nreps;
  r = new MbRandom;
}

std::vector<double> Bootstrap::operator()(Quartet &qrt){
  std::vector<double> weights(3, 0.0);
  std::vector<double> weightSums(3,0.0);
  double overallSum;
  std::vector<int> ix(qrt.seqPtr->nSites, 0);
  std::vector< std::vector<double> > boot_results(reps, std::vector<double>(3, 0.0));
  for(uint i = 0; i < reps; i++){
    randomVector(0, qrt.seqPtr->nSites - 1, ix);
    boot_results[i] = qrt.eval2(ix);
    weightSums[0] += boot_results[i][0];
    weightSums[1] += boot_results[i][1];
    weightSums[2] += boot_results[i][2];
  }
  overallSum = weightSums[0] + weightSums[1] + weightSums[2];
  weights[0] = weightSums[0] / overallSum;
  weights[1] = weightSums[1] / overallSum;
  weights[2] = weightSums[2] / overallSum;
  return weights;
}

void Bootstrap::randomVector(int low, int high, std::vector<int> &vec){
  for(uint i = 0; i < vec.size(); i++){
    vec[i] = r->sampleInteger(low,high);
  }
}
