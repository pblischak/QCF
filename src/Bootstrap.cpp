#include <vector>
#include <string>
#include <iostream>

#include "SeqData.hpp"
#include "Quartet.hpp"
#include "MbRandom.hpp"
#include "Bootstrap.hpp"
#include "qcf.hpp"

Bootstrap::Bootstrap(int &nreps){
  reps = nreps;
  r = new MbRandom;
}

Bootstrap::~Bootstrap(){delete r;}

std::vector<double> Bootstrap::operator()(Quartet &qrt){
  std::vector<double> weights(3, 0.0);
  std::vector< std::vector<double> > boot_results(reps, std::vector<double>(3, 0.0));
  double blah = 0;
  for(uint i = 0; i < reps; i++){
    blah += i;
  }
  return weights;
}

void Bootstrap::randomVector(int low, int high, std::vector<int> &vec){
  for(uint i = 0; i < vec.size(); i++){
    vec[i] = r->sampleInteger(low,high);
  }
}
