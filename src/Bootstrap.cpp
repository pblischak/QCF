#include <vector>
#include <string>
#include <iostream>

#include "SeqData.hpp"
#include "Quartet.hpp"
#include "Bootstrap.hpp"
#include "qcf.hpp"

std::vector<double> Bootstrap::operator()(Quartet &qrt){
  std::vector<double> weights(3, 0.0);
  std::vector< std::vector<double> > boot_results(reps, std::vector<double>(3, 0.0));
  for(uint r = 0; r < reps; r++){

  }

  return weights;
}
