#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <unordered_map>

#include "SeqData.hpp"
#include "QCFData.hpp"
#include "Quartet.hpp"
#include "qcf.hpp"

void SeqData::_readPhylip(){
  std::ifstream phyStream(file_name);
  std::string val1, val2;
  int row = 0;
  if(phyStream.is_open()){
    phyStream >> val1 >> val2;
    nSeqs = atoi(val1.c_str());
    nSites = atoi(val2.c_str());
    dna.resize(nSeqs);
    for(uint i = 0; i < nSeqs; i++){
      dna[i].resize(nSites);
    }
    while(phyStream >> val1 >> val2){
      if(qcfPtr->hap2tax.count(val1) == 0){
        std::cerr << "\nThe haplotype name \'" << val1 << "\' is not present in the given map file." << std::endl;
        std::cerr << "Skipping this gene (" << file_name << ")...\n" << std::endl;
        skip = 1;
      } else {
        haps.push_back(val1);
      }
      for(uint s = 0; s < val2.length(); s++){
        dna[row][s] = _convert(val2[s]);
      }
      row++;
    }
  } else {
    std::cerr << "Problem reading in file " << file_name << ". Skipping..." << std::endl;
    skip = 1;
  }
}

std::vector<Quartet> SeqData::get_quartets(){
  uint H = haps.size();
  std::vector<Quartet> qrts;
  for(uint i = 0; i < H - 3; i++){
    for(uint j = i + 1; j < H - 2; j++){
      for(uint k = j + 1; k < H - 1; k++){
        for(uint l = k + 1; l < H; l++){
          Quartet qrt(i, j, k, l, this);
          qrts.push_back(qrt);
        }
      }
    }
  }
  return qrts;
}
