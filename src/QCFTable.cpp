#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>

#include "qcf.hpp"
#include "QCFData.hpp"
#include "QCFTable.hpp"

QCFTable::QCFTable(QCFData* qcf){
  qcfPtr = qcf;
  nTaxa  = qcfPtr->nTaxa;
  nQrts  = nCk(nTaxa, four);
  std::vector< std::vector<double> > tmp;
  for(uint i = 0; i < nQrts; i++){
    values.push_back(tmp);
  }
  uint idx = 0;
  for(uint i = 0; i < nTaxa - 3; i++){
    for(uint j = i+1; j < nTaxa - 2; j++){
      for(uint k = j+1; k < nTaxa - 1; k++){
        for(uint l = k+1; l < nTaxa; l++){
          QCFEntry lkp(i, j, k, l, idx);
          lookup.push_back(lkp);
          idx++;
        }
      }
    }
  }
}

void QCFTable::addValue(uint i, uint j, uint k, uint l, std::vector<double> val){
  uint tmpIndex = findIndex(i,j,k,l);
  //std::cout << i << "\t" << j << "\t" << k << "\t" << l << std::endl;
  //std::cout << tmpIndex << std::endl;
  values[tmpIndex].push_back(val);
}

uint QCFTable::findIndex(uint i, uint j, uint k, uint l){
  std::vector<uint> tmp = {i,j,k,l};
  uint res = -1;
  for(uint i = 0; i < lookup.size(); i++){
    if(lookup[i].key == tmp){
      res = lookup[i].index;
      break;
    }
  }
  if(res == -1){
    std::cout << "Shit!" << std::endl;
  }
  return res;
}

void QCFTable::write(std::string pfx){
  std::cout << "Did we get here?" << std::endl;
  std::ofstream qcfStream;
  qcfStream.open(pfx+"-qcf.txt", std::ios::out | std::ios::app);
  if(qcfStream.is_open()){
    qcfStream << "\"t1\",\"t2\",\"t3\",\"t4\",\"CF12_34\",\"CF13_24\",\"CF14_23\"" << std::endl;
  } else {
    std::cerr << "ERROR: Could not open outfile: " << pfx << "-qcf.txt.\n" << std::endl;
    exit(EXIT_FAILURE);
  }
  uint tmpIndex;
  for(uint i = 0; i < nTaxa - 3; i++){
    for(uint j = i+1; j < nTaxa - 2; j++){
      for(uint k = j+1; k < nTaxa - 1; k++){
        for(uint l = k+1; l < nTaxa; l++){
          tmpIndex = findIndex(i,j,k,l);
          double cf12_34 = 0.0, cf13_24 = 0.0, cf14_23 = 0.0, cfSum = 0.0;
          for(uint v = 0; v < values[tmpIndex].size(); v++){
            cf12_34 += values[tmpIndex][v][0];
            cf13_24 += values[tmpIndex][v][1];
            cf14_23 += values[tmpIndex][v][2];
          }
          cfSum = cf12_34 + cf13_24 + cf14_23;
          cf12_34 /= cfSum;
          cf13_24 /= cfSum;
          cf14_23 /= cfSum;
          qcfStream << "\"" << qcfPtr->taxa[i] << "\",\"" << qcfPtr->taxa[j] << "\",\""
                    << qcfPtr->taxa[k] << "\",\"" << qcfPtr->taxa[l] << "\","
                    << cf12_34 << "," << cf13_24 << "," << cf14_23 << std::endl;
        }
      }
    }
  }
}
