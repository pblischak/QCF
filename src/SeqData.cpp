#include <algorithm>
#include <cstring>
#include <fstream>
#include <functional>
#include <stdlib.h>
#include <utility>

#include "qcf.hpp"
#include "SeqData.hpp"
#include "QCFData.hpp"
#include "Quartet.hpp"

void SeqData::readPhylip(){
  std::ifstream phyStream(file_name);
  std::string val1, val2;
  uint row = 0;
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
        return; // we'll just skip this locus
      } else {
        haps.push_back(val1);
        seqIndex.insert({val1, row});
      }
      for(uint s = 0; s < val2.length(); s++){
        dna[row][s] = convert(val2[s]);
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
  std::vector<uint> ordered;
  std::vector<Quartet> qrts;
  for(uint i = 0; i < H - 3; i++){
    for(uint j = i + 1; j < H - 2; j++){
      for(uint k = j + 1; k < H - 1; k++){
        for(uint l = k + 1; l < H; l++){
          if(qcfPtr->hap2tax[haps[i]] == qcfPtr->hap2tax[haps[j]] ||
             qcfPtr->hap2tax[haps[i]] == qcfPtr->hap2tax[haps[k]] ||
             qcfPtr->hap2tax[haps[i]] == qcfPtr->hap2tax[haps[l]] ||
             qcfPtr->hap2tax[haps[j]] == qcfPtr->hap2tax[haps[k]] ||
             qcfPtr->hap2tax[haps[j]] == qcfPtr->hap2tax[haps[l]] ||
             qcfPtr->hap2tax[haps[k]] == qcfPtr->hap2tax[haps[l]]){
            continue;
          } else {
            // order haplotypes by their occurence in the taxon map
            ordered = orderHaps(i,j,k,l);
            Quartet qrt(haps[ordered[0]],
                        haps[ordered[1]],
                        haps[ordered[2]],
                        haps[ordered[3]], this);
            qrts.push_back(qrt);
          }
        }
      }
    }
  }
  return qrts;
}

std::vector<uint> SeqData::orderHaps(uint i, uint j, uint k, uint l){
  uint one   = qcfPtr->hap2tax[haps[i]],
       two   = qcfPtr->hap2tax[haps[j]],
       three = qcfPtr->hap2tax[haps[k]],
       four  = qcfPtr->hap2tax[haps[l]];
  std::vector<uint> res(4);
  std::vector< std::pair<uint, uint> > pairs;
  pairs.push_back(std::make_pair(one,i));
  pairs.push_back(std::make_pair(two,j));
  pairs.push_back(std::make_pair(three,k));
  pairs.push_back(std::make_pair(four,l));
  std::sort(pairs.begin(), pairs.end());
  /*std::cerr << "{" << pairs[0].first << "," << haps[pairs[0].second] << " : "
                   << pairs[1].first << "," << haps[pairs[1].second] << " : "
                   << pairs[2].first << "," << haps[pairs[2].second] << " : "
                   << pairs[3].first << "," << haps[pairs[3].second] << "}" << std::endl;*/
  res = {pairs[0].second,
         pairs[1].second,
         pairs[2].second,
         pairs[3].second};
  return res;
}
