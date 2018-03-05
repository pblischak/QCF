#include <assert.h>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <vector>

#include "QCFData.hpp"
#include "SeqData.hpp"
#include "Quartet.hpp"
#include "qcf.hpp"

Quartet::Quartet(std::string i, std::string j, std::string k, std::string l, SeqData* seq){
  A = seq->seqIndex[i], B = seq->seqIndex[j],
  C = seq->seqIndex[k], D = seq->seqIndex[l];
  hapA = i, hapB = j,
  hapC = k, hapD = l;
  /*std::cerr << hapA << ":" << A << "\t"
            << hapB << ":" << B << "\t"
            << hapC << ":" << C << "\t"
            << hapD << ":" << D << std::endl;*/
  seqPtr = seq;
  makeIndexVec();
}

void Quartet::getCountMatrices(std::vector<int> &ix){
  assert(ix.size() == seqPtr->nSites);
  for(std::vector<int>::iterator s = ix.begin();
      s != ix.end(); s++){
    if(seqPtr->dna[A][*s] < 4 &&
       seqPtr->dna[B][*s] < 4 &&
       seqPtr->dna[C][*s] < 4 &&
       seqPtr->dna[D][*s] < 4){
      ABCD[seqPtr->dna[A][*s] * 4 + seqPtr->dna[B][*s]]
          [seqPtr->dna[C][*s] * 4 + seqPtr->dna[D][*s]] += 1.0;
      ACBD[seqPtr->dna[A][*s] * 4 + seqPtr->dna[C][*s]]
          [seqPtr->dna[B][*s] * 4 + seqPtr->dna[D][*s]] += 1.0;
      ADBC[seqPtr->dna[A][*s] * 4 + seqPtr->dna[D][*s]]
          [seqPtr->dna[B][*s] * 4 + seqPtr->dna[C][*s]] += 1.0;
    } else if(!ignore_amb_sites){
      resolved = resolveAmbiguity(seqPtr->dna[A][*s],
                                  seqPtr->dna[B][*s],
                                  seqPtr->dna[C][*s],
                                  seqPtr->dna[D][*s]);
    }
  }
}

bool Quartet::resolveAmbiguity(const int& one,
                               const int& two,
                               const int& three,
                               const int& four){
  /* Check to see if any combination of three of the four taxa are ambiguous. */
  double denom = 0.0;
  if((one >= 4 && two   >= 4 && three >= 4) ||
     (one >= 4 && two   >= 4 && four  >= 4) ||
     (one >= 4 && three >= 4 && four  >= 4) ||
     (two >= 4 && three >= 4 && four  >= 4)){
    return 0;
  } else if(one == 4 || two == 4 || three == 4 || four == 4){ /* Don't allow gaps. */
    return 0;
  } else {
    denom = baseLookup[one].size()   * baseLookup[two].size()
          * baseLookup[three].size() * baseLookup[four].size();
    for(unsigned i = 0; i < baseLookup[one].size(); i++){
      for(unsigned j = 0; j < baseLookup[two].size(); j++){
        for(unsigned k = 0; k < baseLookup[three].size(); k++){
          for(unsigned r = 0; r < baseLookup[four].size(); r++){
            ABCD[baseLookup[one][i]   * 4 + baseLookup[two][j]]
                [baseLookup[three][k] * 4 + baseLookup[four][r]]  += 1.0 / denom;
            ACBD[baseLookup[one][i]   * 4 + baseLookup[three][j]]
                [baseLookup[two][k]   * 4 + baseLookup[four][r]]  += 1.0 / denom;
            ADBC[baseLookup[one][i]   * 4 + baseLookup[four][j]]
                [baseLookup[two][k]   * 4 + baseLookup[three][r]] += 1.0 / denom;
          }
        }
      }
    }
  }
  return 1;
}

std::vector<double> Quartet::eval(std::vector<int> &vec){
  std::vector<double> scores(3,0.0);
  std::vector<double> weights(3,0.0);
  if(vec.size() == 0){
    getCountMatrices(index);
  } else {
    getCountMatrices(vec);
  }
  for(uint i = 0; i < 4; i++){
    for(uint j = 0; j < 4; j++){
      for(uint k = 0; k < 4; k++){
        for(uint l = 0; l < 4; l++){
          /*if(i == j && k != l){
            scores[0] += ABCD[i * 4 + j][k * 4 + l];
            scores[1] += ACBD[i * 4 + j][k * 4 + l];
            scores[2] += ADBC[i * 4 + j][k * 4 + l];
          }*/
          if(ABCD[i * 4 + j][k * 4 + l] == ABCD[i * 4 + j][l * 4 + k] &&
             ABCD[i * 4 + j][k * 4 + l] != ABCD[j * 4 + i][k * 4 + l]){
            scores[0] += 1.0;
          }
          if(ACBD[i * 4 + j][k * 4 + l] == ACBD[i * 4 + j][l * 4 + k] &&
             ACBD[i * 4 + j][k * 4 + l] != ACBD[j * 4 + i][k * 4 + l]){
            scores[1] += 1.0;
          }
          if(ADBC[i * 4 + j][k * 4 + l] == ADBC[i * 4 + j][l * 4 + k] &&
             ADBC[i * 4 + j][k * 4 + l] != ADBC[j * 4 + i][k * 4 + l]){
            scores[2] += 1.0;
          }
        }
      }
    }
  }
  weights = getWeights(scores);
  return weights;
}

std::vector<double> Quartet::eval2(std::vector<int> &vec){
  std::vector<double> scores(3,0.0);
  std::vector<double> weights(3,0.0);
  std::vector<int> ix;
  double mAB = 0.0, mAC = 0.0, mAD = 0.0, mBC = 0.0, mBD = 0.0, mCD = 0.0;
  if(vec.size() == 0){
    ix = index;
  } else {
    ix = vec;
  }
  assert(ix.size() == seqPtr->nSites);
  for(std::vector<int>::iterator s = ix.begin();
      s != ix.end(); s++){
    if(seqPtr->dna[A][*s] < 4 &&
       seqPtr->dna[B][*s] < 4 &&
       seqPtr->dna[C][*s] < 4 &&
       seqPtr->dna[D][*s] < 4){
      if(seqPtr->dna[A][*s] == seqPtr->dna[B][*s]) {mAB += 1.0;}
      if(seqPtr->dna[A][*s] == seqPtr->dna[C][*s]) {mAC += 1.0;}
      if(seqPtr->dna[A][*s] == seqPtr->dna[D][*s]) {mAD += 1.0;}
      if(seqPtr->dna[B][*s] == seqPtr->dna[C][*s]) {mBC += 1.0;}
      if(seqPtr->dna[B][*s] == seqPtr->dna[D][*s]) {mBD += 1.0;}
      if(seqPtr->dna[C][*s] == seqPtr->dna[D][*s]) {mCD += 1.0;}
    } else if(!ignore_amb_sites){
      mAB += resolveAmbiguity2(seqPtr->dna[A][*s], seqPtr->dna[B][*s]);
      mAC += resolveAmbiguity2(seqPtr->dna[A][*s], seqPtr->dna[C][*s]);
      mAD += resolveAmbiguity2(seqPtr->dna[A][*s], seqPtr->dna[D][*s]);
      mBC += resolveAmbiguity2(seqPtr->dna[B][*s], seqPtr->dna[C][*s]);
      mBD += resolveAmbiguity2(seqPtr->dna[B][*s], seqPtr->dna[D][*s]);
      mCD += resolveAmbiguity2(seqPtr->dna[C][*s], seqPtr->dna[D][*s]);
    }
  }
  scores[0] = mAB + mCD - mAC - mAD - mBC - mBD;
  scores[1] = mAC + mBD - mAB - mAD - mBC - mCD;
  scores[2] = mAD + mBC - mAB - mAC - mBD - mCD;
  weights = getWeights(scores);
  return weights;
}

std::vector<double> Quartet::getWeights(std::vector<double> &vec){
  std::vector<double> res(3,0.0);
  assert(vec.size() == 3);
  double maxVal = vec[0];
  int maxIndex  = 0;
  if(vec[0] == vec[1] && vec[0] == vec[2]){
    res = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
  } else if(vec[0] == vec[1] && vec[0] > vec[2]){
    res = {0.5, 0.5, 0.0};
  } else if(vec[0] == vec[2] && vec[0] > vec[1]){
    res = {0.5, 0.0, 0.5};
  } else if(vec[1] == vec[2] && vec[1] > vec[0]){
    res = {0.0, 0.5, 0.5};
  } else {
    for(uint i = 0; i < 3; i++){
      if(vec[i] > maxVal){
        maxIndex = i;
        maxVal   = vec[i];
      }
    }
    res[maxIndex] += 1.0;
  }
  /*std::cerr << seqPtr->qcfPtr->hap2tax[hapA] << ","
            << seqPtr->qcfPtr->hap2tax[hapB] << ","
            << seqPtr->qcfPtr->hap2tax[hapC] << ","
            << seqPtr->qcfPtr->hap2tax[hapD] << ":\t";
  std::cerr << "0: " << vec[0] << "  1: " << vec[1] << "  2: " << vec[2] << "  max: "
            << maxIndex << "," << res[0] << "," << res[1] << ","<< res[2] << std::endl;*/
  return res;
}

void Quartet::makeIndexVec(){
  for(uint s = 0; s < seqPtr->nSites; s++){
    index.push_back(s);
  }
}
