#include <assert.h>
#include <math.h>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "qcf.hpp"
#include "SeqData.hpp"
#include "Quartet.hpp"

Quartet::Quartet(const std::string i, const std::string j,
                 const std::string k, const std::string l,
                 SeqData* seq){
  A = seq->seqIndex[i], B = seq->seqIndex[j],
  C = seq->seqIndex[k], D = seq->seqIndex[l];
  hapA = i, hapB = j,
  hapC = k, hapD = l;
  /*std::cerr << hapA << ":" << A << "\t"
            << hapB << ":" << B << "\t"
            << hapC << ":" << C << "\t"
            << hapD << ":" << D << std::endl;*/
  seqPtr = seq;
  makeIndexVec_();
}

void Quartet::getCountMatrices_(std::vector<int>& ix){
  assert((int) ix.size() == seqPtr->nSites);
  for(std::vector<int>::iterator s = ix.begin();
      s != ix.end(); ++s){
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
      resolved = resolveAmbiguity_(seqPtr->dna[A][*s],
                                   seqPtr->dna[B][*s],
                                   seqPtr->dna[C][*s],
                                   seqPtr->dna[D][*s]);
    }
  }
}

bool Quartet::resolveAmbiguity_(const int one,
                                const int two,
                                const int three,
                                const int four){
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
    for(size_t i = 0; i < baseLookup[one].size(); ++i){
      for(size_t j = 0; j < baseLookup[two].size(); ++j){
        for(size_t k = 0; k < baseLookup[three].size(); ++k){
          for(size_t r = 0; r < baseLookup[four].size(); ++r){
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

double Quartet::resolveAmbiguity3_(const int a, const int b, double mat[4][4]){
  double denom = 0.0;
  if(a == 4 || b == 4){ // Don't allow gaps
    return 0.0;
  } else {
    denom = baseLookup[a].size() * baseLookup[b].size();
    for(size_t i = 0; baseLookup[a].size(); ++i){
      for(size_t j = 0; baseLookup[b].size(); ++j){
        mat[baseLookup[a][i]][baseLookup[b][j]] += 1.0 / denom;
      }
    }
  }
  return 1.0;
}

std::vector<double> Quartet::eval(std::vector<int>& vec){
  std::vector<double> scores(3,0.0);
  std::vector<double> weights(3,0.0);
  if(vec.size() == 0){
    getCountMatrices_(index);
  } else {
    getCountMatrices_(vec);
  }
  for(int i = 0; i < 4; ++i){
    for(int j = 0; j < 4; ++j){
      for(int k = 0; k < 4; ++k){
        for(int l = 0; l < 4; ++l){
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
  weights = getWeights_(scores);
  return weights;
}

std::vector<double> Quartet::eval2(std::vector<int>& vec){
  /*
    This metric is essentially a parsimony score that compares sequences
    based on how similar they are, and gets a tree topology that combines
    similarity and penalizes for dissimilarity. I don't think that there
    is a technical proof that this works, but in practice I've found that
    it performs pretty well for closely related taxa. The GTR simulations
    show that it works well too, but the LogDet transformation is more
    general and so that's what we use (eval3).
  */
  std::vector<double> scores(3,0.0);
  std::vector<double> weights(3,0.0);
  std::vector<int> ix;
  double mAB = 0.0, mAC = 0.0, mAD = 0.0, mBC = 0.0, mBD = 0.0, mCD = 0.0;
  if(vec.size() == 0){
    ix = index;
  } else {
    ix = vec;
  }
  assert((int) ix.size() == seqPtr->nSites);
  for(std::vector<int>::iterator s = ix.begin();
      s != ix.end(); ++s){
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
      mAB += resolveAmbiguity2_(seqPtr->dna[A][*s], seqPtr->dna[B][*s]);
      mAC += resolveAmbiguity2_(seqPtr->dna[A][*s], seqPtr->dna[C][*s]);
      mAD += resolveAmbiguity2_(seqPtr->dna[A][*s], seqPtr->dna[D][*s]);
      mBC += resolveAmbiguity2_(seqPtr->dna[B][*s], seqPtr->dna[C][*s]);
      mBD += resolveAmbiguity2_(seqPtr->dna[B][*s], seqPtr->dna[D][*s]);
      mCD += resolveAmbiguity2_(seqPtr->dna[C][*s], seqPtr->dna[D][*s]);
    }
  }
  scores[0] = mAB + mCD - mAC - mAD - mBC - mBD;
  scores[1] = mAC + mBD - mAB - mAD - mBC - mCD;
  scores[2] = mAD + mBC - mAB - mAC - mBD - mCD;
  weights = getWeights_(scores);
  return weights;
}

std::vector<double> Quartet::eval3(std::vector<int>& vec){
  /*

  */
  std::vector<double> scores(3,0.0);
  std::vector< std::vector<double> > weights {{1.0,0.0,0.0},
                                              {0.0,1.0,0.0},
                                              {0.0,0.0,1.0},
                                              {0.0,0.0,1.0},
                                              {0.0,1.0,0.0},
                                              {1.0,0.0,0.0}};
  std::vector<int> ix;

  // Matrices for storing pairwise sequence differences
  double fAB[4][4] = {}, fAC[4][4] = {},
         fAD[4][4] = {}, fBC[4][4] = {},
         fBD[4][4] = {}, fCD[4][4] = {};

  // Variables for LogDet scores
  double dAB = 0.0, dAC = 0.0, dAD = 0.0, dBC = 0.0, dBD = 0.0, dCD = 0.0;

  // Variables for number of sites used in each pairwise comparison
  // Used for normalizing the F matrices
  double nAB = 0.0, nAC = 0.0, nAD = 0.0, nBC = 0.0, nBD = 0.0, nCD = 0.0;;

  // Variables for neighbor joining scores to determine quartet topology
  double njAB = 0.0, njAC = 0.0, njAD = 0.0, njBC = 0.0, njBD = 0.0, njCD = 0.0;

  // We do this here so that we can pass in a vector with randomized indices
  // for bootstrapping. An empty vector is passed in otherwise, and we just use
  // the internally stored index.
  if(vec.size() == 0){
    ix = index;
  } else {
    ix = vec;
  }
  assert((int) ix.size() == seqPtr->nSites);

  for(std::vector<int>::iterator s = ix.begin(); s != ix.end(); ++s){
    if(seqPtr->dna[A][*s] < 4 &&
       seqPtr->dna[B][*s] < 4 &&
       seqPtr->dna[C][*s] < 4 &&
       seqPtr->dna[D][*s] < 4){
      fAB[seqPtr->dna[A][*s]][seqPtr->dna[B][*s]] += 1.0;
      fAC[seqPtr->dna[A][*s]][seqPtr->dna[C][*s]] += 1.0;
      fAD[seqPtr->dna[A][*s]][seqPtr->dna[D][*s]] += 1.0;
      fBC[seqPtr->dna[B][*s]][seqPtr->dna[C][*s]] += 1.0;
      fBD[seqPtr->dna[B][*s]][seqPtr->dna[D][*s]] += 1.0;
      fCD[seqPtr->dna[C][*s]][seqPtr->dna[D][*s]] += 1.0;
      nAB += 1.0;
      nAC += 1.0;
      nAD += 1.0;
      nBC += 1.0;
      nBD += 1.0;
      nCD += 1.0;
    } else if(!ignore_amb_sites){
      nAB += resolveAmbiguity3_(seqPtr->dna[A][*s], seqPtr->dna[B][*s], fAB);
      nAC += resolveAmbiguity3_(seqPtr->dna[A][*s], seqPtr->dna[C][*s], fAC);
      nAD += resolveAmbiguity3_(seqPtr->dna[A][*s], seqPtr->dna[D][*s], fAD);
      nBC += resolveAmbiguity3_(seqPtr->dna[B][*s], seqPtr->dna[C][*s], fBC);
      nBD += resolveAmbiguity3_(seqPtr->dna[B][*s], seqPtr->dna[D][*s], fBD);
      nCD += resolveAmbiguity3_(seqPtr->dna[C][*s], seqPtr->dna[D][*s], fCD);
    }
  }

  dAB = LogDet_(fAB, nAB);
  dAC = LogDet_(fAC, nAC);
  dAD = LogDet_(fAD, nAD);
  dBC = LogDet_(fBC, nBC);
  dBD = LogDet_(fBD, nBD);
  dCD = LogDet_(fCD, nCD);

  njAB = (2 * dAB) - (dAB + dAC + dAD) - (dAB + dBC + dBD);
  njAC = (2 * dAC) - (dAB + dAC + dAD) - (dAC + dBC + dCD);
  njAD = (2 * dAD) - (dAB + dAC + dAD) - (dAD + dBD + dCD);
  njBC = (2 * dBC) - (dAB + dBC + dBD) - (dAC + dBC + dCD);
  njBD = (2 * dBD) - (dAB + dBC + dBD) - (dAD + dBD + dCD);
  njCD = (2 * dCD) - (dAC + dBC + dCD) - (dAD + dBD + dCD);
  std::vector<double> njVec = {njAB, njAC, njAD, njBC, njBD, njCD};

  double minVal = njVec[0];
  int minIndex = 0;
  std::vector<int> minIndices = {0};
  for(int i = 1; i < 6; ++i){
    if(njVec[i] < minVal){
      minIndex = i;
      minIndices = {i};
      minVal = njVec[i];
    } else if(njVec[i] == minVal || abs(njVec[i] - minVal) < 1.0e-4){
      minIndices.push_back(i);
    }
  }

  // Now we get the weights
  //scores[0] = 0.5 * (dAC + dBD + dAB + dCD);
  //scores[1] = 0.5 * ();
  //scores[0] = njAB;// + njCD;
  //scores[1] = njAC;// + njBD;
  //scores[2] = njAD;// + njBC;
  //weights = getWeights2_(scores);
  std::vector<double> combinedWeights = {0.0,0.0,0.0};
  if(minIndices.size() > 1){
    for(size_t i = 0; i < minIndices.size(); ++i){
      for(int j = 0; j < 3; ++j){
        combinedWeights[j] += weights[minIndices[i]][j];
      }
    }
    for(int k = 0; k < 3; ++k){
      combinedWeights[k] /= (double) minIndices.size();
    }
    return combinedWeights;
  } else {
    return weights[minIndex];
  }
}

double Quartet::LogDet_(double mat[4][4], double norm){
  double normMat[4][4] = {}, res = 0.0;
  // Normalize all of the entries first
  for(int i = 0; i < 4; ++i){
    for(int j = 0; j < 4; j++){
      normMat[i][j] = mat[i][j] / norm;
    }
  }
  res = determinant_(normMat);
  return -1.0 * std::log(res);
}

std::vector<double> Quartet::getWeights_(std::vector<double>& vec){
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
    for(int i = 0; i < 3; ++i){
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

std::vector<double> Quartet::getWeights2_(std::vector<double>& vec){
  std::vector<double> res(3,0.0);
  assert(vec.size() == 3);
  double minVal = vec[0];
  int minIndex  = 0;
  if(vec[0] == vec[1] && vec[0] == vec[2]){
    res = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
  } else if(vec[0] == vec[1] && vec[0] < vec[2]){
    res = {0.5, 0.5, 0.0};
  } else if(vec[0] == vec[2] && vec[0] < vec[1]){
    res = {0.5, 0.0, 0.5};
  } else if(vec[1] == vec[2] && vec[1] < vec[0]){
    res = {0.0, 0.5, 0.5};
  } else {
    for(int i = 0; i < 3; ++i){
      if(vec[i] < minVal){
        minIndex = i;
        minVal   = vec[i];
      }
    }
    res[minIndex] += 1.0;
  }
  /*std::cerr << seqPtr->qcfPtr->hap2tax[hapA] << ","
            << seqPtr->qcfPtr->hap2tax[hapB] << ","
            << seqPtr->qcfPtr->hap2tax[hapC] << ","
            << seqPtr->qcfPtr->hap2tax[hapD] << ":\t";
  std::cerr << "0: " << vec[0] << "  1: " << vec[1] << "  2: " << vec[2] << "  max: "
            << maxIndex << "," << res[0] << "," << res[1] << ","<< res[2] << std::endl;*/
  return res;
}

void Quartet::makeIndexVec_(){
  for(int s = 0; s < seqPtr->nSites; ++s){
    index.push_back(s);
  }
}
