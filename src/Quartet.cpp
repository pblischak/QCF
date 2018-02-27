

void Quartet::getCountMatrix(const int& one,
                             const int& two,
                             const int& three,
                             const int& four,
                             double cp[16][16]){
  for(int s = 0; s < seqPtr->nSites; s++){
    if(seqPtr->dna[one][s]   < 4 &&
       seqPtr->dna[two][s]   < 4 &&
       seqPtr->dna[three][s] < 4 &&
       seqPtr->dna[four][s]  < 4){
      cp[seqPtr->dna[one][s]   * 4 + seqPtr->dna[two][s]]
        [seqPtr->dna[three][s] * 4 + seqPtr->dna[four][s]] += 1.0;
    } else if(!ignore_amb_sites){
      resolved = _resolveAmbiguity(_dnaMatrix[_taxaMap[_outgroup][i]][s],
                                   _dnaMatrix[_taxaMap[p1][j]][s],
                                   _dnaMatrix[_taxaMap[hyb][k]][s],
                                   _dnaMatrix[_taxaMap[p2][r]][s], cp);
    }
  }
}

void Quartet::resolveAmbiguity(const int& one,
                               const int& two,
                               const int& three,
                               const int& four,
                               double cp[16][16]){
  /* Check to see if any combination of three of the four taxa are ambiguous. */
  double denom = 0.0;
  if((one >= 4 && two  >= 4 && three >= 4) ||
     (one >= 4 && two  >= 4 && four  >= 4) ||
     (one >= 4 && three >= 4 && four  >= 4) ||
     (two  >= 4 && three >= 4 && four  >= 4)){
    return 0.0;
  } else if(one == 4 || two == 4 || three == 4 || four == 4){ /* Don't allow gaps. */
    return 0.0;
  } else {
    denom = _baseLookup[one].size() * _baseLookup[two].size() * _baseLookup[three].size() * _baseLookup[four].size();
    for(unsigned i = 0; i < _baseLookup[one].size(); i++){
      for(unsigned j = 0; j < _baseLookup[two].size(); j++){
        for(unsigned k = 0; k < _baseLookup[three].size(); k++){
          for(unsigned r = 0; r < _baseLookup[four].size(); r++){
            cp[_baseLookup[one][i] * 4 + _baseLookup[two][j]]
              [_baseLookup[three][k] * 4 + _baseLookup[four][r]] += 1.0 / denom;
          }
        }
      }
    }
  }
}
