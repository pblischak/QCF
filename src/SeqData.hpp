#ifndef SEQDATA_HPP
#define SEQDATA_HPP

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <assert.h>

class QCFData;
class Quartet;

class SeqData {
public:
  SeqData(const std::string &name,
          QCFData* qcf){file_name = name;
                        qcfPtr = qcf;
                        readPhylip_();}
  ~SeqData(){}
  std::vector<Quartet> getQuartets();
  std::string file_name;
  std::vector< std::vector<int> > dna;
  std::vector<std::string> haps;
  std::unordered_map<std::string, int> seqIndex;
  std::vector<int> orderHaps(const int i, const int j,
                             const int k, const int l);
  QCFData* qcfPtr;
  int nSeqs = -999, nSites = -999;
  bool skip = 0;
private:
  int convert_(const char str);
  //void checkHap_(std::string& hapName);
  void readPhylip_();
};

/* Convert DNA bases to ints as they are read in. */
inline int SeqData::convert_(const char str){
  int baseCode_ = -999;
  switch (str) {
    case 'A': case 'a': baseCode_ = 0;  break;
    case 'G': case 'g': baseCode_ = 1;  break;
    case 'C': case 'c': baseCode_ = 2;  break;
    case 'T': case 't': case 'U': case 'u': baseCode_ = 3;  break;
    case 'M': case 'm': baseCode_ = 5;  break;
    case 'R': case 'r': baseCode_ = 6;  break;
    case 'W': case 'w': baseCode_ = 7;  break;
    case 'S': case 's': baseCode_ = 8;  break;
    case 'Y': case 'y': baseCode_ = 9;  break;
    case 'K': case 'k': baseCode_ = 10;  break;
    case 'B': case 'b': baseCode_ = 11;  break;
    case 'D': case 'd': baseCode_ = 12;  break;
    case 'H': case 'h': baseCode_ = 13;  break;
    case 'V': case 'v': baseCode_ = 14;  break;
    case 'N': case 'n': baseCode_ = 15;  break;
    case '-': baseCode_ = 4; break;
    case '?': baseCode_ = 15; break;
    default : {std::cerr << "\n** ERROR: Unrecognized base character in DNA matrix: \"" << str << "\". **\n" << std::endl;
               exit(EXIT_FAILURE);} break;
  }
  assert(baseCode_ != -999);
  return baseCode_;
}

#endif //SEQDATA_HPP
