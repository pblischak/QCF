#ifndef QUARTET_HPP
#define QUARTET_HPP

#include <vector>
#include <string>

class SeqData;

class Quartet {
public:
  Quartet(const std::string i, const std::string j,
          const std::string k, const std::string l,
          SeqData *seq);
  ~Quartet(){};

  std::vector<double> eval(std::vector<int>& vec);
  std::vector<double> eval2(std::vector<int>& vec);
  int A, B, C, D;
  std::string hapA, hapB, hapC, hapD;
  bool resolved = 0;
  std::vector<int> index;
  std::vector< std::vector<double> > ABCD{16, std::vector<double>(16)};
  std::vector< std::vector<double> > ACBD{16, std::vector<double>(16)};
  std::vector< std::vector<double> > ADBC{16, std::vector<double>(16)};
  const std::vector< std::vector<int> > baseLookup = { /* {A,G,C,T} == {0,1,2,3} */
    {0},
    {1},
    {2},
    {3},
    {4},          /* - == ignore */
    {0, 2},       /* M == A or C */
    {0, 1},       /* R == A or G */
    {0, 3},       /* W == A ot T */
    {1, 2},       /* S == G or C */
    {2, 3},       /* Y == C or T */
    {1, 3},       /* K == G or T */
    {1, 2, 3},    /* B == C or G or T */
    {0, 1, 3},    /* D == A or G or T */
    {0, 2, 3},    /* H == A or C or T */
    {0, 1, 2},    /* V == A or G or C */
    {0, 1, 2, 3}  /* N == A or G or C or T */
  };
  SeqData* seqPtr;

private:
  void getCountMatrices_(std::vector<int>& ix);
  bool resolveAmbiguity_(const int one,
                         const int two,
                         const int three,
                         const int four);
  double resolveAmbiguity2_(const int a, const int b);
  bool ignore_amb_sites = 0;
  std::vector<double> getWeights_(std::vector<double>& vec);
  void makeIndexVec_();
};

inline double Quartet::resolveAmbiguity2_(const int a, const int b){
  double res = 0.0;
  bool match = 0;
  for(size_t i = 0; i < baseLookup[a].size(); ++i){
    for(size_t j = 0; j < baseLookup[b].size(); ++j){
      if(baseLookup[a][i] == baseLookup[b][j]) {match = 1;}
    }
  }
  if(match){
    double size1 = (double) baseLookup[a].size();
    double size2 = (double) baseLookup[b].size();
    res = 1.0 / (size1 * size2);
  } else {
    res = 0.0;
  }
  return res;
}
#endif //QUARTET_HPP
