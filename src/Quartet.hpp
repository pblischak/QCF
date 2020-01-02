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
  std::vector<double> eval3(std::vector<int>& vec);
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
  double LogDet_(double mat[4][4], double norm);
  double determinant_(double mat[4][4]);
  void getCountMatrices_(std::vector<int>& ix);
  bool resolveAmbiguity_(const int one,
                         const int two,
                         const int three,
                         const int four);
  double resolveAmbiguity2_(const int a, const int b);
  double resolveAmbiguity3_(const int a, const int b, double mat[4][4]);
  bool ignore_amb_sites = 0;
  std::vector<double> getWeights_(std::vector<double>& vec);
  std::vector<double> getWeights2_(std::vector<double>& vec);
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

inline double Quartet::determinant_(double mat[4][4]){
  /*
  Tried to write some recursive code based on stuff I saw
  online. Didn't work -- we'll just do it the long way.
  http://www.euclideanspace.com/maths/algebra/matrix/functions/determinant/fourD/index.htm
  */
  return mat[0][3] * mat[1][2] * mat[2][1] * mat[3][0] - mat[0][2] * mat[1][3] * mat[2][1] * mat[3][0] -
         mat[0][3] * mat[1][1] * mat[2][2] * mat[3][0] + mat[0][1] * mat[1][3] * mat[2][2] * mat[3][0] +
         mat[0][2] * mat[1][1] * mat[2][3] * mat[3][0] - mat[0][1] * mat[1][2] * mat[2][3] * mat[3][0] -
         mat[0][3] * mat[1][2] * mat[2][0] * mat[3][1] + mat[0][2] * mat[1][3] * mat[2][0] * mat[3][1] +
         mat[0][3] * mat[1][0] * mat[2][2] * mat[3][1] - mat[0][0] * mat[1][3] * mat[2][2] * mat[3][1] -
         mat[0][2] * mat[1][0] * mat[2][3] * mat[3][1] + mat[0][0] * mat[1][2] * mat[2][3] * mat[3][1] +
         mat[0][3] * mat[1][1] * mat[2][0] * mat[3][2] - mat[0][1] * mat[1][3] * mat[2][0] * mat[3][2] -
         mat[0][3] * mat[1][0] * mat[2][1] * mat[3][2] + mat[0][0] * mat[1][3] * mat[2][1] * mat[3][2] +
         mat[0][1] * mat[1][0] * mat[2][3] * mat[3][2] - mat[0][0] * mat[1][1] * mat[2][3] * mat[3][2] -
         mat[0][2] * mat[1][1] * mat[2][0] * mat[3][3] + mat[0][1] * mat[1][2] * mat[2][0] * mat[3][3] +
         mat[0][2] * mat[1][0] * mat[2][1] * mat[3][3] - mat[0][0] * mat[1][2] * mat[2][1] * mat[3][3] -
         mat[0][1] * mat[1][0] * mat[2][2] * mat[3][3] + mat[0][0] * mat[1][1] * mat[2][2] * mat[3][3];
}

#endif //QUARTET_HPP
