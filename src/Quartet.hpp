#ifndef QUARTET_HPP
#define QUARTET_HPP

class SeqData;

class Quartet {
public:
  Quartet(std::string i, std::string j,
          std::string k, std::string l,
          SeqData *seq){A = i;
                        B = j;
                        C = k;
                        D = l;
                        seqPtr = seq;};
  Quartet(){};
  std::string A, B, C, D;
  double ABCD[16][16] = {{0.0}},
         ACBD[16][16] = {{0.0}},
         ADBC[16][16] = {{0.0}};
  static std::vector<std::vector<int> > baseLookup = { /* {A,G,C,T} == {0,1,2,3} */
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
  bool ignore_amb_sites = 0;
};

#endif //QUARTET_HPP
