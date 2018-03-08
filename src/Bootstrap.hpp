#ifndef BOOTSTRAP_HPP
#define BOOTSTRAP_HPP

#include <vector>
#include "MbRandom.hpp"

class Quartet;

class Bootstrap {
public:
  Bootstrap(const int nReps);
  ~Bootstrap(){delete r;}
  std::vector<double> operator()(Quartet& qrt);
  int reps = -999;
private:
  void randomVector_(const int low, const int high, std::vector<int>& vec);
  MbRandom* r;
};

#endif //BOOTSTRAP_HPP
