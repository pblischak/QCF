#ifndef BOOTSTRAP_HPP
#define BOOTSTRAP_HPP

#include <vector>
#include "MbRandom.hpp"

class Quartet;

class Bootstrap {
public:
  Bootstrap(int &nreps);
  ~Bootstrap(){delete r;}
  std::vector<double> operator()(Quartet &qrt);
  void randomVector(int low, int high, std::vector<int> &vec);
  MbRandom* r;
  int reps = 1;
};

#endif //BOOTSTRAP_HPP
