#ifndef BOOTSTRAP_HPP
#define BOOTSTRAP_HPP

#include <vector>

class Quartet;
class MbRandom;

class Bootstrap {
public:
  Bootstrap(int &nreps);
  ~Bootstrap(){};
  std::vector<double> operator()(Quartet &qrt);
  void randomVector(int low, int high, std::vector<int> &vec);
  MbRandom* r;
  int reps = 1;
};

#endif //BOOTSTRAP_HPP
