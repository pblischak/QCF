#ifndef BOOTSTRAP_HPP
#define BOOTSTRAP_HPP

class Quartet;

class Bootstrap {
public:
  Bootstrap(int &nreps){reps = nreps;}
  ~Bootstrap(){};
  int reps = 1;
  std::vector<double> operator()(Quartet &qrt);
};

#endif //BOOTSTRAP_HPP
