#ifndef QCFTABLE_HPP
#define QCFTABLE_HPP

#include <vector>
#include <string>
/* Helper class: QCFTable
   The table is a 3-dimensional ragged array.
   Each row is a taxon-level quartet indexed by the order
   of the (N choose 4) possible quartets. Then, each row gets
   a vector of triple scores (also a vector) for each gene tree.
   If no bootstrapping is conducted, only one entry will be positive,
   the others will be 0.0. If bootstrapping is conducted, then the scores
   will be the weighted scores across bootstrap replicates.
*/

class QCFData;

class QCFEntry {
public:
  QCFEntry(const int i, const int j, const int k,
           const int l, const int idx){key = {i,j,k,l};
                                       index = idx;}
  ~QCFEntry(){};
  std::vector<int> key;
  int index;
};

class QCFTable {
public:
  QCFTable(QCFData* qcf);
  ~QCFTable(){};
  std::vector<QCFEntry> lookup;
  std::vector< std::vector< std::vector<double> > > values;
  //std::vector< std::vector < std::vector<double> > > qtab;
  int nCk(const int& n, int& k);
  void addValue(const int i, const int j, const int k,
                const int l, const std::vector<double> val);
  void write(const std::string pfx);
  int findIndex(const int i, const int j,
                const int k, const int l);
  int four = 4, nTaxa, nQrts;
  QCFData* qcfPtr;
};

/* Calculate binomial coefficient. */
inline int QCFTable::nCk(const int& n, int& k){
  int res = 1;
  if(k > n - k)
    k = n - k;

  for(long int i = 0; i < k; i++){
    res *= (n - i);
    res /= (i + 1.0);
  }
  return res;
}

#endif //QCFTABLE_HPP
