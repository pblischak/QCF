#ifndef QCFDATA_HPP
#define QCFDATA_HPP

/* Forward declarations of classes. */
class SeqData;
class Quartet;
class QCFTable;

class QCFData {
public:
  QCFData(int c, char* v[]);
  ~QCFData(){}
  std::vector<SeqData> get_seqs();
  int _nCk(const int& n, int& k);
  void map2qcf(Quartet &q);
  std::vector<std::string> taxa, seq_files;
  std::unordered_map<std::string, int> ind2tax;
  std::unique_ptr<QCFTable> qcfs;

private:
  void _parseCommandLine(int ac, char* av[]);
  void _checkCommandLineInput();
  void _parseMap();
  void _parseSeqsFile();
  int _nTaxa = -999, _nHaps = -999, four = 4;
};

/* Calculate binomial coefficient. */
inline int QCFData::_nCk(const int& n, int& k){
  int res = 1;
  if(k > n - k)
    k = n - k;

  for(long int i = 0; i < k; i++){
    res *= (n - i);
    res /= (i + 1.0);
  }
  return res;
}

/* Helper class: QCFTable
   The table is a 3-dimensional ragged array.
   Each row is a taxon-level quartet indexed by the order
   of the (N choose 4) possible quartets. Then, each row gets
   a vector of triple scores (also a vector) for each gene tree.
   If no bootstrapping is conducted, only one entry will be positive,
   the others will be 0.0. If bootstrapping is conducted, then the scores
   will be the weighted scores across bootstrap replicates.
*/
class QCFTable {
public:
  QCFTable(int nqrts){qtab.resize(nqrts);}
  ~QCFTable(){};
  std::vector< std::vector < std::vector<double> > > qtab;
  void add_entry(std::vector<double> &entry, int i, int j, int k, int l);
  void write(std::string &outfile);
};

#endif //QCFDATA_HPP
