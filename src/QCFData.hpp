#ifndef QCFDATA_HPP
#define QCFDATA_HPP

#include <vector>
#include <unordered_map>
#include <string>

/* Forward declarations of classes. */
class SeqData;
class Quartet;

class QCFData {
public:
  QCFData(int c, char* v[]);
  ~QCFData(){}
  std::vector<SeqData> get_seqs();
  void map2qcf(Quartet &q);
  std::vector<std::string> taxa, seqFiles;
  std::unordered_map<std::string, uint> hap2tax;
  std::string infile = "none", mapfile = "none";
  int nTaxa = 0, nHaps = 0, bootReps = -999;
  bool quiet = 0;
  std::string prefix = "out";

private:
  void _parseCommandLine(int ac, char* av[]);
  void _checkCommandLineInput();
  void _parseMap();
  void _parseSeqsFile();
  std::vector<std::string> _splitString(const std::string &s, char d);
};

#endif //QCFDATA_HPP
