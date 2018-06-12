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
  QCFData(const int c, char* v[]);
  ~QCFData(){}
  std::vector<SeqData> getSeqs();
  //void map2qcf(Quartet& q);
  std::vector<std::string> taxa, seqFiles;
  std::unordered_map<std::string, int> hap2tax;
  std::string infile = "none", mapfile = "none";
  int nTaxa = 0, nHaps = 0, bootReps = 0;
  bool rawOutput = 0, quiet = 0;
  std::string prefix = "out";

private:
  void parseCommandLine_(const int ac, char* av[]);
  void checkCommandLineInput_();
  void parseMap_();
  void parseSeqsFile_();
  std::vector<std::string> splitString_(const std::string& s, char d);
};

#endif //QCFDATA_HPP
