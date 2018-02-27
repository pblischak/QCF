#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>

#include "QCFData.hpp"
#include "SeqData.hpp"
#include "Quartet.hpp"

QCFData::QCFData(int c, char* v[]){
  /*

  */
  _parseCommandLine(c, v);
  _checkCommandLineInput();
  _parseSeqsFile();
  _parseMap();
  std::unique_ptr<QCFTable> qcfs(new QCFTable(_nCk(nHaps, four)));
}

void QCFData::_parseCommandLine(int ac, char* av[]){
  /*

  */
  int invalidArgCount = 0;
  std::vector<char*> invalidArgs;
  for(int i = 1; i < ac; i++){
    if(strcmp(av[i], "-i") == 0 || strcmp(av[i], "--infile") == 0){
      infile = av[i + 1];
    } else if(strcmp(av[i], "-m") == 0 || strcmp(av[i], "--map") == 0){
      mapfile = av[i + 1];
    } else if (strcmp(av[i], "-b") == 0 || strcmp(av[i], "--bootstrap") == 0){
      bootReps = atoi(av[i + 1]);
    } else if(strcmp(av[i], "--prefix") == 0){
      prefix = av[i + 1];
    } else if(strcmp(av[i], "-q") == 0 || strcmp(av[i], "--quiet") == 0){
      quiet = 1;
    } else if(av[i][0] == '-'){ /* This checks if it is a flag that isn't valid. */
      invalidArgCount++;
      invalidArgs.push_back(av[i]);
    }
  }

  if(invalidArgCount != 0){
    std::cerr << "\n** ERROR: Unrecognized command line flag(s). **\n" << std::endl;
    for(unsigned i = 0; i < invalidArgs.size(); i++){
      std::cerr << "   " << invalidArgs[i] << std::endl;
    }
    std::cerr << "\nType 'qcf -h' for command line options.\n" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void QCFData::_checkCommandLineInput(){
  /*

  */
  int errorCaught = 0;
  if(strcmp(infile.c_str(), "none") == 0){
    std::cerr << "\nMissing or invalid option for -i [--infile]: " << infile << std::endl;
    errorCaught++;
  }
  if(strcmp(mapfile.c_str(),"none") == 0){
    std::cerr << "\nMissing or invalid option for -m [--map]: " << mapfile << std::endl;
    errorCaught++;
  }
  if(bootReps <= 0 && bootReps != -999){
    std::cerr << "\nInvalid option for number of bootstrap replicates [-b,--bootstrap]: "
              << bootReps << std::endl;
  }
  if(errorCaught > 0){
    std::cerr << "** ERROR: " << errorCaught << " command line option(s) improperly specified. **\n" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void QCFData::_parseSeqsFile(){
  /*

  */
  std::string inValue;
  std::ifstream inStream(infile);
  if(inStream.is_open()){
    while(inStream >> inValue){
      seqFiles.push_back(inValue);
    }
  } else {
    std::cerr << "Could not open input file: " << infile << std::endl;
    exit(EXIT_FAILURE);
  }
}

void QCFData::_parseMap(){
  /*

  */
  int taxonIndex = 0;
  std::string inString;
  std::vector<std::string> str1, str2;
  std::ifstream mapStream(mapfile);
  char d1 = ':', d2 = ',';
  if(mapStream.is_open()){
    while(mapStream >> inString){
      str1 = _splitString(inString, d1);
      if(str1.size() != 2){
        std::cerr << "\nError reading map file on line " << taxonIndex << "." << std::endl;
        std::cerr << "More than one \":\" present.\n" << std::endl;
        exit(EXIT_FAILURE);
      }
      taxa.push_back(str1[0]);
      str2 = _splitString(str1[1], d2);
      if(str2.size() <= 0){
        std::cerr << "\nError reading map file on line " << taxonIndex << "." << std::endl;
        std::cerr << "Haplotype names could not be split (should be comma delimted).\n" << std::endl;
        exit(EXIT_FAILURE);
      }
      for(uint i = 0; i < str2.size(); i++){
        hap2tax.insert({str2[i], taxonIndex});
        nHaps++;
      }
      taxonIndex++;
    }
    nTaxa = taxa.size();
  } else {
    std::cerr << "Could not open map file: " << mapfile << std::endl;
    exit(EXIT_FAILURE);
  }
}

std::vector<std::string> QCFData::_splitString(const std::string &s,
                                               char d){
  /*
  Split an input string (s) based on character delimeter (d).
  Return string tokens as a vector.
  */
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, d)){
    tokens.push_back(token);
  }
  return tokens;
}

std::vector<SeqData> QCFData::get_seqs(){
  /*

  */
  std::vector<SeqData> seqs;
  for(uint s = 0; s < seqFiles.size(); s++){
    SeqData seq(seqFiles[s], this);
    seqs.push_back(seq);
  }
  /*std::string t = "test.txt";
  SeqData seq1(t, this);
  SeqData seq2(t, this);
  SeqData seq3(t, this);
  seqs.push_back(seq1);
  seqs.push_back(seq2);
  seqs.push_back(seq3);*/
  return seqs;
}
