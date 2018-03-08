#include <iostream>
#include <fstream>
#include <cstring>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>

#include "qcf.hpp"
#include "QCFData.hpp"
#include "SeqData.hpp"
#include "Quartet.hpp"

QCFData::QCFData(const int c, char* v[]){
  /*

  */
  parseCommandLine_(c, v);
  checkCommandLineInput_();
  parseSeqsFile_();
  parseMap_();
}

void QCFData::parseCommandLine_(const int ac, char* av[]){
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

void QCFData::checkCommandLineInput_(){
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
  if(bootReps < 0){
    std::cerr << "\nInvalid option for number of bootstrap replicates [-b,--bootstrap]: "
              << bootReps << std::endl;
    errorCaught++;
  }
  if(errorCaught > 0){
    std::cerr << "** ERROR: " << errorCaught << " command line option(s) improperly specified. **\n" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void QCFData::parseSeqsFile_(){
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

void QCFData::parseMap_(){
  /*

  */
  int taxonIndex = 0;
  std::string inString;
  std::vector<std::string> str1, str2;
  std::ifstream mapStream(mapfile);
  char d1 = ':', d2 = ',';
  if(mapStream.is_open()){
    while(mapStream >> inString){
      str1 = splitString_(inString, d1);
      if(str1.size() != 2){
        std::cerr << "\nError reading map file on line " << taxonIndex << "." << std::endl;
        std::cerr << "More than one \":\" present.\n" << std::endl;
        exit(EXIT_FAILURE);
      }
      taxa.push_back(str1[0]);
      str2 = splitString_(str1[1], d2);
      if(str2.size() <= 0){
        std::cerr << "\nError reading map file on line " << taxonIndex << "." << std::endl;
        std::cerr << "Haplotype names could not be split (should be comma delimted).\n" << std::endl;
        exit(EXIT_FAILURE);
      }
      for(int i = 0; i < str2.size(); i++){
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

std::vector<std::string> QCFData::splitString_(const std::string& s,
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

std::vector<SeqData> QCFData::getSeqs(){
  /*

  */
  std::vector<SeqData> seqs;
  for(int s = 0; s < seqFiles.size(); s++){
    SeqData seq(seqFiles[s], this);
    seqs.push_back(seq);
  }
  return seqs;
}
