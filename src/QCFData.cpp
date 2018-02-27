#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <memory>

#include "QCFData.hpp"
#include "SeqData.hpp"
#include "Quartet.hpp"

QCFData::QCFData(int c, char* v[]){
  _parseCommandLine(c, v);
  _checkCommandLineInput();
  ind2tax.reserve(_nHaps);
  std::unique_ptr<QCFTable> qcfs(new QCFTable(_nCk(_nHaps, four)));
  _parseSeqsFile();
  _parseMap();
}

void QCFData::_parseCommandLine(int ac, char* av[]){

}

std::vector<SeqData> QCFData::get_seqs(){
  std::vector<SeqData> seqs;
  std::string t = "test.txt";
  SeqData seq1(t, this);
  SeqData seq2(t, this);
  SeqData seq3(t, this);
  seqs.push_back(seq1);
  seqs.push_back(seq2);
  seqs.push_back(seq3);
  return seqs;
}

void QCFData::_checkCommandLineInput(){

}

void QCFData::_parseSeqsFile(){

}

void QCFData::_parseMap(){

}
