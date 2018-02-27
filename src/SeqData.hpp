#ifndef SEQDATA_HPP
#define SEQDATA_HPP

class QCFData;

class SeqData {
public:
  SeqData(std::string &name,
          QCFData* ptr){file_name = name;
                        qcfd = ptr;}
  ~SeqData(){}
  std::string file_name;
  std::vector< std::vector<int> > dna;
  std::vector<std::string> haps;
  QCFData* qcfd;
private:
  int _convert(char str);
};

/* Convert DNA bases to ints as they are read in. */
inline int SeqData::_convert(char str){
  int _baseCode = 999;
  switch (str) {
    case 'A': case 'a': _baseCode = 0;  break;
    case 'G': case 'g': _baseCode = 1;  break;
    case 'C': case 'c': _baseCode = 2;  break;
    case 'T': case 't': case 'U': case 'u': _baseCode = 3;  break;
    case 'M': case 'm': _baseCode = 5;  break;
    case 'R': case 'r': _baseCode = 6;  break;
    case 'W': case 'w': _baseCode = 7;  break;
    case 'S': case 's': _baseCode = 8;  break;
    case 'Y': case 'y': _baseCode = 9;  break;
    case 'K': case 'k': _baseCode = 10;  break;
    case 'B': case 'b': _baseCode = 11;  break;
    case 'D': case 'd': _baseCode = 12;  break;
    case 'H': case 'h': _baseCode = 13;  break;
    case 'V': case 'v': _baseCode = 14;  break;
    case 'N': case 'n': _baseCode = 15;  break;
    case '-': _baseCode = 4; break;
    case '?': _baseCode = 15; break;
    default : {std::cerr << "\n** ERROR: Unrecognized base character in DNA matrix: \"" << str << "\". **\n" << std::endl;
               exit(EXIT_FAILURE);} break;
  }
  return _baseCode;
}

#endif //SEQDATA_HPP
