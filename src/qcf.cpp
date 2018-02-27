#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <unordered_map>

#include "QCFData.hpp"
#include "SeqData.hpp"
#include "Quartet.hpp"
#include "Bootstrap.hpp"
#include "qcf.hpp"

int main(int argc, char* argv[]){
  /* Initial parsing of command line options to look for -h / --help,
   * -V / --version, or an incorrect number of arguments. */
  if(argc < 2){
    std::cerr << "\n** ERROR: Incorrect number of arguments. **" << std::endl;
    usage();
    exit(EXIT_FAILURE);
  }

  if(strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0){
    usage();
    exit(EXIT_SUCCESS);
  }

  if(strcmp(argv[1],"-V") == 0 || strcmp(argv[1],"--version") == 0){
    std::cerr << "\nThis is QCF version " << Version << " (" << Date << ").\n" << std::endl;
    exit(EXIT_SUCCESS);
  }

  int boot_reps = 4;
  QCFData qcf(argc, argv);
  Bootstrap boot(boot_reps);
  std::vector<int> test(10001);
  boot.randomVector(0,500,test);
  //std::cout << "[";
  for(uint i = 1; i < test.size(); i++){
    std::cout << test[i] << std::endl;
  }
  //std::cout << test[test.size() - 1] << "]" << std::endl << std::endl;
  std::vector<SeqData> seqs = qcf.get_seqs();
  for(uint i = 0; i < seqs.size(); i++){
    //std::cout << seqs[i].qcfd->taxa.size() << std::endl;
  }

  return EXIT_SUCCESS;
}

void usage(){
  std::cerr << "\nUsage: qcf -i <infile> -m <taxon-map> [additional options]" << std::endl
            << "\nInformation options:" << std::endl
            << "  -h [--help]          Prints this help message" << std::endl
            << "  -V [--version]       Prints version information" << std::endl
            << "\nProgram options:" << std::endl
            << "  -i [--infile]        Name of the input file" << std::endl
            << "  -m [--map]           Map of haplotypes to taxa" << std::endl
            << "\nAdditional options:" << std::endl
            << "  -b [--bootstrap]     Perform bootstrap resampling" << std::endl
            << "  --prefix             Append a prefix to the beginning of outfile and logfile\n" << std::endl;
}
