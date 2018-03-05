#include <iostream>
#include <cstring>
#include <vector>

#include "qcf.hpp"
#include "QCFTable.hpp"
#include "QCFData.hpp"
#include "SeqData.hpp"
#include "Quartet.hpp"
#include "Bootstrap.hpp"

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
    std::cerr << "\nThis is QCF version " << VERSION << " (" << DATE << ").\n" << std::endl;
    exit(EXIT_SUCCESS);
  }

  std::vector<int> empty; // An emty vector to pass if no bootstrapping is done.
  bool bootstrap = 0;
  QCFData qcf(argc, argv);
  QCFTable table(&qcf);
  Bootstrap boot(qcf.bootReps);
  std::vector<double> res(3,0.0);
  if(boot.reps > 0){
    bootstrap = 1;
  }
  std::vector<SeqData> seqs = qcf.get_seqs();
  for(uint s = 0; s < seqs.size(); s++){
    if(seqs[s].skip){
      continue;
    }
    std::cout << "Analyzing gene " << s+1 << " (" << seqs[s].haps.size() << "): ";
    std::vector<Quartet> qrts = seqs[s].get_quartets();
    std::cout << qrts.size() << " quartets.\n";
    for(uint q = 0; q < qrts.size(); q++){
      if(bootstrap){
        res = boot(qrts[q]);
      } else {
        res = qrts[q].eval2(empty);
      }
      //std::cout << res[0] << "\t" << res[1] << "\t" << res[2] << std::endl;
      //std::cout << qcf.hap2tax[qrts[q].hapA] << std::endl;
      //std::cout << qcf.hap2tax[qrts[q].hapB] << std::endl;
      //std::cout << qcf.hap2tax[qrts[q].hapC] << std::endl;
      //std::cout << qcf.hap2tax[qrts[q].hapD] << std::endl;
      table.addValue(qcf.hap2tax[qrts[q].hapA],
                     qcf.hap2tax[qrts[q].hapB],
                     qcf.hap2tax[qrts[q].hapC],
                     qcf.hap2tax[qrts[q].hapD], res);
    }
  }
  table.write(qcf.prefix);
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
