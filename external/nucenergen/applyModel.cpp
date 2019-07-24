#include <iostream>
#include <fstream>
#include "io_handler.h"
#include "option2.h"
#include "nonSpec.h"
#include "spatRes.h"

int main  (int argc, char* argv[]) {

  const int DEF_STRAND = 1;
  //const int DEF_N = 2;
  const int DEF_MOD = 1;
  std::ostringstream os;
  os << "usage:\t -f fasta.fna -eps model.eps -out output.table < -mod " << DEF_MOD << " -strand " << DEF_STRAND << ">\n"
     <<"\t\tmodels:\t 1 = position independent (default), 2 = spatially resolved\n"
     << "\t\tstrand:\t -1 = reverse, 0 = both, 1 = forward (default)\n";

  //     << "\t\tmodels: 1=nonSpecific";
  std::string usage(os.str());

  string_vector args;
  args.push_back("f");   // 0
  args.push_back("eps"); // 1
  args.push_back("out"); // 2
  args.push_back("mod"); // 3
  args.push_back("strand"); // 4
  string_vector opts = option2_parseOptions(argc,argv,args);
  if (opts.size()<args.size()) { // if you didn't provide any options, or if, somehow, the parser failed
    std::cout << usage;
    return 1;
  }

  char fastaFile[99];
  sprintf(fastaFile,"%s",opts[0].c_str());
  char epsFile[99];
  sprintf(epsFile,"%s",opts[1].c_str());
  char outFile[99];
  sprintf(outFile,"%s",opts[2].c_str());
  int mod = DEF_MOD;
  if (opts[3].length()) {
    mod = option2_stringToInt(opts[3]);
  } 
  int strand = DEF_STRAND;
  if (opts[4].length()) {
    strand = option2_stringToInt(opts[4]);
  }

  Cio_handler myIo;
  
  char_vector * seq = new char_vector;
  myIo.readFasta(fastaFile,seq);

  Cmodel * model;
  coefficientMap coef; 
  switch(mod) {
  case 1:
    model = new CnonSpec(epsFile,coef);
    break;
  case 2:
    model = new CspatRes(epsFile,coef);
    break;
  default:
    std::cout << usage;
    exit(1);
  }
  
  model->makeEnergy(outFile,*seq,coef,(short)strand);
  return 0;
}


  /*
    // verify readCoef
  std::string s;
  coefficientMap::iterator it1;
  double_vector::iterator it2;
  int i;
  for (it1=coef.begin(); it1!=coef.end(); ++it1) {
    s = it1->first;
    i=0;
    for(it2=it1->second.begin();it2!=it1->second.end();++it2) {
      std::cout << s << ++i << ") " << *it2 << "\n";
    }
  }
  return 0;
  printf("outFile = '%s'\n",outFile);
  printf("length of seq = %d\n",(int)seq->size());
  printf("strand = %d\n",strand);
  */
