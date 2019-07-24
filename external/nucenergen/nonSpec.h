#ifndef C_NONSPEC
#define C_NONSPEC 1

#include "model.h"

#ifndef STRING_VECTOR
#define STRING_VECTOR 1
typedef std::vector<std::string> string_vector;
#endif


class CnonSpec : public Cmodel {
 public:
  CnonSpec();
  CnonSpec(unsigned short N);
  CnonSpec(char epsFile[],coefficientMap &coef);
  ~CnonSpec() {};

  int makeRTable(char outFile[],const char_vector& seq,const double_vector& energy,const short_vector& filter,bool db); // seq = entire sequence
  int makeEnergy(char outFile[],const char_vector& seq,coefficientMap coef,short strand); // given coefficients, make an energy landscape
  
 protected:
  void shiftFrontWords(string_vector &words, char nextLetter); // keep track of which words enter the 'active region'
  void shiftBackWords(string_vector &words, char nextLetter); // keep track of which words leave the 'active region'  

 private:
  int initCounts(const char_vector& seq,char_vector::const_iterator seq0, short_vector_vector &counts,string_vector &forwardWords,string_vector &backWords); // assumes you've checked filter already
  void dumpCounts(short_vector_vector counts);
};

#endif


//int readCoef(char file[],coefficientMap &coef); // sets coef.  see io_handler.h for error code return values

