/* this is an interface to the various models for nucleosome specificity */

#ifndef C_MODEL
#define C_MODEL 1

#include <string>
#include <vector>
#include <map> // note: trying to access an undefined element of a map with [] will _instantiate_ that object if it is not defined!!
//#include <hash_map>

#ifndef DOUBLE_VECTOR
#define DOUBLE_VECTOR 1
typedef std::vector<double> double_vector;
#endif
#ifndef DOUBLE_VECTOR_VECTOR
#define DOUBLE_VECTOR_VECTOR 1
typedef std::vector<double_vector> double_vector_vector;
#endif
#ifndef COEFFICIENT_MAP 
#define COEFFICIENT_MAP 1
typedef std::vector <double_vector_vector> coefficientMap; // coef[AA123] = coef[1][0][122]
//                                                         //                    n length of word
//                                                         //                       w wordToInt           
//                                                         //                          l position of word - 1
//                                                         // store intercept at coef[0][4][0]
#endif
#ifndef CHAR_VECTOR
#define CHAR_VECTOR 1
typedef std::vector<char> char_vector;
#endif
#ifndef CHAR_VECTOR_VECTOR
#define CHAR_VECTOR_VECTOR 1
typedef std::vector<char_vector > char_vector_vector;
#endif
#ifndef SIZE_T_VECTOR // need for word vectors, in case of long words
#define SIZE_T_VECTOR 1
typedef std::vector<size_t> size_t_vector;
#endif
#ifndef SHORT_VECTOR
#define SHORT_VECTOR 1
typedef std::vector<short> short_vector;
#endif
#ifndef SHORT_VECTOR_VECTOR
#define SHORT_VECTOR_VECTOR 1
typedef std::vector<short_vector> short_vector_vector;
#endif
#ifndef TRULES
#define TRULES 1
typedef std::vector<short_vector_vector> Trules;
#endif

/*
  notes on Trules
  Trules[n-1] = all the rules for words of length n
  Trules[n-1][i] = all the rules for word i (where i = wordToInt[i])
  Trules[n-1][i][j] = the specific rule for mapping word i to word j
  so, word 'C' is Trules[0][3] = (-1,-1,-1)
  rules for 'CA' to GA is TRules[1][12][4]=-1
*/

class Cmodel {
 public:
  explicit Cmodel();
  explicit Cmodel(short unsigned N);
  explicit Cmodel(char epsFile[],coefficientMap &coef); // read from file, set coef & N
  virtual ~Cmodel();
  int readCoef(char file[],coefficientMap &coef,unsigned short& N); // sets coef.  see io_handler.h for error code return values
  int makeCombinedRTable(char baseOutFile[],const char_vector_vector& seqs,const double_vector_vector& energies,const short_vector_vector& filters,size_t maxBp); //
  // take a list of sequences, energies, filters, and make R-tables out of all of them.                                               <-  //
  // the tables concatenate the inputs even as they break them up.  The idea is to do a regression on all your data without           <-  //
  // overflowing RAM with R-Tables that are too big.                                                                                  <-  //
  // Args are lists of <TYPE> where TYPE is as for makeRTable                                                                         <-  //

  virtual int makeRTable(char outFile[],const char_vector &seq,const double_vector &energy,const short_vector &filter,bool db)=0; // seq = entire sequence
  virtual int makeEnergy(char outFile[],const char_vector &seq,coefficientMap coef,short strand)=0; // given coefficients, make an energy landscape
  
  //void test();
                                                                                             // strand= (0,1,-1) = (both,F,R)

 protected:
  void defineCharMap();
  void init(short unsigned N);
  Trules makeRules(); // return a 3d array with all the rules up to a given N (see notes above)
  int wordToInt(std::string w); // base 4 map AA = 0, AT = 1, AG = 2,AC=3, TA = 4, and so forth
  std::string intToWord(short n, int i); // base 4 map AA = 0, AT = 1, AG = 2, TA = 4, and so forth note: n is one-based!
  bool threesMatch(int i, int j,short_vector threes); // if i has threes and j matches i's non-threes, then return 1, else 0
  bool hasThrees(short n,int i);
  int findThrees(short n,int i,short_vector &threes); // ret is number of threes in word
  int addRules(short_vector &rule1,short_vector rule2); // vector addition
  void scaleRule(short s,short_vector &rule); // vector scalar multiplication
  int rc(short n,int i); // get reverse compliment of word i
  void fillCoef(coefficientMap &coef); // add in energy values for words including C using rules
  void dumpCoef(coefficientMap coef);

  unsigned short N; // max length of word
  char shortMap[4]; // short to char
  short charMap[128]; // (int)char to short  
  unsigned short OBJ_SIZE;
  short OUTSIDE;  // size of excluded region on the flanks of the nucleosome
  Trules rules;
};

#endif


/* usage notes
given foo,bar : public Cmodel
junk() {
  std::string chooseModel;
  Cmodel * model;
  if (chooseModel.compare("foo")) {
    model = new foo(args); 
  } else {
    model = new bar(args);
  }
  model->doStuff();
}

*/
