#ifndef C_IO_HANDLER
#define C_IO_HANDLER 1

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <vector>
#include "handyDandy.h"

#ifndef DOUBLE_VECTOR
#define DOUBLE_VECTOR 1
typedef std::vector<double> double_vector;
#endif
#ifndef SHORT_VECTOR
#define SHORT_VECTOR 1
typedef std::vector<short> short_vector;
#endif
#ifndef CHAR_VECTOR
#define CHAR_VECTOR 1
typedef std::vector<char> char_vector;
#endif
#ifndef CHAR_VECTOR_VECTOR
#define CHAR_VECTOR_VECTOR 1
typedef std::vector<char_vector > char_vector_vector;
#endif
#ifndef STRING_VECTOR
#define STRING_VECTOR 1
typedef std::vector<std::string > string_vector;
#endif
#ifndef STRING_VECTOR_VECTOR
#define STRING_VECTOR_VECTOR 1
typedef std::vector<string_vector > string_vector_vector;
#endif

#ifndef MAX_LENGTH
#define MAX_LENGTH 50000
#endif

class Cio_handler {

 public:

  Cio_handler();
  /*
    read
  */
  int readOcc(char fileName[],long double prob[], long double occ0[], int &seqLen); // read in occupancy file
  int readOcc(char fileName[],long double prob[], int &seqLen); // as above, disregard occ
  /* // old: int readEn(char fileName[],long double en[], int &seqLen); // read in energy file 
     //^^// note that for most purposes the following line is required:
     //^^//       seqLen = seqLen + OBJ_SIZE -1;  */
  int readEn(char fileName[],double_vector *en); // read in energy file
  int readFilter(char fileName[],short_vector *filt); // read in filter file
  int readFasta(char fileName[],char_vector *seq); // reads only first sequence!!  
  // recognizes lowercase acgt
  // will need something like readMultiFasta(char,char_vector seqs[],string_vector seqNames[])
  //   so that i can make genome-wide measurements
  int readMultiFasta(char fileName[],char_vector_vector *seqs); // reads only first sequence!!  
  int readList(char fileName[],string_vector *elems); // read first element of each line in file
  int readList2(char fileName[],string_vector_vector *elems); // read each element in each line, split by white space
  int readTable(char fileName[],unsigned short index,double_vector *tab,bool header); // read n'th column in table
  int readTable(char fileName[],unsigned short index,double_vector *tab); // just calls readTable(,,,TRUE) 
  
  /*
    write
  */
  int writeOcc(char outFile[], int seqLen, long double prob[], long double pL[]);
  int writeEn(char outFile[], int seqLen, long double en[]);

  /*
    Error codes:
    -1 = failed to open file
    -2 = didn't find header in file
    -3 = found something in file that I didn't understand
    -4 = too short
  */

  /*
    other
  */
  void dumpOcc(int seqLen,long double prob[],long double pL[]); // print occpuancy to stdout
  char charMapper(char c); // interface to shortMap;
  
 private:
  char charMap[128]; // map fasta sequence to short
  // char shortMap[4]; // map short to readable DNA


};

#endif
