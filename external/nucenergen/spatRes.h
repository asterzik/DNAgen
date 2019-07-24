#ifndef C_SPATRES
#define C_SPATRES 1

#include "model.h"
#include <deque>

// spatially resolved model

#ifndef STRING_VECTOR
#define STRING_VECTOR 1
typedef std::vector<std::string> string_vector;
#endif
#ifndef INT_QUEUE
#define INT_QUEUE 1
typedef std::deque<int> int_queue;
#endif
#ifndef INT_QUEUE_VECTOR
#define INT_QUEUE_VECTOR 1
typedef std::vector<int_queue> int_queue_vector;
#endif


class CspatRes : public Cmodel {
 public:
  CspatRes();
  CspatRes(unsigned short N);
  CspatRes(char epsFile[],coefficientMap &coef);
  ~CspatRes() {};

  int makeRTable(char outFile[],const char_vector& seq,const double_vector& energy,const short_vector& filter,bool db); // seq = entire sequence
  int makeEnergy(char outFile[],const char_vector& seq,coefficientMap coef,short strand); // given coefficients, make an energy landscape
  
 protected: 
  int catWord(int w, char c);  // given a word w of length n, return word w.c of length n+1
  int initPastWords(const char_vector& seq,char_vector::const_iterator seq0, int_queue_vector &pastWords);
  void dumpPastWords(int_queue_vector pastWords);

};

#endif


//int readCoef(char file[],coefficientMap &coef); // sets coef.  see io_handler.h for error code return values





/*

  int_queue_vector holds currently active words.  
    once you set the struct up, all you have to do is push on the next letter and unshift the front words

  output should be a queue of rules, that way you only have to update the front.
  

    

*/
