#include "spatRes.h"
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include "handyDandy.h"


CspatRes::CspatRes() {}

CspatRes::CspatRes(unsigned short n) {
  defineCharMap();
  init(n);

  int_queue iq;
  int i;
  for (i=0;i<10;i++) {
    iq.push_back(i+1);
  }
}

CspatRes::CspatRes(char epsFile[], coefficientMap &coef) {

  defineCharMap();

  unsigned short n;
  int ret = readCoef(epsFile,coef,n);

  if (ret!=1) {
    std::cerr << "unable to read eps file " << epsFile << "\n";
    exit(ret);
  }
  init(n);
  fillCoef(coef);


  /*
  dumpCoef(coef);
  exit(1);

  debug catWord(int,char)
  char c = 'C';
  for (int i=0;i<16;i++) {
    int k = catWord(i,c);
    //std::cout << intToWord(2,i) << "." << c << " = " << intToWord(3,k) << "\n";
  }

  exit(0);
  */
}


// seq = entire sequence
int CspatRes::makeRTable(char outFile[],const char_vector& seq,const double_vector& energy,const short_vector& filter,bool db) {

  short strand = 0;

  if (!db && energy.size()!=seq.size()-OBJ_SIZE+1) {
    printf("WARNING: energy size %d doesn't match sequence.  expected %d\n",(int)energy.size(), (int)seq.size()-OBJ_SIZE+1);
  }
  /*
   passing energy and filter as const references, so you have to do this elsewhere to make a database.
 
    else if (db) {
    energy.clear();
    for (int i=1;i<=(int)seq.size()-OBJ_SIZE+1;i++) {
      energy.push_back(i);
    }
    filter.clear();
    filter.resize(seq.size(),1);
  }
  */
  if (seq.size()<OBJ_SIZE) {
    printf("fasta sequence is smaller than object size\n");
    return -4;
  }

  FILE * oFile;  
  oFile = fopen (outFile,"w"); // "w" is over-write
  
  /*
    Make Header 
  */
  std::string s;
  fprintf(oFile,"%14s\t","E");
  std::stringstream st1(std::stringstream::in | std::stringstream::out);
  int n,d,i;
  for (n=1;n<=N;n++) {
    for (d=1;d<=OBJ_SIZE-2*OUTSIDE-n+1;d++) {
      for (i=0;i<pow(4,n);i++) {
	if (hasThrees(n,i)) {
	  continue;
	}
	st1.str("");
	st1 << intToWord(n,i);
	st1 << d;
	s = st1.str();
	st1.clear();
	const char *p = s.c_str();
	fprintf(oFile,"\t%s",p);
      }
    }
  }
  fprintf(oFile,"\n");

  int_queue_vector pastWords; // starting at position i, pastWords[n-1][j] = word (length n) at position i+j+OUTSIDE
  // initialize pastWords
  int_queue iq;
  for (n = 1; n<=N; n++) {
    pastWords.push_back(iq);
  }
  int_queue_vector theseWords;
  int_queue_vector::iterator iqvi;
  int_queue::iterator iqi;
  int_queue::reverse_iterator riqi;

  char_vector::const_iterator seq1,seq0; // seq1 reads the next, seq0 reads the end
  seq0 = seq1 = seq.begin(); // used exclusively to initialize counts. 
  seq0+=OUTSIDE; // the nucleotide at the start of the active region -- used only in initialization
  seq1+=OBJ_SIZE-OUTSIDE-1; // the next nucleotide to enter the active area (ie what's under the nucleotide but not in the 'outside' excluded region
  
  short_vector::const_iterator f; // filter iterator;
  f = filter.begin();
  double_vector::const_iterator en; // energy iterator;
  en = energy.begin();
  int seqN = 1; // unused variable
  bool initFlag=true;  // flag this when the counts variable needs to be refilled
  bool updateFlag=false; // set to false when you've just initialized counts
  bool writeFlag=false; // set to true when you want to write this line
  Trules myRules; // all the rules for this position: myRules[i-1][n-1] = vector of all rules for words length n at position i 
  //              // (eg myRules[3][2][1] = AAG4 rule (3 -> 4, 2-> length 3 word, 1 -> AAG)
  short_vector rule1,rule2; 
  short_vector::iterator svi;
  int ret;
  char c;
  //int iSaySo=0; // garbledina
  /*
    MAIN LOOP: calculate the r table contents and write them to file
  */
  while (en!=energy.end()) {    
    if (seqN%10000==0) {
      printf("%d...\n",seqN);
    }

    if (seq1==seq.end() || f==filter.end()) {
      printf("found end of sequence inside seqN = %d\n",seqN);
      return -3;
    }
    
    if (initFlag && *f==1) { 
      ret = initPastWords(seq,seq0,pastWords);
      if (ret!=1) {
	printf("there was some error (%d) in initCounts\n",ret);
	return ret;
      }
      initFlag=false;
      updateFlag=false;
      writeFlag=true;
    } else if (*f==0) {
      updateFlag=false;
      initFlag=true;
      writeFlag=false;
    }

    if (updateFlag) {
      c = *seq1;
      for (n=N-1;n>0; n--) { 
	pastWords[n].push_back(catWord(pastWords[n-1].back(),c)); // next word is the new letter plus the old length n-1 word 
	pastWords[n].pop_front(); // drop the last word
      }
      pastWords[0].push_back(charMap[(int)c]);
      pastWords[0].pop_front();
    }
    
    if (writeFlag) {
      // calculate rules
      if (strand >= 0) { // both or F

	// print line of table
	if (db) { 
	  fprintf(oFile,"%14d\t",(int)*en);
	} else {
	  fprintf(oFile,"%14.6e\t",*en);
	}
	
	n=1;
	for (iqvi = pastWords.begin();iqvi<pastWords.end();iqvi++,n++) {
	  d=1;
	  for (iqi = iqvi->begin(); iqi<iqvi->end(); iqi++,d++) {
	    rule1 = rules[n-1][*iqi];
	    i=0;
	    //printf ("\tn = %d\n",n);
	    for (svi = rule1.begin();svi<rule1.end();svi++,i++) {
	      if (hasThrees(n,i)) {
		continue;
	      }
	      fprintf(oFile,"%d\t",*svi);
	    }
	  }	  
	}
	fprintf(oFile,"\n");
      } 
      if (strand <= 0) { // both or R
	// reverse compliment

	// print line of table
	if (db) { 
	  fprintf(oFile,"%14d\t",-(int)*en);
	} else {
	  fprintf(oFile,"%14.6e\t",*en);
	}

	n=1;
	for (iqvi = pastWords.begin();iqvi!=pastWords.end();++iqvi) {
	  for (riqi = iqvi->rbegin(); riqi<iqvi->rend(); riqi++) { // iterate backwards
	    rule1 = rules[n-1][rc(n,*riqi)]; // get reverse compliment of word
	    i=0;
	    for (svi = rule1.begin();svi<rule1.end();svi++,i++) {
	      if (hasThrees(n,i)) {
		continue;
	      }
	      fprintf(oFile,"%d\t",*svi);
	    }
	  }
	  n++;
	}	
	fprintf(oFile,"\n");
      }
      updateFlag=true;      
    }
    en++;
    f++;
    seq0++;
    seq1++;
    seqN++; // unused variable
  }
  fclose(oFile);
  
  //printf("seqN = %d\n",seqN);
  return 1;
}
// given coefficients, make an energy landscape
int CspatRes::makeEnergy(char outFile[],const char_vector& seq,coefficientMap coef,short strand) {
  
  FILE * oFile;  
  oFile = fopen (outFile,"w"); // "w" is over-write
  if(ferror(oFile)) {
    std::cerr << "spatRes::makeEnergy failed to open '" << outFile 
	      << "' to write to. error = " << ferror(oFile) << "\n";
  }
  /*
    Make Header 
  */
  fprintf(oFile,"%15s","seqN");
  if (strand==0 || strand==1) {
    fprintf(oFile," %14s","E");
  }
  if (strand==0 || strand==-1) {
    fprintf(oFile," %14s","E_rc");
  }
  fprintf(oFile,"\n");

  int_queue_vector pastWords; // starting at position i, pastWords[n-1][j] = word (length n) at position i+j+OUTSIDE
  // initialize pastWords
  int_queue sq; 
  for (int n = 1; n<=N; n++) {
    pastWords.push_back(sq);
  }
  int_queue_vector::iterator iqvi;
  int_queue::iterator iqi;
  int_queue::reverse_iterator riqi;
  
  char_vector::const_iterator seq1, seq0; // seq1: next word, seq0: first word.  
  
  seq0 = seq1 = seq.begin(); // used exclusively to initialize pastWords. 
  seq0+=OUTSIDE; // the nucleotide at the start of the active region -- used only in initialization
  seq1+=OBJ_SIZE-OUTSIDE-1; // the next nucleotide to enter the active area (ie what's under the nucleotide but not in the 'outside' excluded region

  bool initFlag=true;  // flag this when the pastWords variable needs to be refilled
  bool updateFlag=false; // set to false when you've just initialized pastWords
  bool writeFlag=false; // set to true when you want to write this line
  short_vector_vector myRules; // all the rules for this position: myRules[n-1] = vector of all rules for words length n
  short_vector rule1,rule2; 
  int n,i,ret;
  char c;
  double E;
  for (int seqN=1;seqN<=(int)seq.size()-OBJ_SIZE+1;seqN++) {
    if (seqN%10000==0) {
      printf("%d...\n",seqN);
    }

    if (seq1==seq.end()) {
      printf("found end of sequence inside seqN = %d\n",seqN);
      return -3;
    }
    
    if (initFlag) { 
      ret = initPastWords(seq,seq0,pastWords);
      if (ret!=1) {
	printf("there was some error (%d) in initPastWords\n",ret);
	return ret;
      }
      initFlag=false;
      updateFlag=false;
      writeFlag=true;
    }

    if (updateFlag) {
      //printf("%d) updating...\n",seqN);
      c = *seq1;
      for (n=N-1;n>0; n--) {      
	pastWords[n].push_back(catWord(pastWords[n-1].back(),c)); // next word is the new letter plus the old length n-1 word 
	pastWords[n].pop_front(); // drop the last word
      }
      pastWords[0].push_back(charMap[(int)c]);
      pastWords[0].pop_front();
    }
    
    if (writeFlag) {

      fprintf(oFile,"%15d",seqN);
      
      //int cnt;
      std::string s;
      
      if (strand >= 0) { // both or F
	// make E from rules
	E=coef[0][4][0]; // add intercept 
	n=1;
	for (iqvi = pastWords.begin();iqvi!=pastWords.end();++iqvi) {
	  i=0;
	  for (iqi = iqvi->begin(); 
	       iqi < iqvi->end(); 
	       iqi++) {
	    //printf("looking for coef[%d][%d][%d]\n",n-1,*iqi,i);
	    E += coef[n-1][*iqi][i];
	    i++;
	  }
	  n++;
	}	
	fprintf(oFile," %14e",E);
      }

      // reverse compliment
      if (strand <= 0) {
	E=coef[0][4][0]; // add intercept 
	n=1;
	for (iqvi = pastWords.begin();iqvi!=pastWords.end();++iqvi) {
	  i=OBJ_SIZE-OUTSIDE*2-n;
	  //i=0;
	  for (riqi = iqvi->rbegin(); riqi<iqvi->rend(); riqi++,i--) {
	    E += coef[n-1][*riqi][i];
	  }
	  n++;
	}	
	fprintf(oFile," %14e",E);
      }

      fprintf(oFile,"\n");
      updateFlag=true;
    }
    seq0++;
    seq1++;
  }
  fclose(oFile);

  return 0;
}


// fill up the pastWords vector, which pastWords up the number of all words of all lengths up to N underneath a nucleosome starting at seq1
// assumes you've checked filter already
int CspatRes::initPastWords(const char_vector& seq,char_vector::const_iterator seq1, int_queue_vector &pastWords) {

  // clear and reinitialize pastWords (STL queue has no clear method)
  pastWords.clear();
  int_queue sq; 
  for (int n = 1; n<=N; n++) {
    pastWords.push_back(sq);
  }
  int c,n;
  std::string s;
  string_vector::iterator si;
  for (int i=0;i<OBJ_SIZE-2*OUTSIDE;i++) {
    if (seq1==seq.end()) {
      printf("Got to end of sequence inside initCounts\n");
      return -1;
    }

    c = (int)*seq1;
    for (n=N-1;n>0; n--) {      
      if (i-n<0) { // must be at least n letters in before trying to make add a word of length n
	continue;
      }
      pastWords[n].push_back(catWord(pastWords[n-1].back(),c));
    }
    pastWords[0].push_back(charMap[c]);
    seq1++;    
  }

  return 1;
} 

int CspatRes::catWord(int w, char c) {
  return 4*w+charMap[(int)c];
}

void CspatRes::dumpPastWords(int_queue_vector pastWords) {
  
  int_queue_vector::iterator iqvi;

  short n=1;
  int i;
  std::string w0,w1,w;
  for (iqvi = pastWords.begin();iqvi!=pastWords.end();++iqvi) {
    i=1;
    printf("size of pastWords[%d] = %d\n",n,(int)iqvi->size());
    if (n==1) {
      w0 = "";
      while (!iqvi->empty()) {
	w = intToWord(n,iqvi->front());
	std::cout << w;
	iqvi->pop_front();
      }
    } else {
      while (!iqvi->empty()) {
	w = intToWord(n,iqvi->front());
	std::cout << w << " ";
	iqvi->pop_front();
	if (i%14==0) {
	  printf("\n");
	}
	i++;
      }
    }
    printf("\n");
    n++;
  }
  return;
}


	/*
	  myRules.clear();      
	  for (n=0;n<N;n++) {
	  rule1.clear();
	  rule1.resize((int)pow(4,n+1),0);
	  for (j=0;j<pow(4,n+1);j++) {
	  rule2 = rules[n][j];
	  scaleRule(counts[n][j],rule2);
	  addRules(rule1,rule2);
	  }
	  myRules.push_back(rule1);
	  }
  int_queue::reverse_iterator riqi;
  i=0;
  printf("hio\n");
  for (riqi = iq.rbegin(); riqi<iq.rend(); riqi++) {
    std::cout << "iq[" << i << "] = " << *riqi << "\n";
  }
  exit(0);

	*/
