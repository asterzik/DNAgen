#include "nonSpec.h"
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include "handyDandy.h"


CnonSpec::CnonSpec() {}

CnonSpec::CnonSpec(unsigned short n) {
  defineCharMap();
  init(n);
}

CnonSpec::CnonSpec(char epsFile[], coefficientMap &coef) {

  defineCharMap();
  unsigned short n;
  int ret = readCoef(epsFile,coef,n);
  if (ret!=1) {
    std::cerr << "unable to read eps file " << epsFile << "\n";
    exit(ret);
  }
  printf("found maximum word length %d in coef file %s\n",n,epsFile);
  init(n);
 
}

// seq = entire sequence
int CnonSpec::makeRTable(char outFile[],const char_vector& seq,const double_vector& energy,const short_vector& filter,bool db) {
  
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
  fprintf(oFile,"%14s","E");
  for (int n=1;n<=N;n++) {
    for (int i=0;i<pow(4,n);i++) {
      if (hasThrees(n,i)) {
	continue;
      }
      s = intToWord(n,i);
      s+= "1";
      const char *p = s.c_str();
      fprintf(oFile," %7s",p);
    }
  }
  fprintf(oFile,"\n");


  short_vector_vector counts; // counts[n-1][i] = counts of word i which has length n
  // initialize counts
  short_vector sv; 
  for (int n = 1; n<=N; n++) {
    sv.resize((int)pow(4,n),0);
    counts.push_back(sv);
  }

  string_vector frontWords; // the words at the 'front' (entering) end of the active region
  string_vector backWords;  // "              " 'back' (leaving)   "                      "                            
  string_vector::iterator wi; // iterate through the above
  char_vector::const_iterator seq1,seq0,seqAtStart; // seq1 reads the next, seq0 reads the end
                                              // the idea is that you'll add counts from the end and subtract counts from the beginning  
  
  short_vector::const_iterator f; // filter iterator;
  double_vector::const_iterator en; // energy iterator;

  seqAtStart = seq0 = seq1 = seq.begin(); // used exclusively to initialize counts. 
  f = filter.begin();
  en = energy.begin();
  seqAtStart+=OUTSIDE; // the nucleotide at the start of the active region -- used only in initialization
  seq0+=OUTSIDE+N-1; // the next nucleotide to enter the 'backWords' vector  
                     // for N=3,OUTSIDE=3: XXXACGTACGT...
                     //                pos=i    ^ *seq0==G at position i+OUTSIDE+N-1
  seq1+=OBJ_SIZE-OUTSIDE-1; // the next nucleotide to enter the active area (ie what's under the nucleotide but not in the 'outside' excluded region
  int seqN = 1; // unused variable
  bool initFlag=true;  // flag this when the counts variable needs to be refilled
  bool updateFlag=false; // set to false when you've just initialized counts
  bool writeFlag=false; // set to true when you want to write this line
  short_vector_vector myRules; // all the rules for this position: myRules[n-1] = vector of all rules for words length n
  short_vector rule1,rule2; 
  int n,j,w;
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
      int ret = initCounts(seq,seqAtStart,counts,frontWords,backWords);
      initFlag=false;
      updateFlag=false;
      if (ret!=1) {
	printf("there was some error (%d) in initCounts\n",ret);
	return ret;
      }
      writeFlag=true;
    } else if (*f==0) {
      updateFlag=false;
      initFlag=true;
      writeFlag=false;
    }

    if (updateFlag) {
      //printf("%d) updating...\n",seqN);
      n=0;
      for (wi=backWords.begin();wi!=backWords.end();++wi) {
	w = wordToInt(*wi);
	//std::cout << seqN << ") removing word " << *wi << " from beginning\n";
	counts[n][w]--;
	n++;
      }
      shiftBackWords(backWords,*seq0);
      shiftFrontWords(frontWords,*seq1);
      n=0;
      for (wi=frontWords.begin();wi!=frontWords.end();++wi) {
	//std::cout << seqN << ") adding word " << *wi << " to end\n";
	w = wordToInt(*wi);
	counts[n][w]++;
	n++;
      }
    }
    
    //printf("%d) counts dump: \n",seqN);
    //dumpCounts(counts);
    //printf("\n");
      
    if (writeFlag) {
      // calculate rules
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

      // print line of table
      if (db) { 
	fprintf(oFile,"%14d",(int)*en);
      } else {
	fprintf(oFile,"%14.6e",*en);
      }
      for (n=1;n<=N;n++) {
	for (int i=0;i<pow(4,n);i++) {
	  if (hasThrees(n,i)) {
	    continue;
	  }
	  fprintf(oFile," %7d",myRules[n-1][i]);
	}
      }
      fprintf(oFile,"\n");
	
      // reverse compliment
      myRules.clear();      
      for (n=0;n<N;n++) {
	rule1.clear();
	rule1.resize((int)pow(4,n+1),0);
	for (j=0;j<pow(4,n+1);j++) { // loop over words
	  //k = rc(j);
	  rule2 = rules[n][rc(n+1,j)]; // apply rule for rc(j) scaled by counts(j)
	  scaleRule(counts[n][j],rule2);
	  addRules(rule1,rule2);
	}
	myRules.push_back(rule1);
      }

      // print line of table
      if (db) {
	fprintf(oFile,"%14d",-(int)*en); 
      } else {
	fprintf(oFile,"%14.6e",*en);
      }
      for (n=1;n<=N;n++) {
	for (int i=0;i<pow(4,n);i++) {
	  if (hasThrees(n,i)) {
	    continue;
	  }
	  fprintf(oFile," %7d",myRules[n-1][i]);
	}
      }
      fprintf(oFile,"\n");
      
  
      updateFlag=true;
    }

    en++;
    f++;
    seqAtStart++;
    seq0++;
    seq1++;
    seqN++;
  }
  fclose(oFile);
  
  //printf("seqN = %d\n",seqN);
  return 1;
}
// given coefficients, make an energy landscape
int CnonSpec::makeEnergy(char outFile[],const char_vector& seq,coefficientMap coef,short strand) {
  
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

  short_vector_vector counts; // counts[n-1][i] = counts of word i which has length n
  // initialize counts
  short_vector sv; 
  for (int n = 1; n<=N; n++) {
    sv.resize((int)pow(4,n),0);
    counts.push_back(sv);
  }

  
  string_vector frontWords; // the words at the 'front' (entering) end of the active region
  string_vector backWords;  // "              " 'back' (leaving)   "                      "                            
  string_vector::iterator wi; // iterate through the above
  char_vector::const_iterator seq1,seq0,seqAtStart; // seq1 reads the next, seq0 reads the end
                                              // the idea is that you'll add counts from the end and subtract counts from the beginning  
  
  seqAtStart = seq0 = seq1 = seq.begin(); // used exclusively to initialize counts. 
  seqAtStart+=OUTSIDE; // the nucleotide at the start of the active region -- used only in initialization
  seq0+=OUTSIDE+N-1; // the next nucleotide to enter the 'backWords' vector  
                     // for N=3,OUTSIDE=3: XXXACGTACGT...
                     //                pos=i    ^ *seq0==G at position i+OUTSIDE+N-1
  seq1+=OBJ_SIZE-OUTSIDE-1; // the next nucleotide to enter the active area (ie what's under the nucleotide but not in the 'outside' excluded region

  bool initFlag=true;  // flag this when the counts variable needs to be refilled
  bool updateFlag=false; // set to false when you've just initialized counts
  bool writeFlag=false; // set to true when you want to write this line
  short_vector_vector myRules; // all the rules for this position: myRules[n-1] = vector of all rules for words length n
  short_vector rule1,rule2; 
  int n,j,w;
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
      int ret = initCounts(seq,seqAtStart,counts,frontWords,backWords);
      initFlag=false;
      updateFlag=false;
      if (ret!=1) {
	printf("there was some error (%d) in initCounts\n",ret);
	return ret;
      }
      writeFlag=true;
    }

    if (updateFlag) {
      //printf("%d) updating...\n",seqN);
      n=0;
      for (wi=backWords.begin();wi!=backWords.end();++wi) {
	w = wordToInt(*wi);
	counts[n][w]--;
	n++;
      }
      shiftBackWords(backWords,*seq0);
      shiftFrontWords(frontWords,*seq1);
      n=0;
      for (wi=frontWords.begin();wi!=frontWords.end();++wi) {
	w = wordToInt(*wi);
	counts[n][w]++;
	n++;
      }
    }
    
    if (writeFlag) {
      fprintf(oFile,"%15d",seqN);

      
      double c;
      int cnt;
      std::string s;
      
      if (strand >= 0) { // both or F
	myRules.clear();      
	// calculate rules
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
	
	// make E from rules
	E=coef[0][4][0]; // add intercept
	for (n=1;n<=N;n++) {
	  for (int i=0;i<pow(4,n);i++) {
	    if (hasThrees(n,i)) {
	      continue;
	    }
	    
	    c=coef[n-1][i][0];
	    cnt = myRules[n-1][i];
	    //std::cout << "forward: " << E << " += " << cnt << " * coef[" << s << "][0] = " << cnt << " * " << coef[s][0] <<" = " << cnt * coef[s][0] << "\n";
	    
	    E += myRules[n-1][i] * coef[n-1][i][0];
	  }
	}	
	fprintf(oFile," %14e",E);
      }

      // reverse compliment
      if (strand <= 0) {
	myRules.clear();      
	for (n=0;n<N;n++) {
	  rule1.clear();
	  rule1.resize((int)pow(4,n+1),0);
	  for (j=0;j<pow(4,n+1);j++) { // loop over words
	    //k = rc(j);
	    rule2 = rules[n][rc(n+1,j)]; // apply rule for rc(j) scaled by counts(j)
	    scaleRule(counts[n][j],rule2);
	    addRules(rule1,rule2);
	  }
	  myRules.push_back(rule1);
	}
	
	// make E from rules
	E=coef[0][4][0]; // add intercept
	for (n=1;n<=N;n++) {
	  for (int i=0;i<pow(4,n);i++) {
	    if (hasThrees(n,i)) {
	      continue;
	    }
	    
	    c=coef[n-1][i][0];
	    cnt = myRules[n-1][i];
	    //std::cout << "reverse: " << E << " += " << cnt << " * coef[" << s << "][0] = " << cnt << " * " << coef[s][0] <<" = " << cnt * coef[s][0] << "\n";
	    
	    E += myRules[n-1][i] * coef[n-1][i][0];
	  }
	}	
	fprintf(oFile," %14e",E);
      }

      fprintf(oFile,"\n");
      updateFlag=true;
    }
    seqAtStart++;
    seq0++;
    seq1++;
  }
  fclose(oFile);

  return 0;
}

void CnonSpec::shiftFrontWords(string_vector &words, char next) {
  int n;
  for (n=words.size()-1;n>=1;n--) {
    words[n] = words[n-1]+next;
  }
  words[0] = next;
  return;
}

void CnonSpec::shiftBackWords(string_vector &words, char next) {
  int n;
  for (n=0;n<(int)words.size()-1;n++) {
    words[n] = words[n+1].substr(1,n+1);
  }
  std::string w = words[n].substr(1,n);
  w+=next;
  words[n]=w;
  return;
}


// fill up the counts vector, which counts up the number of all words of all lengths up to N underneath a nucleosome starting at seq1
// assumes you've checked filter already
int CnonSpec::initCounts(const char_vector& seq,char_vector::const_iterator seq1, short_vector_vector &counts,string_vector &frontWords,string_vector &backWords) {

  frontWords.clear();
  backWords.clear();
  counts.clear();
  short_vector sv;
  for (int n = 1; n<=N; n++) {
    sv.resize((int)pow(4,n),0);
    counts.push_back(sv);
  }

  int w,j,n;
  std::string s;
  string_vector::iterator si;
  for (int i=0;i<OBJ_SIZE-2*OUTSIDE;i++) {
    if (seq1==seq.end()) {
      printf("Got to end of sequence inside initCounts\n");
      return -1;
    }
    if (i==0) {
      s = *seq1;
      backWords.push_back(s);
      frontWords.push_back(s);
      w = wordToInt(s);
      counts[0][w]++;
    } else if (i<N) {
      
      s = backWords[i-1] + *seq1;
      backWords.push_back(s);
      frontWords.push_back(s);

      for (j=i-2;j>=0;j--) {
	s = frontWords[j] + *seq1;
	frontWords[j+1] = s;
      }
      frontWords[0] = *seq1;
      n=0;
      for (si=frontWords.begin(); si!=frontWords.end();++si) {
	w = wordToInt(*si);
	counts[n][w]++;
	n++;
      }
      //frontWords[0]+=*seq1;
      /*
	next = a
	fwrds = (g,cg,acg)
	
	want fwrds = (a,ga,cga)
       */
      

    } else { // now front and backwords are filled, you can concentrate on counts
      shiftFrontWords(frontWords,*seq1);
      n=0;
      for (si=frontWords.begin(); si!=frontWords.end();++si) {
	w = wordToInt(*si);
	counts[n][w]++;
	n++;
      }
    }
    seq1++;    
  }

  return 1;
} 

void CnonSpec::dumpCounts(short_vector_vector counts) {
  
  short_vector_vector::iterator svvi;
  short_vector::iterator svi;      
  short n=1;
  int w;
  for (svvi = counts.begin();svvi!=counts.end();++svvi) {
    w=0;
    for (svi=svvi->begin();svi!=svvi->end();++svi) {
      std::cout << intToWord(n,w) << ": " << *svi << "\t";
      w++;
    }
    printf("\n");
    n++;
  }
  return;
}

