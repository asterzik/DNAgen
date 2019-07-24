#include "model.h"
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include "handyDandy.h"

Cmodel::Cmodel() {}

Cmodel::Cmodel(unsigned short n) {
  defineCharMap();
  init(n);
}

// required instantiator for child classes
// create a model according to epsFile, store coeff in coef
Cmodel::Cmodel(char epsFile[], coefficientMap &coef) {
  
  std::cerr << "**************Cmodel(char epsFile[], coefficientMap &coef) should never execute\n";

  /*
  int ret = readCoef(epsFile,coef);
  if (ret!=1) {
    printf("unable to read eps file %s\n",epsFile);
    exit(ret);
  }

  std::string s;
  coefficientMap::iterator it1;
  int n=-1;
  for (it1=coef.begin(); it1!=coef.end(); ++it1) {
    s = it1->first;
    if (s.compare("Intercept") && (int)s.length()>n) {
      n = s.length();
    }
  }

  init(n);
  */
}
Cmodel::~Cmodel() {}

void Cmodel::defineCharMap() {
  /*
  // NOTE: reverse compliments must be separated by two.  rc(word) uses this.

  //  TEMPORARY MAP COLLAPSING ON _A_ INSTEAD OF _C_
  for (int i=0;i<128;i++) {
    switch (i) {
    case 65: // A
    case 97: // a
      charMap[i] = 3;
      break;
    case 67: // C
    case 99: // c
      charMap[i] = 0;
      break;
    case 71:  // G
    case 103: // g
      charMap[i] = 2;
      break;
    case 84:  // T
    case 116: // t
      charMap[i] = 1;
      break;
    default:
      charMap[i] = -1;
    }
  }


  shortMap[0] = 'C';
  shortMap[1] = 'T';
  shortMap[2] = 'G';
  shortMap[3] = 'A';
  // END TMP MAP
  */
  
  // PROTECTED MAP
  for (int i=0;i<128;i++) {
    switch (i) {
    case 65: // A
    case 97: // a
      charMap[i] = 0;
      break;
    case 67: // C
    case 99: // c
      charMap[i] = 3;
      break;
    case 71:  // G
    case 103: // g
      charMap[i] = 1;
      break;
    case 84:  // T
    case 116: // t
      charMap[i] = 2;
      break;
    default:
      charMap[i] = -1;
    }
  }

  shortMap[0] = 'A';
  shortMap[1] = 'G';
  shortMap[2] = 'T';
  shortMap[3] = 'C';
  //END PROTECTED MAP

}


void Cmodel::init(unsigned short n) {

  N=n;

  OBJ_SIZE = 147;
  OUTSIDE = 3;

  rules = makeRules();
  


  //DEBUG
  /*
    // debug vector addition and multiplication
  short_vector v1=rules[0][1];
  short_vector v2=rules[0][2];
  scaleRule(10,v1);
  scaleRule(12,v2);
  int ret = addRules(v1,v2);
  if (ret!=1) {
    printf("couldn't add v1 and v2\n");
    return;
  }
  short_vector::iterator it2;
  int j=0;
  std::cout << "added 10*" << intToWord(1,1) << " to 12*" <<intToWord(1,2)  << ":\n";
  for (it2=v1.begin();it2!=v1.end();++it2) {
    std::cout << "\t" << intToWord(1,j) << "->"<< *it2;
    j++;

  }
  printf("\n");
  */

  /*
  // debug rules
  short_vector_vector svv = rules[1];
  short_vector_vector::iterator it1;
  short_vector::iterator it2;
  int i,j;
  i=0;
  for (it1=svv.begin();it1!=svv.end();++it1) {
    j=0;
    std::cout << "rules for " << intToWord(2,i) << ":";
    for (it2=it1->begin();it2!=it1->end();++it2) {
      std::cout << "\t" << intToWord(2,j) << "->"<< *it2;
      j++;
    }
    printf("\n");
    i++;
  }
  */

  /*
  printf("about to debug reverse compliment\n");
  //debug reverse complement
  std::string w;
  int n,r;
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
      for (int k=0;k<4;k++) {
	w.clear();
	w+=shortMap[i];
	w+=shortMap[j];
	w+=shortMap[k];
	n = wordToInt(w);
	r = rc(3,n);
	std::cout << "rc(" << w << ") = " << intToWord(3,r) << "\n";
      }
    }
  }
  //exit(1);
  */
}



// coef format:
// coef[AA123] = coef[1][0][122]
//                    n length of word
//                       w wordToInt(AA)           
//                          l position of word - 1
// store intercept at coef[0][4][0]
int Cmodel::readCoef(char inFile[],coefficientMap &coef,unsigned short& maxN) {

  coef.clear();

  maxN=0; // maximum word length;
  double intercept;

  std::string line;
  std::ifstream myFile;
  myFile.open(inFile);
  
  if (!myFile.is_open()) {
    printf("failed to open infile %s\n",inFile);
    return -1;
  }
  
  std::stringstream st1(std::stringstream::in | std::stringstream::out);
  string_vector spl; // split line at white space
  string_vector::iterator sv_it; 

  int tokenLoc, epsLoc,i;
  epsLoc=tokenLoc=-1;
   // first find header
  bool noheader = true;
  while (noheader && ! myFile.eof()) {
    getline (myFile,line);
    if (line[0]=='#') {
      continue;
    }    
    if (line.length()==0) {
      continue;
    }
    spl = hand_splitWhiteSpace(line);    
    i=0;
    for (sv_it = spl.begin(); sv_it!=spl.end(); ++sv_it) {
      //std::cout << "comparing spl[" << i << "] to eps :" << (*sv_it).compare(
      if (((*sv_it).compare("eps"))==0) {
	epsLoc = i;
      } else if ((*sv_it).compare("token")==0 || (*sv_it).compare("E")==0) {
	tokenLoc = i;
      }
      i++;
    }
    noheader = false;
  }
  //printf("(epsLoc,tokenLoc) = (%d, %d)\n",epsLoc,tokenLoc);

  if (myFile.eof()) {    
    printf("couldn't find header\n");
    return -2;
  } else if (tokenLoc<0) {
    printf("couldn't find token column\n");
    return -2;
  } else if (epsLoc<0) {
    printf("couldn't find eps column\n");
    return -2;
  }
  unsigned int minElems = 1+(epsLoc>tokenLoc)?epsLoc:tokenLoc; // need at least this many elements per line
  
  double_vector v;
  std::string token,s;
  unsigned short tokenL; // length of token
  int tokenI; // int representing word
  unsigned short n; // which zone/position (eg aa1 -> l=0, aa123 -> l=122)
  double eps;
  while (! myFile.eof() ) {
    getline (myFile,line);
    if (line[0]=='#') {
      continue;
    }
    if (!line.length()) {
      continue;
    }
    spl = hand_splitWhiteSpace(line);
    if (spl.size()<minElems) {
      continue;
    }

    s = spl[tokenLoc];
    st1.str(spl[epsLoc]);
    st1 >> eps;
    if (st1.fail()) {
      std::cerr << "couldn't parse eps from " << st1.str() << ", got error " << st1.fail() << ", badbit state: " << st1.bad() << "\n";
      return -3;
    }
    st1.clear();
    if (s.compare("Intercept")==0) {
      intercept = eps;
      continue;
    } 
    
    n = -1;
    i=0;
    while (!isdigit(s[i])&&i<(int)s.length()) {
      i++;
    }
    token = s.substr(0,i);
    st1.str(s.substr(i,s.length()-i));
    st1 >> n;
    if (st1.fail()) {
      std::cerr << "couldn't get a position from from " << st1.str() << ", got error " << st1.fail() << ", badbit state: " << st1.bad() << "\n";
      return -3;
    }
    if (n<0) {
      std::cerr << "couldn't find a position for token " << s << "\n";
      return -3;
    }
    st1.clear();

    tokenL = token.length();
    if (maxN<tokenL) {
      maxN=tokenL;
    }
    tokenI = wordToInt(token);

    if (coef.size()<tokenL) { 
      coef.resize(tokenL);
    }
    if ((int)coef[tokenL-1].size() < tokenI+1) {
      coef[tokenL-1].resize(tokenI+1);
    }
    if (coef[tokenL-1][tokenI].size()) { // if this vector exists      
      if ((int)coef[tokenL-1][tokenI].size()<n) {
	coef[tokenL-1][tokenI].resize(n,0);
      }      
      coef[tokenL-1][tokenI][n-1]=eps;
    } else {
      v.clear();
      v.resize(n);
      v[n-1] = eps;
      coef[tokenL-1][tokenI]=v;
    }
  } // end read file loop
  myFile.close();

  tokenI = 4;
  if (coef[0].size()<4) {  // we expect coef[0][0,1,2] to be defined and no other in coef[0]
    coef[0].resize(4);
  }
  v.clear();
  v.push_back(intercept);
  coef[0].push_back(v);
  

  return 1;
}

/*
  take a list of sequences, energies, filters, and make R-tables out of all of them.                                      
  the tables concatenate the inputs even as they break them up.  The idea is to do a regression on all your data without  
  overflowing RAM with R-Tables that are too big.                                                                         
  Args are lists of <TYPE> where TYPE is as for makeRTable                                                                
*/
int Cmodel::makeCombinedRTable(char baseOutFile[],const char_vector_vector& seqs,const double_vector_vector& energies,const short_vector_vector& filters,size_t maxBp) {
  
  /*
    method: 
    concatenate all the data in the lists (flatten the vector_vectors into vectors)
    break up the flattened vectors into sections of the desired length and send them to makeRTable()
  */


  if (seqs.size()!=energies.size() || seqs.size()!=filters.size()) {
    printf("makeCombinedRTable tried to combine sets of unequal size\n");
    return -1;
  }
  
  double_vector blanks;
  unsigned short blanksSize = OBJ_SIZE-1;
  blanks.resize(blanksSize,0);

  char_vector * bigSeq = new char_vector;
  double_vector * bigEn = new double_vector;
  short_vector * bigFilter = new short_vector;
  char_vector_vector::const_iterator cvvi;
  double_vector_vector::const_iterator dvvi;
  short_vector_vector::const_iterator svvi;

  for (cvvi=seqs.begin();cvvi<seqs.end();cvvi++) {
    bigSeq->insert(bigSeq->end(),cvvi->begin(),cvvi->end());
  }

  size_t_vector enEnds;
  size_t_vector chrEnds;
  for (dvvi=energies.begin();dvvi<energies.end();dvvi++) {
    bigEn->insert(bigEn->end(),dvvi->begin(),dvvi->end());
    enEnds.push_back(bigEn->size());
    bigEn->insert(bigEn->end(),blanks.begin(),blanks.end());
    chrEnds.push_back(bigEn->size());
  }
  
  for (svvi=filters.begin();svvi<filters.end();svvi++) {
    bigFilter->insert(bigFilter->end(),svvi->begin(),svvi->end());
  }
  
  size_t totalSize = chrEnds.back();

  if (bigFilter->size()!=totalSize) {
    printf("failure to match sizes of energies and filters\n");
    return -1;
  }
  if (bigSeq->size()!=totalSize) {
    printf("failure to match sizes of energies and sequences\n");
    return -1;
  }
  
  // make sure to filter the ends of each chromosome, as the energies from End to End-OBJ_SIZE+1 are junk
  size_t_vector::iterator stvi1,stvi2;
  for (stvi1 = enEnds.begin(),stvi2 = chrEnds.begin(); stvi1 < enEnds.end() && stvi2 < chrEnds.end(); stvi1++,stvi2++) {
    for (size_t i=*stvi1;i < *stvi2; i++) {
      (*bigFilter)[i] = 0;
    }
  }

  size_t segSize = totalSize;
  int nsegs = 1;
  while (segSize>maxBp) {
    nsegs++;
    segSize = totalSize/nsegs;
  }
  segSize++;  //  adding this line seems to prevent the last segment from having length<10
  while (nsegs*segSize<totalSize) { 
    nsegs++; 
  }

  // format file-names
  // given prefix.body.tab, 
  // files to appear as prefix.body.XXXXX-YYYYY.tab 
  //   the number of digits should remain fixed so that they ls better
  short totalDigits = 1+log(totalSize)/log(10);
  string_vector base = hand_split(".",baseOutFile);
  std::string suffix = base.back();
  base.pop_back();

  // now split up the data and make tables!
  char_vector * seq = new char_vector;
  double_vector * en = new double_vector;
  short_vector * filter = new short_vector;

  unsigned int segStart = 1;
  unsigned int segStop = segStart+segSize-1;
  string_vector out;
  std::string numberString;
  std::string outFile = "";
  std::stringstream st1(std::stringstream::in | std::stringstream::out);
  short nDigits;
  char file[199];
  for (int i=0;i<nsegs;i++) {
    // format the output file (lengthy code...)
    numberString = "";
    nDigits = 1+log(segStart)/log(10);    
    while (nDigits<totalDigits) {
      numberString+="0";
      nDigits++;
    }
    st1.str("");
    st1 << segStart;
    numberString+=st1.str();
    numberString+="-";
    nDigits = 1+log(segStop)/log(10);    
    while (nDigits<totalDigits) {
      numberString+="0";
      nDigits++;
    }
    st1.str("");
    st1 << segStop;
    numberString+=st1.str();
    
    out = base;    
    out.push_back(numberString);
    out.push_back(suffix);
    outFile = hand_join(".",out);
    std::cout << outFile << "\n";
    strcpy(file,outFile.c_str());

    // excerpt the data
    seq->clear();
    seq->insert(seq->begin(),bigSeq->begin()+segStart-1,bigSeq->begin()+segStop-1+OBJ_SIZE-1);
    en->clear();
    en->insert(en->begin(),bigEn->begin()+segStart-1,bigEn->begin()+segStop-1);
    filter->clear();
    filter->insert(filter->begin(),bigFilter->begin()+segStart-1,bigFilter->begin()+segStop-1);
    this->makeRTable(file,*seq,*en,*filter,false);
    

    segStart+=segSize;
    segStop+=segSize;
    if (segStop>totalSize) {
      segStop = totalSize;
    }
  }
  delete bigSeq;
  delete bigEn;
  delete bigFilter;
  delete seq;
  delete en;
  delete filter;

  return 0;

}

/***********************//**********//***********/
/*  protected methods ****  protected methods ***/
/***********//*********//************************/

// base 4 map AA = 0, AG = 1, AT = 2,AC=3, GA = 4, and so forth
int Cmodel::wordToInt(std::string w) {
  int n = (int)w.length();
  int ret=0;
  int j;
  for (int i=0;i<n;i++) {
    j = charMap[(int)w[i]];
    if (j==-1) {
      std::cerr << "(int)w[" << i <<"] == " << (int)w[i]<< " is unidentified nucleotide in word " << w << "\n";
      return -1;
    }
    ret+=j*pow(4,n-1-i);
  }
  return ret;
}

std::string Cmodel::intToWord(short n,int w) {
  std::string ret;
  
  int j;
  for (int i=n-1;i>=0;i--) {
    j = w/pow(4,i);
    j = j%4;
    ret+=shortMap[j];
  }

  return ret;
}

// return a 3d array with all the rules up to a given N (see notes below)
Trules Cmodel::makeRules() {

  /*
    notes on Trules
    Trules[n-1] = all the rules for words of length n
    Trules[n-1][i] = all the rules for word i (where i = wordToInt[i])
    Trules[n-1][i][j] = the specific rule for mapping word i to word j
    so, word 'C' is Trules[0][3] = (-1,-1,-1)
    rules for 'CA' to GA is TRules[1][12][4]=-1
    note that rules from XX to YZ where Y and/or Z is a C are always 0!!
      could have used base three to index the words mapped to 
      but that would have been more confusing
  */

  Trules ret;
  short_vector_vector svv;  
  short_vector rule; // rule for a specific word
  short_vector threes; // consider i as base 4.  which digits are 3s ('C's)?
  short n;
  int i,j,cCount;
  //bool jHasThree; // if you find a 3 in j then set flag
  for (n=0;n<N;n++) {
    svv.clear();    
    for (i=0;i<pow(4,n+1);i++) {
      rule.clear();
      threes.clear();
      cCount = findThrees(n,i,threes);
      for (j=0;j<pow(4,n+1);j++) {
	if (hasThrees(n,j)) {
	  rule.push_back(0);
	} else {
	  //printf("(n,i,j) = (%d,%d,%d)\t",n,i,j);
	  //std::cout << "cCount for " << intToWord(n+1,i) << " and " << intToWord(n+1,j) << " is " << cCount << "\n";
	  if (threesMatch(i,j,threes)) {
	    if (cCount%2==0) {
	      rule.push_back(1);
	    } else {
	      rule.push_back(-1);
	    }
	    //svv[i][j]=-1;
	  } else if (i==j) {
	    rule.push_back(1);
	    //svv[i][j]= 1;
	  } else {
	    rule.push_back(0);
	    //svv[i][j]= 0;
	  }
	  // note: if you don't push_back *each time* then the whole thing breaks
	} // end if no flag
      } // end j loop
      svv.push_back(rule);
    } // end i loop
    ret.push_back(svv);
  } // end n loop
  //exit(1);
  return ret;
}

// if i has threes and j matches i's non-threes, then return the number of c's in i; else 0
// for instance:  AAC and AAT "match", and will return 1;
//                ACC and ATG "match", and will return 2;
//                ACC and TTG don't match, return 0;
bool Cmodel::threesMatch(int i, int j,short_vector threes) {
  
  int N = threes.size(); // size of word
  
  short_vector i_w,j_w;
  int x,k;
  for (k=N-1;k>=0;k--) {
    x = i/pow(4,k);
    i_w.push_back(x%4);
    x = j/pow(4,k);
    j_w.push_back(x%4);    
  }

  bool threePresent = false;
  for (k=0;k<N;k++) {
    if (!threes[k]&&i_w[k]!=j_w[k]) { // if two non-C letters don't match, return false
      return false;
    }
    if (threes[k]==1) { // only return true if i actually has threes
      threePresent = true;
    }
  }
  return threePresent;
}

// return true iff there is at least one three
bool Cmodel::hasThrees(short n,int i) {
  short k,x;
  for (k=n;k>=0;k--) {	
    x = i/pow(4,k);
    x = x%4;
    if (x==3) {
      return true;
    } 
  }
  return false;
}

int Cmodel::findThrees(short n,int i,short_vector &threes) {
  short k,x;
  int ret=0; // count of threes in word
  for (k=n;k>=0;k--) {	
    x = i/pow(4,k);
    x = x%4;
    if (x==3) {
      ret++;
      threes.push_back(1);
    } else {
      threes.push_back(0);
    }
  }
  return ret;
}

// vector addition
int Cmodel::addRules(short_vector &rule1,short_vector rule2) {
  if (rule1.size()!=rule2.size()) {
    printf("tried to add together two rules of different sizes\n");
    return -1;
  }
  short_vector::iterator i1;
  short_vector::iterator i2;
  for (i1=rule1.begin(),i2=rule2.begin();i1!=rule1.end();++i1,++i2) {
    (*i1)+=*i2;
  }
  return 1;
}

// vector scalar multiplication
void Cmodel::scaleRule(short s,short_vector &rule) {
  short_vector::iterator i;
  for (i=rule.begin();i!=rule.end();++i) {
    (*i) = (*i) * s;
  }  
}

int  Cmodel::rc(short n, int i) {
  int ret = 0;
  int x;
  //printf("rc(%d):\n",i);
  for (short j=0;j<n;j++) {
    x = i/pow(4,j);
    //printf("\t%dth place (%d)",j,x);
    x = (x+2)%4; // (A,G,T,C) -> (T,C,A,G) shift by 2
    //printf(" becomes (%d)",x);
    x = x*pow(4,n-j-1);
    //printf(", which goes into the %dth place, with value %d\n",n-j,x);
    ret+=x;
  }



  return ret;
}

void Cmodel::fillCoef(coefficientMap &coef) {

  int i,j,x; // i and j are words, x is a position index eg AG123 -> coef[i][x] -> coef[2][122]
  //std::string intToWord(short n, int i); // base 4 map AA = 0, AG = 1, AT = 2, TA = 4, and so forth note: n is one-based!
  std::string s;

 
  /*
    proof of principle
  dv.resize(10,0);
  j=10;
  for (dvi = dv.begin(); dvi < dv.end(); dvi++) {
    *dvi += j;
    j++;
  }
  j=0;
  for (dvi = dv.begin(); dvi < dv.end(); dvi++) {
    std::cout << "dv[" << j << "] = " << *dvi << "\n";	
    j++;
  }

  return;
  */

  size_t * sizes = new size_t[N];
  for (i = 0; i<N;i++) {
    sizes[i] = 0;
  }

  coefficientMap::iterator cmi;
  double_vector_vector::iterator dvvi;
  double_vector::iterator dvi;
  //double_vector dv;
  short_vector rule;
  short n=1;
  for (cmi = coef.begin();cmi<coef.end();cmi++,n++) {
    i=0;
    if (n>1) { // n=1 already resized for intercept at coef[0][4]
      cmi->resize(pow(4,n));
    } 
    for (dvvi = cmi->begin(); dvvi<cmi->end();dvvi++,i++) {
      //for (i = 0;i<pow(4,n);i++) {
      //dv = *dvvi;
      if (sizes[n-1]==0) { 
	// coef will have no vector for words with C's, and we are making those
	// we have to know what size the vector should be
	// rely on fact that C's are defined as coming at the end of this iteration (3's)
	sizes[n-1] = dvvi->size();
      }
      if (! hasThrees(n,i)) { // you only need to fill words with C's (threes)
	continue;
      }
      dvvi->resize(sizes[n-1],0);
      rule = rules[n-1][i];      
      
      x=0; // position index
      //std::cout << "filling " << intToWord(n,i) << "\n";
      for (dvi = dvvi->begin(); dvi < dvvi->end(); dvi++,x++) {	
      //for (dvi = dv.begin(); dvi < dv.end(); dvi++,x++) { 
	
	for (j = 0;j<pow(4,n);j++) {
	  if (hasThrees(n,j)) {
	    continue;
	  }
	  /*
	  if (x==15) {
	    std::cout << "\t + " << rule[j] << " * coef[" << intToWord(n,j) << "][" << x << "]\n";
	  }
	  */

	  *dvi += coef[n-1][j][x] * rule[j];
	}
	//exit(0);
      }
      //coef[n-1][i] = dv;
    }
  }
  
}

void Cmodel::dumpCoef(coefficientMap coef) {

  int i,j;
  //std::string intToWord(short n, int i); // base 3 map AA = 0, AG = 1, AT = 2, GA = 3, and so forth note: n is one-based!
  std::string s;
  double_vector dv;
  double_vector::iterator dvi;
  for (short n=1;n<=N;n++) {
    for (i = 0;i<pow(4,n);i++) { 
      s = intToWord(n,i);
      if ((int)coef[n-1].size()<(i+1)) {
	std::cerr << "can't find coef [" << n-1 << "]->[" << i << "]\n";
	exit(1);
      }
      dv = coef[n-1][i];
      j=0;
      std::cout << s << std::endl;
      for (dvi = dv.begin(); dvi < dv.end(); dvi++) {
	std::cout << "\tcoef[" << n-1 << "]->[" << i << "]->[" << j << "] = " << *dvi << "\n";	
	j++;
	if (0 && j>10) {
	  std::cout << "\t...\n";
	  break;
	}
      }
    }
  }
  return;
}

