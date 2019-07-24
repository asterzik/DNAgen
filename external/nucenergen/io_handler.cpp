#include "io_handler.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include "handyDandy.h"

Cio_handler::Cio_handler() {

  // fill charMap
  for (int i=0;i<128;i++) {
    switch (i) {
    case 65: // A
    case 97: // a
      charMap[i] = 'A';
      break;
    case 67: // C
    case 99: // c
      charMap[i] = 'C';
      break;
    case 71:  // G
    case 103: // g
      charMap[i] = 'G';
      break;
    case 84:  // T
    case 116: // t
      charMap[i] = 'T';
      break;
    default:
      charMap[i] = 'X';
    }
  }

  /*
  // fill shortMap
  shortMap[0] = 'A';
  shortMap[1] = 'G';
  shortMap[2] = 'T';
  shortMap[3] = 'C';
  */
}

/*
  utilities
*/

// second arg is pass by ref; return -1 if fail
// read .occ file
//   probabilities are filled into second arg
//   this version differs because it gets the input occupancy as well;
int Cio_handler::readOcc(char inFile[],long double prob0[], long double occ0[], int &seqLen) {

  std::string line;
  std::ifstream myFile;
  myFile.open(inFile);
  
  if (!myFile.is_open()) {
    printf("failed to open infile %s\n",inFile);
    return -1;
  }

  int len; // length of line
  int i;   // loop over line
  int seqN;         // where in sequence are we
  long double prob; // prob at seqN
  long double occ;  // occ  at seqN
  std::string first;     // convert to seqN, read from file
  std::string second;    // convert to prob, read from file
  std::string third;     // convert to occ,  read from file
  std::stringstream st1(std::stringstream::in | std::stringstream::out);
  int count; // count first or second
  seqLen = 0;

  while (! myFile.eof() ) {

    getline (myFile,line);
    if (line[0]=='#') {
      continue;
    }
    len = line.length();
    if (!len) {
      continue;
    }
    first = "";
    second = "";
    third = "";
    count = 1;
    for (i=0; i<len; i++) {
      if (line[i]!=' ') {
	if (count==1) {
	  first += line[i];
	} else if (count==2) {
	  second +=line[i];
	} else if (count==3) {
	  third +=line[i];
	}
      } else if (count==1 && first.length()) {
	count++;
      } else if (count==2 && second.length()) {
	count++;
      } else if (count==3 && third.length()) {
	count++;
      }
    } 
    //std::cout << "first = \'" << first << "\', second = \'" << second << "\', third \'" << third << "\'\n";
    if (!first.length() || !second.length() || !third.length()) {
      continue;
    } else if (first.compare("seqN")==0) {
      printf("found header\n");
      continue;
    }
    
    st1.str(first);
    st1 >> seqN;
    if (st1.fail()) {
      std::cout << "couldn't parse first string " << st1.str() << "\nQutting\n";
      printf("seqN = %d, prob = %Le\n", seqN,prob);
      return -1;
    }
    st1.clear();
    
    st1.str(second);
    st1 >> prob;
    if (st1.fail()) {
      std::cout << "couldn't parse second string " << st1.str() << "\nQutting\n";
      printf("seqN = %d, prob = %Le\n", seqN,prob);
      return -1;
    }
    st1.clear();
    
    st1.str(third);
    st1 >> occ;
    if (st1.fail()) {
      std::cout << "couldn't parse third string " << st1.str() << "\nQutting\n";
      printf("seqN = %d, prob = %Le\n", seqN,prob);
      return -1;
    }
    st1.clear();
    //printf("got:   seqN = %d, prob = %Le, occ = %Le\n\n",seqN,prob,occ);
    
    if (seqN<MAX_LENGTH && seqN>=0) {
      prob0[seqN-1] = prob;
      occ0[seqN-1] = occ;
    } else {
      std::cout << "trying to fill sequence at position " << seqN << ".  Fail.\n";
      return -1;
    }
    
    if (seqN>seqLen) {
      seqLen = seqN;
    }
  } // end read file loop

  myFile.close();

  return 1;
} // end read occ

// mirrors above, but doesn't bother to record occupancy
int Cio_handler::readOcc(char inFile[],long double prob0[], int &seqLen) {
  long double dummyOcc[MAX_LENGTH];
  return this->readOcc(inFile,prob0,dummyOcc,seqLen);
}

// print occupancy to stdout
void Cio_handler::dumpOcc(int seqLen,long double prob[],long double pL[]) {

  printf("%7s  %12s  %12s\n","seqN","prob","pL");
  for (int i=0; i<seqLen; i++) {
    printf("%7i  %11Le  %11Le\n",i,prob[i],pL[i]);
  }

  return;
}

// write occupancy to file
int Cio_handler::writeOcc(char outFile[], int seqLen, long double prob[], long double pL[]) {
  
  std::ofstream oFile(outFile);
  if (!oFile.is_open()) {
    printf("failed to open file %s\n",outFile);
    return -1;
  }

  
  std::string header = "# Bw =   X.XX, bkgr_model = none # note this is a dummy occupancy file\n# <E(1)> =  XXXXX, var(E(1)) =  xxxxx\n    seqN      Prob1       Occ1\n";
    
  oFile << header;

  char line[81];
  for (int i=0; i<seqLen; i++) {
    sprintf(line,"%8d   %12.8Le   %12.8Le\n",(i+1),prob[i],1-pL[i]);
    oFile << line;
  }
  
  oFile.close();
  return 1;
} // end write occ


// read .en file
// second arg is pass by ref; return -1 if fail  
//// note that for most purposes the following line is required:
////       seqLen = seqLen + OBJ_SIZE -1;  
int Cio_handler::readEn(char inFile[],double_vector * en0) {

  std::string line;
  std::ifstream myFile;
  myFile.open(inFile);  
  if (!myFile.is_open()) {
    printf("failed to open energy file %s\n", inFile);
    return -1;
  }

  int i;   // loop over line
  int seqN;         // where in sequence are we
  int prevSeqN=0;   // where we just were - throw an exception when this and prev are not identical
  double en;        // en at seqN
  int seqNLoc= -1;  // where in line is seqN?
  int enLoc= -1;    // where in line is en?
  std::stringstream st1(std::stringstream::in | std::stringstream::out);
  string_vector spl; // split line at white space
  string_vector::iterator sv_it; 
  
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
      //std::cout << "comparing spl[" << i << "] to seqN :" << (*sv_it).compare(
      if (((*sv_it).compare("seqN"))==0) {
	seqNLoc = i;
      } else if ((*sv_it).compare("en")==0 || (*sv_it).compare("E")==0) {
	enLoc = i;
      }
      i++;
    }
    noheader = false;
  }

  if (myFile.eof()) {    
    printf("couldn't find header\n");
    return -2;
  } else if (enLoc<0) {
    printf("couldn't find energy column\n");
    return -2;
  } else if (seqNLoc<0) {
    printf("couldn't find seqN column\n");
    return -2;
  }
  unsigned int minElems = 1+(seqNLoc>enLoc)?seqNLoc:enLoc; // need at least this many elements per line
  
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

    st1.str(spl[seqNLoc]);
    st1 >> seqN;
    if (st1.fail()) {
      std::cout << "couldn't parse first string " << st1.str() << "\n";
      printf("seqN = %d, en = %e\n", seqN,en);
      return -1;
    }
    st1.clear();
    
    st1.str(spl[enLoc]);
    st1 >> en;
    if (st1.fail()) {
      std::cout << "couldn't parse second string " << st1.str() << "\nQutting\n";
      printf("seqN = %d, en = %e\n", seqN,en);
      return -1;
    }
    st1.clear();
    //printf("goODt:   seqN = %d, en = %Le\n\n",seqN,en);
    
    if (seqN==prevSeqN+1) {
      en0->push_back(en);
      prevSeqN = seqN;
    } else {
      std::cout << "trying to fill sequence at position " << seqN << " with prevSeqN = " << prevSeqN << "  Fail.\n";
      return -1;
    }
    
  } // end read file loop
  myFile.close();

  return 1;
} // end read en

// read .en file
// second arg is pass by ref; return -1 if fail  
//// note that for most purposes the following line is required:
////       seqLen = seqLen + OBJ_SIZE -1;  
int Cio_handler::readFilter(char inFile[],short_vector * filt) {

  std::string line;
  std::ifstream myFile;
  myFile.open(inFile);  
  if (!myFile.is_open()) {
    printf("failed to open energy file %s\n", inFile);
    return -1;
  }

  int i;   // loop over line
  int seqN;         // where in sequence are we
  int prevSeqN=0;   // where we just were - throw an exception when this and prev are not identical
  short f;        // en at seqN
  int seqNLoc= -1;  // where in line is seqN?
  int filtLoc= -1;    // where in line is en?
  std::stringstream st1(std::stringstream::in | std::stringstream::out);
  string_vector spl; // split line at white space
  string_vector::iterator sv_it; 
  
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
      //std::cout << "comparing spl[" << i << "] to seqN :" << (*sv_it).compare(
      if (((*sv_it).compare("seqN"))==0) {
	seqNLoc = i;
      } else if ((*sv_it).compare("filter")==0) {
	filtLoc = i;
      }
      i++;
    }
    noheader = false;
  }

  if (myFile.eof()) {    
    printf("couldn't find header\n");
    return -2;
  } else if (filtLoc<0) {
    printf("couldn't find energy column\n");
    return -2;
  } else if (seqNLoc<0) {
    printf("couldn't find seqN column\n");
    return -2;
  }
  unsigned int minElems = 1+(seqNLoc>filtLoc)?seqNLoc:filtLoc; // need at least this many elements per line
  
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

    st1.str(spl[seqNLoc]);
    st1 >> seqN;
    if (st1.fail()) {
      std::cout << "couldn't parse first string " << st1.str() << "\n";
      printf("seqN = %d, f = %d\n", seqN,f);
      return -1;
    }
    st1.clear();
    
    st1.str(spl[filtLoc]);
    st1 >> f;
    if (st1.fail()) {
      std::cout << "couldn't parse second string " << st1.str() << "\nQutting\n";
      printf("seqN = %d, f = %d\n", seqN,f);
      return -1;
    }
    st1.clear();
    //printf("goODt:   seqN = %d, f = %Le\n\n",seqN,f);
    
    if (seqN==prevSeqN+1) {
      filt->push_back(f);
      prevSeqN = seqN;
    } else {
      std::cout << "trying to fill sequence at position " << seqN << " with prevSeqN = " << prevSeqN << "  Fail.\n";
      return -1;
    }
    
  } // end read file loop
  myFile.close();

  return 1;
} // end read flter


int Cio_handler::readFasta(char inFile[],char_vector *seq) {

  std::ifstream myFile;
  myFile.open(inFile);  
  if (!myFile.is_open()) {
    printf("failed to open fasta file %s\n", inFile);
    return -1;
  }

  std::string line;
  bool foundSequence = false;  // make sure you don't read two chromosomes/sequences in one .fasta file
  unsigned int i,j;
  short s;
  while (! myFile.eof() ) {
    getline (myFile,line);
    if (!line.length()) {
      continue;
    }
    if (line[0]=='#') {
      continue;
    }
    if (line[0]=='>') {
      if (foundSequence) {
	std::cerr << "WARNING: found more than one sequence in fasta file.  Only reading the first one\n";
	break;
      } else {
	continue;
      }
    }
    foundSequence = true;
    for (i=0;i<line.length();i++) {
      j = (unsigned int)line[i];
      if (j>127) { j=0; }; // j is unsigned so it won't be less than zero
      s = charMap[j];
      if (s=='X') {
	printf("***Error: found a strange letter '%c' in sequence\n",line[i]);
	return -3;
      }
      seq->push_back(s);
    }
  }
  return 1;
}

int Cio_handler::readMultiFasta(char inFile[],char_vector_vector *seqs) {
  
  std::ifstream myFile;
  myFile.open(inFile);  
  if (!myFile.is_open()) {
    printf("failed to open fasta file %s\n", inFile);
    return -1;
  }

  char_vector emptyCV;
  

  char_vector_vector::iterator cvvi;
  std::string line;
  //bool foundSequence = false;  // make sure you don't read two chromosomes/sequences in one .fasta file
  unsigned int i,j;
  short s;
  while (! myFile.eof() ) {
    getline (myFile,line);
    if (!line.length()) {
      continue;
    }
    if (line[0]=='#') {
      continue;
    }
    if (line[0]=='>') {
      std::cout << "reading " << line << "...\n";
      if (seqs->size()) {	
	seqs->push_back(emptyCV);	
	cvvi++;
	continue;
      } else {
	seqs->push_back(emptyCV);	
	cvvi=seqs->begin();	
	continue;
      }
    }

    for (i=0;i<line.length();i++) {
      j = (unsigned int)line[i];
      if (j>127) { j=0; }; // j is unsigned so it won't be less than zero
      s = charMap[j];
      if (s=='X') {
	printf("***Error: found a strange letter '%c' in sequence\n",line[i]);
	return -3;
      }
      cvvi->push_back(s);
    }
  }
  return 1;
  
  
  
}

// read first element of each line in file
int Cio_handler::readList(char inFile[],string_vector *elems) {
  std::ifstream myFile;
  myFile.open(inFile);  
  if (!myFile.is_open()) {
    printf("failed to open list file %s\n", inFile);
    return -1;
  }
  if (myFile.fail()) {
    printf("failbit set in readList\n");
    return -1;
  }

  elems->clear();

  std::string line;
  string_vector spl;
  int i=0;
  while (! myFile.eof() ) {
    i++;
    getline (myFile,line);
    if (!line.length()) {
      continue;
    }
    if (line[0]=='#') {
      continue;
    }
    spl = hand_splitWhiteSpace(line);
    if (spl.size()==0) {
      continue;
    }
    elems->push_back(spl[0]);
  }  
  return 1;
}

// read each element in each line, split by white space
int Cio_handler::readList2(char inFile[],string_vector_vector *elems) {
  elems->clear();
  std::ifstream myFile;
  myFile.open(inFile);  
  if (!myFile.is_open()) {
    printf("failed to open list file %s\n", inFile);
    return -1;
  }
  if (myFile.fail()) {
    printf("failbit set in readList\n");
    return -1;
  }

  std::string line;
  int i=0;
  while (! myFile.eof() ) {
    i++;
    getline (myFile,line);
    //std::cout << i << ") '" << line << std::endl;
    if (!line.length()) {
      continue;
    }
    if (line[0]=='#') {
      continue;
    }
    string_vector spl = hand_splitWhiteSpace(line);
    //dumpVector(spl,"\t");
    if (spl.size()==0) {
      continue;
    }
    elems->push_back(spl);
  } 

  return 1;
}

int Cio_handler::readTable(char fileName[],unsigned short index,double_vector *tab) {
  return readTable(fileName,index,tab,true);
}

// read n'th column in table
int Cio_handler::readTable(char fileName[],unsigned short index,double_vector *tab,bool header) {

  if (index==0) {
    printf("readTable ( %s ): index must be greater than 0!\n",fileName);
    return -2;
  }

  std::ifstream myFile;
  myFile.open(fileName);  
  if (!myFile.is_open()) {
    printf("failed to open list file %s\n", fileName);
    return -1;
  }
  if (myFile.fail()) {
    printf("failbit set in readTable (%s)\n",fileName);
    return -1;
  }

  tab->clear();
  
  std::string line;
  string_vector spl;
  if (header) { // find header
    header = false;
    do {
      getline (myFile,line);
      if (!line.length()) {
	continue;
      }
      if (line[0]=='#') {
	continue;
      }
      header = true;
    } while (!header);
  }
  
  std::stringstream st1(std::stringstream::in | std::stringstream::out);

  double d;
  while (! myFile.eof() ) {
    getline (myFile,line);
    if (!line.length()) {
      continue;
    }
    if (line[0]=='#') {
      continue;
    }
    spl = hand_splitWhiteSpace(line);
    if (spl.size()<index) {
      continue;
    }    
    
    st1.str(spl[index-1]);
    st1 >> d;
    if (st1.fail()) {
      std::cout << "readTable ( %s ) failed to parse string " << st1.str() << " (line: " << line << ")\n";
      return -1;
    }
    st1.clear();

    tab->push_back(d);
  }  

  return 1;
}


/*
// format the energy array as a text file
// return -1 if can't open outFile
int Cio_handler::writeEn(char outFile[], int seqLen,long double en[]) {
   
  std::ofstream oFile(outFile);
  if (!oFile.is_open()) {
    printf("failed to open file %s\n",outFile);
    return -1;
  }

  char line[81];
  sprintf(line,"%8s%8s\n","seqN","E");
  oFile << line;
  
  for (int i=0; i<(seqLen-OBJ_SIZE+1); i++) {
    sprintf(line,"%8d%10.3Lf\n",i+1,en[i]);
    oFile << line;
  }
  oFile.close();

  return 1;
}
*/

char Cio_handler::charMapper(char c) {
  return charMap[int(c)];
}
