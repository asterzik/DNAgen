#include "handyDandy.h"

string_vector hand_splitWhiteSpace( std::string s ) {

  std::string elem = "";
  //int count = 0;
  int len = s.length();
  string_vector ret;

  if (!len) {
    return ret;
  }

  for (int i=0; i<len; i++) {
    if (s[i]!=' ') {
      elem+=s[i];
    } else if (elem.length()) {
      ret.push_back(elem);
      elem = "";
    } 
  } 

  if (elem.length()) {
    ret.push_back(elem);
    elem = "";
  } 

  return ret;
}

// as perl split /$delim+/, $target 
string_vector hand_split(std::string delim,std::string target) {

  std::string elem = "";
  //int count = 0;
  int len = target.length();
  string_vector ret;

  if (!len) {
    return ret;
  }

  std::string sub;
  for (int i=0; i<len; i++) {
    sub = target.substr(i,1);
    if (sub.compare(delim)) { // unless sub eq delim
      elem+=sub;
    } else if (elem.length()) {
      ret.push_back(elem);
      elem = "";
    } 
  } 
  if (elem.length()) {
    ret.push_back(elem);
    elem = "";
  } 

  return ret;
}

// as perl join($delim,@arr);
std::string hand_join(std::string delim,string_vector arr) {
  std::string ret = "";
  string_vector::iterator svi;
  for (svi = arr.begin(); svi < arr.end(); svi++) {
    ret+=*svi;
    if (svi+1<arr.end()) {
      ret+=delim;
    }
  }
  return ret;
}
