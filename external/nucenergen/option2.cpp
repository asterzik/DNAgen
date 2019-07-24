#include "option2.h"
#include <cstring>
#include <stdio.h>
#include <sstream>

// arguments:
//  argc, argv, command line arguments
//  usage: instructions to use the script
//  args: a list of the arguments you wish to search for
string_vector option2_parseOptions(int argc, char* argv[],string_vector args) {

  string_vector ret; // if you didn't find a value, set ""
  if (argc < 2) {
    return ret;
  }

  string_vector::iterator svi;
  std::string val;
  std::string found = "this arg was defined without a value"; // eg you set -mode, not -mode X.
  int i;
  //char a[304];
  for (svi=args.begin();svi<args.end();svi++) {
    val = "";
    for (i=1; i<argc; i++) {
      if (strstr(argv[i],"-")==argv[i] && strstr(argv[i],svi->c_str()) && strlen(argv[i]) == svi->length()+1) {
	if (i+1 < argc) {
	  val.replace(0,0,argv[i+1]);  
	  // supposing you wrote -arg1 -arg2 and you searched for arg1, 
	  //   value is set to -arg2
	} else {
	  val = found;
	}
      }
    }
    ret.push_back(val); 
  }
  return ret;
}

int option2_stringToInt(std::string s) {
  std::stringstream st1(std::stringstream::in | std::stringstream::out);
  st1.str(s);
  int ret;
  st1 >> ret;
  return ret;
}
double option2_stringToDouble(std::string s) {
  std::stringstream st1(std::stringstream::in | std::stringstream::out);
  st1.str(s);
  double ret;
  st1 >> ret;
  return ret;
}
