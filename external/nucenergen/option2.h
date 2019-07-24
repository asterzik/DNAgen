#ifndef C_OPTION2
#define C_OPTION2 1

#include <iostream>
#include <string>
#include <vector>

#ifndef STRING_VECTOR
#define STRING_VECTOR 1
typedef std::vector<std::string > string_vector;
#endif

string_vector option2_parseOptions(int argc, char* argv[],string_vector args);
int option2_stringToInt(std::string s);
double option2_stringToDouble(std::string s);

#endif
