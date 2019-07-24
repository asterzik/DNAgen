#ifndef H_HANDY_DANDY
#define H_HANDY_DANDY 1

#include <iostream>
#include <string>
#include <vector>

#ifndef STRING_VECTOR
#define STRING_VECTOR 1
typedef std::vector<std::string > string_vector;
#endif

string_vector hand_splitWhiteSpace( std::string s ); // as perl split / +/, $s 
                                                // returns no empty strings
                                                // leading whitespace is ignored (kind of handy, really)
                                                                                          // handy dandy? 
string_vector hand_split(std::string delim,std::string target); // as perl split /$delim+/, $target 

std::string hand_join(std::string delim,string_vector arr); // as perl join($delim,@arr);

#endif
