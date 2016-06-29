// Class for parsing procedures
// load mgf, load ami masses, load ami deltas

#ifndef PARSER_H
#define PARSER_H

#include "ms1ms2.h"
#include "spectra.h"
#include <vector>
#include <string>
#include <map>

//  functions to split string into vector<string> by delimeter
std::vector<std::string> split(const std::string &s, char delim);

class MParser {
public:
  // large and strange function to read and parse mgf file
  // and return vector<Spectrum> object
  // TODO(dima) we need universal proteomic data parser
  std::vector<Spectrum> loadMGF(std::string filename,
    std::map<char, double> masses, double delta,
    bool annot, bool addions);


  // load MS1 delta values from file to vector<ADelta> object
  std::vector<ADelta> loadDeltas(std::string filename);

  // create default MS1 delta values
  // TODO(dima) may be inline file?
  std::vector<ADelta> defaultDeltas();

  // load amino acid masses
  std::map<char, double> loadMasses(std::string filename);

  // default amino acid masses
  // TODO(dima) may be inline file?
  std::map<char, double> defaultMasses();
};

#endif
