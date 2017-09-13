// Class for parsing procedures
// load mgf, load ami masses, load ami deltas

#include "mparser.h"
#include "ms1ms2.h"
#include "spectra.h"
#include "scoring.h"

#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>

#define mnI 100 // mean intensity
#define minPeaks 10 // min number of peaks in spectra
#define minBYPeaks 5 // min number of by peaks in lib

// functions to split string into vector<string> by delimeter
std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

// large and strange function to read and parse mgf file
// and return vector<Spectrum> object
// TODO(dima) we need universal proteomic data parser
std::vector<Spectrum> MParser::loadMGF(
  const std::string& filename,
  const std::map<char, double>& masses,
  double delta,
  bool annot,
  bool addions) {
  
  std::vector<Spectrum> spv = std::vector<Spectrum>(0);
  int id = 0;
  std::string line;
  Spectrum sp = Spectrum();
  // nocharge bool flag
  bool isNoCharge = true;
  // read mgf
  std::ifstream infile(filename.c_str());
  if (!infile) throw std::invalid_argument("There is no such file with mgf.");
  // go through each line
  while(getline(infile, line)) {
    // TODO(dima) end of line in windows format
    if (!line.empty() && line[line.size() - 1] == '\r')
      line.erase(line.size() - 1);
    if(!line.compare(0, 1, "#") || line[0] == 0) continue;
    if(!line.compare(0, 10, "BEGIN IONS")) {
      id++;
      sp = Spectrum();      
      isNoCharge = true;
    } else {
      std::replace(line.begin(), line.end(), '\t', ' ');
      std::vector<std::string> el = split(line, '=');
      if (el[0].compare(0,5,"TITLE") == 0) {
        sp.set_title(el[1]);
      } else if (!el[0].compare(0, 7, "PEPMASS")) {
        sp.get_ms1p()->set_mz(std::atof(el[1].c_str()));
      } else if (!el[0].compare(0, 6, "CHARGE")) {
        sp.get_ms1p()->set_charge(el[1][0] - '0');
        // charge value exists
        isNoCharge = false;
      } else if (!el[0].compare(0, 3, "SEQ")) { sp.set_seq(el[1]);
      } else if (!el[0].compare(0, 11, "RTINSECONDS")) {
        sp.get_ms1p()->set_rt(std::atof(el[1].c_str()));
      } else if(el.size() == 1) {
        // create ms2 peaks
        do {
          if (!line.empty() && line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);
          if (!line.compare(0, 1, "#") || line[0] == 0) continue;
          std::replace(line.begin(), line.end(), '\t', ' ');
          std::vector<std::string> ms2el = split(line, ' ');
          MS2Peak ms2 = MS2Peak(std::atof(ms2el[0].c_str()));
          ms2.set_intensity(std::atof(ms2el[1].c_str()));
          if (ms2el.size() == 3 && ms2el[2][0] != 0)
            ms2.set_charge(ms2el[2][0] - '0');
          sp.addMS2Peak(ms2);
        } while(getline(infile, line) && line.compare(0,8,"END IONS"));
        
        // annotate if we need it
        if (annot) {
          // get reference spectrum by sequence
          Spectrum refp = sp.getTheorSpectrum(masses, addions);
          // annotate spectrum
          RangeC r1 = refp.getMS2range(true);
          RangeC r2 = sp.getMS2range(false);
          // overlap peaks mz
          olapdata nc = overlap2range(r1, r2, delta);
          for (std::size_t i = 0; i < nc.op1.size(); i++) {
            sp.get_ms2p()->at(nc.op2[i]).set_series(refp.get_ms2()[nc.op1[i]].get_series());
            sp.get_ms2p()->at(nc.op2[i]).set_number(refp.get_ms2()[nc.op1[i]].get_number());
          }
          // replace MS1 mass with theoretical
          sp.recalc_mass(masses);
        }

        // add spectrum to result vector if it pass some filters
        if((sp.get_ms2p()->size() >= minPeaks))  {
          if(!annot || (annot && sp.get_by_size() >= minBYPeaks)) {
            // add charge to title
            std::string sp_tit = sp.get_title();
            sp.set_title(sp_tit + "." + intToString(sp.get_ms1p()->get_charge()));
            spv.push_back(sp);
            // if no charge selected add spectra with several default charges
            // TODO(dima) accurate default charges sets
            if (isNoCharge) {
              sp.set_title(sp_tit + ".3");
              sp.get_ms1p()->set_charge(3);
              spv.push_back(sp);
            }
          }
        }
      }
    }
  }
  
  infile.close();
  return spv;
}

// load MS1 delta values from file to vector<ADelta> object
std::vector<ADelta> MParser::loadDeltas(const std::string& filename) {
  std::vector<ADelta> ms1deltas = std::vector<ADelta>(0);
  // open file with deltas
  std::ifstream infile(filename.c_str());
  if (!infile) throw std::invalid_argument("There is no such file with deltas.");
  std::string line;
  while(std::getline(infile, line)) {
    std::vector<std::string> el = split(line, ' ');
    ADelta d = ADelta(el[0], el[1], std::atof(el[2].c_str()));
    ms1deltas.push_back(d);
  }

  infile.close();
  return ms1deltas;
}

// create default MS1 delta values
std::vector<ADelta> MParser::defaultDeltas() {
  std::vector<ADelta> ms1deltas = std::vector<ADelta>(0);
  
  // manual deltas
  ms1deltas.push_back(ADelta("K", "M", 2.945522));
  ms1deltas.push_back(ADelta("P", "T", 3.994915));
  ms1deltas.push_back(ADelta("Q", "H", 9.000334));
  ms1deltas.push_back(ADelta("S", "P", 10.020736));
  ms1deltas.push_back(ADelta("T", "I", 12.036385));
  ms1deltas.push_back(ADelta("T", "N", 12.995248));
  ms1deltas.push_back(ADelta("V", "I", 14.01565));
  ms1deltas.push_back(ADelta("V", "L", 14.01565));
  ms1deltas.push_back(ADelta("D", "E", 14.01565));
  ms1deltas.push_back(ADelta("G", "A", 14.01565));
  ms1deltas.push_back(ADelta("S", "T", 14.015651));
  ms1deltas.push_back(ADelta("N", "K", 14.052036));
  ms1deltas.push_back(ADelta("L", "Q", 14.974514));
  ms1deltas.push_back(ADelta("I", "K", 15.010899));
  ms1deltas.push_back(ADelta("V", "D", 15.958529));
  ms1deltas.push_back(ADelta("S", "C", 72.977157)); //+57
  ms1deltas.push_back(ADelta("F", "Y", 15.994906));
  ms1deltas.push_back(ADelta("A", "S", 15.994914));
  ms1deltas.push_back(ADelta("P", "L", 16.0313));
  ms1deltas.push_back(ADelta("L", "M", 17.956421));
  ms1deltas.push_back(ADelta("I", "M", 17.956421));
  ms1deltas.push_back(ADelta("H", "R", 19.042199));
  ms1deltas.push_back(ADelta("D", "H", 22.031969));
  ms1deltas.push_back(ADelta("N", "H", 23.015985));
  ms1deltas.push_back(ADelta("L", "H", 23.974848));
  ms1deltas.push_back(ADelta("M", "R", 25.060626));
  ms1deltas.push_back(ADelta("H", "Y", 26.004408));
  ms1deltas.push_back(ADelta("A", "P", 26.01565));
  ms1deltas.push_back(ADelta("S", "L", 26.052036));
  ms1deltas.push_back(ADelta("S", "I", 26.052036));
  ms1deltas.push_back(ADelta("S", "N", 27.010899));
  ms1deltas.push_back(ADelta("T", "K", 27.047284));
  ms1deltas.push_back(ADelta("K", "R", 28.006148));
  ms1deltas.push_back(ADelta("A", "V", 28.0313));
  ms1deltas.push_back(ADelta("Q", "R", 28.042533));
  ms1deltas.push_back(ADelta("V", "E", 29.974179));
  ms1deltas.push_back(ADelta("R", "W", 29.978202));
  ms1deltas.push_back(ADelta("T", "M", 29.992806));
  ms1deltas.push_back(ADelta("G", "S", 30.010564));
  ms1deltas.push_back(ADelta("A", "T", 30.010565));
  ms1deltas.push_back(ADelta("P", "Q", 31.005814));
  ms1deltas.push_back(ADelta("V", "M", 31.972071));
  ms1deltas.push_back(ADelta("L", "F", 33.98435));
  ms1deltas.push_back(ADelta("I", "F", 33.98435));
  ms1deltas.push_back(ADelta("P", "H", 40.006148));
  ms1deltas.push_back(ADelta("G", "V", 42.04695));
  ms1deltas.push_back(ADelta("L", "R", 43.017047));
  ms1deltas.push_back(ADelta("I", "R", 43.017047));
  ms1deltas.push_back(ADelta("A", "D", 43.989829));
  ms1deltas.push_back(ADelta("C", "F", 44.059229));
  ms1deltas.push_back(ADelta("G", "C", 102.987721)); //+57
  ms1deltas.push_back(ADelta("V", "F", 48));
  ms1deltas.push_back(ADelta("D", "Y", 48.036377));
  ms1deltas.push_back(ADelta("N", "Y", 49.020393));
  //ms1deltas.push_back(ADelta("C", "R", 53.091926));
  ms1deltas.push_back(ADelta("T", "R", 55.053432));
  ms1deltas.push_back(ADelta("A", "E", 58.005479));
  ms1deltas.push_back(ADelta("G", "D", 58.005479));
  ms1deltas.push_back(ADelta("P", "R", 59.048347));
  ms1deltas.push_back(ADelta("S", "F", 60.036386));
  //ms1deltas.push_back(ADelta("C", "Y", 60.054135));
  ms1deltas.push_back(ADelta("S", "R", 69.069083));
  ms1deltas.push_back(ADelta("G", "E", 72.021129));
  ms1deltas.push_back(ADelta("L", "W", 72.995249));
  ms1deltas.push_back(ADelta("S", "Y", 76.031292));
  ms1deltas.push_back(ADelta("C", "W", 26.070128)); //-57
  ms1deltas.push_back(ADelta("S", "W", 99.047285));
  ms1deltas.push_back(ADelta("G", "R", 99.079647));
  ms1deltas.push_back(ADelta("G", "W", 129.057849));
  return ms1deltas;
}


// load amino acid masses
std::map<char, double> MParser::loadMasses(const std::string& filename) {
  std::map<char, double> masses = std::map<char, double>();
  // open file with masses
  std::ifstream infile(filename.c_str());
  if (!infile) throw std::invalid_argument("There is no such file with amino acids masses.");
  std::string line;
  while(std::getline(infile, line)) {
    std::vector<std::string> el = split(line, ' ');
    masses.insert(std::pair<char, double>(el[0][0],
          (double)std::atof(el[1].c_str())));
  }
  infile.close();

  return masses;
}

// default amino acid masses
std::map<char, double> MParser::defaultMasses() {
  std::map<char, double> masses = std::map<char, double>();
  // default masses
  masses['A'] = 71.037114;
  masses['R'] = 156.101111;
  masses['N'] = 114.042927;
  masses['D'] = 115.026943;
  //masses['C'] = 103.009185;
  masses['C'] = 160; // +carboxymehyl mod
  masses['E'] = 129.042593;
  masses['Q'] = 128.058578;
  masses['G'] = 57.021464;
  masses['H'] = 137.058912;
  masses['I'] = 113.084064;
  masses['L'] = 113.084064;
  masses['K'] = 128.094963;
  masses['M'] = 131.040485;
  masses['F'] = 147.068414;
  masses['P'] = 97.052764;
  masses['S'] = 87.032028;
  masses['T'] = 101.047679;
  masses['U'] = 150.95363;
  masses['W'] = 186.079313;
  masses['Y'] = 163.06332;
  masses['V'] = 99.068414;

  return masses;
}
