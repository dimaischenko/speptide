// Class for calculations angles between spectra and libs

#ifndef SCORING_H
#define SCORING_H

#include<vector>
#include<string>
#include<map>

#include "ms1ms2.h"
#include "spectra.h"
#include "result.h" 

// structure for compspecBY function
struct compdata {
  int qa, sa, qo, so, o;
  double angle;
};

// int to string
std::string intToString(int i);

class Scoring {
public:
  Scoring() {}

// calculate cos of angle between two vectors
// if one of vector equal zero vector (0, 0, 0, 0) return 0
double sangle(const std::vector<double>& v1, const std::vector<double>& v2);

// calculate distance between two spectra ms2 data
double comp2speci (const RangeC& r1, const RangeC& r2, double delta, double pdiv);

// Calculation of distance between two spectra
// first spectrum(sp1) is query spectrum. second spectrum (sp2) is referent spectrum
// function apply some modification to referent spectrum (e.g. shift b/y ions according to delta)
// select top (mz / pdiv) MS2 peaks from query spectrum by intensities
// intersect them with referent MS2 peaks and calculate spectral angle
// return structure with angle, amounts of intersected peaks
compdata comp2spec (Spectrum sp1, Spectrum sp2, double delta, double pdiv,
  std::vector< std::pair<int, double> > dpos);

// function compare two spectra object my therir ms1 and ms2 data
ResultIdent comp2spvAngle (std::vector<Spectrum> sp1v, std::vector<Spectrum> sp2v,
  double value, bool isPpm, double delta, double pdiv, double ath, bool norm,
  double iConst, char trAlg); 


// function compare two vector<Specta> object my their ms1 and ms2 data
// second vector<Spectrum> considered as a reference (with sequences and annotated peaks)
// we find MS1 deltas between query and referent spectra and check them
// whether they can be candidates for amino acid substitution. MS1 delta beints to {ADeltas}
// next step. for each candidate pair we find all possible amino acid substitution
// positions (by reference sequence) and calculate spectral angle for each realised
// substitution (with shifting referent peaks)
ResultSap comp2spvAngleAap (std::vector<Spectrum> sp1v,
  std::vector<Spectrum> sp2v, double value, bool isPpm, double delta, double pdiv,
  std::vector<ADelta> deltas, std::map<char, double> masses, double ath,
  bool norm, double iConst, char trAlg, double refdiv);
};

#endif
