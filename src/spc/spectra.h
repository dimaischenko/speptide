// Main class is spectrum class for spectrum that consitst of one MS1 Peak
// and some amount of MS2 Peaks. Spectrum has title
// may has sequence.  

#ifndef SPECTRA_H
#define SPECTRA_H

#include <vector>
#include <map>
#include <string>

#include "ms1ms2.h"

// Simple class for amino acid delta
// Class for deltas between amino acids
// contain only POSITIVE delta values (because negative is symmetric)
// can contain multiple deltas for two, three and more amino acids
// for example AB -> CD
// when we have negative delta and want return ami for it
// methods get_ami1(bool plus) and get_ami2(bool plus) swap amino acids
class ADelta {
private:
  std::string ami1; // referent ami group
  std::string ami2; // substitution ami group
  double delta;

public:
  ADelta();
  ADelta(std::string _ami1, std::string _ami2, double _delta);
  // get set block
  double get_delta() const;
  std::string get_ami1(bool plus = true);
  std::string get_ami2(bool plus = true);
};
// operator for comparsion Deltas (order by delta values)
bool operator< (ADelta a, ADelta b);

// function check if delta of masses in vector with deltas
std::vector<ADelta> inRange(double m1, double m2,
  std::vector<ADelta> deltas, double value, bool isPpm);

// calc max delta for vector with ADelta
double calc_max_delta(const std::vector<ADelta>& vdelta);

// function to check consitense element in vector
template< typename T>
bool inV(T x, std::vector<T> v);

// Class for spectrum objects
// contains one MS1Peak object and std::vector of MS2Peak
// in addition contains title and sequence of spectrum
class Spectrum {
private:
  MS1Peak ms1; // MS1Peak
  std::vector<MS2Peak> ms2; // MS2Peaks
  std::string title; // title of spectrum (something like id)
  std::string seq; // sequence (if we have it)
  std::string protein; // something like protein id
  int weight; // number of merged spectra (for spectral library)

public:
  Spectrum();
  Spectrum(MS1Peak _ms1);

  // add ms2 peak to std::vector<MS2Peak> ms2
  void addMS2Peak(MS2Peak ms2p);

  // get all ms2 peaks as Range Collection object 
  RangeC getMS2range(bool query) const;
  RangeC getMS2range(bool query, std::vector<char> charges);

  // set get block
  MS1Peak get_ms1();
  MS1Peak* get_ms1p();
  std::vector<MS2Peak> get_ms2();
  std::vector<MS2Peak>* get_ms2p();
  void set_title(std::string _title);
  void set_seq(std::string _seq);
  void set_protein(std::string _protein);
  void set_weight(int _weight);
  std::string get_title();
  std::string get_seq();
  int get_lseq();
  double get_theor_mass(std::map<char, double> masses);
  void recalc_mass(std::map<char, double> masses);
  // number of ms2 peaks
  int ms2length();
  int get_weight();
  // return number of top (precursor mz / pdiv) ms2 peaks by intensity
  int top_mz_peaks(double pdiv);
  // return number of top (precursor m / pdiv) ms2 peaks by intensity
  int top_m_peaks(double pdiv);
  // number of annotated peaks in spectra
  int get_by_size();
  // return top (precursor mz / pdiv) ms2 peaks by intensity
  std::vector<MS2Peak> get_filter_ms2(double pdiv);
  // retain only top n peaks
  void filter_ms2(std::size_t n);
  // return mean intensity value for all ms2 peaks from spectrum
  double get_mean_intensity();
  // multiply all intensities by constant in spectrum and
  // transform by selected method (iConst + trAlg(mint * intens))
  void inc_intenisies(double mint, char trAlg, double iConst);
  // normalize spectrum intensities in way that mean intensity of
  // spectrum equal mnI constant. if norm is false - then only
  // transform intensities according selected algorithm 
  void normalize(bool norm, char trAlg, double iConst);
  
  // return theoretical Spectrum for sequence
  // create all b/y/B/Y (B = b - NH3, Y = y - NH3, P = b - H2O, I = y - H2O) ions with charge = 1
  // and intensity = 0
  Spectrum getTheorSpectrum(std::map<char, double> masses, bool addions);
  
  // find aa in sequence
  std::vector<int> findAmiSeq(char ami) ;
  
  // find availible delta position in seq for current ADelta
  // some explanation. we have delta, it may be multiple (AB -> CD)
  // and find all availible variants in sequence that may be realisied
  // return matrix (std::vector< std::vector<int> >) with positions in which deltas can appear
  // each row of matrix is one of availible variants
  // each col of matrix corresponds to single ami substitutions
  // for delta AB->CD. first column corresponds to A->C and second to B->D
  std::vector< std::vector<int> > getDeltaPos(ADelta adelta, double ms1delta);

  // Block of shifting methods.
  // We test different approaches and now it is several methods.
  // May be in future i will remove them and left only one.

  // common for all algorithms method to calculate deltas in each position
  // in sequence calculate delta value for each position by array 
  // with <pos, delta> std::pairs
  std::vector<double> deltaPos(std::vector< std::pair<int, double> > dpos);
  // levae only b/y ions and shift them (if ion in needed position)
  // erase other ions
  void shiftBY(std::vector< std::pair<int, double> > dpos);
  // remove all not annotated peaks
  void removeNA();
  // shift only b/y ions
  void shift1(std::vector< std::pair<int, double> > dpos);
  // shift only b/y ions at needed deltas 
  // and other (not annotated) to summary delta
  void shift2(std::vector< std::pair<int, double> > dpos);
  // shift only b/y ions at needed deltas 
  // other (not annotated) duplicate with summary delta
  void shift3(std::vector< std::pair<int, double> > dpos);
  // print spectrum
  void print();
};

// return vector<Range> object with ms1 mz windows from vector<Spectrum> object
// TODO(dima) move to RangeC type
std::vector<Range> getMS1range(std::vector<Spectrum>* spv, bool query);

// Subset from vector of Spectrum by sequences or titles
// checking that value in vector of values
// subset from vector<Spectrum> by titles vector
std::vector<Spectrum> subsetTitleSpectra(std::vector<Spectrum>*spv, 
  std::vector<std::string> titles, bool negate = false);
// subset from vector<Spectrum> by sequences vector
std::vector<Spectrum> subsetSeqSpectra(std::vector<Spectrum>*spv, 
  std::vector<std::string> seqs, bool negate = false);

// Spectral library class
// Spectra library consists of spectra for unique sequences. While creation
// from vector of spectra, spectra for identical sequences merged to one spectrum.
class SpectraLib {
private:
  std::map<std::string, Spectrum> slb; // map of sequence -> spectra

public: 
  SpectraLib();
  // add one spectrum to library
  void add_spectrum(Spectrum sp, double delta);
  // add several spectra to library from vector<Spectrum>
  void add_spectra(std::vector<Spectrum>* spv, double delta);
  // size of spectra library
  int get_size();
  // normalize each ms2 peak intensity to number of spectra for seq
  // TODO(dima) remember that after normalization weights are the same
  void normalize_lib();
  // get vector of spectra
  std::vector<Spectrum> get_spectra_vector();
};

#endif
