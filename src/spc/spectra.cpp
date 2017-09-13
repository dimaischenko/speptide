// Class for spectrum objects
// contains one MS1Peak object and std::vector of MS2Peak
// in addition contains title and sequence of spectrum

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "ms1ms2.h"
#include "spectra.h"

#define mH 1.007825
#define mO 15.99491
#define mN 14.0067

#define mnI 100 // mean intensity

// Simple class for amino acid delta
// Class for deltas between amino acids
// contain only POSITIVE delta values (because negative is symmetric)
// can contain multiple deltas for two, three and more amino acids
// for example AB -> CD
// when we have negative delta and want return ami for it
// methods get_ami1(bool plus) and get_ami2(bool plus) swap amino acids
ADelta::ADelta() {
  ami1 = "";
  ami2 = "";
  delta = 0;
}

ADelta::ADelta(std::string _ami1, std::string _ami2, double _delta) {
  // swap amies when delta is negative
  if (_delta >= 0) {
    ami1 = _ami1;
    ami2 = _ami2;
    delta = _delta;
  } else {
    ami1 = _ami2;
    ami2 = _ami1;
    delta = _delta;
  }
}
// get set block
double ADelta::get_delta() const {
  return delta;
}
std::string ADelta::get_ami1(bool plus) {
  if (plus) return ami1;
  return ami2;
} 
std::string ADelta::get_ami2(bool plus) {
  if (plus) return ami2;
  return ami1;
}
// operator for comparsion Deltas (order by delta values)
bool operator< (ADelta a, ADelta b) {
  return (a.get_delta() < b.get_delta());
}

// function check if delta of masses in vector with deltas
std::vector<ADelta> inRange(double m1, double m2,
  const std::vector<ADelta>& deltas, double value, bool isPpm) {
  double dt = std::abs(m2 - m1);
  std::vector<ADelta> rdeltas = std::vector<ADelta>(0); 

  double mdelta = (isPpm) ? m1 * value * 1e-6 : value; 

  std::size_t i = 0;
  while(dt >= (deltas[i].get_delta() - mdelta) && i < deltas.size()) {
    if (dt <= deltas[i].get_delta() + mdelta) {
      rdeltas.push_back(deltas[i]);
    }
    
    i++;
  }

  return rdeltas;
}

// calc max delta for vector with ADelta
double calc_max_delta(const std::vector<ADelta>& vdelta) {
  double mdelta = 0;
  for (std::size_t i = 0; i < vdelta.size(); i++) {
    double cdelta = vdelta[i].get_delta();
    if (cdelta > mdelta) mdelta = cdelta;
  }
  return mdelta;
}

// function to check consitense element in vector
template< typename T >
bool inV(T x, const std::vector<T>& v) {
  for (int i = 0; i < v.size(); i++) {
    if (x == v[i]) return true;
  }
  return false;
}

// Class for spectrum objects
// contains one MS1Peak object and std::vector of MS2Peak
// in addition contains title and sequence of spectrum
Spectrum::Spectrum() {
  ms1 = MS1Peak();
  ms2 = std::vector<MS2Peak>(0);
  seq = "";
  title = "";
  weight = 1;
  protein = "";
}

Spectrum::Spectrum(const MS1Peak& _ms1) {
  ms1 = MS1Peak(_ms1);
  ms2 = std::vector<MS2Peak>(0);
  seq = "";
  title = "";
  weight = 1;
  protein = "";
}

// add ms2 peak to std::vector<MS2Peak> ms2
void Spectrum::addMS2Peak(const MS2Peak& ms2p) {
  ms2.push_back(ms2p);
}

// get all ms2 peaks as range collection object 
RangeC Spectrum::getMS2range(bool query) const {
  RangeC rc = RangeC();
  for (std::size_t i = 0; i < ms2.size(); i++) {
    Range r = Range((MSPeak)ms2[i], i, query, i);
    rc.add_range(r); // add range to collection
    rc.add_midinter(i); // add counter for each parent peak
    rc.add_mid2first(i); // add parent to first child peak link
  }
  return rc;
}

// get all ms2 peaks with diff charges as range collection object 
RangeC Spectrum::getMS2range(bool query, const std::vector<char>& charges) {
  RangeC rc = RangeC();
  size_t chsize = charges.size();
  for (std::size_t i = 0; i < ms2.size(); i++) {
    for (std::size_t ch = 0; ch < charges.size(); ch++) {
      MSPeak newp = MSPeak(ms2[i]);
      newp.set_charge(charges[ch]);
      Range r = Range(newp, i * chsize + ch, query, i);
      rc.add_range(r);
    }
    rc.add_midinter(i);
    rc.add_mid2first(i * chsize);
  }
  return rc;
}

// set get block
MS1Peak Spectrum::get_ms1() {
  return ms1;
}

MS1Peak* Spectrum::get_ms1p() {
  return &ms1;
}

std::vector<MS2Peak> Spectrum::get_ms2() const {
  return ms2;
}

std::vector<MS2Peak>* Spectrum::get_ms2p() {
  return &ms2;
}

void Spectrum::set_title(std::string _title) {
  title = _title;
}

void Spectrum::set_seq(std::string _seq) {
  seq = _seq;
}

void Spectrum::set_protein(std::string _protein) {
  protein = _protein;
}

void Spectrum::set_weight(int _weight) {
  weight = _weight;
}

std::string Spectrum::get_title() const {
  return title;
}

std::string Spectrum::get_seq() const {
  return seq;
}

int Spectrum::get_lseq() const {
  return seq.size();
}

double Spectrum::get_theor_mass(const std::map<char, double>& masses) const {
  if (seq.size() == 0) return 0;
  double M = 2 * mH + mO;
  for (std::size_t i = 0; i < seq.size(); i++) {
    M += masses.at(seq[i]);
  }
  return M;
}

void Spectrum::recalc_mass(const std::map<char, double>& masses) {
  ms1.set_m(get_theor_mass(masses));
}

// number of ms2 peaks
int Spectrum::ms2length() const {
  return ms2.size();
}

int Spectrum::get_weight() const {
  return weight;
}
// return number of top (precursor mz / pdiv) ms2 peaks by intensity
int Spectrum::top_mz_peaks(double pdiv) {
  // calcl number and compare with number of ms2 peaks
  int npeaks = std::min(floor(ms1.get_mz() / pdiv), floor(ms2.size()));

  return npeaks;
}
// return number of top (precursor m / pdiv) ms2 peaks by intensity
int Spectrum::top_m_peaks(double pdiv) {
  // calcl number and compare with number of ms2 peaks
  int npeaks = std::min(floor(ms1.get_m() / pdiv), floor(ms2.size()));

  return npeaks;
}

// number of annotated peaks in spectra
int Spectrum::get_by_size() const {
  int bysz = 0;
  for (std::size_t i = 0; i < ms2.size(); i++) {
    if(ms2[i].get_series() != 'n')
      bysz++;
  }
  return bysz;
}

// return top (precursor mz / pdiv) ms2 peaks by intensity
std::vector<MS2Peak> Spectrum::get_filter_ms2(double pdiv) {
  // sort peaks by intensity
  std::sort(ms2.begin(), ms2.end());
  // create subset
  std::vector<MS2Peak> ms2f(ms2.begin(), ms2.begin() + top_mz_peaks(pdiv));

  return ms2f;
}

// retain only top n peaks
void Spectrum::filter_ms2(std::size_t n) {
  // sort peaks by intensity
  std::sort(ms2.begin(), ms2.end());
  if (n < ms2.size())
    ms2.resize(n);
}

// return mean intensity value for all ms2 peaks from spectrum
double Spectrum::get_mean_intensity() {
  double mi = 0;
  for (std::size_t i = 0; i < ms2.size(); i++) {
    mi += ms2[i].get_intensity();
  }
  mi /= ms2.size();
  return mi;
}

// multiply all intensities by constant in spectrum and
// transform by selected method (iConst + trAlg(mint * intens)
void Spectrum::inc_intenisies(double mint, char trAlg, double iConst) {
  for (std::size_t i = 0; i < ms2.size(); i++) {
    ms2[i].set_intensity(ms2[i].get_intensity() * mint, trAlg, iConst);
  }
}

// normalize spectrum intensities in way that mean intensity of
// spectrum equal mnI constant. if norm is false - then only
// transform intensities according selected algorighm 
void Spectrum::normalize(bool norm, char trAlg, double iConst) {
  double mint = norm ? mint = mnI / this->get_mean_intensity() : 1;
  this->inc_intenisies(mint, trAlg, iConst);
}

// return theoretical Spectrum for sequence
// create all b/y/B/Y (B = b - NH3, Y = y - NH3, P = b - H2O, I = y - H2O) ions with charge = 1
// and intensity = 0
Spectrum Spectrum::getTheorSpectrum(const std::map<char, double>& masses, bool addions) {
  Spectrum theorsp = Spectrum();
  int pep_size = seq.size();
  double M = 2 * mH + mO;
  // first last B/Y P/I availible sites
  int BYf = pep_size; int BYl = -1;
  int PIf = pep_size; int PIl = -1;

  // ami for B/Y ions (-NH3)
  const char byarr[] = {'R', 'K', 'Q', 'N'};
  std::vector<char> BYami(byarr, byarr + sizeof(byarr) / sizeof(byarr[0])); 
  // ami for P/I ions (-H20)
  const char pyarr[] = {'S', 'T', 'E','D'};
  std::vector<char> PIami(pyarr, pyarr + sizeof(pyarr) / sizeof(pyarr[0]));
  // calc neutural peptide mass and B/Y P/I first last positions
  for (int i = 0; i < pep_size; i++) {
    M += masses.at(seq[i]);
    if (inV(seq[i], BYami)) {
      if (BYf == pep_size) BYf = i;
      BYl = i;
    }
    if (inV(seq[i], PIami)) {
      if (PIf == pep_size) PIf = i;
      PIl = i;
    }
  }

  double m1 = 0;
  for (int i = 0; i < (pep_size - 1); i++) {
    m1 += masses.at(seq[i]);
    // add b y ioins
    MS2Peak bion = MS2Peak(m1 + mH, 0, 1, true, 'b', (char)(i + 1));
    MS2Peak yion = MS2Peak(M - m1 + mH, 0, 1, true, 'y', (char)(pep_size - (i + 1)));
    theorsp.addMS2Peak(bion);
    theorsp.addMS2Peak(yion);

    if (addions) {
      // add b - NH3 (B) and y - NH3 (Y) ions
      MS2Peak Bion = MS2Peak(m1 + mH - 3 * mH - mN, 0, 1, true, 'B', (char)(i+1));
      MS2Peak Yion = MS2Peak(M - m1 + mH - 3 * mH - mN, 0, 1, true, 'Y', (char)(pep_size - (i+1)));
      if (i >= BYf) theorsp.addMS2Peak(Bion);
      if (i <= BYl) theorsp.addMS2Peak(Yion);
      // add b - H20 (P) and y - H20 (I) ions
      MS2Peak Pion = MS2Peak(m1 + mH - 2 * mH - mO, 0, 1, true, 'P', (char)(i+1));
      MS2Peak Iion = MS2Peak(M - m1 + mH - 2 * mH - mO, 0, 1, true, 'I', (char)(pep_size - (i+1)));
      if (i >= PIf) theorsp.addMS2Peak(Pion);
      if (i <= PIl) theorsp.addMS2Peak(Iion);
    }

  }

  return theorsp;
}
// find positions of aminoaccid in sequence
std::vector<int> Spectrum::findAmiSeq(char ami) {
  std::vector<int> pos = std::vector<int>(0);
  for (std::size_t i = 0; i < seq.size(); i++) {
    if (seq[i] == ami) {
      pos.push_back(i+1);
    }
  }
  return pos;
}

// find availible delta position in seq for current ADelta
// some explanation. we have delta, it may be multiple (AB -> CD)
// and find all availible variants in sequence that may be realisied
// return matrix (std::vector< std::vector<int> >) with positions in which deltas can appear
// each row of matrix is one of availible variants
// each col of matrix corresponds to single ami substitutions
// for delta AB->CD. first column corresponds to A->C and second to B->D
std::vector< std::vector<int> > Spectrum::getDeltaPos(ADelta adelta, double ms1delta) {
  std::vector< std::vector<int> > pos = std::vector< std::vector<int> >(0);
  // get amies
  std::string newc = adelta.get_ami1(ms1delta > 0);
  std::string oldc = adelta.get_ami2(ms1delta > 0);
  
  for (std::size_t ai = 0; ai < oldc.size(); ai++) {
    std::vector<int> apos = findAmiSeq(oldc[ai]);
    if (apos.size() == 0) { pos = std::vector< std::vector<int> >(0); break; }
    int poss = pos.size();  
    for(std::size_t aposi = 0; aposi < apos.size(); aposi++) {
      if (ai == 0) {
        std::vector<int> fami = std::vector<int>(0);
        fami.push_back(apos[aposi]);
        pos.push_back(fami);
      } else {
        for (int posi = 0; posi < poss; posi++) {
          if (aposi == 0) {
            pos[posi].push_back(apos[aposi]);
          } else {
          std::vector<int> np(pos[posi]);
          np[np.size() - 1] = apos[aposi];
          pos.push_back(np);
          }
        }
      }
    }
  }

  return pos;
}

// Block of shifting methods.
// We test different approaches and now it is several methods.
// May be in future i will remove them and left only one.

// common for all algorithms method to calculate deltas in each position
// in sequence calculate delta value for each position by array 
// with <pos, delta> pairs
std::vector<double> Spectrum::deltaPos(std::vector< std::pair<int, double> > dpos) {
  // sort deltas by pos
  std::sort(dpos.begin(), dpos.end());
  // get sequqnce length
  int seql = this->get_lseq();
  // calc delta for each position
  std::vector<double> bdeltas = std::vector<double>(0);
  bdeltas.push_back(0);
  double cdelta = 0; std::size_t cp = 0;
  for (int i = 1; i <= seql; i++) {
    if(i > dpos[cp].first && cp < (dpos.size() - 1)) cp++;
    if(i == dpos[cp].first) cdelta += dpos[cp].second;
    bdeltas.push_back(cdelta);
  }

  return bdeltas;
}

// levae only b/y ions and shift them (if ion in needed position)
// erase other ions
void Spectrum::shiftBY(const std::vector< std::pair<int, double> >& dpos) {
  // get deltas for each position
  std::vector<double> bdeltas = this->deltaPos(dpos);
  // summary delta
  double sdelta = bdeltas[this->get_lseq()];
  // shift by. othere delete
  std::vector<MS2Peak>::iterator sp2i = this->ms2.begin();
  // vector with B-type ions
  const char btype[] = {'b', 'B', 'P'};
  std::vector<char> btypev(btype, btype + sizeof(btype) / sizeof(btype[0]));
  const char ytype[] = {'y', 'Y', 'I'};
  std::vector<char> ytypev(ytype, ytype + sizeof(ytype) / sizeof(ytype[0]));
  
  // TODO(dima) remember that we erase not annotated early.
  for (sp2i = this->ms2.begin(); sp2i != this->ms2.end(); ++sp2i) {
    char cur_num = (*sp2i).get_number();
    char cur_series = (*sp2i).get_series();
    if (inV(cur_series,btypev)) {
      (*sp2i).add_mz(-bdeltas[cur_num] / (*sp2i).get_charge());
    } else {
      (*sp2i).add_mz(-(sdelta - bdeltas[this->get_lseq() - cur_num]) / (*sp2i).get_charge());
    }
  }
}

// remove all not annotated peaks from spectrum
void Spectrum::removeNA() {
  std::vector<MS2Peak>::iterator spi = this->ms2.begin();
  for (spi = this->ms2.begin(); spi != this->ms2.end();) {
    if ((*spi).get_series() == 'n') {
      this->ms2.erase(spi);
    } else {
      ++spi;
    }
  }
}

// shift only b/y ions
void Spectrum::shift1(const std::vector< std::pair<int, double> >& dpos) {
  // get deltas for each position
  std::vector<double> bdeltas = this->deltaPos(dpos);
  // summary delta
  double sdelta = bdeltas[this->get_lseq()];
  // shift by. othere delete
  std::vector<MS2Peak>::iterator sp2i = this->ms2.begin();
  // go through each ms2 peak and add new wih +DELTA mz value
  for (sp2i =this->ms2.begin(); sp2i != this->ms2.end();) {
    char cur_series = (*sp2i).get_series();
    if (cur_series == 'b' || cur_series == 'y') {
      char cur_num = (*sp2i).get_number();
      if (cur_series == 'b') {
        (*sp2i).add_mz(-bdeltas[cur_num] / (*sp2i).get_charge());
      } else {
        (*sp2i).add_mz(-(sdelta - bdeltas[this->get_lseq() - cur_num]) / (*sp2i).get_charge());
      }
    }
    ++sp2i;
  }
}

// shift only b/y ions at needed deltas 
// and other (not annotated) to summary delta
void Spectrum::shift2(const std::vector< std::pair<int, double> >& dpos) {
  // get deltas for each position
  std::vector<double> bdeltas = this->deltaPos(dpos);
  // summary delta
  double sdelta = bdeltas[this->get_lseq()];
  // shift by. othere delete
  std::vector<MS2Peak>::iterator sp2i = this->ms2.begin();
  // go through each ms2 peak and add new wih +DELTA mz value
  for (sp2i = this->ms2.begin(); sp2i != this->ms2.end();) {
    char cur_series = (*sp2i).get_series();
    if (cur_series == 'b' || cur_series == 'y') {
      char cur_num = (*sp2i).get_number();
      if (cur_series == 'b') {
        (*sp2i).add_mz(-bdeltas[cur_num] / (*sp2i).get_charge());
      } else {
        (*sp2i).add_mz(-(sdelta - bdeltas[this->get_lseq() - cur_num]) / (*sp2i).get_charge());
      }
    } else {
      (*sp2i).add_mz(-sdelta / (*sp2i).get_charge());
    }
    ++sp2i;
  }
}

// shift only b/y ions at needed deltas 
// other (not annotated) duplicate with summary delta
void Spectrum::shift3(const std::vector< std::pair<int, double> >& dpos) {
  // get deltas for each position
  std::vector<double> bdeltas = this->deltaPos(dpos);
  // summary delta
  double sdelta = bdeltas[this->get_lseq()];
  // shift by. othere delete
  std::vector<MS2Peak>::iterator sp2i = this->ms2.begin();
  // std::vector for new (shifted MS2 peaks)
  std::vector<MS2Peak> dpeaks = std::vector<MS2Peak>(0);
  // go through each ms2 peak and add new wih +DELTA mz value
  for (sp2i = this->ms2.begin(); sp2i != this->ms2.end();) {
    char cur_series = (*sp2i).get_series();
    if (cur_series == 'b' || cur_series == 'y') {
      char cur_num = (*sp2i).get_number();
      if (cur_series == 'b') {
        (*sp2i).add_mz(-bdeltas[cur_num] / (*sp2i).get_charge());
      } else {
        (*sp2i).add_mz(-(sdelta - bdeltas[this->get_lseq() - cur_num]) / (*sp2i).get_charge());
      }
    } else {
      // duplicate ms2 with delta value
      MS2Peak ms2 = MS2Peak((*sp2i));
      ms2.add_mz(-sdelta / ms2.get_charge());
      dpeaks.push_back(ms2);
    }
    ++sp2i;
  }
  // add new peaks to sp2
  this->ms2.insert(this->ms2.begin(), dpeaks.begin(), dpeaks.end());
}

void Spectrum::print() {
  for (std::size_t i = 0; i < ms2.size(); i++) {
    std::cout << title << '\t';
    std::cout << ms2[i].get_mz() << '\t';
    std::cout << (int)ms2[i].get_charge() << '\t';
    std::cout << ms2[i].get_intensity() << '\t';
    std::cout << ms2[i].get_series() << '\t';
    std::cout << (int)ms2[i].get_number() << std::endl;
  }
}

// return vector<Range> object with ms1 mz windows from vector<Spectrum> object
// TODO(dima) move to RangeC type
std::vector<Range> getMS1range(std::vector<Spectrum>* spv, bool query) {
  std::vector<Range> ms1r;
  for (std::size_t i = 0; i < spv->size(); i++) {
    /*for (int ch = 1; ch <= 3; ch++) {
      double cmz = spv->at(i).get_ms1().get_m(ch);
      Range r(cmz * (1 - ppm * 1e-6) , cmz * (1 + ppm * 1e-6), i, pid);
      ms1r.push_back(r);
    }*/
    Range r((MSPeak)spv->at(i).get_ms1(), i, query, i);
    ms1r.push_back(r);
  }

  return ms1r;
}

// Subset from vector of Spectrum by sequences or titles
// subset from vector<Spectrum> by titles vector
std::vector<Spectrum> subsetTitleSpectra(std::vector<Spectrum>*spv, 
  const std::vector<std::string>& titles, bool negate) {
  std::vector<Spectrum> nSpv = std::vector<Spectrum>(0);
  
  for(std::size_t i = 0; i < spv->size(); i++) {
    if (negate ^ inV(spv->at(i).get_title(), titles)) {
      nSpv.push_back(spv->at(i));
    }
  }
  return nSpv;
}
// subset from vector<Spectrum> by sequences vector
std::vector<Spectrum> subsetSeqSpectra(std::vector<Spectrum>*spv, 
  const std::vector<std::string>& seqs, bool negate) {
  std::vector<Spectrum> nSpv = std::vector<Spectrum>(0);
  for(std::size_t i = 0; i < spv->size(); i++) {
    if (negate ^ inV(spv->at(i).get_seq(), seqs)) {
      nSpv.push_back(spv->at(i));
    }
  }
  return nSpv;
}

// Spectral library class
// Spectra library consists of spectra for unique sequences. While creation
// from vector of spectra, spectra for identical sequences merged to one spectrum.
SpectraLib::SpectraLib() {
  slb = std::map<std::string, Spectrum>();
}

// add one spectrum to library
void SpectraLib::add_spectrum(const Spectrum& sp, double delta) {
  // find sequence of new spectra in spectra library
  std::map<std::string, Spectrum>::iterator oldsp = slb.find(sp.get_seq());
  // check for existande spectrum in lib with same sequence
  if(oldsp == slb.end()) {
    slb.insert(std::pair<std::string, Spectrum>(sp.get_seq(), sp));
  } else {
    // get all ms2 peaks for old and new spectra
    std::vector<MS2Peak> s1v = oldsp->second.get_ms2();
    std::vector<MS2Peak> s2v = sp.get_ms2();
    // create ranges with delta error window
    RangeC r1 = oldsp->second.getMS2range(1);
    RangeC r2 = sp.getMS2range(2);
    // overlap peaks mz
    olapdata nc = overlap2range(r1, r2, delta);
    // go through each overlaped peaks and calc mean intensity in lib spectrum
    for (std::size_t i = 0; i < nc.op1.size(); i++) {
      double oldint = s1v[nc.op1[i]].get_intensity();
      double newint = s2v[nc.op2[i]].get_intensity();
      oldsp->second.get_ms2p()->at(nc.op1[i]).set_intensity(oldint + newint);
    }
    // add new peaks to lib spectrum
    for (std::size_t i = 0; i < nc.np2.size(); i++) {
      oldsp->second.addMS2Peak(s2v[nc.np2[i]]);
    }
    // increase weight for spectrum
    oldsp->second.set_weight(oldsp->second.get_weight() + 1);
  }
}

// add several spectra to library from vector<Spectrum>
void SpectraLib::add_spectra(std::vector<Spectrum>* spv, double delta) {
  for (std::size_t i = 0; i < spv->size(); i++) {
    this->add_spectrum(spv->at(i), delta);
  }
}

// size of spectra library
int SpectraLib::get_size() { return slb.size(); }

// normalize each ms2 peak intensity to number of spectra for seq
// TODO(dima) remember that after normalization weights are the same
void SpectraLib::normalize_lib() {
  for(std::map<std::string, Spectrum>::iterator it = slb.begin();it != slb.end(); it++) {
    int ns = it->second.get_weight();
    for (std::size_t k = 0; k < it->second.get_ms2p()->size(); k++) {
      it->second.get_ms2p()->at(k).div_intensity(ns);
    }
  }
}

// get vector of spectra
std::vector<Spectrum> SpectraLib::get_spectra_vector() {
  std::vector<Spectrum> spv;
  for (std::map<std::string, Spectrum>::iterator it = slb.begin(); it != slb.end(); ++it) {
    spv.push_back(it->second);
  }
  return spv;
}


