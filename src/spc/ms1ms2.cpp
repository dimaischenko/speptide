// Parent MSPeak class and
// child MS1 and MS2 peak classes
// and Ranges models (for intersection)

#include "ms1ms2.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

#define mH 1.007825

// Class for MS peak (parent class)
MSPeak::MSPeak(double _mz, double _intensity,
  char _charge, bool _plus) {
  
  mz = _mz;
  intensity = _intensity;
  charge = _charge;
  plus = _plus;
}
MSPeak::MSPeak(const MSPeak& ms1) {
  mz = ms1.mz;
  intensity = ms1.intensity;
  charge = ms1.charge;
  plus = ms1.plus;
}
// set get block
double MSPeak::get_mz() {
  return mz;
}
double MSPeak::get_intensity() const {
  return intensity;
}
char MSPeak::get_charge() {
  return charge;
}
bool MSPeak::is_plus() {
  return plus;
}
double MSPeak::get_m() {
  return ((mz - mH) * charge);
}
void MSPeak::set_charge(char _charge) {
  charge = _charge;
}
void MSPeak::set_mz(double _mz) {
  mz = _mz;
}
void MSPeak::set_m(double _m) {
  mz = (_m / charge) + mH;
}
void MSPeak::set_plus(bool _plus) {
  plus = _plus;
}
void MSPeak::set_intensity(double _intensity) {
  intensity = _intensity;
}
void MSPeak::div_intensity(int idiv) {
  intensity /= idiv;
}

// Class with MS1 peak information
// contain all needed information about MS1 peak
MS1Peak::MS1Peak(double _mz, double _rt, double _intensity, char _charge,
  bool _plus) : MSPeak(_mz, _intensity, _charge, _plus) {
  rt = _rt;
} 

MS1Peak::MS1Peak(const MS1Peak& ms1p) {
  mz = ms1p.mz;
  intensity = ms1p.intensity;
  charge = ms1p.charge;
  plus = ms1p.plus;
  rt = ms1p.rt;
}

double MS1Peak::get_rt() {
  return rt;
}
void MS1Peak::set_rt(double _rt) {
  rt = _rt;
}

// Class for MS2 peak object
// contain all needed information
MS2Peak::MS2Peak(double _mz, double _intensity, int _charge, 
  bool _plus, char _series, char _number) : MSPeak(_mz, _intensity,
                _charge, _plus) {

  series = _series;
  number = _number;
}

MS2Peak::MS2Peak(const MS2Peak &ms2p) {
  mz = ms2p.mz;
  intensity = ms2p.intensity;
  charge = ms2p.charge;
  plus = ms2p.plus;
  series = ms2p.series;
  number = ms2p.number;
}
char MS2Peak::get_series() {
  return series;
}
char MS2Peak::get_number() {
  return number;
}
void MS2Peak::set_series(char _series) {
  series = _series;
}
void MS2Peak::set_number(char _number) {
  number = _number;
}
void MS2Peak::set_intensity(double _intensity) {
  intensity = _intensity;
}
// set intensity and transform
// several methods for intensity transformation
// TODO(dima) some hash table for different algorithms and good names
void MS2Peak::set_intensity(double _intensity, char trAlg, double iConst) {
  switch (trAlg) {
    // SQRT transform
    case 'a':
    intensity = iConst + sqrt(_intensity);
    break;
    // LN transrom
    case 'b':
    intensity = iConst + log(_intensity);
    break;
    // no transform
    case 'c':
    intensity = iConst + _intensity;
    break;
    default:
    intensity = 0;
  }
}
// add mz value to ms2 mz
void MS2Peak::add_mz(double addmz) {
  mz += addmz;
}
// operator for comparsion MS2 peaks (by intensity)
bool operator< (MS2Peak a, MS2Peak b) {
  return (a.get_intensity() > b.get_intensity());
}

// Range class is helpfull class for intersection peaks
// it contain MS peak object its id and mid (parent peak id)
// parent peak can produce several peaks with different charges
// so these peaks have different ids but the same mid
// bool query flag - indicates is it query or subject peak
// we need it for intersection
Range::Range() {
  id = 0;
  mspeak = MSPeak();
  query = false;
  mid = 0;
}
Range::Range(MSPeak _mspeak, int _id, bool _query, int _mid) {
  mspeak = _mspeak;
  id = _id;
  query = _query;
  mid = _mid;
}
// get block
MSPeak Range::get_mspeak() {
  return mspeak;
}
int Range::get_id() {
  return id;
}
int Range::get_mid() {
  return mid;
}
bool Range::isQuery() {
  return query;
}
double Range::get_m() {
  return mspeak.get_m();
}
double Range::get_intensity() const {
  return mspeak.get_intensity();
}
// set block
void Range::set_mspeak(MSPeak _mspeak) {
  mspeak = _mspeak;
}
void Range::set_id(int _id) {
  id = _id;
}
void Range::set_mid(int _mid) {
  mid = _mid;
}
void Range::set_query(bool _query) {
  query = _query;
}
// operator for comparsion Ranges (order by masses)
bool operator< (Range a, Range b) {
  return (a.get_m() < b.get_m());
}
// operator for equal ranges
bool operator== (Range a, Range b) {
  return (a.get_m() == b.get_m());
}

// Range collection class has vector of Ranges
// and their information about parent peaks
RangeC::RangeC() {
  rc = std::vector<Range>(0);
  mid2first = std::vector<int>(0);
  midinter = std::map<int, int>();
}
int RangeC::get_midl() {
  return mid2first.size();
}
int RangeC::get_size() const {
  return rc.size();
}
double RangeC::get_i_intens(int i) const {
  return rc[i].get_intensity();
}
void RangeC::add_mid2first(int _m2f) {
  mid2first.push_back(_m2f);
}
void RangeC::add_range(Range _range) {
  rc.push_back(_range);
}
void RangeC::add_midinter(int _mid) {
  midinter.insert(std::pair<int, int>(_mid, 0));
}
void RangeC::inc_midinter(int _mid) {
  midinter[_mid]++;
}
std::vector<Range> RangeC::get_ranges() const {
  return rc;
}
std::vector<Range>* RangeC::get_rangesp() {
  return &rc;
}
std::vector<int>* RangeC::get_mid2firstp() {
  return &mid2first;
}
int RangeC::get_mid2first_i(int i) const {
  return mid2first[i];
}
std::map<int, int>* RangeC::get_midinter_p() {
  return &midinter;
}
std::map<int, int> RangeC::get_midinter() const {
  return midinter;
}

// Function overlap two vector<Range> objects
// and return structure with information (indexes)
// with overlapped and not overlapped peaks
//TODO(dima) simplify
olapdata overlap2range(const RangeC& rc1, const RangeC& rc2,
  double delta, double ppm) {
  olapdata nc = olapdata();
  
  std::vector<Range> r1 = rc1.get_ranges();
  std::vector<Range> r2 = rc2.get_ranges(); 
  // concat vectors of ranges and sort it
  std::vector<Range> r3;
  r3.reserve(r1.size() + r2.size());
  r3.insert(r3.end(), r1.begin(), r1.end());
  r3.insert(r3.end(), r2.begin(), r2.end());
  std::sort(r3.begin(), r3.end());
  // vector with inter flag
  std::vector<bool> inter(r3.size());
  std::fill(inter.begin(), inter.end(), false);
  // get midinters
  std::map<int, int> rc1m = rc1.get_midinter();
  std::map<int, int> rc2m = rc2.get_midinter();
  // go through each
  for (std::size_t i = 0; i < r3.size(); i++) {
    // skip if intersected earlier (to best match)
    if (inter[i]) continue;
    std::size_t k = i + 1;
    // delta construction (in case if ppm is selected)
    delta = (ppm == 0) ? delta : std::abs(r3[i].get_m() * ppm * 1e-6);
    while(r3[k].get_m() <= (r3[i].get_m() + delta) && k < r3.size()) {
      if (r3[k].isQuery() != r3[i].isQuery()) {
        // set intersection flags
        inter[i] = true;
        inter[k] = true;
        int qindex = i;
        // is query checking
        if (r3[k].isQuery()) {
          qindex = k;
        }

        // get indexes
        int i1 = r3[qindex].get_id();
        int i2 = r3[i + k - qindex].get_id();
        nc.op1.push_back(i1);
        nc.op2.push_back(i2);
        // increase parent peak intersection counters
        rc1m[r3[qindex].get_mid()]++;
        rc2m[r3[i + k - qindex].get_mid()]++;
        // break (because to best match)
        break;
      }

      k++;
    }
    /*
    // check if not overlaps and fill ovelap data structure
    if (!inter[i]) {
      if(r3[i].isQuery()) {
        nc.np1.push_back(r3[i].get_id());
      } else {
        nc.np2.push_back(r3[i].get_id());
      }
    }*/
    // not overlapped data
    
  }
  
  for (std::map<int, int>::iterator it = rc1m.begin();
      it != rc1m.end(); ++it) {
    if (it->second == 0) {
      nc.np1.push_back(rc1.get_mid2first_i(it->first));
    }
  }
  for (std::map<int, int>::iterator it = rc2m.begin();
      it != rc2m.end(); ++it) {
    if (it->second == 0) {
      nc.np2.push_back(rc2.get_mid2first_i(it->first));
    }
  }
  return nc;
}
