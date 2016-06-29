// Class for calcuations angles between spectra and libs

#include<vector>
#include<string>
#include<map>
#include<sstream>
#include<algorithm>
#include<cmath>

#include "ms1ms2.h"
#include "spectra.h"
#include "result.h" 
#include "scoring.h"

// calculate cos of angle between two vectors
// if one of vector equal zero vector (0, 0, 0, 0) return 0
double Scoring::sangle(const std::vector<double>& v1, const std::vector<double>& v2) {
  // minimal vector length
  int vl = std::min(v1.size(), v2.size());
  if (vl == 0) return 0;

  // calc lengths and scalar
  double s1 = 0, s2 = 0, s12 = 0;
  for (int i = 0; i < vl; i++) {
    s1 += v1[i] * v1[i];
    s2 += v2[i] * v2[i];
    s12 += v1[i] * v2[i];
  }

  if ((s1 * s2) == 0) return 0;

  return (s12 / sqrt(s1 * s2));
}

// int to string
std::string intToString(int i)
{
  std::stringstream ss;
  std::string s;
  ss << i;
  s = ss.str();
  return s;
}

// calculate distance between two spectra ms2 data
double Scoring::comp2speci (const RangeC& r1, const RangeC& r2, double delta, double pdiv) {
  double dist = 0;
  // subset top intensity ms2 peaks
  // TODO(dima) subset early
  //sp1.filter_ms2(pdiv);
  //sp2.filter_ms2(pdiv);

  // create ranges with delta error window
  //RangeC r1 = sp1.getMS2range(true);
  //RangeC r2 = sp2.getMS2range(false);

  // overlap peaks mz
  olapdata nc = overlap2range(r1, r2, delta);

  // create vectors with intensities
  std::vector<double> vi1;
  std::vector<double> vi2;
  for (std::size_t i = 0; i < nc.op1.size(); i++) {
    vi1.push_back(r1.get_i_intens(nc.op1[i]));
    vi2.push_back(r2.get_i_intens(nc.op2[i]));
  }
  for (std::size_t i = 0; i < nc.np1.size(); i++) {
    vi1.push_back(r1.get_i_intens(nc.np1[i]));
    vi2.push_back(0);
  }
  for (std::size_t i = 0; i < nc.np2.size(); i++) {
    vi1.push_back(0);
    vi2.push_back(r2.get_i_intens(nc.np2[i]));
  }
  // return dist
  dist = sangle(vi1, vi2);
  return dist;
}

// Calculation of distance between two spectra
// first spectrum(sp1) is query spectrum. second spectrum (sp2) is referent spectrum
// function apply some modification to referent spectrum (e.g. shift b/y ions according to delta)
// select top (mz / pdiv) MS2 peaks from query spectrum by intensities
// intersect them with referent MS2 peaks and calculate spectral angle
// return structure with angle, amounts of intersected peaks
compdata Scoring::comp2spec (Spectrum sp1, Spectrum sp2, double delta, double pdiv,
  std::vector< std::pair<int, double> > dpos) {
  double dist = 0;
  // subset top peaks in query
  //shit referent ions
  sp2.shiftBY(dpos);
  // TODO(dima) subset early
  sp1.filter_ms2(sp2.ms2length() * pdiv);

  // create ranges with delta error window
  // charges vecor
  //std::vector<char> charges = std::vector<char>(0);
  //charges.push_back(1); charges.push_back(2);

  RangeC r1 = sp1.getMS2range(true);
  RangeC r2 = sp2.getMS2range(false);
  // overlap peaks mz
  olapdata nc = overlap2range(r1, r2, delta);
  // create vectors with intensities
  std::vector<double> vi1 = std::vector<double>(0);
  std::vector<double> vi2 = std::vector<double>(0);
  
  double v1 = 0; double v2 = 0; double v12 = 0;
  //TODO (insert angle calculation for increase speed)
  for (std::size_t i = 0; i < nc.op1.size(); i++) {
    //vi1.push_back(r1.get_i_intens(nc.op1[i]));
    //vi2.push_back(r2.get_i_intens(nc.op2[i]));
    double i1 = r1.get_i_intens(nc.op1[i]);
    double i2 = r2.get_i_intens(nc.op2[i]);
    v12 += i1 * i2;
    v1 += i1 * i1;
    v2 += i2 * i2;
  }
  for (std::size_t i = 0; i < nc.np1.size(); i++) {
    //vi1.push_back(r1.get_i_intens(nc.np1[i]));
    //vi2.push_back(0);
    double i1 = r1.get_i_intens(nc.np1[i]);
    v1 += i1 * i1;
  }
  for (std::size_t i = 0; i < nc.np2.size(); i++) {
    //vi1.push_back(0);
    //vi2.push_back(r2.get_i_intens(nc.np2[i]));
    double i2 = r2.get_i_intens(nc.np2[i]);
    v2 += i2 * i2;
  }
  
  // TODO(dima) this is bad part
  // return dist
  //dist = sangle(vi1, vi2);
  dist = (v1 * v2) == 0 ? 0 : v12 / sqrt(v1 * v2); 

  // fill compdata struct
  compdata res;
  res.angle = dist;
  
  res.qa = r1.get_size();
  res.sa = r2.get_size();
  res.o = nc.op1.size();
  /*
  sort(nc.op1.begin(), nc.op1.end()); sort(nc.op2.begin(), nc.op2.end());
  nc.op1.erase(unique(nc.op1.begin(), nc.op1.end()), nc.op1.end());
  nc.op2.erase(unique(nc.op2.begin(), nc.op2.end()), nc.op2.end());
  res.qo = nc.op1.size();
  res.so = nc.op2.size();
  */
  return res;
}


// function compare two spectra object my therir ms1 and ms2 data
ResultIdent Scoring::comp2spvAngle (std::vector<Spectrum> sp1v, std::vector<Spectrum> sp2v,
  double value, bool isPpm, double delta, double pdiv, double ath, bool norm,
  double iConst, char trAlg) {
  std::vector<RangeC> rc1v = std::vector<RangeC>(0);
  // normalize spectra
  for (std::size_t spi = 0; spi < sp1v.size(); spi++) {
    sp1v[spi].normalize(norm, trAlg, iConst);
    sp1v[spi].filter_ms2(sp1v[spi].top_m_peaks(pdiv));
    rc1v.push_back(sp1v[spi].getMS2range(true));
  }
  std::vector<RangeC> rc2v = std::vector<RangeC>(0);
  // normalize spectra
  for (std::size_t spi = 0; spi < sp2v.size(); spi++) {
    sp2v[spi].normalize(norm, trAlg, iConst);
    sp2v[spi].filter_ms2(sp2v[spi].top_m_peaks(pdiv));
    rc2v.push_back(sp2v[spi].getMS2range(false));
  }
  // result object
  ResultIdent nc;
  // crate ranges
  std::vector<Range> r1 = getMS1range(&sp1v, true);
  std::vector<Range> r2 = getMS1range(&sp2v, false);
  // concat vectors of ranges and sort it
  std::vector<Range> r3;
  r3.reserve(r1.size() + r2.size());
  r3.insert(r3.end(), r1.begin(), r1.end());
  r3.insert(r3.end(), r2.begin(), r2.end());
  std::sort(r3.begin(), r3.end());
  for (std::size_t i = 0; i < r3.size() - 1; i++) {
    std::size_t k = i + 1;
    // delta construction (in case if ppm is selected)
    double ms1delta = (isPpm) ? std::abs(r3[i].get_m() * value * 1e-6) : value;
    while(r3[k].get_m() <= (r3[i].get_m() + ms1delta) && k < r3.size()) {
      if (r3[k].isQuery() != r3[i].isQuery()) {
        int qindex = i;
        if (r3[k].isQuery()) {
          qindex = k;
        }
        
        int i1 = r3[qindex].get_id();
        int i2 = r3[i + k - qindex].get_id();
        double ang = comp2speci(rc1v[i1], rc2v[i2], delta, pdiv);
        
        if (ang >= ath) {
          nc.add_result(sp1v[i1].get_title(), sp2v[i2].get_title(),
            sp2v[i2].get_seq(), ang);
        }
      }
      k++;
    }
  }

  return nc;
}


// function compare two vector<Specta> object my their ms1 and ms2 data
// second vector<Spectrum> considered as a reference (with sequences and annotated peaks)
// we find MS1 deltas between query and referent spectra and check them
// whether they can be candidates for amino acid substitution. MS1 delta beints to {ADeltas}
// next step. for each candidate pair we find all possible amino acid substitution
// positions (by reference sequence) and calculate spectral angle for each realised
// substitution (with shifting referent peaks)
ResultSap Scoring::comp2spvAngleAap (std::vector<Spectrum> sp1v,
  std::vector<Spectrum> sp2v, double value, bool isPpm, double delta, double pdiv,
  std::vector<ADelta> deltas, std::map<char, double> masses, double ath,
  bool norm, double iConst, char trAlg, double refdiv) {
  
  // get max delta
  double mdelta = calc_max_delta(deltas); 
  // create RangesC vector in this point to increase speed
  // std::vector<RangeC> rc1v = std::vector<RangeC>(0);
  // normalize spectra
  for (std::size_t spi = 0; spi < sp1v.size(); spi++) {
    sp1v[spi].normalize(norm, trAlg, iConst);
    //sp1v[spi].filter_ms2(sp1v[spi].top_m_peaks(pdiv));
    //rc1v.push_back(sp1v[spi].getMS2range(true));
  }
  // normalize spectra
  for (std::size_t spi = 0; spi < sp2v.size(); spi++) {
    sp2v[spi].normalize(norm, trAlg, iConst);
    // TODO(dima). remember about it. remove not annotated peaks
    sp2v[spi].filter_ms2(sp2v[spi].top_m_peaks(refdiv));
    sp2v[spi].removeNA();
  }

  // result table
  ResultSap nc = ResultSap();
  // crate ranges
  std::vector<Range> r1 = getMS1range(&sp1v, true);
  std::vector<Range> r2 = getMS1range(&sp2v, false);
  // concat vectors of ranges and sort it
  std::vector<Range> r3;
  r3.reserve(r1.size() + r2.size());
  r3.insert(r3.end(), r1.begin(), r1.end());
  r3.insert(r3.end(), r2.begin(), r2.end());
  sort(r3.begin(), r3.end());
  // let's start
  for (std::size_t i = 0; i < r3.size() - 1; i++) {
    std::size_t k = i + 1;
    double ms1delta = r3[k].get_m() - r3[i].get_m();
    while( ms1delta < mdelta && k < r3.size()) {
      ms1delta = r3[k].get_m() - r3[i].get_m();
      if (r3[k].isQuery() != r3[i].isQuery()) {
        int qindex = i;
        if (r3[k].isQuery()) {
          qindex = k;
          ms1delta = -ms1delta;
        }
        // in range check
        std::vector<ADelta> inrd = inRange(r3[qindex].get_m(),
          r3[i + k - qindex].get_m(), deltas, value, isPpm);
        if (inrd.size() == 0) { k++; continue; }
        // indexes
        int i1 = r3[qindex].get_id();
        int i2 = r3[i + k - qindex].get_id();
        // seq
        std::string sp2seq = sp2v.at(i2).get_seq();
        // for each delta
        for (std::size_t dli = 0; dli < inrd.size(); dli++) {
          std::string newc = inrd[dli].get_ami1(ms1delta > 0);
          std::string oldc = inrd[dli].get_ami2(ms1delta > 0);
          // get availible positions for deltas
          std::vector< std::vector<int> > pos = sp2v.at(i2).getDeltaPos(inrd[dli], ms1delta); 
          // all amies in seq
          if (pos.size() == 0 || pos[0].size() != oldc.size()) continue;
          for (std::size_t posi = 0; posi < pos.size(); posi++) {
            std::vector< std::pair<int, double> > dpos = std::vector< std::pair<int, double> >(0);
            for (std::size_t ai = 0; ai < pos[posi].size(); ai++) {
              dpos.push_back(std::pair<int, double>(pos[posi][ai], masses[oldc[ai]] - masses[newc[ai]]));
            }

            compdata compd = comp2spec(sp1v.at(i1), sp2v.at(i2), delta, pdiv, dpos);
            // add result section
            if (compd.angle >= ath) {
              // concat positions to space separated string
              std::string printpos = intToString(dpos[0].first);
              for (std::size_t dpi = 1; dpi < dpos.size(); dpi++) {
                printpos += " ";
                printpos += intToString(dpos[dpi].first);
              }

              nc.add_result(sp1v.at(i1).get_title(), sp2v.at(i2).get_title(),
                sp1v.at(i1).get_ms1().get_m(),
                sp2v.at(i2).get_ms1().get_m(),
                sp2v.at(i2).get_seq(), compd.angle, newc, oldc, printpos,
                compd.qa, compd.sa, compd.o);
            }
          }
        }

      }

      k++;
    }
  }

  return nc;
}

