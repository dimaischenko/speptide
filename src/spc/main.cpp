/*-------------------------------------------------------
          SPEPTIDE

Application for align two spectra sets (two mgf) and find
amino acids substitutions. Calculate spectral angles
and return possible substitutions.

author: Dima Ischenko
email: ischenko.dmitry@gmail.com
version: 0.92
date start: 29.07.2014
current version date: 24.03.2016
---------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <map>
#include <stdexcept>
#include "time.h"

#include "../ini/iniparser.h"
#include "spectra.h"
#include "mparser.h"
#include "result.h"
#include "scoring.h"

using namespace std;

int main (int argc, char *argv[])
{
  try {
  
  // check input arguments and print help
  if(argc != 4) {
    cout << "  speptide version 0.92 by Dima Ischenko (ischenko.dmitry@gmail.com)" << endl;
    cout << endl;
    cout << "  Usage:" << endl;
    cout << "   speptide <exp mgf> <db mgf> <config>" << endl;
    cout << endl;
    cout << "  Inputs:" << endl;
    cout << "   <exp mgf>    (string)  experiment spectra mgf file." << endl;
    cout << "   <db mgf>     (string)  database spectra mgf file (with SEQ tag)." << endl;
    cout << "   <config>     (string)  ini file with params for search." << endl;
    cout << endl;
    cout << "  Output (tab separated format):" << endl;
    cout << "   [exp spectra] [db spectra] [position] [exp ami] [db ami] [db seq] [cos(theta)]" << endl;
    cout << endl;
    cout << "  Read doc/speptide.pdf for detailed information.";
    cout << endl;
    exit(1);
  }

  // set mgf filenames
  string mgf1F = (string)argv[1];
  string mgf2F = (string)argv[2];
  
  // load ini file and set params
  dictionary* ini;
  ini = iniparser_load(argv[3]);
  if (ini == NULL) return -1;

  // annot params
  double byDelta = iniparser_getdouble(ini, "anot:da", 0.5);

  // ident params
  double i_erms1 = iniparser_getdouble(ini, "ident:ms1er", 10);
  bool i_isPpm = iniparser_getboolean(ini, "ident:ms1ppm", true);
  double i_erDa = iniparser_getdouble(ini, "ident:da", 0.5);
  double i_pdiv = iniparser_getdouble(ini, "ident:n", 100);
  char i_trAlg = iniparser_getstring(ini, "ident:trans", (char *)"b")[0];
  double i_ath = iniparser_getdouble(ini, "ident:cos", 0);
  bool i_norm = iniparser_getboolean(ini, "ident:norm", true);
  double i_iConst = iniparser_getdouble(ini, "ident:const", 0);
  double i_ath05 = iniparser_getdouble(ini, "ident:cos05", 0);

  // sap params
  double erms1 = iniparser_getdouble(ini, "sap:ms1er", 10);
  bool isPpm = iniparser_getboolean(ini, "sap:ms1ppm", true);
  double erDa = iniparser_getdouble(ini, "sap:da", 0.5);
  double pdiv = iniparser_getdouble(ini, "sap:n", 40);
  char trAlg = iniparser_getstring(ini, "sap:trans", (char *)"b")[0];
  double modath = iniparser_getdouble(ini, "sap:mcos", 0);
  double ath = iniparser_getdouble(ini, "sap:cos", 0);
  bool norm = iniparser_getboolean(ini, "sap:norm", true);
  double iConst = iniparser_getdouble(ini, "sap:const", 0);
  bool addions = iniparser_getboolean(ini, "sap:addions", false);
  double refdiv = iniparser_getdouble(ini, "sap:refdiv", 1);
  bool fIdent = iniparser_getboolean(ini, "sap:fident", true);
  bool fMod = iniparser_getboolean(ini, "sap:fmod", true);
  
  // free ini
  iniparser_freedict(ini);

  // loader
  MParser msp = MParser();

  // load amino acid masses
  map<char, double> masses = msp.defaultMasses();
  // load vector with deltas and sort
  vector<ADelta> ms1deltas = msp.defaultDeltas();
  sort(ms1deltas.begin(), ms1deltas.end());

  // load and parse mgf data. one set is query(annotation - false)
  // another (second) database with seq tab(annotation - true)
  vector<Spectrum> spv1 = msp.loadMGF(mgf1F, masses, byDelta, false, false);
  vector<Spectrum> spv2 = msp.loadMGF(mgf2F, masses, byDelta, true, addions);
  
  // output double precision
  cout << setprecision(7);
  
  // create calculator class
  Scoring score = Scoring();

  vector<Spectrum> spv1fin = vector<Spectrum>();
  vector<Spectrum> spv2fin = vector<Spectrum>();

  // Alogrithm pipeline 
  if (fIdent) {
    // find identical
    ResultIdent nci = score.comp2spvAngle(spv1, spv2,
        i_erms1, i_isPpm, i_erDa, i_pdiv, i_ath05,
        i_norm, i_iConst, i_trAlg);
    
    // subset not identical from experiment
    vector<Spectrum> spv1ni = subsetTitleSpectra(&spv1,
        nci.get_uniq_title1ang(i_ath), true);
    // subset not identical from spectral library
    spv2fin = subsetSeqSpectra(&spv2,
        nci.get_uniq_seq2(), true);

    if (fMod) {
      // subset identical from library for modifications search
      vector<Spectrum> spv2sb = subsetSeqSpectra(&spv2,
          nci.get_uniq_seq2());
      
      // find modifications
      if (spv1ni.size() > 0 && spv2sb.size() > 0) {
        ResultSap ncmod = score.comp2spvAngleAap(spv1ni, spv2sb,
            erms1, isPpm, erDa, pdiv, ms1deltas, masses, modath,
            norm, iConst, trAlg, refdiv);
        spv1fin = subsetTitleSpectra(&spv1ni,
            ncmod.get_uniq_title1(), true);
      } else {
        spv1fin = spv1ni;
      }
    } else {
      spv1fin = spv1ni;
    }
  } else {
    spv1fin = spv1;
    spv2fin = spv2;
  }
  
  // final sap search
  if (spv1fin.size() > 0 && spv2fin.size() > 0) {
    ResultSap ncfin = score.comp2spvAngleAap(spv1fin, spv2fin,
        erms1, isPpm, erDa, pdiv, ms1deltas, masses, ath,
        norm, iConst, trAlg, refdiv);
    // print final result
    ncfin.print();
  }
  } catch(exception& e) {
    cerr << "Error: " << e.what() << endl;
    exit(1);
  }

  // ola
  return 0;
}
