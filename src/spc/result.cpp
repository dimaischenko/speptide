// Class for algorithms results processing
// Subsetting with thresholds, printing etc

#include "result.h"
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <set>

// Result parent class
// Something like data.frame representation (list of vectors)
Result::Result() {
  title1 = std::vector<std::string>(0);
  title2 = std::vector<std::string>(0);
  seq2 = std::vector<std::string>(0);
  angle = std::vector<double>(0);
  m1 = std::vector<double>(0);
  m2 = std::vector<double>(0);
}

// get unique title1 from results
std::vector<std::string> Result::get_uniq_title1() {
  std::set<std::string> stitle1(title1.begin(), title1.end());
  std::vector<std::string> utitle1(stitle1.begin(), stitle1.end());
  return utitle1;
}

// get unique title1 from results that have angle greater than ath
std::vector<std::string> Result::get_uniq_title1ang(double _ath) {
  std::set<std::string> stitle1;
  for (std::size_t i = 0; i < title1.size(); i++) {
    if(angle[i] >= _ath) stitle1.insert(title1[i]);
  }
  std::vector<std::string> utitle1(stitle1.begin(), stitle1.end());
  return utitle1;
}

// get unique seq2 from results
std::vector<std::string> Result::get_uniq_seq2() {
  std::set<std::string> sseq2(seq2.begin(), seq2.end());
  std::vector<std::string> useq2(sseq2.begin(), sseq2.end());;
  return useq2;
}

// Result ident class
ResultIdent::ResultIdent() : Result() {};

// add new result line
void ResultIdent::add_result(std::string _title1, std::string _title2,
  std::string _seq2, double _angle) {
  this->title1.push_back(_title1);
  this->title2.push_back(_title2);
  this->seq2.push_back(_seq2);
  this->angle.push_back(_angle);
}

// print results
void ResultIdent::print() {
  for(std::size_t i = 0; i < this->title1.size(); i++) {
    std::cout << title1[i] << '\t';
    std::cout << title2[i] << '\t';
    std::cout << seq2[i] << '\t';
    std::cout << angle[i];
    std::cout << std::endl;
  }
}

// Result sap class
ResultSap::ResultSap() : Result() {
  this->ami1 = std::vector<std::string>(0);
  this->ami2 = std::vector<std::string>(0);
  this->pos = std::vector<std::string>(0);
}

// add new result line
void ResultSap::add_result(const std::string& _title1, const std::string& _title2,
    double _m1, double _m2, 
    const std::string& _seq2, double _angle, const std::string& _ami1, const std::string& _ami2,
    const std::string& _pos, int _qa, int _sa, int _o) {
  this->title1.push_back(_title1);
  this->title2.push_back(_title2);
  this->m1.push_back(_m1);
  this->m2.push_back(_m2);
  this->seq2.push_back(_seq2);
  this->angle.push_back(_angle);
  this->ami1.push_back(_ami1);
  this->ami2.push_back(_ami2);
  this->pos.push_back(_pos);
  this->qa.push_back(_qa);
  this->sa.push_back(_sa);
  this->o.push_back(_o);
}

// print results
void ResultSap::print() {
  for(std::size_t i = 0; i < this->title1.size(); i++) {
    std::cout << title1[i] << '\t';
    std::cout << title2[i] << '\t';
    std::cout << m1[i] << '\t';
    std::cout << m2[i] << '\t';
    std::cout << pos[i] << '\t';
    std::cout << ami1[i] << '\t';
    std::cout << ami2[i] << '\t';
    //std::cout << qa[i] << '\t';
    //std::cout << sa[i] << '\t';
    //std::cout << o[i] << '\t';
    std::cout << seq2[i] << '\t';
    std::cout << angle[i];
    std::cout << std::endl;
  }
}
