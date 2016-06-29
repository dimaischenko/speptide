// Class for algorithms results processing
// Subsetting with thresholds, printing etc

#ifndef RESULT_H
#define RESULT_H

#include <vector>
#include <string>

// Result parent class
// Something like data.frame representation (list of vectors)
class Result {
protected:
  // some members for alignment
  std::vector<std::string> title1;
  std::vector<std::string> title2;
  std::vector<std::string> seq2;
  std::vector<double> m1;
  std::vector<double> m2;
  std::vector<double> angle;
  std::vector<int> qa;
  std::vector<int> sa;
  std::vector<int> o;
public:
  Result();

  // get unique ids and sequences from result
  std::vector<std::string> get_uniq_title1();
  std::vector<std::string> get_uniq_seq2();
  std::vector<std::string> get_uniq_title1ang(double ath);

  void print();
};

// Result child class for identical algorithm
class ResultIdent : public Result {
public:
  ResultIdent();

  void add_result(std::string _title1, std::string _title2,
    std::string _seq2, double _angle);

  void print();
};

// Result child class for sap algorithm
class ResultSap : public Result {
private:
  std::vector<std::string> ami1;
  std::vector<std::string> ami2;
  std::vector<std::string> pos;
public:
  ResultSap();

  void add_result(std::string _title1, std::string _title2,
      double _m1, double _m2,
      std::string _seq2, double _angle, std::string _ami1, std::string _ami2,
      std::string _pos, int _qa, int _sa, int _o);

  void print();
};

#endif
