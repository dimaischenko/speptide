// Parent MSPeak class and
// child MS1 and MS2 peak classes
// and Ranges models (for intersection)

#ifndef MS1MS2_H
#define MS1MS2_H

#include <vector>
#include <map>

// Class for MS peak (parent class)
class MSPeak {
protected:
  double mz; // mz ratio
  double intensity; // intensity
  char charge; // charge
  bool plus; // plus or minus charge state
public:
  MSPeak(double _mz = 0, double _intensity = 0,
    char _charge = 1, bool _plus = true);
  MSPeak(const MSPeak& msp);

  // get block
  double get_mz();
  double get_intensity() const;
  char get_charge();
  bool is_plus();
  double get_m();

  // set block
  void set_mz(double _mz);
  void set_m(double _m);
  void set_intensity(double _intensity);
  void set_charge(char _charge);
  void set_plus(bool _plus);

  // methods
  void div_intensity(int idiv);
};

// Class with MS1 peak information
// contain all needed information about MS1 peak
class MS1Peak : public MSPeak {
private:
  double rt; // retention time

public:
  MS1Peak(double _mz = 0, double _rt = 0, double _intensity = 0,
    char _charge = 2, bool _plus = true);
  MS1Peak(const MS1Peak& ms1);

  // set get block
  double get_rt();
  void set_rt(double _rt);
};

// Class for MS2 peak object
// contain all needed information about MS2 ion peak
class MS2Peak : public MSPeak {
private:
  char series; // ion series {b, y, B, Y, P, I ...}
  char number; // ion number {1, 2, 3, ...} in series

public:
  MS2Peak(double _mz = 0, double _intensity = 0, int _charge = 1,
    bool _plus = true, char _series = 'n', char _number = 0);
  MS2Peak(const MS2Peak &ms2p);

  // set get block
  char get_series() const;
  char get_number();
  void set_series(char _series);
  void set_number(char _number);
  void set_intensity(double _intensity);
  // set intensity and transform
  // several methods for intensity transformation
  // TODO(dima) some hash table for different algorithms and good names
  void set_intensity(double _intensity, char trAlg, double iConst);
  // add mz value to ms2 mz
  void add_mz(double addmz);
};

// operator for comparsion MS2 peaks (by intensity)
bool operator< (MS2Peak a, MS2Peak b);

// Range class is helpfull class for intersection peaks
// it contain MS peak object its id and mid (parent peak id)
// parent peak can produce several peaks with different charges
// so these peaks have different ids but the same mid
// bool query flag - indicates is it query or subject peak
// we need it for intersection
class Range {
private:
  // MSPeak object
  MSPeak mspeak;
  // range id
  int id;
  // ms id (for one ms we can create several ranges with diff charges)
  int mid;
  // range type (query or not)
  // use when intersect two vector<Range> to distinguish ranges from
  // different spectrum
  bool query;
public:
  Range();
  Range(MSPeak _mspeak, int _id, bool _query, int _mid);

  // get block
  MSPeak get_mspeak();
  int get_id();
  bool isQuery();
  int get_mid();
  double get_m();
  double get_intensity() const;

  // set block
  void set_mspeak(MSPeak _mspeak);
  void set_id(int _id);
  void set_mid(int _mid);
  void set_query(bool _query);
};

// operator for comparsion Ranges (order by masses)
bool operator< (Range a, Range b);
// operator for equal ranges
bool operator== (Range a, Range b);

// Range collection class has vector of Ranges
// and their information about parent peaks
class RangeC {
private:
  std::vector<Range> rc; // vector with ranges
  std::vector<int> mid2first; // for -i mid first index of child peak
  std::map<int, int> midinter; // map for count intersection for each mid
public:
  RangeC();
  void add_range(Range _range);
  int get_midl();
  int get_size() const;
  double get_i_intens(int i) const;
  void add_mid2first(int _m2f);
  void add_midinter(int _mid);
  void inc_midinter(int _mid);
  std::vector<Range> get_ranges() const;
  std::vector<Range>* get_rangesp();
  std::vector<int>* get_mid2firstp();
  int get_mid2first_i(int i) const;
  std::map<int, int>* get_midinter_p();
  std::map<int, int> get_midinter() const;
};

// structure for overlapping range (overlap2range function) result data
struct olapdata {
  std::vector<int> op1, op2, np1, np2;
};
// Function overlap two vector<Range> objects
// and return structure with information (indexes)
// with overlapped and not overlapped peaks
// TODO(dima) simplify
olapdata overlap2range(const RangeC& r1, const RangeC& r2,
  double delta, double ppm = 0);

#endif
