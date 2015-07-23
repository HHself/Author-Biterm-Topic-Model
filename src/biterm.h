#ifndef _BITERM_H
#define _BITERM_H

#include <cmath>
#include <string>
#include <sstream>
#include <cassert>

using namespace std;

class Biterm {
private:
  int wi;
  int wj;
  int z; // topic assignment
  int a; // author assignment
  
public: 
  Biterm(int w1, int w2): z(-1) {
	wi = min(w1, w2);
	wj = max(w1, w2);
        a  = -1;
  }

  // s format:wi   wj    z    a
  Biterm(string s) {
	istringstream iss(s);
	iss >> wi >> wj >> z >> a;
  }
  
  int get_wi() const {return wi;}
  int get_wj() const {return wj;}
  
  int get_z() const {return z;}
  int get_a() const {return a;}
  void set_z(int k) {z = k;}
  void set_a(int au) {a = au;}
  void reset_z() {z = -1;}
  void reset_a() {a = -1;}
  
  string str() const {
	ostringstream os;
	os << wi << '\t' << wj << '\t' << '\t' << z << '\t' << a;
	return os.str();
  }  
};

#endif
