#ifndef BRANCHVARIABLE_H_
#define BRANCHVARIABLE_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <time.h>


struct branch_var{
	std::string name;
	std::string type;
	std::string associated_hist;
	double value;	
	branch_var(std::string n, std::string t, std::string a) : name(n), type(t), associated_hist(a) {value=0;};	
	virtual void* getValue(){};
};

struct branch_var_i: public branch_var{
	int value_i; 
	branch_var_i(std::string n, std::string t, std::string a) : branch_var(n,t,a) {value_i=0;};	
	void* getValue(){ return &value_i;}
};

struct branch_var_d: public branch_var{
	double value_d; 
	branch_var_d(std::string n, std::string t, std::string a) : branch_var(n,t,a) {value_d=0;};	
	void* getValue(){ return &value_d;}
};

struct branch_var_f: public branch_var{
	float value_f; 
	branch_var_f(std::string n, std::string t, std::string a) : branch_var(n,t,a) {value_f=0;};	
	void* getValue(){ return &value_f;}
};

#endif
