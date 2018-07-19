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

	bool oscillate;
	std::string true_param_name;
	std::string true_L_name;

	//int value_i;
	//int true_value_i;
	//int true_L_i;

	float value_f;
	float true_value_f;
	float true_L_f;

	double value_d;
	double true_value_d;
	double true_L_d;

	branch_var(std::string n, std::string t, std::string a) : name(n), type(t), associated_hist(a) {oscillate=false;};
	virtual void* getValue(){};
	virtual void* getTrueValue(){};
	virtual void* getTrueL(){};

	int setOscillate(bool inbool){ oscillate = inbool;};
	bool getOscillate(){ return oscillate;};

};

/*struct branch_var_i: public branch_var{
	branch_var_i(std::string n, std::string t, std::string a) : branch_var(n,t,a) {value_i=0;true_value_i=0; true_L_i =0;};
	void* getValue(){ return &value_i;}
	void* getTrueValue(){ return &true_value_i;}
	void* getTrueL(){ return &true_L_i;}
};
*/

struct branch_var_d: public branch_var{
	branch_var_d(std::string n, std::string t, std::string a) : branch_var(n,t,a) {value_d=0;true_value_d=0; true_L_d = 0;};
	void* getValue(){ return &value_d;}
	void* getTrueValue(){ return &true_value_d;}
	void* getTrueL(){ return &true_L_d;}
};

struct branch_var_f: public branch_var{
	branch_var_f(std::string n, std::string t, std::string a) : branch_var(n,t,a) {value_f=0;true_value_f=0; true_L_f = 0;};
	void* getValue(){ return &value_f;}
	void* getTrueValue(){ return &true_value_f;}
	void* getTrueL(){ return &true_L_f;}
};

#endif
