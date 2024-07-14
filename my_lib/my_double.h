#pragma once
/*****************************************************************//**
 * \file   my_double.h
 * \brief  my_double class is same as double but has function to compute operation number 
 * 
 * \author 26259
 * \date   January 2023
 *********************************************************************/
#ifdef use_my_double
#include<math.h>
#include<iostream>
#include"my_float.h"
using namespace std;
class my_double_auxiliary_storage {
public:
	static unsigned long long operation_number;
	static unsigned long long compare_number;
};
unsigned long long my_double_auxiliary_storage::operation_number = 0;
unsigned long long my_double_auxiliary_storage::compare_number = 0;

// write it as val POD mode for speed, in release mode this has same speed as double
class my_double
{
private:
	double val;
public:
	static void init() {}		// do nothing

	/* define construction functions */
	my_double() = default;		// this is POD
	my_double(double b) :val(b) {}
	// allow converting my_float to my_double implicitly, but versa should be done explicitly, be careful
	my_double(my_float b) :val(b) {}
	my_double(float b) :val(b) {}
	my_double(unsigned long long b) :val(double(b)) {}
	my_double(int b) :val(b) {}
	my_double(bool b) :val(b) {}

	my_double operator - () const {
		return -val;		// add inverse
	}

	void operator += (const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number
		val += F2.val;
	}
	void operator -= (const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number
		val -= F2.val;
	}
	void operator *= (const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number
		val *= F2.val;
	}
	void operator /= (const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number
		val /= F2.val;
	}

	friend my_double operator + (const my_double& F1, const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1.val + F2.val;
	}
	friend my_double operator - (const my_double& F1, const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1.val - F2.val;						// result is same as '+' in my_double
	}
	friend my_double operator * (const my_double& F1, const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1.val * F2.val;
	}
	friend my_double operator / (const my_double& F1, const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1.val / F2.val;
	}

	friend ostream& operator << (ostream& out, my_double c) {
		out << c.val;
		return out;
	}
	friend istream& operator >> (istream& in, my_double& c) {
		in >> c.val;
		return in;
	}
	friend bool operator == (const my_double& F1, const my_double& F2) {
		return F1.val == F2.val;
	}
	friend bool operator !=  (const my_double& F1, const my_double& F2) {
		return F1.val != F2.val;
	}
	friend bool operator >(const my_double& F1, const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::compare_number++;
#endif // count_compare_number
		return F1.val > F2.val;
	}
	friend bool operator <(const my_double& F1, const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::compare_number++;
#endif // count_compare_number
		return F1.val < F2.val;
	}
	friend bool operator >=(const my_double& F1, const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::compare_number++;
#endif // count_compare_number
		return F1.val >= F2.val;
	}
	friend bool operator <=(const my_double& F1, const my_double& F2) {
#ifdef count_operation_number
		my_double_auxiliary_storage::compare_number++;
#endif // count_compare_number
		return F1.val <= F2.val;
	}
	explicit operator double() const {
		return val;
	}
	explicit operator float() const {
		return (float)val;
	}
	explicit operator int() const {
		return (int)val;
	}
	explicit operator unsigned long long() const {
		return (unsigned long long)val;
	}

	friend my_double sqrt(const my_double& c) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number

		return sqrt(c.val);
	}
	friend my_double exp(const my_double& c) {

#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number

		return exp(c.val);
	}
	friend my_double log(const my_double& c) {
#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number

		return log(c.val);
	}
	friend my_double pow(const my_double& c1, const my_double& c2) {

#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number

		return pow(c1.val, c2.val);
	}
	friend my_double log10(const my_double& c) {

#ifdef count_operation_number
		my_double_auxiliary_storage::operation_number++;
#endif // count_operation_number

		return log10(c.val);
	}
	friend my_double round(const my_double& c) {
		return round(c.val);
	}

	// triangle function is not consider here, you can easily add them
};

#endif // use_my_double
