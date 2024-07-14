#pragma once
/*****************************************************************//**
 * \file   my_float.h
 * \brief  my_float class is same as float but has function to compute operation number 
 * 
 * \author lilili
 * \date   April 2023
 *********************************************************************/

#ifdef use_my_double
#include<math.h>
#include<iostream>
using namespace std;
class my_float_auxiliary_storage {
public:
	static unsigned long long operation_number;
};
unsigned long long my_float_auxiliary_storage::operation_number = 0;

// write it as val POD mode for speed, in release mode this has same speed as float
class my_float
{
private:
	float val;
public:
	static void init() {}		// do nothing

	/* define construction functions */
	my_float() = default;		// this is POD
	my_float(double b) :val(float(b)) {}
	my_float(float b) :val(b) {}
	my_float(unsigned long long b) :val(float(b)) {}
	my_float(int b) :val(float(b)) {}
	my_float(bool b) :val(b) {}

	my_float operator - () const {
		return -val;		// add inverse
	}

	void operator += (const my_float& F2) {
#ifdef count_operation_number
		my_float_auxiliary_storage::operation_number++;
#endif // count_operation_number
		val += F2.val;
	}
	void operator -= (const my_float& F2) {
#ifdef count_operation_number
		my_float_auxiliary_storage::operation_number++;
#endif // count_operation_number
		val -= F2.val;
	}
	void operator *= (const my_float& F2) {
#ifdef count_operation_number
		my_float_auxiliary_storage::operation_number++;
#endif // count_operation_number
		val *= F2.val;
	}
	void operator /= (const my_float& F2) {
#ifdef count_operation_number
		my_float_auxiliary_storage::operation_number++;
#endif // count_operation_number
		val /= F2.val;
	}

	friend my_float operator + (const my_float& F1, const my_float& F2) {
#ifdef count_operation_number
		my_float_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1.val + F2.val;
	}
	friend my_float operator - (const my_float& F1, const my_float& F2) {
#ifdef count_operation_number
		my_float_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1.val - F2.val;						// result is same as '+' in my_float
	}
	friend my_float operator * (const my_float& F1, const my_float& F2) {
#ifdef count_operation_number
		my_float_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1.val * F2.val;
	}
	friend my_float operator / (const my_float& F1, const my_float& F2) {
#ifdef count_operation_number
		my_float_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1.val / F2.val;
	}

	friend ostream& operator << (ostream& out, my_float c) {
		out << c.val;
		return out;
	}
	friend istream& operator >> (istream& in, my_float& c) {
		in >> c.val;
		return in;
	}
	friend bool operator == (const my_float& F1, const my_float& F2) {
		return F1.val == F2.val;
	}
	friend bool operator !=  (const my_float& F1, const my_float& F2) {
		return F1.val != F2.val;
	}
	friend bool operator >(const my_float& F1, const my_float& F2) {
		return F1.val > F2.val;
	}
	friend bool operator <(const my_float& F1, const my_float& F2) {
		return F1.val < F2.val;
	}
	friend bool operator >=(const my_float& F1, const my_float& F2) {
		return F1.val >= F2.val;
	}
	friend bool operator <=(const my_float& F1, const my_float& F2) {
		return F1.val <= F2.val;
	}
	explicit operator double() const {
		return val;
	}
	explicit operator float() const {
		return val;
	}
	explicit operator int() const {
		return (int)val;
	}
	explicit operator unsigned long long() const {
		return (unsigned long long)val;
	}

	friend my_float sqrt(const my_float& c) {
		return sqrt(c.val);
	}
	friend my_float exp(const my_float& c) {
		return exp(c.val);
	}
	friend my_float log(const my_float& c) {
		return log(c.val);
	}
	friend my_float pow(const my_float& c1, const my_float& c2) {
		return pow(c1.val, c2.val);
	}
	friend my_float log10(const my_float& c) {
		return log10(c.val);
	}
	friend my_float round(const my_float& c) {
		return round(c.val);
	}

	// triangle function is not consider here, you can easily add them
};

#endif // use_my_double
