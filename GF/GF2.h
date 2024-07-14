#pragma once
/*****************************************************************//**
 * \file   GF.h
 * \brief  Galois Field(2), and its operation implementation
 * 
 * \author lilili
 * \date   September 2022
 *********************************************************************/
#include<iostream>

#ifdef use_my_double
#include"../my_lib/my_double.h"
#include"../my_lib/my_float.h"
#endif

using namespace std;
class GF2_auxiliary_storage {
public:
	static unsigned long long operation_number;				// including compare operation number if defined 'count_compare_operation_number'
	static unsigned long long GE_bit_plane_number;
	static double GE_bit_plane_norm_number;					// recording the n-k bit operations as one bit-plane nrom operation
	static unsigned long long re_encoding_bit_plane_norm_number;
	static unsigned long long iteration_number;
};

unsigned long long	GF2_auxiliary_storage::operation_number = 0;
unsigned long long	GF2_auxiliary_storage::GE_bit_plane_number = 0;
double				GF2_auxiliary_storage::GE_bit_plane_norm_number = 0;
unsigned long long	GF2_auxiliary_storage::re_encoding_bit_plane_norm_number = 0;
unsigned long long	GF2_auxiliary_storage::iteration_number = 0;

// write it as a POD mode for speed, but in release mode they are almost same speed (10% of time difference)
class GF2
{
private:
	int a;		// bool is only faster than int about 5%, private is ok
public:
	static void init() {}		// do nothing, for consistance
	/* define construction functions */
	GF2() = default;		// this is POD
	GF2(int b) :a(b) {}	
	GF2 operator - () const {
		return (*this);		// add inverse of a is a in GF(2)
	}
	GF2 flip() const {
#ifdef count_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif // count_operation_number

		return a == 0;
	}

	void operator += (GF2 F2) {
		a = a ^ F2.a;
#ifdef count_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif // count_operation_number
	}
	void operator -= (GF2 F2) {
		a = a ^ F2.a;
#ifdef count_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif // count_operation_number
	}
	void operator *= (GF2 F2) {
		a = a & F2.a;
#ifdef count_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif // count_operation_number
	}
	void operator /= (GF2 F2) {
#ifdef count_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif // count_operation_number
		if (F2.a != 0);
		else {
			throw "error of divide 0";
		}
	}

	friend GF2 operator + (GF2 F1, GF2 F2) {
		F1.a = F1.a ^ F2.a;
#ifdef count_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1;
	}
	friend GF2 operator - (GF2 F1, GF2 F2) {
		F1.a = F1.a ^ F2.a;
#ifdef count_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1;						// result is same as '+' in GF2
	}
	friend GF2 operator * (GF2 F1, GF2 F2) {
		F1.a = F1.a & F2.a;
#ifdef count_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1;
	}
	friend GF2 operator / (GF2 F1, GF2 F2) {
#ifdef count_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif // count_operation_number
		if (F2.a != 0) {
			return F1;
		}
		else {
			throw "error of divide 0";
			return F1;
		}
	}	
	friend GF2 pow(GF2 F1, int n) {		// constrain: n>=0
		return F1;
	}

	friend ostream& operator << (ostream& out, GF2 c) {
		out << c.a;
		return out;
	}
	friend istream& operator >> (istream& in, GF2& c) {
		in >> c.a;
		return in;
	}
	friend bool operator == (GF2 F1, GF2 F2) {
#ifdef count_operation_number
#ifdef count_compare_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif
#endif // count_operation_number

		return F1.a == F2.a;
	}
	friend bool operator !=  (GF2 F1, GF2 F2) {
#ifdef count_operation_number
#ifdef count_compare_operation_number
		GF2_auxiliary_storage::operation_number++;
#endif
#endif // count_operation_number

		return F1.a != F2.a;
	}
	friend bool operator >(GF2 F1, GF2 F2) {
		return true;		// we use true to indicate it has no ordered structure
	}
	friend bool operator <(GF2 F1, GF2 F2) {
		return true;
	}
	friend bool operator >=(GF2 F1, GF2 F2) {
		return true;
	}
	friend bool operator <=(GF2 F1, GF2 F2) {
		return true;
	}
#ifdef use_my_double
	explicit operator my_double() const {
		return (my_double)a;
	}
	explicit operator my_float() const {
		return (my_float)a;
	}
#endif
	explicit operator double() const {
		return (double)a;
	}
	explicit operator int() const {
		return (int)a;
	}
};
