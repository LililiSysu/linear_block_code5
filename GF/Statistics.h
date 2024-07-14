/*****************************************************************//**
 * \file   Statistics.h
 * \brief  count the operation number and other statistics for the program
 * 
 * \author 26259
 * \date   August 2023
 *********************************************************************/

#pragma once

#include<iostream>
#include"GF2.h"
#include"GF2e.h"
using namespace std;

#ifdef use_my_double
#include"../my_lib/my_double.h"
#else
#define my_double double
#endif  // use_my_double

class ope_count{
public:
#ifdef use_my_double
	unsigned long long double_ope;
#endif // use_my_double
	unsigned long long GF2_ope;
	unsigned long long GF2e_ope;

	ope_count() = default;

	// output the ope_count
	friend ostream& operator << (ostream& out, const ope_count& a) {
#ifdef use_my_double
		out << "[double]: " << a.double_ope << "\t";
#endif // use_my_double
		out << "[GF2]: " << a.GF2_ope << "\t"
			<< "[GF2e]: " << a.GF2e_ope << "\n";

		return out;
	}

	void count() {
#ifdef use_my_double
		double_ope = my_double_auxiliary_storage::operation_number;
#endif // use_my_double
		GF2_ope = GF2_auxiliary_storage::operation_number;
		GF2e_ope = GF2e_auxiliary_storage::operation_number;
	}

	friend ope_count operator - (const ope_count& a, const ope_count& b) {
		ope_count ans;
#ifdef use_my_double
		ans.double_ope = a.double_ope - b.double_ope;
#endif // use_my_double
		ans.GF2_ope = a.GF2_ope - b.GF2_ope;
		ans.GF2e_ope = a.GF2e_ope - b.GF2e_ope;
		return ans;
	}
};
