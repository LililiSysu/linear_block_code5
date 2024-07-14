#pragma once
/*****************************************************************//**
 * \file   GFq.h
 * \brief  Galois Field(q)
 * 
 * \author lilili
 * \date   September 2022
 *********************************************************************/
#include"../my_lib/my.h"

#ifdef use_my_double
#include"../my_lib/my_double.h"
#include"../my_lib/my_float.h"
#endif

#include<cmath>
#include<vector>
#include<iostream>
using namespace std;

class GFq_auxiliary_storage {
public:
	static vector<int> mul_inv_tab;
	static unsigned long long operation_number;			
		// cannot declear a static number in templete class, hence put it in GFq_auxiliary_storage

	/* build mul_inv_tab */
	/**
	* .define static member function
	* .build the multiply inverse table for division compute
	*
	*/
	static void build_mul_inv(int q) {
		if (my::is_prime(q)) {
			mul_inv_tab.clear();
			mul_inv_tab.resize(q - 1LL);
			for (int i = 1; i < q; ++i) {
				for (int j = 1; j < q; ++j) {
					if ((j * q + 1) % i != 0);
					else {
						mul_inv_tab[i - 1LL] = (j * q + 1) / i;
						break;
					}
				}
			}
		}
		else {
			throw "error in build_mul_inv, q is not prime";
		}
	}
};
vector<int> GFq_auxiliary_storage::mul_inv_tab(0);
unsigned long long GFq_auxiliary_storage::operation_number = 0;

/**
 * .Galois Field<q>
 */
template<int q> class GF {		// no extension field, and q must be prime number
private:
	int a;
	/* multiply inverse table */
	// but we cannot declear a static varible in template class !!
	//static vector<int> mul_inv;

public:
	/**
	 * .must call this function before use of GF<q>, shell of build_mul_inv()
	 *
	 */
	static void init() {
		GFq_auxiliary_storage::build_mul_inv(q);
	}

	/* define construction functions */
	/**
	 * .initialize 'a' as 0
	 *
	 */
	GF() = default;
	/**
	 * .
	 *
	 * \param _a: value in that field
	 * \return
	 */
	GF(int _a) : a(_a) {}

	/**
	 * .-
	 *
	 * \return (- this_class)
	 */
	GF<q> operator - () const {
#ifdef count_operation_number
		GFq_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return q - a;			// add inverse of a is a in GF
	}	  // add inverse of a

	void operator += (GF<q> F2) {
		a = (a + F2.a) % q;
#ifdef count_operation_number
		GFq_auxiliary_storage::operation_number++;
#endif // count_operation_number
	}
	void operator -= (GF<q> F2) {
		a -= F2.a;
		a += a < 0 ? q : 0;
#ifdef count_operation_number
		GFq_auxiliary_storage::operation_number++;
#endif // count_operation_number
	}
	void operator *= (GF<q> F2) {
		a = (a * F2.a) % q;
#ifdef count_operation_number
		GFq_auxiliary_storage::operation_number++;
#endif // count_operation_number
	}
	void operator /= (GF<q> F2) {
#ifdef count_operation_number
		GFq_auxiliary_storage::operation_number++;
#endif // count_operation_number
		if (F2.a != 0) {
			a = (a * GFq_auxiliary_storage::mul_inv_tab[F2.a - 1]) % q;
		}
		else {
			throw "error of divide 0";
		}
	}

	/**
	 * .+
	 *
	 * \param F1: GF class to be added
	 * \param F2: GF class to add
	 * \return 'F1' plus 'F2'
	 */
	friend GF<q> operator + (GF<q> F1, GF<q> F2) {
		F1.a = (F1.a + F2.a) % q;
#ifdef count_operation_number
		GFq_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1;
	}
	/**
	 * .-
	 *
	 * \param F1: GF class to be minus
	 * \param F2: GF class to minus
	 * \return 'F1' minus 'F2'
	 */
	friend GF<q> operator - (GF<q> F1, GF<q> F2) {
		int a = (F1.a - F2.a);
		F1.a = a < 0 ? a + q : a;
#ifdef count_operation_number
		GFq_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1;
	}
	/**
	 * .*
	 *
	 * \param F1: GF class to be multiplied
	 * \param F2: GF class to multiply
	 * \return 'F1' multiply 'F2'
	 */
	friend GF<q> operator * (GF<q> F1, GF<q> F2) {
		F1.a = (F1.a * F2.a) % q;
#ifdef count_operation_number
		GFq_auxiliary_storage::operation_number++;
#endif // count_operation_number
		return F1;
	}
	/**
	 * ./
	 *
	 * \param F1: GF class to be divided
	 * \param F2: GF class to divide
	 * \return 'F1' divided by 'F2'
	 */
	friend GF<q> operator / (GF<q> F1, GF<q> F2) {
#ifdef count_operation_number
		GFq_auxiliary_storage::operation_number++;
#endif // count_operation_number
		if (F2.a != 0) {
			F1.a = (F1.a * GFq_auxiliary_storage::mul_inv_tab[F2.a - 1]) % q;
			return F1;
		}
		else {
			throw "error of divide 0";
			return 0;
		}
	}

	friend GF<q> pow(GF<q> F1, int n) {		// constrain: n>=0
		n = n % (q - 1);
		F1.a = (int) pow(F1.a, n);
		F1.a %= q;
		return F1;
	}

	/* define operator {<<, >>}*/

	/**
	 * .<<
	 *
	 * \param out: something like cout
	 * \param c: the class to output
	 * \return out
	 */
	friend ostream& operator << (ostream& out, GF<q> c) {
		out << c.a;
		return out;
	}
	/**
	 * .>>
	 *
	 * \param in: something like cin
	 * \param c: the class to input
	 * \return in
	 */
	friend istream& operator >> (istream& in, GF<q>& c) {
		in >> c.a;
		return in;
	}

	/* define operator (==, !=) */

	/**
	 * .==
	 *
	 * \param F: class to compare
	 * \return (this_class == F)
	 */
	friend bool operator == (GF<q> F1, GF<q> F2) {
		return F1.a == F2.a;
	}
	/**
	 * .!=
	 *
	 * \param F: class to compare
	 * \return (this_class != F)
	 */
	friend bool operator != (GF<q> F1, GF<q> F2) {
		return F1.a != F2.a;
	}

	/* define compare operator, always return false since no ordered structure in GF */
	friend bool operator >(GF<q> F1, GF<q> F2) {
		return true;
	}
	friend bool operator <(GF<q> F1, GF<q> F2) {
		return true;
	}
	friend bool operator >=(GF<q> F1, GF<q> F2) {
		return true;
	}
	friend bool operator <=(GF<q> F1, GF<q> F2) {
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
