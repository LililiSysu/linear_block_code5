#pragma once
/*****************************************************************//**
 * \file   GF2e.h
 * \brief  Galois Field(2^m), the extension filed
 * 
 * \author lilili
 * \date   September 2022
 *********************************************************************/
#include"GF2.h"
#include"../my_lib/Matrix.h"
#ifdef use_my_double
#include"../my_lib/my_double.h"
#include"../my_lib/my_float.h"
#endif
using namespace std;

class GF2e_auxiliary_storage;

class polynomial_GF2 {
	// can be replaced by 'polinomial<GF2>' if including "polinomial.h" and "GF2.h"
	// some function is not present as in 'polinomial<GF2>', recommand to use 'polinomial<GF2>'
private:
	/**
	 * .coefficient of a GF2 polynomial, such that coeff {1,0,0,1} indicates 1+0x+0x^2+1x^3
	 */
	Matrix<GF2> coeff;
public:
	/* define construction functions */

	polynomial_GF2() {}
	/**
	 * .this construction function give a way to initialize with length
	 * 
	 * \param len: length of this polynomial
	 * \param type:'0' are usually used to initialize a polynomial
	 */
	polynomial_GF2(int len, char type_0_1_i_b_N) :coeff(Matrix<GF2>(1, len, type_0_1_i_b_N)) {}
	polynomial_GF2(const Matrix<GF2>& _coeff) :coeff(_coeff) {}

	// Access the individual elements
	GF2& operator()(const unsigned& pos_ind) {
		return coeff(pos_ind);
	}
	const GF2& operator()(const unsigned& pos_ind) const {
		return coeff(pos_ind);
	}

	/**
	 * .format the coeff to certain length, that is to append 0s
	 * 
	 * \param len: the length
	 */
	void format_len(int len = 1) {
		// extend to num coefficents with 0s
		int len_odd = simplify() + 1;
		if (len_odd >= len)  return;
		else {
			Matrix<GF2> tmp(1, len - len_odd, '0');
			coeff = coeff.combine_right(tmp);
			return;
		}
	}
	/**
	 * .set coeff.c = index of end non 0
	 * 
	 * \return the index of the last 1, for example, {0,0,1,1,0} will return 3
	 */
	inline int simplify() {
		return coeff.cut_end_0();
	}
	/**
	 * .
	 * 
	 * \return the index of the last 1, for example, {0,0,1,1,0} will return 3
	 */
	inline int ind_of_end_non_0() const {
		return coeff.ind_of_end_non_0();
	}
	/**
	 * .
	 * 
	 * \param len_ind: coeff with index less or equal than 'len_ind' will be fetched
	 * \return the former part of coeff
	 */
	inline polynomial_GF2 get_former(int len_ind) const {
		// get the former "num" coefficient		
		return polynomial_GF2(coeff.get_part(0, 0, 0, len_ind));
	}
	/**
	 * . append 0 to the left of coeff
	 * 
	 * \param num: the number of 0s to append
	 * \return the result 'polynomial_GF2'
	 */
	void shift_right(int num) {
		coeff.shift_right(num);
	}
	/**
	 * .
	 * 
	 * \param num
	 * \return 
	 */
	polynomial_GF2 get_shift_right(int num) const{
		return coeff.get_shift_right(num);
	}
	inline int size() const {
		return coeff.size();
	}
	inline Matrix<GF2> get_coeff() const {
		return coeff;
	}

	/* define operator {=, +, -, *, /} */
	polynomial_GF2& operator = (const polynomial_GF2& p) {
		if (this != &p) {	// judge if set class for itself
			coeff = p.coeff;
		}
		return *this;
	}
	polynomial_GF2& operator = (const Matrix<GF2>& _coeff) {
		coeff = _coeff;
		return *this;
	}
	polynomial_GF2 operator - () const {
		return -coeff;
	}		// add inverse of coeff is itself

	void operator += (const polynomial_GF2& p2) {
		int len1 = size();
		int len2 = p2.size();
		if (len1 >= len2) {
			for (int i = 0; i < len2; ++i) {
				coeff(i) += p2.coeff(i);
			}
		}
		else {
			coeff.resize(1, len2, true);
			for (int i = 0; i < len1; ++i) {
				coeff(i) += p2.coeff(i);
			}
			for (int j = len1; j < len2; ++j) {
				coeff(j) = p2.coeff(j);
			}
		}
	}
	void operator -= (const polynomial_GF2& p2) {
		int len1 = size();
		int len2 = p2.size();
		if (len1 >= len2) {
			for (int i = 0; i < len2; ++i) {
				coeff(i) -= p2.coeff(i);
			}
		}
		else {
			coeff.resize(1, len2, true);
			for (int i = 0; i < len1; ++i) {
				coeff(i) -= p2.coeff(i);
			}
			for (int j = len1; j < len2; ++j) {
				coeff(j) = -p2.coeff(j);
			}
		}
	}
	void operator *= (const polynomial_GF2& p2) {			// this cannot be written as in-place operation
		(*this) = (*this) * p2;
	}	
	void operator /= (const polynomial_GF2& p2) {
		// return the quotient only, p1/p2
		//polynomial_GF2 result;
		//Matrix<GF2> coeff1 = p1.coeff;
		Matrix<GF2> coeff2 = p2.coeff;
		simplify();
		coeff2.cut_end_0();
		if (coeff2.size() == 1 && coeff2(0) == 0)
			throw "divide by 0";
		int len1 = coeff.size();
		int len2 = coeff2.size();

		if (len1 >= len2) {
			Matrix<GF2> result_coeff(1, len1 - len2 + 1, '0');

			while (len1 >= len2) {
				result_coeff(len1 - len2) = 1;
				for (int k = 0; k < len2; ++k) {
					coeff(len1 - 1 - k) = coeff(len1 - 1 - k) - coeff2(len2 - 1 - k);
				}
				len1 = coeff.ind_of_end_non_0() + 1;
			}
			simplify();
			//coeff2.simplify();
			if (coeff == coeff2) {
				result_coeff(0) = 1;
			}
			coeff = result_coeff;
		}
		else {
			coeff = Matrix<GF2>(1, 1, '0');
		}
	}
	void operator %= (const polynomial_GF2& p2) {
		// return the remainder only, p1%p2
		//polynomial_GF2 result;
		Matrix<GF2> coeff2 = p2.coeff;
		simplify();
		coeff2.cut_end_0();		// force the polynomial be no end 0s is important
		if (coeff2.size() == 1 && coeff2(0) == 0)
			throw "divide by 0";
		int len1 = coeff.size();
		int len2 = coeff2.size();

		while (len1 >= len2) {
			for (int k = 0; k < len2; ++k) {
				coeff(len1 - 1 - k) = coeff(len1 - 1 - k) - coeff2(len2 - 1 - k);
			}
			len1 = coeff.ind_of_end_non_0() + 1;
		}

		simplify();
	}

	friend polynomial_GF2 operator + (const polynomial_GF2& p1, const polynomial_GF2& p2) {

		int len1 = p1.size();
		int len2 = p2.size();
		if (len1 == len2) {
			return p1.coeff + p2.coeff;
		}
		else if (len1 > len2) {
			polynomial_GF2 result(len1, '0');
			for (int i = 0; i < len2; ++i) {
				result(i) = p1.coeff(i) + p2.coeff(i);
			}
			for (int j = len2; j < len1; ++j) {
				result(j) = p1.coeff(j);
			}
			return result;
		}
		else {
			polynomial_GF2 result(len2, '0');
			for (int i = 0; i < len1; ++i) {
				result(i) = p1.coeff(i) + p2.coeff(i);
			}
			for (int j = len1; j < len2; ++j) {
				result(j) = p1.coeff(j);
			}
			return result;
		}
	}
	friend polynomial_GF2 operator - (const polynomial_GF2& p1, const polynomial_GF2& p2) {
		int len1 = p1.size();
		int len2 = p2.size();
		if (len1 == len2) {
			return p1.coeff - p2.coeff;
		}
		else if (len1 > len2) {
			polynomial_GF2 result(len1, '0');
			for (int i = 0; i < len2; ++i) {
				result(i) = p1.coeff(i) - p2.coeff(i);
			}
			for (int j = len2; j < len1; ++j) {
				result(j) = p1.coeff(j);
			}
			return result;
		}
		else {
			polynomial_GF2 result(len2, '0');
			for (int i = 0; i < len1; ++i) {
				result(i) = p1.coeff(i) - p2.coeff(i);
			}
			for (int j = len1; j < len2; ++j) {
				result(j) = -p2.coeff(j);
			}
			return result;
		}
	}
	friend polynomial_GF2 operator * (const polynomial_GF2& p1, const polynomial_GF2& p2) {
		int len1 = p1.size();
		int len2 = p2.size();
		polynomial_GF2 result(Matrix<GF2>(1, len1 + len2 - 1, '0'));
		for (int j = 0; j < len2; ++j) {
			if (p2.coeff(j) == 0);
			else {
				for (int i = 0; i < len1; ++i) {
					result.coeff(i + j)+= p1.coeff(i) * p2.coeff(j);
				}
			}
		}
		result.simplify();
		return result;
	}
	friend polynomial_GF2 operator / (const polynomial_GF2& p1, const polynomial_GF2& p2) {
		// return the quotient only, p1/p2
		//polynomial_GF2 result;
		Matrix<GF2> coeff1 = p1.coeff;
		Matrix<GF2> coeff2 = p2.coeff;
		coeff1.cut_end_0();
		coeff2.cut_end_0();
		if (coeff2.size() == 1 && coeff2(0) == 0)
			throw "divide by 0";
		int len1 = coeff1.size();
		int len2 = coeff2.size();

		if (len1 >= len2) {
			Matrix<GF2> result_coeff(1, len1 - len2 + 1, '0');

			while (len1 >= len2) {
				result_coeff(len1 - len2)= 1;
				for (int k = 0; k < len2; ++k) {
					coeff1(len1 - 1 - k) -= coeff2(len2 - 1 - k);
				}
				len1 = coeff1.ind_of_end_non_0() + 1;
			}
			coeff1.cut_end_0();
			//coeff2.simplify();
			if (coeff1 == coeff2) {
				result_coeff(0)= 1;
			}
			return result_coeff;
			//result.coeff = result_coeff;
			//return result;
		}
		else {
			return Matrix<GF2>(1, 1, '0');
			//return result;
		}
	}
	friend polynomial_GF2 operator % (const polynomial_GF2& p1, const polynomial_GF2& p2) {
		// return the remainder only, p1%p2
		//polynomial_GF2 result;
		Matrix<GF2> coeff1 = p1.coeff;
		Matrix<GF2> coeff2 = p2.coeff;
		coeff1.cut_end_0();
		coeff2.cut_end_0();
		if (coeff2.size() == 1 && coeff2(0) == 0)
			throw "divide by 0";
		int len1 = coeff1.size();
		int len2 = coeff2.size();

		while (len1 >= len2) {
			for (int k = 0; k < len2; ++k) {
				coeff1(len1 - 1 - k) -= coeff2(len2 - 1 - k);
			}
			len1 = coeff1.ind_of_end_non_0() + 1;
		}

		coeff1.cut_end_0();

		//result.coeff = coeff1;
		return coeff1;
	}

	/* define operator {<<,>>} */
	friend ostream& operator << (ostream& out, const polynomial_GF2& p) {
		out << p.coeff;
		return out;
	}
	friend istream& operator >> (istream& in, polynomial_GF2& p) {
		in >> p.coeff;
		return in;
	}

	/* define operator {==,!=} */
	friend bool operator == (const polynomial_GF2& p1, const polynomial_GF2& p2) {
		return p1.coeff == p2.coeff;
	}
	friend bool operator != (const polynomial_GF2& p1, const polynomial_GF2& p2) {
		return p1.coeff != p2.coeff;
	}

	/* class to store calculation tables for GF2e */
	friend GF2e_auxiliary_storage;
};

class GF2e_auxiliary_storage {
private:
	// this function should not be visit by anyone
	static void generate_polynomial_table_and_alpha_table(const vector<polynomial_GF2>& addition_table) {	
		// the two table are only valid for GF2e, pay attention

		// using the result of addition table
		unsigned len = (unsigned)addition_table.size();
		int poly_len = addition_table[0].coeff.size();
		polynomial_table.resize(len, 0);
		alpha_table.resize(len, 0);
		for (unsigned i = 0; i < len; ++i) {
			int poly_index = 0;
			for (int j = 0; j < poly_len; ++j) {
				poly_index <<= 1;
				poly_index = poly_index + (addition_table[i].coeff(j) == 1);
			}
			polynomial_table[i] = poly_index;
			alpha_table[poly_index] = i;
		}
	}

public:

	// primitive polynomial's coefficient, starting from x^0 to x^m
	static vector<vector<GF2>> p_polynomials;

	// GF table for switching between alpha^n to a polynomial
	static vector<int> polynomial_table;	// like 'addition_table' but store polynomial index

	// GF table for switching between a polynomail to alpha^n
	static vector<int> alpha_table;

	// number of operations in GF2e
	static unsigned long long operation_number;

	static unsigned long long add_number;

	static unsigned long long mul_number;

	/**
	 * .generate 'polynomial_table' and 'alpha_table' for addition calculation in GF2e
	 *
	 * \param m: to indicate GF(2^m)
	 */
	static void generate_compute_table(int m) {			// only for GF2 extension field

		// GF table for switching between alpha^n to a polynomial
		vector<polynomial_GF2> addition_table;		// automatically let this destory

		// 2^m polynomials corresponding to 0, 1, alpha^1 to alpha^{2^m-2}
		polynomial_GF2 equ(Matrix<GF2>(1, m + 1, p_polynomials[m - 2]));
		int len2 = equ.ind_of_end_non_0();
		polynomial_GF2 equ_left(equ.get_former(len2 - 1));
		int order = 1 << m;
		addition_table.push_back(polynomial_GF2(m, '0'));	// for all 0 element
		for (int i = 1; i < order; ++i) {
			if (i <= m) {
				addition_table.push_back(polynomial_GF2(m, '0'));
				addition_table.back()(i - 1) = 1;
			}
			else {
				polynomial_GF2 tmp = addition_table.back().get_shift_right(1);
				if (tmp.ind_of_end_non_0() != len2);
				else {
					tmp = tmp.get_former(len2 - 1) + equ_left; // this is only valid for GF2 extension field
				}
				tmp.format_len(m);
				addition_table.push_back(tmp);
			}
		}
		generate_polynomial_table_and_alpha_table(addition_table);
		q = (1 << m);

		// for new GF2e
		vector<int> tmp_store_alpha_table = alpha_table;
		for (int i = 1, imax = (int)alpha_table.size(); i < imax; ++i) {
			polynomial_table[i - 1] = my::rev_bits(polynomial_table[i], m);		// of size q-1;
			alpha_table[my::rev_bits(i, m)] = tmp_store_alpha_table[i] - 1;
		}
		polynomial_table.pop_back();
		alpha_table[0] = -1;		// just a dummy element, please do not visit alpha_table[0], be careful!

		q_minus_1 = q - 1;
	}

	/* size of GF2e */
	static int q;

	/* a new variable for GF2e multiplication */
	static int q_minus_1;

#ifdef need_trace_table
	static vector<GF2> trace_table;
	static vector<int> GF2e_that_trace_to_0;
	static vector<int> GF2e_that_trace_to_1;		// this should be removed, saving space
#endif
};

unsigned long long GF2e_auxiliary_storage::operation_number = 0;
unsigned long long GF2e_auxiliary_storage::add_number = 0;
unsigned long long GF2e_auxiliary_storage::mul_number = 0;

vector<vector<GF2>> GF2e_auxiliary_storage::p_polynomials = vector<vector<GF2>>(
	{					// starting from m=2 to m=24
		{1,1,1},											// 2
		{1,1,0,1},											// 3
		{1,1,0,0,1},										// 4
		{1,0,1,0,0,1},										// 5
		{1,1,0,0,0,0,1},									// 6
		{1,0,0,1,0,0,0,1},									// 7
		{1,0,1,1,1,0,0,0,1},								// 8
		{1,0,0,0,1,0,0,0,0,1},								// 9
		{1,0,0,1,0,0,0,0,0,0,1},							// 10
		{1,0,1,0,0,0,0,0,0,0,0,1},							// 11
		{1,1,0,0,1,0,1,0,0,0,0,0,1},						// 12
		{1,1,0,1,1,0,0,0,0,0,0,0,0,1},						// 13
		{1,1,0,0,0,0,1,0,0,0,1,0,0,0,1},					// 14
		{1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},					// 15
		{1,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1},				// 16
		{1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},				// 17
		{1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1},			// 18
		{1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1},			// 19
		{1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},		// 20
		{1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},		// 21
		{1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},	// 22
		{1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},	// 23
		{1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}	// 24
	}
);
vector<int> GF2e_auxiliary_storage::polynomial_table(0);
vector<int> GF2e_auxiliary_storage::alpha_table(0);
int GF2e_auxiliary_storage::q = 0;
int GF2e_auxiliary_storage::q_minus_1 = 0;

#ifdef need_trace_table
vector<GF2> GF2e_auxiliary_storage::trace_table(0);
vector<int> GF2e_auxiliary_storage::GF2e_that_trace_to_0(0);
vector<int> GF2e_auxiliary_storage::GF2e_that_trace_to_1(0);
#endif

template<int m> class GF2e {		// extension field GF(2^m), and m must be greater than 0
private:
	int a;	// ranging from 0 to (1<<m)-1
	/* table for conversion between 'int' to 'polynomial_GF2' */
	/**
	 * .but we cannot declear a static varible in template class !! Now, we can!
	 * use GF2e_auxiliary_storage::addition_table and
	 * GF2e_auxiliary_storage::polynomial_table and GF2e_auxiliary_storage::alpha_table instead
	 */

	// trace function to get an GF2
	GF2 trace_inside() const {
		GF2e<m> ans(*this);
		for (int i = 1; i < m; ++i) {
			ans += pow((*this), 1 << i);
		}
		GF2 ret = (int)ans;
		return ret;
	}

public:

	/**
	 * .must call this function before any use of GF2e
	 *
	 */
	static void init() {
		GF2e_auxiliary_storage::generate_compute_table(m);

		// initialize trace table

#ifdef need_trace_table
		GF2e<m> tmp;
		GF2e_auxiliary_storage::trace_table.clear();
		GF2e_auxiliary_storage::GF2e_that_trace_to_0.clear();
		GF2e_auxiliary_storage::GF2e_that_trace_to_1.clear();
		for (int i = 0; i < GF2e_auxiliary_storage::q; ++i) {
			tmp = i;
			GF2 trace_result = tmp.trace_inside();
			GF2e_auxiliary_storage::trace_table.push_back(trace_result);
			if (trace_result == 0) {
				GF2e_auxiliary_storage::GF2e_that_trace_to_0.push_back(i);
			}
			else{
				GF2e_auxiliary_storage::GF2e_that_trace_to_1.push_back(i);
			}
		}
#endif
	}

	/* define construction functions */
	/**
	 * .initialize 'a' as 0, will not judge if generate_addition_table(m) has done, for speed
	 *
	 */
	GF2e() = default;
	/**
	 * .will not judge if integer is in field
	 *
	 * \param _a: value in that field, the bits indicates polynomial coefficients in GF2
	 */
	GF2e(int _a) :a(_a) {};
	/**
	 * . initialize with sub field
	 */
	GF2e(GF2 _a) :a(_a) {};

	// we omitted the = operator, using the function generated by the compiler by default

	// set by alpha^_a
	void set_by_alpha_power(int _a) {
		// if _a can be any value, please reprocess _a with
		//	_a %= GF2e_auxiliary_storage::q_minus_1;

		// _a in {-1} U [0,q-2] representing alpha^_a
		if (_a >= 0) {
			a = GF2e_auxiliary_storage::polynomial_table[_a];											// alpha^0 = 1
		}
		else {
			a = GF2e_auxiliary_storage::polynomial_table[GF2e_auxiliary_storage::q_minus_1 + _a];		// alpha^-1 = alpha^(q-2)
		}
	}

	// trace function to get an GF2
	inline GF2 trace() const {
#ifdef need_trace_table
		return GF2e_auxiliary_storage::trace_table[a];
#else
		return trace_inside();
#endif		
	}

	inline int alpha_power() const {
		return GF2e_auxiliary_storage::alpha_table[a];
	}

	/**
	 * .-
	 *
	 * \return (- this_class)
	 */
	GF2e<m> operator - () const {
		return (*this);
	}	  // add inverse of a

	void operator += (const GF2e<m>& F2) {
		// the following code are only valid for GF2e
		a ^= F2.a;

#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::add_number++;
#endif // count_add_mul_number
	}
	void operator -= (const GF2e<m>& F2) {
		// the following code are only valid for GF2e
		a ^= F2.a;

#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::add_number++;
#endif // count_add_mul_number
	}
	void operator *= (const GF2e<m>& F2) {
		if (a != 0 && F2.a != 0) {

			int tmp = GF2e_auxiliary_storage::alpha_table[a] + GF2e_auxiliary_storage::alpha_table[F2.a];
			tmp -= tmp >= GF2e_auxiliary_storage::q_minus_1 ? GF2e_auxiliary_storage::q_minus_1 : 0;
			a = GF2e_auxiliary_storage::polynomial_table[tmp];
		}
		else {
			a = 0;
		}

#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::mul_number++;
#endif // count_add_mul_number
	}
	void operator /= (const GF2e<m>& F2) {
#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::mul_number++;
#endif // count_add_mul_number
		if (a != 0 && F2.a != 0) {
			int tmp = GF2e_auxiliary_storage::alpha_table[a] - GF2e_auxiliary_storage::alpha_table[F2.a];
			tmp += tmp < 0 ? GF2e_auxiliary_storage::q_minus_1 : 0;
			a = GF2e_auxiliary_storage::polynomial_table[tmp];
		}
		else if (F2.a == 0) {
			throw "error of divide 0";
		}
	}

	/**
	 * .+
	 *
	 * \param F1: GF2e<m> class to be added
	 * \param F2: GF2e<m> class to add
	 * \return 'F1' plus 'F2'
	 */
	friend GF2e<m> operator + (const GF2e<m>& F1, const GF2e<m>& F2) {
		// the following code are only valid for GF2e

#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::add_number++;
#endif // count_add_mul_number
		return F1.a ^ F2.a;
	}
	/**
	 * .-
	 *
	 * \param F1: GF2e<m> class to be added
	 * \param F2: GF2e<m> class to add
	 * \return 'F1' plus 'F2'
	 */
	friend GF2e<m> operator - (const GF2e<m>& F1, const GF2e<m>& F2) {
		// the following code are only valid for GF2e

#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::add_number++;
#endif // count_add_mul_number
		return F1.a ^ F2.a;
	}
	/**
	 * .*
	 *
	 * \param F1: GF2e<m> class to be multiplied
	 * \param F2: GF2e<m> class to multiply
	 * \return 'F1' multiply 'F2'
	 */
	friend GF2e<m> operator * (const GF2e<m>& F1, const GF2e<m>& F2) {

#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::mul_number++;
#endif // count_add_mul_number

		if (F1.a != 0 && F2.a != 0) {
			int tmp = GF2e_auxiliary_storage::alpha_table[F1.a] + GF2e_auxiliary_storage::alpha_table[F2.a];
			tmp -= tmp >= GF2e_auxiliary_storage::q_minus_1 ? GF2e_auxiliary_storage::q_minus_1 : 0;
			return GF2e_auxiliary_storage::polynomial_table[tmp];
		}
		else {
			return 0;
		}
	}
	/**
	 * ./
	 *
	 * \param F1: GF2e<m> class to be divided
	 * \param F2: GF2e<m> class to divide
	 * \return 'F1' divided by 'F2'
	 */
	friend GF2e<m> operator / (const GF2e<m>& F1, const GF2e<m>& F2) {

#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::mul_number++;
#endif // count_add_mul_number

		if (F1.a != 0 && F2.a != 0) {
			int tmp = GF2e_auxiliary_storage::alpha_table[F1.a] - GF2e_auxiliary_storage::alpha_table[F2.a];
			tmp += tmp < 0 ? GF2e_auxiliary_storage::q_minus_1 : 0;
			return GF2e_auxiliary_storage::polynomial_table[tmp];
		}
		else if (F2.a == 0) {
			throw "error of divide 0";
		}
		return 0;
	}

	/**
	 * .*
	 *
	 * \param n: number of F2 add together
	 * \param F2: GF2e<m> class to accumulate
	 * \return F2+F2+...+F2 (n of them)
	 */
	friend GF2e<m> operator * (const int& n, const GF2e<m>& F2) {
		// n F2 add together
#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::add_number++;
#endif // count_add_mul_number
		return (n & 1) ? F2 : 0;
	}

	friend GF2e<m> pow(const GF2e<m>& F1, const int& n) {				// constrain: n>=0

#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::mul_number++;
#endif // count_add_mul_number

		if (F1.a != 0) {
			return GF2e_auxiliary_storage::polynomial_table[\
				(GF2e_auxiliary_storage::alpha_table[F1.a] * n) % GF2e_auxiliary_storage::q_minus_1];
		}
		else {
			return 0;
		}
	}

	/* define operator {<<, >>}*/

	/**
	 * .<<
	 *
	 * \param out: something like cout
	 * \param c: the class to output
	 * \return out
	 */
	friend ostream& operator << (ostream& out, const GF2e<m>& F1) {
		if (F1.a != 0) {
			out << GF2e_auxiliary_storage::alpha_table[F1.a];
		}
		else {
			out << '.';			// this symbol represents 0 element of finite field
		}
		return out;
	}
	/**
	 * .>>
	 *
	 * \param in: something like cin
	 * \param c: the class to input
	 * \return in
	 */
	friend istream& operator >> (istream& in, GF2e<m>& F1) {
		// not finished
		int tmp;
		in >> tmp;
		if (tmp == -1) {
			// use -1 to indacate 0 in GF2e<m>, as 0 represents alpha^0 = 1
			F1.a = 0;
		}
		else {
			F1.set_by_alpha_power(tmp);
		}

		return in;
	}

	/* define operator (==, !=) */

	/**
	 * .==
	 *
	 * \param F: class to compare
	 * \return (this_class == F)
	 */
	friend bool operator == (const GF2e<m>& F1, const GF2e<m>& F2) {

#ifdef count_compare_operation_number
#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::add_number++;
#endif // count_add_mul_number
#endif // count_compare_operation_number

		return F1.a == F2.a;
	}
	/**
	 * .!=
	 *
	 * \param F: class to compare
	 * \return (this_class != F)
	 */
	friend bool operator != (const GF2e<m>& F1, const GF2e<m>& F2) {

#ifdef count_compare_operation_number
#ifdef count_operation_number
		GF2e_auxiliary_storage::operation_number++;
#endif // count_operation_number

#ifdef count_add_mul_number
		GF2e_auxiliary_storage::add_number++;
#endif // count_add_mul_number
#endif // count_compare_operation_number

		return F1.a != F2.a;
	}

	/* define compare operator, always return false since no ordered structure in GF */
	friend bool operator >(const GF2e<m>& F1, const GF2e<m>& F2) {
		return true;
	}
	friend bool operator <(const GF2e<m>& F1, const GF2e<m>& F2) {
		return true;
	}
	friend bool operator >=(const GF2e<m>& F1, const GF2e<m>& F2) {
		return true;
	}
	friend bool operator <=(const GF2e<m>& F1, const GF2e<m>& F2) {
		return true;
	}
#ifdef use_my_double
	explicit operator my_double() const {
		return (my_double)a;
	}
	explicit operator my_float() const {
		return (my_float)a;
	}
	explicit operator GF2() const {
		return a != 0;
	}
#endif
	explicit operator double() const {
		return (double)a;
	}
	explicit operator int() const {
		return (int)a;
	}

	/**
	 * .turn a Matrix on field T to GF2, where each element is replaces by a column vector on GF2
	 */
	static Matrix<GF2> to_bits_row_extention(const Matrix<GF2e<m>>& M) {
		int r = M.row();
		int c = M.col();
		Matrix<GF2> result(r * m, c);
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
				int sym = (int)M(i, j);
				//cout << "sym=" << sym << endl;
				int result_r_ind = i * m;
				for (int k = 0; k < m; ++k) {
					result(result_r_ind + k, j) = sym & 1;
					sym >>= 1;
				}
			}
		}
		return result;
	}
};

