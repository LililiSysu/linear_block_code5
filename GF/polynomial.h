#pragma once
/*****************************************************************//**
 * \file   polynomial.h
 * \brief  polynomial class for any type of coefficients
 * 
 * \author lilili
 * \date   September 2022
 *********************************************************************/
#include"../my_lib/Matrix.h"
using namespace std;

template<class T>
class polynomial {
private:
	/**
	 * .coefficient of a T polynomial, such that coeff {1,0,0,1} indicates 1+0x+0x^2+1x^3
	 * also valid for GF2e element such that coeff {a^3, a^1 ,0 ,1} indicates a^3 + a^1 x + 0 x^2 + x^3
	 */
	Matrix<T> coeff;		// we prevent storing ending zeors in coeff
public:
	/* define construction functions */

	// to compose it into own written class Matrix, we have to define a one parameter constructor
	polynomial(const T& element = 0) : coeff(T(element)) {}; // it has at least one coefficient, representing the polynomial is zero
	/**
	 * .this construction function give a way to initialize with length
	 *
	 * \param len: length of this polynomial
	 * \param type:'0' are usually used to initialize a polynomial
	 */
	polynomial(int len, char type_0_1_i_b_N) :coeff(Matrix<T>(1, len, type_0_1_i_b_N)) {}
	polynomial(const Matrix<T>& _coeff) :coeff(_coeff) {}

	// Access the individual elements
	T& operator()(const unsigned& pos_ind) {
		return coeff(pos_ind);
	}
	const T& operator()(const unsigned& pos_ind) const {
		return coeff(pos_ind);
	}

	/**
	 * .format the coeff to certain length, that is to append 0s, do not use it during polynomial computation
	 *
	 * \param len: the length
	 */
	void format_len(int len = 1) {
		// extend to num coefficents with 0s
		int len_old = simplify() + 1;
		if (len_old >= len)  return;
		else {
			coeff.resize(1, len, true);
			coeff.set_part(0, len_old, 0, len - 1, T(0));
			return;
		}
	}
	void rev() {		// this should be fixed using inherit property of C++
		coeff.rev();
	}

	/**
	 * .set coeff.c = index of end non 0
	 *
	 * \return the index of the last 1, for example, {0,0,1,1,0} will return 3
	 */
	inline int simplify() {
		if (T() < T())					// if it is not order structure
			return coeff.cut_end_0();
		else
			return coeff.cut_end_small();
	}
	/**
	 * .
	 *
	 * \return the index of the last 0, for example, {0,0,1,1,0} will return 3, {0,0,0,0} will return 0
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
	inline polynomial<T> get_former(int len_ind) const {
		// get the former "num" coefficient		
		return polynomial<T>(coeff.get_part_ind(0, 0, 0, len_ind));
	}
	inline Matrix<T> get_coeff() const {
		return coeff;
	}
	inline int size() const {
		return coeff.size();
	}
	inline const T* ptr() const {
		return coeff.ptr();
	}
	/**
	 * . append 0 to the left of coeff
	 *
	 * \param num: the number of 0s to append
	 * \return the result 'polynomial'
	 */
	void shift_right(int num) {
		coeff.shift_right(num);
	}
	polynomial<T> get_shift_right(int num) const {
		return coeff.get_shift_right(num);
	}
	/**
	 * .
	 * 
	 * \return the derivative of polynomial
	 */
	polynomial<T> get_derivative() const {
		int len = coeff.ind_of_end_non_0();
		if (len != 0) {
			polynomial<T> result(len, '0');	// be careful of empty polynomial
			for (int i = 0; i < len; ++i) {
				//T coeff_element = coeff(i + 1);
				//T coeff_element_result = coeff_element;
				//for (int j = 0; j < i; ++j) {				// write w * a= a+a+...+a (w a add together) for computation on finite field
				//	coeff_element_result += coeff_element;
				//}
				//result(i)= coeff_element_result;
				result(i) = (i + 1) * coeff(i + 1);
			}
			return result;
		}
		else {
			return polynomial<T>(1, '0');
		}
	}

	/* define operator {=, +, -, *, /} */
	//polynomial<T>& operator = (const polynomial<T>& p) {
	//	if (this != &p) {	// judge if set class for itself
	//		coeff = p.coeff;
	//	}
	//	return *this;
	//}
	//polynomial<T>& operator = (const vector<T>& _coeff){
	//	coeff = Matrix<T>(_coeff);
	//	return *this;
	//}
	//polynomial<T>& operator = (const Matrix<T>& _coeff){
	//	coeff = _coeff;
	//	return *this;
	//}
	polynomial<T> operator - () const {
		return -coeff;
	}		// add inverse of coeff is itself

	// we assume that the polynomials, during computation, is simplified
	void operator += (const polynomial<T>& p2) {
		int len1 = size();
		int len2 = p2.size();
		if (len1 >= len2) {
			for (int i = 0; i < len2; ++i) {
				coeff(i) += p2(i);
			}
		}
		else {
			coeff.resize(1, len2, true);
			for (int i = 0; i < len1; ++i) {
				coeff(i) += p2(i);
			}
			for (int j = len1; j < len2; ++j) {
				coeff(j) = p2(j);
			}
		}
		simplify();
	}
	void operator -= (const polynomial<T>& p2) {
		int len1 = size();
		int len2 = p2.size();
		if (len1 >= len2) {
			for (int i = 0; i < len2; ++i) {
				coeff(i) -= p2(i);
			}
		}
		else {
			coeff.resize(1, len2, true);
			for (int i = 0; i < len1; ++i) {
				coeff(i) -= p2(i);
			}
			for (int j = len1; j < len2; ++j) {
				coeff(j) = -p2(j);
			}
		}
		simplify();
	}
	void operator *= (const polynomial<T>& p2) {// this cannot be written as in-place operation
		int len1 = size();
		int len2 = p2.size();
		coeff.resize(1, len1 + len2 - 1, true);

		T store;
		for (int k = len1 + len2 - 2; k >= 0; --k) {
			// k is the sum of processing degree of (*this) and p2
			int lower_bound = my::max(0, k - len2 + 1);
			store = 0;
			for (int i = my::min(len1 - 1, k); i >= lower_bound; --i) {
				// i is the degree of proceessing (*this)
				store += coeff(i) * p2(k - i);
			}
			coeff(k) = store;
		}
		simplify();
	}
	void operator *= (const T& num) {
		if (num != T(0)) {
			int len = size();
			for (int i = 0; i < len; ++i) {
				coeff(i) *= num;
			}
		}
		else {
			coeff = T(0);
		}
		
	}
	void operator /= (const polynomial<T>& p2) {
		// return the quotient only, p1/p2
		(*this) = (*this) / p2;
	}
	void operator /= (const T& num) {
		int len = size();
		T one_over_num = 1 / num;
		for (int i = 0; i < len; ++i) {
			coeff(i) *= one_over_num;
		}
	}
	void operator %= (const polynomial<T>& p2) {
		// return the remainder only, p1%p2
		int len2 = p2.size();
		if (len2 == 1 && p2(0) == 0)
			throw "divide by 0";

		int len1 = size();

		while (len1 >= len2 && !is_approach_0()) {
			T num = (*this)(len1 - 1) / p2(len2 - 1);
			for (int k = 0; k < len2; ++k) {
				(*this)(len1 - 1 - k) -= num * p2(len2 - 1 - k);
			}
			simplify();
			len1 = size();
		}
	}

	friend polynomial<T> operator + (const polynomial<T>& p1, const polynomial<T>& p2) {		
		int len1 = p1.size();
		int len2 = p2.size();
		if (len1 >= len2) {
			polynomial<T> result(len1, '0');
			for (int i = 0; i < len2; ++i) {
				result(i)= p1(i) + p2(i);
			}
			for (int j = len2; j < len1; ++j) {
				result(j)= p1(j);
			}
			result.simplify();
			return result;
		}
		else {
			polynomial<T> result(len2, '0');
			for (int i = 0; i < len1; ++i) {
				result(i)= p1(i) + p2(i);
			}
			for (int j = len1; j < len2; ++j) {
				result(j)= p1(j);
			}
			result.simplify();
			return result;
		}
	}
	friend polynomial<T> operator - (const polynomial<T>& p1, const polynomial<T>& p2) {
		int len1 = p1.size();
		int len2 = p2.size();
		if (len1 >= len2) {
			polynomial<T> result(len1, '0');
			for (int i = 0; i < len2; ++i) {
				result(i) = p1(i) - p2(i);
			}
			for (int j = len2; j < len1; ++j) {
				result(j) = p1(j);
			}
			result.simplify();
			return result;
		}
		else {
			polynomial<T> result(len2, '0');
			for (int i = 0; i < len1; ++i) {
				result(i) = p1(i) - p2(i);
			}
			for (int j = len1; j < len2; ++j) {
				result(j) = -p2(j);		// note that for non-GF2(GF2e), - symbol is necessary
			}
			result.simplify();
			return result;
		}
	}
	friend polynomial<T> operator * (const polynomial<T>& p1, const polynomial<T>& p2) {
		int len1 = p1.size();
		int len2 = p2.size();
		polynomial<T> result(Matrix<T>(1, len1 + len2 - 1, '0'));
		for (int j = 0; j < len2; ++j) {
			if (p2(j) != 0){
				for (int i = 0; i < len1; ++i) {
					result(i + j)+= p1(i) * p2(j);
				}
			}
		}
		result.simplify();
		return result;
	}

	friend polynomial<T> operator * (const T& num, polynomial<T> p) {
		if (num != 0) {
			int len = p.size();
			for (int i = 0; i < len; ++i) {
				p(i) *= num ;
			}
			return p;
		}
		else {
			return T(0);
		}
		
	}
	friend polynomial<T> operator * (const polynomial<T>& p, const T& num) {
		return num * p;		// switching property of *
	}

	friend polynomial<T> operator / (polynomial<T> p1, const polynomial<T>& p2) {
		// return the quotient only, p1/p2
		int len2 = p2.size();
		if (len2 == 1 && p2(0) == 0)
			throw "divide by 0";

		int len1 = p1.size();

		if (len1 >= len2) {
			polynomial<T> result(len1 - len2 + 1, '0');

			while (len1 >= len2 && !p1.is_approach_0()) {
				T num = p1(len1 - 1) / p2(len2 - 1);
				result(len1 - len2)= num;
				for (int k = 0; k < len2; ++k) {
					p1(len1 - 1 - k) -= num * p2(len2 - 1 - k);
				}

				p1.simplify();
				len1 = p1.size();
			}
			return result;
		}
		else {
			return T(0);
		}
	}
	friend polynomial<T> operator / (const polynomial<T>& p, const T& num) {
		if (num != 0) {
			return (1 / num) * p;
		}
		else {
			throw "divide by 0";
		}
		
	}

	friend polynomial<T> operator % (polynomial<T> p1, const polynomial<T>& p2) {
		// return the remainder only, p1%p2
		int len2 = p2.size();
		if (len2 == 1 && p2(0) == 0)
			throw "divide by 0";

		int len1 = p1.size();

		while (len1 >= len2 && !p1.is_approach_0()) {
			T num = p1(len1 - 1) / p2(len2 - 1);
			for (int k = 0; k < len2; ++k) {
				p1(len1 - 1 - k) -= num * p2(len2 - 1 - k);
			}
			p1.simplify();
			len1 = p1.size();
		}
		return p1;
	}

	/* define operator {<<,>>} */
	friend ostream& operator << (ostream& out, const polynomial<T>& p) {
		out << p.coeff;
		return out;
	}
	friend istream& operator >> (istream& in, const polynomial<T>& p) {
		in >> p.coeff;
		return in;
	}

	/* define operator {==,!=} */
	friend bool operator == (const polynomial<T>& p1, const polynomial<T>& p2) {
		return p1.coeff == p2.coeff;
	}
	friend bool operator != (const polynomial<T>& p1, const polynomial<T>& p2) {
		return p1.coeff != p2.coeff;
	}

	/* define evaluate function */
	/**
	 * . compute the polynomial given x in field T
	 * 
	 * \param x: the varible to plug into the polynomial
	 * \return the evaluate result
	 */
	T evaluate(const T& x) const {
		int coeff_len = coeff.size();
		T result = coeff(coeff_len - 1);
		for (int i = coeff_len - 2; i >= 0; --i) {
			result *= x;
			result += coeff(i);
		}
		return result;
	}
	/**
	 * .find roots for the polynomial, by searching all elements in {0,1,...,field_cardinality-1}
	 * 
	 * \return a Matrix<T> containing all roots
	 */
	Matrix<T> find_roots(int field_cardinality) const {
		Matrix<T> result(1, 0, 'v');

		// searching
		for (int i = 0; i < field_cardinality; ++i) {
			T can = T(i);
			if (evaluate(can) == 0) {
				result.push_back(can);
			}
		}
		return result;
	}
	
	/* 1-dimensional Hasse derivative */
	T Hasse_D(int r, const T& alpha) const {
		if (alpha == 0) {
			return coeff(0, r);
		}
		T result = 0;
		int co = coeff.col();
		for (int j = r; j < co; ++j) {		// over columns
			result += my::n_choose_k(j, r) * coeff(0, j) * pow(alpha, j - r);
		}
		return result;
	}

	/* for matrix polynomial computation, we introduce abs */
	friend bool operator < (const polynomial<T>& p1, const polynomial<T>& p2) {
		int p1s = p1.size();
		int p2s = p2.size();
		if (p1s < p2s) {			// polynomial with little size are assumed larger
			return false;
		}
		else if (p1s == p2s) {		// treat polynomials as element in T
			return my::abs(p1(p1s - 1)) < my::abs(p2(p2s - 1));
		}
		else {
			return true;
		}
	}
	// note that this is not simply < || ==
	friend bool operator <= (const polynomial<T>& p1, const polynomial<T>& p2) {
		int p1s = p1.size();
		int p2s = p2.size();
		if (p1s < p2s) {			// polynomial with little size are assumed larger
			return false;
		}
		else if (p1s == p2s) {		// treat polynomials as element in T
			return my::abs(p1(p1s - 1)) <= my::abs(p2(p2s - 1));
		}
		else {
			return true;
		}
	}
	friend bool operator > (const polynomial<T>& p1, const polynomial<T>& p2) {
		int p1s = p1.size();
		int p2s = p2.size();
		if (p1s < p2s) {		// this will make 0 greater than all 1 x-ordered polynomials, not ture 
			return true;
		}
		else if (p1s == p2s) {
			return my::abs(p1(p1s - 1)) > my::abs(p2(p2s - 1));
		}
		else {
			return false;
		}
	}
	friend bool operator >= (const polynomial<T>& p1, const polynomial<T>& p2) {
		int p1s = p1.size();
		int p2s = p2.size();
		if (p1s < p2s) {
			return true;
		}
		else if (p1s == p2s) {
			return my::abs(p1(p1s - 1)) >= my::abs(p2(p2s - 1));
		}
		else {
			return false;
		}
	}
	bool is_approach_0() const{
		return size() == 1 && (double)my::abs(coeff(0)) < my::zero_approximation;
	}

	explicit operator double() const {
		return (double)coeff(0);		// may incur problem when the first term of coeff is zero
	}
};

class Matrix_polynomials{
public:
	template<class T>
	static void to(const Matrix<T>& src, Matrix<polynomial<T>>& dst) {
		// make sure that src and dst have the same size
		int sr = src.row();
		int sc = src.col();
		if (sr == dst.row() && sc == dst.col()) {
			for (int i = 0; i < sr; ++i) {
				for (int j = 0; j < sc; ++j) {
					dst(i, j) = src(i, j);
				}
			}
		}
		else {
			cout << "size error, to_polynomial failed" << endl;
			cout << "src.row()=" << sr << ",\tsrc.col()=" << sc << endl;
			cout << "dst.row()=" << dst.row() << ",\tdst.col()=" << dst.col() << endl;
		}
	}

	template<class T>
	static void from(const Matrix<polynomial<T>>& src, const T& x, Matrix<T>& dst) {
		// make sure that src and dst have the same size
		int sr = src.row();
		int sc = src.col();
		if (sr == dst.row() && sc == dst.col()) {
			for (int i = 0; i < sr; ++i) {
				for (int j = 0; j < sc; ++j) {
					dst(i, j) = src(i, j).evaluate(x);
				}
			}
		}
		else {
			cout << "size error, to_polynomial failed" << endl;
			cout << "src.row()=" << sr << ",\tsrc.col()=" << sc << endl;
			cout << "dst.row()=" << dst.row() << ",\tdst.col()=" << dst.col() << endl;
		}
	}
};
