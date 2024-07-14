/*****************************************************************//**
 * \file	Matrix.h
 * \brief	matrix implementation, must include "my.h"
 * \see		we cannot seperate implementation and declear under templete function
 *			so we should put them in one file, this make this head file long, PAY ATTENTION
 * 
 * \todo	(1) the class implementation has too much space declaration and destruction,
 *				try to define new interface that the 'ans' parameter is pass through reference.		// donot do it, a waste of time
 *				the effect of this to the run time of program is still to be validate,
 *				since in release mode the optimization is not well understood. 
 * 
 * 
 * \author	lilili
 * \date	October 2022
 *********************************************************************/
#pragma once

#include<iomanip>
#include<vector>		// already contain a lot
#include<unordered_map>
#include<iostream>
#include<math.h>
#include"my.h"
#include"Complex.h"		// this can be not included, for GF use as convenient
using namespace std;

#ifdef Complex_class
#else
#define Complex_class
class Complex {
public:
	static const double J;
};
const double Complex::J = 0;
#endif // Complex_class

//#define see_eig_iteration_num
#define see_eig_warning

// unless you know the consequence exactly, please do not use memcpy and memset,
// if so, the "=" operator will not called in for elements,
// since matrix may contain pointer, release the pointer unintendedly will cause problem

//Matrix class, if you have problem with ordered structure, i.e., operator '<' not defined, please use std standard class vector
// or if you can not trun int to type T, use std standard class vector also
template<class T>
class Matrix
{
protected:
	int r;
	int c;
	T* e;
	// a int varible indicating decleared size, and push_back, pop_back function to make it a vector like
	int capacity;		// but this will induce out of index error not been found, this is special

	// Access the individual elements for function of class, optional
	inline T& entry(const unsigned& pos_ind) {
		return e[pos_ind];
	}
	inline const T& entry(const unsigned& pos_ind) const {
		return e[pos_ind];
	}
	inline T& entry(const unsigned& row_ind, const unsigned& col_ind) {
		return e[row_ind * c + col_ind];
	}
	inline const T& entry(const unsigned& row_ind, const unsigned& col_ind) const {
		return e[row_ind * c + col_ind];
	}

public:

	//constructors and destructor
	Matrix() :r(1), c(0), e(NULL), capacity(0){}
	// accept one parameter and construct a matrix with size 1 by 1, for convinience only
	Matrix(const T& element):r(1),c(1) {
		e = new T[1];
		e[0] = element;
		capacity = 1;
	}
	Matrix(int row, int column) :r(row), c(column) {
		capacity = r * c;
		e = new T[capacity];
	}
	Matrix(int row, int column, T* elements) :r(row), c(column) {
		capacity = r * c;
		e = new T[capacity];
		for (int i = 0; i < capacity; ++i) {
			e[i] = elements[i];
		}
	}	
	/**
	 * .initialize matrix
	 *
	 * \param row
	 * \param column
	 * \param v: the vector to copy, NOTE if using {} to replace 'v', then {}
	 * with 1 element like {5}, is not allowed and this will set matrix with 5 0s
	 */
	Matrix(int row, int column, const vector<T>& v) :r(row), c(column) {
		// not that this function is not valid for v={1}, that
		// is v set with {} and has 1 zero element !!!!!!!!
		// really creazy !!
		capacity = r * c;
		e = new T[capacity];
		for (int i = 0, imax = my::min((int)v.size(), capacity); i < imax; ++i) {
			// make sure not overflow
			e[i] = v[i];
		}
	}
	Matrix(int row, int column, char type_0_1_i_b_N_v) :r(row), c(column) {
		capacity = r * c;
		e = new T[capacity];

		switch (type_0_1_i_b_N_v) {
		case '0':							/* zero */
			for (int i = 0; i < capacity; i++)
				e[i] = T(0);
			break;
		case '1':							/* one */
			for (int i = 0; i < capacity; i++)
				e[i] = T(1);
			break;
		case 'i':							/* identity */
			for (int i = 0; i < capacity; i++)
				e[i] = T(0);
			for (int i = 0, min_rc = my::min(r, c), ind = 0; i < min_rc; ++i, ind += c + 1)
				e[ind] = T(1);
			break;
		case 'b':							/* binary: uniform distribution of 0 and 1 */
			for (int i = 0; i < capacity; i++)
				e[i] = T(my::rand_01());
			break;
		case 'N':							/* Natrual number in sequence start from 0 */
			for (int i = 0; i < capacity; i++)
				e[i] = T(i);
			break;
		case 'v':							/* initialize it as an empty row vector */
			r = 1;		// it must be a row vector
			c = 0;		// let the column be 0, starting position for push_back
			break;
		default:
			cout << "warning: Matrix type undefined. type_0_1_i_b_N_v = " << type_0_1_i_b_N_v << endl;
		}
	}
	Matrix(const T& start, const T& diff, const T& end, char type_d){
		if (type_d != 'd') {
			// not giving an identifying message
			c = 0;
			r = 0;
			e = NULL;
			capacity = 0;
			return;
		}

		r = 1;		// create a row vector
		c = (int)(my::abs(end - start) / my::abs(diff)) + 1;
		e = new T[c];
		capacity = c;

		e[0] = start;
		for (int i = 1; i < c; ++i) {
			e[i] = e[i - 1] + diff;
		}
	}
	// to do: the constructor function with (int,int,function pointer) to call the function to initialize every element in Matrix
	Matrix(int row, int column, T(*generate_element)(void)) :r(row), c(column) {
		capacity = r * c;
		e = new T[capacity];
		for (int i = 0; i < capacity; ++i) {
			e[i] = generate_element();
		}
	}

	Matrix(const Matrix<T>& A) {
		r = A.r;
		c = A.c;
		capacity = A.capacity;
		e = new T[capacity];
		for (int i = 0, imax = size(); i < imax; i++) {
			e[i] = A.e[i];
		}
	}
	Matrix<T>& operator = (const Matrix<T>& A) {
		if (this != &A) {
			// here we should compare the declear_size of A
			if (capacity >= A.capacity);		// if original size is enough, donot need to delete and new			
			else {
				delete[] e;
				capacity = A.capacity;			// make the declear size equal
				e = new T[capacity];
			}

			r = A.r;
			c = A.c;
			for (int i = 0, imax = size(); i < imax; ++i) {
				e[i] = A.e[i];
			}
		}
		return *this;
	}
	~Matrix() {
		// before this, destructor of e[i] has already been called, hence it will be okay to destroy all the element, no memory waste
		delete[] e;
	}

	Matrix<T>& cp_valid(const Matrix<T>& A) {
		if (this != &A) {
			r = A.r;
			c = A.c;
			int new_size = r * c;

			// here we should compare the size of A
			if (capacity >= new_size);		// if original size is enough, donot need to delete and new			
			else {
				delete[] e;
				capacity = new_size;		// make the declear size equal to the size of A
				e = new T[capacity];
			}

			for (int i = 0; i < new_size; ++i) {
				e[i] = A.e[i];
			}
		}
		return *this;
	}
	void clear_space() {				// release all the declared memory, just for inteface consistance with vector 
		r = 1;
		c = 0;
		delete[] e;
		e = NULL;
		capacity = 0;		// clear out the space, same as destructor except for keep the variable's name
	}
	inline void clear() {
		r = 1;
		c = 0;
	}
	inline bool empty() {
		return c == 0 || r == 0;
	}
	
	// set all element in the matrix be 0
	void reset(const T& elements = T()) {
		for (int i = 0, imax = size(); i < imax; ++i) {
			e[i] = elements;
		}
	}

	// vector function, make sure r=1
	void push_back(const T& element) {
		if (c < capacity);
		else if (capacity != 0) {
			// space is not enough, double the current space
			T* tmp = e;
			e = new T[capacity * 2];
			for (int i = 0; i < capacity; ++i) {
				e[i] = tmp[i];
			}

			capacity *= 2;			// this is not recommand for speed crucial program, set a interuption here to check
			delete[] tmp;	
			// this do not call destructor of the elements inside temp, hence not release memory if elements are not trivially destructible
		}
		else {
			// declear space of 10 element at the first time

			delete[] e;		// this should be deleted, since e is initiallized with 0 size of array, taking spaces
			e = new T[10];
			capacity = 10;
		}
		e[c] = element;
		c++;
	}
	inline void pop_back() {
		c--;
		if (c >= 0);
		else {
			cout << "warning: under flow of pop_back(), setting c be 0" << endl;
			c = 0;
		}
	}
	inline T& back() {
		if (c > 0)
			return e[c - 1];
		else {
			cout << "warning: under flow of back(), no element at back" << endl;
			return e[0];
		}
	}
	inline const T& back() const{
		if (c > 0)
			return e[c - 1];
		else {
			cout << "warning: under flow of back(), no element at back" << endl;
			return e[0];
		}
	}

	// Access the individual elements, we define both [] operator and () operator, for convinience
	inline T& operator[](const unsigned& pos_ind) {
		return e[pos_ind];
	}
	inline const T& operator[](const unsigned& pos_ind) const {
		return e[pos_ind];
	}
	inline T& operator()(const unsigned& pos_ind) {
		return e[pos_ind];
	}
	inline const T& operator()(const unsigned& pos_ind) const {
		return e[pos_ind];
	}
	// in high performance vital codes, consider to replace the following two () functions with the above two [] functions
	inline T& operator()(const unsigned& row_ind, const unsigned& col_ind) {
		return e[row_ind * c + col_ind];
	}
	inline const T& operator()(const unsigned& row_ind, const unsigned& col_ind) const {
		return e[row_ind * c + col_ind];
	}
	
	inline pair<int, int> locate(const unsigned& pos_ind) const {
		pair<int, int> ans;
		ans.first = pos_ind / c;
		ans.second = pos_ind - ans.first * c;
		return ans;
	}

	// return the pointer e
	inline const T* ptr() const{
		return e;
	}
	// return the number of row
	inline int row() const {
		return r;
	}
	// return the number of column
	inline int col() const {
		return c;
	}
	// return the number of all element
	inline int size() const {
		return r * c;
	}
	//return the [row or column] vector of all elements
	Matrix<T> vectorize(bool to_row = true) const {
		int sss = size();
		Matrix<T> ans(1, sss);
		if (to_row);
		else {
			// to a column vector
			ans.r = sss;
			ans.c = 1;
		}
		for (int i = 0; i < sss; ++i) {
			ans.e[i] = e[i];
		}
		return ans;
	}
	//return the a certain part in matrix, -1 for row_end_ind indicates r-1, and -1 for column_end_ind indicates c-1
	Matrix<T> get_part(int row_start_ind = -1, int column_start_ind = -1, int row_end_ind = -1, int column_end_ind = -1) const {

		row_start_ind = row_start_ind == -1 ? r - 1 : row_start_ind;
		column_start_ind = column_start_ind == -1 ? c - 1 : column_start_ind;

		row_end_ind = row_end_ind == -1 ? r - 1 : row_end_ind;
		column_end_ind = column_end_ind == -1 ? c - 1 : column_end_ind;

		if (row_start_ind == 0 && column_start_ind == 0 && row_end_ind == r - 1 && column_end_ind == c - 1) {
			// all get
			return (*this);
		}

		Matrix<T> result(row_end_ind - row_start_ind + 1, column_end_ind - column_start_ind + 1);
		int k = 0;
		int row_pos = row_start_ind * c;
		for (int i = row_start_ind; i <= row_end_ind; i++) {
			for (int j = column_start_ind; j <= column_end_ind; j++) {
				result.e[k] = e[row_pos + j];
				k++;
			}
			row_pos += c;
		}
		return result;
	}
	// the inplace operation of get_part, please allocate enough size for result, for inplace operation, note that &result != this
	void get_part(int row_start_ind, int column_start_ind, int row_end_ind, int column_end_ind, Matrix<T>& result) const {

		if (row_start_ind == 0 && column_start_ind == 0 && row_end_ind == r - 1 && column_end_ind == c - 1) {
			// all get
			result = (*this);
			return;
		}

		result.resize(row_end_ind - row_start_ind + 1, column_end_ind - column_start_ind + 1, false);
		int k = 0;
		int row_pos = row_start_ind * c;
		for (int i = row_start_ind; i <= row_end_ind; i++) {
			for (int j = column_start_ind; j <= column_end_ind; j++) {
				result.e[k] = e[row_pos + j];
				k++;
			}
			row_pos += c;
		}
	}
	
	//return the row vector of diagonal elements
	Matrix<T> get_diag(bool ret_col = true) const {
		if (r <= c) {
			Matrix<T> result(r, 1);
			if (ret_col);
			else
				result.resize(1, r, false);		// ret_col = false, then return a row vector

			for (int i = 0; i < r; i++)
				result(i) = e[i * c + i];
			return result;
		}
		else {
			Matrix<T> result(c, 1);
			if (ret_col);
			else
				result.resize(1, c, false);

			for (int j = 0; j < c; j++)
				result(j) = e[j * c + j];
			return result;
		}
	}
	//return the row vector in row 'row_ind'
	Matrix<T> get_row(int row_ind) const {
		Matrix<T> result(1, c);
		int row_pos = row_ind * c;
		for (int i = 0; i < c; ++i) {
			result.e[i] = e[row_pos + i];
		}
		return result;
	}
	//return the column vector in column 'column_ind'
	Matrix<T> get_col(int column_ind) const {
		Matrix<T> result(r, 1);
		for (int i = 0; i < r; i++) {
			result(i) = e[column_ind];
			column_ind += c;
		}
		return result;
	}
	// return the rows in the Matrix, indicating by row_ind
	Matrix<T> get_rows(const Matrix<int>& row_ind) const {
		int row_ind_size = row_ind.size();
		if (row_ind_size == r) {		// all get
			return (*this);
		}

		Matrix<T> result(row_ind_size, c);
		for (int i = 0; i < row_ind_size; i++) {
			for (int j = 0; j < c; ++j) {
				result(i, j) = (*this)(row_ind(i), j);
			}
		}
		return result;
	}
	// return the colums in the Matrix, indicating by column_ind
	Matrix<T> get_cols(const Matrix<int>& column_ind) const {
		int column_ind_size = column_ind.size(); 
		if (column_ind_size == c) {		// all get
			return (*this);
		}

		Matrix<T> result(r, column_ind_size);
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < column_ind_size; ++j) {
				result(i, j) = (*this)(i, column_ind(j));
			}
		}
		return result;
	}
	// inplace operation for get_cols
	/**
	 * .\param row_end_ind: rows from index [0, row_end_ind] will be copied, by default it will copy all rows
	 */
	void get_cols(const Matrix<int>& column_ind, Matrix<T>& result, int row_end_ind = -1) const {
		int column_ind_size = column_ind.size();
		row_end_ind = row_end_ind == -1 ? r - 1 : row_end_ind;

		if (column_ind_size == c) {		// all get
			(*this).get_part(0, 0, row_end_ind, c, result);
			return;
		}

		result.resize(row_end_ind + 1, column_ind_size, false);
		for (int i = 0; i <= row_end_ind; i++) {
			for (int j = 0; j < column_ind_size; ++j) {
				result(i, j) = (*this)(i, column_ind(j));
			}
		}
	}

	/**
	 * .put the element back as m_ind indicates, inverse of permute
	 *
	 * \param m_ind: the element of position i will occupy position m_ind(i)
	 * e.g., {a,b,c}.permute_back({2,0,1}) gives {b,c,a}
	 */
	void permute_back(const Matrix<int>& m_ind) {
		int len = size();
		int m_ind_len = m_ind.size();
		if (m_ind_len <= len && m_ind_len != 0) {
			Matrix<T> result(*this);
			for (int i = 0; i < m_ind_len; ++i) {
				result(m_ind(i)) = (*this)(i);		// partial permute back for first part
			}
			(*this) = result;
		}
		else {
			cout << "m_ind_len greater than matrix size, permute_back failed" << endl;
			cout << "m_ind_len = " << m_ind_len << ", matrix size = " << len << endl;
		}
	}

	/**
	 * .put the element as m_ind indicates, inverse of permute_back
	 *
	 * \param m_ind: the element of position m_ind(i) will occupy position i
	 * e.g., {a,b,c}.permute({2,0,1}) gives {c,a,b}
	 */
	void permute(const Matrix<int>& m_ind) {
		int len = size();
		int m_ind_len = m_ind.size();
		if (m_ind_len <= len) {
			Matrix<T> result(1, m_ind_len);
			for (int i = 0; i < m_ind_len; ++i) {
				result(i) = (*this)(m_ind(i));		// partial permute for the first part
			}
			for (int i = 0; i < m_ind_len; ++i) {
				e[i] = result.e[i];
			}
		}
		else {
			cout << "m_ind_len greater than matrix size, permute failed" << endl;
			cout << "m_ind_len = " << m_ind_len << ", matrix size = " << len << endl;
		}
		/*for (int i = m_ind_len; i < len; ++i) {
			result(i) = (*this)(i);				
		}*/
	}

	/**
	 * .put the column elements as m_ind indicates, inverse of permute_back
	 *
	 * \param m_ind: the element of position m_ind(i) will occupy position i
	 * e.g., {a,b,c}.permute({2,0,1}) gives {c,a,b}
	 */
	void permute_col(const Matrix<int>& m_ind) {
		int m_ind_len = m_ind.size();
		if (m_ind_len <= c) {
			Matrix<T> result(*this);
			for (int i = 0; i < m_ind_len; ++i) {
				if (i != m_ind(i)) {				// not permuted column can be skipped
					for (int j = 0; j < r; ++j) {
						result(j, i) = (*this)(j, m_ind(i));	// partial permute for the first part
					}
				}
			}
			(*this) = result;
		}
		else {
			cout << "m_ind_len greater than matrix column number, permute_col failed" << endl;
			cout << "m_ind_len = " << m_ind_len << ", matrix column number = " << c << endl;
		}
	}
	/**
	 * .put the column elements as m_ind indicates, inverse of permute_back
	 *
	 * \param m_ind: the element of position m_ind(i) will occupy position i
	 */
	void permute_row(const Matrix<int>& m_ind) {
		int m_ind_len = m_ind.size();
		if (m_ind_len <= r) {
			Matrix<T> result(*this);
			for (int i = 0; i < m_ind_len; ++i) {
				if (i != m_ind(i)) {				// not permuted column can be skipped
					for (int j = 0; j < c; ++j) {
						result(i, j) = (*this)(m_ind(i), j);	// partial permute for the first part
					}
				}
			}
			(*this) = result;
		}
		else {
			cout << "m_ind_len greater than matrix row number, permute_row failed" << endl;
			cout << "m_ind_len = " << m_ind_len << ", matrix row number = " << r << endl;
		}
	}
	void permute_col_back(const Matrix<int>& m_ind) {
		Matrix<int> permute_back_record(1, m_ind.size(), 'N');
		permute_back_record.permute_back(m_ind);
		permute_col(permute_back_record);			// that is simple compound method
	}
	void permute_row_back(const Matrix<int>& m_ind) {	// not test yet
		Matrix<int> permute_back_record(1, m_ind.size(), 'N');
		permute_back_record.permute_back(m_ind);
		permute_row(permute_back_record);			// that is simple compound method
	}
	/**
	 * .randomly permute the elements of Matrix
	 */
	void permute_rand() {
		int sss = size();
		Matrix<int> m_ind(1, sss, 'N');
		m_ind = m_ind.get_random_element(sss);
		permute(m_ind);
	}

	/**
	 * .return matrix with each element take abs from (*this) matrix
	 */
	Matrix<T> get_abs() const {
		Matrix<T> result(r, c);

		for (int i = 0, imax = size(); i < imax; ++i) {
			result.e[i] = my::abs(e[i]);
		}
		return result;
	}
	/**
	 * .get num elements from (*this) matrix, not overlapping
	 */
	Matrix<T> get_random_element(int num = 1) const{
		int sss = size();
		if (num <= sss) {
			Matrix<T> tmp(*this);		// to prevent change Matrix (*this)
			Matrix<T> result(1, num, 'v');
			for (int i = 0; i < num; ++i) {
				//int choose_ind = my::rand_int_adv(0, sss - 1);
				int choose_ind = rand() % sss;	// not using the advanced random generator to keep the seed consist
				result.push_back(tmp(choose_ind));
				tmp.switch_ele(choose_ind, sss - 1);
				tmp.pop_back();
				sss--;
			}
			return result;
		}
		else {
			cout << "(get_random_element) warning: num > size()" << endl;
			cout <<"num = " << num << ", size() = " << size() << endl;
			return Matrix<T>();
		}
	}
	/**
	 * .if the inner product of v and each row of the (*this) matrix is zero, return true, else return false
	 */
	bool check_inner_product(const Matrix<T>& v) const {
		if (c == 0) {
			return true;
		}
		else if (v.col() == c) {
			int row_pos = 0;
			for (int i = 0; i < r; ++i) {
				T inner_product = 0;
				for (int j = 0; j < c; ++j) {
					inner_product += v(j) * e[row_pos + j];
				}
				if (inner_product == 0);
				else{
					return false;
				}
				row_pos += c;
			}
			return true;
		}
		else {
			cout << "warning: column not mach. v.col()=" << v.col() << ", c=" << c << endl;
			return false;
		}
	}
	/**
	 * .return the index where the elements of (*this) is non zero
	 */
	Matrix<int> get_non_0_ind() const {
		int sss = size();
		Matrix<int> ans(1, sss, 'v');
		for (int i = 0; i < sss; ++i) {
			if (e[i] != 0) {
				ans.push_back(i);
			}
		}
		return ans;
	}
	/**
	 * .'ans' be the index where the elements of (*this) is non zero
	 */
	void diff_ind_inplace(const Matrix<T>& ref, Matrix<int>& ans) const {
		int sss = size();
		int refs = ref.size();

		ans.resize(1, 0, false);
		if (refs == 0) {		// consider it as all zero vector
			for (int i = 0; i < sss; ++i) {
				if (e[i] != 0) {
					ans.push_back(i);
				}
			}
		}
		else {		// compare the first number of min{refs,sss} to generate ans
			int min_size = my::min(refs, sss);
			for (int i = 0; i < min_size; ++i) {
				if (e[i] != ref(i)) {
					ans.push_back(i);
				}
			}
		}
	}
	/**
	 * . to be fixed over Viterbi class
	 */
	bool check_inner_product_4_GF2(const Matrix<T>& v, const Matrix<T>& ref, \
		Matrix<int>& non_0_ind_aux, const Matrix<T>& ref_result = Matrix<T>()) const {
		if (c == 0) {
			return true;
		}

		bool ref_result_unset = ref_result.size() == 0;		
		v.diff_ind_inplace(ref, non_0_ind_aux);

		int nsss = non_0_ind_aux.size();
		int row_pos = 0;
		for (int i = 0; i < r; ++i) {
			T inner_product = 0;
			for (int j = 0; j < nsss; ++j) {
				inner_product += e[row_pos + non_0_ind_aux(j)];
			}
			if ((ref_result_unset && inner_product == 0) \
				|| (!ref_result_unset && inner_product == ref_result(i)));
			else {
				return false;
			}
			row_pos += c;
		}
		return true;
	}

	/**
	 * . function the specified for inner porduct of a 'vector' over a Matrix, can be used to Matrix over Matrix
	 * 
	 * \param B: the Matrix to be multiplied, saving the step of transpose B
	 * \return (*this)*B.Transpose() (actually the B.Hermitian() is more useful, not implement here)
	 */
	Matrix<T> multiply_transpose_of(const Matrix<T>& B) const {
		if (c == B.c) {
			Matrix<T> ans(r, B.r, '0');
			int row_pos;
			int row_pos2;

			for (int t = 0; t < c; t++) {
				row_pos = 0;
				row_pos2 = 0;
				for (int i = 0; i < r; i++) {	// the compiler will parallize it, using release mode

					//theoretically will have parallize degree A.r * B.c

					for (int j = 0; j < B.r; j++) {
						ans.e[row_pos + j] += e[row_pos2 + t] * B(j, t);
					}
					row_pos2 += c;
					row_pos += B.r;
				}
			}
			return ans;
		}
		else {
			cout << "size error (c == B.c), multiply_transpose_of failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			return Matrix<T>();
		}
	}
	/**
	 * .sum of columns indexed by 'column_ind' has the 
	 */
	Matrix<T> sum_columns(const Matrix<int>& column_ind) const {

		Matrix<T> ans(1, r, '0');

		// check out each row of H
		int column_num = column_ind.size();
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < column_num; ++j) {
				ans(i) += (*this)(i, column_ind(j));
			}
		}

		return ans;
	}
	/**
	 * .if sum of columns indexed by 'column_ind' == ref, return true, else return false
	 */
	bool sum_columns_is_ref(const Matrix<int>& column_ind, const Matrix<T>& ref = Matrix<T>()) const {

		// check out each row of H
		int column_num = column_ind.size();
		for (int i = 0; i < r; ++i) {
			T check = T(0);
			for (int j = 0; j < column_num; ++j) {
				check += (*this)(i, column_ind(j));
			}
			if (check != ref(i)) {
				return false;
			}
		}
		return true;
	}

	// the next 6 functions return the element with [max/min] [my::abs] value

	// these fucntions are not suitable for GF matrix as element has no order, pay attention

	// return max element and its index
	int max_ele_ind() const {
		int sss = size();
		if (sss > 0) {
			int ans = 0;
			T max_can = e[0];
			T can;
			bool update;
			for (int i = 1; i < sss; ++i) {
				can = e[i];
				update = max_can < can;
				ans = update ? i : ans;
				max_can = update ? can : max_can;
			}
			return ans;
		}
		else {
			cout << "warning: an empty matrix calling max_ele_ind()" << endl;
			cout << "return 0" << endl;
			return 0;
		}
	}
	int min_ele_ind() const {
		int sss = size();
		if (sss > 0) {
			int ans = 0;
			T min_can = e[0];
			T can;
			bool update;
			for (int i = 1; i < sss; ++i) {
				can = e[i];
				update = min_can > can;
				ans = update ? i : ans;
				min_can = update ? can : min_can;
			}
			return ans;
		}
		else {
			cout << "warning: an empty matrix calling min_ele_ind()" << endl;
			cout << "return 0" << endl;
			return 0;
		}
	}

	// return max abs element and its index
	int max_abs_ele_ind() const {
		int sss = size();
		if (sss > 0) {
			int ans = 0;
			T max_can = my::abs(e[0]);
			T can;
			bool update;
			for (int i = 1; i < sss; ++i) {
				can = my::abs(e[i]);
				update = max_can < can;
				ans = update ? i : ans;
				max_can = update ? can : max_can;
			}
			return ans;
		}
		else {
			cout << "warning: an empty matrix calling max_ele_ind()" << endl;
			cout << "return 0" << endl;
			return 0;
		}
	}
	int min_abs_ele_ind() const {
		int sss = size();
		if (sss > 0) {
			int ans = 0;
			T min_can = e[0];
			T can;
			bool update;
			for (int i = 1; i < sss; ++i) {
				can = my::abs(e[i]);
				update = min_can > can;
				ans = update ? i : ans;
				min_can = update ? can : min_can;
			}
			return ans;
		}
		else {
			cout << "warning: an empty matrix calling min_ele_ind()" << endl;
			cout << "return 0" << endl;
			return 0;
		}
	}

	//return the max element in matrix
	T max_ele() const {
		return e[max_ele_ind()];
	}
	T max_abs_ele() const {
		return e[max_abs_ele_ind()];
	}

	//return the min element in matrix
	T min_ele() const {
		return e[min_ele_ind()];
	}
	T min_abs_ele() const {
		return e[min_abs_ele_ind()];
	}

	// the next 6 function is for special use

	// return the max element and its position (position_ind) in matrix start from the element with positoin (start_position_ind)
	T max_abs_ele_start_from(int start_position_ind, int& position_ind) const {

		T max = e[start_position_ind];
		T max_abs = my::abs(max);
		T temp_abs;
		position_ind = start_position_ind;
		int sss = size();
		for (int j = start_position_ind + 1; j < sss; j++) {
			temp_abs = my::abs(e[j]);
			if (max_abs >= temp_abs);
			else {
				max = e[j];
				max_abs = temp_abs;
				position_ind = j;
			}
		}
		return max;
	}
	// this function is not valid for GF since my::abs() and < operator is not valid for GF2
	T max_abs_col_ele_start_from(int column_ind, int row_start_ind, int& row_ind) const {
		int start_pos = row_start_ind * c + column_ind;
		T max = e[start_pos];
		T max_abs = my::abs(max);
		T temp_abs;
		row_ind = row_start_ind;
		for (int i = row_start_ind + 1; i < r; i++) {
			start_pos += c;
			temp_abs = my::abs(e[start_pos]);
			if (max_abs < temp_abs) {
				max = e[start_pos];
				max_abs = temp_abs;
				row_ind = i;
			}
		}
		return max;
	}

	// return the max element and its position (position_ind) in matrix end from the element with positoin (end_position_ind)
	T max_abs_ele_end_from(int end_position_ind, int& position_ind) const {

		T max = e[end_position_ind];
		T max_abs = my::abs(max);
		T temp_abs;
		position_ind = end_position_ind;
		for (int j = end_position_ind - 1; j >= 0; j--) {
			temp_abs = my::abs(e[j]);
			if (max_abs >= temp_abs);
			else {
				max = e[j];
				max_abs = temp_abs;
				position_ind = j;
			}
		}
		return max;
	}
	// this function is not valid for GF since my::abs() and < operator is not valid for GF2
	T max_abs_col_ele_end_from(int column_ind, int row_end_ind, int& row_ind) const {
		int start_pos = row_end_ind * c + column_ind;
		T max = e[start_pos];
		T max_abs = my::abs(max);
		T temp_abs;
		row_ind = row_end_ind;
		for (int i = row_end_ind - 1; i >= 0; i--) {
			start_pos -= c;
			temp_abs = my::abs(e[start_pos]);
			if (max_abs < temp_abs) {
				max = e[start_pos];
				max_abs = temp_abs;
				row_ind = i;
			}
		}
		return max;
	}

	// use this function instead of the last one, do not need partial pivoting in GF since no miscalculation in GF, pay attention for this!
	T non_0_col_ele_start_from(int column_ind, int row_start_ind, int& row_ind) const {
		// NOTE: this is changed for GF2 

		int start_pos = row_start_ind * c + column_ind;
		T max = T(0);
		for (int i = row_start_ind; i < r; i++) {
			max = e[start_pos];
			row_ind = i;
			if (max != T(0)) {
				return max;
			}
			start_pos += c;
		}
		return max;
	}
	// use this function instead of the last one, do not need partial pivoting in GF since no miscalculation in GF, pay attention for this!
	T non_0_col_ele_end_from(int column_ind, int row_end_ind, int& row_ind) const {
		// NOTE: this is changed for GF2 

		int start_pos = row_end_ind * c + column_ind;
		T max = 0;
		for (int i = row_end_ind; i >= 0; i--) {
			max = e[start_pos];
			row_ind = i;
			if (max != 0) {
				return max;
			}
			start_pos -= c;
		}
		return max;
	}

	// for matrix polynomials only
	T max_col_ele_start_from_4_matrix_polynomials(int column_ind, int row_start_ind, int& row_ind) const {
		int start_pos = row_start_ind * c + column_ind;
		T max = e[start_pos];
		row_ind = row_start_ind;
		for (int i = row_start_ind + 1; i < r; i++) {
			start_pos += c;
			// should prevent max being less than my::zero_approximation, already comparing abs value
			if ((max.is_approach_0() ? false : e[start_pos] <= max) \
				|| e[start_pos].is_approach_0());
			else {
				max = e[start_pos];
				row_ind = i;
			}


		}
		return max;
	}

	// judge functions
	bool isDiagonal() const {
		if (r == c) {
			for (int i = 0; i < r; i++)
				for (int j = 0; j < c; j++)
					if (i != j)
						if ((double) my::abs(e[i * c + j]) < my::zero_approximation);			// wider condition, approach any num < my::zero_approximation be 0, last value: 1e-12
						else return false;
			return true;
		}
		else
			return false;
	}
	bool isSymmetric() const {
		if (r == c) {
			for (int i = 0; i < r; i++)
				for (int j = i + 1; j < c; j++)
					if (e[i * c + j] == e[j * c + i]);
					else return false;
			return true;
		}
		else
			return false;
	}
	bool isHermitian() const {
		if (r == c) {
			for (int i = 0; i < r; i++)
				for (int j = i; j < c; j++)		// diagonal elements need to be real
					if (e[i * c + j] == my::conj(e[j * c + i]));
					else return false;
			return true;
		}
		else
			return false;
	}
	bool isVector() const {
		return (r == 1 || c == 1);
	}
	bool isDeficient() const {
		// if deficient, row deficient or column deficient
		int _rank = rank();
		int min_r_c = r < c ? r : c;
		return _rank < min_r_c;
		// will not give which one, but tell whether deficient
	}
	bool isPositiveDefinite() const {
		if (r == c) {
			// to be done
			return true;
		}
		else {
			cout << "invalid size r!=c for isPositiveDefinite()" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return false;
		}
	}
	bool isSemiPositiveDefinite() const {
		if (r == c) {
			// to be done
			return true;
		}
		else {
			cout << "invalid size r!=c for isPositiveDefinite()" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return false;
		}
	}
	bool isZero() const {
		int sss = size();
		for (int i = 0; i < sss; ++i) {
			if (e[i] == 0);
			else {
				return false;
			}
		}
		return true;
	}

	// set function(change the matrix)

	// set size without changing matrix element
	void set_row(int row_ind, const Matrix<T>& elements) {
		int row_pos = row_ind * c;
		for (int i = 0; i < c; ++i) {
			e[row_pos + i] = elements.e[i];
		}
	}
	void set_col(int column_ind, const Matrix<T>& elements) {
		for (int i = 0; i < r; i++) {
			e[column_ind] = elements(i);
			column_ind += c;
		}
	}
	void set_part(int row_start_ind, int column_start_ind, const Matrix<T>& elements) {
		int er = elements.row();
		int ec = elements.col();
		int start_ind = row_start_ind * c + column_start_ind;
		for (int i = 0; i < er; ++i) {
			for (int j = 0; j < ec; ++j) {
				e[start_ind + j] = elements(i, j);
			}
			start_ind += c;
		}
	}
	void set_part(int row_start_ind, int column_start_ind, \
		int row_end_ind = -1, int column_end_ind = -1, const T& element = 0) {

		row_end_ind = row_end_ind == -1 ? r - 1 : row_end_ind;
		column_end_ind = column_end_ind == -1 ? c - 1 : column_end_ind;

		int row_ind = row_start_ind * c;
		for (int i = row_start_ind; i <= row_end_ind; ++i) {
			for (int j = column_start_ind; j <= column_end_ind; ++j) {
				e[row_ind + j] = element;
			}
			row_ind += c;
		}
	}
	/**
	 * . resize the matrix
	 * 
	 * \param row: row number of the matrix
	 * \param column: column number of the matrix
	 * \param keep_original: if keep the original element if a new space decleared, by default it keeps
	 */
	void resize(int row, int column, bool keep_original = true) {
		int size_2 = row * column;
		if (capacity >= size_2);		// if original size is enough, do not create a new array
		else if (!keep_original) {
			delete[] e;
			e = new T[size_2];
			capacity = size_2;		// capacity never shrink
		}
		else {
			int size_1 = size();
			T* temp = e;
			e = new T[size_2];
			for (int i = 0; i < size_1; ++i) {
				e[i] = temp[i];
			}
			delete[] temp;
			capacity = size_2;		// capacity never shrink
		}
		r = row;
		c = column;
	}
	void scale(const T& element) {	// a number times a matrix
		int size_all = size();
		for (int i = 0; i < size_all; i++)
			e[i] *= element;
	}
	void add_diag(const T& element) {
		int diag_pos = 0;
		if (r <= c)
			for (int i = 0; i < r; i++) {
				e[diag_pos] += element;
				diag_pos += (c + 1);
			}
		else
			for (int j = 0; j < c; j++) {
				e[diag_pos] += element;
				diag_pos += (c + 1);
			}
	}
	/**
	 * .
	 * 
	 * \param column_ind: indexes of column to be set
	 * \param element: the left part columns of elements will be plug into columns of (*this), according to column_ind
	 */
	void set_cols(const Matrix<int>& column_ind, const Matrix<T>& elements, bool duplicate_first_column_of_elements = false) {
		/**
		 * e.g.
		 *	(*this) =
		 	 1             0             0             0             0
             0             1             0             0             0
             0             0             1             0             0
             0             0             0             1             0
             0             0             0             0             1
		 * 
		 *  column_ind =
		     3			   4
		 
		 *  elements =
		     9             1             2
             3             4             5
             6             7             8
		 *
		 *  (*this) after = 
			 1             0             0             9             1
			 0             1             0             3             4
			 0             0             1             6             7
			 0             0             0             1             0
			 0             0             0             0             1
		 * */

		int cs = column_ind.size();
		if (cs != 0) {
			int ec = elements.col();
			int er = elements.row();
			int row_start = 0;
			int row_start2 = 0;
			if (!duplicate_first_column_of_elements) {
				for (int i = 0; i < er; i++) {
					for (int k = 0; k < cs; ++k) {
						e[row_start + column_ind(k)] = elements(row_start2 + k);// fetch the first 'cs' columns of element
					}
					row_start += c;
					row_start2 += ec;
				}
			}
			else {
				for (int i = 0; i < er; i++) {
					for (int k = 0; k < cs; ++k) {
						e[row_start + column_ind(k)] = elements(row_start2);	// fetch the first column of element only
					}
					row_start += c;
					row_start2 += ec;
				}
			}
		}		
	}
	// set this into an row_size * row_size identity Matrix
	void set_identity(int row_size, int col_size) {
		resize(row_size, col_size, false);
		int size_all = size();
		for (int i = 0; i < size_all; ++i) {
			e[i] = 0;
		}

		int min_rc = r < c ? r : c;
		for (int i = 0, ind = 0; i < min_rc; ++i, ind += c + 1) {
			e[ind] = 1;
		}
	}
	void set_natural() {
		for (int i = 0, imax = size(); i < imax; ++i) {
			e[i] = i;
		}
	}

	inline void switch_ele(int position1_ind, int position2_ind) {
		//if (position_1 == position_2) return;
		T temp = e[position1_ind];
		e[position1_ind] = e[position2_ind];
		e[position2_ind] = temp;
	}
	void switch_row(int row1_ind, int row2_ind) {
		if (row1_ind != row2_ind) {
			row1_ind = row1_ind * c;
			row2_ind = row2_ind * c;
			for (int j = 0; j < c; j++) {
				switch_ele(row1_ind, row2_ind);
				row1_ind++;
				row2_ind++;
			}
		}
	}
	void switch_col(int column1_ind, int column2_ind) {
		if (column1_ind != column2_ind) {
			for (int i = 0; i < r; i++) {
				switch_ele(column1_ind, column2_ind);
				column1_ind += c;
				column2_ind += c;
			}
		}
	}

	void unitize() {	// this is only valid for float and double
		int size_all = size();
		T sqrt_sum_abs_2 = 0;		// use double make it strange
		for (int i = 0; i < size_all; i++) {
			sqrt_sum_abs_2 += my::conj(e[i]) * e[i];
		}
		sqrt_sum_abs_2 = sqrt(sqrt_sum_abs_2);
		if (sqrt_sum_abs_2 >= my::zero_approximation) {
			for (int i = 0; i < size_all; i++) {
				e[i] /= sqrt_sum_abs_2;
			}
		}
	}
	void unitize_row(int row_ind) {		// this is only valid for float and double
		row_ind = row_ind * c;
		T sqrt_sum_abs_2 = 0;
		for (int i = 0; i < c; i++) {
			sqrt_sum_abs_2 += my::conj(e[row_ind + i]) * e[row_ind + i];
		}
		sqrt_sum_abs_2 = sqrt(sqrt_sum_abs_2);
		if (sqrt_sum_abs_2 >= my::zero_approximation) {
			for (int i = 0; i < c; i++) {
				e[row_ind + i] /= sqrt_sum_abs_2;
			}
		}
	}
	void unitize_col(int column_ind = -1) {	// this is only valid for float and double
		if (column_ind != -1) {
			int temp_col = column_ind;
			T sqrt_sum_abs_2 = 0;
			for (int j = 0; j < r; j++) {
				sqrt_sum_abs_2 += my::conj(e[temp_col]) * e[temp_col];
				temp_col += c;
			}
			sqrt_sum_abs_2 = sqrt(sqrt_sum_abs_2);		// this is improtant, you have fogotten it before
			if (sqrt_sum_abs_2 >= my::zero_approximation)
				for (int j = 0; j < r; j++) {
					e[column_ind] /= sqrt_sum_abs_2;
					column_ind += c;
				}
		}
		else {
			// unitize every column if column_ind is not given
			for (int i = 0; i < c; ++i) {
				unitize_col(i);
			}
		}
	}
	/**
	 * .reverse the matrix
	 *
	 */
	void rev() {
		int len = size() - 1;
		int len_over_2 = len / 2;
		for (int i = 0; i <= len_over_2; ++i) {
			switch_ele(i, len - i);
		}
	}

	/**
	 * .sort eleemnts from small to big, < as default, quick sort is the best choice
	 * 
	 * \return 
	 */
	void sort(char symbol_lt_gt = '<') {
		if (symbol_lt_gt == '<') {
			quick_sort_recur_lt(0, size() - 1);
		}
		else{
			quick_sort_recur_gt(0, size() - 1);
		}
	}
	void quick_sort_recur_lt(int start_ind, int end_ind) {
		if (start_ind >= end_ind)	return;
		int center = my::rand_int(start_ind, end_ind);
		switch_ele(start_ind, center);
		center = start_ind;
		int i = start_ind + 1;
		int j = end_ind;
		while (i != j) {
			if (e[i] <= e[center]) i++;
			else {
				switch_ele(i, j);
				j--;
			}
		}
		if (e[i] <= e[center]) {
			switch_ele(i, center);
			center = i;
		}
		else {
			switch_ele(i - 1, center);
			center = i - 1;
		}
		quick_sort_recur_lt(start_ind, center - 1);
		quick_sort_recur_lt(center + 1, end_ind);
	}
	void quick_sort_recur_gt(int start_ind, int end_ind) {
		if (start_ind >= end_ind)	return;
		int center = my::rand_int(start_ind, end_ind);
		switch_ele(start_ind, center);
		center = start_ind;
		int i = start_ind + 1;
		int j = end_ind;
		while (i != j) {
			if (e[i] >= e[center]) i++;
			else {
				switch_ele(i, j);
				j--;
			}
		}
		if (e[i] >= e[center]) {
			switch_ele(i, center);
			center = i;
		}
		else {
			switch_ele(i - 1, center);
			center = i - 1;
		}
		quick_sort_recur_gt(start_ind, center - 1);
		quick_sort_recur_gt(center + 1, end_ind);
	}

	/**
	 * .sort eleemnts from small to big, < as default, quick sort is the best choice
	 * 
	 * \return original index for the corresponding elements, the matrix can be recovered with this->permute_back();
	 */
	Matrix<int> sort_with_ind(char symbol_lt_gt = '<') {
		Matrix<int> result(1, size(), 'N');
		if (symbol_lt_gt == '<') {
			quick_sort_recur_lt_with_ind(0, size() - 1, result);
		}
		else {
			quick_sort_recur_gt_with_ind(0, size() - 1, result);
		}
		return result;
	}
	void quick_sort_recur_lt_with_ind(int start_ind, int end_ind, Matrix<int>& m_ind) {
		if (start_ind >= end_ind)	return;
		int center = my::rand_int(start_ind, end_ind);
		switch_ele(start_ind, center);
		m_ind.switch_ele(start_ind, center);
		center = start_ind;
		int i = start_ind + 1;
		int j = end_ind;
		while (i != j) {
			if (e[i] <= e[center]) i++;
			else {
				switch_ele(i, j);
				m_ind.switch_ele(i, j);
				j--;
			}
		}
		if (e[i] <= e[center]) {
			switch_ele(i, center);
			m_ind.switch_ele(i, center);
			center = i;
		}
		else {
			switch_ele(i - 1, center);
			m_ind.switch_ele(i - 1, center);
			center = i - 1;
		}
		quick_sort_recur_lt_with_ind(start_ind, center - 1, m_ind);
		quick_sort_recur_lt_with_ind(center + 1, end_ind, m_ind);
	}
	void quick_sort_recur_gt_with_ind(int start_ind, int end_ind, Matrix<int>& m_ind) {
		if (start_ind >= end_ind)	return;
		int center = my::rand_int(start_ind, end_ind);
		switch_ele(start_ind, center);
		m_ind.switch_ele(start_ind, center);
		center = start_ind;
		int i = start_ind + 1;
		int j = end_ind;
		while (i != j) {
			if (e[i] >= e[center]) i++;
			else {
				switch_ele(i, j);
				m_ind.switch_ele(i, j);
				j--;
			}
		}
		if (e[i] >= e[center]) {
			switch_ele(i, center);
			m_ind.switch_ele(i, center);
			center = i;
		}
		else {
			switch_ele(i - 1, center);
			m_ind.switch_ele(i - 1, center);
			center = i - 1;
		}
		quick_sort_recur_gt_with_ind(start_ind, center - 1, m_ind);
		quick_sort_recur_gt_with_ind(center + 1, end_ind, m_ind);
	}

	/**
	 * .before calling this, make sure both v1 and v2 is sorted from small to big
	 * this function only consider merge sort from small to big
	 * total_size and v1_start_ind is dedicated for Viterbi_optimized
	 */
	void merge_sort_lt(const Matrix<T>& v1, const Matrix<T>& v2, int total_size, int v1_start_ind = 0) {
		int v1_size = v1.size();
		int v2_size = v2.size();
		resize(1, total_size, false);
		int v1_ind = v1_start_ind;
		int v2_ind = 0;
		int i = 0;
		while (i < total_size && v1_ind < v1_size && v2_ind < v2_size) {
			if (v1(v1_ind) <= v2(v2_ind)) {
				e[i] = v1(v1_ind);
				v1_ind++;
			}
			else {
				e[i] = v2(v2_ind);
				v2_ind++;
			}
			++i;
		}
		if (v1_ind == v1_size) {
			while (i < total_size && v2_ind < v2_size) {
				e[i] = v2(v2_ind);
				v2_ind++;
				++i;
			}
		}
		else {		// v2_ind == v2_size
			while (i < total_size && v1_ind < v1_size) {
				e[i] = v1(v1_ind);
				v1_ind++;
				++i;
			}
		}
	}

	/**
	 * .binary search only for sorted Matrix, return the index of lower bound of upper bound of val, i.e., 
	 * {1,2,4,6,6,7} is the sorted Matrix, to search the lower bound
	 *	search 3, return 1, 
	 *	search 4, return 2, 
	 *	search 6, may return 3 or 4
	 *  search 8, return 5
	 * 
	 * \param val: the value of element
	 * \param upper_or_lower_bound: 'L' -> get lower bound's index, 'U' -> get upper bound's index
	 * \return the index closest to val, if multiple result, randomly return anyone
	 */
	int binary_search(T val, bool is_lower_bound = true) const {
		// to be finished
		// use the binary search method
		int start = 0;
		int end = size() - 1;

		//int mid = (start + end) >> 1;
		// binary search over v, O(log(end)) comparison and program jump
		while (start < end) {
			int mid = (start + end) >> 1;
			if (e[mid] < val) {
				start = mid + 1;
			}
			else {
				end = mid;
			}
		}
		if (is_lower_bound) {
			return e[start] <= val ? start : start - 1;
		}
		else {
			return e[start] >= val ? start : start + 1;
		}

	}

	/**
	 * .
	 *
	 * \return ret/2 = number of rows shifted back, ret%2 = 0 even valid permutation, = 1 odd valid permutation
	 */
	int row_transformation_to_up_triangle() {
		int permute_row_ind;
		int col_pos;
		int col_pos2;
		T temp;
		int ic, ir = 0;
		int r_shift = 0;
		int sign = 0;
		bool is_ordered_structure = !(T() < T());
		for (int i = 0; i + r_shift < r && i < c; i++) {
			//cout << "*this" << *this;
			ir = i + r_shift;
			ic = i;

			if (is_ordered_structure) {
				temp = max_abs_col_ele_start_from(ic, ir, permute_row_ind);	// for computation accuracy

				// Note: approach my::zero_approximation as 0
				if ((double) my::abs(temp) >= my::zero_approximation);							// this will affect the rank
				else {
					r_shift--;
					continue;
				}
			}
			else {
				temp = non_0_col_ele_start_from(ic, ir, permute_row_ind);	// specifically for GF

				if (temp != T(0));
				else {
					r_shift--;
					continue;
				}
			}

			//cout << "ir=" << ir << "\t" << "ic=" << ic << "\t" << "permute_row_ind=" << permute_row_ind << endl;
			if (ir != permute_row_ind) {
				switch_row(ir, permute_row_ind);
				sign = sign ^ 1;
			}

			col_pos = ir * c + ic;
			col_pos2 = col_pos;
			for (int j = ir + 1; j < r; j++) {
				col_pos2 += c;
				if (e[col_pos2] != T(0)) {		// in Hessenberg matrix this can skip, and faster
					temp = e[col_pos2] / e[col_pos];
					for (int k = ic; k < c; ++k) {
						(*this)(j, k) -= (*this)(ir, k) * temp;
					}
				}
			}
		}
		//cout << "*this" << *this;
		return sign + (-r_shift) * 2;
	}

	/**
	 * .only for matrix_polynomials
	 *
	 * \return ret/2 = number of rows shifted back, ret%2 = 0 even valid permutation, = 1 odd valid permutation
	 */
	int row_transformation_to_up_triangle_4_matrix_polynomials() {
		int permute_row_ind;
		int col_pos;
		int col_pos2;
		T temp;
		int ic, ir = 0;
		int r_shift = 0;
		int sign = 0;
		bool is_ordered_structure = !(T() < T());
		for (int i = 0; i + r_shift < r /* we consider the last line */ && i < c; i++) {
			//cout << "*this" << *this;
			ir = i + r_shift;
			ic = i;

			while (true) {
				//cout << "*this" << *this;

				if (is_ordered_structure) {
					// problem
					temp = max_col_ele_start_from_4_matrix_polynomials(ic, ir, permute_row_ind);		// for computation accuracy

					// Note: approach my::zero_approximation as 0
					if (!temp.is_approach_0());					// this will affect the rank
					else {
						r_shift--;		// this indicates that matrix singular, for determinant calcualtion, should return 0
						break;
					}
				}
				else {
					temp = non_0_col_ele_start_from(ic, ir, permute_row_ind);	// specifically for GF

					if (temp != T(0));
					else {
						r_shift--;
						break;
					}
				}

				//cout << "ir=" << ir << "\t" << "ic=" << ic << "\t" << "permute_row_ind=" << permute_row_ind << endl;
				if (ir != permute_row_ind) {
					switch_row(ir, permute_row_ind);
					sign = sign ^ 1;
				}

				//cout << "*this switch" << *this;


				col_pos = ir * c + ic;
				col_pos2 = col_pos;
				bool is_ir_below_zero = true;
				for (int j = ir + 1; j < r; j++) {
					col_pos2 += c;
					if (e[col_pos2] != T(0)) {		// in Hessenberg matrix this can skip, and faster

						//cout << "e[col_pos2]" << e[col_pos2];
						//cout << "e[col_pos]" << e[col_pos];
						//cout << " e[col_pos](e[col_pos].size() - 1) = " << e[col_pos](e[col_pos].size() - 1) << endl;

						temp = e[col_pos2] / e[col_pos];
						//cout << "temp" << temp;

						for (int k = ic; k < c; ++k) {
							/*cout << "k=" << k << endl;
							cout << "(*this)(" << ir << ", " << k << ")" << (*this)(ir, k);
							cout << "(*this)(" << ir << ", " << k << ") * temp" << (*this)(ir, k) * temp;
							cout << "(*this)(" << j << ", " << k << ")" << (*this)(j, k);
							cout << "minus" << (*this)(j, k) - (*this)(ir, k) * temp;*/

							(*this)(j, k) -= (*this)(ir, k) * temp;


						}
						is_ir_below_zero = is_ir_below_zero && (*this)(j, ic).is_approach_0();
					}
				}

				//cout << "*this end" << *this;

				if (!is_ir_below_zero);
				else {
					break;
				}
			}
		}
		//cout << "*this" << *this;
		return sign + (-r_shift) * 2;
	}

	/**
	 * .
	 * \return the recorded permutation
	 */
	Matrix<int> col_permute_to_full_rank_on_left() {
		Matrix<int> permute_record(1, c, 'N');
		bool is_ordered_structure = !(T() < T());

		// this step is after row transformation to up triangle
		if (!is_ordered_structure) {
			for (int i = 0; i < r; ++i) {
				for (int j = i; j < c; ++j) {	// remember to use strict equal only
					if ((*this)(i, j) != 0) {		// strictly equals to 0, okay with GF2
						if (i == j);
						else {
							switch_col(i, j);
							permute_record.switch_ele(i, j);
							//cout << "test=" << (*this);
						}
						break;
					}
				}
			}
		}
		else {
			for (int i = 0; i < r; ++i) {
				for (int j = i; j < c; ++j) {	// remember to use strict equal only
					if ((double)my::abs((*this)(i, j)) >= my::zero_approximation) {		// strictly equals to 0, okay with GF2
						if (i == j);
						else {
							switch_col(i, j);
							permute_record.switch_ele(i, j);
							//cout << "test=" << (*this);
						}
						break;
					}
				}
			}
		}
		return permute_record;
	}
	void row_transformation_left_up_triangle_to_identity() {
		// this function should be called after col_permute_to_full_rank_on_left
		for (int i = r - 1; i >= 0; --i) {	// diag ind
			// always (*this)(i, i) != 0
			if ((*this)(i, i) != T(1)) {			// scale row, not enter in GF2
				T factor = (*this)(i, i);
				for (int j = r; j < c; ++j) {	// col ind, knowing that col ind from i+1 to r-1 is 0 already
					(*this)(i, j) /= factor;
				}
				(*this)(i, i) = T(1);
			}
			for (int j = 0; j < i; ++j) {		// row ind
				if ((*this)(j, i) != T(0)) {		// in Tridiagonal matrix this can skip, and faster
					// row minus
					T factor = (*this)(j, i);
					for (int k = r; k < c; ++k) {	// col ind
						(*this)(j, k) -= factor * (*this)(i, k);
					}
					(*this)(j, i) = T(0);
				}
				//cout << "G in =" << (*this);
			}
		}
	}

	/**
	 * .
	 *
	 * \return ret/2 = number of rows shifted forward, ret%2 = 0 even valid permutation, = 1 odd valid permutation
	 */
	int row_transformation_to_low_triangle() {

		int permute_row_ind;
		int col_pos;
		int col_pos2;
		T temp;
		int ic = c - 1, ir = r - 1;
		int r_shift = 0;
		int sign = 0;
		bool is_ordered_structure = !(T() < T());
		for (; ir >= 0 && ic >= 0; ic--) {
			//cout << "*this" << *this;

			ir = ic - c + r + r_shift;

			if (is_ordered_structure) {
				temp = max_abs_col_ele_end_from(ic, ir, permute_row_ind);	// for computation accuracy

				// Note: approach my::zero_approximation as 0
				if ((double)my::abs(temp) >= my::zero_approximation);					// this will affect the rank
				else {
					r_shift++;
					continue;
				}
			}
			else {
				temp = non_0_col_ele_end_from(ic, ir, permute_row_ind);		// specifically for GF

				if (temp != T(0));
				else {
					r_shift++;
					continue;
				}
			}

			//cout << "ir=" << ir << "\t" << "ic=" << ic << "\t" << "permute_row_ind=" << permute_row_ind << endl;
			if (ir != permute_row_ind) {
				switch_row(ir, permute_row_ind);
				sign = sign ^ 1;
			}

			col_pos = ir * c + ic;
			col_pos2 = col_pos;
			for (int j = ir - 1; j >= 0; j--) {
				col_pos2 -= c;
				if (e[col_pos2] != 0) {		// in Hessenberg matrix this can skip, and faster
					temp = e[col_pos2] / e[col_pos];
					for (int k = ic; k >= 0; --k) {
						(*this)(j, k) -= (*this)(ir, k) * temp;
					}
				}
			}
		}
		//cout << "*this" << *this;
		return sign + r_shift * 2;
	}

	/**
	 * .
	 * \return the record of permute
	 */
	Matrix<int> col_permute_to_full_rank_on_right() {
		// this step is after row transformation to low triangle

		Matrix<int> permute_record(1, c, 'N');
		bool is_ordered_structure = !(T() < T());

		if (!is_ordered_structure) {
			for (int i = r - 1; i >= 0; --i) {
				for (int j = i - r + c; j >= 0; --j) { 	// remember to use strict equal only
					if ((*this)(i, j) != 0) {	// strictly equals to 0, okay with GF2
						if (i == j - c + r);
						else {
							switch_col(j, i - r + c);
							permute_record.switch_ele(j, i - r + c);
							//cout << "test=" << (*this);
						}
						break;
					}
				}
			}
		}
		else {
			for (int i = r - 1; i >= 0; --i) {
				for (int j = i - r + c; j >= 0; --j) {	// remember to use strict equal only
					if ((double)my::abs((*this)(i, j)) >= my::zero_approximation) {		// strictly equals to 0, okay with GF2
						if (i == j - c + r);
						else {
							switch_col(j, i - r + c);
							permute_record.switch_ele(j, i - r + c);
							//cout << "test=" << (*this);
						}
						break;
					}
				}
			}
		}

		return permute_record;
	}
	void row_transformation_right_low_triangle_to_identity() {
		// this function should be called after col_permute_to_full_rank_on_left
		for (int i = 0; i < r; ++i) {	// row ind
			int col_ind = c - r + i;
			// always (*this)(i, col_ind) != 0
			if ((*this)(i, col_ind) != 1) {
				T factor = (*this)(i, col_ind);
				for (int j = c - r - 1; j >= 0; --j) {	// col ind, knowing that col ind from c-r to i-1 is 0 already
					(*this)(i, j) /= factor;
				}
				(*this)(i, col_ind) = 1;
			}
			for (int j = r-1; j > i; --j) {		// row ind
				if ((*this)(j, col_ind) != 0) {		// in Tridiagonal matrix this can skip, and faster
					// row minus
					T factor = (*this)(j, col_ind);
					for (int k = c - r - 1; k >= 0; --k) {	// col ind
						(*this)(j, k) -= factor * (*this)(i, k);
					}
					(*this)(j, col_ind) = 0;
				}
				//cout << "G in =" << (*this);
			}
		}
	}

	/**
	 * .		Gaussian Elimination to make identity as left as possible
	 *			we donot need to consider row permutation record, or singular case
	 *			(*this) should be fat matrix and full rank of rows
	 * 
	 * \param	perform GE to have echelon form, this save only 1% of the complexity, time cost is almost same
	 *			doesn't make a difference, discarding, write it just to be in correspondence with the algorithm description in the book
	 * 
	 * \return columns' permutation record
	 */
	Matrix<int> GE_left_identity_4_GF2_echelon() {		
		
		// record the permutation
		Matrix<int> ans(1, r, 'v');

		//cout << "column proceed = " << 0 << ", density() = " << density() << endl;

		// row transformation to turn the matrix into an up-triangle Matrix
		for (int row_now = 0, column_now = 0; row_now < r && column_now < c; ++row_now, ++column_now) {

			if ((*this)(row_now, column_now) == 0) {
				// find the row below row_now that have entry 1
				int nz_row_below;
				for (nz_row_below = row_now + 1; nz_row_below < r && (*this)(nz_row_below, column_now) == 0; ++nz_row_below);

				if (nz_row_below != r) {
					// in this case we can permute nz_row_below with row_now to get entry 1 in position (row_now, column_now)
					switch_row(nz_row_below, row_now);
				}
				else {
					// in this case the column column_now is all zero and we have to seek another column
					row_now--;			// keep row_now unchanged and move to the next column
					continue;
				}
			}

			// here, we can use row_now to eleminate other rows of the Matrix		
			ans.push_back(column_now);
			//cout << "ans" << ans;

			for (int i = row_now + 1; i < r; ++i) {
				if ((*this)(i, column_now) == 1) {

					(*this)(i, column_now) = 0;		// this should be count as a xor operation !!
					for (int j = column_now + 1; j < c; ++j) {
						(*this)(i, j) += (*this)(row_now, j);
					}
				}
			}

			// check out the density change after the row elimination
			//cout << "column proceed = " << column_now + 1 << ", density() = " << density(row_now + 1, column_now + 1) << endl;
		}
		//cout << "*this" << *this;

		// row transformation to generate an echelon form matrix
		int next_column = ans.back();
		int start_col_to_eliminate = next_column + 1;
		Matrix<int> more_to_eliminate(1, start_col_to_eliminate - r, 'v');		// the size is accurate.

		for (int ans_ind = ans.size() - 1; ans_ind > 0; --ans_ind) {
			int column_now = next_column;

			// eliminate the rows above ans_ind
			for (int i = 0; i < ans_ind; ++i) {
				if ((*this)(i, column_now) == 1) {
					(*this)(i, column_now) = 0;
					// eliminate the cols start from start_col_to_eliminate to the end
					for (int j = start_col_to_eliminate; j < c; ++j) {
						(*this)(i, j) += (*this)(ans_ind, j);
					}
					for (int j = 0, jmax = more_to_eliminate.size(); j < jmax; ++j) {
						(*this)(i, more_to_eliminate(j)) += (*this)(ans_ind, more_to_eliminate(j));	
							// for less operation, not for time saving, only less than 1% of operation can be saved!
					}
				}
			}

			next_column = ans(ans_ind - 1);
			for (int i = column_now - 1; i > next_column; --i) {
				more_to_eliminate.push_back(i);
			}
		}

		// ans should be in echelon form		

		// column permute the Matrix to be full rank on left

		// we do not have decreasing reliability shortly after k, this does not affect optimal condition for OSD, be careful
		Matrix<int> permute_record(1, c, 'N');
		for (int i = 0; i < r; ++i) {
			int j;
			for (j = i; j < c && (*this)(i, j) == 0; ++j);	// remember to use strict equal only
			
			// strictly equals to 0, okay with GF2

			if (i == j);		// no switch is needed
			else{
				switch_col(i, j);
				permute_record.switch_ele(i, j);
			}
		}

		// row transformation to turn the matrix's left part to identity

		return permute_record;
	}

	/**
	 * .perform Gaussian Elimination to make identity as left as possible
	 * we donot need to consider row permutation record, or singular case
	 * (*this) should be fat matrix and full rank of rows
	 *
	 * \return columns' permutation record
	 */
	Matrix<int> GE_left_identity_4_GF2() {

		// row transformation to turn the matrix into an up-triangle Matrix
		for (int row_now = 0, column_now = 0; row_now < r && column_now < c; ++row_now, ++column_now) {

			if ((*this)(row_now, column_now) == 0) {
				// find the row below row_now that have entry 1
				int nz_row_below;
				for (nz_row_below = row_now + 1; nz_row_below < r && (*this)(nz_row_below, column_now) == 0; ++nz_row_below);

				if (nz_row_below != r) {
					// in this case we can permute nz_row_below with row_now to get entry 1 in position (row_now, column_now)
					switch_row(nz_row_below, row_now);
				}
				else {
					// in this case the column column_now is all zero and we have to seek another column
					row_now--;			// keep row_now unchanged and move to the next column
					continue;
				}
			}

			// here, we can use row_now to eleminate other rows of the Matrix		

			for (int i = row_now + 1; i < r; ++i) {
				if ((*this)(i, column_now) == 1) {

					(*this)(i, column_now) = 0;
					for (int j = column_now + 1; j < c; ++j) {
						(*this)(i, j) += (*this)(row_now, j);
					}
				}
			}
		}

		// column permute the Matrix to be full rank on left

		// we do not have decreasing reliability shortly after k, this does affect optimal condition for OSD, be careful
		Matrix<int> permute_record(1, c, 'N');
		for (int i = 0; i < r; ++i) {
			int j;
			for (j = i; j < c && (*this)(i, j) == 0; ++j);	// remember to use strict equal only

			// strictly equals to 0, okay with GF2

			if (i == j);		// no switch is needed
			else {
				switch_col(i, j);
				permute_record.switch_ele(i, j);
			}
		}

		// row transformation to turn the matrix's left part to identity

		// start from the column r-1
		for (int column_now = r - 1; column_now >= 0; --column_now){	
			for (int i = column_now - 1; i >= 0; --i) {
				if ((*this)(i, column_now) == 1) {

					(*this)(i, column_now) = 0;
					for (int j = r; j < c; ++j) {
						// we only need to eleminate the right part of the matrix
						(*this)(i, j) += (*this)(column_now, j);
					}
				}				
			}
		}

		return permute_record;
	}

	/**
	 * .perform Gaussian Elimination to make identity as left as possible,
	 * take in the last weight-1 column with 1 in the last row, as a part of identity
	 * we donot need to consider row permutation record, or singular case
	 * (*this) should be fat matrix and full rank of rows
	 *
	 * .Be careful to use this function, make sure the last column is independent to all the other columns
	 * 
	 * \return columns' permutation record
	 */
	Matrix<int> GE_left_identity_4_GF2_with_last_w1_col() {
		// please never use that function. Have a problem

		// row transformation to turn the matrix into an up-triangle Matrix
		for (int row_now = 0, column_now = 0; row_now < r && column_now < c; ++row_now, ++column_now) {

			if (entry(row_now, column_now) == 0) {
				// find the row below row_now that have entry 1
				int nz_row_below;
				for (nz_row_below = row_now + 1; nz_row_below < r && entry(nz_row_below, column_now) == 0; ++nz_row_below);

				if (nz_row_below != r) {
					// in this case we can permute nz_row_below with row_now to get entry 1 in position (row_now, column_now)
					switch_row(nz_row_below, row_now);
				}
				else {
					// in this case the column column_now is all zero and we have to seek another column
					row_now--;			// keep row_now unchanged and move to the next column
					continue;
				}
			}

			// here, we can use row_now to eleminate other rows of the Matrix		

			for (int i = row_now + 1; i < r; ++i) {
				if (entry(i, column_now) == 1) {

					entry(i, column_now) = 0;
					for (int j = column_now + 1; j < c; ++j) {
						entry(i, j) += entry(row_now, j);
					}
				}
			}
		}

		// column permute the Matrix to be full rank on left

		// we do not have decreasing reliability shortly after k, this does not affect optimal condition for OSD, be careful
		Matrix<int> permute_record(1, c, 'N');
		for (int i = 0, imax = r - 1; i < imax; ++i) {		// here we consider the first r-1 columns and the last column, which is different
			int j;
			for (j = i; j < c && entry(i, j) == 0; ++j);	// remember to use strict equal only

			// strictly equals to 0, okay with GF2

			if (i == j);		// no switch is needed
			else {
				switch_col(i, j);
				permute_record.switch_ele(i, j);
			}
		}

		// permute the last column to the position r-1
		switch_col(c - 1, r - 1);
		permute_record.switch_ele(c - 1, r - 1);

		//cout << "(*this) (middle)" << (*this) << endl;

		// row transformation to turn the matrix's left part to identity

		// start from the column r-1, this is different, remember that
		for (int column_now = r - 1; column_now >= 0; --column_now) {
			for (int i = column_now - 1; i >= 0; --i) {
				if (entry(i, column_now) == 1) {

					entry(i, column_now) = 0;
					for (int j = r; j < c; ++j) {
						// we only need to eleminate the right part of the matrix
						entry(i, j) += entry(column_now, j);
					}
				}
			}
		}

		return permute_record;
	}

	/**
	 * .perform Gaussian-Jordan Elimination to make identity as left as possible
	 * we donot need to consider row permutation record, or singular case
	 * (*this) should be fat matrix and full rank of rows
	 *
	 * \param num_of_column_processed: record number of columns traversed, which is the number of sequential steps
	 * 
	 * \return columns' permutation record
	 */
	Matrix<int> GJE_left_identity_4_GF2(int& num_of_column_processed) {

		// row transformation to turn the matrix into an up-triangle Matrix
		int column_now = 0;
		for (int row_now = 0; row_now < r && column_now < c; ++row_now, ++column_now) {

			if ((*this)(row_now, column_now) == 0) {
				// find the row below row_now that have entry 1
				int nz_row_below;
				for (nz_row_below = row_now + 1; nz_row_below < r && (*this)(nz_row_below, column_now) == 0; ++nz_row_below);

				if (nz_row_below != r) {
					// in this case we can permute nz_row_below with row_now to get entry 1 in position (row_now, column_now)
					switch_row(nz_row_below, row_now);
				}
				else {
					// in this case the column column_now is all zero and we have to seek another column
					row_now--;			// keep row_now unchanged and move to the next column
					continue;
				}
			}

			// here, we can use row_now to eleminate other rows of the Matrix		

			// eliminate rows above and below the pivot row
			for (int i = 0; i < r; ++i) {
				if (i != row_now && (*this)(i, column_now) == 1) {

					for (int j = column_now; j < c; ++j) {
						(*this)(i, j) += (*this)(row_now, j);
					}
				}
			}
		}
		num_of_column_processed = column_now;

		// column permute the Matrix to be full rank on left

		// we do not have decreasing reliability shortly after k, this does affect optimal condition for OSD, be careful
		Matrix<int> permute_record(1, c, 'N');
		for (int i = 0; i < r; ++i) {
			int j;
			for (j = i; j < c && (*this)(i, j) == 0; ++j);	// remember to use strict equal only

			// strictly equals to 0, okay with GF2

			if (i == j);		// no switch is needed
			else {
				switch_col(i, j);
				permute_record.switch_ele(i, j);
			}
		}

		return permute_record;
	}

	/**
	 * .count the density of matrix of the right down corner, i.e., (row_start_ind to end, col_start_ind to end)
	 * 
	 * \param row_start_ind
	 * \param col_start_ind
	 * \return the density
	 */
	double density(int row_start_ind = 0, int col_start_ind = 0) {
		int num_of_1 = 0;
		for (int i = row_start_ind; i < r; ++i) {
			for (int j = col_start_ind; j < c; ++j) {
				if (entry(i, j) != 0) {
					num_of_1++;
				}
			}
		}
		return double(num_of_1) / (r - row_start_ind) / (c - col_start_ind);
	}

	/**
	 * .		shift the columns of matrix circularly to the right
	 *			maybe we can consider pass a reference as parameter, reducing memory declaration and destruction
	 * 
	 * \param	shift_value: the number of columns to shift, positive stands for shift right, negtive stands for shift left 
	 *				e.g., (*this) is [a,b,c,d], and shift_value is 1, the output matrix is [d,a,b,c]
	 *				e.g., (*this) is [a,b,c,d], and shift_value is -1, the output matrix is [b,c,d,a]
	 *			
	 *			ans: the returned matrix
	 * 
	 */
	void col_shift_right_cir(int shift_value, Matrix<T>& ans) const{
		if (&ans != this) {
			ans.resize(r, c, false);
			int col_start_ind = 0;
			shift_value = shift_value % c;		// restrict shift_value in (-c, c) preventing overflow of col_ind

			if (shift_value >= 0) {		// shift right circularly
				for (int i = 0; i < r; ++i) {

					// copy the columns before shift_value
					for (int j_ans = col_start_ind, j_ans_max = col_start_ind + shift_value, j_this = col_start_ind + c - shift_value; \
						j_ans < j_ans_max; ++j_ans, ++j_this) {

						ans.e[j_ans] = e[j_this];		// I believe this is the best program without parallelization
						// never use c function memcpy in c++
					}

					// copy the columns after shift_value
					for (int j_ans = shift_value + col_start_ind, j_ans_max = c + col_start_ind, j_this = col_start_ind; \
						j_ans < j_ans_max; ++j_ans, ++j_this) {

						ans.e[j_ans] = e[j_this];
					}

					col_start_ind += c;
				}
			}
			else {					// shift left circularly
				shift_value = -shift_value;
				for (int i = 0; i < r; ++i) {

					// copy the columns before shift_value
					for (int j_ans = col_start_ind + c - shift_value, j_ans_max = col_start_ind + c, j_this = col_start_ind; \
						j_ans < j_ans_max; ++j_ans, ++j_this) {

						ans.e[j_ans] = e[j_this];
					}

					// copy the columns after shift_value
					for (int j_ans = col_start_ind, j_this = col_start_ind + shift_value, j_this_max = col_start_ind + c; \
						j_this < j_this_max; ++j_ans, ++j_this) {

						ans.e[j_ans] = e[j_this];
					}

					col_start_ind += c;
				}
			}
		}
		else {		// it is strange that by setting &ans == this, the function changed (*this), but const is a property of this function
			cout << "(col_shift_right_cir) warning: please don't use (*this) as ans, this will make this function not const" << endl;
		}
	}

	void col_ind_mul_cir(int mul_num, Matrix<T>& ans) const {
		if (&ans != this) {
			ans.resize(r, c, false);
			int col_start_ind = 0;
			
			for (int i = 0; i < r; ++i) {
				for (int j = 0; j < c; ++j) {
					int col_ind = j * mul_num;
					col_ind = col_ind % c;

					ans.e[col_start_ind + col_ind] = e[col_start_ind + j];
				}
				col_start_ind += c;
			}
		}
		else {		// it is strange that by setting &ans == this, the function changed (*this), but const is a property of this function
			cout << "(col_ind_mul_cir) warning: please don't use (*this) as ans, this will make this function not const" << endl;
		}
	}

	// following functions do not change the matrix

	// new function for polynomials, principle: eleminate friend, the functions should be simplified, keep less interface to avoide confusion

	// with num>0 shift right and num<0 shift left
	void shift_right(int num) {
		if (r == 1) {		// r == 0 is no longer allowed, but c == 0 is allowed
			if (num > 0) {
				int c_new = c + num;
				if (capacity < c_new) {
					T* temp = e;
					e = new T[c_new];

					for (int i = c_new - 1; i >= num; --i) {
						e[i] = temp[i - num];
					}
					for (int i = num - 1; i >= 0; --i) {
						e[i] = T(0);
					}
					delete[] temp;
					c = c_new;
					capacity = c_new;		// capacity never shrink
				}
				else {
					for (int i = c_new - 1; i >= num; --i) {
						e[i] = e[i - num];
					}
					for (int i = num - 1; i >= 0; --i) {
						e[i] = T(0);
					}
					c = c_new;
				}
			}
			else if (num == 0);
			else {
				for (int i = 0; i < c; ++i) {
					e[i] = e[i - num];
				}
				c += num;
				//get_part(0, -num, 0, c - 1);
			}
		}
		else {
			cout << "size error, shift failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
		}
	}
	Matrix<T> get_shift_right(int num) const{
		if (r == 1) {		// r == 0 is no longer allowed, but c == 0 is allowed
			if (num > 0) {
				Matrix<T> result(1, c + num);
				for (int i = 0; i < num; ++i) {
					result.e[i] = 0;
				}
				for (int i = 0; i < c; ++i) {
					result.e[num + i] = e[i];
				}
				return result;
			}
			else if (num == 0) {
				return *this;
			}
			else {
				return get_part(0, -num, 0, c - 1);
			}
		}
		else {
			cout << "size error, shift failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return Matrix<T>();
		}
	}
	/**
	 * .
	 * 
	 * \return index of end non 0, ret >=0, if array all 0, set c=1 and return 0
	 */
	int cut_end_0() {		// return same as 'ind_of_end_non_0()'
		if (r == 1) {		// the decleared size is not changed and space can be used in the future
			int i = c - 1;
			for (; i > 0; i--)
				if (e[i] == 0);		// will not judge e[0]
				else break;

			// prevent moving the array, since it shrinks

			c = i + 1;		// c is at least 1, if array all 0, c=1
			return i;		// return value is at least 0
		}
		else {
			cout << "size error, cut_end_0 failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return 0;
		}
	}
	/**
	 * .
	 * 
	 * \return >=0, if array all 0, set c=1 and return 0
	 */
	int ind_of_end_non_0() const {
		if (r == 1) {
			int i = c - 1;
			for (; i > 0; i--)
				if (e[i] == 0);		// will not judge e[0]
				else break;
			return i;				// return value is at least 0
		}
		else {
			cout << "size error, ind_of_end_non_0 failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return 0;
		}
	}
	/**
	 * .
	 *
	 * \return index of end > my::zero_approximation for ordered structure, ret >=0, if array all 0, set c=1 and return 0
	 */
	int cut_end_small() {		// return similar as 'ind_of_end_non_0()'
		if (r == 1) {		// the decleared size is not changed and space can be used in the future
			int i = c - 1;
			for (; i > 0; i--)
				if ((double)my::abs(e[i]) < my::zero_approximation);		// will not judge e[0]
				else break;

			// prevent moving the array, since it shrinks

			c = i + 1;		// c is at least 1, if array all 0, c=1
			return i;		// return value is at least 0
		}
		else {
			cout << "size error, cut_end_0 failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return 0;
		}
	}

	//compute functions(return a new matrix if can)
	T norm_1() const {
		// only valid for vector, return abs sum of each element
		int size_all = size();
		T ans = T(0);
		for (int i = 0; i < size_all; ++i) {
			ans += my::abs(e[i]);
		}
		return ans;
	}
	T norm_2() const {
		if (isVector()) {
			int size_all = size();
			T sum_abs_2 = T(0);
			for (int i = 0; i < size_all; i++)
				sum_abs_2 += my::conj(e[i]) * e[i];
			return sqrt(sum_abs_2);
		}
		else {
			/*Matrix<T> AH_A = Hermit() * (*this);
			Matrix<T> eig_values = AH_A.eig_val();
			return sqrt(eig_values(0));*/
			return T(0);			// to be solved
		}
	}
	T det() const {
		if (r == c) {
			Matrix<T> U(*this);
			int r_shift = U.row_transformation_to_up_triangle();
			if (r_shift / 2 == 0) {
				int diag_pos = 0;
				T ans = T(1);
				for (int i = 0; i < r; i++) {
					ans *= U(diag_pos);
					diag_pos += (c + 1);
				}
				return (r_shift % 2 == 0 ? ans : -ans);
			}
			else {
				return T(0);
			}
		}
		else {
			cout << "size error, [det] failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return T(0);
		}
	}
	T det_4_matrix_polynomials() const {
		// designed for matrix polynomials, which is very powerful
		if (r == c) {
			Matrix<T> U(*this);
			int r_shift = U.row_transformation_to_up_triangle_4_matrix_polynomials();
			//cout << "U" << U;
			if (r_shift / 2 == 0) {
				int diag_pos = 0;
				T ans = T(1);
				for (int i = 0; i < r; i++) {
					ans *= U(diag_pos);
					diag_pos += (c + 1);
				}
				return (r_shift % 2 == 0 ? ans : -ans);
			}
			else {
				return T(0);
			}
		}
		else {
			cout << "size error, [det_4_matrix_polynomials] failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return T(0);
		}
	}
	int rank() const {
		Matrix<T> U(*this);
		U.row_transformation_to_up_triangle();
		//cout << "U" << U;
		int defficient_num = 0;
		bool is_ordered_structure = !(T() < T());

		if (is_ordered_structure) {
			for (int i = (r < c ? r : c) - 1; i >= 0; --i) {
				for (int j = i; j < c; ++j) {
					if ((double) my::abs(U(i, j)) >= my::zero_approximation) {
						return i + 1;
					}
				}
			}
		}
		else {
			for (int i = (r < c ? r : c) - 1; i >= 0; --i) {
				for (int j = i; j < c; ++j) {
					if (U(i, j) != 0) {
						return i + 1;
					}
				}
			}
		}
		return 0;
		//return (U.c > U.r ? U.r : U.c) - (r_shift / 2);
	}
	T condition_num() const {
		int rk = rank();
		if (rk == 0) {
			return (T)NAN;
		}
		if (rk < r && rk < c) {
			return (T)INFINITY;
		}
		Matrix<T> AH_A = Hermit() * (*this);
		Matrix<T> eig_values = AH_A.eig_val();
		return sqrt(eig_values(0)) / sqrt(eig_values(c - 1));	// max singular value over min singular value
	}
	T norm_F() const {
		int sss = size();
		T result = 0;
		for (int i = 0; i < sss; ++i) {
			result += my::conj(e[i]) * e[i];
		}
		return sqrt(result);
	}

	// if A and B are vectors, return A^H * B, if A and B are matrices, return trace of A'*B or A*B'
	T dot_product(const Matrix<T>& B) const {
		int size_all = size();
		if (isVector() && B.isVector() && size_all == B.size()) {
			// case of vector dot product
			T ans = T(0);
			for (int i = 0; i < size_all; i++)
				ans += my::conj(e[i]) * B.e[i];
			return ans;
		}
		else if (r == B.r) {
			// case of matrix dot product
			return (Hermit() * B).trace();
		}
		else if (c == B.c) {
			return ((*this) * B.Hermit()).trace();
		}
		else {
			cout << "size error, dot_product failed" << endl;
			cout << "r=" << row() << "c=" << col() << endl;
			cout << ",\tB.row()=" << B.row() << ",\tB.col()=" << B.col() << endl;
			return T(0);
		}
	}
	// convolution
	Matrix<T> conv(const Matrix<T>& B) const {
		int len1 = size();
		int len2 = B.size();
		Matrix<T> result(1, len1 + len2 - 1, '0');
		for (int j = 0; j < len2; ++j) {
			if (B(j) != 0) {
				for (int i = 0; i < len1; ++i) {
					result(i + j) += (*this)(i) * B(j);
				}
			}
		}
		return result;		// may contain ending 0
	}
	Matrix<T> inv() const {
		if (r == c) {
			Matrix<T> ext = combine_right(Matrix<T>(r, c, 'i'));
			//cout << "ext(1)" << ext;

			/* row_transformation_to_up_triangle */
			int r_shift = ext.row_transformation_to_up_triangle();
			//cout << "ext(2)" << ext;

			if (r_shift / 2 == 0);
			else {
				cout << "Matrix singular, inv error" << endl;
				cout << "*this" << *this;
				return ext.get_part(0, c, r - 1, ext.c - 1);
			}

			/* row transformation left up triangle to identity */
			ext.row_transformation_left_up_triangle_to_identity();
			//cout << "ext(3)" << ext;
			return ext.get_part(0, c, r - 1, ext.c - 1);
		}
		else {
			Matrix<T> H = Hermit();
			//will perform pseudo inverse
			if (r > c) {
				// find left inverse
				return (H * (*this)).inv() * H;
			}
			else {
				// find right inverse
				return H * ((*this) * H).inv();
			}
		}
	}
	// use Gaussian elemination to compute left multiply of inv(), same as \ operator in Matlab, i.e., return (*this)\B
	Matrix<T> inv_then_multiply(const Matrix<T>& B) const {
		if (r == c) {
			Matrix<T> ext = combine_right(B);
			//cout << "ext" << ext;

			/* row_transformation_to_up_triangle */
			int r_shift = ext.row_transformation_to_up_triangle();
			//cout << "ext" << ext;
			if (r_shift / 2 == 0);
			else {
				cout << "Matrix singular, inv error" << endl;
				cout << "*this" << *this;
				return ext.get_part(0, c, r - 1, ext.c - 1);
			}

			/* row transformation left up triangle to identity */
			ext.row_transformation_left_up_triangle_to_identity();
			//cout << "ext" << ext;
			return ext.get_part(0, c, r - 1, ext.c - 1);
		}
		else {
			Matrix<T> H = Hermit();
			//will perform pseudo inverse
			if (r > c) {
				// find left inverse
				return (H * (*this)).inv_then_multiply(H)* B;
			}
			else {
				// find right inverse
				return H * ((*this) * H).inv_then_multiply(B);
			}
		}
	}

	// return trace of a matrix
	T trace() const {
		int min_rc = r < c ? r : c;
		T ans = 0;
		for (int i = 0; i < min_rc; ++i) {
			ans += (*this)(i, i);
		}
		return ans;
	}

	// Schimidt orthogonalization
	Matrix<T> Schimidt_orthogonalization(bool by_column = true, \
		double orthogonal_precision = my::zero_approximation * my::zero_approximation) const {
		// if row equals column, and not deficient, this return an orthogonal matrix and it doesn't matter to deal with row or column vector

		// but if deficient, deal with column vector as default and represent complementary of column   space as zero column vector

		// consider throwing away zero vector? then output may have size less than original matrix, to be done, but i think it problematic

		Matrix<T> beta = (*this) / max_abs_ele();		// first normalize the max element of matrix be 1

		//cout << "max_abs_ele() = " << max_abs_ele() << endl;
		//cout << "beta" << beta;									// problem in dealing with complex numbers


		if (r > c || (r == c && by_column)) {

			for (int i = 0; i < c; ++i) {
				for (int j = 0; j < i; ++j) {
					// compute sum[j](*this)(col j)^T*(*this)(col i)
					T tmp = 0;
					for (int k = 0; k < r; ++k) {
						tmp += my::conj(beta(k, j)) * beta(k, i);
					}
					for (int k = 0; k < r; ++k) {
						beta(k, i) -= tmp * beta(k, j);
					}
				}
				// compute column norm
				T col_norm2 = 0;
				for (int j = 0; j < r; ++j) {
					col_norm2 += my::conj(beta(j, i)) * beta(j, i);
				}

				//cout << "col_norm2 = " << col_norm2 << endl;

				if (col_norm2 >= orthogonal_precision) {		// square of my::zero_approximation
					col_norm2 = sqrt(col_norm2);	// now be norm2
					for (int j = 0; j < r; ++j) {
						beta(j, i) /= col_norm2;
					}
				}
				else {
					//cout << "column: " << i << " is zero" << endl;
					//cout << "col_norm2=" << col_norm2 << endl;
					for (int j = 0; j < r; ++j) {
						beta(j, i) = 0;		// denote it as zero vector, meaning rank deficient
					}
				}
			}
			return beta;
		}
		else {

			for (int i = 0; i < r; ++i) {
				for (int j = 0; j < i; ++j) {
					// compute sum[j](*this)(row j)^T*(*this)(row i)
					T tmp = 0;
					for (int k = 0; k < c; ++k) {
						tmp += my::conj(beta(j, k)) * beta(i, k);
					}
					for (int k = 0; k < c; ++k) {
						beta(i, k) -= tmp * beta(j, k);
					}
				}

				// compute row norm
				T row_norm2 = 0;

				row_norm2 = 0;
				for (int j = 0; j < c; ++j) {
					row_norm2 += my::conj(beta(i, j)) * beta(i, j);
				}

				//cout << "row_norm2 = " << row_norm2 << endl;

				if (row_norm2 >= orthogonal_precision) {
					row_norm2 = sqrt(row_norm2);
					for (int j = 0; j < c; ++j) {
						beta(i, j) /= row_norm2;
					}
				}
				else {
					for (int j = 0; j < c; ++j) {
						beta(i, j) = 0;		// denote it as zero vector, meaning rank deficient
					}
				}
			}

			return beta;
		}
	}
	// find orthogonal complementary with only one 1 in row/col vector, this may add nothing
	Matrix<T> add_w1_orthogonal_complementary() const {
		// before call this, make sure (*this) is orthogonal, you can call Schimidt_orthogonalization() to do that
		if (r > c) {
			// recording weight 1 ind
			Matrix<int> w1_ind(1, r);
			int w1_actual_size = 0;
			for (int i = 0; i < r; ++i) {
				T row_abs_sum = 0;
				for (int j = 0; j < c && (double) row_abs_sum < my::zero_approximation; ++j) {
					row_abs_sum += my::abs((*this)(i, j));
				}
				if ((double) row_abs_sum >= my::zero_approximation);		// not orthogonal
				else {
					// orthogonal
					w1_ind(w1_actual_size) = i;
					w1_actual_size++;		// just like push_back
				}
			}
			if (w1_actual_size == 0) {
				return (*this);
			}
			else {
				// add those weigth 1 vector back into the original matrix
				Matrix<T> add_back_vec(r, w1_actual_size, '0');
				for (int i = 0; i < w1_actual_size; ++i) {
					add_back_vec(w1_ind(i), i) = 1;
				}
				return combine_right(add_back_vec);
			}
		}
		else if (r < c) {
			// recording weight 1 ind
			Matrix<int> w1_ind(1, c);
			int w1_actual_size = 0;
			for (int j = 0; j < c; ++j) {
				T col_abs_sum = 0;
				for (int i = 0; i < r && (double) col_abs_sum < my::zero_approximation; ++i) {
					col_abs_sum += my::abs((*this)(i, j));		// this is stricter than row-by-row orthogonal
				}
				if ((double) col_abs_sum >= my::zero_approximation);		// not orthogonal
				else {
					// orthogonal
					w1_ind(w1_actual_size) = j;
					w1_actual_size++;		// just like push_back
				}
			}
			if (w1_actual_size == 0) {
				return (*this);
			}
			else {
				// add those weigth 1 vector back into the original matrix
				Matrix<T> add_back_vec(w1_actual_size, c, '0');
				for (int i = 0; i < w1_actual_size; ++i) {
					add_back_vec(i, w1_ind(i)) = 1;
				}
				return combine_down(add_back_vec);
			}
		}
		else {
			return (*this);
		}
	}
	// start with a random vector, use Schimidt orthogonalization to get orthogonal complementary vector
	Matrix<T> add_rand_orthogonal_complementary() const {
		// before call this, make sure (*this) is orthogonal, you can call Schimidt_orthogonalization() to do that
		if (r > c) {
			Matrix<T> beta(r, r);
			for (int i = 0; i < r; ++i) {
				for (int j = 0; j < c; ++j) {
					beta(i, j) = (*this)(i, j);			
				}
				for (int j = c; j < r; ++j) {
					beta(i, j) = my::rand_u();	// set it random, then use method of Schimidt orthogonalization
				}
			}

			// compute column norm
			T col_norm2 = 0;

			for (int i = c; i < r; ++i) {
				for (int j = 0; j < i; ++j) {
					// compute sum[j](*this)(col j)^H*(*this)(col i)
					T tmp = 0;
					for (int k = 0; k < r; ++k) {
						tmp += my::conj(beta(k, j)) * beta(k, i);
					}
					for (int k = 0; k < r; ++k) {
						beta(k, i) -= tmp * beta(k, j);
					}
				}
				col_norm2 = 0;
				for (int j = 0; j < r; ++j) {
					col_norm2 += my::conj(beta(j, i)) * beta(j, i);
				}
				col_norm2 = sqrt(col_norm2);
				for (int j = 0; j < r; ++j) {
					beta(j, i) /= col_norm2;
				}
			}

			return beta;
		}
		else if (r < c) {
			Matrix<T> beta(c, c);
			for (int i = 0; i < r; ++i) {
				for (int j = 0; j < c; ++j) {
					beta(i, j) = (*this)(i, j);
				}
			}
			for (int i = r; i < c; ++i) {
				for (int j = 0; j < c; ++j) {
					beta(i, j) = my::rand_u();	// set it random, then use method of Schimidt orthogonalization
				}
			}
			// compute column norm
			T row_norm2 = 0;

			for (int i = r; i < c; ++i) {
				for (int j = 0; j < i; ++j) {
					// compute sum[j](*this)(col j)^H*(*this)(col i)
					T tmp = 0;
					for (int k = 0; k < c; ++k) {
						tmp += my::conj(beta(j, k)) * beta(i, k);
					}
					for (int k = 0; k < c; ++k) {
						beta(i, k) -= tmp * beta(j, k);
					}
				}
				row_norm2 = 0;
				for (int j = 0; j < c; ++j) {
					row_norm2 += my::conj(beta(i, j)) * beta(i, j);
				}
				row_norm2 = sqrt(row_norm2);
				for (int j = 0; j < c; ++j) {
					beta(i, j) /= row_norm2;
				}
			}

			return beta;
		}
		else {
			return (*this);		// considering (*this) is already orthogonal
		}
	}
	// find orthogonal complementary in general
	Matrix<T> add_orthogonal_complementary() const {
		// before call this, make sure (*this) is orthogonal, you can call Schimidt_orthogonalization() to do that

		// cascade above tow function for simplicity
		return add_w1_orthogonal_complementary().add_rand_orthogonal_complementary();
	}
	
	Matrix<T> erase_cols(Matrix<int> column_ind, bool is_sorted = false) const {

		int cs = column_ind.size();
		int col_left = c - cs;
		if (col_left == 0) {
			return Matrix<T>();
		}

		if (is_sorted == false) {
			column_ind.sort('<');
		}
		Matrix<T> ans(r, col_left);
		int column_ind_ind = 0;
		for (int j = 0; j < c; ++j) {
			if (column_ind_ind < cs && j == column_ind(column_ind_ind)) {
				column_ind_ind++;
			}
			else {
				// a valid column to copy
				for (int i = 0; i < r; ++i) {
					ans(i, j - column_ind_ind) = (*this)(i, j);
				}
			}
		}
		return ans;
	}
	Matrix<T> erase_rows(Matrix<int> row_ind, bool is_sorted = false) const {

		int rs = row_ind.size();
		int row_left = r - rs;
		if (row_left == 0) {
			return Matrix<T>();
		}

		if (is_sorted == false) {
			row_ind.sort('<');
		}
		Matrix<T> ans(row_left, c);

		int row_ind_ind = 0;
		for (int i = 0; i < r; ++i) {
			// remember judging if row_ind_ind is valid in row_ind, very painful teach, time consuming
			if (row_ind_ind < rs && i == row_ind(row_ind_ind)) {

				row_ind_ind++;
			}
			else {
				// a valid row to copy

				for (int j = 0; j < c; ++j) {
					ans(i - row_ind_ind, j) = (*this)(i, j);
				}
			}
		}
		return ans;
	}

	Matrix<T> insert_cols(const Matrix<T>& columns_inserted, const Matrix<int>& column_ind) const {
		int cs_temp = column_ind.size();
		if (cs_temp != 0) {
			int col_size_all = c + cs_temp;
			Matrix<T> ans(r, col_size_all);
			int col_shift_circularly_tmp = 0;
			for (int j = 0; j < col_size_all; ++j) {
				if (col_shift_circularly_tmp >= cs_temp || j != column_ind(col_shift_circularly_tmp)) {
					for (int i = 0; i < r; ++i) {
						ans(i, j) = (*this)(i, j - col_shift_circularly_tmp);
					}
				}
				else {
					for (int i = 0; i < r; ++i) {
						ans(i, j) = columns_inserted(i, col_shift_circularly_tmp);
					}
					col_shift_circularly_tmp++;
				}
			}
			return ans;
		}
		else {
			return *this;
		}
	}
	Matrix<T> insert_rows(const Matrix<T>& rows_inserted, const Matrix<int>& row_ind) const {
		int rs_tmp = row_ind.size();
		if (rs_tmp != 0) {
			int row_size_all = r + rs_tmp;
			Matrix<T> ans(row_size_all, c);
			int row_shift = 0;
			for (int i = 0; i < row_size_all; ++i) {
				if (row_shift >= rs_tmp || i != row_ind(row_shift)) {
					for (int j = 0; j < c; ++j) {
						ans(i, j) = (*this)(i - row_shift, j);
					}
				}
				else {
					for (int j = 0; j < c; ++j) {
						ans(i, j) = rows_inserted(row_shift, j);
					}
					row_shift++;
				}
			}
			return ans;
		}
		else {
			return *this;		// nothing inserted
		}		
	}
	Matrix<T> place_diag() const {		
		int sss = size();
		Matrix<T> ans(sss, sss, '0');
		for (int i = 0; i < sss; ++i) {
			ans(i, i) = (*this)(i);
		}
		return ans;
	}

	// return the index of zero col or row
	Matrix<int> find_zero_col()const {
		Matrix<int> ans(1, c);		// assign ram in advance
		int ans_ind = 0;
		for (int j = 0; j < c; ++j) {
			int i = 0;
			for (; i < r && (*this)(i, j) == 0; ++i);
			if (i != r);
			else{
				ans(ans_ind) = j;
				ans_ind++;
			}
		}
		ans.resize(1, ans_ind);
		return ans;
	}
	Matrix<int> find_zero_row()const {
		Matrix<int> ans(1, r);		// assign ram in advance
		int ans_ind = 0;
		for (int j = 0; j < r; ++j) {
			int i = 0;
			for (; i < c && (*this)(j, i) == 0; ++i);
			if (i != c);
			else {
				ans(ans_ind) = j;
				ans_ind++;
			}
		}
		ans.resize(1, ans_ind, true);
		return ans;
	}

	/**
	 * .Cholesky factorization,
	 * (*this) must be symmetric positive definite matrix, be careful
	 * 
	 * \param L: (*this) = L*L.Transpose(), L is low triangular
	 */
	void chol(Matrix<T>& L) const {
		L = (* this);
		for (int k = 0; k < r - 1; ++k) {		// row ind
			for (int i = k + 1; i < r; ++i) {	// col ind
				T ratio= L(k, i) / L(k, k);			// using symmetric of L
				for (int j = i; j < r; ++j) {	// col ind
					L(i, j) = L(i, j) - ratio * L(k, j); // computed 1 / 6 * (r ^ 3 - r) times
				}
			}
		}
		for (int i = 0; i < r; ++i) {
			L(i, i) = sqrt(L(i, i));	// adjust diag
			for (int j = i + 1; j < r; ++j) {
				L(j, i) = L(i, j) / L(i, i); // cause of split out diag
				L(i, j) = 0; // make L low triangular
			}
		}
	}

	/**
	 * .lanczos iteration to transform a matrix into tridiagonal matrix,
	 * 	(*this) must be Hermitian matrix
	 *  it will be better to still apply arnoldi algorithm for accuracy, be careful
	 *
	 * \param Tridiagonal:	result tridagonal matrix, satisfying Q.transpose()*(*this)*Q=Tridiagonal
	 * \param Q:			orthogonal matrix output
	 * \param iter:			iteration will decide accuracy of ritz values on eigen value, default: row of matrix
	 */
	void lanczos(Matrix<T>& Tridiagonal, Matrix<T>& Q, int iter = -1) const {
		if (r == c) {
			iter = iter == -1 ? r : iter;		// r==c
			Matrix<T> uk(r, 1);
			Q = Matrix<T>(r, iter + 1);
			Tridiagonal = Matrix<T>(iter + 1, iter + 1, '0');
			Matrix<T> x0(r, 1);
			for (int i = 0; i < r; ++i) {
				x0(i) = T(1);
			}
			Matrix<T> qj(r, 1);
			x0.unitize();
			Q.set_col(0, x0);

			// for case k==0
			int k = 0; 
			uk = (*this) * Q.get_col(k);
			int j = 0;
			qj = Q.get_col(j);
			Tridiagonal(j, k) = qj.dot_product(uk);
			uk -= Tridiagonal(j, k) * qj;				// projection

			Tridiagonal(k + 1, k) = uk.norm_2();
			if ((double) my::abs(Tridiagonal(k + 1, k)) < my::zero_approximation) {		
					// we approach my::zero_approximation be 0, all over the class of matrix
				Tridiagonal = Tridiagonal.get_part(0, 0, 0, 0);
				Q = Q.get_part(0, 0, r - 1, 0);
				return;
			}
			Q.set_col(k + 1, uk / Tridiagonal(k + 1, k));

			// for case k>=0
			k++;
			for (; k < iter; ++k) {
				uk = (*this) * Q.get_col(k);
				Tridiagonal(k - 1, k) = Tridiagonal(k, k - 1);		// symmetric

				uk -= Tridiagonal(k - 1, k) * Q.get_col(k - 1);		// big problem when matrix get large
				j = k;
				qj = Q.get_col(j);
				Tridiagonal(j, k) = qj.dot_product(uk);
				uk -= Tridiagonal(j, k) * qj;				// projection

				Tridiagonal(k + 1, k) = uk.norm_2();
				if ((double) my::abs(Tridiagonal(k + 1, k)) < my::zero_approximation) {		
						// we approach my::zero_approximation be 0, all over the class of matrix
					break;
				}
				Q.set_col(k + 1, uk / Tridiagonal(k + 1, k));
			}

			k = k == iter ? k : k + 1;
			Tridiagonal = Tridiagonal.get_part(0, 0, k - 1, k - 1);
			Q = Q.get_part(0, 0, r - 1, k - 1);
		}
		else {
			cout << "matrix not squared, lanczos iteration failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
		}
	}

	/**
	 * .Arnoldi iteration to transform a matrix into a Hessenberg matrix 
	 * 	(*this) must be square matrix
	 * 
	 * \param Hessenberg:	result Hessenberg matrix, satisfying Q.transpose()*(*this)*Q=Hessenberg
	 * \param Q:			orthogonal matrix output
	 * \param iter:			iteration will decide accuracy of ritz values on eigen value, default: row of matrix
	 */
	void arnoldi(Matrix<T>& Hessenberg, Matrix<T>& Q, int iter = -1) const {		// very unstable
		if (r == c) {
			iter = iter == -1 ? r : iter;		// r==c
			Matrix<T> uk(r, 1);
			Q = Matrix<T>(r, iter + 1);
			Hessenberg = Matrix<T>(iter + 1, iter + 1, '0');
			Matrix<T> x0(r, 1);
			for (int i = 0; i < r; ++i) {
				x0(i) = T(2 * my::rand_u() - 1);
			}
			Matrix<T> qj(r, 1);
			x0.unitize();
			Q.set_col(0, x0);
			int k = 0;
			for (; k < iter; ++k) {
				uk = (*this) * Q.get_col(k);
				for (int j = 0; j <= k; ++j) {
					qj = Q.get_col(j);
					Hessenberg(j, k) = qj.dot_product(uk);
					uk -= Hessenberg(j, k) * qj;				// projection
				}

				Hessenberg(k + 1, k) = uk.norm_2();
				if ((double) my::abs(Hessenberg(k + 1, k)) < my::zero_approximation) {		
						// we approach my::zero_approximation be 0, all over the class of matrix
					break;
				}
				Q.set_col(k + 1, uk / Hessenberg(k + 1, k));
			}

			k = k == iter ? k : k + 1;
			Hessenberg = Hessenberg.get_part(0, 0, k - 1, k - 1);
			Q = Q.get_part(0, 0, r - 1, k - 1);
		}
		else {
			cout << "matrix not squared, lanczos iteration failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
		}
	}

	/**
	 * .transform a hessenberg matrix to tridiagonlize matrix by simply set 0
	 * 
	 */
	void tridiagonlize() {
		for (int i = 0; i < r; ++i) {
			for (int j = i + 2; j < c; ++j) {
				(*this)(i, j) = 0;
			}
		}
	}

	Matrix<T> multiply_result_Hessenberg(const Matrix<T>& B) const{
		if (c == B.r) {
			Matrix<T> ans(r, B.c, '0');
			int row_pos;
			int row_pos2;

			int tmp_j = 0;
			for (int t = 0; t < c; t++) {
				row_pos = 0;
				row_pos2 = 0;

				// i==0
				for (int j = 0; j < B.c; j++) {		// ans be Hessenberg
					ans.e[row_pos + j] += e[row_pos2 + t] * B.e[j + tmp_j];
				}
				row_pos2 += c;
				row_pos += B.c;

				for (int i = 1; i < r; i++) {	// the compiler will parallize it, using release mode

					//theoretically will have parallize degree r * B.c

					for (int j = i - 1; j < B.c; j++) {		// ans be Hessenberg
						ans.e[row_pos + j] += e[row_pos2 + t] * B.e[j + tmp_j];
					}
					row_pos2 += c;
					row_pos += B.c;
				}
				tmp_j += B.c;
			}
			return ans;
		}
		else {
			cout << "size error, multiply failed" << endl;
			cout << "r=" << r << ",\tA.c=" << c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			return Matrix<T>();
		}
	}

	Matrix<T> multiply_result_Tridiagonal(const Matrix<T>& B) const {
		if (c == B.r) {


			Matrix<T> ans(r, B.c, '0');

			if (r == 1) {
				// complete that
				for (int j = 0; j < c; ++j) {
					for (int i = 0; i < B.c; ++i) {
					
						ans(i) += (*this)(0, j) * B(j, i);
					}
				}
				return ans;
			}

			int row_pos;
			int row_pos2;

			int tmp_j = 0;
			for (int t = 0; t < c; t++) {
				row_pos = 0;
				row_pos2 = 0;

				// i==0, ans be Tridiagonal, only 2 to compute
				int j = 0;
				ans.e[row_pos + j] += e[row_pos2 + t] * B.e[j + tmp_j];
				j++;

				ans.e[row_pos + j] += e[row_pos2 + t] * B.e[j + tmp_j];
				row_pos2 += c;
				row_pos += B.c;

				for (int i = 1; i < r - 1; i++) {	// the compiler will parallize it, using release mode

					//theoretically will have parallize degree r * B.c

					// ans be Tridiagonal, only 3 to compute
					j = i - 1;
					ans.e[row_pos + j] += e[row_pos2 + t] * B.e[j + tmp_j];
					j++;
					ans.e[row_pos + j] += e[row_pos2 + t] * B.e[j + tmp_j];
					j++;
					ans.e[row_pos + j] += e[row_pos2 + t] * B.e[j + tmp_j];

					row_pos2 += c;
					row_pos += B.c;
				}


				// i==r-1, ans be Tridiagonal, only 2 to compute
				j = r - 2;
				ans.e[row_pos + j] += e[row_pos2 + t] * B.e[j + tmp_j];
				j++;
				ans.e[row_pos + j] += e[row_pos2 + t] * B.e[j + tmp_j];
				row_pos2 += c;
				row_pos += B.c;

				tmp_j += B.c;
			}
			return ans;
		}
		else {
			cout << "size error, multiply failed" << endl;
			cout << "r=" << r << ",\tA.c=" << c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			return Matrix<T>();
		}
	}

	/**
	 * .you have to call arnoldi or lanczos iteration first to turn the matrix into Hessenberg or Tridiagonal_Hermitian
	 *	Here we consider ritz_vec is needed, as the parameter indicates
	 * 
	 * \param ritz_val: the returned ritz value
	 * \param ritz_vec: the returned ritz vaecot
	 * \param is_Hermitian: if false, matrix will be consider as Hessenberg matrix, else, tridiagonal matrix
	 * \param iteration: max iteration to run
	 * \param variance_ratio: if variance of ritz_val less than this value, stop iteration
	 * \return ritz vector if need_ritz_vec is true
	 */
	void ritz_Hessenberg(Matrix<T>& ritz_val, Matrix<T>& ritz_vec, bool is_Hermitian = false, \
		int iteration = 1500, double ritz_precision = my::zero_approximation) {
		// this function is similar to eig_val, but i believe only this function is reliable

		// (*this) must be Hessenberg matrix, or specially, tridiagonal matrix, and [positive definite ?]

		// perform QR iteration, the output may be not correct, need 'eig_shift_inv'

		// this may converge very slowly, a problem

		// need to find a stop criteria, for the case of perturbation on eigen value
		
		if (size() == 1) {			// kick out the case of 1*1 matrix
			ritz_val.resize(1, 1, false);
			ritz_val(0) = e[0];

			ritz_vec.resize(1, 1, false);
			ritz_vec(1) = T(1);
			return;
		}

		if (r == c) {				// we donot seperate r and c in this if clause
			Matrix<T> A(*this);		// prevent changing (*this)
			Matrix<T> Q, R;
			int n = r;
			ritz_val.resize(1, n, false);
			double sqrt_ritz_precision = sqrt(ritz_precision);

			Matrix<T> Q_now(n, n, 'i');
			ritz_vec.set_identity(n, n);

			// dealear enough size to use in loop
			Matrix<T> sub_A(n - 1, n - 1);
			Matrix<int> valid_col(1, n);
			Matrix<int> invalid_col(1, n);
			Matrix<T> A_valid_part(n, n);
			Matrix<T> upT(n, n);

			int ritz_iteration = 0;
			for (; ritz_iteration < iteration; ritz_iteration++) {
				
				//cout << "ritz_iteration = " << ritz_iteration << endl;
				//cout << "A" << A;

				// use this, this is fast
				T a = A(n - 2, n - 2);
				T b = A(n - 2, n - 1);
				T c = A(n - 1, n - 2);
				T d = A(n - 1, n - 1);
				T in_sqrt = (a - d) * (a - d) + 4 * b * c;
				T shift = 0;
				//cout << "shift = " << shift << endl;

				T sqrt_ans = sqrt(in_sqrt);
				if (!isnan((double)sqrt_ans)) {
					//cout << "in_sqrt = " << in_sqrt << endl;
					//cout << "sqrt(in_sqrt) = " << sqrt(in_sqrt) << endl;
					T can_1 = (a + d - sqrt_ans) / 2;
					T can_2 = (a + d + sqrt_ans) / 2;
					if (my::abs(d - can_1) > my::abs(d - can_2)) {
						shift = can_2;		// eigen values of bottom right matrix closer to most bottom right element
					}
					else {
						shift = can_1;
					}
				}
				//cout << "shift = " << shift << endl;

				A.add_diag(-shift);
				//cout << "A_add_shift" << A;
				if (is_Hermitian) {		// this is only for tridiagonal Hermitian matrix
					A.givens_Tridiagonal(Q, R);
					A = R.multiply_result_Tridiagonal(Q);	// O(r^2) faster
				}
				else {
					A.givens_Hessenberg(Q, R);
					A = R.multiply_result_Hessenberg(Q);	// O(r^3) still

					/*A.QR(Q, R);
					A = R * Q;*/
				}

				Q_now = Q_now * Q;							// record Q

				//cout << "Q" << Q;
				//cout << "Q_now" << Q_now;

				A.add_diag(shift);
				int k = n;				// k be the rows number of A
				while (n > 1 && A(n - 1, n - 2) < ritz_precision) {
					ritz_val(n - 1) = A(n - 1, n - 1);
					n--;
				}
				if (n != k) {							// size of A changed, now simplify it, n and k is not changed under this if clause
					do {
						T ritz_now = A(k - 1, k - 1);	// any ritz value close to this is seem as in one group, adding multiplicity
						int n_prime = k - 2;			// find out number of multiplicity 
						for (; n_prime >= n; --n_prime) {
							if (my::abs(ritz_now - A(n_prime, n_prime)) <= sqrt_ritz_precision) {
								// under this case, we consider the eigenvalue are identical
							}
							else {
								break;		// stop counting identical eigen values
							}
						}
						n_prime++;			// with usual case, no eigenvalue are same, n_prime will be k-1

						if (!is_Hermitian) {
							A.get_part(0, 0, n_prime - 1, n_prime - 1, sub_A);				// the left-upper part of A
							//cout << "sub_A" << sub_A;
							sub_A.add_diag(-ritz_val(k - 1));		// as the induction says
							//cout << "A" << A;
							valid_col.resize(1, 0, false);				// record the eigen value without eigenvector
							valid_col.push_back(n_prime);
							invalid_col.resize(1, 0, false);				// record the eigen value without corresponding eigen vectors
							for (int t = n_prime + 1; t < k; ++t) {
								int w = n_prime;
								for (; w < t; ++w) {
									if (A(w, t) == 0);				// normal entry
									else {
										invalid_col.push_back(t);	// the column that can not be eleminated, producing null in ritz vector
										break;
									}
								}
								if (w == t) {
									valid_col.push_back(t);			// normal columns recorded
								}
							}
							A.get_cols(valid_col, A_valid_part, n_prime - 1);				// remain the valid columns, only n_prime rows

							//cout << "A_valid_part" << A_valid_part;
							//cout << "sub_A: add diag" << sub_A;
							upT.set_identity(k, k);					// upT is up triangle Matrix, ensuring last columns weight 1 at diagonal
							upT.set_cols(valid_col, -sub_A.inv_then_multiply(A_valid_part));

							//cout << "upT" << upT;
							//cout << "upT.inv() * A * upT()" << upT.inv() * A * upT;		// the designed similarity transformation

							//ritz_vec = ritz_vec * (Q_now * upT).combine_diag(Matrix<T>(r - k, r - k, 'i'));
							ritz_vec.set_part(0, 0, ritz_vec.get_part(0, 0, r - 1, k - 1) * (Q_now * upT));	// same as the last line

							//cout << "ritz_vec" << ritz_vec;
							if (invalid_col.size() == 0);
							else {						// if has invalid column, set the corresponding eigen vector to a valid one
								ritz_vec.set_cols(invalid_col, ritz_vec.get_col(n_prime), true);			// this is very rare
							}
						}
						else {
							// Hermitian Matrix is sure to be decomposible and, no transformation to the last columns is needed
							ritz_vec.set_part(0, 0, ritz_vec.get_part(0, 0, r - 1, k - 1) * Q_now);
						}

						//cout << "ritz_val" << ritz_val;

						Q_now.set_identity(n_prime, n_prime);
						k = n_prime;
					} while (n < k);		// looping back is very rare

					A = A.get_part(0, 0, n - 1, n - 1);		// decompose A in to a smaller Matrix
					if (n == 1) {
						ritz_val(0) = A(0, 0);
						ritz_iteration++;
						break;
					}
				}
				//cout << "A_back" << A;
			}

			//cout << "ritz_vec: unitize before" << ritz_vec;
			ritz_vec.unitize_col();
			// comparing to Matlab, it still need a rotation in complex domain, 
			// to force the first or second element in the vector be real

			//cout << "ritz_vec: unitize after" << ritz_vec;

#ifdef see_eig_warning
			if (ritz_iteration == iteration) {		// this is nearly possible since the convregence is fast now
				cout << "warning: ritz values are perturbative, ritz iteration=" << ritz_iteration << endl;
				cout << "A" << A;
			}
#endif // see_eig_warning


#ifdef see_eig_iteration_num
			cout << "ritz_val_iteration=" << ritz_iteration << endl;		// keep an eye on this
#endif


			// arrange the eigen_Values selection sort, dosen't matter
			int max_ind;
			for (int i = 0; i < r; ++i) {	// only need to sort A.r eigen value
				ritz_val.max_abs_ele_start_from(i, max_ind);
				ritz_val.switch_ele(i, max_ind);
				ritz_vec.switch_col(i, max_ind);
			}
		}
		else {
			cout << "size error, get ritz_val failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
		}
	}
	/**
	 * .you have to call arnoldi or lanczos iteration first to turn the matrix into Hessenberg or Tridiagonal_Hermitian
	 *	Here we consider ritz vector is not needed
	 * \param ritz_val: the returned ritz value
	 * \param is_Hermitian: if false, matrix will be consider as Hessenberg matrix, else, tridiagonal matrix
	 * \param iteration: max iteration to run
	 * \param variance_ratio: if variance of ritz_val less than this value, stop iteration
	 * \return ritz vector if need_ritz_vec is true
	 */
	void ritz_Hessenberg(Matrix<T>& ritz_val, bool is_Hermitian = false, \
		int iteration = 1500, double ritz_precision = my::zero_approximation) {
		// this function is similar to eig_val, but i believe only this function is reliable

		// (*this) must be Hessenberg matrix, or specially, tridiagonal matrix, and [positive definite ?]

		// perform QR iteration, the output may be not correct, need 'eig_shift_inv'

		// this may converge very slowly, a problem

		// need to find a stop criteria, for the case of perturbation on eigen value

		if (size() == 1) {			// kick out the case of 1*1 matrix
			ritz_val.resize(1, 1, false);
			ritz_val(0) = e[0];
			return;
		}

		if (r == c) {				// we donot seperate r and c in this if clause
			Matrix<T> A(*this);		// prevent changing (*this)
			Matrix<T> Q, R;
			int n = r;
			ritz_val.resize(1, n, false);

			int ritz_iteration = 0;
			for (; ritz_iteration < iteration; ritz_iteration++) {

				//cout << "ritz_iteration = " << ritz_iteration << endl;
				//cout << "A" << A;

				// use this, this is fast
				T a = A(n - 2, n - 2);
				T b = A(n - 2, n - 1);
				T c = A(n - 1, n - 2);
				T d = A(n - 1, n - 1);
				T in_sqrt = (a - d) * (a - d) + 4 * b * c;
				T shift = 0;
				//cout << "shift = " << shift << endl;

				T sqrt_ans = sqrt(in_sqrt);
				if (!isnan((double)sqrt_ans)) {
					//cout << "in_sqrt = " << in_sqrt << endl;
					//cout << "sqrt(in_sqrt) = " << sqrt(in_sqrt) << endl;
					T can_1 = (a + d - sqrt_ans) / 2;
					T can_2 = (a + d + sqrt_ans) / 2;
					if (my::abs(d - can_1) > my::abs(d - can_2)) {
						shift = can_2;		// eigen values of bottom right matrix closer to most bottom right element
					}
					else {
						shift = can_1;
					}
				}
				//cout << "shift = " << shift << endl;

				A.add_diag(-shift);
				//cout << "A_add_shift" << A;
				if (is_Hermitian) {		// this is only for tridiagonal Hermitian matrix
					A.givens_Tridiagonal(Q, R);
					A = R.multiply_result_Tridiagonal(Q);	// O(r^2) faster
				}
				else {
					A.givens_Hessenberg(Q, R);
					A = R.multiply_result_Hessenberg(Q);	// O(r^3) still

					/*A.QR(Q, R);
					A = R * Q;*/
				}

				//cout << "Q" << Q;

				A.add_diag(shift);
				int k = n;				// k be the rows number of A
				while (n > 1 && A(n - 1, n - 2) < ritz_precision) {
					ritz_val(n - 1) = A(n - 1, n - 1);
					n--;
				}
				if (n != k) {			// size of A changed, now simplify it, n and k is not changed under this if clause					
					A = A.get_part(0, 0, n - 1, n - 1);	// decompose A in to a smaller Matrix
					if (n == 1) {
						ritz_val(0) = A(0, 0);
						ritz_iteration++;
						break;
					}
				}
				//cout << "A_back" << A;
			}

#ifdef see_eig_warning
			if (ritz_iteration == iteration) {		// this is nearly possible since the convregence is fast now
				cout << "warning: ritz values are perturbative, ritz iteration=" << ritz_iteration << endl;
				cout << "A" << A;
			}
#endif // see_eig_warning


#ifdef see_eig_iteration_num
			cout << "ritz_val_iteration=" << ritz_iteration << endl;		// keep an eye on this
#endif

			// arrange the eigen_Values selection sort, dosen't matter
			int max_ind;
			for (int i = 0; i < r; ++i) {	// only need to sort A.r eigen value
				ritz_val.max_abs_ele_start_from(i, max_ind);
				ritz_val.switch_ele(i, max_ind);				
			}
		}
		else {
			cout << "size error, get ritz_val failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
		}
	}

	// before calling the function, make sure (*this) is non-defficient (by QR,RQ opeartion). or unexpected behavior will occur
	// adding arnodi iteration before calling ritz_inner, we do not consider finding partly ritz value, since it is unstable
	// but if you do want to, try by passing parameter iter over arnoldi, but the resulted ritz value is very unstable
	// in this case, we assume ritz vectors are needed, as the parameter indicates
	void ritz_non_defficient(Matrix<T>& ritz_val, Matrix<T>& ritz_vec, bool is_Hermitian = false, \
		int iteration = 1500, double ritz_precision = my::zero_approximation) const {

		Matrix<T> Hessenberg, Q;
		arnoldi(Hessenberg, Q);
		//lanczos(Hessenberg, Q);		// if is_Hermitian, use lanczos will be faster but unstable, we choose arnoldi whatever

		//cout << "Hessenberg" << Hessenberg;
		Hessenberg.ritz_Hessenberg(ritz_val, ritz_vec, is_Hermitian, iteration, ritz_precision);
		
		ritz_vec = Q * ritz_vec;
	}
	// in this case we consider ritz vector not needed, as parameter indicates
	void ritz_non_defficient(Matrix<T>& ritz_val, bool is_Hermitian = false, \
		int iteration = 1500, double ritz_precision = my::zero_approximation) const {

		Matrix<T> Hessenberg, Q;
		arnoldi(Hessenberg, Q);
		//lanczos(Hessenberg, Q);		// if is_Hermitian, use lanczos will be faster but unstable, we choose arnoldi whatever

		//cout << "Hessenberg" << Hessenberg;
		Hessenberg.ritz_Hessenberg(ritz_val, is_Hermitian, iteration, ritz_precision);
	}

	/**
	 * .find ritz values and vectors, for universe case of Matrix
	 * 
	 * \param ritz_val: the returned ritz value
	 * \param ritz_vec: the returned ritz vector
	 * \param is_Hermitian: whether the Matrix is Hermitian
	 * \param num_of_ritz_need: -2 indicates only non-zero ritz value needed, -1 indicates all ritz value need,
	 *  positive number indicates exact number of ritz value needed
	 * \param ritz_iter: max iteration to run
	 * \param ritz_precision: if variance of ritz_val less than this value, stop iteration
	 * \return ritz vector if need_ritz_vec is true
	 */
	 void ritz(Matrix<T>& ritz_val, Matrix<T>& ritz_vec, bool is_Hermitian = false, \
		int num_of_ritz_need = -2, int ritz_iter = 1500, double ritz_precision = my::zero_approximation) const{

		if (r != c) {
			cout << "size error, [ritz] failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return;
		}
		if (num_of_ritz_need == 0 || num_of_ritz_need < -2) {
			cout << "do not need ritz value and vector, return empty" << endl;
			cout << "num_of_ritz_need = " << num_of_ritz_need;
			return;
		}

		// condider case for not Hermitian first
		Matrix<T> Q, R;
		QR(Q, R);
		if (R.size() == 0) {
			return;	// no need to consider further for zero rank matrix
		}

		Matrix<T> A = R * Q;		// use first iteration to throw zero space of (*this), to be optimized
		A.ritz_non_defficient(ritz_val, ritz_vec, is_Hermitian, ritz_iter, ritz_precision);
		int rk = ritz_val.col();

		ritz_vec = Q * ritz_vec;		// this will be the ritz vector of (*this)
		
		if (num_of_ritz_need == -2 || (rk == c && num_of_ritz_need == -1)) {
			// do not need to consider the null space
		}
		else if (num_of_ritz_need == -1 || num_of_ritz_need > rk) {
			if (num_of_ritz_need > c) {
				cout << "Warning: the number of ritz value you required is out of Matrix size" << endl;
				cout << "num_of_ritz_need = " << num_of_ritz_need << "\t c = " << c << endl;
				cout << "returning the max number of eigenvalue and eigenvector as possible" << endl;
				num_of_ritz_need = c;
			}
			num_of_ritz_need = num_of_ritz_need == -1 ? c : num_of_ritz_need;
			ritz_val = ritz_val.combine_right(Matrix<T>(1, num_of_ritz_need - rk, '0'));
			Matrix<T> row_ortho = Schimidt_orthogonalization(false);
			//cout << "row_ortho" << row_ortho;
			Matrix<int> zero_row_ind = row_ortho.find_zero_row();
			//cout << "zero_row_ind" << zero_row_ind;
			row_ortho = row_ortho.erase_rows(zero_row_ind);
			//cout << "cut_row" << cut_row;
			row_ortho = row_ortho.add_orthogonal_complementary();
			ritz_vec = ritz_vec.combine_right(row_ortho.get_part(rk, 0, num_of_ritz_need - 1, c - 1).Hermit());// change from Transpose
		}
		else{		// (num_of_ritz_need <= rk)
			ritz_val.resize(1, num_of_ritz_need, true);		// setting true or false doesn't matter, since it's shrinking
			ritz_vec = ritz_vec.get_part(0, 0, -1, num_of_ritz_need - 1);
		}
		//cout << "ritz_Values: after" << ritz_Values;								// check if it can contains same ritz values
	}
	/**
	 * .find ritz values only as the parameters indicate
	 *
	 */
	void ritz(Matrix<T>& ritz_val, bool is_Hermitian = false, \
		int num_of_ritz_need = -2, int ritz_iter = 1500, double ritz_precision = my::zero_approximation) const {

		if (r != c) {
			cout << "size error, [ritz] failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			return;
		}
		if (num_of_ritz_need == 0 || num_of_ritz_need < -2) {
			cout << "do not need ritz value and vector, return empty" << endl;
			cout << "num_of_ritz_need = " << num_of_ritz_need;
			return;
		}

		// condider case for not Hermitian first
		Matrix<T> Q, R;
		QR(Q, R);
		if (R.size() == 0) {
			return;	// no need to consider further for zero rank matrix
		}

		Matrix<T> A = R * Q;		// use first iteration to throw zero space of (*this), to be optimized
		A.ritz_non_defficient(ritz_val, is_Hermitian, ritz_iter, ritz_precision);
		int rk = ritz_val.col();

		if (num_of_ritz_need == -2 || (rk == c && num_of_ritz_need == -1)) {
			// do not need to consider the null space
		}
		else if (num_of_ritz_need == -1 || num_of_ritz_need > rk) {
			if (num_of_ritz_need > c) {
				cout << "Warning: the number of ritz value you required is out of Matrix size" << endl;
				cout << "num_of_ritz_need = " << num_of_ritz_need << "\t c = " << c << endl;
				cout << "returning the max number of eigenvalue and eigenvector as possible" << endl;
				num_of_ritz_need = c;
			}
			num_of_ritz_need = num_of_ritz_need == -1 ? c : num_of_ritz_need;
			ritz_val = ritz_val.combine_right(Matrix<T>(1, num_of_ritz_need - rk, '0'));
		}
		else {		// (num_of_ritz_need <= rk)
			ritz_val.resize(1, num_of_ritz_need);			
		}
		//cout << "ritz_Values: after" << ritz_Values;								// check if it can contains same ritz values

	}

	// givens rotation to get QR from squared Hessenberg matrix
	void givens_Hessenberg(Matrix<T>& Q, Matrix<T>& R) const {		// a full Q will be retruned, to reduce size consider QR function
		if (r == c) {
			R = (*this);
			Q.set_identity(r, c);
			for (int i = 0; i < c - 1; ++i) {
				int j = i + 1;		// the row to eleminate

				// this if clause is with small probability, donot consider in complexity, always skip
				if ((double) my::abs(R(i, i)) < my::zero_approximation) {
					// permute a non-zero row to ith row
					
					// since it is Hessenberg matrix, permute the next row is ok
					if ((double) my::abs(R(j, i)) < my::zero_approximation) {
						// in this case we cannot find a valid perputation
						continue;		// skip this column directly
					}
					else {		// permute this row with next row, keep diagonal position non-zero
						for (int k = i; k < c; ++k) {
							T a = R(i, k);
							R(i, k) = R(j, k);
							R(j, k) = a;

							T q = Q(i, k);
							Q(i, k) = Q(j, k);
							Q(j, k) = q;
						}
					}
				}

				T a = R(i, i);
				T e = R(j, i);

				T sq = sqrt(a * my::conj(a) + e * my::conj(e));
				T co = a / sq;
				T si = e / sq;
				T co_conj = my::conj(co);
				T si_conj = my::conj(si);

				// R=G*R
				// Q=Q*G^H

				// multiplication of a givens matrix can be optimized as follow
				R(i, i) = co_conj * a + si_conj * e;
				R(j, i) = 0;
				for (int k = i + 1; k < c; ++k) {		// observing that R is uptriangle
					a = R(i, k);
					e = R(j, k);
					R(i, k) = co_conj * a + si_conj * e;
					R(j, k) = co * e - si * a;
				}

				for (int k = 0; k <= j; ++k) {		// observing that Q is Hessenberg
					a = Q(k, i);					// complexity: O(r*r)
					e = Q(k, j);
					Q(k, i) = a * co + e * si;
					Q(k, j) = e * co_conj - a * si_conj;
				}
			}
			
			// for the last column
			int i = c - 1;
			if ((double)my::abs(R(i, i)) >= my::zero_approximation) {
				T a = R(i, i);
				T sq = my::abs(a);
				T co = a / sq;
			 
				// multiplication of a givens matrix can be optimized as follow
				R(i, i) = sq;
				for (int k = 0; k < c; ++k) {		// observing that Q is Hessenberg, complexity: O(r*r)
					Q(k, i) *= co;
				}
			}
		}
		else {
			cout << "size error, get eigen_val failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
		}
	}

	// givens rotation to get QR from squared Tridiagonal matrix
	void givens_Tridiagonal(Matrix<T>& Q, Matrix<T>& R) const {
		// consider the diagonal element of R is non-zero
		if (r == c) {
			R = (*this);
			Q = Matrix<T>(r, c, 'i');
			for (int i = 0; i < c - 1; ++i) {
				int j = i + 1;		// the row to eleminate

				T a = R(i, i);
				T e = R(j, i);

				// this part is not fixed for Complex matrix !! to be done
				T sq = sqrt(a * my::conj(a) + e * my::conj(e));
				T co = a / sq;
				T si = e / sq;
				T co_conj = my::conj(co);
				T si_conj = my::conj(si);

				// R=G*R
				// Q=G*Q

				// multiplication of a givens matrix can be optimized as follow
				R(i, i) = co_conj * a + si_conj * e;
				R(j, i) = 0;
				for (int k = i + 1; k < r && k < i + 3; ++k) {		// observing that R is Tridiagonal
					a = R(i, k);
					e = R(j, k);

					R(i, k) = co_conj * a + si_conj * e;
					R(j, k) = co * e - si * a;

					/*R(i, k) = co * a + si * e;
					R(j, k) = -si * a + co * e;*/
				}

				for (int k = 0; k <= j; ++k) {		// observing that Q is Hessenberg
					a = Q(k, i);		// complexity: O(r*r)
					e = Q(k, j);

					Q(k, i) = a * co + e * si;
					Q(k, j) = e * co_conj - a * si_conj;

					/*Q(k, i) = co * a + si * e;
					Q(k, j) = -si * a + co * e;*/
				}
			}

			// for the last column
			int i = c - 1;
			if ((double)my::abs(R(i, i)) >= my::zero_approximation) {
				T a = R(i, i);
				T sq = my::abs(a);
				T co = a / sq;

				// multiplication of a givens matrix can be optimized as follow
				R(i, i) = sq;
				for (int k = 0; k < c; ++k) {		// observing that Q is Hessenberg, complexity: O(r*r)
					Q(k, i) *= co;
				}
			}
		}
		else {
			cout << "size error, get eigen_val failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
		}
	}

	// QR decomposition
	void QR(Matrix<T>& Q, Matrix<T>& R, bool need_all_Q = false, \
		double QR_precision = my::zero_approximation * my::zero_approximation) const {
		// do we need the full form of QR factorization? consider throwing away 0 column vectors in Q and corresponding row vector in R
		if (r >= c) {
			// QR fractorization is only useful for matrix with row >= column, note that matrix should have full column rank

			// how to deal with deficient matrix
			Q = *this;
			R = Matrix<T>(c, c, '0');
			Matrix<int> zero_col_ind(1, c);		// allocate enough space in advance
			int zero_col_ind_ind = 0;
			int row_pos;
			for (int i = 0; i < c; i++) {		// col ind
				row_pos = i * R.c;

				T col_i_norm2 = 0;
				for (int j = 0; j < r; ++j) {		// row ind
					col_i_norm2 += my::conj(Q(j, i)) * Q(j, i);
				}
				if ((double) col_i_norm2 >= QR_precision) {			// square of my::zero_approximation
					R.e[row_pos + i] = sqrt(col_i_norm2);	//get R[i][i]

					for (int j = 0; j < r; ++j) {			// row ind
						Q(j, i) /= R.e[row_pos + i];
					}

					for (int k = i + 1; k < c; k++) {		// the varible l(a letter) and 1(a number) is uneasy to distinguish, use k
						T qi_dot_qk = 0;
						for (int j = 0; j < r; ++j) {
							qi_dot_qk += my::conj(Q(j, i)) * Q(j, k);		// complexity O(c * c * r)
						}
						R.e[row_pos + k] = qi_dot_qk;					//R[i][k]=q[i]^T*q[k]

						for (int j = 0; j < r; ++j) {
							Q(j, k) -= R.e[row_pos + k] * Q(j, i);	//q[k]=q[k]-R[i][k] * q[i]
						}
					}
				}
				else {
					zero_col_ind(zero_col_ind_ind) = i;
					zero_col_ind_ind++;
					R.e[row_pos + i] = 0;
					for (int j = 0; j < r; ++j) {		// row ind
						Q(j, i) = 0;
					}
					for (int k = i + 1; k < c; k++) {
						R.e[row_pos + k] = 0;					//R[i][k]=q[i]^T*q[k]
					}
				}
			}
			zero_col_ind.resize(1, zero_col_ind_ind, true);		// seting ture will make no difference
			if (zero_col_ind.size() == 0);
			else {
				Q = Q.erase_cols(zero_col_ind);
				if (!need_all_Q) {
					R = R.erase_rows(zero_col_ind);
				}
				else {
					Matrix<T> Q_add_orthogonal_complementary = Q.add_orthogonal_complementary();
					Q = Q.insert_cols(Q_add_orthogonal_complementary.get_part(0, Q.col(), Q.row() - 1, c - 1), zero_col_ind);
				}
			}
		}
		else {
			// r<c in this case Q has size r*r and R has size r*n, which is not useful, except for getting orthogonal column base
			Matrix<T> Left = (*this).get_part(0, 0, r - 1, r - 1);
			Matrix<T> Right = (*this).get_part(0, r, r - 1, c - 1);
			Matrix<T> R_Left;
			Left.QR(Q, R_Left, need_all_Q);
			R = R_Left.combine_right(Q.Hermit() * Right);		// change from Transpose
		}
	}
	
	/**
	 * . computing svd decomposition, U, Sigma, VH are all needed, (*this) = U * Sigma * VH
	 * 
	 * \param U
	 * \param Sigma
	 * \param VH 
	 * \param num_of_singular_need: -2 indicates only non-zero singular value needed, -1 indicates all singular value need,
	 *  positive number indicates exact number of singular value needed
	 */
	void svd(Matrix<T>& U, Matrix<T>& Sigma, Matrix<T>& VH, \
		int num_of_singular_need = -2) const {

		if (num_of_singular_need == 0 || num_of_singular_need < -2) {
			cout << "do not need singular value and vector, return empty" << endl;
			cout << "num_of_singular_need = " << num_of_singular_need;
			return;
		}

		if (r < c) {
			// a fat matrix, compute its transpose
			Matrix<T> tmp;
			Hermit().svd(tmp, Sigma, VH, num_of_singular_need);
			U = VH.Hermit();		// change from Transpose
			VH = tmp.Hermit();
			return;
		}

		
		// correcting this computation, for element of Sigma>=0
		//int rk = rank();
		Matrix<T> Q, R;		// do QR to reduce size
		QR(Q, R);

		Matrix<T> RH = R.Hermit();
		// squared condition number, consider fixing it, but using svd make the QR iteration less, for large entries
		Matrix<T> R_RH = R * RH;
		//cout << "R_RH" << R_RH;
		Matrix<T> sigma_2;
		//cout << "A_AH.rank()" << A_AH.rank() << endl;

		// update to ritz, more faster, note that R_RH is hermitian non-deficient matrix for sure
		R_RH.ritz_non_defficient(sigma_2, U, true);
		int rk = sigma_2.col();

		Sigma.set_identity(rk, rk);				// squared matrix
		for (int i = 0; i < rk; ++i) {
			Sigma(i, i) = 1 / sqrt(sigma_2(i));	// sigma_inv for now
		}

		VH = (RH * U * Sigma).Hermit();

		for (int i = 0; i < rk; ++i) {
			Sigma(i, i) = 1 / Sigma(i, i);		// back to sigma
		}
		U = Q * U;

		if (num_of_singular_need == -2 || (rk == c && num_of_singular_need == -1));
		else if (num_of_singular_need == -1 || num_of_singular_need > rk) {
			if (num_of_singular_need > c) {
				cout << "Warning: the number of singular value you required is out of Matrix size" << endl;
				cout << "num_of_ritz_need = " << num_of_singular_need << "\t c = " << c << endl;
				cout << "returning the max number of eigenvalue and eigenvector as possible" << endl;
				num_of_singular_need = c;
			}
			num_of_singular_need = num_of_singular_need == -1 ? c : num_of_singular_need;
			Sigma = Sigma.combine_diag(Matrix<T>(num_of_singular_need - rk, num_of_singular_need - rk, '0'));

			U = U.add_orthogonal_complementary().get_part(0, 0, -1, num_of_singular_need - 1);		// can be optimized
			VH = VH.add_orthogonal_complementary().get_part(0, 0, num_of_singular_need - 1, -1);
		}
		else {		// (num_of_ritz_need <= rk)
			Sigma = Sigma.get_part(0, 0, num_of_singular_need - 1, num_of_singular_need - 1);
			U = U.get_part(0, 0, -1, num_of_singular_need - 1);
			VH = VH.get_part(0, 0, num_of_singular_need - 1, -1);
		}
	}
	// find svd values only
	void svd(Matrix<T>& Sigma, int num_of_singular_need = -2) const {

		if (num_of_singular_need == 0 || num_of_singular_need < -2) {
			cout << "do not need singular value and vector, return empty" << endl;
			cout << "num_of_singular_need = " << num_of_singular_need;
			return;
		}

		if (r < c) {
			// a fat matrix, compute its transpose
			Hermit().svd(Sigma, num_of_singular_need);
			return;
		}


		// correcting this computation, for element of Sigma>=0
		//int rk = rank();
		Matrix<T> Q, R;		// do QR to reduce size
		QR(Q, R);

		Matrix<T> RH = R.Hermit();
		// squared condition number, consider fixing it, but using svd make the QR iteration less, for large entries
		Matrix<T> R_RH = R * RH;
		//cout << "R_RH" << R_RH;

		// update to ritz, more faster, note that R_RH is hermitian non-deficient matrix for sure
		R_RH.ritz_non_defficient(Sigma, true);
		int rk = Sigma.col();

		for (int i = 0; i < rk; ++i) {
			Sigma(i) = sqrt(Sigma(i));
		}

		if (num_of_singular_need == -2 || (rk == c && num_of_singular_need == -1));
		else if (num_of_singular_need == -1 || num_of_singular_need > rk) {
			if (num_of_singular_need > c) {
				cout << "Warning: the number of singular value you required is out of Matrix size" << endl;
				cout << "num_of_ritz_need = " << num_of_singular_need << "\t c = " << c << endl;
				cout << "returning the max number of eigenvalue and eigenvector as possible" << endl;
				num_of_singular_need = c;
			}
			num_of_singular_need = num_of_singular_need == -1 ? c : num_of_singular_need;
			Sigma = Sigma.combine_right(Matrix<T>(1, num_of_singular_need - rk, '0'));
		}
		else {		// (num_of_ritz_need <= rk)
			Sigma = Sigma.get_part(0, 0, 0, num_of_singular_need - 1);
		}
	}

	//extend functions (return a new matrix)
	Matrix<T> Transpose() const {
		if (!isVector()) {
			Matrix<T> ans(c, r);
			int row_pos;
			int col_pos;
			for (int i = 0; i < r; i++) {
				row_pos = i * c;
				col_pos = i;
				for (int j = 0; j < c; j++) {
					ans.e[col_pos] = e[row_pos + j];
					col_pos += ans.c;
				}
			}
			return ans;
		}
		else {
			Matrix<T> ans = (*this);
			ans.resize(c, r, true);		// if is a vector, resize it is all okay
			return ans;
		}
	}
	Matrix<T> Hermit() const {
		Matrix<T> ans(c, r);
		int row_pos;
		int col_pos;
		for (int i = 0; i < r; i++) {
			row_pos = i * c;
			col_pos = i;
			for (int j = 0; j < c; j++) {
				ans.e[col_pos] =  my::conj(e[row_pos + j]);
				col_pos += ans.c;
			}
		}
		return ans;
	}
	/*
	   *this B
	   C     D
	 */
	Matrix<T> combine(const Matrix<T>& B, const Matrix<T>& C, const Matrix<T>& D) const {
		if (r == B.r && c == C.c && B.c == D.c && C.r == D.r) {
			int row_all = r + C.r;
			int column_all = c + B.c;
			int row_pos;
			Matrix<T> ans(row_all, column_all);
			for (int i = 0; i < row_all; i++) {
				row_pos = i * column_all;
				for (int j = 0; j < column_all; j++) {
					if (i < r && j < c)
						ans.e[row_pos + j] = e[i * c + j];
					else if (i < r && j >= c)
						ans.e[row_pos + j] = B.e[i * B.c + j - c];
					else if (i >= r && j < c)
						ans.e[row_pos + j] = C.e[(i - r) * C.c + j];
					else if (i >= r && j >= c)
						ans.e[row_pos + j] = D.e[(i - r) * D.c + j - c];
				}
			}
			return ans;
		}
		else {
			cout << "size error, combine failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			cout << "C.r=" << C.r << ",\tC.c=" << C.c << endl;
			cout << "D.r=" << D.r << ",\tD.c=" << D.c << endl;
			return Matrix<T>();
		}
	}
	/*
	   *this 0
	   0     B
	 */
	Matrix<T> combine_diag(const Matrix<T>& B) const {
		if (B.size() != 0) {
			int row_all = r + B.r;
			int column_all = c + B.c;
			int row_pos;
			Matrix<T> ans(row_all, column_all);
			for (int i = 0; i < row_all; i++) {
				row_pos = i * column_all;
				for (int j = 0; j < column_all; j++) {
					if (i < r && j < c)
						ans.e[row_pos + j] = e[i * c + j];
					else if (i >= r && j >= c)
						ans.e[row_pos + j] = B.e[(i - r) * B.c + j - c];
					else
						ans.e[row_pos + j] = T(0);
				}
			}
			return ans;
		}
		else {
			return (*this);
		}
	}
	/*
	   *this B
	 */
	Matrix<T> combine_right(const Matrix<T>& B) const {
		if (B.c == 0)
			return (*this);

		if (r == B.r) {
			int row_all = r;
			int column_all = c + B.c;
			int row_pos;
			Matrix<T> ans(row_all, column_all);
			for (int i = 0; i < row_all; i++) {
				row_pos = i * column_all;
				for (int j = 0; j < column_all; j++) {
					if (j < c)
						ans.e[row_pos + j] = e[i * c + j];
					else
						ans.e[row_pos + j] = B.e[i * B.c + j - c];
				}
			}
			return ans;
		}
		else {
			cout << "size error, combine_right failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			return Matrix<T>();
		}
	}
	/*
	   *this
	   B
	 */
	Matrix<T> combine_down(const Matrix<T>& B) const {
		if (B.r == 0) {
			return (*this);
		}

		if (c == B.c) {
			int row_all = r + B.r;
			int column_all = c;
			int row_pos;
			Matrix<T> ans(row_all, column_all);
			for (int i = 0; i < row_all; i++) {
				row_pos = i * column_all;
				for (int j = 0; j < column_all; j++) {
					if (i < r)
						ans.e[row_pos + j] = e[i * c + j];
					else
						ans.e[row_pos + j] = B.e[(i - r) * B.c + j];
				}
			}
			return ans;
			
		}
		else {
			cout << "size error, combine_down failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			return Matrix<T>();
		}
	}

	/**
	 * .Hamming distance is the number of different elements between (*this) and B
	 */
	int Hamming_distance(const Matrix<T> B) const {
		int this_size = r * c;
		int B_size = B.size();
		if (this_size == B_size) {
			int result = 0;
			for (int i = 0; i < this_size; ++i) {
				result += (*this)(i) != B(i);
			}
			return result;
		}
		else {
			cout << "size error, Hamming_distance fail" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			return 0;
		}
	}
	/**
	 * .Hamming weight is the number of non-zero elements
	 */
	int Hamming_weight() const {
		int this_size = r * c;
		int result = 0;
		for (int i = 0; i < this_size; ++i) {
			result += (*this)(i) != 0;
		}
		return result;
	}

	Matrix<T> operator - () const {
		Matrix<T> ans(r, c);
		int sss = size();
		for (int i = 0; i < sss; i++)
			ans.e[i] = -e[i];
		return ans;
	}
	void operator += (const Matrix<T>& B) {
		if (r == B.r && c == B.c) {
			int sss = size();
			for (int i = 0; i < sss; i++)
				e[i] += B.e[i];
		}
		else {
			cout << "size error, += failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
		}
	}
	void operator -= (const Matrix<T>& B) {
		if (r == B.r && c == B.c) {
			int sss = size();
			for (int i = 0; i < sss; i++)
				e[i] -= B.e[i];
		}
		else {
			cout << "size error, += failed" << endl;
			cout << "r=" << r << ",\tc=" << c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
		}
	}
	void operator *= (const Matrix<T>& B) {
		*this = (*this) * B;		// this cannot be written as in-place operation
	}
	void operator *= (const T& num) {
		int sss = size();
		for (int i = 0; i < sss; i++)
			e[i] *= num;
	}
	void operator /= (const Matrix<T>& B) {
		*this = (*this) * B.inv();	// this cannot be written as in-place operation
	}
	void operator /= (const T& num) {
		int sss = size();
		for (int i = 0; i < sss; i++)
			e[i] /= num;
	}

	// all friends of Matrix, saving IO
	friend Matrix<T> operator + (const Matrix<T>& A, const Matrix<T>& B) {
		if (A.r == B.r && A.c == B.c) {
			Matrix<T> ans(A.r, A.c);
			int sss = A.r * A.c;
			for (int i = 0; i < sss; i++)
				ans.e[i] = A.e[i] + B.e[i];
			return ans;
		}
		else {
			cout << "size error, add failed" << endl;
			cout << "A.r=" << A.r << ",\tA.c=" << A.c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			return Matrix<T>();
		}
	}
	friend Matrix<T> operator - (const Matrix<T>& A, const Matrix<T>& B) {
		if (A.r == B.r && A.c == B.c) {
			Matrix<T> ans(A.r, A.c);
			int sss = A.r * A.c;
			for (int i = 0; i < sss; i++)
				ans.e[i] = A.e[i] - B.e[i];
			return ans;
		}
		else {
			cout << "size error, minus failed" << endl;
			cout << "A.r=" << A.r << ",\tA.c=" << A.c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			return Matrix<T>();
		}
	}
	// the following is optimized for parallization
	friend Matrix<T> operator * (const Matrix<T>& A, const Matrix<T>& B) {
		if (A.c == B.r) {
			Matrix<T> ans(A.r, B.c, '0');
			int row_pos;
			int row_pos2;

			int tmp_j = 0;
			for (int t = 0; t < A.c; t++) {
				row_pos = 0; 
				row_pos2 = 0;
				for (int i = 0; i < A.r; i++) {	// the compiler will parallize it, using release mode
					
					//theoretically will have parallize degree A.r * B.c

					for (int j = 0; j < B.c; j++) {
						ans.e[row_pos + j] += A.e[row_pos2 + t] * B.e[j + tmp_j];
					}
					row_pos2 += A.c;
					row_pos += B.c;
				}
				tmp_j += B.c;
			}
			return ans;
		}
		else {
			cout << "size error, multiply failed" << endl;
			cout << "A.r=" << A.r << ",\tA.c=" << A.c << endl;
			cout << "B.r=" << B.r << ",\tB.c=" << B.c << endl;
			return Matrix<T>();
		}
	}
	friend Matrix<T> operator * (const T& num, const Matrix<T>& A) {
		Matrix<T> ans(A.r, A.c);
		int sss = A.size();
		for (int i = 0; i < sss; i++)
			ans.e[i] = A.e[i] * num;
		return ans;
	}
	friend Matrix<T> operator * (const Matrix<T>& A, const T& num) {
		return num * A;
	}
	friend Matrix<T> operator / (const Matrix<T>& A, const Matrix<T>& B) {
		return A * B.inv();
	}
	friend Matrix<T> operator / (const T& num, const Matrix<T>& A) {
		return num * A.inv();
	}
	friend Matrix<T> operator / (const Matrix<T>& A, const T& num) {
		return (1 / num) * A;
	}
	friend bool operator == (const Matrix<T>& A, const Matrix<T>& B) {
		if (A.r == B.r && A.c == B.c) {
			int size_all = A.r * A.c;
			for (int i = 0; i < size_all; i++)
				if (A.e[i] == B.e[i]);
				else
					return false;
			return true;
		}
		else
			return false;
	}
	friend bool operator != (const Matrix<T>& A, const Matrix<T>& B) {
		return !(A == B);
	}
	friend bool operator >(const Matrix<T>& A, const Matrix<T>& B) {
		return false;
	}
	friend bool operator <(const Matrix<T>& A, const Matrix<T>& B) {
		return false;
	}
	friend bool operator >=(const Matrix<T>& A, const Matrix<T>& B) {
		return false;
	}
	friend bool operator <=(const Matrix<T>& A, const Matrix<T>& B) {
		return false;
	}

	// output a matrix
	friend ostream& operator << (ostream& out, const Matrix<T>& A) {
		out << " (size: " << A.r << " * " << A.c << ")" << endl;
		out << setprecision(4);
		int row_pos;
		bool is_ordered_structure = !(T() < T());

		for (int i = 0; i < A.r; i++) {
			row_pos = i * A.c;
			for (int j = 0; j < A.c; j++) {
				if (is_ordered_structure /*&& A.e[row_pos + j] >= my::zero_approximation*/)
					out << setw(14) << A.e[row_pos + j];
				/*else if(A.is_ordered_structure)
					out<< setw(14) << 0.000;*/
					// approach my::zero_approximation be 0 at output, may cause confusion
				else {
					//out << setw(6) << A.e[row_pos + j];	// for GF display, no decimal point

					// GF2 display
					out << setw(3) << A.e[row_pos + j];

					// for display over lenovo notebook
					/*if ((j + 1) % 26 == 0) {
						out << " ";
					}*/
				}
			}
			out << endl;
		}
		out << setprecision(6);
		return out;
	}
	// input 2 enter to finish the matrix
	friend istream& operator >> (istream& in, Matrix<T>& A) {
		vector<T> v;
		char c;
		T num = T();

		//input process
		int input_column = 0;
		int input_row = 0;
		while ((c = getchar()) != '\n') {
			ungetc(c, stdin);
			input_column = 0;
			while ((c = getchar()) != '\n')
			{
				ungetc(c, stdin);
				in >> num;
				v.push_back(num);
				input_column++;
			}
			input_row++;
		}
		A = Matrix<T>(input_row, input_column, v);
		return in;
	}

	void print() const {
		int row_strat_ind;

		for (int i = 0; i < r; i++) {
			row_strat_ind = i * c;
			for (int j = 0; j < c; j++) {
				if (T() < T()) {
					cout << e[row_strat_ind + j] << ",";		// for GF2 output
				}
				else {
					cout << e[row_strat_ind + j] << ",";
				}
			}
			cout << endl;
		}
		cout << endl;
	}

	/* --------- re-write the GE and GJE: left identity ------------- */

	/**
	 * . this step is after row transformation to up triangle
	 *	 permute the column of (*this) to get an left up triangle matrix,
	 *   not disturbing the original column order as much as possible
	 *
	 * \return
	 */
	void col_permute_up_triangle_to_full_rank_left_v2(Matrix<int>& permute_record) {
		// column permute the Matrix to be full rank on left
		bool is_ordered_structure = !(T() < T());

		Matrix<int> second_permutation(1, c);
		Matrix<bool> MRIPs_mark(1, c, '0');
		for (int i = 0; i < r; ++i) {
			int j;
			if (is_ordered_structure == false) {
				for (j = i; j < c && (*this)(i, j) == 0; ++j);	// remember to use strict equal only for finite fields
			}
			else {
				for (j = i; j < c && (double)my::abs((*this)(i, j)) < my::zero_approximation; ++j);	// for real or complex field
			}

			// strictly equals to 0, okay with GF2

			MRIPs_mark[j] = 1;
		}


		////cout << "MRIPs_mark" << MRIPs_mark << endl;

		//// This ensures both the MRIP part and LRP part are ordered in decreasing reliability

		int index_cnt = 0;
		for (int i = 0; i < c; ++i) {
			if (MRIPs_mark[i] == 1) {
				second_permutation[index_cnt] = i;
				index_cnt++;
				if (index_cnt == r)
					break;
			}
		}
		for (int i = 0; i < c; ++i) {
			if (MRIPs_mark[i] == 0) {
				second_permutation[index_cnt] = i;
				index_cnt++;
			}
		}

		//cout << "second_permutation" << second_permutation;

		permute_col(second_permutation);
		permute_record.permute(second_permutation);
	}

	/**
	 * .perform Gaussian-Jordan Elimination to make identity as left as possible
	 * we donot need to consider row permutation record, or singular case
	 * (*this) should be fat matrix and full rank of rows
	 *
	 * \param permute_record: permutation will be mirrored into permute_record
	 * \param num_of_column_processed: record number of columns traversed, which is the number of sequential steps
	 *
	 * \return Marix after Gaussian-Jordan Elimination
	 */
	void GJE_4_GF2_left(Matrix<int>& permute_record) {

		// row transformation to turn the matrix into an up-triangle Matrix
		int column_now = 0;
		for (int row_now = 0; row_now < r && column_now < c; ++row_now, ++column_now) {
			GF2_auxiliary_storage::iteration_number++;

			if ((*this)(row_now, column_now) == 0) {
				// find the row below row_now that have entry 1
				int nz_row_below;
				for (nz_row_below = row_now + 1; nz_row_below < r && (*this)(nz_row_below, column_now) == 0; ++nz_row_below);

				if (nz_row_below != r) {
					// in this case we can permute nz_row_below with row_now to get entry 1 in position (row_now, column_now)
					switch_row(nz_row_below, row_now);
				}
				else {
					// in this case the column column_now is all zero and we have to seek another column
					row_now--;			// keep row_now unchanged and move to the next column
					continue;
				}
			}

			// here, we can use row_now to eleminate other rows of the Matrix		

			// eliminate rows above and below the pivot row
			for (int i = 0; i < r; ++i) {
				if (i != row_now && (*this)(i, column_now) == 1) {

					GF2_auxiliary_storage::GE_bit_plane_number++;
					GF2_auxiliary_storage::GE_bit_plane_norm_number += (double)(c - column_now) / (c - r);

					int i_row_start_ind = i * c;
					int row_now_row_start_ind = row_now * c;					// optimize only in cubic for loop

					for (int j = column_now; j < c; ++j) {
						(*this)[i_row_start_ind + j] += (*this)[row_now_row_start_ind + j];
					}
				}
			}
		}
		
		col_permute_up_triangle_to_full_rank_left_v2(permute_record);

	}

	/**
	 * .this cause more operations than non-v2 version, for examing all rows
	 *
	 * \return ret/2 = number of rows shifted back, ret%2 = 0 even valid permutation, = 1 odd valid permutation
	 */
	int row_transformation_to_up_triangle_v2() {
		int permute_row_ind;
		int col_pos;
		int col_pos2;
		T temp;
		int ic, ir = 0;
		int r_shift = 0;
		int sign = 0;
		bool is_ordered_structure = !(T() < T());
		for (int i = 0; i + r_shift < r && i < c; i++) {
			//cout << "*this" << *this;
			ir = i + r_shift;
			ic = i;

			if (is_ordered_structure) {
				temp = max_abs_col_ele_start_from(ic, ir, permute_row_ind);	// for computation accuracy

				// Note: approach my::zero_approximation as 0
				if ((double)my::abs(temp) >= my::zero_approximation);							// this will affect the rank
				else {
					r_shift--;
					continue;
				}
			}
			else {
				temp = non_0_col_ele_start_from(ic, ir, permute_row_ind);	// specifically for GF

				if (temp != T(0));
				else {
					r_shift--;
					continue;
				}
			}

			//cout << "ir=" << ir << "\t" << "ic=" << ic << "\t" << "permute_row_ind=" << permute_row_ind << endl;
			if (ir != permute_row_ind) {
				switch_row(ir, permute_row_ind);
				sign = sign ^ 1;
			}

			col_pos = ir * c + ic;
			col_pos2 = col_pos;
			for (int j = ir + 1; j < r; j++) {
				col_pos2 += c;

				//if (e[col_pos2] != T(0)) {		// in Hessenberg matrix this can skip, and faster

				temp = e[col_pos2] / e[col_pos];
				for (int k = ic; k < c; ++k) {
					(*this)(j, k) -= (*this)(ir, k) * temp;
				}

				//}

			}
		}
		//cout << "*this" << *this;
		return sign + (-r_shift) * 2;
	}

	/**
	 * .this cause more operations than non-v2 version, for examing all rows
	 */
	void row_transformation_left_up_triangle_to_identity_v2() {
		// this function should be called after col_permute_to_full_rank_on_left
		for (int i = r - 1; i >= 0; --i) {	// diag ind
			// always (*this)(i, i) != 0
			if ((*this)(i, i) != T(1)) {			// scale row, not enter in GF2
				T factor = (*this)(i, i);
				for (int j = r; j < c; ++j) {	// col ind, knowing that col ind from i+1 to r-1 is 0 already
					(*this)(i, j) /= factor;
				}
				(*this)(i, i) /= factor;
			}
			for (int j = 0; j < i; ++j) {		// row ind

				if ((*this)(j, i) != T(0)) {		// in Tridiagonal matrix this can skip, and faster
					
					// row minus
					T factor = (*this)(j, i);		// bull-shift for matching with jiabao's complexity
					for (int k = r; k < c; ++k) {	// col ind
						(*this)(j, k) -= factor * (*this)(i, k);
					}
					(*this)(j, i) -= factor;

				}
				//cout << "G in =" << (*this);
			}
		}
	}

	/**
	 * .Gaussing elimination on (*this) to transform it in [I|P] form, 
	 *  with column permutation recorded as 'permute_record'
	 * 
	 * \return This is an inplace function
	 */
	void GE_left_identity_match_jiabao(Matrix<int>& permute_record) {
		row_transformation_to_up_triangle_v2();

		col_permute_up_triangle_to_full_rank_left_v2(permute_record);

		row_transformation_left_up_triangle_to_identity_v2();
	}

	/* --------- re-write the GE and GJE: right identity ------------- */

	/**
	 * . this step is after row transformation to up triangle
	 *	 permute the column of (*this) to get an left up triangle matrix,
	 *   not disturbing the original column order as much as possible
	 *
	 * \return
	 */
	void col_permute_low_triangle_to_full_rank_right_v2(Matrix<int>& permute_record) {

		// column permute the Matrix to be full rank on right

		Matrix<int> second_permutation(1, c);
		Matrix<bool> MRIPs_mark(1, c, '0');
		for (int i = r - 1; i >= 0; --i) {
			int j;
			for (j = c - (r - i); j >= 0 && (*this)(i, j) == 0; --j);	// remember to use strict equal only

			// strictly equals to 0, okay with GF2

			MRIPs_mark[j] = 1;
		}

		//cout << "MRIPs_mark" << MRIPs_mark << endl;

		//// This ensures both the MRIP part and LRP part are ordered in decreasing reliability

		int index_cnt = 0;
		int k = c - r;
		for (int i = 0; i < c; ++i) {
			if (MRIPs_mark[i] == 0) {
				second_permutation[index_cnt] = i;
				index_cnt++;
				if (index_cnt == k)
					break;
			}
		}
		for (int i = 0; i < c; ++i) {
			if (MRIPs_mark[i] == 1) {
				second_permutation[index_cnt] = i;
				index_cnt++;
			}
		}

		//cout << "second_permutation" << second_permutation;

		permute_col(second_permutation);
		permute_record.permute(second_permutation);
	}

	/**
	 * .perform Gaussian-Jordan Elimination to make identity as left as possible
	 * we donot need to consider row permutation record, or singular case
	 * (*this) should be fat matrix and full rank of rows
	 *
	 * \param permute_record: permutation will be mirrored into permute_record
	 * \param num_of_column_processed: record number of columns traversed, which is the number of sequential steps
	 *
	 * \return Marix after Gaussian-Jordan Elimination
	 */
	void GJE_4_GF2_right(Matrix<int>& permute_record) {

		// row transformation to turn the matrix into an up-triangle Matrix
		int column_now = c - 1;
		for (int row_now = r - 1; row_now >= 0 && column_now >= 0; --row_now, --column_now) {
			GF2_auxiliary_storage::iteration_number++;

			if ((*this)(row_now, column_now) == 0) {
				// find the row above row_now that have entry 1
				int nz_row_above;
				for (nz_row_above = row_now - 1; nz_row_above >= 0 && (*this)(nz_row_above, column_now) == 0; --nz_row_above);

				if (nz_row_above != -1) {
					// in this case we can permute nz_row_above with row_now to get entry 1 in position (row_now, column_now)
					switch_row(nz_row_above, row_now);
				}
				else {
					// in this case the column column_now is all zero and we have to seek another column
					row_now++;			// keep row_now unchanged and move to the next column
					continue;
				}
			}

			// here, we can use row_now to eleminate other rows of the Matrix		

			// eliminate rows above and below the pivot row
			for (int i = 0; i < r; ++i) {
				if (i != row_now && (*this)(i, column_now) == 1) {

					GF2_auxiliary_storage::GE_bit_plane_number++;
					GF2_auxiliary_storage::GE_bit_plane_norm_number += (column_now + 1.0) / r;

					int i_row_start_ind = i * c;
					int row_now_row_start_ind = row_now * c;					// optimize only in cubic for loop

					for (int j = column_now; j >= 0; --j) {
						(*this)[i_row_start_ind + j] += (*this)[row_now_row_start_ind + j];
					}
				}
			}
		}
		//num_of_column_processed = column_now;

		//cout << "(*this)" << (*this);

		col_permute_low_triangle_to_full_rank_right_v2(permute_record);
	}

	/**
	 * .perform Gaussian-Jordan Elimination to make identity as left as possible
	 *	we donot need to consider row permutation record, or singular case
	 *	(*this) should be fat matrix and full rank of rows
	 *	worse than 'GJE_4_GF2_right' incuring *, and checking all rows
	 *
	 * \param permute_record: permutation will be mirrored into permute_record
	 * \param num_of_column_processed: record number of columns traversed, which is the number of sequential steps
	 *
	 * \return Marix after Gaussian-Jordan Elimination
	 */
	void GJE_right_v2(Matrix<int>& permute_record) {

		// row transformation to turn the matrix into an up-triangle Matrix
		int column_now = c - 1;
		for (int row_now = r - 1; row_now >= 0 && column_now >= 0; --row_now, --column_now) {

			if ((*this)(row_now, column_now) == 0) {
				// find the row above row_now that have entry 1
				int nz_row_above;
				for (nz_row_above = row_now - 1; nz_row_above >= 0 && (*this)(nz_row_above, column_now) == 0; --nz_row_above);

				if (nz_row_above != -1) {
					// in this case we can permute nz_row_above with row_now to get entry 1 in position (row_now, column_now)
					switch_row(nz_row_above, row_now);
				}
				else {
					// in this case the column column_now is all zero and we have to seek another column
					row_now++;			// keep row_now unchanged and move to the next column
					continue;
				}
			}

			// here, we can use row_now to eleminate other rows of the Matrix		

			// eliminate rows above and below the pivot row
			for (int i = 0; i < r; ++i) {
				//if (i != row_now && (*this)(i, column_now) == 1) {
				//	int i_row_start_ind = i * c;
				//	int row_now_row_start_ind = row_now * c;					// optimized
				//	for (int j = column_now; j >= 0; --j) {
				//		(*this)[i_row_start_ind + j] -= (*this)[row_now_row_start_ind + j];
				//	}
				//}

				if (i != row_now) {
					T factor = (*this)(i, column_now);
					int i_row_start_ind = i * c;
					int row_now_row_start_ind = row_now * c;					// dis-optimize
					for (int j = column_now; j >= 0; --j) {
						(*this)[i_row_start_ind + j] -= (*this)[row_now_row_start_ind + j] * factor;
					}
				}
			}
		}
		//num_of_column_processed = column_now;

		//cout << "(*this)" << (*this);

		col_permute_low_triangle_to_full_rank_right_v2(permute_record);
	}

	/**
	 * .this cause more operations than non-v2 version, for examing all rows
	 *
	 * \return ret/2 = number of rows shifted forward, ret%2 = 0 even valid permutation, = 1 odd valid permutation
	 */
	int row_transformation_to_low_triangle_v2() {

		int permute_row_ind;
		int col_pos;
		int col_pos2;
		T temp;
		int ic = c - 1, ir = r - 1;
		int r_shift = 0;
		int sign = 0;
		bool is_ordered_structure = !(T() < T());
		for (; ir >= 0 && ic >= 0; ic--) {
			//cout << "*this" << *this;

			ir = ic - c + r + r_shift;

			if (is_ordered_structure) {
				temp = max_abs_col_ele_end_from(ic, ir, permute_row_ind);	// for computation accuracy

				// Note: approach my::zero_approximation as 0
				if ((double)my::abs(temp) >= my::zero_approximation);					// this will affect the rank
				else {
					r_shift++;
					continue;
				}
			}
			else {
				temp = non_0_col_ele_end_from(ic, ir, permute_row_ind);		// specifically for GF

				if (temp != T(0));
				else {
					r_shift++;
					continue;
				}
			}

			//cout << "ir=" << ir << "\t" << "ic=" << ic << "\t" << "permute_row_ind=" << permute_row_ind << endl;
			if (ir != permute_row_ind) {
				switch_row(ir, permute_row_ind);
				sign = sign ^ 1;
			}

			col_pos = ir * c + ic;
			col_pos2 = col_pos;
			for (int j = ir - 1; j >= 0; j--) {
				col_pos2 -= c;
				//if (e[col_pos2] != 0) {		// in Hessenberg matrix this can skip, and faster

				temp = e[col_pos2] / e[col_pos];
				for (int k = ic; k >= 0; --k) {
					(*this)(j, k) -= (*this)(ir, k) * temp;
				}

				//}
			}
		}
		//cout << "*this" << *this;
		return sign + r_shift * 2;
	}

	/**
	 * .this cause more operations than non-v2 version, for examing all rows
	 */
	void row_transformation_right_low_triangle_to_identity_v2() {
		// this function should be called after col_permute_to_full_rank_on_left
		for (int i = 0; i < r; ++i) {	// row ind
			int col_ind = c - r + i;
			// always (*this)(i, col_ind) != 0
			if ((*this)(i, col_ind) != 1) {
				T factor = (*this)(i, col_ind);
				for (int j = c - r - 1; j >= 0; --j) {	// col ind, knowing that col ind from c-r to i-1 is 0 already
					(*this)(i, j) /= factor;
				}
				(*this)(i, col_ind) /= factor;
			}
			for (int j = r - 1; j > i; --j) {		// row ind
				if ((*this)(j, col_ind) != 0) {		// in Tridiagonal matrix this can skip, and faster
					
					// row minus
					T factor = (*this)(j, col_ind);	// this wirting is bull-shit, just for matching jiabao's complexity
					for (int k = c - r - 1; k >= 0; --k) {	// col ind
						(*this)(j, k) -= factor * (*this)(i, k);
					}
					(*this)(j, col_ind) -= factor;

				}
				//cout << "G in =" << (*this);
			}
		}
	}

	/**
	 * .Gaussing elimination on (*this) to transform it in [I|P] form,
	 *  with column permutation recorded as 'permute_record'
	 *
	 * \return This is an inplace function
	 */
	void GE_right_identity_match_jiabao(Matrix<int>& permute_record) {
		row_transformation_to_low_triangle_v2();

		col_permute_low_triangle_to_full_rank_right_v2(permute_record);

		row_transformation_right_low_triangle_to_identity_v2();
	}

	/* --------- extra functions ------------- */

	void erase_element(int ind) {
		// the matrix must be transformed to a vector, i.e., r = 1

		for (int i = ind + 1; i < c; ++i) {
			e[i - 1] = e[i];
		}
		if (c > 0) {
			c--;		// pop back
		}
	}

	Matrix<T> erase_all_0_rows() const {

		Matrix<int> all_zero_row_ind(1, r, 'v');
		for (int i = 0; i < r; ++i) {
			bool is_all_zero_row = true;
			for (int j = 0; j < c; ++j) {
				if ((*this)(i, j) != 0) {
					is_all_zero_row = false;
					break;
				}
			}

			if (is_all_zero_row == true) {
				all_zero_row_ind.push_back(i);
			}
		}
		//cout << "all_zero_row_ind" << all_zero_row_ind;

		Matrix<T> ans = erase_rows(all_zero_row_ind, true);
		return ans;
	}
};

/**
 * .the sorted vector, with each element contaning index int sorted from small to big, and the value
 * this class is sutible for vector that donot need to change after its definieion
 * 
 * this is a little complicated, not easy to use, it's the problem of function interface
 */
template<class T>
class Sorted_vector {
public:
	vector<pair<int, T>> v;
	Sorted_vector(int len = 0) {
		// you have to assign enough space for the sorted v at the very begining
		v = vector<pair<int, T>>(len);		// the len cannot be zero
	}
	void assign_ind(int index, int key) {
		// please assign the element in the vector in an ordered way
		// i.e., if position_1 < position_2, then index_1 < index_2
		v[index].first = key;

		//cout << "v[" << index << "].first=" << v[index].first << endl;

		// we donot recommand change the class during the use
		// since it will likely voilate the structure or induce high complexity
		// so you have to defined it at the begining, and later visit it only
	}
	int size() const {
		return (int)v.size();
	}
	// access the individual element, given key, return value
	T& operator[](const int& key) {
		// for a non existing key, it may return some value, please make sure the key exists
		
		// use the binary search method
		int start = 0;
		int end = (int)v.size() - 1;

		// binary search over v, but time consuming
		while (start < end) {
			int mid = (start + end) >> 1;
			if (v[mid].first < key) {
				start = mid + 1;
			}
			else {
				end = mid;
			}
		}
		return v[start].second;
		
	}
	// return the index corresponding to key, -1 indecates key not found
	int find(const int key) const {
		// use the binary search method
		int start = 0;
		int end = (int)v.size() - 1;

		//int mid = (start + end) >> 1;
		// binary search over v, O(log(end)) comparison and program jump
		while (start < end) {		// the standard binary search template
			int mid = (start + end) >> 1;
			if (v[mid].first < key) {
				start = mid + 1;
			}
			else {
				end = mid;
			}
		}
		return v[start].first == key ? start : -1;
	}
	// key corresponding to index
	const int& key(const int& index) const {
		return v[index].first;
	}
	// value corresponding to index
	T& val(const int& index) {		// val can be changed
		return v[index].second;
	}
	friend ostream& operator << (ostream& out, const Sorted_vector<T>& sv) {
		int svs = (int)sv.v.size();
		out << endl;
		for (int i = 0; i < svs; ++i) {
			out << i << ": " << sv.v[i].first << " -> " << sv.v[i].second << endl;
		}
		return out;
	}
};

/**
 * .flexible Matrix, i.e., Matrix's rows may have different column size
 */
template<class T>
class Matrix_flex_col
{
protected:
	int r;
	int* c;		// size of r + 1, starting index for each row, the last index is size of the Matrix
	T* e;		// we donot allow change the size of the Matrix during processing

	// note that copy constructor and operator = is not allowed
	Matrix_flex_col(const Matrix_flex_col<T>& A) {}
	void operator = (const Matrix_flex_col<T>& A) {}
public:

	//constructors and destructor,
	Matrix_flex_col():r(0) {
		c = NULL;
		e = NULL;
	}

	// please be very careful of defining column_size, which is column size of each row
	Matrix_flex_col(int _r, const Matrix<int>& column_size) {
		if (_r == column_size.size()) {
			r = _r;
			c = new int[r + 1];
			c[0] = 0;
			for (int i = 1; i <= r; ++i) {
				c[i] = c[i - 1] + column_size(i - 1);
			}
			e = new T[c[r]];
		}
		else {
			cout << "error in initializing Matrix_flexible_column, for _r != column_size.size()" << endl;
			cout << "_r = " << _r << "\tcolumn_size.size() = " << column_size.size() << endl;
			r = 0;
			c = NULL;
			e = NULL;
		}
	}

	~Matrix_flex_col() {
		// before this, destructor of e[i] has already been called, 
		// hence it will be okay to destroy all the element, no memory waste
		delete[] c;
		delete[] e;
	}

	// set size, not recommand to use repeatedly
	void resize(const Matrix<int>& column_size) {

		delete[] c;
		delete[] e;		// you can delete first, no original data kept

		r = column_size.size();
		c = new int[r + 1];
		c[0] = 0;
		for (int i = 1; i <= r; ++i) {
			c[i] = c[i - 1] + column_size(i - 1);
		}
		e = new T[c[r]];
	}

	// set all element in the matrix by default
	void reset(const T& elements = T()) {
		int sss = c[r];
		for (int i = 0; i < sss; ++i) {
			e[i] = elements;
		}
	}

	// Access the individual elements, we define both [] operator and () operator, for convinience
	inline T& operator[](const unsigned& pos_ind) {
		return e[pos_ind];
	}
	inline const T& operator[](const unsigned& pos_ind) const {
		return e[pos_ind];
	}
	inline T& operator()(const unsigned& pos_ind) {
		return e[pos_ind];
	}
	inline const T& operator()(const unsigned& pos_ind) const {
		return e[pos_ind];
	}
	inline T& operator()(const unsigned& row_ind, const unsigned& col_ind) {
		return e[c[row_ind] + col_ind];
	}
	inline const T& operator()(const unsigned& row_ind, const unsigned& col_ind) const {
		return e[c[row_ind] + col_ind];
	}

	// return the number of row
	inline int row() const {
		return r;
	}
	// return the number of column at specific row
	inline int col(int r_ind) const {
		return c[r_ind + 1] - c[r_ind];
	}
	// return the number of all element
	inline int size() const {
		return c[r];
	}

	// otuput
	friend ostream& operator << (ostream& out, const Matrix_flex_col<T>& A) {
		out << " (row size: " << A.r << ")" << endl;
		out << setprecision(4);
		for (int i = 0; i < A.r; i++) {
			out << " (row: " << i << ", size: " << A.c[i + 1] - A.c[i] << ")" << endl;
			for (int j = A.c[i]; j < A.c[i + 1]; j++) {
				out << setw(14) << A.e[j];
			}
			out << endl;
		}
		out << setprecision(6);
		return out;
	}
};

template<class T>
class Heap_max :public Matrix<T>
{
public:
	void ascend_to_father(int child_ind) {
		while (child_ind != 0) {
			int father_ind = (child_ind - 1) >> 1;
			if (this->e[child_ind] > this->e[father_ind]) {
				this->switch_ele(child_ind, father_ind);
				child_ind = father_ind;
			}
			else {
				return;
			}
		}
	}
	void build() {
		int half_size = this->size() >> 1;
		for (int father_ind = half_size - 1; father_ind >= 0; --father_ind) {
			descend_to_child(father_ind);
		}
	}

	/**
	 *
	 * \return if descending happened
	 */
	bool descend_to_child(int father_ind) {
		int max_ind = this->size() - 1;
		int can_child_ind = (father_ind << 1) + 1;
		int orig_father_ind = father_ind;
		while (can_child_ind < max_ind) {
			// if the father_ind has 2 children, enter while loop

			can_child_ind += this->e[can_child_ind] < this->e[can_child_ind + 1];
			if (this->e[father_ind] < this->e[can_child_ind]) {
				this->switch_ele(father_ind, can_child_ind);
				father_ind = can_child_ind;
				can_child_ind = (father_ind << 1) + 1;
			}
			else {
				// father stop descending
				return orig_father_ind != father_ind;
			}
		}

		if (can_child_ind == max_ind && this->e[father_ind] < this->e[can_child_ind]) {
			// have one child only
			this->switch_ele(father_ind, can_child_ind);
			return true;
		}
		return orig_father_ind != father_ind;
	}
	void pop(int pop_ind = 0) {
		int last_ind = this->size() - 1;
		if (pop_ind < last_ind) {
			this->switch_ele(last_ind, pop_ind);
			this->pop_back();
			bool is_descended = descend_to_child(pop_ind);
			if (is_descended);
			else {
				ascend_to_father(pop_ind);
			}
		}
		else {
			// poping back donot change the property of heap
			this->pop_back();
		}
	}
	void pop_leaf(int pop_ind) {
		int last_ind = this->size() - 1;
		if (pop_ind < last_ind) {
			this->switch_ele(last_ind, pop_ind);
			this->pop_back();
			ascend_to_father(pop_ind);
		}
		else {
			this->pop_back();
		}
	}
	void push(const T& val) {
		this->push_back(val);
		ascend_to_father(this->size() - 1);
	}
	void replace_top(const T& val) {
		this->e[0] = val;
		descend_to_child(0);
	}
	int bottom() const {
		// find the index of min element on the heap
		int sss = this->size();
		int ans = sss >> 1;
		for (int i = ans + 1; i < sss; ++i) {
			ans = this->e[i] < this->e[ans] ? i : ans;
		}
		return ans;
	}
	bool is_heap() const {
		int sss = this->size();
		if (sss == 0)	return true;

		int father_ind = 0;
		while ((father_ind << 1) + 2 < sss) {
			// this father has two child
			if (this->e[father_ind] < this->e[(father_ind << 1) + 1] || this->e[father_ind] < this->e[(father_ind << 1) + 2]){
				return false;
			}
			father_ind++;
		}
		if ((father_ind << 1) + 1 < sss && this->e[father_ind] < this->e[(father_ind << 1) + 1]) {
			return false;
		}
		return true;

	}
	inline T top() const {
		return this->e[0];
	}
	void replace(int ind, const T& val) {
		if (val < this->e[ind]) {		// special for max heap
			this->e[ind] = val;
			descend_to_child(ind);
		}
		else {
			this->e[ind] = val;
			ascend_to_father(ind);
		}
	}
};

template<class T>
class Heap_min :public Matrix<T>
{
public:
	void ascend_to_father(int child_ind) {
		while (child_ind != 0) {
			int father_ind = (child_ind - 1) >> 1;
			if (this->e[child_ind] < this->e[father_ind]) {
				this->switch_ele(child_ind, father_ind);
				child_ind = father_ind;
			}
			else {
				return;
			}
		}
	}
	void build() {
		int half_size = this->size() >> 1;
		for (int father_ind = half_size - 1; father_ind >= 0; --father_ind) {
			descend_to_child(father_ind);
		}
	}

	/**
	 *
	 * \return if descending happened
	 */
	bool descend_to_child(int father_ind) {
		int max_ind = this->size() - 1;
		int can_child_ind = (father_ind << 1) + 1;
		int orig_father_ind = father_ind;
		while (can_child_ind < max_ind) {
			// if the father_ind has 2 children, enter while loop

			can_child_ind += this->e[can_child_ind] > this->e[can_child_ind + 1];
			if (this->e[father_ind] > this->e[can_child_ind]) {
				this->switch_ele(father_ind, can_child_ind);
				father_ind = can_child_ind;
				can_child_ind = (father_ind << 1) + 1;
			}
			else {
				// father stop descending
				return orig_father_ind != father_ind;
			}
		}

		if (can_child_ind == max_ind && this->e[father_ind] > this->e[can_child_ind]) {
			// have one child only
			this->switch_ele(father_ind, can_child_ind);
			return true;
		}
		return orig_father_ind != father_ind;
	}
	void pop(int pop_ind = 0) {
		this->switch_ele(this->size() - 1, pop_ind);
		this->pop_back();
		bool is_descended = descend_to_child(pop_ind);
		if (is_descended);
		else {
			ascend_to_father(pop_ind);
		}
	}
	void pop_leaf(int pop_ind) {
		int last_ind = this->size() - 1;
		if (pop_ind < last_ind) {
			this->switch_ele(last_ind, pop_ind);
			this->pop_back();
			ascend_to_father(pop_ind);
		}
		else {
			this->pop_back();
		}
	}
	void replace_top(const T& val) {
		this->e[0] = val;
		descend_to_child(0);
	}
	void push(const T& val) {
		this->push_back(val);
		ascend_to_father(this->size() - 1);
	}
	int bottom() const {
		// find the index of max element on the heap
		int sss = this->size();
		int ans = sss >> 1;
		for (int i = ans + 1; i < sss; ++i) {
			ans = this->e[i] > this->e[ans] ? i : ans;
		}
		return ans;
	}
	bool is_heap() const {
		int sss = this->size();
		if (sss == 0)	return true;

		int father_ind = 0;
		while ((father_ind << 1) + 2 < sss) {
			// this father has two child
			if (this->e[father_ind] > this->e[(father_ind << 1) + 1] || this->e[father_ind] > this->e[(father_ind << 1) + 2]){
				return false;
			}
			father_ind++;
		}
		if ((father_ind << 1) + 1 < sss && this->e[father_ind] > this->e[(father_ind << 1) + 1]) {
			return false;
		}
		return true;

	}
	inline T top() const {
		return this->e[0];
	}
	void replace(int ind, const T& val) {
		if (val > this->e[ind]) {		// special for min heap
			this->e[ind] = val;
			descend_to_child(ind);
		}
		else {
			this->e[ind] = val;
			ascend_to_father(ind);
		}
	}
};

// a singly-linked List in continous memory, and you can only remove contents once the list is built
template<class T>
class sList : public Matrix<T>
{
protected:
	int* next_ind;
	int start_ind;
	int current_ind;		// actually single linked list is enough for GE_PSM data structure.
public:

	void resize(int row, int column, bool keep_original = false) {		// you have to call resize to iniitalize next_ind
		Matrix<T>::resize(row, column, keep_original);

		delete[] next_ind;
		next_ind = new int[Matrix<T>::size()];
	}
	void build() {
		
		int sss = this->size();
		// flat the matrix into a row vector
		Matrix<T>::resize(1, sss, true);

		// we only build non-empty matrix into list
		if (sss > 0) {

			// build the list link relation, from start to end connected one by one
			for (int i = 1; i < sss; ++i) {
				next_ind[i - 1] = i;
			}
			next_ind[sss - 1] = -1;				// the tail of list has no next node
			start_ind = 0;
			current_ind = 0;					// the current pointer always initialized to the head node index
		}
	}
	void move_start() {
		current_ind = start_ind;
	}
	void move_next() {
		if (current_ind != -1)
			current_ind = next_ind[current_ind];
	}
	void remove_start() {
		start_ind = next_ind[start_ind];
		this->c--;
	}
	void remove_next() {
		if (current_ind != -1 && next_ind[current_ind] != -1) {
			// current index is not the end of list, and also its next, then we can remove
			next_ind[current_ind] = next_ind[next_ind[current_ind]];
			this->c--;
		}
	}

	inline int get_start_ind() {
		return start_ind;
	}
	inline int get_current_ind() {
		return current_ind;
	}
	inline int get_next_ind() {
		if (current_ind != -1)
			return next_ind[current_ind];
		else
			return -1;
	}

	inline T& start_val() {
		return this->e[start_ind];
	}
	inline const T& start_val() const {
		return this->e[start_ind];
	}
	inline T& current_val() {
		return this->e[current_ind];
	}
	inline const T& current_val() const {
		return this->e[current_ind];
	}

	~sList() {
		delete[] next_ind;
	}
	friend ostream& operator << (ostream& out, const sList<T>& A) {
		out << " (size: " << A.c << ")" << endl;
		out << setprecision(4);
		bool is_ordered_structure = !(T() < T());
		for (int i = A.start_ind; i != -1; i = A.next_ind[i]) {
			if (is_ordered_structure /*&& A.e[row_pos + j] >= my::zero_approximation*/)
				out << setw(14) << A.e[i];
			/*else if(A.is_ordered_structure)
				out<< setw(14) << 0.000;*/
				// approach my::zero_approximation be 0 at output, may cause confusion
			else {
				//out << setw(6) << A.e[row_pos + j];	// for GF display, no decimal point

				// GF2 display
				out << setw(3) << A.e[i];

				// for display over lenovo notebook
				/*if ((j + 1) % 26 == 0) {
					out << " ";
				}*/
			}
		}
		out << endl;
		out << setprecision(6);
		return out;
	}
};


// a singly-linked List in continous memory, and you can only remove contents once the list is built
template<class T>
class dList : public Matrix<T>
{
protected:
	int* next_ind;
	int* former_ind;
	int start_ind;
	int current_ind;		// actually single linked list is enough for GE_PSM data structure.
public:
	// challange: add end_ind, move_former, insert, all in O(1) with enough space decleared, holes allowed

	void resize(int row, int column, bool keep_original = false) {		// you have to call resize to iniitalize next_ind
		Matrix<T>::resize(row, column, keep_original);
		int sss = this->size();

		delete[] next_ind;
		next_ind = new int[sss];

		delete[] former_ind;
		former_ind = new int[sss];		// simply delete the old one and build a new one, donot compare capacity here.
	}
	void build() {

		int sss = this->size();
		// flat the matrix into a row vector
		Matrix<T>::resize(1, sss, true);

		// we only build non-empty matrix into list
		if (sss > 0) {

			// build the list link relation, from start to end connected one by one
			former_ind[0] = -1;					// the head of list has no former node

			for (int i = 1; i < sss; ++i) {
				next_ind[i - 1] = i;
				former_ind[i] = i - 1;
			}
			next_ind[sss - 1] = -1;				// the tail of list has no next node
			start_ind = 0;
			current_ind = 0;					// the current pointer always initialized to the head node index
		}
	}
	void move_start() {
		current_ind = start_ind;
	}
	void move_next() {
		if (current_ind != -1)
			current_ind = next_ind[current_ind];
	}
	void remove() {
		if (current_ind != -1) {
			// current index is not the end of list, then we can remove

			if (former_ind[current_ind] != -1) {
				// link the former node to the next node
				next_ind[former_ind[current_ind]] = next_ind[current_ind];
			}
			else {
				// the head node is removed, update start_ind
				start_ind = next_ind[current_ind];
			}

			if (next_ind[current_ind] != -1) {
				// link the next node to the former node
				former_ind[next_ind[current_ind]] = former_ind[current_ind];
			}

			// move current_ind to next node's index
			current_ind = next_ind[current_ind];
			this->c--;
		}
	}	
	inline int get_start_ind() {
		return start_ind;
	}
	inline int get_current_ind() {
		return current_ind;
	}

	~dList() {
		delete[] next_ind;
		delete[] former_ind;
	}
	friend ostream& operator << (ostream& out, const dList<T>& A) {
		out << " (size: " << A.c << ")" << endl;
		out << setprecision(4);
		bool is_ordered_structure = !(T() < T());
		for (int i = A.start_ind; i != -1; i = A.next_ind[i]) {
			if (is_ordered_structure /*&& A.e[row_pos + j] >= my::zero_approximation*/)
				out << setw(14) << A.e[i];
			/*else if(A.is_ordered_structure)
				out<< setw(14) << 0.000;*/
				// approach my::zero_approximation be 0 at output, may cause confusion
			else {
				//out << setw(6) << A.e[row_pos + j];	// for GF display, no decimal point

				// GF2 display
				out << setw(3) << A.e[i];

				// for display over lenovo notebook
				/*if ((j + 1) % 26 == 0) {
					out << " ";
				}*/
			}
		}
		out << endl;
		out << setprecision(6);
		return out;
	}
};