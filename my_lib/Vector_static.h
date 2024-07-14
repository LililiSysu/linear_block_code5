#pragma once
/*****************************************************************//**
 * \file	Vector_static.h
 * \brief	Self-defined class of vector, for a bunch of function of my need,
 *			but the function interface is similar to class Matrix in Matrix.h.
 *			The calss my_vector is static, not size extension allowed.
 *			Please allocate enough space at the compile state of the program.
 *			This make the porgram simpler and faster.
 * 
 * \todo	(1) implement Mat_s and all the functions associated with Matrix
 * 
 *			(2) implement deQ_s by inheriting class Vec_s
 * 
 *			(3) and implement the insert, move_former, move_end for dList_s, all in O(1) time cost, 
 *				adjust the inteface in consistance with those in Matrix.h
 * 
 * 
 * \author	26259
 * \date	July 2023
 *********************************************************************/

#include"my.h"
#include<iostream>
#include<iomanip>
#include<vector>
using namespace std;


// this class is just for storage and fetching, not for computing, we assume the type donot have any opearator defined yet.
template<class T, int cap>
class Vec_s
{
protected:
	T e[cap];		// will call constructor of e automatically
	int dim;
public:

	//constructors
	Vec_s() :dim(0), e() {};
	// accept one parameter and construct a matrix with size 1 by 1, for convinience only
	Vec_s(int _dim) :dim(_dim), e() {}
	/**
	 * .
	 * \param v:	the vector to copy, NOTE if using {} to replace 'v', then {}
	 *				with 1 element like {5}, is not allowed and this will set matrix with 5 0s
	 */
	Vec_s(const vector<T>& v){
		// not that this function is not valid for v={1}, that is v set with {} and has 1 zero element

		resize((int)v.size());
		for (int i = 0; i < dim; ++i) {
			// make sure not overflow
			e[i] = v[i];
		}
	}	

	// set all element in the matrix be 0
	inline void reset(const T& elements = T()) {
		for (int i = 0; i < dim; ++i) {
			e[i] = elements;
		}
	}
	inline int size() const {
		return dim;
	}
	inline int space() const {
		return cap;		// you have to declear size for each my_vector
	}
	inline bool empty() const {
		return dim == 0;
	}

	// interface same as vector
	inline void push_back(const T& element) {
#ifdef DEBUG
		if (dim < cap) {
#endif
			e[dim] = element;
			dim++;
#ifdef DEBUG
		}
		else {
			cout << "warning: [push_back] overflow in [my_vector]" << endl;
		}
#endif
	}
	inline void pop_back() {
#ifdef DEBUG
		if (dim > 0) {
#endif
			dim--;
#ifdef DEBUG
		}
		else {
			cout << "warning: [pop_back] underflow in [my_vector]" << endl;
		}
#endif
	}
	inline T& back() {
#ifdef DEBUG
		if (dim > 0) {
#endif
			return e[dim - 1];
#ifdef DEBUG
		}
		else {
			cout << "warning: [back] underflow in [my_vector], no element at back" << endl;
			return e[0];
		}
#endif
	}
	inline const T& back() const {
#ifdef DEBUG
		if (dim > 0) {
#endif
			return e[dim - 1];
#ifdef DEBUG
		}
		else {
			cout << "warning: [back] underflow in [my_vector], no element at back" << endl;
			return e[0];
		}
#endif
	}

	// reset size of my_vector, very important function
	inline void resize(int _dim) {
		dim = my::min(cap, _dim);		// if _dim > cap, something unexpected will happen
	}
	inline void clear() {
		dim = 0;
	}
	inline void switch_ele(int position1_ind, int position2_ind) {
		//if (position_1 == position_2) return;
		T temp = e[position1_ind];
		e[position1_ind] = e[position2_ind];
		e[position2_ind] = temp;
	}
	/**
	 * .reverse the Vec_s
	 *
	 */
	void rev() {
		int len = size() - 1;
		int len_over_2 = len / 2;
		for (int i = 0; i <= len_over_2; ++i) {
			switch_ele(i, len - i);
		}
	}
	
	// Donot use 'ans' to call the following function, if you have to, decleare a new Vec_s<T, cap> to store 'ans'
	// you have to make sure ans_static not used in the future, since it is corrupted here
	
	// the inplace operation of keeping part of (*this)
	void keep_part(int start_ind, int end_ind) {

		// view -1 as the last index
		start_ind = start_ind == -1 ? dim - 1 : start_ind;
		end_ind = end_ind == -1 ? dim - 1 : end_ind;

		resize(end_ind - start_ind + 1);
		if (start_ind != 0) {
			for (int j = start_ind, k = 0; k < dim; j++, k++) {
				e[k] = e[j];		// the copy is free of interact, hence &ans_static == this will be okay
			}
		}
	}

	/**
	 * .				put the element back as v_ind indicates, inverse of permute
	 *
	 * \param v_ind:	the element of position i will occupy position v_ind[i]
	 *					e.g., {a,b,c}.permute_back({2,0,1}) gives {b,c,a}
	 * 
	 */
	void permute_back(const Vec_s<int, cap>& v_ind) {

		Vec_s<T, cap> ans_static(*this);	// copying the original part
		for (int i = 0; i < dim; ++i) {
			e[v_ind[i]] = ans_static[i];			// partial permute back for first part
		}
	}
	/**
	 * .				put the element as v_ind indicates, inverse of permute_back
	 *
	 * \param v_ind:	the element of position v_ind[i] will occupy position i
	 *					e.g., {a,b,c}.permute({2,0,1}) gives {c,a,b}
	 * 
	 */
	void permute(const Vec_s<int, cap>& v_ind) {
		Vec_s<T, cap> ans_static(*this);
		for (int i = 0; i < dim; ++i) {
			e[i] = ans_static[v_ind[i]];		// partial permute for the first part
		}
	}
	/**
	 * .randomly permute the elements of Vec_s
	 */
	void permute_rand() {
		Vec_s<int, cap> v_ind(dim);
		for (int i = 0; i < dim; ++i) {
			v_ind[i] = i;
		}
		v_ind.keep_random_element(dim);
		permute(v_ind);
	}
	/**
	 * .keep num elements from (*this) Vec_s, not overlapping
	 */
	void keep_random_element(int num) {
		int sss = size();
		if (num <= sss) {
			resize(0);
			for (int i = 0; i < num; ++i) {
				int choose_ind = my::rand_int_adv(i, sss - 1);
				switch_ele(dim, choose_ind);
				dim++;
			}
		}
		else {
			cout << "(get_random_element) warning: num > size()" << endl;
			cout << "num = " << num << ", size() = " << size() << endl;
		}
	}
	/**
	 * .		shift the columns of matrix circularly to the right
	 *			maybe we can consider pass a reference as parameter, reducing memory declaration and destruction
	 *
	 * \param	shift_value: the number of columns to shift, positive stands for shift right, negtive stands for shift left
	 *				e.g., (*this) is [a,b,c,d], and shift_value is 1, the result vector is [d,a,b,c]
	 *				e.g., (*this) is [a,b,c,d], and shift_value is -1, the result vector is [b,c,d,a]
	 *
	 */
	void shift_right_cir(int shift_value) {
		Vec_s<T, cap> ans_static(*this);

		if (shift_value >= 0) {		// shift right circularly

			// copy the columns before shift_value
			for (int j_ans = 0, j_ans_max = shift_value, j_this = dim - shift_value; j_ans < j_ans_max; ++j_ans, ++j_this) {

				e[j_ans] = ans_static[j_this];		// I believe this is the best program without parallelization
				// never use c function memcpy in c++
			}

			// copy the columns after shift_value
			for (int j_ans = shift_value, j_ans_max = dim, j_this = 0; j_ans < j_ans_max; ++j_ans, ++j_this) {

				e[j_ans] = ans_static[j_this];
			}
		}
		else {					// shift left circularly
			shift_value = -shift_value;

			// copy the columns before shift_value
			for (int j_ans = dim - shift_value, j_ans_max = dim, j_this = 0; j_ans < j_ans_max; ++j_ans, ++j_this) {

				e[j_ans] = ans_static[j_this];
			}

			// copy the columns after shift_value
			for (int j_ans = 0, j_this = shift_value, j_this_max = dim;	j_this < j_this_max; ++j_ans, ++j_this) {

				e[j_ans] = ans_static[j_this];
			}
		}
	}

	// Access the individual elements, we define both [] only
	inline T& operator[](const unsigned& pos_ind) {
		return e[pos_ind];
	}
	inline const T& operator[](const unsigned& pos_ind) const {
		return e[pos_ind];
	}
	// output a vector, if << called, it will require operator << defined for T
	friend ostream& operator << (ostream& out, const Vec_s<T, cap>& v) {
		out << " (size: " << v.dim << ")" << endl;
		out << setprecision(4);

		for (int j = 0; j < v.dim; j++) {
			out << setw(14) << v[j];
		}

		out << setprecision(6);
		out << endl;
		return out;
	}

	/**
	 * .sort eleemnts from small to big, < as default, quick sort is the best choice
	 *
	 */
	void sort(char symbol_lt_gt = '<') {
		if (symbol_lt_gt == '<') {
			quick_sort_recur_lt(0, this->size() - 1);
		}
		else {
			quick_sort_recur_gt(0, this->size() - 1);
		}
	}
	void quick_sort_recur_lt(int start_ind, int end_ind) {
		if (start_ind >= end_ind)	return;
		int center = my::rand_int(start_ind, end_ind);
		this->switch_ele(start_ind, center);
		center = start_ind;
		int i = start_ind + 1;
		int j = end_ind;
		while (i != j) {
			if (this->e[i] <= this->e[center]) i++;
			else {
				this->switch_ele(i, j);
				j--;
			}
		}
		if (this->e[i] <= this->e[center]) {
			this->switch_ele(i, center);
			center = i;
		}
		else {
			this->switch_ele(i - 1, center);
			center = i - 1;
		}
		quick_sort_recur_lt(start_ind, center - 1);
		quick_sort_recur_lt(center + 1, end_ind);
	}
	void quick_sort_recur_gt(int start_ind, int end_ind) {
		if (start_ind >= end_ind)	return;
		int center = my::rand_int(start_ind, end_ind);
		this->switch_ele(start_ind, center);
		center = start_ind;
		int i = start_ind + 1;
		int j = end_ind;
		while (i != j) {
			if (this->e[i] >= this->e[center]) i++;
			else {
				this->switch_ele(i, j);
				j--;
			}
		}
		if (this->e[i] >= this->e[center]) {
			this->switch_ele(i, center);
			center = i;
		}
		else {
			this->switch_ele(i - 1, center);
			center = i - 1;
		}
		quick_sort_recur_gt(start_ind, center - 1);
		quick_sort_recur_gt(center + 1, end_ind);
	}

	/**
	 * .sort eleemnts from small to big, < as default, quick sort is the best choice
	 *
	 * \param v_ind: original index for the corresponding elements, the matrix can be recovered with this->permute_back(v_ind);
	 *				 we can still pass Vec_s_cmp<int,cap> as v_ind on this function
	 */
	void sort_with_ind(Vec_s<int, cap>& v_ind, char symbol_lt_gt = '<') {
		// initialize v_ind
		v_ind.resize(this->dim);
		for (int i = 0; i < this->dim; ++i) {
			v_ind[i] = i;
		}

		if (symbol_lt_gt == '<') {
			quick_sort_recur_lt_with_ind(0, this->size() - 1, v_ind);
		}
		else {
			quick_sort_recur_gt_with_ind(0, this->size() - 1, v_ind);
		}
	}
	void quick_sort_recur_lt_with_ind(int start_ind, int end_ind, Vec_s<int, cap>& v_ind) {
		if (start_ind >= end_ind)	return;
		int center = my::rand_int(start_ind, end_ind);
		this->switch_ele(start_ind, center);
		v_ind.switch_ele(start_ind, center);
		center = start_ind;
		int i = start_ind + 1;
		int j = end_ind;
		while (i != j) {
			if (this->e[i] <= this->e[center]) i++;
			else {
				this->switch_ele(i, j);
				v_ind.switch_ele(i, j);
				j--;
			}
		}
		if (this->e[i] <= this->e[center]) {
			this->switch_ele(i, center);
			v_ind.switch_ele(i, center);
			center = i;
		}
		else {
			this->switch_ele(i - 1, center);
			v_ind.switch_ele(i - 1, center);
			center = i - 1;
		}
		quick_sort_recur_lt_with_ind(start_ind, center - 1, v_ind);
		quick_sort_recur_lt_with_ind(center + 1, end_ind, v_ind);
	}
	void quick_sort_recur_gt_with_ind(int start_ind, int end_ind, Vec_s<int, cap>& v_ind) {
		if (start_ind >= end_ind)	return;
		int center = my::rand_int(start_ind, end_ind);
		this->switch_ele(start_ind, center);
		v_ind.switch_ele(start_ind, center);
		center = start_ind;
		int i = start_ind + 1;
		int j = end_ind;
		while (i != j) {
			if (this->e[i] >= this->e[center]) i++;
			else {
				this->switch_ele(i, j);
				v_ind.switch_ele(i, j);
				j--;
			}
		}
		if (this->e[i] >= this->e[center]) {
			this->switch_ele(i, center);
			v_ind.switch_ele(i, center);
			center = i;
		}
		else {
			this->switch_ele(i - 1, center);
			v_ind.switch_ele(i - 1, center);
			center = i - 1;
		}
		quick_sort_recur_gt_with_ind(start_ind, center - 1, v_ind);
		quick_sort_recur_gt_with_ind(center + 1, end_ind, v_ind);
	}

	// return max element and its index
	int max_ele_ind() const {
		int sss = this->size();
		if (sss > 0) {
			int ans = 0;
			T max_can = this->e[0];
			T can;
			bool update;
			for (int i = 1; i < sss; ++i) {
				can = this->e[i];
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
		int sss = this->size();
		if (sss > 0) {
			int ans = 0;
			T min_can = this->e[0];
			T can;
			bool update;
			for (int i = 1; i < sss; ++i) {
				can = this->e[i];
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
};


// Heap with max at top, with static storage
template<class T, int cap>
class Heap_max_s :public Vec_s<T, cap>
{
public:

	// redefine the default constructor from father to son, this is a must
	Heap_max_s() = default;
	Heap_max_s(int _dim) : Vec_s<T, cap>(_dim) {}
	Heap_max_s(const vector<T>& v) : Vec_s<T, cap>(v) {}

	// copy constructor to set father as son
	Heap_max_s(const Vec_s<T, cap>& father) :Vec_s<T, cap>(father) {}

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
			if (this->e[father_ind] >= this->e[(father_ind << 1) + 1] && this->e[father_ind] >= this->e[(father_ind << 1) + 2]);
			else {
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
		if (val <= this->e[ind]) {		// special for max heap
			this->e[ind] = val;
			descend_to_child(ind);
		}
		else {
			this->e[ind] = val;
			ascend_to_father(ind);
		}
	}
};


// Heap with min at top, with static storage
template<class T, int cap>
class Heap_min_s :public Vec_s<T, cap>
{
public:

	// redefine the default constructor from father to son, this is a must
	Heap_min_s() = default;
	Heap_min_s(int _dim) : Vec_s<T, cap>(_dim) {}
	Heap_min_s(const vector<T>& v) : Vec_s<T, cap>(v) {}

	// copy constructor to set father as son
	Heap_min_s(const Vec_s<T, cap>& father) :Vec_s<T, cap>(father) {}		// great

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
	void push(const T& val) {
		this->push_back(val);
		ascend_to_father(this->size() - 1);
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
			if (this->e[father_ind] <= this->e[(father_ind << 1) + 1] && this->e[father_ind] <= this->e[(father_ind << 1) + 2]);
			else {
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
};


// a singly-linked List in continous memory, and you can only remove contents once the list is built
template<class T, int cap>
class sList_s : public Vec_s<T, cap>
{
protected:
	int next_ind[cap];
	int start_ind;
	int current_ind;		// actually single linked list is enough for GE_PSM data structure.

public:
	// redefine the opearator since it cannot be inherited

	// output a vector, this should after build of list
	friend ostream& operator << (ostream& out, const sList_s<T, cap>& sL) {
		out << " (size: " << sL.dim << ")" << endl;
		out << setprecision(4);

		for (int i = sL.start_ind; i != -1; i = sL.next_ind[i]) {
			out << setw(14) << sL[i];
		}

		out << setprecision(6);
		out << endl;
		return out;
	}

	// redefine the default constructor from father to son, this is a must
	sList_s() = default;
	sList_s(int _dim) : Vec_s<T, cap>(_dim), next_ind(), start_ind(), current_ind() {}
	sList_s(const vector<T>& v) : Vec_s<T, cap>(v), next_ind(), start_ind(), current_ind() {}

	// copy constructor to set father as son
	sList_s(const Vec_s<T, cap>& father) :Vec_s<T, cap>(father), next_ind(), start_ind(), current_ind() {}	// great

	bool build() {			// return if build is success

		// we only build non-empty vector into list
		if (this->dim > 0) {

			// build the list link relation

			for (int i = 1; i < this->dim; ++i) {
				next_ind[i - 1] = i;
			}
			next_ind[this->dim - 1] = -1;			// the tail of list has no next node
			start_ind = 0;
			current_ind = 0;						// the current pointer always initialized to the head node index

			return true;
		}
		else {
			return false;
		}
	}
	void move_start() {
		current_ind = start_ind;
	}
	bool move_next() {		// return if move is success
		if (current_ind != -1) {
			current_ind = next_ind[current_ind];
			return true;
		}
		else {
			return false;
		}
	}
	bool remove_start() {	// return if remove is success
		if (start_ind != -1) {
			start_ind = next_ind[start_ind];
			this->dim--;
			return true;
		}
		else {
			return false;
		}
	}
	bool remove_next() {	// return if remove is success
		if (current_ind != -1 && next_ind[current_ind] != -1) {
			// current index is not the end of list, and also its next, then we can remove
			next_ind[current_ind] = next_ind[next_ind[current_ind]];
			this->dim--;
			return true;
		}
		else {
			return false;
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

};

// a doubly-linked List in continous memory, and you can only remove contents once the list is built
template<class T, int cap>
class dList_s : public Vec_s<T, cap>
{
	// to do: how to represent an empty list?

protected:
	int next_ind[cap];
	int former_ind[cap];

	int start_ind;
	int current_ind;		// current index will equals to -1 
							// while at this note current index can not be move_next or move_former, 
							// consider add end node and start node to the current index
							// to do: use a node for list containing element T
	int end_ind;

	Vec_s<int,cap> free_ind;

public:

	// output a vector, this should after build of list
	friend ostream& operator << (ostream& out, const dList_s<T, cap>& dL) {
		out << " (size: " << dL.dim << ")" << endl;
		out << setprecision(4);

		// -1 is an invalid index indicates that the node not existed
		for (int i = dL.start_ind; i != -1; i = dL.next_ind[i]) {
			out << setw(14) << dL[i];
		}

		out << setprecision(6);
		out << endl;
		return out;
	}

	// redefine the default constructor from father to son, this is a must
	dList_s() = default;
	dList_s(int _dim) : Vec_s<T, cap>(_dim), next_ind(), former_ind(), start_ind(), current_ind(), end_ind() {}
	dList_s(const vector<T>& v) : Vec_s<T, cap>(v), next_ind(), former_ind(), start_ind(), current_ind(), end_ind() {}

	// copy constructor to set father as son
	dList_s(const Vec_s<T, cap>& father) :Vec_s<T, cap>(father), next_ind(), former_ind(), start_ind(), current_ind(), end_ind() {}	// great

	// challange: add insert, all in O(1) with enough space decleared, holes allowed

	void build() {	// return if build is success

		// we only build non-empty vector into list, also we can build empty list to be inserted
		if (this->dim > 0) {

			// build the list link relation
			former_ind[0] = -1;						// the head of list has no former node

			for (int i = 1; i < this->dim; ++i) {
				next_ind[i - 1] = i;
				former_ind[i] = i - 1;
			}

			next_ind[this->dim - 1] = -1;			// the tail of list has no next node

			start_ind = 0;
			current_ind = 0;						// the current pointer always initialized to the head node index
			end_ind = this->dim - 1;

			// fill the free_ind, with index not occupied with elements
			for (int i = this->dim; i < cap; ++i) {
				free_ind.push_back(i);
			}
		}
	}
	void move_start() {
		current_ind = start_ind;
	}
	bool move_next() {		// return if move is successful
		if (current_ind != -1) {
			current_ind = next_ind[current_ind];
			return true;
		}
		else {
			return false;
		}
	}
	void move_former() {	// return if move is successful
		if (current_ind != -1) {
			current_ind = former_ind[current_ind];
			return true;
		}
		else {
			return false;
		}
	}
	void move_end() {
		current_ind = end_ind;
	}
	bool remove() {			// return if remove is successful
		if (current_ind != -1 && this->dim > 0) {
			// current index is not the end of list, then we can remove

			if (current_ind != start_ind) {
				// current ind is not the start of dList, link the former node to the next node
				next_ind[former_ind[current_ind]] = next_ind[current_ind];
			}
			else {
				// the start node is removed, update start_ind
				start_ind = next_ind[current_ind];
			}

			if (current_ind != end_ind) {
				// current ind is not the end of dList, link the next node to the former node
				former_ind[next_ind[current_ind]] = former_ind[current_ind];
			}
			else {
				// the end node is removed, update end_ind
				end_ind = former_ind[current_ind];
			}

			// add the index removed to free_ind
			free_ind.push_back(current_ind);

			// move current_ind to next node's index, this is a default setting to be remembered in dList
			current_ind = next_ind[current_ind];

			// decrease the size of dList
			this->dim--;

			return true;
		}
		else {
			return false;
		}
	}
	bool insert_next(const T& val) {
		if (!free_ind.empty()) {
			// we can only insert if there is free space

			if (!Vec_s<T, cap>::empty()) {
				// insert to an non-empty list

				if (current_ind != -1) {
					// current index is valid

					int allocated_ind = free_ind.back();
					free_ind.pop_back();

					// insert the val into the next index of current index in the dList
					this->e[allocated_ind] = val;
					next_ind[allocated_ind] = next_ind[current_ind];
					former_ind[allocated_ind] = current_ind;

					if (current_ind != end_ind) {
						// current_ind is not the end of dList, we have to form the backward linkage
						former_ind[next_ind[current_ind]] = allocated_ind;
					}
					else {
						// insert val to the end of dList, we only need to update end_ind
						end_ind = allocated_ind;
					}
					// form the forward linkage
					next_ind[current_ind] = allocated_ind;

					this->dim++;
					return true;
				}
				else {
					// current index is invalid
					return false;
				}
			}
			else {
				// insert to an empty list

				int allocated_ind = free_ind.back();
				free_ind.pop_back();

				// we don't have node in the list, current_ind must be -1

				this->e[allocated_ind] = val;
				next_ind[allocated_ind] = -1;
				former_ind[allocated_ind] = -1;

				start_ind = allocated_ind;
				current_ind = allocated_ind;			// the current pointer always initialized to the head node index
				end_ind = allocated_ind;

				this->dim++;
				return true;
			}
		}
		else {
			// no free space to insert
			return false;
		}
	}
	bool insert_former(const T& val) {
		if (!free_ind.empty()) {
			// we can only insert if there is free space

			if (!Vec_s<T, cap>::empty()) {
				// insert to an non-empty list

				if (current_ind != -1) {
					// current index is valid

					int allocated_ind = free_ind.back();
					free_ind.pop_back();

					// insert the val into the next index of current index in the dList
					this->e[allocated_ind] = val;
					former_ind[allocated_ind] = former_ind[current_ind];
					next_ind[allocated_ind] = current_ind;

					if (current_ind != start_ind) {
						// current_ind is not the start of dList, we have to form the forward linkage
						next_ind[former_ind[current_ind]] = allocated_ind;
					}
					else {
						// insert val to the start of dList, we only need to update start_ind
						start_ind = allocated_ind;
					}
					// form the backward linkage
					former_ind[current_ind] = allocated_ind;

					this->dim++;
					return true;
				}
				else {
					// current index is invalid
					return false;
				}
			}
			else {
				// insert to an empty list

				int allocated_ind = free_ind.back();
				free_ind.pop_back();

				// we don't have node in the list, current_ind must be -1

				this->e[allocated_ind] = val;
				next_ind[allocated_ind] = -1;
				former_ind[allocated_ind] = -1;

				start_ind = allocated_ind;
				current_ind = allocated_ind;			// the current pointer always initialized to the head node index
				end_ind = allocated_ind;

				this->dim++;
				return true;
			}
		}
		else {
			// no free space to insert
			return false;
		}
	}
	void clear() {
		// this will induce an empty list, we should define the state of an empty list
		move_start();
		while (current_ind != -1) {
			free_ind.push_back(current_ind);		// set all current index into the free index
			move_next();
		}
		start_ind = -1;								// indicate that no a single node in the list
		current_ind = -1;
		end_ind = -1;

		Vec_s<T, cap>::clear();
	}

	inline int get_start_ind() {
		return start_ind;
	}
	inline int get_current_ind() {
		return current_ind;
	}
	inline int get_end_ind() {
		return end_ind;
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
	inline T& end_val() {
		return this->e[end_ind];
	}
	inline const T& end_val() const {
		return this->e[end_ind];
	}
};

