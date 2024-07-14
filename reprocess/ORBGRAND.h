/*****************************************************************//**
 * \file   ORBGRAND.h
 * \brief  class of ordered reliablility bit GRAND, reliability values fixed to [1,2, .. n]
 * 
 * \author 26259
 * \date   April 2024
 *********************************************************************/

#pragma once


#include"OSD.h"

class Test_vector_ORBGRAND {
protected:
	//int max_weight;						// the maximum allowed logistic weights (LW)

	Matrix<int> candidate;					// logistic weights (LW) unused
	Matrix<int> cur;						// current generating test vector
	int ans_ind_cur;						// current generating index of 'ans'

	/**
	 * .generate test vectors recursively, wiht their logistic weights equals to target, and add them back to 'ans'
	 *
	 * \param start:	starting position of 'candidate'
	 * \param target:	required LW
	 * \param cur:		candidate test vector
	 */
	void gen_recur(int start, int target, Matrix<int>& cur) {

		if (target == 0) {
			//ans[ans_ind_cur] = cur;

			// only copy the valid part of 'cur'
			ans[ans_ind_cur].cp_valid(cur);				// saving 2M (19 M ->17 M) for decoding (63,45)BCH code, case of using 4e5 test vectors

			ans_ind_cur++;
			return;
		}

		if (target < 0 || start >= n || target < candidate[start] || ans_ind_cur == max_num) {
			return;
		}

		// flipping position of start
		cur.push_back(candidate[start] - 1);
		gen_recur(start + 1, target - candidate[start], cur);
		cur.pop_back();					// un-flip position of start

		// not flipping position of start 
		gen_recur(start + 1, target, cur);
	}

public:
	int max_num;							// the maximum allowed number of test vectors
	int n;									// length of code
	Matrix<Matrix<int>> ans;				// test vector results

	void gen(int _max_num, int _n) {

		n = _n; 
		max_num = _max_num;

		// initialize candidate
		candidate.resize(1, n);
		for (int i = 0; i < n; ++i) {
			candidate(i) = i + 1;
		}
		//cout << "candidate" << candidate;

		ans.resize(1, max_num);

		ans_ind_cur = 0;
		for (int i = 1; i <= n; ++i) {				// traverse all LW from small to large
			gen_recur(0, i, cur);

			if (ans_ind_cur == max_num - 1) {
				return;
			}
		}
	}
};

class ORBGRAND {
protected:

	int n;
	int k;
	int n_minus_k;
	Matrix<GF2> H;
	Test_vector_ORBGRAND TV_generator;
	int max_check_num;

	/* variables during decoding */

	Matrix<GF2> y;						// Hard-decisions of received vector
	Matrix<my_double> r_abs;			// Received vector in absolute value
	Matrix<int> permutation_first;		// Permutation that ensures decreasing reliability
	Matrix<GF2> y_bar;					// 'y_Gs' without 'permutation_second'
	Matrix<my_double> r_abs_bar;		// 'r_abs' sorted in REAL decreasing reliability, 'r_abs_Gs' without 'permutation_second'

	Matrix<GF2> syndrome;				// result of y * H^T
	Matrix<int> TEP_unpermuted;			// the permute-back version of 'TEP_generator.now' 

public:

	int type;
	Matrix<GF2> c_hat;

	int TEP_num;
	bool is_early_termination;

	ORBGRAND(const Matrix<GF2>& _H, int _max_check_num){
		n = _H.col();
		k = _H.col() - _H.row();
		n_minus_k = _H.row();
		max_check_num = _max_check_num;

		TV_generator.gen(max_check_num - 1, n);

		H = _H;
		syndrome.resize(1, n_minus_k);
		TEP_unpermuted.resize(1, n);
		r_abs.resize(1, n);
		y.resize(1, n);

		type = 2;						// the most simple type
		c_hat.resize(1, n);

		printf("ORBGRAND(%.0E)\n", (double)max_check_num);

		//cout << "TV_generator.ans" << TV_generator.ans;
	}

	void solve(const Matrix<my_double>& r) {

		for (int i = 0; i < n; ++i) {
			r_abs[i] = my::abs(r[i]);
			y[i] = r[i] > 0 ? 0 : 1;
		}
		permutation_first = Matrix<int>(1, n, 'N');								// Vector_ext::natual<n>(permutation_first);
		r_abs_bar = r_abs;
		r_abs_bar.quick_sort_recur_lt_with_ind(0, n - 1, permutation_first);	// r_abs.sort_gt_with_record(permutation_first, r_abs_bar);
		y_bar = y;
		y_bar.permute(permutation_first);										// y.permute(permutation_first, y_bar);

		// initialize the output variables
		TEP_num = 1;
		is_early_termination = false;

		//cout << "y" << y;

		compute_syndrome();
		if (syndrome.isZero() == true) {
			c_hat = y;

			is_early_termination = true;
			return;
		}


		while (TEP_num < max_check_num) {

			int TV_generator_ind = TEP_num - 1;

			if (is_permuted_codeword(TV_generator.ans[TV_generator_ind]) == true) {
				c_hat = y_bar;

				// adding the TEP
				int flip_num = TV_generator.ans[TV_generator_ind].size();
				for (int i = 0; i < flip_num; ++i) {
					c_hat(TV_generator.ans[TV_generator_ind][i]) += 1;
				}

				// permute back
				c_hat.permute_back(permutation_first);

				is_early_termination = true;
				return;
			}
			TEP_num++;
		}

		// will reach here
		c_hat.reset(0);
		return;
	}

	void compute_syndrome() {
		syndrome.reset(0);

		// using y
		for (int j = 0; j < n; ++j) {
			if (y(j) == 1) {
				// adding column-i to the syndrome
				for (int i = 0; i < n_minus_k; ++i) {
					syndrome(i) += H(i, j);
				}
			}
		}
	}

	bool is_permuted_codeword(const Matrix<int>& TEP_now) {
		int ts = TEP_now.size();
		TEP_unpermuted.resize(1, ts);
		for (int i = 0; i < ts; ++i) {
			TEP_unpermuted(i) = permutation_first(TEP_now(i));
		}


		for (int i = 0; i < n_minus_k; ++i) {
			GF2 syn_tmp = syndrome(i);
			for (int j = 0; j < ts; ++j) {
				syn_tmp += H(i, TEP_unpermuted(j));
			}
			if (syn_tmp != 0) {
				return false;
			}
		}
		return true;
	}
};
