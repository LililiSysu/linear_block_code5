/*****************************************************************//**
 * \file   SGRAND.h
 * \brief  Soft-GRAND as an ML decoder
 * 
 * \author 26259
 * \date   April 2024
 *********************************************************************/

#pragma once

#include"OSD.h"

struct reliability_and_TV_ind {
	my_double reliability;
	int TV_ind;
	bool operator < (const reliability_and_TV_ind& A) const {
		return reliability < A.reliability;		// creating a heap
	}
	bool operator > (const reliability_and_TV_ind& A) const {
		return reliability > A.reliability;		// creating a heap
	}
	reliability_and_TV_ind(my_double _reliability = 0, int _TV_ind = 0) {
		reliability = _reliability;
		TV_ind = _TV_ind;
	}
};

/**
 * . the difficulty of this algorithm is the tree search
 */

class SGRAND {
protected:

	int n;
	int k;
	int n_minus_k;
	Matrix<GF2> H;
	int max_check_num;


	/* variables during decoding */

	Matrix<GF2> y;						// Hard-decisions of received vector
	Matrix<my_double> r_abs;			// Received vector in absolute value
	Matrix<int> permutation_first;		// Permutation that ensures decreasing reliability
	Matrix<GF2> y_bar;					// 'y_Gs' without 'permutation_second'
	Matrix<my_double> r_abs_bar;		// 'r_abs' sorted in REAL decreasing reliability, 'r_abs_Gs' without 'permutation_second'

	Matrix<GF2> syndrome;				// result of y * H^T
	Matrix<int> TEP_unpermuted;			// the permute-back version of 'TEP_generator.now' 

	/* test vectors and selection structure */
	Matrix<Matrix<int>> TVs;
	Heap_min<reliability_and_TV_ind> hp;

public:

	int type;
	Matrix<GF2> c_hat;

	int TEP_num;
	bool is_early_termination;

	SGRAND(const Matrix<GF2>& _H, int _max_check_num) {
		n = _H.col();
		k = _H.col() - _H.row();
		n_minus_k = _H.row();
		max_check_num = _max_check_num;


		H = _H;
		syndrome.resize(1, n_minus_k);
		TEP_unpermuted.resize(1, n);
		r_abs.resize(1, n);
		y.resize(1, n);

		type = 3;						// the most simple type
		c_hat.resize(1, n);

		printf("SGRAND(%.0E)\n", (double)max_check_num);

		//cout << "TV_generator.ans" << TV_generator.ans;

		// declear enough space for TV in SGRAND
		TVs.resize(1, max_check_num);
		TVs.reset(Matrix<int>(1, n));	// it should be resize to 0 every new vector used

		hp.resize(1, max_check_num);	// it should be resize to 0 in every decoding event
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

		// initialize TVs and pq
		TVs[0].resize(1, 1);
		TVs[0][0] = 0;
		int TVs_range = 0;

		hp.resize(1,0);
		hp.push(reliability_and_TV_ind(r_abs_bar(0), 0));

		while (TEP_num < max_check_num) {

			int cur_TV_ind = hp.top().TV_ind;
			my_double cur_TV_reliability = hp.top().reliability;

			if (is_permuted_codeword(TVs[cur_TV_ind]) == true) {
				c_hat = y_bar;

				// adding the TEP
				int flip_num = TVs[cur_TV_ind].size();
				for (int i = 0; i < flip_num; ++i) {
					c_hat(TVs[cur_TV_ind][i]) += 1;
				}

				// permute back
				c_hat.permute_back(permutation_first);

				is_early_termination = true;
				return;
			}

			// update TVs and pq
			int cur_TV_last = TVs[cur_TV_ind].back();
			if (cur_TV_last != n - 1) {
				// we can generate a new TV
				my_double rel_add_pos = cur_TV_reliability + r_abs_bar(cur_TV_last + 1);
				my_double rel_move_pos = rel_add_pos - r_abs_bar(cur_TV_last);

				TVs_range++;
				TVs[TVs_range].cp_valid(TVs[cur_TV_ind]);	
				TVs[TVs_range].push_back(cur_TV_last + 1);	// genrate TV add_pos by appending to TVs, use current TV

				TVs[cur_TV_ind].back()++;			// generate the new TV move_pos in place

				// update hp to include new TVs' reliability and index
				hp.replace_top(reliability_and_TV_ind(rel_move_pos, cur_TV_ind));		// inplace operation of including TV move_pos
				hp.push(reliability_and_TV_ind(rel_add_pos, TVs_range));
			}
			else {
				hp.pop();		// by excluding index, we will not reach excluded TV
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
