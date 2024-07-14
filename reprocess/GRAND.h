/*****************************************************************//**
 * \file   GRAND.h
 * \brief  the original GRAND algorithm
 * 
 * \author 26259
 * \date   April 2024
 *********************************************************************/

#pragma once

#include"OSD.h"

class GRAND{
protected:

	int n;
	int k;
	int n_minus_k;
	Matrix<GF2> H;
	OSD_TEP TEP_generator;
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

	GRAND(const Matrix<GF2>& _H, int _max_check_num): TEP_generator(_H.col()) {
		n = _H.col();
		k = _H.col() - _H.row();
		n_minus_k = _H.row();
		max_check_num = _max_check_num;

		H = _H;
		syndrome.resize(1, n_minus_k); 
		TEP_unpermuted.resize(1, n);
		r_abs.resize(1, n);
		y.resize(1, n);

		type = 1;						// the most simple type
		c_hat.resize(1, n);

		printf("GRAND(%.0E)\n", (double)max_check_num);
	}

	void solve(const Matrix<my_double>& r) {

		for (int i = 0; i < n; ++i) {
			r_abs[i] = my::abs(r[i]);
			y[i] = r[i] > 0 ? 0 : 1;
		}
		permutation_first = Matrix<int>(1, n, 'N');								// Vector_ext::natual<n>(permutation_first);
		r_abs_bar = r_abs;
		r_abs_bar.quick_sort_recur_gt_with_ind(0, n - 1, permutation_first);	// r_abs.sort_gt_with_record(permutation_first, r_abs_bar);
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


		int order = 1;
		while (true) {
			TEP_generator.set_weight(order);
			while (TEP_generator.now[0] > 0) {

				if (TEP_num == max_check_num) {
					c_hat.reset(0);
					return;
				}

				if (is_permuted_codeword(TEP_generator.now) == true) {
					c_hat = y_bar;

					// adding the TEP
					for (int i = 0; i < order; ++i) {
						c_hat(TEP_generator.now(i)) += 1;
					}

					// permute back
					c_hat.permute_back(permutation_first);

					is_early_termination = true;
					return;
				}
				TEP_num++;
				TEP_generator.next();
			}
			order++;
		}

		// will never reach here
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
