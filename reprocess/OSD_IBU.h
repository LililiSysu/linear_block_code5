#pragma once
/*****************************************************************//**
 * \file   OSD_IBU.h
 * \brief  OSD utilizing IBU
 *
 * \author 26259
 * \date   January 2024
 *********************************************************************/

#include"IBU.h"
#include"OSD.h"

 // to be optimized
class OSD_IBUc {
protected:
	int n;
	int k;
	int min_Hamming_distance;
	int s;

	// variables during decoding
	Matrix<GF2> G_sys_fixed;		// The systematic genrator matrix for the code
	Matrix<GF2> y;						// Hard-decisions of received vector
	Matrix<GF2> y_Gs;					// Permuted hard-decisions of received vector
	Matrix<GF2> c0_Gs;					// First permuted codeword candidate
	Matrix<GF2> cb_Gs;					// The best permuted codeowrd candidate
	Matrix<int> nat;						// [0,1, ... ,n-1] for constant
	Matrix<int> inv_permutation_Gs_record;	// inverse permutation of ibu.r_abs_Gs_record
	Matrix<int> permutation_Gs_to_bar;	// Permutation that turn  'cb_Gs' into sequence with decreasing reliability
	my_double w_min;							// Minimum sort Hamming weight that is recorded through re-encoding
	OSD_TEP TEP_generator;				// TEP 'TEP_generator.now' will be generated sequentially
	Matrix<GF2> ce_Gs;					// Permuted codeowrd candidate generated by TEP 'TEP_generator.now'
	my_double w_ce_Gs;						// Sort Hamming weight of 'ce_Gs'
	int error_num;								// Number of different bits between 'y_Gs' and 'ce_Gs' ( 'c0_Gs' when in phase-0 )

	my_double G_threshold;					// threshold for ML criterion
	my_double G_OSD_threshold;				// soft Hamming weight lower bound for OSD phase-i re-processing
	Matrix<my_double> G_threshold_component;	// store r_abs_bar values adding in that constitute G_threshold

	Matrix<GF2> y_bar;					// 'y_Gs' without 'permutation_second'
	Matrix<GF2> cb_bar;					// 'cb_Gs' without 'permutation_second'

public:
	IBUc ibuc;
	int order;
	Matrix<GF2> c_hat;					// Optimal codeword for output
	int type;

	// simulation statistics
	int TEP_num;						// number of TEP used for re-encoding
	bool is_early_termination;			// if ture, the decoding is terminated by stopping condition

	OSD_IBUc(const Matrix<GF2>& _G, int _min_Hamming_distance, \
		int _s, int _order, bool is_cycle_used = true): ibuc(_G.col(), _G.row(), _s), TEP_generator(_G.col(), _G.row()) {

		n = _G.col();
		k = _G.row();
		s = _s;
		min_Hamming_distance = _min_Hamming_distance;


		G_sys_fixed = _G;
		order = _order;
		y.resize(1, n);
		c0_Gs.resize(1, n);
		ce_Gs.resize(1, n);
		nat = Matrix<int>(1, n, 'N');
		//Vector_ext::natual(nat);

		// initialize to arbitary values
		w_min = 0;
		w_ce_Gs = 0;
		error_num = 0;
		G_threshold = 0;
		G_OSD_threshold = 0;
		ibuc.set_cycle_used(is_cycle_used);
		if (is_cycle_used == true) {
			cout << "OSD_IBUc(" << _order << ")" << endl;
			type = 3;
		}
		else {
			cout << "OSD_IBU(" << _order << ")" << endl;
			type = 2;
		}

		TEP_num = 0;
		is_early_termination = false;
	}

	// the decoding result is 'c_hat', which can be fetched from the class
	void solve(const Matrix<my_double>& r) {
		for (int i = 0; i < n; ++i) {
			y[i] = r[i] > 0 ? 0 : 1;
		}

		ibuc.solve(G_sys_fixed, r);
		//cout << "ibuc.Gs" << ibuc.Gs;
		inv_permutation_Gs_record = nat;
		inv_permutation_Gs_record.permute_back(ibuc.r_abs_Gs_record);
		permutation_Gs_to_bar = inv_permutation_Gs_record;
		permutation_Gs_to_bar.permute(ibuc.r_abs_bar_record);

		y_bar = y;
		y_bar.permute(ibuc.r_abs_bar_record);
		y_Gs = y;
		y_Gs.permute(ibuc.r_abs_Gs_record);

		//y_Gs.get_part<k>(0, k - 1, yB_Gs);
		//cout << "y_Gs" << y_Gs;
		//cout << "yB_Gs" << yB_Gs;
		
		TEP_num = 0;
		is_early_termination = false;

		//c0_Gs = yB_Gs * ibuc.Gs;
		first_re_encoding();
		//cout << "c0_Gs" << c0_Gs;
		w_min = soft_Hamming_weight(c0_Gs);
		//cout << "w_min = " << w_min << endl;

		cb_Gs = c0_Gs;
		//cout << "cb_Gs" << cb_Gs;

		if (check_ML_codeword_general() == true) {		// in order-0, this is the same as the OSD specific type
			c_hat = cb_Gs;
			c_hat.permute_back(ibuc.r_abs_Gs_record);
			is_early_termination = true;
			return;
		}

		//cout << "c0_Gs" << c0_Gs;
		//cout << "d_min = " << d_min << endl;

		G_OSD_threshold = 0;
		for (int i = 1; i <= order; ++i) {
			TEP_generator.set_weight(i);

			if (check_ML_codeword_OSD_new_phase(i) == true) {
				c_hat = cb_Gs;
				c_hat.permute_back(ibuc.r_abs_Gs_record);
				is_early_termination = true;
				return;
			}

			while (TEP_generator.now[0] >= 0) {

				//cout << "TEP_generator.now" << TEP_generator.now;

				re_encoding(TEP_generator.now);

				//cout << "ce_Gs" << ce_Gs;
				w_ce_Gs = soft_Hamming_weight(ce_Gs);
				//cout << "w_ce_Gs = " << w_ce_Gs << endl;

				if (w_ce_Gs < w_min) {
					w_min = w_ce_Gs;
					cb_Gs = ce_Gs;

					/*if (check_ML_codeword_general() == true) {
						c_hat = cb_Gs.permute_back(permutation_all);
						return c_hat;
					}*/

					if (check_ML_codeword_OSD_specific(i) == true) {
						c_hat = cb_Gs;
						c_hat.permute_back(ibuc.r_abs_Gs_record);
						is_early_termination = true;
						return;
					}
				}
				TEP_generator.next();
			}
		}

		//cout << "cb_Gs" << cb_Gs;

		c_hat = cb_Gs;
		c_hat.permute_back(ibuc.r_abs_Gs_record);
		return;
	}

	void first_re_encoding() {
		// re-encode 'yB_Gs' into 'c0_Gs', i.e., c0_Gs = yB_Gs * ibuc.Gs;

		// reset to all-0 vector
		c0_Gs.reset(0);

		for (int j = 0; j < k; ++j) {

			if (y_Gs[j] == 1) {		// consider to write it as 'yB_Gs[j] + 1 == 0' for counting one GF2 operation

				// copy systematic part
				c0_Gs[j] = y_Gs[j];

				// re-encode parity-check part
				int row_start_ind = j * n;
				for (int p = k; p < n; ++p) {
					c0_Gs[p] += ibuc.Gs[row_start_ind + p];		// add accumulately, this cause my second problem, mark down here
				}
			}

		}

		TEP_num++;
	}

	void re_encoding(const Matrix<int>& TEP_now) {
		// re-encode to ce_Gs

		// initialize ans
		for (int i = 0; i < n; ++i) {
			ce_Gs[i] = c0_Gs[i];
		}

		for (int j = 0; j < TEP_now.size(); ++j) {
			// decide which bit to flip

			//cout << "ibuc.permute_record_MRIP" << ibuc.permute_record_MRIP;
			int flip_ind = ibuc.permute_record_MRIP[TEP_now[j]];	// flip the least reliable bit first in the MRIPs

			// flip systematic part
			ce_Gs[flip_ind] += 1;

			// re-encode parity-check part
			int row_start_ind = flip_ind * n;
			for (int p = k; p < n; ++p) {
				ce_Gs[p] += ibuc.Gs[row_start_ind + p];		// add accumulately, this cause my second problem, mark down here
			}
		}

		TEP_num++;
	}

	my_double soft_Hamming_weight(const Matrix<GF2>& c_Gs) {
		my_double ans = 0;
		error_num = 0;
		for (int i = 0; i < n; ++i) {
			if (c_Gs[i] != y_Gs[i]) {
				ans += ibuc.r_abs_Gs[i];
				error_num++;
			}
		}
		return ans;
	}

	bool check_ML_codeword_general() {
		// Will always check cb_Gs		
		// 'G_threshold(c_Gs,d)' is only needed to compute

		// Only need 'D0_add_size' reliabilities adding up to form the 'G_threshold'
		int D0_add_size = min_Hamming_distance - error_num;

		// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		cb_bar = cb_Gs;
		cb_bar.permute(permutation_Gs_to_bar);

		G_threshold = 0;
		G_threshold_component.clear();

		for (int i = n - 1; i >= 0 && D0_add_size > 0; --i) {

			if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
				G_threshold += ibuc.r_abs_bar[i];
				D0_add_size--;

				if (w_min < G_threshold) return true;
				G_threshold_component.push_back(ibuc.r_abs_bar[i]);
			}
		}
		//cout << "G_threshold = " << G_threshold << endl;

		return false;		// ML criterion for general codeword

	}

	bool check_ML_codeword_OSD_specific(int TEP_weight) {
		// Will always check cb_Gs		
		// 'G_threshold(c_Gs,d)' is only needed to compute

		G_threshold = G_OSD_threshold;
		/*for (int i = 0; i < TEP_weight; ++i) {
			G_threshold += ibuc.r_abs_Gs_B_bar[k - 1 - i];
		}*/

		// Only need 'D0_add_size' reliabilities adding up to 'G_threshold' further form the entire 'G_threshold'
		int D0_add_size = min_Hamming_distance - error_num - TEP_weight;

		// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		cb_bar = cb_Gs;
		cb_bar.permute(permutation_Gs_to_bar);

		G_threshold_component.clear();

		for (int i = n - 1; i >= 0 && D0_add_size > 0; --i) {

			if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
				G_threshold += ibuc.r_abs_bar[i];
				D0_add_size--;

				if (w_min < G_threshold) return true;
				G_threshold_component.push_back(ibuc.r_abs_bar[i]);
			}
		}
		return false;		// ML criterion for general codeword
	}

	bool check_ML_codeword_OSD_new_phase(int phase_num) {

		int flip_ind = ibuc.permute_record_MRIP[k - phase_num];	// the phase_num-th least reliable bit first in the MRIPs

		// update 'G_threshold' for a new judgement
		if (!G_threshold_component.empty()) {

			G_threshold += ibuc.r_abs_Gs[flip_ind] - G_threshold_component.back();
			G_threshold_component.pop_back();
		}
		else {
			G_threshold += ibuc.r_abs_Gs[flip_ind];
		}

		// judge the current best codeword if it is ML
		if (w_min < G_threshold) {
			return true;
		}
		else {
			// update 'G_OSD_threshold' for the following codeword candidates
			G_OSD_threshold += ibuc.r_abs_Gs[flip_ind];
			return false;
		}
	}

	void call_IBUc_only(const Matrix<my_double>& r) {
		ibuc.solve(G_sys_fixed, r);
	}
};

class OSD_IBUc_dual {
protected:
	int n;
	int n_minus_k;
	int min_Hamming_distance;
	int s;

	// variables during decoding
	Matrix<GF2> G_sys_fixed;		// The systematic genrator matrix for the code
	Matrix<GF2> y;						// Hard-decisions of received vector
	Matrix<GF2> y_Gs;					// Permuted hard-decisions of received vector
	Matrix<GF2> c0_Gs;					// First permuted codeword candidate
	Matrix<GF2> cb_Gs;					// The best permuted codeowrd candidate
	Matrix<int> nat;						// [0,1, ... ,n-1] for constant
	Matrix<int> inv_permutation_Gs_record;	// inverse permutation of ibu.r_abs_Gs_record
	Matrix<int> permutation_Gs_to_bar;	// Permutation that turn  'cb_Gs' into sequence with decreasing reliability
	my_double w_min;							// Minimum sort Hamming weight that is recorded through re-encoding
	OSD_TEP TEP_generator;				// TEP 'TEP_generator.now' will be generated sequentially
	Matrix<GF2> ce_Gs;					// Permuted codeowrd candidate generated by TEP 'TEP_generator.now'
	my_double w_ce_Gs;						// Sort Hamming weight of 'ce_Gs'
	int error_num;								// Number of different bits between 'y_Gs' and 'ce_Gs' ( 'c0_Gs' when in phase-0 )

	my_double G_threshold;					// threshold for ML criterion
	my_double G_OSD_threshold;				// soft Hamming weight lower bound for OSD phase-i re-processing
	Matrix<my_double> G_threshold_component;	// store r_abs_bar values adding in that constitute G_threshold

	Matrix<GF2> y_bar;					// 'y_Gs' without 'permutation_second'
	Matrix<GF2> cb_bar;					// 'cb_Gs' without 'permutation_second'

public:
	IBUc ibuc;
	int order;
	Matrix<GF2> c_hat;					// Optimal codeword for output
	int type;

	// simulation statistics
	int TEP_num;						// number of TEP used for re-encoding
	bool is_early_termination;			// if ture, the decoding is terminated by stopping condition

	OSD_IBUc_dual(const Matrix<GF2>& _G, int _min_Hamming_distance, \
		int _s, int _order, bool is_cycle_used = true) : ibuc(_G.col(), _G.row(), _s), TEP_generator(_G.col(), _G.col() - _G.row()) {

		n = _G.col();
		n_minus_k = _G.row();
		s = _s;
		min_Hamming_distance = _min_Hamming_distance;

		G_sys_fixed = _G;
		order = _order;
		y.resize(1, n);
		c0_Gs.resize(1, n);
		ce_Gs.resize(1, n);
		nat = Matrix<int>(1, n, 'N');

		// initialize to arbitary values
		w_min = 0;
		w_ce_Gs = 0;
		error_num = 0;
		G_threshold = 0;
		G_OSD_threshold = 0;

		ibuc.set_dual();
		ibuc.set_cycle_used(is_cycle_used);
		

		if (is_cycle_used == true) {
			cout << "OSD_IBUc_dual(" << _order << ")" << endl;
			type = 6;
		}
		else {
			cout << "OSD_IBU_dual(" << _order << ")" << endl;
			type = 5;
		}

		TEP_num = 0;
		is_early_termination = false;
	}

	// the decoding result is 'c_hat', which can be fetched from the class
	void solve(const Matrix<my_double>& r) {
		for (int i = 0; i < n; ++i) {
			y[i] = r[i] > 0 ? 0 : 1;
		}

		ibuc.solve(G_sys_fixed, r);
		//cout << "ibuc.Gs" << ibuc.Gs;
		inv_permutation_Gs_record = nat;
		inv_permutation_Gs_record.permute_back(ibuc.r_abs_Gs_record);
		permutation_Gs_to_bar = inv_permutation_Gs_record;
		permutation_Gs_to_bar.permute(ibuc.r_abs_bar_record);

		y_bar = y;
		y_bar.permute(ibuc.r_abs_bar_record);
		y_Gs = y;
		y_Gs.permute(ibuc.r_abs_Gs_record);

		//y_Gs.get_part<k>(0, k - 1, yB_Gs);
		//cout << "y_Gs" << y_Gs;
		//cout << "yB_Gs" << yB_Gs;


		TEP_num = 0;
		is_early_termination = false;

		//c0_Gs = yB_Gs * ibuc.Gs;
		//cout << "ibuc.Gs" << ibuc.Gs;
		first_re_encoding();
		//cout << "c0_Gs" << endl;
		//c0_Gs.print();
		w_min = soft_Hamming_weight(c0_Gs);
		//cout << "w_min = " << w_min << endl;

		cb_Gs = c0_Gs;
		//cout << "cb_Gs" << cb_Gs;

		if (check_ML_codeword_general() == true) {		// in order-0, this is the same as the OSD specific type
			c_hat = cb_Gs;
			c_hat.permute_back(ibuc.r_abs_Gs_record);
			is_early_termination = true;
			return;
		}

		//cout << "c0_Gs" << c0_Gs;
		//cout << "d_min = " << d_min << endl;

		G_OSD_threshold = 0;
		for (int i = 1; i <= order; ++i) {
			TEP_generator.set_weight(i);

			if (check_ML_codeword_OSD_new_phase(i) == true) {
				c_hat = cb_Gs;
				c_hat.permute_back(ibuc.r_abs_Gs_record);
				is_early_termination = true;
				return;
			}

			while (TEP_generator.now[0] >= 0) {

				//cout << "TEP_generator.now" << TEP_generator.now;

				re_encoding(TEP_generator.now);

				//cout << "ce_Gs" << endl;
				//ce_Gs.print();
				w_ce_Gs = soft_Hamming_weight(ce_Gs);
				//cout << "w_ce_Gs = " << w_ce_Gs << endl;

				if (w_ce_Gs < w_min) {
					w_min = w_ce_Gs;
					cb_Gs = ce_Gs;

					/*if (check_ML_codeword_general() == true) {
						c_hat = cb_Gs.permute_back(permutation_all);
						return c_hat;
					}*/

					if (check_ML_codeword_OSD_specific(i) == true) {
						c_hat = cb_Gs;
						c_hat.permute_back(ibuc.r_abs_Gs_record);
						is_early_termination = true;
						return;
					}
				}
				TEP_generator.next();
			}
		}

		//cout << "cb_Gs" << cb_Gs;

		c_hat = cb_Gs;
		c_hat.permute_back(ibuc.r_abs_Gs_record);
		return;
	}

	void first_re_encoding() {
		// re-encode 'yB_Hs' into 'c0_Hs', i.e., c0_Hs = yB_Hs * Hs;

		// reset to all-0 vector
		c0_Gs.reset(0);

		for (int j = n_minus_k; j < n; ++j) {

			if (y_Gs[j] == 1) {		// consider to write it as 'yB_Hs[j] + 1 == 0' for counting one GF2 operation

				// copy systematic part
				c0_Gs[j] = y_Gs[j];

				// re-encode parity-check part
				int col_traverse_ind = j;
				for (int p = 0; p < n_minus_k; ++p) {
					c0_Gs[p] += ibuc.Gs[col_traverse_ind];		// add accumulately, this cause my second problem, mark down here
					col_traverse_ind += n;
				}
			}

		}

		TEP_num++;
	}

	void re_encoding(const Matrix<int>& TEP_now) {
		// re-encode to ce_Gs

		// initialize ans
		for (int i = 0; i < n; ++i) {
			ce_Gs[i] = c0_Gs[i];
		}

		for (int j = 0; j < TEP_now.size(); ++j) {
			// decide which bit to flip

			//cout << "ibuc.permute_record_LRP" << ibuc.permute_record_LRP;
			int flip_ind = ibuc.permute_record_LRP[n - n_minus_k - TEP_now[j] - 1];		// flip the least reliable bit first in the MRIPs
			//cout << "flip_ind = " << flip_ind << endl;

			// flip systematic part
			ce_Gs[flip_ind] += 1;

			// re-encode parity-check part
			int col_traverse_ind = flip_ind;
			for (int p = 0; p < n_minus_k; ++p) {
				ce_Gs[p] += ibuc.Gs[col_traverse_ind];		// add accumulately, this cause my second problem, mark down here
				col_traverse_ind += n;
			}
		}

		TEP_num++;
	}

	my_double soft_Hamming_weight(const Matrix<GF2>& c_Gs) {
		my_double ans = 0;
		error_num = 0;
		for (int i = 0; i < n; ++i) {
			if (c_Gs[i] != y_Gs[i]) {
				ans += ibuc.r_abs_Gs[i];
				error_num++;
			}
		}
		return -ans;
	}

	bool check_ML_codeword_general() {
		// Will always check cb_Gs		
		// 'G_threshold(c_Gs,d)' is only needed to compute

		// Only need 'D0_add_size' reliabilities adding up to form the 'G_threshold'
		int D0_add_size = min_Hamming_distance - error_num;

		// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		cb_bar = cb_Gs;
		cb_bar.permute(permutation_Gs_to_bar);

		G_threshold = 0;
		G_threshold_component.clear();

		for (int i = 0; i < n && D0_add_size > 0; ++i) {

			if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
				G_threshold += ibuc.r_abs_bar[i];
				D0_add_size--;

				if (w_min < -G_threshold) return true;
				G_threshold_component.push_back(ibuc.r_abs_bar[i]);
			}
		}
		//cout << "G_threshold = " << G_threshold << endl;

		return false;		// ML criterion for general codeword

	}

	bool check_ML_codeword_OSD_specific(int TEP_weight) {
		// Will always check cb_Gs		
		// 'G_threshold(c_Gs,d)' is only needed to compute

		G_threshold = G_OSD_threshold;
		/*for (int i = 0; i < TEP_weight; ++i) {
			G_threshold += ibuc.r_abs_Gs_B_bar[k - 1 - i];
		}*/

		// Only need 'D0_add_size' reliabilities adding up to 'G_threshold' further form the entire 'G_threshold'
		int D0_add_size = min_Hamming_distance - error_num - TEP_weight;

		// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		cb_bar = cb_Gs;
		cb_bar.permute(permutation_Gs_to_bar);

		G_threshold_component.clear();

		for (int i = 0; i < n && D0_add_size > 0; ++i) {

			if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
				G_threshold += ibuc.r_abs_bar[i];
				D0_add_size--;

				if (w_min < -G_threshold) return true;
				G_threshold_component.push_back(ibuc.r_abs_bar[i]);
			}
		}
		return false;		// ML criterion for general codeword
	}

	bool check_ML_codeword_OSD_new_phase(int phase_num) {

		int flip_ind = ibuc.permute_record_LRP[phase_num - 1];	// the phase_num-th least reliable bit first in the MRIPs

		// update 'G_threshold' for a new judgement
		if (!G_threshold_component.empty()) {

			G_threshold += ibuc.r_abs_Gs[flip_ind] - G_threshold_component.back();
			G_threshold_component.pop_back();
		}
		else {
			G_threshold += ibuc.r_abs_Gs[flip_ind];
		}

		// judge the current best codeword if it is ML
		if (w_min < -G_threshold) {
			return true;
		}
		else {
			// update 'G_OSD_threshold' for the following codeword candidates
			G_OSD_threshold += ibuc.r_abs_Gs[flip_ind];
			return false;
		}
	}

	void call_IBUc_only(const Matrix<my_double>& r) {
		ibuc.solve(G_sys_fixed, r);
	}
};
