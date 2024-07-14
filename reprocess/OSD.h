#pragma once
/*****************************************************************//**
 * \file   OSD.h
 * \brief  Order Statistic Decoding class
 * 
 * \author 26259
 * \date   January 2024
 *********************************************************************/

#include"../my_lib/my_double.h"
#include"../GF/GF2.h"
#include"../my_lib/Matrix.h"

const int max_OSD_order = 10;

class OSD_TEP {
public:
	Matrix<int> now;		// supporting max order of 10
	int k;
	int w;

	OSD_TEP(int _k) {
		k = _k;
		w = 1;
		now.resize(1, max_OSD_order);
	}


	void set_weight(int _w) {
		w = _w;

		// initialize the TEP to be the most unreliable pattern
		now.resize(1, w);
		for (int i = 0; i < w; ++i) {
			now[i] = (k - w + i);
		}
	}

	/**
	 * . take k = 4, w = 2 as an example.
	 *   This function generate the next TEP (the variable now) in such a sequence
	 *
				 2             3		(the first TEP should be set manually)
				 1             3
				 1             2
				 0             3
				 0             2
				 0             1
	 *
	 * this is for a good OSD TEP ordering
	 *
	 */
	void next() {

		// minus one, counting carry bit
		int q = w - 1;
		while (q > 0) {
			if (now[q] - 1 == now[q - 1]) {
				// minus one to the former bit and this bit to be changed
			}
			else {
				break;
			}
			--q;
		}

		now[q]--;

		q++;
		// we must have q >= 0
		while (q < w) {
			now[q] = k - w + q;
			++q;
		}
	}
};

class OSD_TEP_seg {
public:
	int seg_1_len;
	int seg_2_len;
	int seg_1_w;
	int seg_2_w;

	Matrix<int> now;
	int len_all;
	int w_all;

	OSD_TEP_seg(int _seg_1_len, int _seg_2_len) :seg_1_len(_seg_1_len), seg_2_len(_seg_2_len) {
		now.resize(1, 2 * max_OSD_order);

		len_all = seg_1_len + seg_2_len;

		seg_1_w = 0;
		seg_2_w = 0;
		w_all = 0;
	}

	void set_weight(int _seg_1_w, int _seg_2_w) {
		seg_1_w = _seg_1_w;
		seg_2_w = _seg_2_w;
		w_all = seg_1_w + seg_2_w;
		now.resize(1, seg_1_w + seg_2_w);

		for (int i = 0; i < seg_1_w; ++i) {
			now[i] = (seg_1_len - seg_1_w + i);
		}

		for (int i = 0; i < seg_2_w; ++i) {
			now[seg_1_w + i] = (len_all - seg_2_w + i);
		}
	}

	void seg_1_next() {

		// minus one, counting carry bit
		int q = seg_1_w - 1;
		while (q > 0) {
			if (now[q] - 1 == now[q - 1]) {
				// minus one to the former bit and this bit to be changed
			}
			else {
				break;
			}
			--q;
		}

		if (q < 0) {
			now[0] = -1;		// first order is 0, stop generating TEPs
			return;
		}

		now[q]--;

		q++;
		// we must have q >= 0
		while (q < seg_1_w) {
			now[q] = seg_1_len - seg_1_w + q;
			++q;
		}
	}

	/**
	 * . take _seg_1_len = 3, _seg_2_len = 2, _seg_1_w = 1, seg_2_w = 1,  as an example.
	 *   This function generate the next TEP (the variable now) in such a sequence
	 *
				 2             4			(the first TEP should be set manually)
				 2             3
				 1			   4
				 1			   3			
				 0			   4
				 0			   3
	 *
	 * this is for a good OSD TEP ordering
	 *
	 */
	void next() {

		// minus one, counting carry bit

		int q = w_all - 1;
		while (q > seg_1_w) {
			if (now[q] - 1 == now[q - 1]) {
				// minus one to the former bit and this bit to be changed
			}
			else {
				break;
			}
			--q;
		}


		if (now[q] - 1 < seg_1_len) {
			// reset seg_2 TEP
			for (int i = 0; i < seg_2_w; ++i) {
				now[seg_1_w + i] = (len_all - seg_2_w + i);
			}

			// update seg_1 TEP
			seg_1_next();
		}
		else {
			now[q]--;
			// continue to update seg_2 TEP
			q++;
			// we must have q >= 0
			while (q < w_all) {
				now[q] = len_all - seg_2_w + q;
				++q;
			}

		}
	}
};

class OSD_base {
public:
	virtual void call_GE_only(const Matrix<my_double>& r) {};
	virtual void solve(const Matrix<my_double>& r) {};

	Matrix<GF2> c_hat;					// Optimal codeword for output
	int type;

	// simulation statistics
	int TEP_num;						// number of TEP used for re-encoding
	int float_TEP_num;					// number of TEP used for correlation distance calculation
	bool is_early_termination;			// if ture, the decoding is terminated by stopping condition

	// interface
	my_double sigma;					// denote the received signal noise, keep for order-skipping condition
	int ET_type;						// 0 -> no stopping criterion
										// 1 -> general ML stopping criterion
										// 2 -> OSD specific ML stopping criterion
										// 3 -> OS stopping cirterion
										// 4 -> Hamming distance criterion, threshold = d / 2 = t
										// 5 -> PNC + PSC
										// 6 -> FSD + PSC
										// 7 -> OS + PSC

	int tau_PSC;						// PSC termination parameter, discarding re-encoded codeword
	my_double beta;						// FSD termination parameter, order-skipping
};

class OSD_v4: public OSD_base {
protected:
	int n;
	int k;
	int min_Hamming_distance;
	int error_correction_capacity;			// min_Hamming_distance / 2

	// variables during decoding
	Matrix<int> r_abs_Gs_record;			// Permutation for obtaining MRIPs
	Matrix<GF2> y;						// Hard-decisions of received vector
	Matrix<my_double> r_abs;				// Received vector in absolute value
	Matrix<my_double> r_Gs;			// Permuted received vector
	Matrix<my_double> r_abs_Gs;		// Permuted received vector in absolute value
	Matrix<GF2> Gs;					// Systematic generator matrix corresponding to MRIPs
	Matrix<GF2> Hs;					// Systematic generator matrix corresponding to MRIPs
	Matrix<GF2> y_Gs;					// Permuted hard-decisions of received vector
	Matrix<GF2> c0_Gs;					// First permuted codeword candidate
	Matrix<GF2> cb_Gs;					// The best permuted codeowrd candidate
	my_double w_min;							// Minimum sort Hamming weight that is recorded through re-encoding
	OSD_TEP TEP_generator;				// TEP 'TEP_generator.now' will be generated sequentially
	Matrix<GF2> ce_Gs;					// Permuted codeowrd candidate generated by TEP 'TEP_generator.now'
	my_double w_ce_Gs;						// Sort Hamming weight of 'ce_Gs'
	int error_num;								// Number of different bits between 'y_Gs' and 'ce_Gs' ( 'c0_Gs' when in phase-0 )

	my_double G_threshold;					// threshold for ML criterion
	my_double G_OSD_threshold;				// soft Hamming weight lower bound for OSD phase-i re-processing
	Matrix<my_double> G_threshold_component;	// store r_abs_bar values adding in that constitute G_threshold

	Matrix<int> permutation_first;		// Permutation that ensures decreasing reliability
	Matrix<int> permutation_second;		// Permutation that ensures the Gs systematic on left, disrupts decreasing reliability
	Matrix<int> permutation_all;			// Permutation that combines 'permutation_first' and 'permutation_second'
	Matrix<my_double> r_abs_bar;			// 'r_abs' sorted in REAL decreasing reliability, 'r_abs_Gs' without 'permutation_second'
	Matrix<GF2> y_bar;					// 'y_Gs' without 'permutation_second'
	Matrix<GF2> cb_bar;					// 'cb_Gs' without 'permutation_second'

public:
	Matrix<GF2> G;
	Matrix<GF2> H;
	int order;
	bool by_generator;



	OSD_v4(const Matrix<GF2>& _G, const Matrix<GF2>& _H, int _d, int _order \
	) : TEP_generator(_G.row() <= _H.row() ? _G.row() : _H.col() - _H.row()) {

		by_generator = _G.row() <= _H.row();
		n = by_generator ? _G.col() : _H.col();

		type = 1;
		cout << "OSD_v4(" << _order << ")" << endl;

		if (by_generator == true) {
			G = _G;
			k = G.row();
		}
		else {
			H = _H;
			k = n - H.row();
			Gs = Matrix<GF2>(k, n, 'i');
		}

		min_Hamming_distance = _d;
		error_correction_capacity = _d / 2;
		order = _order;
		r_abs.resize(1, n);
		y.resize(1, n);
		c0_Gs.resize(1, n);
		ce_Gs.resize(1, n);

		// initialize to arbitary values
		w_min = 0;
		w_ce_Gs = 0;
		error_num = 0;
		G_threshold = 0;
		G_OSD_threshold = 0;

		TEP_num = 0;
		float_TEP_num = 0;
		is_early_termination = false;

		sigma = 0;				// reset arbitary, need to change with the channel SNR
		ET_type = 3;			// the most efficient and non-significant performance loss ET

		tau_PSC = n + 1;		// by default, not skipping any re-encoded codeword
		beta = 0;				// by default, FSD = PNC
	}

	void call_GE_only(const Matrix<my_double>& r) {

		//for (int i = 0; i < n; ++i) {
		//	r_abs[i] = my::abs(r[i]);
		//	y[i] = r[i] > 0 ? 0 : 1;
		//}
		//permutation_first = Matrix<int>(1, n, 'N');								// Vector_ext::natual<n>(permutation_first);
		//r_abs_bar = r_abs;
		//r_abs_bar.quick_sort_recur_gt_with_ind(0, n - 1, permutation_first);	// r_abs.sort_gt_with_record(permutation_first, r_abs_bar);
		//y_bar = y;
		//y_bar.permute(permutation_first);										// y.permute(permutation_first, y_bar);
		//Gs = G;
		//Gs.permute_col(permutation_first);
		//// G.permute_column(permutation_first, Gs);
		//permutation_second = Matrix<int>(1, n, 'N'); // Vector_ext::natual<n>(permutation_second);
		//Gs.GJE_4_GF2_left(permutation_second);
		////cout << "Gs" << Gs;

		for (int i = 0; i < n; ++i) {
			r_abs[i] = my::abs(r[i]);
			y[i] = r[i] > 0 ? 0 : 1;
		}
		permutation_first = Matrix<int>(1, n, 'N');								// Vector_ext::natual<n>(permutation_first);
		r_abs_bar = r_abs;
		r_abs_bar.quick_sort_recur_gt_with_ind(0, n - 1, permutation_first);	// r_abs.sort_gt_with_record(permutation_first, r_abs_bar);
		y_bar = y;
		y_bar.permute(permutation_first);										// y.permute(permutation_first, y_bar);

		if (by_generator == true) {
			Gs = G;
			Gs.permute_col(permutation_first);
			// G.permute_column(permutation_first, Gs);

			permutation_second = Matrix<int>(1, n, 'N'); // Vector_ext::natual<n>(permutation_second);
			Gs.GJE_4_GF2_left(permutation_second);

			//Gs.GE_left_identity_match_jiabao(permutation_second);

			//cout << "Gs" << Gs;
		}
		else {
			Hs = H;
			Hs.permute_col(permutation_first);

			permutation_second = Matrix<int>(1, n, 'N'); // Vector_ext::natual<n>(permutation_second);
			Hs.GJE_4_GF2_right(permutation_second);

			//Hs.GE_right_identity_match_jiabao(permutation_second);

			// set Gs = [I_k | P^T]
			for (int i = 0; i < k; ++i) {
				for (int j = k; j < n; ++j) {
					Gs(i, j) = Hs(j - k, i);
				}
			}
			//cout << "Gs" << Gs;
		}

		permutation_all = permutation_first;
		permutation_all.permute(permutation_second);	// combine the two permutation, forming a permutation to obtain MRIPs
		//cout << "permutation_all" << permutation_all;

		r_Gs = r;
		r_Gs.permute(permutation_all);			// remember this is a new permutation, all vectors should be re-freshed
		r_abs_Gs = r_abs;
		r_abs_Gs.permute(permutation_all);	// absence of this line cause my first problem, mark down here
		y_Gs = y;
		y_Gs.permute(permutation_all);
		//cout << "r_Gs" << r_Gs;
	}

	// the decoding result is 'c_hat', which can be fetched from the class
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

		if (by_generator == true) {
			Gs = G;
			Gs.permute_col(permutation_first);
			// G.permute_column(permutation_first, Gs);

			permutation_second = Matrix<int>(1, n, 'N'); // Vector_ext::natual<n>(permutation_second);
			
			Gs.GJE_4_GF2_left(permutation_second);

			//Gs.GE_left_identity_match_jiabao(permutation_second);

			//cout << "Gs" << Gs;
		}
		else {
			Hs = H;
			Hs.permute_col(permutation_first);

			permutation_second = Matrix<int>(1, n, 'N'); // Vector_ext::natual<n>(permutation_second);
			Hs.GJE_4_GF2_right(permutation_second);
			
			//Hs.GJE_right_v2(permutation_second);
			
			//Hs.GE_right_identity_match_jiabao(permutation_second);


			// set Gs = [I_k | P^T]
			for (int i = 0; i < k; ++i) {
				for (int j = k; j < n; ++j) {
					Gs(i, j) = Hs(j - k, i);
				}
			}
			//cout << "Gs" << Gs;
		}

		permutation_all = permutation_first;
		permutation_all.permute(permutation_second);	// combine the two permutation, forming a permutation to obtain MRIPs
		//cout << "permutation_all" << permutation_all;

		r_Gs = r;
		r_Gs.permute(permutation_all);			// remember this is a new permutation, all vectors should be re-freshed
		r_abs_Gs = r_abs;
		r_abs_Gs.permute(permutation_all);	// absence of this line cause my first problem, mark down here
		y_Gs = y;
		y_Gs.permute(permutation_all);
		//cout << "r_Gs" << r_Gs;


		//cout << "my_double_auxiliary_storage::operation_number (1) = " << my_double_auxiliary_storage::operation_number << endl;

		//y_Gs.get_part<k>(0, k - 1, yB_Gs);
		//cout << "y_Gs" << y_Gs;
		//cout << "yB_Gs" << yB_Gs;

		TEP_num = 0;
		float_TEP_num = 0;
		is_early_termination = false;

		//c0_Gs = yB_Gs * Gs;
		first_re_encoding();
		//cout << "c0_Gs" << c0_Gs;
		w_min = soft_Hamming_weight(c0_Gs);
		//cout << "w_min = " << w_min << endl;

		cb_Gs = c0_Gs;
		//cout << "cb_Gs" << cb_Gs;

		//cout << "my_double_auxiliary_storage::operation_number (2) = " << my_double_auxiliary_storage::operation_number << endl;

		if (phase_0_termination() == true) {
			return;
		}

		//cout << "my_double_auxiliary_storage::operation_number (3) = " << my_double_auxiliary_storage::operation_number << endl;

		//cout << "c0_Gs" << c0_Gs;
		//cout << "d_min = " << d_min << endl;

		//cout << "G_threshold_component" << G_threshold_component;
		//cout << "G_threshold = " << G_threshold << endl;


		for (int i = 1; i <= order; ++i) {
			TEP_generator.set_weight(i);

			if (new_order_termination(i) == true) {
				return;
			}

			while (TEP_generator.now[0] >= 0) {

				//cout << "TEP_generator.now" << TEP_generator.now;

				re_encoding(TEP_generator.now);

				if (judge_discard_codeword() == true) {
					TEP_generator.next();
					continue;
				}

				//cout << "ce_Gs" << ce_Gs;
				w_ce_Gs = soft_Hamming_weight(ce_Gs);
				//cout << "w_ce_Gs = " << w_ce_Gs << endl;

				if (w_ce_Gs < w_min) {
					w_min = w_ce_Gs;
					cb_Gs = ce_Gs;

					if (optimal_termination(i) == true) {
						return;
					}
				}
				TEP_generator.next();
			}
		}

		//cout << "my_double_auxiliary_storage::operation_number (4) = " << my_double_auxiliary_storage::operation_number << endl;

		//cout << "cb_Gs" << cb_Gs;
		c_hat = cb_Gs;
		c_hat.permute_back(permutation_all);
		return;
	}

	void first_re_encoding() {
		// re-encode 'yB_Gs' into 'c0_Gs', i.e., c0_Gs = yB_Gs * Gs;

		// reset to all-0 vector
		c0_Gs.reset(0);

		for (int j = 0; j < k; ++j) {

			if (y_Gs[j] == 1) {		// consider to write it as 'yB_Gs[j] + 1 == 0' for counting one GF2 operation

				// copy systematic part
				c0_Gs[j] = y_Gs[j];

				// re-encode parity-check part

				GF2_auxiliary_storage::re_encoding_bit_plane_norm_number++;

				int row_start_ind = j * n;
				for (int p = k; p < n; ++p) {
					c0_Gs[p] += Gs[row_start_ind + p];		// add accumulately, this cause my second problem, mark down here
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
			// flip systematic part
			ce_Gs[TEP_now[j]] += 1;

			// re-encode parity-check part

			GF2_auxiliary_storage::re_encoding_bit_plane_norm_number++;

			int row_start_ind = TEP_now[j] * n;
			for (int p = k; p < n; ++p) {
				ce_Gs[p] += Gs[row_start_ind + p];		// add accumulately, this cause my second problem, mark down here
			}
		}

		TEP_num++;
	}

	my_double soft_Hamming_weight(const Matrix<GF2>& c_Gs) {
		my_double ans = 0;
		error_num = 0;
		for (int i = 0; i < n; ++i) {
			if (c_Gs[i] != y_Gs[i]) {
				ans += r_abs_Gs[i];
				error_num++;
			}
		}
		float_TEP_num++;
		return ans;
	}

	void compute_Hamming_weight(const Matrix<GF2>& c_Gs) {
		error_num = 0;
		for (int i = 0; i < n; ++i) {
			if (c_Gs[i] != y_Gs[i]) {
				error_num++;		// compute Hamming distance of current codeword and hard decision vector
			}
		}
	}

	bool check_ML_codeword_general() {
		// Will always check cb_Gs		
		// 'G_threshold(c_Gs,d)' is only needed to compute


		// Only need 'D0_add_size' reliabilities adding up to form the 'G_threshold'
		int D0_add_size = min_Hamming_distance - error_num;

		// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		cb_bar = cb_Gs;
		cb_bar.permute_back(permutation_second);

		G_threshold = 0;

		for (int i = n - 1; i >= 0 && D0_add_size > 0; --i) {

			if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
				G_threshold += r_abs_bar[i];
				D0_add_size--;

				//if (w_min < G_threshold) return true;		// in align with Guo's Program
			}
		}
		return (w_min < G_threshold);
		//return false;		// ML criterion for general codeword
	}

	bool check_ML_codeword_OSD_specific(int TEP_weight) {
		// Will always check cb_Gs		
		// 'G_threshold(c_Gs,d)' is only needed to compute

		G_threshold = G_OSD_threshold;
		/*for (int i = 0; i < TEP_weight; ++i) {
			G_threshold += r_abs_Gs[k - 1 - i];
		}*/

		// Only need 'D0_add_size' reliabilities adding up to 'G_threshold' further form the entire 'G_threshold'
		int D0_add_size = min_Hamming_distance - error_num - TEP_weight;

		// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		cb_bar = cb_Gs;
		cb_bar.permute_back(permutation_second);

		G_threshold_component.clear();

		for (int i = n - 1; i >= 0 && D0_add_size > 0; --i) {

			if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
				G_threshold += r_abs_bar[i];
				D0_add_size--;

				//if (w_min < G_threshold) return true;		// in align with Guo
				G_threshold_component.push_back(r_abs_bar[i]);
			}
		}
		return (w_min < G_threshold);
		//return false;		// ML criterion for general codeword
	}

	bool check_ML_codeword_OSD_new_phase(int phase_num) {

		//cout << "r_abs_Gs" << r_abs_Gs;
		//cout << "r_abs_Gs[k - phase_num] = " << r_abs_Gs[k - phase_num] << endl;

		// update 'G_threshold' for a new judgement
		if (!G_threshold_component.empty()) {

			G_threshold += r_abs_Gs[k - phase_num] - G_threshold_component.back();
			G_threshold_component.pop_back();
		}
		else {
			G_threshold += r_abs_Gs[k - phase_num];
		}

		// judge the current best codeword if it is ML
		if (w_min < G_threshold) {
			return true;
		}
		else {
			// update 'G_OSD_threshold' for the following codeword candidates
			G_OSD_threshold += r_abs_Gs[k - phase_num];
			return false;
		}
	}

	void G_threshold_OS_OSD_ini() {
		G_threshold = 0;

		my_double scale_by_sigma = 2 / (sigma * sigma);
		for (int i = k; i < n; ++i) {
			G_threshold += r_abs_Gs(i) / (1 + exp(scale_by_sigma * r_abs_Gs(i)));
		}
	}

	void G_threshold_OS_OSD_increase(int phase_num) {
		G_threshold += r_abs_Gs(k - phase_num);
	}

	void early_termination() {
		c_hat = cb_Gs;
		c_hat.permute_back(permutation_all);
		is_early_termination = true;
		return;
	}

	bool phase_0_termination() {
		if (ET_type == 1) {
			if (check_ML_codeword_general() == true) {		// in order-0, this is the same as the OSD specific type			
				early_termination();
				return true;
			}
		}
		else if (ET_type == 2) {
			G_OSD_threshold = 0;

			if (check_ML_codeword_OSD_specific(0) == true) {		// in order-0, this is the same as the OSD specific type			
				early_termination();
				return true;
			}
		}
		else if (ET_type == 3) {
			G_threshold_OS_OSD_ini();
			if (w_min < G_threshold) {		// use order-skipping ET in order-0				
				early_termination();
				return true;
			}
		}
		else if (ET_type == 4) {
			if (error_num <= error_correction_capacity) {
				early_termination();
				return true;
			}
		}
		else if (ET_type == 5) {
			G_threshold = 0;
			if (error_num <= error_correction_capacity) {
				early_termination();
				return true;
			}
		}
		else if (ET_type == 6) {
			G_threshold = beta;				// beta is an emperical parameter
			if (error_num <= error_correction_capacity) {
				early_termination();
				return true;
			}
		}
		else if (ET_type == 7) {
			if (error_num <= error_correction_capacity) {
				early_termination();
				return true;
			}
			G_threshold_OS_OSD_ini();
		}

		return false;
	}

	bool new_order_termination(int phase_num) {
		if (ET_type == 2) {
			if (check_ML_codeword_OSD_new_phase(phase_num) == true) {
				early_termination();
				//cout << "my_double_auxiliary_storage::operation_number (4) = " << my_double_auxiliary_storage::operation_number << endl;
				//cout << "GF2e_auxiliary_storage::operation_number (4) = " << GF2e_auxiliary_storage::operation_number << endl;
				return true;
			}
		}
		else if (ET_type == 3) {
			G_threshold_OS_OSD_increase(phase_num);
			if (w_min < G_threshold) {
				early_termination();
				//cout << "my_double_auxiliary_storage::operation_number (4) = " << my_double_auxiliary_storage::operation_number << endl;
				//cout << "GF2e_auxiliary_storage::operation_number (4) = " << GF2e_auxiliary_storage::operation_number << endl;
				return true;
			}
		}
		else if (ET_type == 5 || ET_type == 6 || ET_type == 7) {
			G_threshold_OS_OSD_increase(phase_num);
			if (w_min < G_threshold) {
				early_termination();
				return true;
			}
		}
		return false;
	}

	bool optimal_termination(int phase_num) {
		if (ET_type == 1) {
			if (check_ML_codeword_general() == true) {
				early_termination();
				return true;
			}
		}
		else if (ET_type == 2) {
			if (check_ML_codeword_OSD_specific(phase_num) == true) {		// it actually is an ML criterion even for LLOSD		
				early_termination();
				//cout << "my_double_auxiliary_storage::operation_number (4) = " << my_double_auxiliary_storage::operation_number << endl;
				//cout << "GF2e_auxiliary_storage::operation_number (4) = " << GF2e_auxiliary_storage::operation_number << endl;
				return true;
			}
		}
		else if (ET_type == 3) {
			if (w_min < G_threshold) {
				early_termination();
				//cout << "my_double_auxiliary_storage::operation_number (4) = " << my_double_auxiliary_storage::operation_number << endl;
				//cout << "GF2e_auxiliary_storage::operation_number (4) = " << GF2e_auxiliary_storage::operation_number << endl;
				return true;
			}
		}
		else if (ET_type == 4) {
			if (error_num <= error_correction_capacity) {
				early_termination();
				return true;
			}
		}
		else if (ET_type == 5 || ET_type == 6 || ET_type == 7) {
			if (error_num <= error_correction_capacity) {
				early_termination();
				return true;
			}
		}

		return false;
	}

	bool judge_discard_codeword() {
		if (ET_type == 5 || ET_type == 6 || ET_type == 7) {
			compute_Hamming_weight(ce_Gs);
			return error_num >= tau_PSC;
		}
		return false;
	}
};


// the below class can be omitted

class OSD_v4_dual {
protected:
	int n;
	int n_minus_k;
	int min_Hamming_distance;

	// variables during decoding
	Matrix<int> r_abs_Hs_record;			// Permutation for obtaining MRIPs
	Matrix<GF2> y;						// Hard-decisions of received vector
	Matrix<my_double> r_abs;				// Received vector in absolute value
	Matrix<my_double> r_Hs;			// Permuted received vector
	Matrix<my_double> r_abs_Hs;		// Permuted received vector in absolute value
	Matrix<GF2> Hs;					// Systematic generator matrix corresponding to MRIPs
	Matrix<GF2> y_Hs;					// Permuted hard-decisions of received vector
	Matrix<GF2> c0_Hs;					// First permuted codeword candidate
	Matrix<GF2> cb_Hs;					// The best permuted codeowrd candidate
	my_double w_min;							// Minimum sort Hamming weight that is recorded through re-encoding
	OSD_TEP TEP_generator;				// TEP 'TEP_generator.now' will be generated sequentially
	Matrix<GF2> ce_Hs;					// Permuted codeowrd candidate generated by TEP 'TEP_generator.now'
	my_double w_ce_Hs;						// Sort Hamming weight of 'ce_Hs'
	int error_num;								// Number of different bits between 'y_Hs' and 'ce_Hs' ( 'c0_Hs' when in phase-0 )

	my_double G_threshold;					// threshold for ML criterion
	my_double G_OSD_threshold;				// soft Hamming weight lower bound for OSD phase-i re-processing
	Matrix<my_double> G_threshold_component;	// store r_abs_bar values adding in that constitute G_threshold

	Matrix<int> permutation_first;		// Permutation that ensures decreasing reliability
	Matrix<int> permutation_second;		// Permutation that ensures the Hs systematic on left, disrupts decreasing reliability
	Matrix<int> permutation_all;			// Permutation that combines 'permutation_first' and 'permutation_second'
	Matrix<my_double> r_abs_bar;			// 'r_abs' sorted in REAL decreasing reliability, 'r_abs_Hs' without 'permutation_second'
	Matrix<GF2> y_bar;					// 'y_Hs' without 'permutation_second'
	Matrix<GF2> cb_bar;					// 'cb_Hs' without 'permutation_second'

public:
	Matrix<GF2> H;
	int order;
	Matrix<GF2> c_hat;					// Optimal codeword for output
	int type;

	// simulation statistics
	int TEP_num;						// number of TEP used for re-encoding
	bool is_early_termination;			// if ture, the decoding is terminated by stopping condition

	OSD_v4_dual(const Matrix<GF2>& _H, int _min_Hamming_distance, int _order) : TEP_generator(_H.col() - _H.row()) {
		cout << "OSD_dual(" << _order << ")" << endl;
		type = 4;
		n = _H.col();
		n_minus_k = _H.row();		// actually the number of redundant bits
		min_Hamming_distance = _min_Hamming_distance;

		H = _H;
		order = _order;
		r_abs.resize(1, n);
		y.resize(1, n);
		c0_Hs.resize(1, n);
		ce_Hs.resize(1, n);

		// initialize to arbitary values
		w_min = 0;
		w_ce_Hs = 0;
		error_num = 0;
		G_threshold = 0;
		G_OSD_threshold = 0;

		TEP_num = 0;
		is_early_termination = false;
	}

	// the decoding result is 'c_hat', which can be fetched from the class
	void solve(const Matrix<my_double>& r) {
		for (int i = 0; i < n; ++i) {
			r_abs[i] = -my::abs(r[i]);		// changing sign, therefore, sort with gt results in LRPs first
			y[i] = r[i] > 0 ? 0 : 1;
		}
		permutation_first = Matrix<int>(1, n, 'N');
		r_abs_bar = r_abs;
		r_abs_bar.quick_sort_recur_gt_with_ind(0, n - 1, permutation_first);
		//r_abs.sort_gt_with_record(permutation_first, r_abs_bar);
		y_bar = y;
		y_bar.permute(permutation_first);

		Hs = H;
		Hs.permute_col(permutation_first);

		permutation_second = Matrix<int>(1, n, 'N');
		//cout << "permutation_second = " << permutation_second;
		//Vector_ext::natual<n>(permutation_second);
		Hs.GJE_4_GF2_left(permutation_second);
		//cout << "Hs" << Hs;

		permutation_all = permutation_first;
		permutation_all.permute(permutation_second);	// combine the two permutation, forming a permutation to obtain MRIPs
		//cout << "permutation_all" << permutation_all;

		r_Hs = r;
		r_Hs.permute(permutation_all);			// remember this is a new permutation, all vectors should be re-freshed
		r_abs_Hs = r_abs;
		r_abs_Hs.permute(permutation_all);	// absence of this line cause my first problem, mark down here
		y_Hs = y;
		y_Hs.permute(permutation_all);
		//cout << "r_Hs" << r_Hs;

		//y_Hs.get_part<k>(0, k - 1, yB_Hs);
		//cout << "y_Hs" << y_Hs;
		//cout << "yB_Hs" << yB_Hs;

		TEP_num = 0;
		is_early_termination = false;

		//c0_Hs = yB_Hs * Hs;
		first_re_encoding();
		//cout << "c0_Hs" << c0_Hs;
		w_min = soft_Hamming_weight(c0_Hs);
		//cout << "w_min = " << w_min << endl;

		cb_Hs = c0_Hs;
		//cout << "cb_Hs" << cb_Hs;

		if (check_ML_codeword_general() == true) {		// in order-0, this is the same as the OSD specific type

			//cout << "G_threshold_component" << G_threshold_component;

			c_hat = cb_Hs;
			c_hat.permute_back(permutation_all);
			is_early_termination = true;
			return;
		}
		//cout << "G_threshold_component" << G_threshold_component;
		//cout << "G_threshold = " << G_threshold << endl;

		//cout << "c0_Hs" << c0_Hs;
		//cout << "d_min = " << d_min << endl;


		G_OSD_threshold = 0;
		for (int i = 1; i <= order; ++i) {
			TEP_generator.set_weight(i);
			if (check_ML_codeword_OSD_new_phase(i) == true) {
				c_hat = cb_Hs;
				c_hat.permute_back(permutation_all);
				is_early_termination = true;
				return;
			}

			//cout << "G_OSD_threshold = " << G_OSD_threshold << endl;

		/*for (int i = 1; i <= order; ++i) {
			TEP_generator.set_weight(i);*/
			while (TEP_generator.now[0] >= 0) {

				//cout << "TEP_generator.now" << TEP_generator.now;

				re_encoding(TEP_generator.now);

				//cout << "ce_Hs" << ce_Hs;
				w_ce_Hs = soft_Hamming_weight(ce_Hs);
				//cout << "w_ce_Hs = " << w_ce_Hs << endl;

				if (w_ce_Hs < w_min) {
					w_min = w_ce_Hs;
					cb_Hs = ce_Hs;

					/*if (check_ML_codeword_general() == true) {
						c_hat = cb_Hs.permute_back(permutation_all);
						return c_hat;
					}*/

					if (check_ML_codeword_OSD_specific(i) == true) {
						c_hat = cb_Hs;
						c_hat.permute_back(permutation_all);
						is_early_termination = true;
						return;
					}
				}
				TEP_generator.next();
			}
		}

		//cout << "cb_Hs" << cb_Hs;

		c_hat = cb_Hs;
		c_hat.permute_back(permutation_all);
		return;
	}

	void first_re_encoding() {
		// re-encode 'yB_Hs' into 'c0_Hs', i.e., c0_Hs = yB_Hs * Hs;

		// reset to all-0 vector
		c0_Hs.reset(0);

		for (int j = n_minus_k; j < n; ++j) {

			if (y_Hs[j] == 1) {		// consider to write it as 'yB_Hs[j] + 1 == 0' for counting one GF2 operation

				// copy systematic part
				c0_Hs[j] = y_Hs[j];

				// re-encode parity-check part
				int col_traverse_ind = j;
				for (int p = 0; p < n_minus_k; ++p) {
					c0_Hs[p] += Hs[col_traverse_ind];		// add accumulately, this cause my second problem, mark down here
					col_traverse_ind += n;
				}
			}

		}

		TEP_num++;
	}

	void re_encoding(const Matrix<int>& TEP_now) {
		// re-encode to ce_Hs

		// initialize ans
		for (int i = 0; i < n; ++i) {
			ce_Hs[i] = c0_Hs[i];
		}

		for (int j = 0; j < TEP_now.size(); ++j) {
			// flip systematic part
			int col_ind = n - 1 - TEP_now[j];
			ce_Hs[col_ind] += 1;

			// re-encode parity-check part
			int col_traverse_ind = col_ind;
			for (int p = 0; p < n_minus_k; ++p) {
				ce_Hs[p] += Hs[col_traverse_ind];		// add accumulately, this cause my second problem, mark down here
				col_traverse_ind += n;
			}
		}

		TEP_num++;
	}

	my_double soft_Hamming_weight(const Matrix<GF2>& c_Hs) {
		my_double ans = 0;
		error_num = 0;
		for (int i = 0; i < n; ++i) {
			if (c_Hs[i] != y_Hs[i]) {
				ans += r_abs_Hs[i];
				error_num++;
			}
		}
		return -ans;
	}

	bool check_ML_codeword_general() {
		// Will always check cb_Hs		
		// 'G_threshold(c_Hs,d)' is only needed to compute


		// Only need 'D0_add_size' reliabilities adding up to form the 'G_threshold'
		int D0_add_size = min_Hamming_distance - error_num;

		// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		cb_bar = cb_Hs;
		cb_bar.permute_back(permutation_second);

		G_threshold = 0;
		G_threshold_component.clear();

		for (int i = 0; i < n && D0_add_size > 0; ++i) {

			if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
				G_threshold += r_abs_bar[i];
				D0_add_size--;

				if (w_min < -G_threshold) return true;
				G_threshold_component.push_back(r_abs_bar[i]);
			}
		}
		return false;		// ML criterion for general codeword

		//my_double G_threshold = 0;

		//// Only need 'D0_add_size' reliabilities adding up to form the 'G_threshold'
		//int D0_add_size = min_Hamming_distance - error_num;

		//// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		//cb_Hs.permute_back(permutation_second, cb_bar);

		//for (int i = 0; i < n  && D0_add_size > 0; ++i) {

		//	if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
		//		G_threshold += r_abs_bar[i];
		//		D0_add_size--;
		//	}
		//}
		//return w_min < -G_threshold;		// ML criterion for general codeword
	}

	bool check_ML_codeword_OSD_specific(int TEP_weight) {
		// Will always check cb_Hs		
		// 'G_threshold(c_Hs,d)' is only needed to compute

		G_threshold = G_OSD_threshold;
		//for (int i = 0; i < TEP_weight; ++i) {
		//	G_threshold += r_abs_Hs[n_minus_k + i];
		//}

		// Only need 'D0_add_size' reliabilities adding up to 'G_threshold' further form the entire 'G_threshold'
		int D0_add_size = min_Hamming_distance - error_num - TEP_weight;

		// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		cb_bar = cb_Hs;
		cb_bar.permute_back(permutation_second);

		G_threshold_component.clear();

		for (int i = 0; i < n && D0_add_size > 0; ++i) {

			if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
				G_threshold += r_abs_bar[i];
				D0_add_size--;

				if (w_min < -G_threshold) return true;
				G_threshold_component.push_back(r_abs_bar[i]);
			}
		}
		return false;		// ML criterion for general codeword

		//my_double G_threshold = 0;
		//for (int i = 0; i < TEP_weight; ++i) {
		//	G_threshold += r_abs_Hs[n_minus_k + i];
		//}

		//// Only need 'D0_add_size' reliabilities adding up to 'G_threshold' further form the entire 'G_threshold'
		//int D0_add_size = min_Hamming_distance - error_num - TEP_weight;

		//// This re-arrange the codeword candidate in decreasing reliability. If not, will not guarantee an ML codeword
		//cb_Hs.permute_back(permutation_second, cb_bar);

		//for (int i = 0; i < n && D0_add_size > 0; ++i) {

		//	if (cb_bar[i] == y_bar[i]) {		// 'D0' set, traverse starting from the least reliable position
		//		G_threshold += r_abs_bar[i];
		//		D0_add_size--;
		//	}
		//}
		//return w_min < -G_threshold;		// ML criterion for general codeword
	}

	bool check_ML_codeword_OSD_new_phase(int phase_num) {

		//cout << "r_abs_Hs" << r_abs_Hs;
		//cout << "r_abs_Hs[n_minus_k + phase_num - 1] = " << r_abs_Hs[n_minus_k + phase_num - 1] << endl;

		// update 'G_threshold' for a new judgement
		if (!G_threshold_component.empty()) {
			
			G_threshold += r_abs_Hs[n_minus_k + phase_num - 1] - G_threshold_component.back();
			G_threshold_component.pop_back();
		}
		else {
			G_threshold += r_abs_Hs[n_minus_k + phase_num - 1];
		}

		//cout << "G_threshold_component 2" << G_threshold_component;
		//cout << "G_threshold 2 = " << G_threshold << endl;

		// judge the current best codeword if it is ML
		if (w_min < -G_threshold) {
			return true;
			//return false;
		}
		else {
			// update 'G_OSD_threshold' for the following codeword candidates
			G_OSD_threshold += r_abs_Hs[n_minus_k + phase_num - 1];
			return false;
		}
	}

	void call_GE_only(const Matrix<my_double>& r) {
		for (int i = 0; i < n; ++i) {
			r_abs[i] = -my::abs(r[i]);		// changing sign, therefore, sort with gt results in LRPs first
			y[i] = r[i] > 0 ? 0 : 1;
		}
		permutation_first = Matrix<int>(1, n, 'N');
		r_abs_bar = r_abs;
		r_abs_bar.quick_sort_recur_gt_with_ind(0, n - 1, permutation_first);
		//r_abs.sort_gt_with_record(permutation_first, r_abs_bar);
		y_bar = y;
		y_bar.permute(permutation_first);

		Hs = H;
		Hs.permute_col(permutation_first);

		permutation_second = Matrix<int>(1, n, 'N');
		//Vector_ext::natual<n>(permutation_second);
		Hs.GJE_4_GF2_left(permutation_second);
		//cout << "Hs" << Hs;
	}
};

