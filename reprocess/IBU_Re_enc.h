#pragma once

/*****************************************************************//**
 * \file   IBUc_Re_enc.h
 * \brief  While Updating generator matrices, re-encode the TEPs to generate codeword candidates.
 *		   Select codeword candidate using minimum correlation distance rule and half Hamming distance rule
 *
 * \author 26259
 * \date   July 2024
 *********************************************************************/

#include"reprocess_common.h"

// highly similar to IBU, 's' is max power of U that not equal to identity Matrix, 
// where U: i->(2i) mod n. For BCH codes, s=log2(n+1).

class IBUc_Re_enc {		// delete this, to be done
protected:

	/* ------ Inner variables ----------- */
	int n;
	int k;
	int s;

	Matrix<my_double> r_abs_bar_B_temp;
	Matrix<int> r_abs_bar_B_temp_record;

	/* ------ B ----------- */

	sList<my_double> r_abs_Gs_B_sList;	// sList version, for 'IBU_core', use a dummy head would be simpler
	Matrix<my_double> r_abs_Gs_P;		// vector version, for 'IBU_core_v2', erasing element take O(n) but not a matter
	Matrix<int> mu_candidate;			// vector version, for 'IBU_core_v2'

	//Vector<my_double, k> r_abs_Gs_B_temp;
	Matrix<int> r_abs_Gs_B_temp_record;			// elements added to the original basis
	Matrix<int> r_abs_Gs_B_temp_record_partial_sort;

	Matrix<int> r_abs_bar_record_inv;
	Matrix<int> blank;							// a temporary vector, used to fill MRIPs to generate sorted MRIPs
	Matrix<GF2> occupation_mark;
	Matrix<int> partial_rank_mark;

	/* ------ Permutation selection ----------- */

	Matrix<int> U_pow;
	Matrix<GF2> MRPs_occupy;
	Matrix<int> overlap_max;
	Matrix<int> shift_max;

	int s_in_use;
	bool is_dual;
	bool is_cycle_used;

public:

	/* ------ Prerequisite ----------- */
	//Vector<int, n> pre_permute_record;			// IBU use a 'natual' as default, since 'Gs' un-permuted

	/* ------ Output ----------- */
	Matrix<GF2> Gs;
	Matrix<my_double> r_abs;
	Matrix<my_double> r_abs_bar;
	Matrix<my_double> r_abs_Gs;

	Matrix<int> r_abs_bar_record;
	Matrix<int> r_abs_Gs_record;					// 'Gs' column swap record
	Matrix<int> permute_record_MRIP;				// TEP flip order, and early termination add order

	Matrix<int> r_abs_Gs_P_temp_record;
	Matrix<int>  r_abs_Gs_P_temp_record_partial_sort;
	Matrix<int> r_abs_Gs_P_bar_record;
	Matrix<int> permute_record_LRP;					// TEP flip order, and early termination add order, using IBUc dual

	Matrix<my_double> r_abs_Gs_B_bar;
	Matrix<int> r_abs_Gs_B_bar_record;				// the MRIPs with gt sort, for early termination

	// new added for re-encoding during the IBU

	Matrix<GF2> y;						// Hard-decisions of received vector
	Matrix<GF2> y_Gs;					// Permuted hard-decisions of received vector
	Matrix<GF2> c0_Gs;					// First permuted codeword candidate
	Matrix<GF2> cb_Gs;					// The best permuted codeowrd candidate

	Matrix<GF2> codeword_candidates;				// store the present re-encoding codeword candidate
	Matrix<GF2> best_codeword_candidate;
	my_double min_correlation_distance;

	IBUc_Re_enc(int _n, int _k, int _s) {
		n = _n;
		k = _k;
		s = _s;

		r_abs.resize(1, n);
		r_abs_Gs_record.resize(1, n);
		r_abs_bar_B_temp_record.resize(1, n);

		r_abs_bar_record_inv.resize(1, n);
		blank.resize(1, n);
		occupation_mark.resize(1, n);
		partial_rank_mark.resize(1, n);

		r_abs_Gs_B_temp_record.resize(1, k);
		r_abs_Gs_B_temp_record_partial_sort.resize(1, k);
		r_abs_Gs_B_bar_record.resize(1, k);
		permute_record_MRIP.resize(1, k);					// should be put to the initializing function

		r_abs_Gs_P_temp_record.resize(1, n - k);
		r_abs_Gs_P_temp_record_partial_sort.resize(1, n - k);
		r_abs_Gs_P_bar_record.resize(1, n - k);
		permute_record_LRP.resize(1, n - k);

		// U^i stored in Vector of Vector
		U_pow.resize(s, n);
		MRPs_occupy.resize(s, n + k);
		for (int i = 0; i < n; ++i) {
			U_pow[i] = i;
		}
		for (int i = 1; i < s; ++i) {
			int i_row_start_ind = i * n;
			int i_minus_one_row_start_ind = i_row_start_ind - n;

			for (int j = 0; j < n; ++j) {
				U_pow[i_row_start_ind + j] = (U_pow[i_minus_one_row_start_ind + j] * 2) % n;
			}
		}

		//cout << "U_pow" << U_pow;

		//cout << "MRPs_occupy" << MRPs_occupy;

		overlap_max.resize(1, s);
		shift_max.resize(1, s);

		s_in_use = s;
		is_dual = false;
		is_cycle_used = true;
	}
	void reset_s(int _s) {
		s_in_use = _s;
	}
	void set_dual() {
		is_dual = true;
	}
	void set_cycle_used(bool _is_cycle_used) {
		is_cycle_used = _is_cycle_used;
	}

	void solve(const Matrix<GF2>& Gs_ini, const Matrix<my_double>& r) {
		// determine the 'pre_permute_record' for IBU to solve

		Gs = Gs_ini;

		if (is_dual == false) {
			for (int i = 0; i < n; ++i) {
				r_abs[i] = my::abs(r[i]);
			}
		}
		else {
			for (int i = 0; i < n; ++i) {
				r_abs[i] = -my::abs(r[i]);		// find LRB for dual code, use negative values
			}
		}

		//cout << "r_abs" << r_abs;
		//cout << "pre_permute_record" << pre_permute_record;

		r_abs_bar = r_abs;
		r_abs_bar_record = Matrix<int>(1, n, 'N');

		//cout << "r_abs_bar_record" << r_abs_bar_record;

		r_abs_bar.quick_sort_recur_gt_with_ind(0, n - 1, r_abs_bar_record);

		//cout << "r_abs_bar" << r_abs_bar;
		//cout << "r_abs_bar_record" << r_abs_bar_record;

		// write the sorted index into blank
		for (int i = 0; i < n; ++i) {
			r_abs_bar_record_inv[r_abs_bar_record[i]] = i;
		}
		//cout << "r_abs_bar_record_inv" << r_abs_bar_record_inv;

		if (is_cycle_used) {
			determine_best_pre_permutation();				// the best pre-permutation is generated as 'r_abs_Gs_record'
		}
		else {
			r_abs_Gs_record = Matrix<int>(1, n, 'N');		// Vector_ext::natual<n>();
		}

		//Vector<int, n> pre_permute_record;							// to be updated
		//Vector_ext::natual(pre_permute_record);

		r_abs_bar_B_temp = r_abs;
		r_abs_bar_B_temp.permute(r_abs_Gs_record);			// permute with selected pre-permutation that lead to good basis

		// generate the re-encoding codeword for order-1 re-encoding of the current matrix
		for (int i = 0; i < n; ++i) {
			y[i] = r[i] > 0 ? 0 : 1;
		}
		y_Gs = y;
		y_Gs.permute(r_abs_Gs_record);						// ***** new *****

		//partial_gt_sort_k_and_n_k();
		//cout << "r_abs_Gs_record" << r_abs_Gs_record;
		partial_gt_sort_split(r_abs_Gs_record, r_abs_bar_B_temp_record);
			// 'r_abs_bar_B_temp_record' is the permutation of split sort

		//cout << "r_abs_bar_B_temp_record" << r_abs_bar_B_temp_record;

		// after permutation, 
		// the first 'k' part and the last 'n-k' part of 'r_abs_bar_B_temp' are sorted in decending order (gt order)
		r_abs_bar_B_temp.permute(r_abs_bar_B_temp_record);
		// this is permute back, not permute, wasting me half an hour

		//cout << "r_abs_bar_B_temp" << r_abs_bar_B_temp;
		//cout << "r_abs_bar_B_temp_record" << r_abs_bar_B_temp_record;


		// IBU process

		IBU_core_v2();

		// get 'r_abs_Gs'
		r_abs_Gs = r_abs;
		r_abs_Gs.permute(r_abs_Gs_record);

		// sort MRB for TEP re-encoding, Interface as 'permute_record_MRIP'

		if (is_dual == false) {
			sort_MRB_for_TEP_re_encoding();
		}
		else {
			sort_LRB_for_TEP_re_encoding();
		}
	}

	void determine_best_pre_permutation() {
		// output as 'r_abs_Gs_record'

		MRPs_occupy.reset(0);
		overlap_max.reset(0);

		// compute the overlap between candidate basis and MRPs
		for (int i = 0; i < s_in_use; ++i) {
			int i_row_start_ind = i * (n + k);

			// seperate identity and shifting, Mark the MRPs in each permutation class
			int s_minus_i_row_start_ind = (i == 0) ? 0 : (s - i) * n;
			for (int j = 0; j < k; ++j) {
				MRPs_occupy[i_row_start_ind + U_pow[s_minus_i_row_start_ind + r_abs_bar_record[j]]] = 1;
			}

			for (int j = 0; j < k; ++j) {
				if (MRPs_occupy[i_row_start_ind + j] == 1) {

					// cyclic copy the first k-1 position onto the tail of 'MRPs_occupy'
					MRPs_occupy[i_row_start_ind + n + j] = 1;

					// count the first 'k' overlap on MRPs for each permutation class
					overlap_max[i]++;

				}
			}
			shift_max[i] = 0;

			int overlap_can = overlap_max[i];
			//bool fomer_plus = false;

			// count the cyclic shifting overlap and keep the max
			int j_minus_k_store = 0;
			for (int j = k; j < n + k - 1; ++j, ++j_minus_k_store) {

				overlap_can += (int)MRPs_occupy[i_row_start_ind + j] - (int)MRPs_occupy[i_row_start_ind + j_minus_k_store];
				if (overlap_max[i] < overlap_can) {

					overlap_max[i] = overlap_can;
					shift_max[i] = j_minus_k_store + 1;
				}
			}
		}

		//cout << "MRPs_occupy" << MRPs_occupy; 
		//cout << "overlap_max" << overlap_max;
		//cout << "shift_max" << shift_max;				// to be simplified to record only one max overlap and max shift

		// determine the max overlap and max shift

		// find the best permutation U^i T^j
		int permutation_opt_ind = 0;
		int overlap_max_final = overlap_max[0];
		for (int i = 1; i < s_in_use; ++i) {
			if (overlap_max_final < overlap_max[i]) {
				permutation_opt_ind = i;
				overlap_max_final = overlap_max[i];		// can be simplified
			}
		}

		// determine the pre-permutation_record

		int permutation_opt_ind_row_start_ind = permutation_opt_ind * n;
		int shift_bias = permutation_opt_ind_row_start_ind + shift_max[permutation_opt_ind];
		int tail_copy = n - shift_max[permutation_opt_ind];

		for (int i = 0; i < tail_copy; ++i) {
			r_abs_Gs_record[i] = U_pow[shift_bias + i];
		}
		for (int i = tail_copy, j = 0; i < n; ++i, ++j) {
			r_abs_Gs_record[i] = U_pow[permutation_opt_ind_row_start_ind + j];
		}
		//pre_permute_record.shift_right_cir(-shift_max[permutation_opt_ind]);

		//cout << "r_abs_Gs_record" << r_abs_Gs_record;
	}

	void partial_gt_sort_split(const Matrix<int>& input, Matrix<int>& output) {

		// Input as 'r_abs_Gs_record' and output as 'r_abs_bar_B_temp_record'

		// sort the first 'k' column indices of Gs, as well as the last 'n - k' column indices 

		// read the sorting order of first 'k' column indices of Gs and write occupation_mark
		occupation_mark.reset(0);
		for (int i = 0; i < k; ++i) {
			occupation_mark[r_abs_bar_record_inv[input[i]]] = 1;
		}
		//cout << "occupation_mark" << occupation_mark;

		// bring down the 'r_abs_Gs_B_record' into a 'k' length permutation
		int lt_k_num = 0;

		// bring down the 'r_abs_Gs_B_record' into a 'k' to 'n' permutation
		int ge_k_lt_n_num = k;

		for (int i = 0; i < n; ++i) {
			if (occupation_mark[i] == 1) {
				partial_rank_mark[i] = lt_k_num;
				lt_k_num++;
			}
			else {
				partial_rank_mark[i] = ge_k_lt_n_num;
				ge_k_lt_n_num++;
			}
		}
		//cout << "occupation_mark" << occupation_mark;

		// get 'r_abs_Gs_B_gt_record' and 'r_abs_Gs_P_gt_record', all in one 'r_abs_bar_B_temp_record'
		for (int i = 0; i < n; ++i) {
			//B_P_sort[i] = occupation_mark[blank[r_abs_Gs_record[i]]];
			output[partial_rank_mark[r_abs_bar_record_inv[input[i]]]] = i;	// find inverse permutation
		}
		//cout << "r_abs_bar_B_temp_record" << r_abs_bar_B_temp_record;

		// --------------------------------------
	}

	void partial_gt_sort_first_k(const Matrix<int>& input, Matrix<int>& output) {

		// Input as 'r_abs_Gs_record' and output as 'r_abs_bar_B_temp_record'

		// sort the first 'k' column indices of Gs

		// read the sorting order of first 'k' column indices of Gs and write occupation_mark
		occupation_mark.reset(0);
		for (int i = 0; i < k; ++i) {
			occupation_mark[r_abs_bar_record_inv[input[i]]] = 1;
		}
		//cout << "occupation_mark" << occupation_mark;

		// bring down the 'r_abs_Gs_B_record' into a 'k' length permutation
		int lt_k_num = 0;

		for (int i = 0; i < n; ++i) {
			if (occupation_mark[i] == 1) {
				partial_rank_mark[i] = lt_k_num;
				lt_k_num++;
			}
		}
		//cout << "occupation_mark" << occupation_mark;

		// get 'r_abs_Gs_B_gt_record' and 'r_abs_Gs_P_gt_record', all in one 'r_abs_bar_B_temp_record'
		for (int i = 0; i < k; ++i) {
			//B_P_sort[i] = occupation_mark[blank[r_abs_Gs_record[i]]];
			output[partial_rank_mark[r_abs_bar_record_inv[input[i]]]] = i;	// find inverse permutation
		}
		//cout << "r_abs_bar_B_temp_record" << r_abs_bar_B_temp_record;

		// --------------------------------------
	}

	void partial_gt_sort_last_n_k(const Matrix<int>& input, Matrix<int>& output) {

		// Input as 'r_abs_Gs_record' and output as 'r_abs_bar_B_temp_record'

		// sort the first 'k' column indices of Gs

		// read the sorting order of first 'k' column indices of Gs and write occupation_mark
		occupation_mark.reset(0);
		for (int i = 0; i < n - k; ++i) {
			occupation_mark[r_abs_bar_record_inv[input[i]]] = 1;
		}
		//cout << "occupation_mark" << occupation_mark;

		// bring down the 'r_abs_Gs_B_record' into a 'k' length permutation
		int lt_k_num = 0;

		for (int i = 0; i < n; ++i) {
			if (occupation_mark[i] == 1) {
				partial_rank_mark[i] = lt_k_num;
				lt_k_num++;
			}
		}
		//cout << "occupation_mark" << occupation_mark;

		// get 'r_abs_Gs_B_gt_record' and 'r_abs_Gs_P_gt_record', all in one 'r_abs_bar_B_temp_record'
		for (int i = 0; i < n - k; ++i) {
			//B_P_sort[i] = occupation_mark[blank[r_abs_Gs_record[i]]];
			output[partial_rank_mark[r_abs_bar_record_inv[input[i]]]] = i;	// find inverse permutation
		}
		//cout << "r_abs_bar_B_temp_record" << r_abs_bar_B_temp_record;

		// --------------------------------------
	}

	// not used function, to be deleted
	void partial_gt_sort_k_and_n_k() {

		// Input as 'r_abs_Gs_record' and output as 'r_abs_bar_B_temp_record'

		// sort the first 'k' column indices of Gs

		// write the sorted index into blank
		for (int i = 0; i < n; ++i) {
			blank[r_abs_bar_record[i]] = i;
		}
		//cout << "blank" << blank;

		// read the sorting order of first 'k' column indices of Gs and write occupation_mark
		occupation_mark.reset(0);
		for (int i = 0; i < k; ++i) {
			occupation_mark[blank[r_abs_Gs_record[i]]] = 1;
		}
		//cout << "occupation_mark" << occupation_mark;

		// bring down the 'r_abs_Gs_B_record' into a 'k' length permutation
		int lt_k_num = 0;

		// bring down the 'r_abs_Gs_B_record' into a 'k' to 'n' permutation
		int ge_k_lt_n_num = k;

		for (int i = 0; i < n; ++i) {
			if (occupation_mark[i] == 1) {
				partial_rank_mark[i] = lt_k_num;
				lt_k_num++;
			}
			else {
				partial_rank_mark[i] = ge_k_lt_n_num;
				ge_k_lt_n_num++;
			}
		}
		//cout << "occupation_mark" << occupation_mark;

		// get 'r_abs_Gs_B_gt_record' and 'r_abs_Gs_P_gt_record', all in one 'r_abs_bar_B_temp_record'
		for (int i = 0; i < n; ++i) {
			//B_P_sort[i] = occupation_mark[blank[r_abs_Gs_record[i]]];
			r_abs_bar_B_temp_record[partial_rank_mark[blank[r_abs_Gs_record[i]]]] = i;	// find inverse permutation
		}
		//cout << "r_abs_bar_B_temp_record" << r_abs_bar_B_temp_record;

		// --------------------------------------
	}

	void IBU_core_v2() {
		// initialize 'r_abs_Gs_P'
		int n_minus_k = n - k;

		// initialize 'mu_candidate'
		mu_candidate.resize(1, n_minus_k);
		for (int i = 0; i < n_minus_k; ++i) {
			mu_candidate[i] = n - 1 - i;
		}

		int lambda = k - 1;		// ranging from 0,1, ... , k-1

		while (lambda>=0 && mu_candidate.empty()==false && r_abs_bar_B_temp[lambda] < r_abs_bar_B_temp[mu_candidate.back()])
		{
			GF2_auxiliary_storage::iteration_number++;
			int lambda_real = r_abs_bar_B_temp_record[lambda];
				// corresponding to the partial order 'r_abs_bar_B_temp', range unchanged

			int mu_candidate_tail_ind = mu_candidate.size() - 1;
			int mu = mu_candidate[mu_candidate_tail_ind];			// ranging from k, k+1, ... ,n-1
			int mu_real = r_abs_bar_B_temp_record[mu];
				// corresponding to the partial order 'r_abs_bar_B_temp', range unchanged
				
			// identifying mu by making sure Gs(lambda_real, mu_real) == 1
			while (Gs(lambda_real, mu_real) == 0 && mu_candidate_tail_ind != 0) {

				// search for the next mu
				mu_candidate_tail_ind--;
				mu = mu_candidate[mu_candidate_tail_ind];
				mu_real = r_abs_bar_B_temp_record[mu];			// update mu_real
			}

			if (Gs(lambda_real, mu_real) == 1) {
				// mu is found

				if (r_abs_bar_B_temp[lambda] < r_abs_bar_B_temp[mu]) {
					// do basis update

					// Record permutation
					r_abs_Gs_record.switch_ele(lambda_real, mu_real);

					// update Gs
					int lambda_row_start_ind = lambda_real * n;				// optimize for key code
					Gs[lambda_row_start_ind + mu_real] = 0;

					for (int i = 0; i < k; ++i) {
						int i_row_start_ind = i * n;
						if (Gs[i_row_start_ind + mu_real] == 1) {	// if (Gs(i,mu)==1)

							// this is not needed for now
							GF2_auxiliary_storage::GE_bit_plane_number++;
							GF2_auxiliary_storage::GE_bit_plane_norm_number += is_dual == false ? 1 : (double)(n - k) / k;

							for (int j = k; j < n; ++j) {
								// add row-lamda of P to row-i of P, where P is parity submatrix of 'Gs'
								Gs[i_row_start_ind + j] += Gs[lambda_row_start_ind + j];
							}
						}
					}
					Gs[lambda_row_start_ind + mu_real] = 1;
				}

				// Update non-MRB candidates by removing MRB element
				lambda--;
			}

			// Update MRB candidates by removing non-MRB element
			mu_candidate.erase_element(mu_candidate_tail_ind);
		}
	}

	void sort_MRB_for_TEP_re_encoding() {
		// input as 'r_abs_Gs_B_temp_record' and output as 'permute_record_MRIP'

		//// merge 'r_abs_Gs_B_sList' and 'r_abs_Gs_B_temp' into partial lt sort order as 'r_abs_Gs_B_temp' and 'r_abs_Gs_B_temp_record'
		//// this form the MRB, a merge sort can turn 'r_abs_Gs_B_temp' in full order, for OSD TEP flipping positions
		//// it can also be finished with sorted result of 'r_abs_bar'		

		for (int i = 0; i < k; ++i) {
			r_abs_Gs_B_temp_record[i] = r_abs_Gs_record[i];		// copy the MRB part indices for sorting
		}

		partial_gt_sort_first_k(r_abs_Gs_B_temp_record, r_abs_Gs_B_temp_record_partial_sort);

		// determine the order of MRIPs in Gs

		r_abs_Gs_B_bar_record = r_abs_Gs_B_temp_record;
		r_abs_Gs_B_bar_record.permute(r_abs_Gs_B_temp_record_partial_sort);
		//cout << "r_abs_Gs_B_bar_record (new)" << r_abs_Gs_B_temp_record;

		// fill in the blank
		for (int i = 0; i < k; ++i) {
			blank[r_abs_Gs_record[i]] = i;		// unsorted Gs indices, read the MRB
		}
		// read the blank and form 'permute_record_MRIP'
		for (int i = 0; i < k; ++i) {
			// 'permute_record_MRIP' are the gt sorted MRB indices, ranging from 0 to k-1, corresponds to Gs
			permute_record_MRIP[i] = blank[r_abs_Gs_B_bar_record[i]];
		}

		//cout << "permute_record_MRIP" << permute_record_MRIP;
	}

	void sort_LRB_for_TEP_re_encoding() {

		for (int i = 0; i < n - k; ++i) {
			r_abs_Gs_P_temp_record[i] = r_abs_Gs_record[i + k];		// copy the MRB part indices for sorting
		}

		partial_gt_sort_last_n_k(r_abs_Gs_P_temp_record, r_abs_Gs_P_temp_record_partial_sort);

		// determine the order of MRIPs in Gs

		r_abs_Gs_P_bar_record = r_abs_Gs_P_temp_record;
		r_abs_Gs_P_bar_record.permute(r_abs_Gs_P_temp_record_partial_sort);
		//cout << "r_abs_Gs_B_bar_record (new)" << r_abs_Gs_B_temp_record;

		// fill in the blank
		for (int i = 0; i < n - k; ++i) {
			blank[r_abs_Gs_record[i + k]] = i;		// unsorted Gs indices, read the MRB
		}
		// read the blank and form 'permute_record_MRIP'
		for (int i = 0; i < n - k; ++i) {
			// 'permute_record_MRIP' are the gt sorted MRB indices, ranging from 0 to k-1, corresponds to Gs
			permute_record_LRP[i] = k + blank[r_abs_Gs_P_bar_record[i]];
		}

		//cout << "permute_record_LRP" << permute_record_LRP;
	}
};