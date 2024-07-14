#pragma once
/*****************************************************************//**
 * \file   LC_OSD_r.h
 * \brief  local constrain OSD for r<k
 * 
 * \author lilili
 * \date   April 2023
 *********************************************************************/
#include"reprocess_common.h"
#include"../Viterbi/Viterbi_advanced.h"
using namespace std;

class LC_OSD_r {
public:

	Matrix<GF2> PM;		// we restrict ourself to consider PM only
	Matrix<GF2> right_systematic_PM;
	Matrix<GF2> upper_PM;
	Viterbi_unordered_map vit;
	int d_min;
	int n;
	int k;
	int r;
	int delta;
	unsigned long long total_used_list_num;		// counting used list num
	unsigned long long ML_used_list_num;		
		// while this is the minimum list number required, we are preparing to achieve that
	int times_reach_max_list_num;
	Matrix<GF2> prob_pattern;
	my_double sigma;

	LC_OSD_r(const Matrix<GF2>& _PM, int _d_min, int _delta) :PM(_PM), d_min(_d_min), delta(_delta) {
			 
		n = PM.col();
		k = n - PM.row();
		r = n - k;
		vit.resize(k + delta, delta); 
		total_used_list_num = 0;
		ML_used_list_num = 0;
		times_reach_max_list_num = 0;
		sigma = 1;						// we must correctly set sigma before decoding
	}

	// updated to ref encode, to be fixed
	Matrix<GF2> encode_by_right_systematic_PM(const Matrix<GF2>& information_bit) const{
		Matrix<GF2> ans(1, n);		// unit test to be done
		int k_prime = k + delta;
		int r_prime = r - delta;

		for (int i = k_prime; i < n; ++i) {	// parity check bit
			ans(i) = 0;
		}

		for (int i = 0; i < k_prime; ++i) {	// information bit, systematically encoded
			ans(i) = information_bit(i);
			if (ans(i) != 0) {
				for (int j = 0; j < r_prime; ++j) {	
					// parity check bit, add the column of Patrity check matrix
					ans(j + k_prime) += right_systematic_PM(j, i);
				}
			}
		}
		return ans;
	}

	// the input and output will be ref_codeword 
	void ref_encode_by_right_systematic_PM(Matrix<GF2>& ref_codeword, const Matrix<int>& flip_ind) {

		int k_prime = k + delta;
		int r_prime = r - delta;

		int flip_ind_size = flip_ind.size();
		//cout << "flip_ind_size = " << flip_ind_size << endl;
		for (int i = 0; i < flip_ind_size; ++i) {	// information bit, systematically encoded
			int actual_col = flip_ind(i);		// reverse the flip ind, starting from the least reliable opsition in MRIP
			ref_codeword(actual_col) += 1;
			for (int j = 0; j < r_prime; ++j) {	// parity check bit, add the column of Patrity check matrix
				ref_codeword(j + k_prime) += right_systematic_PM(j, actual_col);
			}
		}
	}

	/**
	 * .get the most reliable independent positions
	 *
	 * \param G: generator matrix from code instance, the process change it to systematic G1
	 * \param r_abs: reliability of each received bit
	 * \return the permutation result
	 */
	Matrix<int> get_MRIP(Matrix<my_double>& r_abs) {
		// include producing systematic G through parity check matrix

		// step 1: turn r to abs(r) and arrange r_abs from big to small

		/*cout << "G_or_H=" << G_or_H;*/
		Matrix<int> Pi = r_abs.sort_with_ind('>');
		//r_abs.permute_back(Pi);		// this process return the original un-permuted r_abs

		right_systematic_PM = PM;

		right_systematic_PM.permute_col(Pi);
		//cout << "Pi=" << Pi;

		// step2: Gaussian cancellation to find independent positions
		right_systematic_PM.row_transformation_to_low_triangle();

		//cout << "right_systematic_PM=" << right_systematic_PM;
		//cout << "Pi=" << Pi;

		Matrix<int> nPi = right_systematic_PM.col_permute_to_full_rank_on_right();
		Pi.permute(nPi);
		r_abs.permute(nPi);

		//cout << "right_systematic_PM=" << right_systematic_PM;
		//cout << "Pi=" << Pi;
		right_systematic_PM.row_transformation_right_low_triangle_to_identity();

		//cout << "right_systematic_PM=" << right_systematic_PM;
		//cout << "Pi=" << Pi;

		nPi = Matrix<int>(1, k, 'N');
		for (int i = 1; i < k; ++i) {		// for k element in r_abs
			int j = i;
			while (j != 0 && r_abs(j) > r_abs(j - 1)) {		// bubble sort
				r_abs.switch_ele(j, j - 1);
				nPi.switch_ele(j, j - 1);
				//G_or_H.switch_col(j, j - 1);
				j--;
			}
		}
		Pi.permute(nPi);							// partial permute
		right_systematic_PM.permute_col(nPi);

		//cout << "right_systematic_PM=" << right_systematic_PM;
		//cout << "Pi=" << Pi;

		// by now the first k col of G are most reliable independent positions
		// that is the first k element as index in r, are most reliable independent positions
		// Pi is the record of perbutation

		return Pi;
	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by flipping all possible <=order bits in k most reliable independent position
	 *
	 * \param recv: received codeword, after BPSK_demodulation
	 * \param beta: for stopping rule, smaller beta, less constrain, less candidate codeword, earlier termination
	 * \param max_list_num: the max list number for list Verterbi decoding
	 * \return decoded result, u
	 */
	Matrix<GF2> decode_v(Matrix<my_double> recv, int max_list_num) {

#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		int double_ope_num_before = my_double_auxiliary_storage::operation_number;
		int double_ope_num_after;
#endif // use_my_double

		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		bool isby_G = true;		// fix for compute dependent column

		Matrix<my_double> r_abs = recv.get_abs();
		Matrix<int> Pi = get_MRIP(r_abs);
		//cout << "Pi=" << Pi; cout << "r=" << r;
		recv.permute(Pi);
		//cout << "recv=" << recv;
		//Matrix<GF2> hdr = BPSK::demodulation(recv);
		//cout << "hdr=" << hdr;

		//cout << "right_systematic_PM" << right_systematic_PM;
		
		// partition the Parity Matrix, taking out the upper delta rows
		upper_PM = right_systematic_PM.get_part(0, 0, delta - 1, k + delta - 1);
		//cout << "upper_PM" << upper_PM;				// we have to cope with all zero columns

		// shrink the right_systematic_PM
		right_systematic_PM = right_systematic_PM.get_part(delta, 0, -1, -1);
		//cout << "right_systematic_PM" << right_systematic_PM;

		vit.change_PM(upper_PM);

		Matrix<GF2> best_dv = vit.decode_v_4_LC_OSD_once(recv, *this, max_list_num);

		best_dv.permute_back(Pi);

#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return best_dv;
	}

	void set_R_pattern(const Matrix<GF2>& _prob_pattern) {
		// this method is not greater than computing expectation
		prob_pattern = _prob_pattern;
	}
	void set_sigma(my_double _sigma) {
		sigma = _sigma;					// we must correctly set sigma before decoding
	}
};

template<class T>
Matrix<GF2> Viterbi_unordered_map::decode_v_4_LC_OSD_once(Matrix<T> r_or_hdr, LC_OSD_r& lc_osd_r, int max_list_num) {

#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

	int double_ope_num_before = my_double_auxiliary_storage::operation_number;
	int double_ope_num_after;
#endif // use_my_double

	int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
	int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

	Matrix<my_double> r_L = r_or_hdr.get_part(0, 0, 0, lc_osd_r.k + lc_osd_r.delta - 1);
	Matrix<my_double> r_R = r_or_hdr.get_part(0, lc_osd_r.k + lc_osd_r.delta, 0, -1);
	//cout << "r_L" << r_L;
	decode_v_once_inner(r_L, max_list_num);

	relative_metric.resize(1, max_list_num, false);

	can_v_hat = list_v_hat.get_row(0);
	can_v_hat.permute_back(PM_permutation);
	//cout << "can_v_hat" << can_v_hat;

	// set for listing
	//Matrix<GF2> hdr = BPSK::demodulation(r_or_hdr);


	/*cout << "r_or_hdr" << r_or_hdr;
	cout << "r_L" << r_L;
	cout << "hdr" << hdr;*/

	// it seems problematic
	/*my_double gamma_Z = -Measure::Euclidean_distance(r_part, BPSK::modulation(hdr.get_part(0, 0, 0, lc_osd_r.k + lc_osd_r.delta - 1)));
	my_double gamma_ZL = -Measure::Euclidean_distance(r_or_hdr.get_part(0, lc_osd_r.k + lc_osd_r.delta, 0, -1), \
			BPSK::modulation(hdr.get_part(0, lc_osd_r.k + lc_osd_r.delta, 0, -1)));*/

	//my_double d_Z = Measure::Euclidean_distance(r_or_hdr, BPSK::modulation(hdr));
	//my_double d_ZL = Measure::Euclidean_distance(r_L, BPSK::modulation(hdr.get_part(0, 0, 0, lc_osd_r.k + lc_osd_r.delta - 1)));

	int r_size = r_R.size();
	Matrix<my_double> r_R_abs = r_R.get_abs();
	my_double r_abs_minus_1_squared = 0;
	my_double tmp_calculation;
	for (int i = 0; i < r_size; ++i) {
		tmp_calculation = r_R_abs(i) - 1;
		r_abs_minus_1_squared += tmp_calculation * tmp_calculation;
	}
	//cout << "r_abs_minus_1_squared = " << r_abs_minus_1_squared << endl;

	// sometimes use d_ZR to ensure ML can greatly improve performance
	
	//Matrix<GF2> hdr_R = BPSK::demodulation(r_R);
	/*if (lc_osd_r.prob_pattern.size() != 0) {
		hdr_R = hdr_R + lc_osd_r.prob_pattern;
	}*/
	//my_double d_ZR = Measure::Euclidean_distance(r_R, hdr_R, r_abs_minus_1_squared);

	// replace it to E{d_ZR}, reference: 'A Low-Complexity Ordered Statistic Decoding of Short Block Codes'
	
	my_double d_ZR = 0;
	my_double scale_by_sigma = 2 / (lc_osd_r.sigma * lc_osd_r.sigma);
	for (int i = 0; i < r_size; ++i) {
		d_ZR += r_R_abs(i) / (1 + exp(scale_by_sigma * r_R_abs(i)));
	}
	d_ZR *= 4;
	d_ZR += r_abs_minus_1_squared;

	//cout << "d_ZR = " << d_ZR << endl;

	//cout << "beta_gamma_ZR = " << beta_gamma_ZR << endl;
	//cout << "beta_minus_1_gamma_ZL = " << beta_minus_1_gamma_ZL << endl;
	
	//cout << "r_part" << r_part;

	// first not flip any position, for the first set we do demand in v space
	Matrix<GF2> best_dv = lc_osd_r.encode_by_right_systematic_PM(can_v_hat);
	int best_dv_list_num = 1;
	Matrix<GF2> first_dv = best_dv;
	//cout << "best_dv: 0" << best_dv << endl;

	// best_dv after permute back must in v space
	my_double lambda_min = absolute_metric + absolute_zero_col_mertric + \
		Measure::Euclidean_distance(r_R, best_dv.get_part(0, lc_osd_r.k + lc_osd_r.delta, 0, -1), r_abs_minus_1_squared);
	//my_double lambda_min = Measure::Euclidean_distance(r_or_hdr, BPSK::modulation(best_dv));
	//cout << "lambda_min = " << lambda_min << endl;
	
	//my_double gamma_C_opt = -lambda_min;
	//cout << "beta_gamma_C_opt = " << beta_gamma_C_opt << endl;
	my_double d_copt = lambda_min;

	/* assume that if PM is full rank, then there is no need to use unused_PM */

	Matrix<GF2> dv;		// for ref encode
	//my_double d_CL_old = Measure::Euclidean_distance(r_L, BPSK::modulation(can_v_hat));
	//cout << "d_CL_old = " << d_CL_old << endl;

	Matrix<int> flip_ind(1, n_total);
	int used_list_num = 1;
	if (max_list_num > 1) {
		relative_metric(0) = 0;
		//can_opt = priority_queue<metric_point>();
		list_proceed = 0;
		diverge_time = n - 1;
		can_opt_set.clear();
		add_zero_col_pos_2_can_opt_set(max_list_num - 1);

		do {
			next_subopt_v_set(max_list_num - used_list_num);

			//my_double d_CL = Measure::Euclidean_distance(r_L, BPSK::modulation(can_v_hat)); 
			my_double d_CL = relative_metric(used_list_num) + absolute_metric + absolute_zero_col_mertric;

			// this is ML		
			if (d_copt < d_CL + d_ZR) {
				/*cout << "break here, i=" << i << endl;
				cout << "d_copt = " << d_copt << endl;
				cout << "d_CL = " << d_CL << endl;*/

				lc_osd_r.total_used_list_num += (unsigned long long)used_list_num + 1;
				lc_osd_r.ML_used_list_num += best_dv_list_num;
				return best_dv;
			}

			can_v_hat = list_v_hat.get_row(used_list_num);
			can_v_hat.permute_back(PM_permutation);

			//cout << "d_CL = " << d_CL << endl;
			//if (d_CL < d_CL_old) {
			//	// we assume that by list Viterbi, we have d_CL increase as list output
			//	cout << "error of list Viterbi" << endl;
			//	cout << "d_CL = " << d_CL << endl;
			//	cout << "d_CL_old = " << d_CL_old << endl;
			//}
			//d_CL_old = d_CL;

			//dv = lc_osd_r.encode_by_right_systematic_PM(can_v_hat);
			 
			dv = first_dv;
			can_v_hat.diff_ind_inplace(dv, flip_ind);
			lc_osd_r.ref_encode_by_right_systematic_PM(dv, flip_ind);		// alright

			//my_double gamma_Cr = -Measure::Euclidean_distance(r_part, BPSK::modulation(dv));
			//my_double d_C = lambda;

			/*cout << "-------- i = " << i << "----------" << endl;
			cout << "lambda = " << lambda << endl;
			cout << "beta_gamma_C_opt = " << beta_gamma_C_opt << endl;
			cout << "gamma_Cr = " << gamma_Cr << endl;
			cout << "beta_gamma_C_opt + beta_gamma_ZR + beta_minus_1_gamma_ZL = "
				<< beta_gamma_C_opt + beta_gamma_ZR + beta_minus_1_gamma_ZL << endl;*/

			//if (gamma_Cr + gamma_ZL < beta * (gamma_C_opt + gamma_Z)) {
			//	//cout << "break here" << endl;
			//	return best_dv;
			//}

			/*if (d_C -d_Z > beta * (d_copt - d_Z)) {
				cout << "break here, i=" << i << endl;
				lc_osd_r.total_used_list_num += (unsigned long long)i + 1;
				return best_dv;
			}*/

			// d_CL + d_ZR is the lower bound of the Euclidean distance of unexplored codewords
			/*if (lambda < d_CL + d_ZR - 1e-6) {
				cout << "error detected" << endl;
				cout << "lambda = " << lambda << endl;
				cout << "d_CL + d_ZR = " << d_CL + d_ZR << endl;
			}*/

			my_double lambda = d_CL + \
				Measure::Euclidean_distance(r_R, dv.get_part(0, lc_osd_r.k + lc_osd_r.delta, 0, -1), r_abs_minus_1_squared);
			if (lambda >= lambda_min);
			else {
				// update best_dv to this erasure decoding result
				best_dv = dv;
				lambda_min = lambda;
				d_copt = lambda_min;
				best_dv_list_num = used_list_num + 1;
			}
			/* assume that if PM is full rank, then there is no need to use unused_PM*/
			
			used_list_num++;
		} while (used_list_num < max_list_num);
	}

	//cout << "relative_metric" << relative_metric;

#ifdef RUN_MSG
	/*cout << "can_opt.size() = " << can_opt.size() << endl;
	while (!can_opt.empty()) {
		cout << can_opt.top();
		can_opt.pop();
	}
	cout << endl;*/

	//cout << "list_v_hat" << list_v_hat;
	cout << "valid_list_ind" << valid_list_ind;
	cout << "relative_metric" << relative_metric;

#ifdef count_operation_number
#ifdef use_my_double

	double_ope_num_after = my_double_auxiliary_storage::operation_number;
	cout << "(Viterbi_unordered_map) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

	GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
	cout << "(Viterbi_unordered_map) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

	/*if (valid_list_ind(0) != 0) {
		cout << "select list change the result" << endl;
	}*/
	//cout << "valid_list_ind" << valid_list_ind;

	lc_osd_r.total_used_list_num += max_list_num;
	lc_osd_r.ML_used_list_num += best_dv_list_num;
	lc_osd_r.times_reach_max_list_num++;
	return best_dv;
}
