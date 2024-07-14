#pragma once
/*****************************************************************//**
 * \file   LL_OSD_Viterbi.h
 * \brief  LL_OSD by bit operation over PM
 *
 * \author lilili
 * \date   April 2023
 *********************************************************************/

#include"reprocess_common.h"
#include"../Viterbi/Viterbi_advanced.h"
#include"../code/RS.h"
#include"../code/BCH.h"
#include"LL_OSD.h"

template<int m, int t, int k_prime, int selected_row_num>
class LL_OSD_Viterbi {
public:

	RS<m, k_prime> rs;
	BCH<m, t> bch;
	Matrix<GF2> extend_generator;
	Viterbi_unordered_map vit;
	int d;
	int n;
	int k;
	int r;
	int r_prime;
	int num_invalid_list;
	unsigned long long total_used_list_num;		// counting used list num
	int num_early_stop;

	// sigma is not included, the stopping criteria is ML criteria, the ave_list_num is much greater than ML_use_list_num
	// can be fix in the future but not meaningful. please see the implementation of LL_OSD_Hybrid_flip_Viterbi.

	LL_OSD_Viterbi() {
		d = bch.get_d();
		n = bch.get_n();
		k = bch.get_k();
		r = n - k;
		r_prime = d - 1;
		num_invalid_list = 0;
		total_used_list_num = 0;
		num_early_stop = 0;

		vit.resize(k_prime, selected_row_num);
		//cout << "k_prime = " << k_prime << endl;
	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by serial list Viterbi algorithm to find candidate codeword
	 *
	 * \param recv: received codeword, after BPSK_demodulation
	 * \param rs: the mother code rs
	 * \param alpha: for stopping rule, smaller alpha, less constrain, less candidate codeword, earlier termination
	 * \param max_list_num: the max list number for list Verterbi decoding
	 * \return decoded result, u
	 */
	Matrix<GF2> decode_v(const Matrix<my_double>& recv, int max_list_num) {

#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		int double_ope_num_before = my_double_auxiliary_storage::operation_number;
		int double_ope_num_after;
#endif // use_my_double

		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;

		int GF2e_ope_num_before = GF2e_auxiliary_storage::operation_number;
		int GF2e_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		//cout << "r" << r;
		Matrix<my_double> r_abs = recv.get_abs();
		Matrix<int> permute_ind = r_abs.sort_with_ind('>');		// from big to small
		//cout << "recv" << recv;
		//cout << "r_abs" << r_abs;

		//cout << "permute_ind" << permute_ind;

		/*cout << "GF2e_auxiliary_storage::operation_number: before generate_systematic_generator_any_pos_best_r= " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/
			// rs generator matrix from most reliable information set of bch received codeword
		rs.generate_systematic_generator_any_pos_best(permute_ind.get_part(0, 0, 0, k_prime - 1));

		/*cout << "GF2e_auxiliary_storage::operation_number: after generate_systematic_generator_any_pos_best_r= " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/

			//cout << "rs.generator_M_systematic_any_pos" << rs.generator_M_systematic_any_pos;
			//cout << "rs.information_set_ind" << rs.information_set_ind;
			// the non-systematic part of rs generator matrix, and transposed, n * r_prime matrix
		Matrix<GF2e<m>> parity_part = rs.generator_M_systematic_any_pos.get_cols(rs.redundancy_set_ind).Transpose();
		//cout << "parity_part" << parity_part;
		// extend GF2e<m> to GF2 and apply to each element of parity_part, generate a (m*r_prime) * n matrix
		Matrix<GF2> extend_parity_part = GF2e<m>::to_bits_row_extention(parity_part);
		//cout << "extend_parity_part" << extend_parity_part;
		// throw away the rows having ending 1 at last r_prime columns
		Matrix<int> erase_row_ind(1, r_prime);
		for (int i = 0; i < r_prime; ++i) {
			erase_row_ind(i) = i * m;
		}
		//cout << "erase_row_ind" << erase_row_ind;
		Matrix<GF2> sub_parity_matrix = extend_parity_part.erase_rows(erase_row_ind);
		//cout << "sub_parity_matrix" << sub_parity_matrix;

		sub_parity_matrix = sub_parity_matrix.erase_all_0_rows();
		//cout << "sub_parity_matrix (erased rows)" << sub_parity_matrix;

		// randomly choose (k_prime - k) line
		Matrix<int> row_ind(1, sub_parity_matrix.row(), 'N');
		Matrix<int> chosen_row_ind = row_ind.get_random_element(selected_row_num);
		//cout << "chosen_row_ind" << chosen_row_ind;

		// by size of (k_prime - k) * k_prime
		Matrix<GF2> chosen_parity_check_matrix = sub_parity_matrix.get_rows(chosen_row_ind);
		//cout << "choosen_parity_check_matrix" << chosen_parity_check_matrix;

		vit.change_PM(chosen_parity_check_matrix);

		// select the valid codeword from list

		//Matrix<GF2> backup_not_processed_parity_check_matrix = sub_parity_matrix.erase_rows(chosen_row_ind);
		Matrix<GF2> not_processed_parity_check_matrix = sub_parity_matrix.erase_rows(chosen_row_ind);
		//cout << "not_processed_parity_check_matrix" << not_processed_parity_check_matrix;

		// compute the inner product of vector in viterbi_list and each row of not_processed_parity_check_matrix
		// if one result is non zero, then the codeword is invalid and should be discarded

		vit.change_unused_PM(not_processed_parity_check_matrix);

		//cout << "rs.information_set_ind" << rs.information_set_ind;	

		//cout << "viterbi_list" << viterbi_list;
		extend_generator = extend_parity_part.get_rows(erase_row_ind).Transpose();
		//cout << "extend_generator" << extend_generator << endl;

		Matrix<GF2> v_hat = vit.decode_v_4_LL_OSD_once(recv, *this, max_list_num);

		if (v_hat.size() == 0) {
			// there is no valid list, decode fail
			//cout << "there is no valid list, i=" << i << endl;
			v_hat = Matrix<GF2>(1, n, '0');
			num_invalid_list++;
		}


#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(LL_OSD_viterbi) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(LL_OSD_viterbi) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << "(LL_OSD_viterbi) GF2e_ope_num = " << GF2e_ope_num_after - GF2e_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return v_hat;
	}
};

template<class T, class LL_OSD_type>
Matrix<GF2> Viterbi_unordered_map::decode_v_4_LL_OSD_once(Matrix<T> recv, LL_OSD_type& ll_osd_Viterbi, int max_list_num) {

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

	//cout << "recv_k_prime_part" << recv_k_prime_part;
	//cout << "n = " << n << endl;

	Matrix<T> recv_k_prime_part = recv.get_cols(ll_osd_Viterbi.rs.information_set_ind);		// r_L
	Matrix<T> recv_redundancy_part = recv.get_cols(ll_osd_Viterbi.rs.redundancy_set_ind);	// r_R

	decode_v_once_inner(recv_k_prime_part, max_list_num);

	//cout << "r_part" << r_part;

	relative_metric.resize(1, max_list_num, false);

	// first not flip any position, for the first set we do demand in v space
	Matrix<GF2> best_dv;
	my_double lambda_min = invalid_large_metric;
	my_double d_copt = lambda_min;
	can_v_hat = list_v_hat.get_row(0);

	//cout << "can_v_hat" << can_v_hat;
	//cout << "unused_PM" << unused_PM;

	// set for listing
	Matrix<GF2> hdr_R = BPSK::demodulation(recv_redundancy_part);
	int r_size = recv_redundancy_part.size();
	Matrix<my_double> recv_redundancy_part_abs = recv_redundancy_part.get_abs();
	my_double r_abs_minus_1_squared = 0;
	my_double tmp_calculation;
	for (int i = 0; i < r_size; ++i) {
		tmp_calculation = recv_redundancy_part_abs(i) - 1;
		r_abs_minus_1_squared += tmp_calculation * tmp_calculation;
	}

	//my_double d_Z = Measure::Euclidean_distance(recv, BPSK::modulation(hdr));
	/*my_double n_plus_sum_r_squared_minus_2_r_abs = r_size;
	for (int i = 0; i < r_size; ++i) {
		n_plus_sum_r_squared_minus_2_r_abs += recv_redundancy_part(i) * recv_redundancy_part(i) - 2 * my::abs(recv_redundancy_part(i));
	}*/
	my_double d_ZR = Measure::Euclidean_distance(recv_redundancy_part, hdr_R, r_abs_minus_1_squared);
	Matrix<GF2> v_redundancy;

	can_first_v_hat.resize(1, 0, false);
	if (unused_PM.check_inner_product_4_GF2(can_v_hat, can_first_v_hat, non_0_ind_aux)) {
		can_first_v_hat = can_v_hat;
		can_v_hat.permute_back(PM_permutation);
		v_redundancy = can_v_hat * ll_osd_Viterbi.extend_generator;
		best_dv = can_v_hat.insert_cols(v_redundancy, ll_osd_Viterbi.rs.redundancy_set_ind);
		lambda_min = absolute_metric + absolute_zero_col_mertric + \
			Measure::Euclidean_distance(recv_redundancy_part, v_redundancy, r_abs_minus_1_squared);
		d_copt = lambda_min;
	}

	//cout << "best_dv: 0" << best_dv << endl;

	// best_dv after permute back must in v space
	//cout << "lambda_min = " << lambda_min << endl;

	//my_double gamma_C_opt = -lambda_min;
	//cout << "alpha_gamma_C_opt = " << alpha_gamma_C_opt << endl;

	/* assume that if PM is full rank, then there is no need to use unused_PM */

	Matrix<GF2> dv;
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

			my_double d_CL = relative_metric(used_list_num) + absolute_metric + absolute_zero_col_mertric;
			if (d_copt < d_CL + d_ZR) {

				ll_osd_Viterbi.total_used_list_num += (unsigned long long)used_list_num + 1;
				ll_osd_Viterbi.num_early_stop++;
				return best_dv;
			}
			can_v_hat = list_v_hat.get_row(used_list_num);

			if (unused_PM.check_inner_product_4_GF2(can_v_hat, can_first_v_hat, non_0_ind_aux)) {

				if (can_first_v_hat.size() == 0) {
					can_first_v_hat = can_v_hat;
				}
				can_v_hat.permute_back(PM_permutation);

				// ref decoding to be developed, but with hybrid decoding, this doesn't matter
				v_redundancy = can_v_hat * ll_osd_Viterbi.extend_generator;
				dv = can_v_hat.insert_cols(v_redundancy, ll_osd_Viterbi.rs.redundancy_set_ind);


				//my_double gamma_Cr = -Measure::Euclidean_distance(r_part, BPSK::modulation(dv));
				//my_double lambda = Measure::Euclidean_distance(recv, BPSK::modulation(dv));
				my_double lambda = d_CL + Measure::Euclidean_distance(recv_redundancy_part, \
					v_redundancy, r_abs_minus_1_squared);
				//my_double d_C = lambda;

				/*cout << "-------- i = " << i << "----------" << endl;
				cout << "d_C = " << d_C << endl;
				cout << "d_copt = " << d_copt << endl;
				cout << "d_C - d_Z = " << d_C - d_Z << " ?> alpha * (d_copt - d_Z) = " << alpha * (d_copt - d_Z) << endl;*/

				//if (gamma_Cr + gamma_ZL < alpha * (gamma_C_opt + gamma_Z)) {
				//	//cout << "break here" << endl;
				//	return best_dv;
				//}

				//if (d_C - d_Z > beta * (d_copt - d_Z)) {
				//	//cout << "break here, i = " << i << endl;
				//	ll_osd_Viterbi.total_used_list_num += used_list_num + 1;
				//	return best_dv;
				//}

				if (lambda >= lambda_min);
				else {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					d_copt = lambda_min;

					//cout << "_______________" << endl;
					//cout << "used_list_num = " << used_list_num << endl;
					//cout << "best_dv" << best_dv;
					//cout << "lambda_min = " << lambda_min << endl;
					//cout << "_______________" << endl;
				}

				/* assume that if PM is full rank, then there is no need to use unused_PM*/
			}
			used_list_num++;
		} while (used_list_num < max_list_num);
		//cout << "relative_metric" << relative_metric;
	}


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

	ll_osd_Viterbi.total_used_list_num += max_list_num;
	return best_dv;
}


template<int m, int t, int k_prime, int selected_row_num>
class LL_OSD_Hybrid_flip_Viterbi {
public:

	BCH<m, t> bch;
	RS<m, k_prime> rs;

	Matrix<my_double> recv_permuted;
	Matrix<GF2> sub_parity_matrix;
	Matrix<GF2> extend_generator;		// note the extend generator is adjusted to follow LL_OSD_Viterbi
	Viterbi_unordered_map vit;

	Matrix<int> permuted_ind;
	Matrix<int> permuted_ind_left;
	Matrix<int> permuted_ind_right;
	Matrix<int> left_sort_ind;
	Matrix<int> right_sort_ind;
	Matrix<GF2> codeword_right;
	int d;
	int n;
	int k;
	int r;
	int r_prime;
	int num_invalid_list;
	int num_early_stop;
	unsigned long long total_used_list_num;		// counting used list num
	bool is_ML;

	unsigned long long ML_used_list_num;		// this is the minimum list number required
	int times_reach_max_list_num;
	my_double sigma;

	LL_OSD_Hybrid_flip_Viterbi() {
		d = bch.get_d();
		n = bch.get_n();
		k = bch.get_k();
		r = n - k;
		r_prime = d - 1;
		num_invalid_list = 0;
		num_early_stop = 0;
		total_used_list_num = 0;
		codeword_right.resize(1, r_prime, false);
		is_ML = false;

		ML_used_list_num = 0;
		times_reach_max_list_num = 0;
		sigma = 1;

		vit.resize(k_prime, selected_row_num);
		//cout << "k_prime = " << k_prime << endl;
	}

	Matrix<GF2> encode_by_right_systematic_PM(const Matrix<GF2>& information_bit) {
		Matrix<GF2> ans(1, n);		// unit test to be done
		codeword_right.reset(0);

		for (int i = 0; i < k_prime; ++i) {	// information bit, systematically encoded
			ans(i) = information_bit(i);
			if (ans(i) != 0) {
				for (int j = 0; j < r_prime; ++j) {	// parity check bit, add the column of Patrity check matrix
					codeword_right(j) += extend_generator(i, j);
				}
			}
		}

		codeword_right.permute_back(right_sort_ind);
		for (int j = 0; j < r_prime; ++j) {
			ans(k_prime + j) = codeword_right(j);
		}

		return ans;
	}

	// the input and output will be ref_codeword 
	void ref_encode_by_right_systematic_PM(Matrix<GF2>& ref_codeword, const Matrix<int>& filp_ind) {
		codeword_right.reset(0);
		//cout << "codeword_right" << codeword_right;
		int flip_ind_size = filp_ind.size();
		//cout << "flip_ind_size = " << flip_ind_size << endl;
		for (int i = 0; i < flip_ind_size; ++i) {	// information bit, systematically encoded
			ref_codeword(filp_ind(i)) += 1;
			for (int j = 0; j < r_prime; ++j) {	// parity check bit, add the column of Patrity check matrix
				codeword_right(j) += extend_generator(filp_ind(i), j);
			}
		}
		/*cout << "extend_generator" << extend_generator;
		cout << "filp_ind" << filp_ind;
		cout << "codeword_right" << codeword_right;
		cout << "right_sort_ind" << right_sort_ind;*/

		codeword_right.permute_back(right_sort_ind);
		for (int j = 0; j < r_prime; ++j) {	// parity check bit, add the column of Patrity check matrix
			ref_codeword(k_prime + j) += codeword_right(j);
		}
	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by flipping all possible <=order bits in k most reliable independent position
	 *
	 * \param recv: received codeword, after BPSK_demodulation
	 * \param rs: the mother code rs
	 * \param alpha: for stopping rule, smaller alpha, less constrain, less candidate codeword, earlier termination
	 * \param max_list_num: the max list number for list Verterbi decoding
	 * \return decoded result, u
	 */
	Matrix<GF2> decode_v(const Matrix<my_double>& recv, int max_list_num, int order = 1) {

#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		int double_ope_num_before = my_double_auxiliary_storage::operation_number;
		int double_ope_num_after;
#endif // use_my_double

		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;

		int GF2e_ope_num_before = GF2e_auxiliary_storage::operation_number;
		int GF2e_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		// the parity Matrix should be sorted, and the function to generate the entire codeword should be added

		//cout << "r" << r;
		Matrix<my_double> r_abs = recv.get_abs();
		permuted_ind = r_abs.sort_with_ind('>');		// from big to small
		//cout << "recv" << recv;
		//cout << "r_abs" << r_abs;

		//cout << "permuted_ind" << permuted_ind;
		recv_permuted = recv;
		recv_permuted.permute(permuted_ind);
		//cout << "recv_permuted" << recv_permuted;

		permuted_ind_left = permuted_ind.get_part(0, 0, 0, k_prime - 1);
		left_sort_ind = permuted_ind_left.sort_with_ind('<');

		permuted_ind_right = permuted_ind.get_part(0, k_prime, 0, -1);
		//cout << "permuted_ind_right" << permuted_ind_right;

		right_sort_ind = permuted_ind_right.sort_with_ind('<');
		//cout << "right_sort_ind" << right_sort_ind;

		/*cout << "GF2e_auxiliary_storage::operation_number: before generate_systematic_generator_any_pos_best_r= " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/
			// rs generator matrix from most reliable information set of bch received codeword
		rs.generate_systematic_generator_any_pos_best(permuted_ind_left, true);

		/*cout << "GF2e_auxiliary_storage::operation_number: after generate_systematic_generator_any_pos_best_r= " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/

			//cout << "rs.generator_M_systematic_any_pos" << rs.generator_M_systematic_any_pos;
			//cout << "rs.information_set_ind" << rs.information_set_ind;

			// the non-systematic part of rs generator matrix, and transposed, n * r_prime matrix
		Matrix<GF2e<m>> parity_part = rs.generator_M_systematic_any_pos.get_cols(rs.redundancy_set_ind).Transpose();
		//cout << "parity_part" << parity_part;
		//cout << "left_sort_ind" << left_sort_ind;
		parity_part.permute_col_back(left_sort_ind);		// sort back according to recv_permuted
		//cout << "parity_part" << parity_part;

		// extend GF2e<m> to GF2 and apply to each element of parity_part, generate a (m*r_prime) * n matrix
		Matrix<GF2> extend_parity_part = GF2e<m>::to_bits_row_extention(parity_part);
		//cout << "extend_parity_part" << extend_parity_part;
		// throw away the rows having ending 1 at last r_prime columns
		Matrix<int> erase_row_ind(1, r_prime);
		for (int i = 0; i < r_prime; ++i) {
			erase_row_ind(i) = i * m;
		}
		//cout << "erase_row_ind" << erase_row_ind;
		sub_parity_matrix = extend_parity_part.erase_rows(erase_row_ind);
		//cout << "sub_parity_matrix" << sub_parity_matrix << endl;

		// randomly choose (k_prime - k) line
		Matrix<int> row_ind(1, sub_parity_matrix.row(), 'N');
		Matrix<int> chosen_row_ind = row_ind.get_random_element(selected_row_num);
		//cout << "chosen_row_ind" << chosen_row_ind;

		// by size of (k_prime - k) * k_prime
		Matrix<GF2> chosen_parity_check_matrix = sub_parity_matrix.get_rows(chosen_row_ind);
		//cout << "choosen_parity_check_matrix" << chosen_parity_check_matrix;

		vit.change_PM(chosen_parity_check_matrix);

		// select the valid codeword from list

		//Matrix<GF2> backup_not_processed_parity_check_matrix = sub_parity_matrix.erase_rows(chosen_row_ind);
		Matrix<GF2> not_processed_parity_check_matrix = sub_parity_matrix.erase_rows(chosen_row_ind);
		//cout << "not_processed_parity_check_matrix" << not_processed_parity_check_matrix;

		// compute the inner product of vector in viterbi_list and each row of not_processed_parity_check_matrix
		// if one result is non zero, then the codeword is invalid and should be discarded

		vit.change_unused_PM(not_processed_parity_check_matrix);

		extend_generator = extend_parity_part.get_rows(erase_row_ind).Transpose();
		//cout << "extend_generator" << extend_generator << endl;

		int invalid_record = 0;
		Matrix<GF2> v_hat = flip_to_find(order);
		if (v_hat.size() == 0) {
			// there is no valid list, decode fail
			//cout << "there is no valid list, i=" << i << endl;
			v_hat = Matrix<GF2>(1, n, '0');
			invalid_record++;
		}
		v_hat.permute_back(permuted_ind);
		if (is_ML) {
			num_early_stop++;
			return v_hat;
		}

		Matrix<GF2> v_hat2 = vit.decode_v_4_LL_OSD_Hybrid_once(recv_permuted, *this, max_list_num);

		if (v_hat2.size() == 0) {
			// there is no valid list, decode fail
			//cout << "there is no valid list, i=" << i << endl;
			v_hat2 = Matrix<GF2>(1, n, '0');
			invalid_record++;
		}
		v_hat2.permute_back(permuted_ind);
		num_invalid_list += invalid_record == 2;

		// hybrid them
		my_double measure_llosd = Measure::correlation_discrepancy_v(recv, v_hat);
		my_double measure_llosdv = Measure::correlation_discrepancy_v(recv, v_hat2);

#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(LL_OSD) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(LL_OSD) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << "(LL_OSD) GF2e_ope_num = " << GF2e_ope_num_after - GF2e_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return measure_llosd < measure_llosdv ? v_hat : v_hat2;
	}

	Matrix<GF2> flip_to_find(int order) {

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

		is_ML = false;
		Matrix<my_double> recv_abs = recv_permuted.get_abs();

		//cout << "r=" << r;
		Matrix<GF2> hdr = BPSK::demodulation(recv_permuted);
		//cout << "hdr=" << hdr;
		int error_num;

		//cout << "k_prime = " << k_prime << endl;
		Matrix<GF2> can_v_hat = hdr.get_part(0, 0, 0, k_prime - 1);
		Matrix<GF2> can_first_v_hat;
		Matrix<GF2> best_dv;
		Matrix<GF2> dv_order_0_left = can_v_hat;
		Matrix<GF2> dv_order_0_all = encode_by_right_systematic_PM(can_v_hat);
		//cout << "dv_order_0_all" << dv_order_0_all;
		my_double lambda_min = 3e10;			// a super large number to be updated

		Matrix<int> non_0_ind_aux(1, k_prime);
		total_used_list_num++;
		ML_used_list_num++;
		if (sub_parity_matrix.check_inner_product_4_GF2(can_v_hat, can_first_v_hat, non_0_ind_aux)) {
			can_first_v_hat = can_v_hat;
			//cout << "can_v_hat" << can_v_hat;
			best_dv = dv_order_0_all;
			//cout << "best_dv" << best_dv;

			// best_dv after permute back must in v space
			lambda_min = Measure::correlation_discrepancy_v(recv_abs, hdr, best_dv);
			error_num = Measure::error_num;
			if (optimum_condition_OSD::is_optimum(d, k_prime, best_dv, hdr, recv_abs, -1, lambda_min, error_num)) {
				//cout << "early stop, o=0" << endl;
				is_ML = true;
				return best_dv;
			}
		}

		// first not flip any position, for the first set we do demand in v space
		//cout << "best_dv: 0" << best_dv;

		int current_test_n = 0;
		Matrix<GF2> dv(1, k_prime);
		Matrix<GF2> dv_all(1, n);
		Matrix<int> real_flipped_mat(1, order, 'v');
		for (int o = 1; o <= order; ++o) {
			Matrix<int> flipped_mat(1, o, 'N'); // the order of flipped mat should be fixed
			do {
				dv = dv_order_0_left;
				real_flipped_mat.resize(1, 0, false);
				int flip_ind_size = flipped_mat.size();
				for (int i = 0; i < flip_ind_size; ++i) {
					real_flipped_mat.push_back(k_prime - 1 - flipped_mat(i));
				}
				for (int i = 0; i < flip_ind_size; ++i) {
					dv(real_flipped_mat(i)) += 1;		// add dv's information part out of ref_encode_by_right_systematic_PM
				}

				total_used_list_num++;
				ML_used_list_num++;
				if (sub_parity_matrix.check_inner_product_4_GF2(dv, can_first_v_hat, non_0_ind_aux)) {
					if (can_first_v_hat.size() == 0) {
						can_first_v_hat = dv;
					}

					dv_all = dv_order_0_all;
					ref_encode_by_right_systematic_PM(dv_all, real_flipped_mat);
					//cout << "dv_all" << dv_all << endl;

					my_double lambda = Measure::correlation_discrepancy_v(recv_abs, hdr, dv);
					error_num = Measure::error_num;
					if (lambda >= lambda_min);
					else {
						// update best_dv to this erasure decoding result
						best_dv = dv_all;
						lambda_min = lambda;
						//cout << "dv=" << dv; cout << "lambda_min=" << lambda << endl;
						if (optimum_condition_OSD::is_optimum(d, k_prime, best_dv, hdr, recv_abs, o - 1, lambda_min, error_num)) {
							//cout << "early stop, o=" << o << endl;
							is_ML = true;
							return best_dv;
						}
					}
				}
			} while (flip_TEP::next(flipped_mat, k_prime, o));
		}
		//cout << "best_dv" << best_dv;


#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(OSD_r) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(OSD_r) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return best_dv;
	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by hybrid of flipping all possible 1 bits in k_prime most reliable independent position, and serial list Viterbi algorithm
	 *
	 * \param recv: received codeword, after BPSK_demodulation
	 * \param rs: the mother code rs
	 * \param alpha: for stopping rule, smaller alpha, less constrain, less candidate codeword, earlier termination
	 * \param max_list_num: the max list number for list Verterbi decoding
	 * \return decoded result, u
	 */
	 //Matrix<GF2> decode_v(const Matrix<my_double>& recv, my_double beta, int max_list_num) {
	 //	Matrix<GF2> best_dv = llosd.decode_v(recv, 1);
	 //	if (llosd.is_ML) {
	 //		return best_dv;
	 //	}
	 //	Matrix<GF2> best_dv_vit = llosdv.decode_v(recv, beta, max_list_num);
	 //	// hybrid them
	 //	my_double measure_llosd = Measure::correlation_discrepancy_v(recv, best_dv);
	 //	my_double measure_llosdv = Measure::correlation_discrepancy_v(recv, best_dv_vit);
	 //	return measure_llosd < measure_llosdv ? best_dv : best_dv_vit;		// problematic
	 //}

	void set_sigma(my_double _sigma) {
		sigma = _sigma;
	}
};

template<class T, class LL_OSD_Hybrid_type>
Matrix<GF2> Viterbi_unordered_map::decode_v_4_LL_OSD_Hybrid_once(Matrix<T> recv, \
	LL_OSD_Hybrid_type& ll_osd_Hybrid_Viterbi, int max_list_num) {

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

	//cout << "r_L" << r_L;
	//cout << "n = " << n << endl;

	int k_prime = ll_osd_Hybrid_Viterbi.rs.get_k();
	Matrix<my_double> recv_k_prime_part = recv.get_part(0, 0, 0, k_prime - 1);		// r_L
	Matrix<my_double> recv_redundancy_part = recv.get_part(0, k_prime, 0, -1);		// r_R

	decode_v_once_inner(recv_k_prime_part, max_list_num);

	//cout << "r_part" << r_part;

	relative_metric.resize(1, max_list_num, false);

	// set for listing
	//Matrix<my_double> r_L = recv_k_prime_part.get_part(0, 0, 0, lc_osd_r.k + lc_osd_r.delta - 1);
	Matrix<GF2> hdr_R = BPSK::demodulation(recv_redundancy_part);

	// it seems problematic
	/*my_double gamma_Z = -Measure::Euclidean_distance(r_part, BPSK::modulation(hdr.get_part(0, 0, 0, lc_osd_r.k + lc_osd_r.delta - 1)));
	my_double gamma_ZL = -Measure::Euclidean_distance(recv_k_prime_part.get_part(0, lc_osd_r.k + lc_osd_r.delta, 0, -1), \
			BPSK::modulation(hdr.get_part(0, lc_osd_r.k + lc_osd_r.delta, 0, -1)));*/

	int r_size = recv_redundancy_part.size();
	Matrix<my_double> recv_redundancy_part_abs = recv_redundancy_part.get_abs();
	my_double r_abs_minus_1_squared = 0;
	my_double tmp_calculation;
	for (int i = 0; i < r_size; ++i) {
		tmp_calculation = recv_redundancy_part_abs(i) - 1;
		r_abs_minus_1_squared += tmp_calculation * tmp_calculation;
	}

	/*my_double n_plus_sum_r_squared_minus_2_r_abs = r_size;
	for (int i = 0; i < r_size; ++i) {
		n_plus_sum_r_squared_minus_2_r_abs += recv_redundancy_part(i) * recv_redundancy_part(i) - 2 * my::abs(recv_redundancy_part(i));
	}*/

	my_double d_ZR = 0;
	my_double scale_by_sigma = 2 / (ll_osd_Hybrid_Viterbi.sigma * ll_osd_Hybrid_Viterbi.sigma);
	for (int i = 0; i < r_size; ++i) {
		d_ZR += recv_redundancy_part_abs(i) / (1 + exp(scale_by_sigma * recv_redundancy_part_abs(i)));
	}
	d_ZR *= 4;
	d_ZR += r_abs_minus_1_squared;

	//my_double d_ZR = Measure::Euclidean_distance(recv_redundancy_part, hdr_R, r_abs_minus_1_squared);
	//my_double d_ZL = Measure::Euclidean_distance(r_L, BPSK::modulation(hdr.get_part(0, 0, 0, lc_osd_r.k + lc_osd_r.delta - 1)));

	//cout << "d_Z = " << d_Z << endl;
	//cout << "alpha_minus_1_gamma_ZL = " << alpha_minus_1_gamma_ZL << endl;

	// first not flip any position, for the first set we do demand in v space
	Matrix<GF2> best_dv;
	Matrix<GF2> first_dv;
	my_double lambda_min = invalid_large_metric;
	my_double d_copt = lambda_min;
	can_v_hat = list_v_hat.get_row(0);
	can_first_v_hat = can_v_hat;
	Matrix<GF2> ref_result = can_v_hat.multiply_transpose_of(unused_PM);		// this can be simplified, to be done
	int best_dv_list_num = max_list_num;
	//cout << "can_v_hat" << can_v_hat;
	//cout << "unused_PM" << unused_PM;

	can_v_hat.permute_back(PM_permutation);
	first_dv = ll_osd_Hybrid_Viterbi.encode_by_right_systematic_PM(can_v_hat);

	if (ref_result.isZero()) {
		// a valid codeword
		best_dv = first_dv;
		best_dv_list_num = 1;
		//best_dv = can_v_hat.insert_cols(can_v_hat * ll_osd_Hybrid_Viterbi.extend_generator, ll_osd_Hybrid_Viterbi.rs.redundancy_set_ind);
		lambda_min = absolute_metric + absolute_zero_col_mertric + \
			Measure::Euclidean_distance(recv_redundancy_part, best_dv.get_part(0, k_prime, 0, -1), r_abs_minus_1_squared);
		d_copt = lambda_min;
	}

	//cout << "best_dv: 0" << best_dv << endl;

	// best_dv after permute back must in v space
	//cout << "lambda_min = " << lambda_min << endl;

	//my_double gamma_C_opt = -lambda_min;
	//cout << "alpha_gamma_C_opt = " << alpha_gamma_C_opt << endl;

	/* assume that if PM is full rank, then there is no need to use unused_PM */

	Matrix<GF2> dv;
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
			my_double d_CL = relative_metric(used_list_num) + absolute_metric + absolute_zero_col_mertric;
			// this is ML		
			if (d_copt < d_CL + d_ZR) {
				/*cout << "break here, i=" << i << endl;
				cout << "d_copt = " << d_copt << endl;
				cout << "d_CL = " << d_CL << endl;*/

				ll_osd_Hybrid_Viterbi.total_used_list_num += (unsigned long long)used_list_num + 1;
				ll_osd_Hybrid_Viterbi.ML_used_list_num += best_dv_list_num;
				return best_dv;
			}

			can_v_hat = list_v_hat.get_row(used_list_num);

			if (unused_PM.check_inner_product_4_GF2(can_v_hat, can_first_v_hat, non_0_ind_aux, ref_result)) {

				can_v_hat.permute_back(PM_permutation);

				dv = first_dv;
				can_v_hat.diff_ind_inplace(dv, flip_ind);
				ll_osd_Hybrid_Viterbi.ref_encode_by_right_systematic_PM(dv, flip_ind);

				//dv = \
				// can_v_hat.insert_cols(can_v_hat * ll_osd_Hybrid_Viterbi.extend_generator, ll_osd_Hybrid_Viterbi.rs.redundancy_set_ind);

				//my_double gamma_Cr = -Measure::Euclidean_distance(r_part, BPSK::modulation(dv));
				//my_double lambda = Measure::Euclidean_distance(recv, BPSK::modulation(dv));
				my_double lambda = d_CL + Measure::Euclidean_distance(recv_redundancy_part, \
					dv.get_part(0, k_prime, 0, -1), r_abs_minus_1_squared);
				//my_double d_C = lambda;

				/*cout << "-------- i = " << i << "----------" << endl;
				cout << "d_C = " << d_C << endl;
				cout << "d_copt = " << d_copt << endl;
				cout << "d_C - d_Z = " << d_C - d_Z << " ?> alpha * (d_copt - d_Z) = " << alpha * (d_copt - d_Z) << endl;*/

				//if (gamma_Cr + gamma_ZL < alpha * (gamma_C_opt + gamma_Z)) {
				//	//cout << "break here" << endl;
				//	return best_dv;
				//}

				//my_double beta = 3;
				//if (d_C - d_Z > beta * (d_copt - d_Z)) {
				//	//cout << "break here, i = " << i << endl;
				//	ll_osd_Hybrid_Viterbi.total_used_list_num += i + 1;
				//	return best_dv;
				//}

				if (lambda >= lambda_min);
				else {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					d_copt = lambda_min;
					best_dv_list_num = used_list_num + 1;
				}

				/* assume that if PM is full rank, then there is no need to use unused_PM*/
			}
			used_list_num++;
		} while (used_list_num < max_list_num);
		//cout << "relative_metric" << relative_metric;
	}


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

	ll_osd_Hybrid_Viterbi.total_used_list_num += max_list_num;
	ll_osd_Hybrid_Viterbi.ML_used_list_num += best_dv_list_num;
	ll_osd_Hybrid_Viterbi.times_reach_max_list_num++;
	return best_dv;
}
