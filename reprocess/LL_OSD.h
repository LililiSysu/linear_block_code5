#pragma once
/*****************************************************************//**
 * \file   LL_OSD.h
 * \brief  LL_OSD by bit operation over PM
 *
 * \author lilili
 * \date   April 2023
 *********************************************************************/

#include"reprocess_common.h"
#include"../code/RS.h"
#include"../code/BCH.h"

// this is problematic, discard it
template<int m, int t, int k_prime>
class LL_OSD {
public:

	BCH<m, t> bch; 
	RS<m, k_prime> rs;

	Matrix<my_double> recv_permuted;
	Matrix<GF2> sub_parity_matrix;
	Matrix<GF2> extend_generator;

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
	unsigned long long total_used_list_num;		// counting used list num, function to be developed
	bool is_ML;

	LL_OSD() {
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
			int actual_col = k_prime - 1 - filp_ind(i);
			ref_codeword(actual_col) += 1;
			for (int j = 0; j < r_prime; ++j) {	// parity check bit, add the column of Patrity check matrix
				codeword_right(j) += extend_generator(actual_col, j);
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
	Matrix<GF2> decode_v(const Matrix<my_double>& recv, int order) {

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

		permuted_ind_right = permuted_ind.get_part(0, k_prime, 0, - 1);
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

		extend_generator = extend_parity_part.get_rows(erase_row_ind).Transpose();
		//cout << "extend_generator" << extend_generator << endl;

		Matrix<GF2> v_hat = flip_to_find(order);

		if (v_hat.size() == 0) {
			// there is no valid list, decode fail
			//cout << "there is no valid list, i=" << i << endl;
			v_hat = Matrix<GF2>(1, n, '0');
			num_invalid_list++;
		}

		v_hat.permute_back(permuted_ind);



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

		return v_hat;
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
				num_early_stop++;
				return best_dv;
			}
		}

		// first not flip any position, for the first set we do demand in v space
		//cout << "best_dv: 0" << best_dv;

		int current_test_n = 0;
		Matrix<GF2> dv(1, k_prime);
		Matrix<GF2> dv_all(1, n);
		for (int o = 1; o <= order; ++o) {
			Matrix<int> flipped_mat(1, o, 'N'); // the order of flipped mat should be fixed
			do {
				dv = dv_order_0_left;
				int flip_ind_size = flipped_mat.size();
				for (int i = 0; i < flip_ind_size; ++i) {
					dv(k_prime - 1 - flipped_mat(i)) += 1;		// add dv's information part out of ref_encode_by_right_systematic_PM
				}

				total_used_list_num++;
				if (sub_parity_matrix.check_inner_product_4_GF2(dv, can_first_v_hat, non_0_ind_aux)) {
					if (can_first_v_hat.size() == 0) {
						can_first_v_hat = dv;
					}

					dv_all = dv_order_0_all;
					ref_encode_by_right_systematic_PM(dv_all, flipped_mat);
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
							num_early_stop++;
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

};

class LLOSD_genius_aid{
private:
	int k;
	int k_prime;
	int first_order;
	int second_order;

public:
	LLOSD_genius_aid(int _k, int _k_prime, int _first_order, int _second_order): \
	k(_k),k_prime(_k_prime),first_order(_first_order), second_order(_second_order) {}

	/**
	 * .
	 * 
	 * \param v: the transmitted codeword
	 * \param r: the received sequence
	 * \return: if the error pattern falls into seg(first_order,second_order), return true, else return false
	 */
	bool is_decode_successful(Matrix<GF2> v, Matrix<my_double> r) {
		Matrix<my_double> r_abs = r.get_abs();
		Matrix<int> sort_permutation = r_abs.sort_with_ind('>');
		v.permute(sort_permutation);
		r.permute(sort_permutation);

		// conuting the error of k position
		int first_error = 0;
		for (int i = 0; i < k; ++i) {
			first_error += !((r(i) >= 0 && v(i) == 0) || (r(i) < 0 && v(i) == 1));
		}

		// conuting the error of k_prime position
		int second_error = 0;
		for (int i = k; i < k_prime; ++i) {
			second_error += !((r(i) >= 0 && v(i) == 0) || (r(i) < 0 && v(i) == 1));
		}

		if (first_error <= first_order && second_error <= second_order) {
			return true;
		}
		else {
			return false;
		}
	}
};
