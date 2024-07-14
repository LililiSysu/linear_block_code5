/*****************************************************************//**
 * \file   test_bin.h
 * \brief  test bin function to be discarded
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"test_common.h"

class test_bin {

public:

#ifdef need_trace_table
	static void generate_Generator_of_dual_RS_test() {
		const int m = 5;
		const int k_prime = 29;
		const int n = (1 << m) - 1;
		const int k_prime_dual = n - k_prime;
		const int d = n - k_prime + 1;
		const int t = d / 2;
		typedef GF2e<m> ty;
		ty::init();

		Matrix<int> total_ind(1, n, 'N');
		cout << "total_ind (int)" << total_ind;
		//Matrix<int> info_set_ind = total_ind.get_random_element(k_prime_dual);

		Matrix<int> info_set_ind(1, k_prime_dual, { 0,1 });
		cout << "info_set_ind (int)" << info_set_ind;

		RS<m, k_prime_dual> rs_4_dual;
		rs_4_dual.generate_systematic_generator_any_pos_dual(info_set_ind, false);
		cout << "rs_4_dual.generator_M_systematic_any_pos_dual (GF2e)" << rs_4_dual.generator_M_systematic_any_pos_dual;

		// verify this is the systematic generator of dual RS(n,k_prime), okay
		RS<m, k_prime> rs;
		cout << "rs.get_generator_matrix().multiply_transpose_of(rs_4_dual.generator_M_systematic_any_pos_dual) (GF2e)" << \
			rs.get_generator_matrix().multiply_transpose_of(rs_4_dual.generator_M_systematic_any_pos_dual);

		// use trace function to get sub-generator of (7,3) BCH dual
		Matrix<GF2> Q_sub_dual(k_prime_dual, n);
		for (int i = 0, imax = Q_sub_dual.size(); i < imax; ++i) {
			Q_sub_dual(i) = rs_4_dual.generator_M_systematic_any_pos_dual(i).trace();
		}
		cout << "Q_sub_dual (GF2)" << Q_sub_dual;

		// verify this is the sub-generator of (7,3) BCH dual
		BCH<m, t> bch;
		cout << "bch.get_generator_matrix().multiply_transpose_of(Q_sub_dual) (GF2)" << \
			bch.get_generator_matrix().multiply_transpose_of(Q_sub_dual);
		int k = bch.get_k();
		int k_dual = n - k;

		Matrix<ty> u_dual(1, k_prime_dual, '0');
		Matrix<ty> v_dual(1, n, '0');
		int k_gap = k_dual - k_prime_dual;
		Matrix<GF2> v_dual_trace(k_gap, n, '0');

		for (int kk = 0; kk < k_gap; ++kk) {
			for (int i = int(GF2e_auxiliary_storage::GF2e_that_trace_to_0.size() - 1); i >= 0; --i) {
				cout << endl << "-------- i = " << i << " -----------" << endl;

				u_dual(0) = GF2e_auxiliary_storage::GF2e_that_trace_to_0[i];
				v_dual = u_dual * rs_4_dual.generator_M_systematic_any_pos_dual;		// encode result
				cout << "v_dual (GF2e)" << v_dual;

				for (int j = 0, jmax = v_dual.size(); j < jmax; ++j) {
					v_dual_trace(j) = v_dual(j).trace();
				}
				cout << "v_dual_trace (GF2)" << v_dual_trace;

				if (v_dual_trace.isZero());
				else {
					break;
				}

				if (i == 0) {
					// we should count the probability, to verify
					cout << "generating v_dual_trace failed" << endl;
				}
			}
		}


		// we can adding the 'v_dual_trace' into the Q_sub_dual to form Q_dual
		Matrix<GF2> Q_dual(k_dual, n);		// to be fixed to k_dual
		Q_dual.set_part(0, 0, Q_sub_dual);
		Q_dual.set_part(k_prime_dual, 0, v_dual_trace);
		cout << "Q_dual (GF2)" << Q_dual;

		// verify this is the generator of (7,3) BCH dual
		cout << "bch.get_generator_matrix().multiply_transpose_of(Q_dual) (GF2)" << \
			bch.get_generator_matrix().multiply_transpose_of(Q_dual);
	}
	static void count_probability_of_generating_zero_vec() {
		const int m = 5;
		const int k_prime = 29;
		const int n = (1 << m) - 1;
		const int k_prime_dual = n - k_prime;
		const int d = n - k_prime + 1;
		const int t = d / 2;
		typedef GF2e<m> ty;
		ty::init();

		Matrix<int> total_ind(1, n, 'N');
		cout << "total_ind (int)" << total_ind;
		//Matrix<int> info_set_ind = total_ind.get_random_element(k_prime_dual);

		Matrix<int> info_set_ind_array = Matrix_common::generating_all_n_choose_k_pattern(n, k_prime_dual);
		Matrix<int> num_using_vectors(1, 10, '0');


		RS<m, k_prime_dual> rs_4_dual;
		RS<m, k_prime> rs;
		BCH<m, t> bch;
		int k = bch.get_k();
		int k_dual = n - k;

		for (int test_ind = 0, test_ind_max = info_set_ind_array.row(); test_ind < test_ind_max; ++test_ind) {

			cout << endl << "--------- test_ind = " << test_ind << "------------" << endl;

			Matrix<int> info_set_ind = info_set_ind_array.get_row(test_ind);
			cout << "info_set_ind (int)" << info_set_ind;

			rs_4_dual.generate_systematic_generator_any_pos_dual(info_set_ind, false);
			//cout << "rs_4_dual.generator_M_systematic_any_pos_dual (GF2e)" << rs_4_dual.generator_M_systematic_any_pos_dual;

			// verify this is the systematic generator of dual RS(n,k_prime), okay
			cout << "rs.get_generator_matrix().multiply_transpose_of(rs_4_dual.generator_M_systematic_any_pos_dual).isZero() = " << \
				rs.get_generator_matrix().multiply_transpose_of(rs_4_dual.generator_M_systematic_any_pos_dual).isZero() << endl;

			// use trace function to get sub-generator of (7,3) BCH dual
			Matrix<GF2> Q_sub_dual(k_prime_dual, n);
			for (int i = 0, imax = Q_sub_dual.size(); i < imax; ++i) {
				Q_sub_dual(i) = rs_4_dual.generator_M_systematic_any_pos_dual(i).trace();
			}
			//cout << "Q_sub_dual (GF2)" << Q_sub_dual;

			// verify this is the sub-generator of (7,3) BCH dual
			cout << "bch.get_generator_matrix().multiply_transpose_of(Q_sub_dual).isZero() = " << \
				bch.get_generator_matrix().multiply_transpose_of(Q_sub_dual).isZero() << endl;

			Matrix<ty> u_dual(1, k_prime_dual, '0');
			Matrix<ty> v_dual(1, n, '0');
			Matrix<GF2> v_dual_trace(1, n, '0');

			for (int i = int(GF2e_auxiliary_storage::GF2e_that_trace_to_0.size() - 1); i >= 0; --i) {
				//cout << endl << "-------- i = " << i << " -----------" << endl;

				u_dual(0) = GF2e_auxiliary_storage::GF2e_that_trace_to_0[i];
				//u_dual(1) = GF2e_auxiliary_storage::GF2e_that_trace_to_0.back();
				v_dual = u_dual * rs_4_dual.generator_M_systematic_any_pos_dual;		// encode result
				//cout << "v_dual (GF2e)" << v_dual;

				for (int i = 0, imax = v_dual.size(); i < imax; ++i) {
					v_dual_trace(i) = v_dual(i).trace();
				}
				//cout << "v_dual_trace (GF2)" << v_dual_trace;

				if (v_dual_trace.isZero());
				else {
					cout << "using vectors = " << GF2e_auxiliary_storage::GF2e_that_trace_to_0.size() - i << endl;
					num_using_vectors(int(GF2e_auxiliary_storage::GF2e_that_trace_to_0.size() - i))++;		// prevent warning
					break;		// we add just one more vector here
				}

				if (i == 0) {
					// we should count the probability, to verify
					cout << "generating v_dual_trace failed" << endl;
				}
			}

			// we can adding the 'v_dual_trace' into the Q_sub_dual to form Q_dual
			Matrix<GF2> Q_dual(k_prime_dual + 1, n);		// to be changed to k_dual
			Q_dual.set_part(0, 0, Q_sub_dual);
			Q_dual.set_part(k_prime_dual, 0, v_dual_trace);
			//cout << "Q_dual (GF2)" << Q_dual;

			// verify this is the generator of (7,3) BCH dual
			cout << "bch.get_generator_matrix().multiply_transpose_of(Q_dual).isZero() = " << \
				bch.get_generator_matrix().multiply_transpose_of(Q_dual).isZero() << endl;
		}

		cout << "------------------------" << endl;
		cout << "num_using_vectors" << num_using_vectors;
		cout << "------------------------" << endl;
	}
	static void trace_generate_linear_combination_test() {
		const int m = 5;
		typedef GF2e<m> ty;
		ty::init();

		for (int i = 0; i < GF2e_auxiliary_storage::q; ++i) {
			cout << i << ": " << GF2e_auxiliary_storage::trace_table[i] << endl;
		}
	}
#endif

	static void new_generation_of_G_residual_test_2() {
		const int m = 5;
		const int k_prime = 25;
		const int n = (1 << m) - 1;
		const int k_prime_dual = n - k_prime;
		const int d = n - k_prime + 1;
		const int t = d / 2;
		typedef GF2e<m> ty;
		ty::init();

		//Matrix<int> total_ind(1, n, 'N');
		//cout << "total_ind (int)" << total_ind;
		//Matrix<int> info_set_ind = total_ind.get_random_element(k_prime_dual);

		//Matrix<int> info_set_ind_array = Matrix_common::generating_all_n_choose_k_pattern(n, k_prime_dual);
		Matrix<int> info_set_ind_array(1, k_prime_dual, 'N');
		Matrix<int> num_using_vectors(1, 10, '0');

		RS<m, k_prime_dual> rs_4_dual;
		RS<m, k_prime> rs;
		BCH<m, t> bch;
		int k = bch.get_k();
		int k_dual = n - k;

		int external_considered_info_bit = (k_dual - k_prime_dual) / 2;
		//int more_considered_info_bit = n - k_dual;

		Matrix<ty> testing_coordinate(1, m - 1, { 2, 4, 9, 16 });	// these corodinate generates all zero bits in k_prime_dual info bits

		subfield_subcode_sys_partial_generator<m> spg(n, k_prime_dual, k_dual, external_considered_info_bit, testing_coordinate);

		for (int test_ind = 0, test_ind_max = info_set_ind_array.row(); test_ind < test_ind_max; ++test_ind) {

			cout << endl << "--------- test_ind = " << test_ind << "------------" << endl;

			Matrix<int> info_set_ind = info_set_ind_array.get_row(test_ind);
			cout << "info_set_ind (int)" << info_set_ind;

			rs_4_dual.generate_systematic_generator_any_pos_dual(info_set_ind, false);
			//cout << "rs_4_dual.generator_M_systematic_any_pos_dual (GF2e)" << rs_4_dual.generator_M_systematic_any_pos_dual;

			// verify this is the systematic generator of dual RS(n,k_prime), okay
			cout << "rs.get_generator_matrix().multiply_transpose_of(rs_4_dual.generator_M_systematic_any_pos_dual).isZero() = " << \
				rs.get_generator_matrix().multiply_transpose_of(rs_4_dual.generator_M_systematic_any_pos_dual).isZero() << endl;

			// use trace function to get sub-generator of (7,3) BCH dual
			Matrix<GF2> Q_sub_dual(k_prime_dual, n);
			for (int i = 0, imax = Q_sub_dual.size(); i < imax; ++i) {
				Q_sub_dual(i) = rs_4_dual.generator_M_systematic_any_pos_dual(i).trace();
				//that is encoding of (1,0,...,0), ... ,(0,0,...,1)
			}
			//cout << "Q_sub_dual (GF2)" << Q_sub_dual;

			// verify this is the sub-generator of (7,3) BCH dual
			cout << "bch.get_generator_matrix().multiply_transpose_of(Q_sub_dual).isZero() = " << \
				bch.get_generator_matrix().multiply_transpose_of(Q_sub_dual).isZero() << endl;

			spg.init(rs_4_dual.generator_M_systematic_any_pos_dual);

			spg.get_generator();

			cout << "Q_sub_dual" << Q_sub_dual;

			cout << "spg.generator_info_part" << spg.generator_info_part;
			cout << "spg.encode_index_and_xor_index" << spg.xor_index;
			cout << "spg.generator_ans" << spg.generator_ans;

			// we can adding the 'v_dual_trace' into the Q_sub_dual to form Q_dual
			Matrix<GF2> Q_dual(k_dual, n, '0');		// to be changed to k_dual
			Q_dual.set_part(0, 0, Q_sub_dual);
			Q_dual.set_part(k_prime_dual, 0, spg.generator_ans);
			cout << "Q_dual (GF2)" << Q_dual;

			// verify this is the generator of (7,3) BCH dual
			cout << "bch.get_generator_matrix().multiply_transpose_of(Q_dual).isZero() = " << \
				bch.get_generator_matrix().multiply_transpose_of(Q_dual).isZero() << endl;
		}
	}
	static void GE_left_identity_4_GF2_with_last_w1_col_test() {
		// never use it, having problem
		const int m = 4;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 2;
		eBCH<m, t> ebch;
		ebch.print_info();

		Matrix<GF2> G = ebch.get_generator_matrix();
		cout << "G (orig)" << G;

		Matrix<int> permute_record = G.GE_left_identity_4_GF2_with_last_w1_col();
		cout << "G (left iden)" << G;

		G.permute_col_back(permute_record);
		cout << "G (col back)" << G;

		Matrix<GF2> H = ebch.get_parity_matrix();
		cout << "G.multiply_transpose_of(H).isZero() = " << G.multiply_transpose_of(H).isZero() << endl;

		// check out the shift property of G
		int n = G.col();
		Matrix<int> natural(1, n - 1, 'N');
		cout << "natural (orig)" << natural;

		// the first n-1 columns of G can be shift, not violating that G is a generator matrix
		Matrix<int> natural_shift;
		natural.col_shift_right_cir(2, natural_shift);
		cout << "natural_shift" << natural_shift;

		G.permute_col(natural_shift);
		cout << "G (partly shift)" << G;
		cout << "G.multiply_transpose_of(H).isZero() = " << G.multiply_transpose_of(H).isZero() << endl;		// right

	}
};
