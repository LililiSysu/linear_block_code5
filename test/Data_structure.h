/*****************************************************************//**
 * \file   Data_structure.h
 * \brief  test data structure in this file
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"test_common.h"

struct no_cmp_tst {
	char c;
};

class test_Data_sturcture {
public:

	static void heap_test() {
		Heap_max<int> A;
		A.resize(1, 20, false);
		for (int i = 0; i < 20; ++i) {
			A[i] = my::rand_int_adv(-8, 34);
		}
		cout << "A" << A;
		A.build();
		cout << "A" << A;
		cout << "A.is_heap() = " << A.is_heap() << endl;
		A.push(12);
		A.push(10);
		A.push(15);
		cout << "A: before" << A;
		A.push(4);
		cout << "A: after" << A;
		A.push(-3);
		A.push(-2);
		A.push(0);
		cout << "A" << A;
		cout << "A.is_heap() = " << A.is_heap() << endl;

		for (int i = 0; i < 6; ++i) {
			int A_bottom_ind = A.bottom();
			cout << "A bottom = " << A(A_bottom_ind) << endl;
			A.pop(A_bottom_ind);
		}

		while (A.size() != 0) {
			cout << "A top = " << A(0) << endl;
			A.pop();
		}

		cout << "A.size() = " << A.size() << endl;
	}
	static void vector_test() {
		vector<int> v = { 53 };
		for (int i = 0, vs = (int)v.size(); i < vs; ++i) {
			cout << v[i] << ", ";
		}
		cout << endl;
	}
	static void flip_TEP_diff_test() {
		Matrix<int> last_flipped_mat;
		int o = 3;
		int k = 15;
		Matrix<int> flipped_mat(1, o, 'N');
		Matrix<int> fmd(1, 2 * o);
		int num = 0;
		do {
			cout << "-------------" << num << "------------" << endl;
			flip_TEP::diff(last_flipped_mat, flipped_mat, fmd);
			last_flipped_mat = flipped_mat;
			cout << "flipped_mat" << flipped_mat;
			cout << "fmd" << fmd;

			num++;
		} while (flip_TEP::next(flipped_mat, k, o));
	}
	static void data_ope_test() {
		/*const int GF2_len = 5;
		const int int_size = 1 << GF2_len;
		for (int i = 0; i < int_size; ++i) {
			Matrix<GF2> v = data_ope::i2GF2_natual(i, GF2_len);
			cout << v;
		}*/

		// assume the code of eBCH(128,64,22)
		// delta = 8, Lm = 8192
		// the length of right part is 128-64-8=56
		// assume Pe_ML = 7.4e-3 @ SNR of 2 dB

		int len = 56;
		my_double Pe_ML = 7.43e-3;
		Matrix<GF2> p = flip_TEP::get_estimated_flip_pattern(len, Pe_ML);
		cout << "p" << p;

		/*int len = 5;
		int pow2_len = 1 << len;
		for (int i = 0; i < pow2_len; ++i) {
			Matrix<GF2> flip_pattern = data_ope::i2GF2_osd(i, len);
			cout << "flip_pattern" << flip_pattern;
		}*/
	}
	static void choose_order_test() {
		Matrix<my_double> y(1, 4, { 0.27, 0.19, 0.12, 0.06 });
		Matrix<GF2> OSD_order(16, 4,
			{
				0,0,0,0,
				0,0,0,1,
				0,0,1,0,
				0,1,0,0,
				1,0,0,0,
				0,0,1,1,
				0,1,0,1,
				1,0,0,1,
				0,1,1,0,
				1,0,1,0,
				1,1,0,0,
				0,1,1,1,
				1,0,1,1,
				1,1,0,1,
				1,1,1,0,
				1,1,1,1
			});

		Matrix<GF2> natual_order(16, 4,
			{
				0,0,0,0,
				0,0,0,1,
				0,0,1,0,
				0,0,1,1,
				0,1,0,0,
				0,1,0,1,
				0,1,1,0,
				0,1,1,1,
				1,0,0,0,
				1,0,0,1,
				1,0,1,0,
				1,0,1,1,
				1,1,0,0,
				1,1,0,1,
				1,1,1,0,
				1,1,1,1,
			});

		// compute the Euclidean distance of each ordering
		Matrix<my_double> OSD_order_ds(1, 16, '0');
		Matrix<my_double> natual_order_ds(1, 16, '0');
		for (int i = 0; i < 16; ++i) {
			for (int j = 0; j < 4; ++j) {
				my_double tmp = my_double(BPSK::modulation(OSD_order(i, j))) - y(j);
				OSD_order_ds(i) += tmp * tmp;

				tmp = my_double(BPSK::modulation(natual_order(i, j))) - y(j);
				natual_order_ds(i) += tmp * tmp;
			}
		}
		cout << "OSD_order_ds" << OSD_order_ds;		// it will be better using osd order
		cout << "natual_order_ds" << natual_order_ds;
	}
	static void Euclidean_distance_test() {
		const int n = 30;
		Matrix<my_double> r(1, n);
		Matrix<GF2> v(1, n);
		for (int i = 0; i < n; ++i) {
			r(i) = my::rand_u_adv();
			v(i) = my::rand_01_adv();
		}
		cout << "r" << r;
		cout << "v" << v;


		cout << "Measure::correlation_discrepancy_v(r, v) = " << Measure::correlation_discrepancy_v(r, v) << endl;
		cout << "error_num = " << Measure::error_num << endl;
		cout << "Measure::Euclidean_distance(r, BPSK::modulation(v)) = " << Measure::Euclidean_distance(r, BPSK::modulation(v)) << endl;
		cout << "error_num = " << Measure::error_num << endl;
		int r_size = r.size();
		my_double n_plus_sum_r_squared_minus_2_r_abs = r_size;
		for (int i = 0; i < r_size; ++i) {
			n_plus_sum_r_squared_minus_2_r_abs += r(i) * r(i) - 2 * my::abs(r(i));
		}
		cout << "Measure::Euclidean_distance(r, v, sum_r_squared_plus_n_minus_2_sum_r_abs) = " << \
			Measure::Euclidean_distance(r, v, n_plus_sum_r_squared_minus_2_r_abs) << endl;
		cout << "error_num = " << Measure::error_num << endl;


	}
	static void P2_matrix_test() {
		Matrix<GF2> P2(30, 57, {
			1,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,0,0,0,0,1,1,0,1,0,0,1,1,0,1,1,1,1,1,0,0,0,0,0,1,1,1,0,1,1,0,1,
			1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,0,1,0,0,1,1,1,1,1,1,0,1,1,1,0,1,0,0,1,0,0,0,0,0,1,1,1,0,1,0,1,0,0,0,0,0,
			1,1,0,1,1,1,1,1,1,1,0,1,0,0,1,1,0,0,0,1,0,1,0,1,1,1,1,0,0,1,0,1,1,0,1,0,1,1,1,0,0,1,0,1,0,0,1,1,0,1,0,1,1,1,1,1,1,
			0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,1,0,0,0,1,0,1,1,0,1,1,0,0,0,0,0,1,
			0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,0,0,0,1,1,1,0,0,1,0,0,1,1,0,1,0,1,1,1,1,1,1,0,1,0,1,0,1,0,0,1,1,1,1,1,1,1,1,0,1,0,1,
			1,1,1,0,1,1,1,1,0,0,0,0,1,0,1,1,1,1,0,0,1,1,0,1,1,0,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,0,1,0,0,1,1,0,1,1,0,1,0,1,1,0,
			0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,1,0,1,1,0,0,0,1,1,0,1,0,0,0,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,1,0,
			0,1,1,0,0,0,1,1,0,1,0,0,0,0,1,1,0,0,1,0,0,1,0,0,1,1,1,0,0,0,1,0,1,0,1,0,0,0,1,1,1,0,1,1,0,0,0,0,1,0,1,0,1,0,0,0,1,
			0,1,1,1,0,1,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,1,1,1,0,0,0,0,0,1,0,1,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,1,1,1,1,0,
			0,1,0,1,1,0,1,1,1,0,0,1,1,0,1,1,0,1,0,0,1,1,1,0,1,0,0,1,1,1,0,0,1,0,0,0,0,1,0,1,0,0,1,0,1,1,1,0,1,1,1,0,1,1,0,1,1,
			0,1,0,1,0,0,0,1,0,1,1,0,0,0,1,0,1,1,0,1,1,1,1,1,1,0,1,1,0,1,0,1,0,0,0,1,0,0,0,1,1,1,1,0,1,1,0,1,0,0,0,1,0,0,0,1,0,
			1,1,1,0,0,0,0,0,0,1,0,0,1,1,0,1,1,1,1,0,1,1,1,1,0,0,1,1,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,1,1,0,1,0,1,
			0,1,1,1,0,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,0,1,1,1,0,0,0,1,0,1,0,0,0,1,1,1,0,1,0,0,1,0,1,1,1,1,1,1,1,1,1,1,0,
			1,1,0,0,1,0,1,1,1,1,1,0,0,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,0,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,1,0,1,1,
			1,1,1,1,1,0,1,0,0,1,1,1,1,1,1,0,1,1,1,1,0,1,0,1,0,0,1,0,0,0,0,1,1,1,1,1,1,0,0,0,1,1,1,1,0,1,0,1,1,0,1,1,0,0,0,1,1,
			0,0,1,0,0,1,1,1,0,1,0,1,0,1,0,0,1,1,0,0,0,0,1,1,1,1,1,1,0,1,0,0,0,1,1,0,1,1,1,1,1,0,1,1,1,0,1,1,0,1,1,0,0,0,1,1,0,
			0,0,1,1,0,0,0,1,1,0,0,0,1,1,0,1,0,1,0,1,0,0,1,1,1,1,0,0,1,1,1,0,0,1,1,1,1,0,1,1,1,1,1,0,1,1,0,1,0,0,0,1,1,0,1,1,0,
			1,0,0,1,0,1,1,1,1,1,0,0,1,1,1,0,1,1,0,1,1,1,0,1,0,0,1,0,0,1,0,0,1,0,1,1,0,0,0,1,1,0,0,1,0,1,1,0,1,1,0,0,0,0,0,0,0,
			0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,1,1,0,0,1,0,0,0,1,1,0,0,0,0,0,0,
			0,1,1,0,1,0,1,0,1,1,1,0,0,0,0,1,1,0,0,1,0,0,1,1,1,0,1,0,0,1,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,0,1,1,1,1,1,0,1,
			0,1,1,1,1,1,1,0,1,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,1,0,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,
			0,1,1,0,0,1,1,0,0,1,1,1,0,1,1,0,0,0,0,0,1,0,1,0,1,1,0,1,0,1,1,1,1,1,0,1,1,0,0,0,0,0,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,
			1,0,0,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,0,1,1,0,1,0,1,1,1,0,1,0,1,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,0,0,
			1,0,0,1,0,1,0,0,1,0,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,0,0,1,0,1,0,0,1,1,0,1,1,0,1,0,
			1,0,0,1,0,0,0,0,0,1,0,1,1,0,1,1,1,1,1,0,1,1,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,0,1,0,0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,
			0,0,1,1,0,0,1,0,1,0,1,0,1,0,1,1,0,1,1,0,0,1,1,0,1,1,1,1,1,0,1,1,1,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,0,0,
			1,0,1,0,1,1,1,1,1,0,1,0,1,1,1,1,0,0,1,0,1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0,1,0,1,0,
			1,1,0,1,1,0,0,1,0,0,0,1,0,0,1,1,1,0,1,0,1,1,1,1,1,1,0,1,0,0,1,0,0,1,0,1,0,1,0,1,0,0,1,1,0,1,0,0,0,0,1,1,0,0,0,0,1,
			0,0,1,1,0,1,0,0,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,0,1,1,1,0,0,1,0,1,1,1,1,1,1,0,1,1,0,0,
			1,1,0,1,1,1,0,0,0,1,0,0,1,1,0,0,1,0,1,1,1,0,0,1,1,0,1,1,0,1,0,1,1,1,0,1,1,0,1,1,1,0,0,0,1,0,0,1,0,0,1,1,1,0,1,0,1,
			});

		cout << "P2" << P2;

		// just for playing, to see how to choose P2
		Matrix<GF2> sub_parity_matrix_copy = P2;		// the original P2

		sub_parity_matrix_copy.row_transformation_to_low_triangle();
		cout << "sub_parity_matrix_copy" << sub_parity_matrix_copy;

		sub_parity_matrix_copy = sub_parity_matrix_copy.get_part(18, 0, -1, -1);


		Matrix<int> PM_permutation_test = sub_parity_matrix_copy.col_permute_to_full_rank_on_right();
		cout << "sub_parity_matrix_copy" << sub_parity_matrix_copy;

		sub_parity_matrix_copy.row_transformation_right_low_triangle_to_identity();
		cout << "sub_parity_matrix_copy" << sub_parity_matrix_copy;

		sub_parity_matrix_copy.row_transformation_to_up_triangle();
		cout << "sub_parity_matrix_copy" << sub_parity_matrix_copy;

		PM_permutation_test.permute(find_opt_PM::column_permute(sub_parity_matrix_copy));
		cout << "sub_parity_matrix_copy" << sub_parity_matrix_copy;
		cout << "PM_permutation_test" << PM_permutation_test;

		// no solution
	}
	static void RS_Generator_Lagrange_interpolation_test() {
		const int m = 5;
		const int k = 27;
		typedef GF2e<m> ty;
		ty::init();
		cout << "GF2e_auxiliary_storage::operation_number: = " << GF2e_auxiliary_storage::operation_number << endl;

		RS<m, k> rs;
		rs.print_info();
		int n = rs.get_n();
		int d = rs.get_d();

		Matrix<int> info_choose(1, (1 << m) - 1, 'N');
		Matrix<int> info_set = info_choose.get_random_element(k);
		//Matrix<int> info_set(1, 7, { 2,1,4,14,10,9,13 });
		//Matrix<int> info_set(1, 5, { 2,1,4,14,10 });
		//info_set.sort();
		cout << "info_set" << info_set;

		unsigned long long before, after, before_add, before_mul, after_add, after_mul;

		before = GF2e_auxiliary_storage::operation_number;
		before_add = GF2e_auxiliary_storage::add_number;
		before_mul = GF2e_auxiliary_storage::mul_number;

		rs.generate_systematic_generator_any_pos_best(info_set);
		cout << "rs" << rs.generator_M_systematic_any_pos << endl;

		after = GF2e_auxiliary_storage::operation_number;
		after_add = GF2e_auxiliary_storage::add_number;
		after_mul = GF2e_auxiliary_storage::mul_number;


		cout << "computation cost after polynomial generation 3 = " << after - before << endl;
		cout << "add cost after polynomial generation 3 = " << after_add - before_add << endl;
		cout << "mul cost after polynomial generation 3 = " << after_mul - before_mul << endl;
	}
	static void matrix_mul_transpose_test() {
		Matrix<int> m(4, 10, my::rand_int);
		Matrix<int> M(6, 10, my::rand_int);
		Matrix<int> real_result = m * M.Transpose();
		Matrix<int> test_result = m.multiply_transpose_of(M);
		cout << "real_result" << real_result;
		cout << "test_result" << test_result;

		cout << "should be zero" << real_result - test_result;
	}
	static void my_compute_test() {
		/*cout << "53 order 6 = " << my::n_choose_k(53, 6) + my::n_choose_k(53, 5) + my::n_choose_k(53, 4) + \
			my::n_choose_k(53, 3) + my::n_choose_k(53, 2) + my::n_choose_k(53, 1) + 1 << endl;*/

			/*cout << "Seg(2,4) = " << (my::n_choose_k_long(36, 0) + my::n_choose_k_long(36, 1) + my::n_choose_k_long(36, 2)) * \
				(my::n_choose_k_long(17, 0) + my::n_choose_k_long(17, 1) + my::n_choose_k_long(17, 2) + \
					my::n_choose_k_long(17, 3) + my::n_choose_k_long(17, 4)) << endl;*/

		int k = 43;
		int k_prime = 99;
		int order_first = 2;
		int order_second = 8;

		int dk = k_prime - k;
		unsigned long long TEP_num_first = 0;
		for (int i = 0; i <= order_first; ++i) {
			TEP_num_first += my::n_choose_k_long(k, i);
		}

		unsigned long long TEP_num_second = 0;
		for (int i = 0; i <= order_second; ++i) {
			TEP_num_second += my::n_choose_k_long(dk, i);
		}

		cout << "TEP_num_first = " << TEP_num_first << endl;
		cout << "TEP_num_second = " << TEP_num_second << endl;

		unsigned long long TEP_num_total = TEP_num_first * TEP_num_second;
		cout << "TEP_num_total = " << TEP_num_total << endl;
	}
	static void BCH_generate_test() {
		const int m = 6;
		const int t = 4;

		ope_count start, end;

		start.count();
		GF2e<m>::init();
		end.count();
		cout << "GF2e init cost\n\t" << end - start;

		cout << "---------" << endl;

		start.count();
		BCH<m, t> bch;			// for eBCH code initialization, a step of matrix inversion cause the large amout of operation
		end.count();
		cout << "BCH init cost\n\t" << end - start;

		cout << "---------" << endl;
		bch.print_info();

		//cout << "bch.get_parity_matrix()" << bch.get_parity_matrix();
	}
	static void nnsBCH_generate_test() {
		const int m = 6;
		const int t = 3;

		ope_count start, end;

		start.count();
		GF2e<m>::init();
		end.count();
		cout << "GF2e init cost\n\t" << end - start;

		cout << "---------" << endl;

		start.count();
		nnsBCH<m, t> nnsbch;			// for eBCH code initialization, a step of matrix inversion cause the large amout of operation
		end.count();
		cout << "nnsBCH init cost\n\t" << end - start;

		cout << "---------" << endl;
		nnsbch.print_info();

		//cout << "nnsbch.get_generator_matrix()" << nnsbch.get_generator_matrix();
		//cout << "nnsbch.get_parity_matrix()" << nnsbch.get_parity_matrix();
		//cout << "zero test" << nnsbch.get_parity_matrix().multiply_transpose_of(nnsbch.get_generator_matrix());
	}
	static void npBCH_generate_test() {
		const int m = 8;
		const int t = 2;
		const int b = 5;		// performance is very worse

		ope_count start, end;

		start.count();
		GF2e<m>::init();
		end.count();
		cout << "GF2e init cost\n\t" << end - start;

		cout << "---------" << endl;

		start.count();
		npBCH<m, t, b> npbch;			// for eBCH code initialization, a step of matrix inversion cause the large amout of operation
		end.count();
		cout << "npBCH init cost\n\t" << end - start;

		cout << "---------" << endl;
		npbch.print_info();

		/*cout << "npbch.get_generator_matrix()" << npbch.get_generator_matrix();
		cout << "npbch.get_parity_matrix()" << npbch.get_parity_matrix();
		cout << "zero test" << npbch.get_parity_matrix().multiply_transpose_of(npbch.get_generator_matrix());*/
	}

	static void my_vector_test() {
		const int cap = 1000;
		for (int i = 0; i < 100000; ++i) {		// no memory leakage 
			Vec_s<double, cap> mv;
			mv.push_back(1);
			mv.push_back(my::pi);
			mv.push_back(exp(1));
			mv.push_back(2 / float(9));

			cout << "----i  = " << i << "----" << endl;
			cout << "mv" << mv;
		}

	}
	static void dList_test() {
		dList<int> A;
		A.resize(1, 10);
		A.set_natural();
		A.build();

		cout << "A" << A;

		A.move_start();
		A.remove();

		cout << "A" << A;

		A.move_start();
		A.remove();

		cout << "A" << A;

		A.move_start();

		A.move_next();
		A.move_next();

		A.remove();

		cout << "A" << A;

		A.move_start();

		A.move_next();
		A.move_next();

		A.remove();

		cout << "A" << A;

		A.move_start();
		A.remove();

		cout << "A" << A;

		A.move_start();
		A.move_next();
		A.move_next();
		A.move_next();
		A.move_next();
		A.remove();

		cout << "A" << A;			// good right now

		A.move_start();
		A.move_next();
		A.remove();

		cout << "A" << A;			// good right now
	}
	static void for_loop_test() {
		int j = 3;
		j += 5;
		for (int i = 0; i < 10; ++i, j = 2 * i, j--) {
			// verify that comma operator in for loop is done by sequential
			cout << "i=" << i << ", " << "j=" << j << endl;
		}
	}

	static void triple_operator_test() {
		int tmp[3] = { 1,2,3 };
		//vector<int> tmp = { 1,2,3 };
		int ind = -10;
		int k = (ind >= 0 && ind <= 2) ? tmp[ind] : -1;         // seems alright, but acutally in error
		cout << "k = " << k << endl;
		int o = tmp[ind];                                       // invalid index retrive
		cout << "o = " << o << endl;
	}
	static void Vec_s_test() {
		no_cmp_tst dd;

		Vec_s<no_cmp_tst, 3> w;     // actually okay with that
		dd.c = '4';
		w.push_back(dd);
		dd.c = 'e';
		w.push_back(dd);
		dd.c = 'f';
		w.push_back(dd);
		//cout << "w" << w; // if we donot call <<, it will be fine for class Vec_sc for type even with out operator << defined

		//w.sort();         // if we donot call sort function, it will be fine to use class Vec_sc for type even without comparator defined
	   // cout << "w" << w;


		typedef Complex ty;

		Vec_s<ty, 10> v(10);
		for (int i = 0; i < 10; ++i) {
			v[i] = i;
		}
		cout << "v (orig)" << v;

		v.permute_rand();
		cout << "v (permute_rand)" << v;

		Vec_s<int, 10> v_ind;
		v.sort_with_ind(v_ind);

		cout << "v (sort)" << v;
		cout << "v_ind (sort_ind)" << v_ind;

		// test copy constructor
		v = vector<ty>({ 3, 2, 5, 2, 6 });          // a good way to set by hand
		cout << "v (Vector_static)" << v;

		Vec_s<ty, 10> p;
		cout << "p (default initializer)" << p;

		Vec_s<ty, 10> g(6);
		cout << "g (initialize size)" << g;

		v.shift_right_cir(3);
		cout << "v" << v;
	}
	static void Heap_max_s_test() {
		// can also test Heap_min

		// test heap function
		const int cap = 10;
		Heap_max_s<double, 10> hms({ 1,5,4,33,6,1 });
		cout << "hms" << hms;
		hms.build();
		cout << "hms" << hms;

		while (hms.size() != 0) {
			double top = hms.top();
			cout << "top = " << top << endl;
			hms.pop();      // the pop function is for heap
		}

		// test heap initializer
		const int sss = 6;
		Heap_min_s<Complex, sss> hc;
		hc = Vec_s<Complex, sss>(0);              // this is allowed for child of child initialization
		for (int i = 0; i < sss; ++i) {
			hc.push_back(Complex(my::rand_u(), my::rand_ga()));
		}
		cout << "hc" << hc;
		hc.build();
		cout << "hc" << hc;

		while (hc.size() != 0) {
			Complex top = hc.top();
			cout << "top = " << top << endl;
			hc.pop();      // the pop function is for heap
		}
	}
	static void sList_s_test() {
		const int cap = 10;
		sList_s<int, cap> A(cap);
		for (int i = 0; i < cap; ++i) {
			A[i] = i;
		}

		A.build();

		cout << "A" << A;

		A.remove_start();

		cout << "A" << A;

		A.remove_start();

		cout << "A" << A;

		A.move_start();
		A.move_next();
		A.remove_next();

		cout << "A" << A;

		A.move_start();
		A.move_next();
		A.remove_next();

		cout << "A" << A;

		A.remove_start();

		cout << "A" << A;

		A.move_start();
		A.move_next();
		A.move_next();
		A.move_next();
		A.remove_next();

		cout << "A" << A;			// good right now

		A.move_start();
		A.remove_next();

		cout << "A" << A;			// good right now
	}
	static void dList_s_test() {
		const int cap = 10;
		dList_s<int, cap> A(cap);
		for (int i = 0; i < cap; ++i) {
			A[i] = i;
		}

		A.build();

		cout << "A" << A;

		A.move_start();
		A.remove();

		cout << "A" << A;

		A.move_start();
		A.remove();

		cout << "A" << A;

		A.move_start();

		A.move_next();
		A.move_next();

		A.remove();

		cout << "A" << A;

		A.move_start();

		A.move_next();
		A.move_next();

		A.remove();

		cout << "A" << A;

		A.move_start();
		A.remove();

		cout << "A" << A;

		A.move_start();
		A.move_next();
		A.move_next();
		A.move_next();
		A.move_next();
		A.remove();

		cout << "A" << A;			// good right now

		A.move_start();
		A.move_next();
		A.remove();

		cout << "A" << A;			// good right now
	}
	static void dList_s_test2() {
		const int cap = 11;
		const int half_cap = cap / 2;
		dList_s<int, cap> A(half_cap);
		for (int i = 0; i < half_cap; ++i) {
			A[i] = i;
		}

		A.build();

		cout << "A" << A;

		A.move_start();
		A.remove();

		cout << "A" << A;

		A.move_next();
		A.insert_next(6);

		cout << "A" << A;

		cout << "A.current_val() = " << A.current_val() << endl;

		A.insert_former(34);

		cout << "A" << A;

		A.move_start();
		A.insert_former(-19);
		A.insert_former(-25);
		A.insert_former(-17);

		cout << "A" << A;

		A.move_end();
		A.remove();
		A.move_end();
		A.remove();

		cout << "A" << A;

		A.move_end();
		A.insert_former(40);
		A.insert_next(21);

		cout << "A" << A;
		cout << "A.size() = " << A.size() << endl;

		A.insert_next(15);
		bool is_success = A.insert_former(9);

		cout << "A" << A;
		cout << "is_success = " << is_success << endl;

		cout << "A.size() = " << A.size() << endl;

		A.clear();
		cout << "A" << A;

		for (int i = 0; i < 5; ++i) {
			A.insert_next(i);
		}
		cout << "A" << A;

		A.clear();
		cout << "A" << A;

		for (int i = 0; i < 5; ++i) {
			A.insert_former(i);
		}
		cout << "A" << A;

		A.move_start();
		while (A.size() != 0) {
			A.remove();
		}
		cout << "A" << A;

		for (int i = 0; i < 5; ++i) {
			A.insert_former(i);
		}
		cout << "A" << A;
	}
};