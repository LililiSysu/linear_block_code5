/*****************************************************************//**
 * \file   NewGE_in_OSD.h
 * \brief  test for prestore systematic matrix and update method, including the 
 *		   iterative basis update (IBU) and iterative basis update utilizing cyclic property (IBUc)
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once
#include"test_common.h"

 /**
  * .this class is for new idea testing
  *
  */
class test_NewGE_in_OSD {
public:

	static void GF2e_inv_test() {
		const int m = 4;
		typedef GF2e<m> ty;
		ty::init();

		ty tmp;
		ty tmp2;
		for (int i = 0; i < GF2e_auxiliary_storage::q_minus_1; ++i) {
			tmp.set_by_alpha_power(i);
			tmp2.set_by_alpha_power(-i);
			cout << "tmp = " << tmp << ",\t tmp2 = " << tmp2 << "\t tmp*tmp2 = " << tmp * tmp2 << endl;
		}
	}
	static void generating_all_n_choose_k_pattern_test() {
		int n = 6;
		int k = 3;

		Matrix<int> n_choose_k_pattern = Matrix_common::generating_all_n_choose_k_pattern(n, k);
		cout << "n_choose_k_pattern" << n_choose_k_pattern;

		cout << "---------------------" << endl;
		Matrix<int> n_choose_k_pattern_rev = Matrix_common::generating_all_n_choose_k_pattern_rev(n, k);
		cout << "n_choose_k_pattern_rev" << n_choose_k_pattern_rev;		// looking good

		cout << "---------------------" << endl;
		Matrix<int> TEP(1, 3, { 3,4,5 });
		while (TEP(0) != -1) {
			cout << "TEP" << TEP;
			Matrix_common::OSD_next_TEP(TEP, 6);						// done
		}

		cout << "---------------------" << endl;
		Matrix<int> TEP_2(1, 3, { 0,1,2 });								// which is silly and unright!
		cout << "TEP_2" << TEP_2;
		for (int i = 0; i < 19; ++i) {
			flip_TEP::next(TEP_2, 6, 3);
			cout << "TEP_2" << TEP_2;
		}
	}
	static void new_generation_of_G_residual_test() {
		const int m = 5;
		const int k_prime = 29;
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
		Matrix<int> info_set_ind_array(1, k_prime_dual, { 0,1 });
		Matrix<int> num_using_vectors(1, 10, '0');

		RS<m, k_prime_dual> rs_4_dual;
		RS<m, k_prime> rs;
		BCH<m, t> bch;
		int k = bch.get_k();
		int k_dual = n - k;

		int considered_info_bit = 2 * (k_dual - k_prime_dual);
		//int more_considered_info_bit = n - k_dual;

		Matrix<ty> testing_coordinate(1, m - 1, { 2, 4, 9, 16 });	// these corodinate generates all zero bits in k_prime_dual info bits

		tree_store_method<m> tree_space(considered_info_bit, k_prime_dual, testing_coordinate);

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

			tree_space.init(\
				rs_4_dual.generator_M_systematic_any_pos_dual.get_part(0, k_prime_dual, -1, k_prime_dual + considered_info_bit - 1));

			tree_space.get_generator();

			cout << "Q_sub_dual" << Q_sub_dual;

			cout << "tree_space.codewords" << tree_space.codewords;
			Matrix<GF2> cw = tree_space.codewords;
			cw.row_transformation_to_up_triangle();
			cout << "cw" << cw;

			// we can adding the 'v_dual_trace' into the Q_sub_dual to form Q_dual
			Matrix<GF2> Q_dual(k_prime_dual + tree_space.codewords.row(), n, '0');		// to be changed to k_dual
			Q_dual.set_part(0, 0, Q_sub_dual);
			//Q_dual.set_part(k_prime_dual, 0, tree_space.codewords);			
			//cout << "Q_dual (GF2)" << Q_dual;

			// verify this is the generator of (7,3) BCH dual
			cout << "bch.get_generator_matrix().multiply_transpose_of(Q_dual).isZero() = " << \
				bch.get_generator_matrix().multiply_transpose_of(Q_dual).isZero() << endl;
		}

		cout << "------------------------" << endl;
		cout << "num_using_vectors" << num_using_vectors;
		cout << "------------------------" << endl;
	}
	static void binary_tree_testing() {
		for (int i = 0; i < 100000000; ++i) {
			if (i % 100 == 0)
				cout << "---------- i = " << i << "---------" << endl;
			binary_tree* head = new binary_tree();
			head->left = new binary_tree();
			head->right = new binary_tree();
			head->left->left = new binary_tree();
			binary_tree::destroy_recur(head);
		}
	}

	static void new_gramma_test() {
		int a = 4;
		for (int i = 0; i < 3; ++i) {
			a++;
		}
		cout << "a = " << a << endl;
	}
	static void compare_generator_method() {
		const int m = 7;
		const int t = 10;
		const int d = (t << 1) + 1;
		const int n = (1 << m) - 1;
		const int k_prime = n - d + 1;
		typedef GF2e<m> ty;
		ty::init();

		RS<m, k_prime> rs;
		rs.print_info();
		unsigned long long before, after, before_mul, after_mul, before_add, after_add;

		// we generate the reliability order randomly
		Matrix<my_double> recv_abs(1, n);
		for (int i = 0; i < n; ++i) {
			recv_abs(i) = my::rand_u_adv();
		}
		Matrix<int> sorted_pos = recv_abs.sort_with_ind('>');
		cout << "sorted_pos" << sorted_pos;
		Matrix<int> info_set_k_prime = sorted_pos.get_part(0, 0, 0, k_prime - 1);

		// use Lagrange interpolation to generate generator of RS code
		before = GF2e_auxiliary_storage::operation_number;
		before_mul = GF2e_auxiliary_storage::mul_number;
		before_add = GF2e_auxiliary_storage::add_number;
		rs.generate_systematic_generator_any_pos_best(info_set_k_prime);
		cout << "rs" << rs.generator_M_systematic_any_pos << endl;

		after = GF2e_auxiliary_storage::operation_number;
		after_mul = GF2e_auxiliary_storage::mul_number;
		after_add = GF2e_auxiliary_storage::add_number;
		cout << "computation cost (add & mul) = " << after - before << endl;
		cout << "computation cost (mul) = " << after_mul - before_mul << endl;
		cout << "computation cost (add) = " << after_add - before_add << endl;

		// use Gaussian elimination to generate generator of BCH code
		BCH<m, t> bch;
		int k = bch.get_k();
		Matrix<int> natual_order(1, n, 'N');

		for (int test_ind = 0; test_ind < 10; test_ind++) {

			Matrix<GF2> H = bch.get_parity_matrix();
			for (int i = 0; i < n; ++i) {
				recv_abs(i) = my::rand_u_adv();
			}
			Matrix<int> sorted_pos = recv_abs.sort_with_ind('>');
			H.permute_col(sorted_pos);

			before = GF2_auxiliary_storage::operation_number;
			H.row_transformation_to_low_triangle();
			Matrix<int> permute_record = H.col_permute_to_full_rank_on_right();
			H.row_transformation_right_low_triangle_to_identity();
			cout << "permute_record - natual_order" << permute_record - natual_order;
			cout << "permute_record" << permute_record;
			cout << "H" << H;

			after = GF2_auxiliary_storage::operation_number;
			cout << "computation cost (add & mul) = " << after - before << endl;
			// it shows that the computation cost of Gaussian elimination is less than that of Lagrange intepolation
		}

	}
	static void PreStored_Matrix_test() {
		const int m = 4;
		const int t = 2;
		const int d = (t << 1) + 1;
		const int n = (1 << m) - 1;
		const int k_prime = n - d + 1;
		typedef GF2e<m> ty;
		ty::init();

		BCH<m, t> bch;
		bch.print_info();
		int k = bch.get_k();
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		cout << "G" << G;

		int G_partition_num = 3;
		PreStored_Matrix psm(G, G_partition_num);
		PreStored_Matrix_red psm_red(G, G_partition_num);

		// to test whether it preserve a systematic generator matrix
		//int set_num = psm.G_set.size();
		//for (int i = 0; i < set_num; ++i) {
		//	cout << "------------" << endl;
		//	cout << "psm.G_set(i)" << psm.G_set(i);
		//	cout << "psm.info_set_matrix(i)" << psm.permute_set(i);

		//	// permute G_set(i) back and check if it can pass the parity check matrix of bch code
		//	Matrix<GF2> G_tmp = psm.G_set(i);
		//	G_tmp.permute_col_back(psm.permute_set(i));
		//	cout << "G_tmp" << G_tmp;
		//	cout << " G_tmp.multiply_transpose_of(H)" << G_tmp.multiply_transpose_of(H);

		//	// perfectly, they all pass the test
		//}

		cout << "psm.sorted_info_set" << psm.sorted_info_set;
		Matrix<int> rand_sequence(1, n, 'N');
		rand_sequence.permute_rand();

		unsigned long long before, after;
		before = GF2_auxiliary_storage::operation_number;
		psm.get_MRIP_sys_G(rand_sequence);
		after = GF2_auxiliary_storage::operation_number;
		cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl;

		Matrix<GF2> G_tmp = psm.G_target;
		G_tmp.permute_col_back(psm.permute_target);
		cout << "G_tmp" << G_tmp;
		cout << " G_tmp.multiply_transpose_of(H)" << G_tmp.multiply_transpose_of(H);
		Matrix<int> info_set_sorted = psm.permute_target.get_part(0, 0, 0, k - 1);
		info_set_sorted.sort('<');
		cout << "//////////////////////////////////////////////////////////////////////////////////////////////////" << endl;

		before = GF2_auxiliary_storage::operation_number;
		psm_red.get_MRIP_sys_G(rand_sequence);
		after = GF2_auxiliary_storage::operation_number;
		cout << "------------------ computation cost (psm)(new) = " << after - before << " ---------------------" << endl;

		Matrix<GF2> G_tmp2 = psm.G_target;
		G_tmp2.permute_col_back(psm.permute_target);
		cout << "G_tmp2" << G_tmp2;
		cout << " G_tmp2.multiply_transpose_of(H)" << G_tmp2.multiply_transpose_of(H);
		Matrix<int> info_set_sorted_2 = psm.permute_target.get_part(0, 0, 0, k - 1);
		info_set_sorted_2.sort('<');
		cout << "//////////////////////////////////////////////////////////////////////////////////////////////////" << endl;

		// count the computation cost of Gaussian elimination

		before = GF2_auxiliary_storage::operation_number;
		Matrix<GF2> G_GE = G;
		G_GE.permute_col(rand_sequence);
		G_GE.row_transformation_to_up_triangle();
		Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
		G_GE.row_transformation_left_up_triangle_to_identity();
		after = GF2_auxiliary_storage::operation_number;
		cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
		Matrix<int> permute_GE_real = rand_sequence;
		permute_GE_real.permute(permute_GE);
		cout << "permute_GE_real" << permute_GE_real;
		cout << "G_GE" << G_GE;
		G_GE.permute_col_back(permute_GE_real);
		cout << "G_GE.multiply_transpose_of(H)" << G_GE.multiply_transpose_of(H);

		Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
		info_set_sorted_GE.sort('<');

		// compare the info_set
		cout << "info_set_sorted - info_set_sorted_GE" << info_set_sorted - info_set_sorted_GE;
		cout << "info_set_sorted_2 - info_set_sorted_GE" << info_set_sorted_2 - info_set_sorted_GE;
	}
	static void PreStored_Matrix_test_multi_trial() {
		const int m = 7;
		const int t = 10;
		typedef GF2e<m> ty;
		ty::init();

		eBCH<m, t> ebch;
		ebch.print_info();
		int k = ebch.get_k();
		int n = ebch.get_n();
		int d = ebch.get_d();
		Matrix<GF2> G = ebch.get_generator_matrix();
		Matrix<GF2> H = ebch.get_parity_matrix();
		cout << "G" << G;

		int G_partition_num = 2;
		PreStored_Matrix psm(G, G_partition_num);
		PreStored_Matrix psm_red(G, G_partition_num);

		// to test whether it preserve a systematic generator matrix
		//int set_num = psm.G_set.size();
		//for (int i = 0; i < set_num; ++i) {
		//	cout << "------------" << endl;
		//	cout << "psm.G_set(i)" << psm.G_set(i);
		//	cout << "psm.info_set_matrix(i)" << psm.permute_set(i);

		//	// permute G_set(i) back and check if it can pass the parity check matrix of bch code
		//	Matrix<GF2> G_tmp = psm.G_set(i);
		//	G_tmp.permute_col_back(psm.permute_set(i));
		//	cout << "G_tmp" << G_tmp;
		//	cout << " G_tmp.multiply_transpose_of(H)" << G_tmp.multiply_transpose_of(H);

		//	// perfectly, they all pass the test
		//}

		cout << "psm.sorted_info_set" << psm.sorted_info_set;

		int test_total = 40000;
		my_double computation_cost_psm_orig_ave = 0;
		my_double computation_cost_psm_new_ave = 0;
		my_double computation_cost_GE_ave = 0;

		unsigned long long computation_cost_psm_orig_max = 0;
		unsigned long long computation_cost_psm_new_max = 0;
		unsigned long long computation_cost_GE_max = 0;

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			//cout << "------------------------------------\n" << "test_ind=" << test_ind << "\n---------------------------------\n";

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			unsigned long long before, after;
			before = GF2_auxiliary_storage::operation_number;
			psm.get_MRIP_sys_G(rand_sequence);
			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl; 
			computation_cost_psm_orig_ave += after - before;
			computation_cost_psm_orig_max = my::max(computation_cost_psm_orig_max, after - before);

			Matrix<GF2> G_tmp = psm.G_target;
			//cout << "G_tmp" << G_tmp;
			G_tmp.permute_col_back(psm.permute_target);
			if (G_tmp.multiply_transpose_of(H).isZero() == false) {
				cout << "(psm)(orig) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted = psm.permute_target.get_part(0, 0, 0, k - 1);
			info_set_sorted.sort('<');

			before = GF2_auxiliary_storage::operation_number;
			psm_red.get_MRIP_sys_G(rand_sequence);
			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(new) = " << after - before << " ---------------------" << endl;
			computation_cost_psm_new_ave += after - before;
			computation_cost_psm_new_max = my::max(computation_cost_psm_new_max, after - before);

			Matrix<GF2> G_tmp_red = psm_red.G_target;
			//cout << "G_tmp2" << G_tmp2;
			G_tmp_red.permute_col_back(psm_red.permute_target);
			if (G_tmp_red.multiply_transpose_of(H).isZero() == false) {
				cout << "(psm)(new) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_red = psm_red.permute_target.get_part(0, 0, 0, k - 1);
			info_set_sorted_red.sort('<');

			// count the computation cost of Gaussian elimination

			before = GF2_auxiliary_storage::operation_number;
			Matrix<GF2> G_GE = G;
			G_GE.permute_col(rand_sequence);

			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2();		// new method

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
			computation_cost_GE_ave += after - before;
			computation_cost_GE_max = my::max(computation_cost_GE_max, after - before);

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
			//cout << "permute_GE_real" << permute_GE_real;
			//cout << "G_GE" << G_GE;
			G_GE.permute_col_back(permute_GE_real);
			if (G_GE.multiply_transpose_of(H).isZero() == false) {
				cout << "(GE) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
			info_set_sorted_GE.sort('<');

			//cout << "-----------" << endl;

			// compare the info_set

			//cout << "(info_set_sorted - info_set_sorted_GE).isZero() = " << (info_set_sorted - info_set_sorted_GE).isZero() << endl;
			//cout << "(info_set_sorted_2 - info_set_sorted_GE).isZero() = " << (info_set_sorted_2 - info_set_sorted_GE).isZero() << endl;
			if ((info_set_sorted - info_set_sorted_GE).isZero() == false) {
				cout << "(psm)(orig) not a MRB, test_ind = " << test_ind << endl;
			}
			if ((info_set_sorted_red - info_set_sorted_GE).isZero() == false) {
				cout << "(psm)(new) not a MRB, test_ind = " << test_ind << endl;
			}
		}

		computation_cost_psm_orig_ave /= test_total;
		computation_cost_psm_new_ave /= test_total;
		computation_cost_GE_ave /= test_total;

		cout << "computation_cost_psm_orig_ave = " << computation_cost_psm_orig_ave << endl;
		cout << "computation_cost_psm_new_ave = " << computation_cost_psm_new_ave << endl;
		cout << "computation_cost_GE_ave = " << computation_cost_GE_ave << endl;

		cout << "-----------" << endl;

		cout << "computation_cost_psm_orig_max = " << computation_cost_psm_orig_max << endl;
		cout << "computation_cost_psm_new_max = " << computation_cost_psm_new_max << endl;
		cout << "computation_cost_GE_max = " << computation_cost_GE_max << endl;
	}
	static void PreStored_Matrix_test_time() {
		const int m = 7;
		const int t = 10;
		typedef GF2e<m> ty;
		ty::init();

		eBCH<m, t> ebch;
		ebch.print_info();
		int k = ebch.get_k();
		int n = ebch.get_n();
		int d = ebch.get_d();
		Matrix<GF2> G = ebch.get_generator_matrix();
		//cout << "G" << G;

		int G_partition_num = 3;
		PreStored_Matrix psm(G, G_partition_num);

		// to test whether it preserve a systematic generator matrix
		//int set_num = psm.G_set.size();
		//for (int i = 0; i < set_num; ++i) {
		//	cout << "------------" << endl;
		//	cout << "psm.G_set(i)" << psm.G_set(i);
		//	cout << "psm.info_set_matrix(i)" << psm.permute_set(i);

		//	// permute G_set(i) back and check if it can pass the parity check matrix of bch code
		//	Matrix<GF2> G_tmp = psm.G_set(i);
		//	G_tmp.permute_col_back(psm.permute_set(i));
		//	cout << "G_tmp" << G_tmp;
		//	cout << " G_tmp.multiply_transpose_of(H)" << G_tmp.multiply_transpose_of(H);

		//	// perfectly, they all pass the test
		//}

		//cout << "psm.sorted_info_set" << psm.sorted_info_set;

		int test_total = 40000;
		my_double computation_cost_psm_orig_ave = 0;
		my_double computation_cost_psm_new_ave = 0;
		my_double computation_cost_GE_ave = 0;

		unsigned long long computation_cost_psm_orig_max = 0;
		unsigned long long computation_cost_psm_new_max = 0;
		unsigned long long computation_cost_GE_max = 0;

		clock_t start, end;
		start = clock();
		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();
			psm.get_MRIP_sys_G(rand_sequence);
		}
		end = clock();
		double time_consume_psm_orig = ((double)end - start) / CLOCKS_PER_SEC / (double)test_total;


		cout << scientific << setprecision(2);
		cout << "time consume (psm) (oirg) = " << time_consume_psm_orig << "s/iteration" << endl;

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void PreStored_Matrix_red_test_time() {
		const int m = 7;
		const int t = 10;
		typedef GF2e<m> ty;
		ty::init();

		eBCH<m, t> ebch;
		ebch.print_info();
		int k = ebch.get_k();
		int n = ebch.get_n();
		int d = ebch.get_d();
		Matrix<GF2> G = ebch.get_generator_matrix();
		//cout << "G" << G;

		int G_partition_num = 3;

		int test_total = 40000;
		my_double computation_cost_psm_orig_ave = 0;
		my_double computation_cost_psm_new_ave = 0;
		my_double computation_cost_GE_ave = 0;

		unsigned long long computation_cost_psm_orig_max = 0;
		unsigned long long computation_cost_psm_new_max = 0;
		unsigned long long computation_cost_GE_max = 0;

		clock_t start, end;

		PreStored_Matrix_red psm_red(G, G_partition_num);
		start = clock();
		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();
			psm_red.get_MRIP_sys_G(rand_sequence);
		}
		end = clock();
		double time_consume_psm_new = ((double)end - start) / CLOCKS_PER_SEC / (double)test_total;

		cout << scientific << setprecision(2);
		cout << "time consume (psm) (new) = " << time_consume_psm_new << "s/iteration" << endl;

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void GE_test_time() {
		const int m = 7;
		const int t = 10;
		typedef GF2e<m> ty;
		ty::init();

		eBCH<m, t> ebch;
		ebch.print_info();
		int k = ebch.get_k();
		int n = ebch.get_n();
		int d = ebch.get_d();
		Matrix<GF2> G = ebch.get_generator_matrix();
		//cout << "G" << G;


		int test_total = 40000;
		my_double computation_cost_psm_orig_ave = 0;
		my_double computation_cost_psm_new_ave = 0;
		my_double computation_cost_GE_ave = 0;

		unsigned long long computation_cost_psm_orig_max = 0;
		unsigned long long computation_cost_psm_new_max = 0;
		unsigned long long computation_cost_GE_max = 0;

		clock_t start, end;

		start = clock();
		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			Matrix<GF2> G_GE = G;
			G_GE.permute_col(rand_sequence);
			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2_echelon();		// simplified method for GE

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
		}
		end = clock();
		double time_consume_GE = ((double)end - start) / CLOCKS_PER_SEC / (double)test_total;

		cout << scientific << setprecision(2);
		cout << "time consume (GE) = " << time_consume_GE << "s/iteration" << endl;

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void GE_left_identity_4_GF2_test() {

		my::random_engine.seed((unsigned int)std::time(nullptr));

		//Matrix<GF2> A(6, 6, {
		//	0,0,0,1,0,1,
		//	1,1,0,1,1,1,
		//	0,0,1,0,1,0,
		//	1,0,1,1,1,1,
		//	0,0,0,1,0,1,
		//	0,1,1,0,0,1
		//});

		/*Matrix<GF2> A(64, 128);
		for (int i = 0, imax = A.size(); i < imax; ++i) {
			A(i) = my::rand_01_adv();
		}*/

		const int m = 7;
		const int t = 10;
		typedef GF2e<m> ty;
		ty::init();

		eBCH<m, t> ebch;
		int k = ebch.get_k();
		Matrix<GF2> A = ebch.get_generator_matrix();
		Matrix<GF2> H = ebch.get_parity_matrix();

		// permute A's column randomly
		Matrix<int> col_permutation(1, ebch.get_n(), 'N');
		col_permutation.permute_rand();
		A.permute_col(col_permutation);

		//cout << "A" << A;

		unsigned long long before, after;

		Matrix<GF2> A_new_fun_test = A;
		before = GF2_auxiliary_storage::operation_number;

		Matrix<int> permute_record_new = A_new_fun_test.GE_left_identity_4_GF2();

		after = GF2_auxiliary_storage::operation_number;
		//cout << "A_new_fun_test" << A_new_fun_test;
		unsigned long long specific_cost = after - before;
		cout << "computation number for specific method = " << specific_cost << endl;

		// compare it with standard GE
		Matrix<GF2> A_standard_fun_test = A;
		before = GF2_auxiliary_storage::operation_number;

		A_standard_fun_test.row_transformation_to_up_triangle();
		Matrix<int> permute_record_standard = A_standard_fun_test.col_permute_to_full_rank_on_left();
		A_standard_fun_test.row_transformation_left_up_triangle_to_identity();

		after = GF2_auxiliary_storage::operation_number;
		//cout << "A_standard_fun_test" << A_standard_fun_test;
		unsigned long long standard_cost = after - before;
		cout << "computation number for standard method = " << standard_cost << endl;

		// count computation cost percentage
		cout << "specific_cost / standard_cost = " << (double)specific_cost / standard_cost << endl;

		cout << "----" << endl;
		cout << "(A_new_fun_test - A_standard_fun_test).isZero() = "
			<< (A_new_fun_test - A_standard_fun_test).isZero() << endl;

		cout << "(permute_record_new - permute_record_standard).get_part(0,0,0,k-1).isZero() = "
			<< (permute_record_new - permute_record_standard).get_part(0, 0, 0, k - 1).isZero() << endl;

		A_new_fun_test.permute_col_back(permute_record_new);
		A_new_fun_test.permute_col_back(col_permutation);

		A_standard_fun_test.permute_col_back(permute_record_standard);
		A_standard_fun_test.permute_col_back(col_permutation);

		cout << "A_new_fun_test.multiply_transpose_of(H).isZero() = " << A_new_fun_test.multiply_transpose_of(H).isZero() << endl;
		cout << "A_standard_fun_test.multiply_transpose_of(H).isZero() = " << A_standard_fun_test.multiply_transpose_of(H).isZero() << endl;
	}
	static void cycle_code_generator_test() {
		const int m = 3;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 1;
		BCH<m, t> bch;		// ebch code is not a cyclic code, we should break it into a cyclic bch code and a parity check position
		bch.print_info();

		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();

		cout << "G" << G;
		cout << "H" << H;

		cout << "G.multiply_transpose_of(H).isZero() = " << G.multiply_transpose_of(H).isZero() << endl;

		// test for cycle property of G
		Matrix<GF2> G_shift;
		G.col_shift_right_cir(-1, G_shift);

		cout << "G_shift" << G_shift;
		cout << "G_shift.multiply_transpose_of(H).isZero() = " << G_shift.multiply_transpose_of(H).isZero() << endl;

		// test for cycle property of H
		Matrix<GF2> H_shift;
		H.col_shift_right_cir(-3, H_shift);

		cout << "H_shift" << H_shift;
		cout << "H_shift.multiply_transpose_of(G).isZero() = " << H_shift.multiply_transpose_of(G).isZero() << endl;

		Matrix<GF2> G_sys = G;
		G_sys.GE_left_identity_4_GF2_echelon();

		// test for ind*2 property of G
		Matrix<GF2> G_ind_mul;
		G_sys.col_ind_mul_cir(2, G_ind_mul);

		cout << "G_ind_mul" << G_ind_mul;
		cout << "G_ind_mul.multiply_transpose_of(H).isZero() = " << G_ind_mul.multiply_transpose_of(H).isZero() << endl;

		Matrix<GF2> G_sys_new = G;
		int sss = G_sys_new.col();
		Matrix<int> m_ind(1, sss, 'N');
		my::set_seed_adv(22);
		m_ind = m_ind.get_random_element(sss);
		G_sys_new.permute_col(m_ind);				// I believe permute randomly in column to get new info set is not so meaningful

		Matrix<int> second_permute = G_sys_new.GE_left_identity_4_GF2_echelon();
		cout << "m_ind" << m_ind;
		cout << "second_permute" << second_permute;

		cout << "G_sys_new" << G_sys_new;
		m_ind.permute(second_permute);
		G_sys_new.permute_col_back(m_ind);

		cout << "G_sys_new (permute_col_back)" << G_sys_new;
		cout << "G_sys_new.multiply_transpose_of(H).isZero() = " << G_sys_new.multiply_transpose_of(H).isZero() << endl;

		G_sys_new.col_shift_right_cir(-6, G_shift);
		cout << "G_shift" << G_shift;
		cout << "G_shift.multiply_transpose_of(H).isZero() = " << G_shift.multiply_transpose_of(H).isZero() << endl;

		G_sys_new.col_ind_mul_cir(2, G_ind_mul);
		// mul 2 or 4 will give the info set that has travelled by not multiplied, as you see the more matrix we added
		// the less new the info set can be generated
		cout << "G_ind_mul" << G_ind_mul;
		cout << "G_ind_mul.multiply_transpose_of(H).isZero() = " << G_ind_mul.multiply_transpose_of(H).isZero() << endl;
	}
	static void density_of_BCH_generator_matrix() {
		const int m = 7;
		const int t = 14;
		GF2e<m>::init();
		BCH<m, t> bch;
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		cout << "-----G--------" << endl;
		//cout << "G" << G;
		int k = G.row();
		int n = G.col();

		// count density of 1 in G
		cout << "density_of_1 (G) = " << G.density() << endl;

		Matrix<GF2> G_sys = G;
		G_sys.GE_left_identity_4_GF2();
		cout << "-----G_sys--------" << endl;
		//cout << "G_sys" << G_sys;

		// count density of 1 in G_sys
		cout << "density_of_1 (G_sys) = " << G_sys.density() << endl;
		cout << "density_of_1 (G_sys parity part) = " << G_sys.density(0, k) << endl;

		// find the minimum density of rows

		Matrix<int> num_of_1_in_row(1, k, '0');
		for (int i = 0; i < k; ++i) {
			for (int j = 0; j < n; ++j) {
				if (G_sys(i, j) == 1) {
					num_of_1_in_row(i)++;
				}
			}
		}

		Matrix<double> density_of_row(1, k);
		for (int i = 0; i < k; ++i) {
			density_of_row(i) = (double)num_of_1_in_row(i) / n;
		}

		//cout << "num_of_1_in_row" << num_of_1_in_row;
		//cout << "density_of_row" << density_of_row;

		// generate the generator matrix by the row of minimum density
		int row_min_density_ind = num_of_1_in_row.min_ele_ind();
		Matrix<GF2> G_low_density(k, n, '0');
		for (int i = 0; i < k; ++i) {
			for (int j = 0; j < n; ++j) {
				G_low_density(i, (j + i) % n) = G_sys(row_min_density_ind, (row_min_density_ind + j) % n);
			}
		}
		cout << "-----G_low_density--------" << endl;
		//cout << "G_low_density" << G_low_density;
		cout << "G_low_density.multiply_transpose_of(H).isZero() = " << G_low_density.multiply_transpose_of(H).isZero() << endl;

		// count density of 1 in G_sys		
		cout << "density_of_1 (G_low_density) = " << G_low_density.density() << endl;

		// perform GE to all 3 matrix and find which induce the lowest average operation number
		unsigned long long G_ope = 0, G_sys_ope = 0, G_low_density_ope = 0;
		Matrix<GF2> G_ans;
		ope_count start, end;
		int test_num = 1;
		for (int i = 0; i < test_num; ++i) {
			// random permute the column of generator matrix to imitate a received sequence reliability sorting
			Matrix<int> natual(1, n, 'N');
			natual.permute_rand();

			start.count();
			G_ans = G;
			G_ans.permute_col(natual);
			G_ans.GE_left_identity_4_GF2_echelon();
			end.count();

			G_ope += end.GF2_ope - start.GF2_ope;

			/*start.count();
			G_ans = G_sys;
			G_ans.permute_col(natual);
			G_ans.GE_left_identity_4_GF2_echelon();
			end.count();

			G_sys_ope += end.GF2_ope - start.GF2_ope;

			start.count();
			G_ans = G_low_density;
			G_ans.permute_col(natual);
			G_ans.GE_left_identity_4_GF2_echelon();
			end.count();

			G_low_density_ope += end.GF2_ope - start.GF2_ope;*/
		}

		// taking average of operation numner
		double G_ope_ave, G_sys_ope_ave, G_low_density_ope_ave;
		G_ope_ave = (double)G_ope / test_num;
		G_sys_ope_ave = (double)G_sys_ope / test_num;
		G_low_density_ope_ave = (double)G_low_density_ope / test_num;

		cout << "-----ope_count-----" << endl;
		cout << "G_ope_ave = " << G_ope_ave << endl;
		cout << "G_sys_ope_ave = " << G_sys_ope_ave << endl;
		cout << "G_low_density_ope_ave = " << G_low_density_ope_ave << endl;
	}
	static void density_of_BCH_generator_matrix_2() {
		const int m = 7;
		const int t = 31;
		GF2e<m>::init();
		BCH<m, t> bch;
		bch.print_info();

		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		cout << "-----G--------" << endl;
		//cout << "G" << G;
		int k = G.row();
		int n = G.col();

		// count density of 1 in G
		cout << "density_of_1 (G) = " << G.density() << endl;

		Matrix<GF2> G_sys = G;
		G_sys.GE_left_identity_4_GF2();
		cout << "-----G_sys--------" << endl;
		//cout << "G_sys" << G_sys;

		// count density of 1 in G_sys
		cout << "density_of_1 (G_sys) = " << G_sys.density() << endl;
		cout << "density_of_1 (G_sys parity part) = " << G_sys.density(0, k) << endl;

	}
	static void check_out_eBCH_test() {
		const int m = 5;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 3;
		eBCH<m, t> ebch;
		ebch.print_info();

		Matrix<GF2> G = ebch.get_generator_matrix();
		Matrix<GF2> H = ebch.get_parity_matrix();

		cout << "G" << G;
		cout << "H" << H;

		cout << "G.multiply_transpose_of(H).isZero() = " << G.multiply_transpose_of(H).isZero() << endl;

		int k = G.row();
		int n = G.col();
		int r = H.row();

		// test for cycle property of G
		Matrix<GF2> G_shift;
		G.get_part(0, 0, k - 1, n - 2).col_shift_right_cir(10, G_shift);
		G_shift = G_shift.combine_right(G.get_col(n - 1));

		cout << "G_shift" << G_shift;
		cout << "G_shift.multiply_transpose_of(H).isZero() = " << G_shift.multiply_transpose_of(H).isZero() << endl;


		// test for cycle property of H
		Matrix<GF2> H_shift;
		H.get_part(0, 0, r - 2, n - 2).col_shift_right_cir(5, H_shift);
		H_shift = H_shift.combine(Matrix<GF2>(r - 1, 1, '0'), H.get_part(-1, 0, -1, n - 2), (GF2)1);

		cout << "H_shift" << H_shift;
		cout << "H_shift.multiply_transpose_of(H).isZero() = " << H_shift.multiply_transpose_of(G).isZero() << endl;
	}

	static void num_test() {
		cout << (double)127 / 64 << endl;
	}
	static void number_sequence_test() {
		double p = 0.181102;
		double sum_of_p = 0;
		int k = 85;
		for (int i = 1; i < k; ++i) {
			cout << "i = " << i << ", p = " << p << endl;
			sum_of_p += p;
			p = p * (1 - p) * (1 + 2 * p);
		}
		cout << "sum_of_p = " << sum_of_p << endl;
		cout << "sum_of_p * k = " << sum_of_p * k << endl;
	}
	static void GE_worst_case_ope_estimation_test() {
		int n = 63;
		int k = 18;
		cout << "XOR ope estimation = " << double(1) / 2 * k * (k - 1) * (n - double(2) / 3 * k - double(1) / 3) << endl;
	}
	static void BU_worst_case_ope_estimation_test() {
		int n = 63;
		int k = 18;
		cout << "XOR ope estimation = " << double(1) / 2 * k * k * (n - 2 * k + k * k / double(n)) << endl;
	}
	static void PreStored_Matrix_red_recursive_partition_test_ini() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 10;
		eBCH<m, t> bch;
		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();

		int partition_level = 2;		// do not be larger than 3, the storage size grows fast
		PreStored_Matrix_red_recursive_partition psm_red_recur_part(G, partition_level);
	}
	static void PreStored_Matrix_red_recursive_partition_test_multi_trial() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 10;
		BCH<m, t> bch;
		bch.print_info();

		//Matrix<GF2> G = bch.get_generator_matrix();
		//Matrix<GF2> H = bch.get_parity_matrix();

		Matrix<GF2> G = bch.get_parity_matrix();		// inspect on dual code
		Matrix<GF2> H = bch.get_generator_matrix();

		int n = bch.get_n();
		int k = G.row();								// no problem

		int partition_level = 1;		// do not be larger than 3, the storage size grows fast
		PreStored_Matrix_red_recursive_partition psm_red_recur_part(G, partition_level);

		int test_total = 400;
		my_double computation_cost_psm_red_recur_part_ave = 0;
		my_double computation_cost_GE_ave = 0;

		unsigned long long computation_cost_psm_red_recur_part_max = 0;
		unsigned long long computation_cost_GE_max = 0;

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			//cout << "------------------------------------\n" << "test_ind=" << test_ind << "\n---------------------------------\n";

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			unsigned long long before, after;
			before = GF2_auxiliary_storage::operation_number;
			psm_red_recur_part.get_MRIP_sys_G(rand_sequence);
			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl; 
			computation_cost_psm_red_recur_part_ave += after - before;
			computation_cost_psm_red_recur_part_max = my::max(computation_cost_psm_red_recur_part_max, after - before);

			Matrix<GF2> G_tmp = psm_red_recur_part.G_target;
			//cout << "G_tmp" << G_tmp;
			G_tmp.permute_col_back(psm_red_recur_part.permute_target);
			if (G_tmp.multiply_transpose_of(H).isZero() == false) {
				cout << "(psm_red_recur_part) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted = psm_red_recur_part.permute_target.get_part(0, 0, 0, k - 1);
			info_set_sorted.sort('<');

			// count the computation cost of Gaussian elimination

			before = GF2_auxiliary_storage::operation_number;
			Matrix<GF2> G_GE = G;
			G_GE.permute_col(rand_sequence);

			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2();		// new method

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
			computation_cost_GE_ave += after - before;
			computation_cost_GE_max = my::max(computation_cost_GE_max, after - before);

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
			//cout << "permute_GE_real" << permute_GE_real;
			//cout << "G_GE" << G_GE;
			G_GE.permute_col_back(permute_GE_real);
			if (G_GE.multiply_transpose_of(H).isZero() == false) {
				cout << "(GE) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
			info_set_sorted_GE.sort('<');

			//cout << "-----------" << endl;

			// compare the info_set

			//cout << "info_set_sorted - info_set_sorted_GE" << info_set_sorted - info_set_sorted_GE << endl;
			if ((info_set_sorted - info_set_sorted_GE).isZero() == false) {
				cout << "(psm_red_recur_part) not a MRB, test_ind = " << test_ind << endl;
			}
		}

		computation_cost_psm_red_recur_part_ave /= test_total;
		computation_cost_GE_ave /= test_total;

		cout << "computation_cost_psm_red_recur_part_ave = " << computation_cost_psm_red_recur_part_ave << endl;
		cout << "computation_cost_GE_ave = " << computation_cost_GE_ave << endl;

		cout << "-----------" << endl;

		cout << "computation_cost_psm_red_recur_part_max = " << computation_cost_psm_red_recur_part_max << endl;
		cout << "computation_cost_GE_max = " << computation_cost_GE_max << endl;
	}
	static void n_ary_incremental_scan_test() {
		Matrix<int> ans = Matrix_common::n_ary_incremental_scan(2, 1);
		cout << "ans" << ans;
	}
	static void n_randomly_divide_into_k_and_res_test() {
		Matrix<int> k_part, n_res;
		Matrix_common::n_randomly_divide_into_k_and_res(10, 8, k_part, n_res);
		cout << "k_part" << k_part;
		cout << "n_res" << n_res;
	}

	static void PreStored_Matrix_red_recursive_partition_test_time() {
		const int m = 7;
		const int t = 10;
		typedef GF2e<m> ty;
		ty::init();

		eBCH<m, t> ebch;
		ebch.print_info();
		int k = ebch.get_k();
		int n = ebch.get_n();
		int d = ebch.get_d();
		Matrix<GF2> G = ebch.get_generator_matrix();
		//cout << "G" << G;

		int G_partition_num = 3;

		int test_total = 40000;
		my_double computation_cost_psm_orig_ave = 0;
		my_double computation_cost_psm_new_ave = 0;
		my_double computation_cost_GE_ave = 0;

		unsigned long long computation_cost_psm_orig_max = 0;
		unsigned long long computation_cost_psm_new_max = 0;
		unsigned long long computation_cost_GE_max = 0;

		clock_t start, end;

		PreStored_Matrix_red_recursive_partition psm_red_recur_part(G, G_partition_num);
		start = clock();
		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();
			psm_red_recur_part.get_MRIP_sys_G(rand_sequence);
		}
		end = clock();
		double time_consume_psm_new = ((double)end - start) / CLOCKS_PER_SEC / (double)test_total;

		cout << scientific << setprecision(2);
		cout << "time consume (psm) (new) = " << time_consume_psm_new << "s/iteration" << endl;

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void PreStored_Matrix_red_loose_bin_test_ini() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 10;
		BCH<m, t> bch;
		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> G = bch.get_parity_matrix();
		Matrix<GF2> H = bch.get_generator_matrix();

		int partition_level = 2;		// do not be larger than 3, the storage size grows fast
		PreStored_Matrix_red_loose_bin psm_red_lb(G, partition_level);
	}
	static void PreStored_Matrix_red_loose_bin_test_multi_trial() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 10;
		BCH<m, t> bch;
		bch.print_info();

		//Matrix<GF2> G = bch.get_generator_matrix();
		//Matrix<GF2> H = bch.get_parity_matrix();

		Matrix<GF2> G = bch.get_parity_matrix();		// inspect on dual code
		Matrix<GF2> H = bch.get_generator_matrix();

		int n = bch.get_n();
		int k = G.row();								// no problem

		int partition_level = 2;		// do not be larger than 3, the storage size grows fast
		PreStored_Matrix_red_loose_bin psm_red_lb(G, partition_level);

		int test_total = 6;
		Matrix<my_double> computation_cost_psm(1, test_total, 'v');
		Matrix<my_double> computation_cost_GE(1, test_total, 'v');

		Matrix<GF2> G_psm(k, n);
		Matrix<int> permute_psm(1, k);		// the result from 'PreStored_Matrix_red_loose_bin'
		Matrix<int> max_common_psm(1, test_total, 'v');

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			//cout << "------------------------------------\n" << "test_ind=" << test_ind << "\n---------------------------------\n";

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			unsigned long long before, after;
			before = GF2_auxiliary_storage::operation_number;

			int max_common_psm_tmp = psm_red_lb.get_MRIP_sys_G(rand_sequence, G_psm, permute_psm);
			max_common_psm.push_back(max_common_psm_tmp);

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl; 
			computation_cost_psm.push_back(after - before);

			//cout << "G_target" << G_target;
			G_psm.permute_col_back(permute_psm);
			if (G_psm.multiply_transpose_of(H).isZero() == false) {
				cout << "(psm_red_recur_part) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted = permute_psm.get_part(0, 0, 0, k - 1);
			info_set_sorted.sort('<');

			// count the computation cost of Gaussian elimination

			before = GF2_auxiliary_storage::operation_number;
			Matrix<GF2> G_GE = G;
			G_GE.permute_col(rand_sequence);

			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2();		// new method

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
			computation_cost_GE.push_back(after - before);

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
			//cout << "permute_GE_real" << permute_GE_real;
			//cout << "G_GE" << G_GE;
			G_GE.permute_col_back(permute_GE_real);
			if (G_GE.multiply_transpose_of(H).isZero() == false) {
				cout << "(GE) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
			info_set_sorted_GE.sort('<');

			//cout << "-----------" << endl;

			// compare the info_set

			//cout << "info_set_sorted - info_set_sorted_GE" << info_set_sorted - info_set_sorted_GE << endl;
			if ((info_set_sorted - info_set_sorted_GE).isZero() == false) {
				cout << "(psm_red_recur_part) not a MRB, test_ind = " << test_ind << endl;
			}
		}

		cout << "-----------" << endl;
		cout << "max_common_psm" << max_common_psm;

		cout << "-----------" << endl;
		cout << "computation_cost_psm (loose bin)" << computation_cost_psm;
		cout << "computation_cost_GE (loose bin)" << computation_cost_GE;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_ave (loose bin) = " << Matrix_common::ave(computation_cost_psm) << endl;
		cout << "computation_cost_GE_ave (loose bin) = " << Matrix_common::ave(computation_cost_GE) << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_max (loose bin) = " << Matrix_common::max(computation_cost_psm) << endl;
		cout << "computation_cost_GE_max (loose bin)  = " << Matrix_common::max(computation_cost_GE) << endl;
	}

	// with partition_pattern_number = 1, this is IBU in the paper
	static void PreStored_Matrix_red_equal_dist_test_ini() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 10;
		BCH<m, t> bch;
		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> G = bch.get_parity_matrix();
		Matrix<GF2> H = bch.get_generator_matrix();

		/*Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();*/

		int partition_pattern_num = 10;
		PreStored_Matrix_red_equal_dist psm_red_ed(G, partition_pattern_num);
	}
	static void PreStored_Matrix_red_equal_dist_test_multi_trial() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 3;
		BCH<m, t> bch;
		bch.print_info();

		/*Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();*/

		Matrix<GF2> G = bch.get_parity_matrix();		// inspect on dual code
		Matrix<GF2> H = bch.get_generator_matrix();

		int n = bch.get_n();
		int k = G.row();								// no problem

		int partition_pattern_num = 1;
		PreStored_Matrix_red_equal_dist psm_red_ed(G, partition_pattern_num);

		int test_total = 100000;						// run 10^5 times
		Matrix<my_double> computation_cost_psm(1, test_total, 'v');
		Matrix<my_double> computation_cost_GE(1, test_total, 'v');

		Matrix<GF2> G_psm(k, n);
		Matrix<int> permute_psm(1, k);		// the result from 'PreStored_Matrix_red_loose_bin'
		Matrix<int> max_common_psm(1, test_total, 'v');

		//my::set_seed_adv(20);

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			//cout << "------------------------------------\n" << "test_ind=" << test_ind << "\n---------------------------------\n";

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			unsigned long long before, after;
			before = GF2_auxiliary_storage::operation_number;

			int max_common_psm_tmp = psm_red_ed.get_MRIP_sys_G(rand_sequence, G_psm, permute_psm);
			max_common_psm.push_back(max_common_psm_tmp);

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl; 
			computation_cost_psm.push_back(after - before);

			//cout << "G_target" << G_target;
			G_psm.permute_col_back(permute_psm);
			if (G_psm.multiply_transpose_of(H).isZero() == false) {
				cout << "(psm_red_recur_part) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted = permute_psm.get_part(0, 0, 0, k - 1);
			//cout << "info_set_sorted (before sorted)" << info_set_sorted;
			info_set_sorted.sort('<');

			// count the computation cost of Gaussian elimination

			before = GF2_auxiliary_storage::operation_number;
			Matrix<GF2> G_GE = G;
			//Matrix<GF2> G_GE = psm_red_ed.G_set(0);					// use the left-identity matrix here, NOTICE THAT
			G_GE.permute_col(rand_sequence);

			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			//Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2();			// new method
			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2_echelon();		// echelon method, 1% less of XOR operation number

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
			computation_cost_GE.push_back(after - before);

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
			//cout << "permute_GE_real" << permute_GE_real;
			//cout << "G_GE" << G_GE;
			G_GE.permute_col_back(permute_GE_real);
			if (G_GE.multiply_transpose_of(H).isZero() == false) {
				cout << "(GE) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
			//cout << "info_set_sorted_GE (before sorted)" << info_set_sorted_GE;
			info_set_sorted_GE.sort('<');

			//cout << "-----------" << endl;

			// compare the info_set

			//cout << "info_set_sorted - info_set_sorted_GE" << info_set_sorted - info_set_sorted_GE << endl;
			if ((info_set_sorted - info_set_sorted_GE).isZero() == false) {
				cout << "(psm_red_recur_part) not a MRB, test_ind = " << test_ind << endl;
			}
		}

		//cout << "-----------" << endl;
		//cout << "max_common_psm" << max_common_psm;

		//cout << "-----------" << endl;
		//cout << "computation_cost_psm (equal dist)" << computation_cost_psm;
		//cout << "computation_cost_GE (equal dist)" << computation_cost_GE;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_ave (equal dist) = " << Matrix_common::ave(computation_cost_psm) << endl;
		cout << "computation_cost_GE_ave (equal dist) = " << Matrix_common::ave(computation_cost_GE) << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_max (equal dist) = " << Matrix_common::max(computation_cost_psm) << endl;
		cout << "computation_cost_GE_max (equal dist)  = " << Matrix_common::max(computation_cost_GE) << endl;
	}

	// compare to IBU, fill the n coordinate indices by shift of k info
	static void PreStored_Matrix_red_equal_dist_fill_test_ini() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 10;
		BCH<m, t> bch;
		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();

		//Matrix<GF2> G = bch.get_parity_matrix();
		//Matrix<GF2> H = bch.get_generator_matrix();

		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();

		int partition_pattern_num = 4;
		PreStored_Matrix_red_equal_dist_fill psm_red_ed_fill(G, partition_pattern_num);
	}
	static void PreStored_Matrix_red_equal_dist_fill_test_multi_trial() {
		const int m = 6;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 5;
		BCH<m, t> bch;
		bch.print_info();

		//Matrix<GF2> G = bch.get_generator_matrix();
		//Matrix<GF2> H = bch.get_parity_matrix();

		Matrix<GF2> G = bch.get_parity_matrix();	// inspect on dual code
		Matrix<GF2> H = bch.get_generator_matrix();

		int n = bch.get_n();
		int k = G.row();							// no problem

		int partition_pattern_num = 4;				// 15 at top for (127,64,22) codes, 7 at top for (63,45,7) codes 
		PreStored_Matrix_red_equal_dist_fill psm_red_ed_fill(G, partition_pattern_num);

		int test_total = 400;
		Matrix<my_double> computation_cost_psm(1, test_total, 'v');
		Matrix<my_double> computation_cost_GE(1, test_total, 'v');

		Matrix<GF2> G_psm(k, n);
		Matrix<int> permute_psm(1, k);		// the result from 'PreStored_Matrix_red_loose_bin'
		Matrix<int> max_common_psm(1, test_total, 'v');

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			//cout << "------------------------------------\n" << "test_ind=" << test_ind << "\n---------------------------------\n";

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			unsigned long long before, after;
			before = GF2_auxiliary_storage::operation_number;

			int max_common_psm_tmp = psm_red_ed_fill.get_MRIP_sys_G(rand_sequence, G_psm, permute_psm);
			max_common_psm.push_back(max_common_psm_tmp);

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl; 
			computation_cost_psm.push_back(after - before);

			//cout << "G_target" << G_target;
			G_psm.permute_col_back(permute_psm);
			if (G_psm.multiply_transpose_of(H).isZero() == false) {
				cout << "(psm_red_recur_part) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted = permute_psm.get_part(0, 0, 0, k - 1);
			info_set_sorted.sort('<');

			// count the computation cost of Gaussian elimination

			before = GF2_auxiliary_storage::operation_number;
			Matrix<GF2> G_GE = G;
			G_GE.permute_col(rand_sequence);

			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2();		// new method

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
			computation_cost_GE.push_back(after - before);

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
			//cout << "permute_GE_real" << permute_GE_real;
			//cout << "G_GE" << G_GE;
			G_GE.permute_col_back(permute_GE_real);
			if (G_GE.multiply_transpose_of(H).isZero() == false) {
				cout << "(GE) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
			info_set_sorted_GE.sort('<');

			//cout << "-----------" << endl;

			// compare the info_set

			//cout << "info_set_sorted - info_set_sorted_GE" << info_set_sorted - info_set_sorted_GE << endl;
			if ((info_set_sorted - info_set_sorted_GE).isZero() == false) {
				cout << "(psm_red_recur_part) not a MRB, test_ind = " << test_ind << endl;
			}
		}

		//cout << "-----------" << endl;
		//cout << "max_common_psm" << max_common_psm;

		//cout << "-----------" << endl;
		//cout << "computation_cost_psm (equal dist)" << computation_cost_psm;
		//cout << "computation_cost_GE (equal dist)" << computation_cost_GE;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_ave (equal dist) = " << Matrix_common::ave(computation_cost_psm) << endl;
		cout << "computation_cost_GE_ave (equal dist) = " << Matrix_common::ave(computation_cost_GE) << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_max (equal dist) = " << Matrix_common::max(computation_cost_psm) << endl;
		cout << "computation_cost_GE_max (equal dist)  = " << Matrix_common::max(computation_cost_GE) << endl;
	}

	// only using shift property of cyclic code, information set shift by {0,1,...,n-1}. Not using multiplication by {2^0,2^1,...,2^s}
	static void PreStored_Matrix_red_equal_dist_cycle_test_ini() {
		const int m = 6;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 3;
		BCH<m, t> bch;
		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> G = bch.get_parity_matrix();
		Matrix<GF2> H = bch.get_generator_matrix();

		int  partition_pattern_num = 2;
		PreStored_Matrix_red_equal_dist_cycle psm_red_ed_ext_cyc(G, partition_pattern_num);
	}
	static void PreStored_Matrix_red_equal_dist_cycle_test_multi_trial() {
		const int m = 6;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 3;
		BCH<m, t> bch;
		bch.print_info();

		/*Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();*/

		Matrix<GF2> G = bch.get_parity_matrix();		// inspect on dual code
		Matrix<GF2> H = bch.get_generator_matrix();

		int n = bch.get_n();
		int k = G.row();								// no problem

		int partition_pattern_num = 4;
		PreStored_Matrix_red_equal_dist_cycle psm_red_ed_cyc(G, partition_pattern_num);

		int test_total = 400;
		Matrix<my_double> computation_cost_psm(1, test_total, 'v');
		Matrix<my_double> computation_cost_GE(1, test_total, 'v');

		Matrix<GF2> G_psm(k, n);
		Matrix<int> permute_psm(1, k);		// the result from 'PreStored_Matrix_red_loose_bin'
		Matrix<int> max_common_psm(1, test_total, 'v');

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			//cout << "------------------------------------\n" << "test_ind=" << test_ind << "\n---------------------------------\n";

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			unsigned long long before, after;
			before = GF2_auxiliary_storage::operation_number;

			int max_common_psm_tmp = psm_red_ed_cyc.get_MRIP_sys_G(rand_sequence, G_psm, permute_psm);
			max_common_psm.push_back(max_common_psm_tmp);

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl; 
			computation_cost_psm.push_back(after - before);

			//cout << "G_target" << G_target;
			G_psm.permute_col_back(permute_psm);
			if (G_psm.multiply_transpose_of(H).isZero() == false) {
				cout << "(psm_red_recur_part) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted = permute_psm.get_part(0, 0, 0, k - 1);
			info_set_sorted.sort('<');

			// count the computation cost of Gaussian elimination

			before = GF2_auxiliary_storage::operation_number;
			Matrix<GF2> G_GE = G;
			G_GE.permute_col(rand_sequence);

			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			//Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2();		// new method
			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2_echelon();		// new method

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
			computation_cost_GE.push_back(after - before);

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
			//cout << "permute_GE_real" << permute_GE_real;
			//cout << "G_GE" << G_GE;
			G_GE.permute_col_back(permute_GE_real);
			if (G_GE.multiply_transpose_of(H).isZero() == false) {
				cout << "(GE) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
			info_set_sorted_GE.sort('<');

			//cout << "-----------" << endl;

			// compare the info_set

			//cout << "info_set_sorted - info_set_sorted_GE" << info_set_sorted - info_set_sorted_GE << endl;
			if ((info_set_sorted - info_set_sorted_GE).isZero() == false) {
				cout << "(psm_red_recur_part) not a MRB, test_ind = " << test_ind << endl;
			}
		}

		//cout << "-----------" << endl;
		//cout << "max_common_psm" << max_common_psm;

		//cout << "-----------" << endl;
		//cout << "computation_cost_psm (equal dist cyc)" << computation_cost_psm;
		//cout << "computation_cost_GE (equal dist cyc)" << computation_cost_GE;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_ave (equal dist cyc) = " << Matrix_common::ave(computation_cost_psm) << endl;
		cout << "computation_cost_GE_ave (equal dist cyc) = " << Matrix_common::ave(computation_cost_GE) << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_max (equal dist cyc) = " << Matrix_common::max(computation_cost_psm) << endl;
		cout << "computation_cost_GE_max (equal dist cyc)  = " << Matrix_common::max(computation_cost_GE) << endl;
	}

	// extend shift property to fit extended cyclic code. no multiplication property used.
	static void PreStored_Matrix_red_equal_dist_extend_cycle_test_ini() {
		const int m = 4;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 3;
		eBCH<m, t> ebch;
		ebch.print_info();
		int n = ebch.get_n();
		int k = ebch.get_k();

		Matrix<GF2> G = ebch.get_generator_matrix();
		Matrix<GF2> H = ebch.get_parity_matrix();
		bool is_generator = true;

		//Matrix<GF2> G = ebch.get_parity_matrix();		// inspect on dual code
		//Matrix<GF2> H = ebch.get_generator_matrix();
		//bool is_generator = false;

		cout << "G" << G;

		int  partition_pattern_num = 2;
		PreStored_Matrix_red_equal_dist_extend_cycle psm_red_ed_ext_cyc_G(G, is_generator, partition_pattern_num);
	}
	static void PreStored_Matrix_red_equal_dist_extend_cycle_test_multi_trial() {
		const int m = 6;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 5;
		eBCH<m, t> ebch;
		ebch.print_info();

		/*Matrix<GF2> G = ebch.get_generator_matrix();
		Matrix<GF2> H = ebch.get_parity_matrix();
		bool is_generator = true;*/

		Matrix<GF2> G = ebch.get_parity_matrix();		// inspect on dual code
		Matrix<GF2> H = ebch.get_generator_matrix();
		bool is_generator = false;

		//cout << "G" << G;

		int n = ebch.get_n();
		int k = G.row();								// no problem

		int partition_pattern_num = 8;
		PreStored_Matrix_red_equal_dist_extend_cycle psm_red_ed_ext_cyc(G, is_generator, partition_pattern_num);
		// less than 1 percent complexity different for setting the correct is_generator

		int test_total = 400;
		Matrix<my_double> computation_cost_psm(1, test_total, 'v');
		Matrix<my_double> computation_cost_GE(1, test_total, 'v');

		Matrix<GF2> G_psm(k, n);
		Matrix<int> permute_psm(1, n);		// the result from 'PreStored_Matrix_red_loose_bin'
		Matrix<int> max_common_psm(1, test_total, 'v');

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			//cout << "------------------------------------\n" << "test_ind=" << test_ind << "\n---------------------------------\n";

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			//cout << "rand_sequence" << rand_sequence << endl;

			unsigned long long before, after;
			before = GF2_auxiliary_storage::operation_number;

			int max_common_psm_tmp = psm_red_ed_ext_cyc.get_MRIP_sys_G(rand_sequence, G_psm, permute_psm);
			max_common_psm.push_back(max_common_psm_tmp);

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl; 
			computation_cost_psm.push_back(after - before);

			//cout << "G_target" << G_target;
			G_psm.permute_col_back(permute_psm);
			if (G_psm.multiply_transpose_of(H).isZero() == false) {
				cout << "(psm_red_recur_part) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted = permute_psm.get_part(0, 0, 0, k - 1);
			info_set_sorted.sort('<');

			// count the computation cost of Gaussian elimination

			before = GF2_auxiliary_storage::operation_number;
			Matrix<GF2> G_GE = G;
			G_GE.permute_col(rand_sequence);

			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2();		// new method

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
			computation_cost_GE.push_back(after - before);

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
			//cout << "permute_GE_real" << permute_GE_real;
			//cout << "G_GE" << G_GE;
			G_GE.permute_col_back(permute_GE_real);
			if (G_GE.multiply_transpose_of(H).isZero() == false) {
				cout << "(GE) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
			info_set_sorted_GE.sort('<');

			//cout << "-----------" << endl;

			// compare the info_set

			//cout << "info_set_sorted" << info_set_sorted;
			//cout << "info_set_sorted_GE" << info_set_sorted_GE;

			if ((info_set_sorted - info_set_sorted_GE).isZero() == false) {
				cout << "(psm_red_recur_part) not a MRB, test_ind = " << test_ind << endl;
			}
		}

		//cout << "-----------" << endl;
		//cout << "max_common_psm" << max_common_psm;

		//cout << "-----------" << endl;
		//cout << "computation_cost_psm (equal dist cyc)" << computation_cost_psm;
		//cout << "computation_cost_GE (equal dist cyc)" << computation_cost_GE;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_ave (equal dist cyc) = " << Matrix_common::ave(computation_cost_psm) << endl;
		cout << "computation_cost_GE_ave (equal dist cyc) = " << Matrix_common::ave(computation_cost_GE) << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_max (equal dist cyc) = " << Matrix_common::max(computation_cost_psm) << endl;
		cout << "computation_cost_GE_max (equal dist cyc)  = " << Matrix_common::max(computation_cost_GE) << endl;
	}

	// IBUc in the paper, both shift and multiplication are used
	static void PreStored_Matrix_red_cycle_ini_test() {
		const int m = 4;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 1;
		BCH<m, t> bch;
		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();
		int red = n - k;

		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();

		PreStored_Matrix_red_cycle psm_red_cyc(H);

		Matrix<int> reliability_rand(1, n, 'N');
		reliability_rand.permute_rand();
		cout << "reliability_rand" << reliability_rand;

		psm_red_cyc.select_G_set_continous_count(reliability_rand);
	}
	static void PreStored_Matrix_red_cycle_test_multi_trial() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 10;
		BCH<m, t> bch;
		bch.print_info();

		//Matrix<GF2> G = bch.get_generator_matrix();
		//Matrix<GF2> H = bch.get_parity_matrix();

		Matrix<GF2> G = bch.get_parity_matrix();		// inspect on dual code
		Matrix<GF2> H = bch.get_generator_matrix();

		/*cout << "G" << G;
		cout << "H" << H;*/

		int n = bch.get_n();
		int k = G.row();								// no problem

		PreStored_Matrix_red_cycle psm_red_cyc(G);

		int test_total = 400000;
		//Matrix<my_double> computation_cost_psm(1, test_total, 'v');
		//Matrix<my_double> computation_cost_GE(1, test_total, 'v');

		my_double sum_ope_num_psm = 0;
		my_double sum_ope_num_GE = 0;
		my_double max_ope_num_psm = 0;
		my_double max_ope_num_GE = 0;

		Matrix<GF2> G_psm(k, n);
		Matrix<int> permute_psm(1, k);		// the result from 'PreStored_Matrix_red_loose_bin'
		//Matrix<int> max_common_psm(1, test_total, 'v');

		//my::set_seed_adv(7);

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			//cout << "------------------------------------\n" << "test_ind=" << test_ind << "\n---------------------------------\n";

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			//rand_sequence = Matrix<int>(1, n, { 2,  0,  5,  1,  6,  3,  4 });

			unsigned long long before, after, ope_num_tmp;
			before = GF2_auxiliary_storage::operation_number;

			int max_common_psm_tmp = psm_red_cyc.get_MRIP_sys_G(rand_sequence, G_psm, permute_psm);
			//max_common_psm.push_back(max_common_psm_tmp);

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl;
			//computation_cost_psm.push_back(after - before);
			ope_num_tmp = after - before;
			sum_ope_num_psm += ope_num_tmp;
			max_ope_num_psm = max_ope_num_psm < ope_num_tmp ? ope_num_tmp : max_ope_num_psm;

			//cout << "psm_red_cyc.common_pos_num" << psm_red_cyc.common_pos_num;
			//cout << "psm_red_cyc.common_pos_max_shift" << psm_red_cyc.common_pos_max_shift;

			//cout << "G_target" << G_target;
			G_psm.permute_col_back(permute_psm);
			if (G_psm.multiply_transpose_of(H).isZero() == false) {
				cout << "***************" << endl;
				cout << "(psm_red_cyc) unexpected result, test_ind = " << test_ind << endl;

				cout << "psm_red_cyc.common_pos_num" << psm_red_cyc.common_pos_num;
				cout << "psm_red_cyc.common_pos_max_shift" << psm_red_cyc.common_pos_max_shift;

				cout << "rand_sequence" << rand_sequence;

				//cout << "computation_cost_psm.back() = " << computation_cost_psm.back() << endl;
			}
			Matrix<int> info_set_sorted = permute_psm.get_part(0, 0, 0, k - 1);
			info_set_sorted.sort('<');

			// count the computation cost of Gaussian elimination

			before = GF2_auxiliary_storage::operation_number;
			Matrix<GF2> G_GE = G;
			G_GE.permute_col(rand_sequence);

			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			//Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2();		// new method
			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2_echelon();		// new method

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
			//computation_cost_GE.push_back(after - before);
			ope_num_tmp = after - before;
			sum_ope_num_GE += ope_num_tmp;
			max_ope_num_GE = max_ope_num_GE < ope_num_tmp ? ope_num_tmp : max_ope_num_GE;

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
			//cout << "permute_GE_real" << permute_GE_real;
			//cout << "G_GE" << G_GE;
			G_GE.permute_col_back(permute_GE_real);
			if (G_GE.multiply_transpose_of(H).isZero() == false) {
				cout << "(GE) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
			info_set_sorted_GE.sort('<');

			//cout << "-----------" << endl;

			// compare the info_set

			//cout << "info_set_sorted - info_set_sorted_GE" << info_set_sorted - info_set_sorted_GE << endl;
			if ((info_set_sorted - info_set_sorted_GE).isZero() == false) {
				cout << "(psm_red_cyc) not a MRB, test_ind = " << test_ind << endl;
			}
		}

		//cout << "-----------" << endl;
		//cout << "max_common_psm" << max_common_psm;

		//cout << "-----------" << endl;
		//cout << "computation_cost_psm (equal dist cyc)" << computation_cost_psm;
		//cout << "computation_cost_GE (equal dist cyc)" << computation_cost_GE;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_ave (equal dist cyc) = " << sum_ope_num_psm / test_total << endl;
		cout << "computation_cost_GE_ave (equal dist cyc) = " << sum_ope_num_GE / test_total << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_psm_max (equal dist cyc) = " << max_ope_num_psm << endl;
		cout << "computation_cost_GE_max (equal dist cyc)  = " << max_ope_num_GE << endl;
	}

	// general test for the GE, the IBU and the IBUc
	static void GE_IBU_IBUc_test_all() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 10;
		BCH<m, t> bch;
		bch.print_info();

		Matrix<GF2> G, H;
		my_double rate = bch.get_k() / (my_double)bch.get_n();
		if (false /*rate > 0.5*/) {
			G = bch.get_parity_matrix();		// inspect on dual code
			H = bch.get_generator_matrix();
		}
		else {
			G = bch.get_generator_matrix();		// inspect on original code
			H = bch.get_parity_matrix();
		}		

		/*cout << "G" << G;
		cout << "H" << H;*/

		int n = bch.get_n();
		int k = G.row();								// no problem

		int partition_pattern_num = 1;
		PreStored_Matrix_red_equal_dist IBU(G, partition_pattern_num);		// the IBU

		PreStored_Matrix_red_cycle IBUc(G);		// the IBUc

		int test_total = 100000;
#if complexity_distribution
		Matrix<my_double> computation_cost_GE(1, test_total, 'v');
		Matrix<my_double> computation_cost_IBU(1, test_total, 'v');
		Matrix<my_double> computation_cost_IBUc(1, test_total, 'v');
#endif

		my_double sum_ope_num_IBU = 0;
		my_double sum_ope_num_IBUc = 0;
		my_double sum_ope_num_GE = 0;
		my_double max_ope_num_IBU = 0;
		my_double max_ope_num_IBUc = 0;
		my_double max_ope_num_GE = 0;

		Matrix<GF2> G_IBU(k, n);
		Matrix<int> permute_IBU(1, k);			// the result of IBU

		Matrix<GF2> G_IBUc(k, n);
		Matrix<int> permute_IBUc(1, k);			// the result of IBUc

		int worst_case_common_num_IBU = 0;
		int worst_case_common_num_IBUc = 0;		// anything is okay
		int min_common_num_IBU_all_test = n;
		int min_common_num_IBUc_all_test = n;	// must be an upper bound
		int update_num_IBU_sum = 0;
		int update_num_IBUc_sum = 0;


		//my::set_seed_adv(7);

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			//cout << "------------------------------------\n" << "test_ind=" << test_ind << "\n---------------------------------\n";


			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();

			unsigned long long before, after, ope_num_tmp;
			before = GF2_auxiliary_storage::operation_number;

			int max_common_num_IBU = IBU.get_MRIP_sys_G(rand_sequence, G_IBU, permute_IBU);		// incorrect parameter, discarded
			update_num_IBU_sum += k - max_common_num_IBU;

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl; 
			ope_num_tmp = after - before;
#if complexity_distribution
			computation_cost_IBU.push_back(ope_num_tmp);
#endif
			sum_ope_num_IBU += ope_num_tmp;
			worst_case_common_num_IBU = max_ope_num_IBU < ope_num_tmp ? max_common_num_IBU : worst_case_common_num_IBU;
			min_common_num_IBU_all_test = \
				min_common_num_IBU_all_test < max_common_num_IBU ? min_common_num_IBU_all_test: max_common_num_IBU;
			max_ope_num_IBU = max_ope_num_IBU < ope_num_tmp ? ope_num_tmp : max_ope_num_IBU;

			//cout << "G_target" << G_target;
			G_IBU.permute_col_back(permute_IBU);
			if (G_IBU.multiply_transpose_of(H).isZero() == false) {
				cout << "(IBU) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_IBU = permute_IBU.get_part(0, 0, 0, k - 1);
			//cout << "info_set_sorted (before sorted)" << info_set_sorted;
			info_set_sorted_IBU.sort('<');


			//rand_sequence = Matrix<int>(1, n, { 2,  0,  5,  1,  6,  3,  4 });

			before = GF2_auxiliary_storage::operation_number;

			int max_common_num_IBUc = IBUc.get_MRIP_sys_G(rand_sequence, G_IBUc, permute_IBUc);
			update_num_IBUc_sum += k - max_common_num_IBUc;
			//max_common_psm.push_back(max_common_psm_IBUc);

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (psm)(orig) = " << after - before << " ---------------------" << endl;
			ope_num_tmp = after - before;
#if complexity_distribution
			computation_cost_IBUc.push_back(ope_num_tmp);
#endif
			sum_ope_num_IBUc += ope_num_tmp;
			worst_case_common_num_IBUc = max_ope_num_IBUc < ope_num_tmp ? max_common_num_IBUc : worst_case_common_num_IBUc; 
			min_common_num_IBUc_all_test = \
				min_common_num_IBUc_all_test < max_common_num_IBUc ? min_common_num_IBUc_all_test: max_common_num_IBUc;		// find the min
			max_ope_num_IBUc = max_ope_num_IBUc < ope_num_tmp ? ope_num_tmp : max_ope_num_IBUc;

			//cout << "IBUc.common_pos_num" << IBUc.common_pos_num;
			//cout << "IBUc.common_pos_max_shift" << IBUc.common_pos_max_shift;

			//cout << "G_target" << G_target;
			G_IBUc.permute_col_back(permute_IBUc);
			if (G_IBUc.multiply_transpose_of(H).isZero() == false) {
				cout << "***************" << endl;
				cout << "(IBUc) unexpected result, test_ind = " << test_ind << endl;

				cout << "IBUc.common_pos_num" << IBUc.common_pos_num;
				cout << "IBUc.common_pos_max_shift" << IBUc.common_pos_max_shift;

				cout << "rand_sequence" << rand_sequence;

				//cout << "computation_cost_psm.back() = " << computation_cost_psm.back() << endl;
			}
			Matrix<int> info_set_sorted_IBUc = permute_IBUc.get_part(0, 0, 0, k - 1);
			info_set_sorted_IBUc.sort('<');

			// count the computation cost of Gaussian elimination

			before = GF2_auxiliary_storage::operation_number;
			Matrix<GF2> G_GE = G;
			G_GE.permute_col(rand_sequence);

			/*G_GE.row_transformation_to_up_triangle();
			Matrix<int> permute_GE = G_GE.col_permute_to_full_rank_on_left();
			G_GE.row_transformation_left_up_triangle_to_identity();*/

			//Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2();		// new method
			Matrix<int> permute_GE = G_GE.GE_left_identity_4_GF2_echelon();		// new method

			after = GF2_auxiliary_storage::operation_number;
			//cout << "------------------ computation cost (GE) = " << after - before << " ---------------------" << endl;
			ope_num_tmp = after - before;
#if complexity_distribution
			computation_cost_GE.push_back(ope_num_tmp);
#endif
			sum_ope_num_GE += ope_num_tmp;
			max_ope_num_GE = max_ope_num_GE < ope_num_tmp ? ope_num_tmp : max_ope_num_GE;

			Matrix<int> permute_GE_real = rand_sequence;
			permute_GE_real.permute(permute_GE);
			//cout << "permute_GE_real" << permute_GE_real;
			//cout << "G_GE" << G_GE;

			G_GE.permute_col_back(permute_GE_real);
			if (G_GE.multiply_transpose_of(H).isZero() == false) {
				cout << "(GE) unexpected result, test_ind = " << test_ind << endl;
			}
			Matrix<int> info_set_sorted_GE = permute_GE_real.get_part(0, 0, 0, k - 1);
			info_set_sorted_GE.sort('<');

			//cout << "-----------" << endl;

			// compare the info_set

			//cout << "info_set_sorted - info_set_sorted_GE" << info_set_sorted - info_set_sorted_GE << endl;
			if ((info_set_sorted_IBU - info_set_sorted_GE).isZero() == false) {
				cout << "(IBU) not a MRB, test_ind = " << test_ind << endl;
			}
			if ((info_set_sorted_IBUc - info_set_sorted_GE).isZero() == false) {
				cout << "(IBUc) not a MRB, test_ind = " << test_ind << endl;
			}
		}

		double ave_case_update_num_IBU = update_num_IBU_sum / (double)test_total;
		double ave_case_update_num_IBUc = update_num_IBUc_sum / (double)test_total;

		int worst_case_update_num_IBU = k - worst_case_common_num_IBU;
		int worst_case_update_num_IBUc = k - worst_case_common_num_IBUc;

		int max_update_num_IBU_all_test = k - min_common_num_IBU_all_test;
		int max_update_num_IBUc_all_test = k - min_common_num_IBUc_all_test;

		//cout << "-----------" << endl;
		//cout << "max_common_psm" << max_common_psm;

		//cout << "-----------" << endl;
		//cout << "computation_cost_psm (equal dist cyc)" << computation_cost_psm;
		//cout << "computation_cost_GE (equal dist cyc)" << computation_cost_GE;

		cout << "-----------" << endl;
		cout << "computation_cost_GE_ave = " << sum_ope_num_GE / test_total << endl;
		cout << "computation_cost_GE_max = " << max_ope_num_GE << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_IBU_ave = " << sum_ope_num_IBU / test_total << endl;
		cout << "computation_cost_IBU_max = " << max_ope_num_IBU << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_IBUc_ave = " << sum_ope_num_IBUc / test_total << endl;
		cout << "computation_cost_IBUc_max = " << max_ope_num_IBUc << endl;


		cout << "-----------" << endl;
		cout << "ave_case_update_num_IBU = " << ave_case_update_num_IBU << endl;
		cout << "worst_case_update_num_IBU = " << worst_case_update_num_IBU << endl;
		cout << "max_update_num_IBU_all_test = " << max_update_num_IBU_all_test << endl;

		cout << "-----------" << endl;
		cout << "ave_case_update_num_IBUc = " << ave_case_update_num_IBUc << endl;
		cout << "worst_case_update_num_IBUc = " << worst_case_update_num_IBUc << endl;
		cout << "max_update_num_IBUc_all_test = " << max_update_num_IBUc_all_test << endl;


		cout << "-----------" << endl;

		// output to file
		
		// set cout into file
		char file_name[550] = { 0 };		// set it large enough
		sprintf_s(file_name, 550, "test_on_GE_IBU_IBUc_n_%d_k_%d.txt", bch.get_n(), bch.get_k());
		FILE* stream1;
		freopen_s(&stream1, file_name, "w", stdout);

		bch.print_info();

		cout << "-----------" << endl;
		cout << "computation_cost_GE_ave = " << sum_ope_num_GE / test_total << endl;
		cout << "computation_cost_GE_max = " << max_ope_num_GE << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_IBU_ave = " << sum_ope_num_IBU / test_total << endl;
		cout << "computation_cost_IBU_max = " << max_ope_num_IBU << endl;

		cout << "-----------" << endl;
		cout << "computation_cost_IBUc_ave = " << sum_ope_num_IBUc / test_total << endl;
		cout << "computation_cost_IBUc_max = " << max_ope_num_IBUc << endl;


		cout << "-----------" << endl;
		cout << "ave_case_update_num_IBU = " << ave_case_update_num_IBU << endl;
		cout << "worst_case_update_num_IBU = " << worst_case_update_num_IBU << endl;
		cout << "max_update_num_IBU_all_test = " << max_update_num_IBU_all_test << endl;

		cout << "-----------" << endl;
		cout << "ave_case_update_num_IBUc = " << ave_case_update_num_IBUc << endl;
		cout << "worst_case_update_num_IBUc = " << worst_case_update_num_IBUc << endl;
		cout << "max_update_num_IBUc_all_test = " << max_update_num_IBUc_all_test << endl;

		cout << "-----------" << endl;

		fclose(stdout);

		system("pause");

#if complexity_distribution
		ofstream GE_complexity_distribution("GE_complexity_distribution.txt");
		for (int i = 0; i < test_total; ++i) {
			GE_complexity_distribution << computation_cost_GE(i) << '\t';
		}
		ofstream IBU_complexity_distribution("IBU_complexity_distribution.txt");
		for (int i = 0; i < test_total; ++i) {
			IBU_complexity_distribution << computation_cost_IBU(i) << '\t';
		}
		ofstream IBUc_complexity_distribution("IBUc_complexity_distribution.txt");
		for (int i = 0; i < test_total; ++i) {
			IBUc_complexity_distribution << computation_cost_IBUc(i) << '\t';
		}
#endif
	}
	static void GJE_test() {
		const int m = 7;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 10;
		BCH<m, t> bch;
		bch.print_info();
		int n = bch.get_n();

		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<int> rand_n(1, n, 'N');

		unsigned long long GF2_ope_before, GF2_ope_after;
		unsigned long long GE_GF2_ope_sum = 0, GE_GF2_ope_wst = 0, GJE_GF2_ope_sum = 0, GJE_GF2_ope_wst = 0;
		int GJE_column_proceed_sum = 0, GJE_column_proceed_wst = 0;

		int test_times = 100000;
		for (int i = 0; i < test_times; ++i) {

			//cout << "i = " << i << endl;
			rand_n.permute_rand();
			G.permute_col(rand_n);

			/////////////////////////////////////////////////////////////////////////////

			GF2_ope_before = GF2_auxiliary_storage::operation_number;

			Matrix<GF2> Gs_GE = G;
			Matrix<int> permute_GE = Gs_GE.GE_left_identity_4_GF2_echelon();

			GF2_ope_after = GF2_auxiliary_storage::operation_number;
			unsigned long long GE_GF2_ope = GF2_ope_after - GF2_ope_before;

			GE_GF2_ope_sum += GE_GF2_ope;
			GE_GF2_ope_wst = GE_GF2_ope_wst < GE_GF2_ope ? GE_GF2_ope : GE_GF2_ope_wst;

			/////////////////////////////////////////////////////////////////////////////

			GF2_ope_before = GF2_auxiliary_storage::operation_number;

			Matrix<GF2> Gs_GJE = G;
			int GJE_column_proceed;
			Matrix<int> permute_GJE = Gs_GJE.GJE_left_identity_4_GF2(GJE_column_proceed);

			GF2_ope_after = GF2_auxiliary_storage::operation_number;
			unsigned long long GJE_GF2_ope = GF2_ope_after - GF2_ope_before;

			GJE_GF2_ope_sum += GJE_GF2_ope;
			GJE_GF2_ope_wst = GJE_GF2_ope_wst < GJE_GF2_ope ? GJE_GF2_ope : GJE_GF2_ope_wst;
			GJE_column_proceed_sum += GJE_column_proceed;
			GJE_column_proceed_wst = GJE_column_proceed_wst < GJE_column_proceed ? GJE_column_proceed : GJE_column_proceed_wst;


			/////////////////////////////////////////////////////////////////////////////

			if ((Gs_GE - Gs_GJE).isZero() && (permute_GE - permute_GJE).isZero());
			else {
				cout << "GJE in error" << endl;
			}
		}
		cout << "GE_GF2_ope_ave = " << GE_GF2_ope_sum / (double)test_times << endl;
		cout << "GE_GF2_ope_wst = " << GE_GF2_ope_wst << endl;

		cout << "GJE_GF2_ope_ave = " << GJE_GF2_ope_sum / (double)test_times << endl;
		cout << "GJE_GF2_ope_wst = " << GJE_GF2_ope_wst << endl;

		cout << "GJE_column_proceed_ave = " << GJE_column_proceed_sum / (double)test_times << endl;
		cout << "GJE_column_proceed_wst = " << GJE_column_proceed_wst << endl;
	}
	static void GE_IBU_IBUc_test_all_2() {
		const int m = 5;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 7;
		BCH<m, t> bch;
		bch.print_info();
		Matrix<GF2> G = bch.get_generator_matrix();							// only test on Generator matrix
		int n = bch.get_n();
		int k = G.row();

		int partition_pattern_num = 1;
		PreStored_Matrix_red_equal_dist IBU(G, partition_pattern_num);		// the IBU
		PreStored_Matrix_red_cycle IBUc(G);									// the IBUc

		int test_total = 100000;

#if complexity_distribution
		Matrix<my_double> computation_cost_GE(1, test_total, 'v');
		Matrix<my_double> computation_cost_IBU(1, test_total, 'v');
		Matrix<my_double> computation_cost_IBUc(1, test_total, 'v');
#endif

		unsigned long long GJE_sum_seq_steps = 0, GJE_wst_seq_steps = 0, GJE_sum_ope_num = 0, GJE_wst_ope_num = 0;
		unsigned long long IBU_sum_seq_steps = 0, IBU_wst_seq_steps = 0, IBU_sum_ope_num = 0, IBU_wst_ope_num = 0;
		unsigned long long IBUc_sum_seq_steps = 0, IBUc_wst_seq_steps = 0, IBUc_sum_ope_num = 0, IBUc_wst_ope_num = 0;

		Matrix<GF2> G_GE(k, n); 
		Matrix<int> permute_GE(1, n);			// the result of IBU

		Matrix<GF2> G_IBU(k, n);
		Matrix<int> permute_IBU(1, n);			// the result of IBU

		Matrix<GF2> G_IBUc(k, n);
		Matrix<int> permute_IBUc(1, n);			// the result of IBUc


		//my::set_seed_adv(7);

		for (int test_ind = 0; test_ind < test_total; ++test_ind) {

			///////////////////////////////////////////////////////////////////////////// setup

			Matrix<int> rand_sequence(1, n, 'N');
			rand_sequence.permute_rand();
			unsigned long long before, after, GJE_ope_num, IBU_ope_num, IBUc_ope_num;

			///////////////////////////////////////////////////////////////////////////// GE


			before = GF2_auxiliary_storage::operation_number;
			G_GE = G;
			G_GE.permute_col(rand_sequence);
			permute_GE = rand_sequence;
			int GJE_seq_steps;
			permute_GE.permute(G_GE.GJE_left_identity_4_GF2(GJE_seq_steps));
			after = GF2_auxiliary_storage::operation_number;
			GJE_ope_num = after - before;

			GJE_sum_seq_steps += GJE_seq_steps;
			GJE_wst_seq_steps = my::max((unsigned long long) GJE_seq_steps, GJE_wst_seq_steps);
			GJE_sum_ope_num += GJE_ope_num;
			GJE_wst_ope_num = my::max(GJE_ope_num, GJE_wst_ope_num);


#if complexity_distribution
			computation_cost_GE.push_back(ope_num_tmp);
#endif

			///////////////////////////////////////////////////////////////////////////// IBU
			
			before = GF2_auxiliary_storage::operation_number;
			unsigned long long IBU_seq_steps = IBU.get_MRIP_sys_G(rand_sequence, G_IBU, permute_IBU);
			after = GF2_auxiliary_storage::operation_number;
			IBU_ope_num = after - before;

			IBU_sum_seq_steps += IBU_seq_steps;
			IBU_wst_seq_steps = my::max(IBU_seq_steps, IBU_wst_seq_steps);
			IBU_sum_ope_num += IBU_ope_num;
			IBU_wst_ope_num = my::max(IBU_ope_num, IBU_wst_ope_num);

#if complexity_distribution
			computation_cost_IBU.push_back(ope_num_tmp);
#endif

			///////////////////////////////////////////////////////////////////////////// IBUc

			before = GF2_auxiliary_storage::operation_number;
			unsigned long long IBUc_seq_steps = IBUc.get_MRIP_sys_G(rand_sequence, G_IBUc, permute_IBUc);		
			after = GF2_auxiliary_storage::operation_number;
			IBUc_ope_num = after - before;

			IBUc_sum_seq_steps += IBUc_seq_steps;
			IBUc_wst_seq_steps = my::max(IBUc_seq_steps, IBUc_wst_seq_steps);
			IBUc_sum_ope_num += IBUc_ope_num;
			IBUc_wst_ope_num = my::max(IBUc_ope_num, IBUc_wst_ope_num);

#if complexity_distribution
			computation_cost_IBUc.push_back(ope_num_tmp);
#endif

			/////////////////////////////////////////////////////////////////////////////  compare the info_set

			Matrix<int> info_set_sorted_GE = permute_GE.get_part(0, 0, 0, k - 1);
			Matrix<int> info_set_sorted_IBU = permute_IBU.get_part(0, 0, 0, k - 1);
			Matrix<int> info_set_sorted_IBUc = permute_IBUc.get_part(0, 0, 0, k - 1);
			info_set_sorted_GE.sort('<');
			info_set_sorted_IBU.sort('<');
			info_set_sorted_IBUc.sort('<');
			//cout << "info_set_sorted - info_set_sorted_GE" << info_set_sorted - info_set_sorted_GE << endl;
			if ((info_set_sorted_IBU - info_set_sorted_GE).isZero() == false) {
				cout << "(IBU) not obtaining a MRB, test_ind = " << test_ind << endl;
			}
			if ((info_set_sorted_IBUc - info_set_sorted_GE).isZero() == false) {
				cout << "(IBUc) not obtaining a MRB, test_ind = " << test_ind << endl;
			}
		}

		double GJE_ave_seq_steps = GJE_sum_seq_steps / (double)test_total;
		double GJE_ave_ope_num = GJE_sum_ope_num / (double)test_total;

		double IBU_ave_seq_steps = IBU_sum_seq_steps / (double)test_total;
		double IBU_ave_ope_num = IBU_sum_ope_num / (double)test_total;

		double IBUc_ave_seq_steps = IBUc_sum_seq_steps / (double)test_total;
		double IBUc_ave_ope_num = IBUc_sum_ope_num / (double)test_total;


		cout << "-----------" << endl;
		cout << "GJE_ave_seq_steps = " << GJE_ave_seq_steps << endl;
		cout << "IBU_ave_seq_steps = " << IBU_ave_seq_steps << endl;
		cout << "IBUc_ave_seq_steps = " << IBUc_ave_seq_steps << endl;
		cout << "GJE_wst_seq_steps = " << GJE_wst_seq_steps << endl;
		cout << "IBU_wst_seq_steps = " << IBU_wst_seq_steps << endl;
		cout << "IBUc_wst_seq_steps = " << IBUc_wst_seq_steps << endl;

		cout << "-----------" << endl;
		cout << "GJE_ave_ope_num = " << GJE_ave_ope_num << endl;
		cout << "IBU_ave_ope_num = " << IBU_ave_ope_num << endl;
		cout << "IBUc_ave_ope_num = " << IBUc_ave_ope_num << endl;
		cout << "GJE_wst_ope_num = " << GJE_wst_ope_num << endl;
		cout << "IBU_wst_ope_num = " << IBU_wst_ope_num << endl;
		cout << "IBUc_wst_ope_num = " << IBUc_wst_ope_num << endl;

		cout << "-----------" << endl;

		// output to file

		// set cout into file
		char file_name[550] = { 0 };		// set it large enough
		sprintf_s(file_name, 550, "test_on_GE_IBU_IBUc_n_%d_k_%d.txt", bch.get_n(), bch.get_k());
		FILE* stream1;
		freopen_s(&stream1, file_name, "w", stdout);

		bch.print_info();

		cout << "-----------" << endl;
		cout << "GJE_ave_seq_steps = " << GJE_ave_seq_steps << endl;
		cout << "GJE_wst_seq_steps = " << GJE_wst_seq_steps << endl;
		cout << "IBU_ave_seq_steps = " << IBU_ave_seq_steps << endl;
		cout << "IBU_wst_seq_steps = " << IBU_wst_seq_steps << endl;
		cout << "IBUc_ave_seq_steps = " << IBUc_ave_seq_steps << endl;
		cout << "IBUc_wst_seq_steps = " << IBUc_wst_seq_steps << endl;

		cout << "-----------" << endl;
		cout << "GJE_ave_ope_num = " << GJE_ave_ope_num << endl;
		cout << "GJE_wst_ope_num = " << GJE_wst_ope_num << endl;
		cout << "IBU_ave_ope_num = " << IBU_ave_ope_num << endl;
		cout << "IBU_wst_ope_num = " << IBU_wst_ope_num << endl;
		cout << "IBUc_ave_ope_num = " << IBUc_ave_ope_num << endl;
		cout << "IBUc_wst_ope_num = " << IBUc_wst_ope_num << endl;

		cout << "-----------" << endl;

		fclose(stdout);

#if complexity_distribution
		ofstream GE_complexity_distribution("GE_complexity_distribution.txt");
		for (int i = 0; i < test_total; ++i) {
			GE_complexity_distribution << computation_cost_GE(i) << '\t';
		}
		ofstream IBU_complexity_distribution("IBU_complexity_distribution.txt");
		for (int i = 0; i < test_total; ++i) {
			IBU_complexity_distribution << computation_cost_IBU(i) << '\t';
		}
		ofstream IBUc_complexity_distribution("IBUc_complexity_distribution.txt");
		for (int i = 0; i < test_total; ++i) {
			IBUc_complexity_distribution << computation_cost_IBUc(i) << '\t';
		}
#endif

		system("pause");
	}

	static void OSD_v2_simulation(bool& is_first_time, my_double SNR_dB = 2.0) {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 5000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 4;
		const int t = 2;
		GF2e<m>::init();
		BCH<m, t> bch;			// here we choose eBCH
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> PM = bch.get_parity_matrix();
		//bch.print_info();

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//cout << "v" << v;

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		int order = 2;
		OSD_v2 osd(PM, false, d_min);		// have problem in high SNR

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;

		if (is_first_time) {

			bch.print_info();
			cout << "max_error_frame          = " << max_error_frame << endl;
			cout << "simulation type          = OSD-" << order << endl;

			cout << "\n------------------------------------------------------------------------------------------------------------"
				<< "--------------------------------------------------------------------------------" << endl;


			cout << setw(7) << "SNR_dB" << setw(20) << "simulation_times" << setw(20) << "error_rate" \
				<< setw(20) << "ave_used_list" << setw(20) << "early_stop_rate";

#ifdef count_operation_number
#ifdef use_my_double

			cout << setw(20) << "double_ope_num";
#endif // use_my_double

			cout << setw(20) << "GF2_ope_num" << setw(20) << "GF2e_ope_num";
#endif // count_operation_number

			cout << setw(20) << "time_consume(s)" << endl;

			is_first_time = false;
		}

		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
		unsigned long long double_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;

		unsigned long long GF2e_ope_num_before = GF2e_auxiliary_storage::operation_number;
		unsigned long long GF2e_ope_num_after;
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> recv = AWGN::pass_standard(c);
			//Matrix<my_double> recv = AWGN::pass(c);
			Matrix<GF2> v_hat = osd.decode_v(recv, order);

			if (v == v_hat);
			else {
				// record the error frame
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
			//cout << "simulation_times = " << i + 1 << '\r';
		}
		end = clock();

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << setw(7) << fixed << setprecision(1) << SNR_dB;
		cout << setw(20) << simulation_times;
		cout << scientific << setprecision(2);
		cout << setw(20) << error_frame / (double)simulation_times;

		cout << setw(20) << osd.total_used_list_num / (double)simulation_times;
		cout << setw(20) << osd.num_early_stop / (double)simulation_times;

#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << setw(20) << (double_ope_num_after - double_ope_num_before) / (double)simulation_times;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << setw(20) << (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << setw(20) << (GF2e_ope_num_after - GF2e_ope_num_before) / (double)simulation_times;
#endif // count_operation_number

		cout << setw(20) << time_consume - pass_AWGN_time << endl;

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);

	}
	static void OSD_v2_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "OSD_r_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(0, 1, 7, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			OSD_v2_simulation(is_first_time, test_SNR(i));
		}
	}
	static void OSD_v2_fixed_received_vector() {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;

		// encoding parameter
		const int m = 4;
		const int t = 2;
		GF2e<m>::init();
		BCH<m, t> bch;			// here we choose eBCH
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> GM = bch.get_generator_matrix();
		Matrix<GF2> PM = bch.get_parity_matrix();
		bch.print_info();
		cout << "GM" << GM;

		Matrix<GF2> v(1,15,{ 1,   1,   0,   0,   0,   1,   1,   0,   1,   1,   1,   1,   1,   0,   0, });
		cout << "v" << v;

		Matrix<my_double> c = BPSK::modulation(v);

		// channel fixed
		int simulation_times = 1;
		/*
		Matrix<my_double> recv_matrix(simulation_times, n, {
			-2.03026,      -1.42874,      0.291925,     -0.676346,      0.366056,     -0.456418,      -1.59104,       0.57337,      -2.17165,     -0.419984,      0.063066,      0.260627,      -1.62616,       1.70782,      0.798275,
			-1.86435,    -0.0836922,       1.22064,      0.844791,     0.0188393,     -0.351043,      -1.00643,   -0.00488477,      -1.24717,     -0.975994,     0.0511008,      -1.11552,     -0.633307,     -0.251753,      0.443881,
			-1.71704,    -0.0355986,      0.953004,       1.03202,       1.27793,     -0.130623,      -0.81576,       1.57263,      -1.14507,     -0.151133,     -0.500191,         -0.54,      0.550899,     -0.165649,       1.15493,
			-1.61437,     -0.155342,       1.20602,       1.21009,    -0.0847072,      -2.10529,      -1.22242,     -0.218381,      -1.82477,      -1.15878,      0.125826,     -0.588704,      -1.20404,     -0.215511,      0.504907,
			0.355955,      -1.06582,       1.19873,      0.141901,       1.00561,      -1.52328,      -1.58794,       1.02118,     -0.972271,     -0.674702,    -0.0482323,     -0.660635,     -0.645956,     -0.618223,     -0.480048,
			-1.02733,     -0.631098,      0.296092,     -0.300616,       0.13171,      0.379461,      -1.42243,       1.02651,     -0.967755,     -0.764566,      -1.32126,     -0.528388,      0.252366,      0.555595,      0.764855,
			-0.331572,      -1.04161,      0.593228,      0.386934,      -0.61471,      -0.83762,     -0.927975,      0.900959,     -0.158708,      -1.41073,      -1.02084,     -0.420716,      0.917463,       1.77284,      0.771059,
			-1.33356,     -0.529523,       1.61871,     -0.064144,      0.103081,       0.64732,      -1.05777,       1.72363,     -0.993881,     -0.664848,     -0.437472,     -0.994339,     -0.403136,        1.1621,       1.40201,
			-0.360228,      -1.14777,       1.15819,     -0.409536,       1.77654,      -2.01085,      -0.23443,       1.38709,     -0.857713,      0.561468,      -1.03599,      -1.89366,     -0.118135,       1.70687,       1.51544,
			-1.02542,     -0.168361,       1.54907,      0.656995,       1.44681,       0.58474,      -1.70288,      0.158409,     -0.353804,      0.277997,      -0.37322,      -1.05205,      -1.56184,    -0.0357174,       1.38492,
			0.177387,      0.296973,      0.691871,      0.630922,      0.591638,     -0.814624,      -1.43985,      0.920406,     -0.513152,      -1.22324,      -1.52547,      -1.22079,      0.105395,      0.615392,     -0.148018,
			-0.160218,     -0.689221,       2.29261,    -0.0270581,      0.954183,     -0.885226,      -0.99962,     0.0774115,      -1.70996,      -0.76205,      0.317057,      0.658935,      -1.01952,       1.04976,      0.887684,
			-0.735279,      0.400444,      0.110947,      0.681083,       1.07874,    0.00768535,      -1.39625,      0.218291,     -0.188676,      0.190814,     -0.427557,      -1.45897,      -1.57618,       1.31273,       1.25178,
			-0.537445,      -1.68673,       1.29394,      0.855285,      0.758461,      -1.38197,      0.719674,       0.86607,     -0.879969,    -0.0956519,     -0.457685,      -1.05948,      0.461018,      0.877282,     -0.160599,
			0.0689851,     -0.356127,       1.74339,      0.165252,        1.3097,      -0.14737,      0.307117,         1.566,      -1.27838,      -1.62825,       -1.1263,     -0.119203,     -0.942066,      0.277937,       1.11303,
			-0.454128,      -1.64019,     0.0794384,      0.910759,       1.40077,     -0.746058,     -0.421587,    -0.0778837,       0.55997,    -0.0847262,     -0.157805,      -0.76714,     -0.903166,      0.962477,      0.530563,
			-1.58855,    -0.0025439,       1.84188,       1.47399,      0.924899,     -0.534172,     0.0491871,      0.212942,      0.584881,      0.193296,      -1.23542,      -1.20221,      -1.28081,      0.477928,      0.592583,
			-1.11297,        0.1984,      0.296542,     -0.460359,       1.03005,      -1.51536,      -1.52131,       1.87348,     -0.521967,      -2.23603,      -0.26068,      -1.29145,      -0.70812,       1.77035,     0.0297999,
			-0.706762,      -1.23744,      0.199301,      0.621001,      0.865233,     -0.592504,      -1.55651,      0.397432,    -0.0722575,      -1.21823,     -0.757339,      0.235976,      -1.26564,       1.85329,     -0.743313,
			-0.29025,      0.390852,      0.864206,       1.57607,       1.85232,     -0.613634,     0.0749959,       0.98904,     -0.349902,      -1.26836,    -0.0943165,      -1.21245,      -1.42427,       0.23051,       1.26262,
			-0.223272,     -0.667113,      0.797805,      0.704619,     -0.467731,     -0.435472,      -1.87405,      0.507653,       0.11009,     -0.235535,    -0.0371184,     -0.339786,     0.0447813,      0.782643,       1.13992,
			-0.880174,      -0.10528,      0.055733,      0.994779,       1.67835,        -1.587,      0.800349,       1.96891,      0.639732,     -0.104254,     -0.725549,      -1.41419,     -0.906423,       1.04531,       1.50612,
			-0.681804,     -0.208054,      -1.63016,        1.9236,       1.49046,     -0.448473,      -1.66501,      0.229209,      -1.23122,     0.0926658,     -0.660581,     -0.193753,     -0.485145,      0.808703,      0.632284,
			-0.552159,     -0.934482,      0.621354,      0.935002,       1.34724,      0.154053,     -0.921944,       1.12312,      -1.71679,       0.98508,      -1.71881,       1.05893,     -0.982705,       1.88416,      0.250396,
			-0.20448,      0.445696,       1.70458,       1.01871,       1.21593,     -0.947265,      -1.26421,       1.01227,     -0.259719,     -0.490467,      -1.15994,     -0.602008,     -0.261318,       2.27996,     0.0131333,
			-0.599254,     -0.086338,       1.03219,      -1.20987,     -0.163552,       0.16579,      -1.03645,       1.76537,       -1.2118,      -0.99932,      -1.28603,      -1.50709,     -0.589108,       1.50766,      0.141206,
			-0.00459127,     -0.395807,       1.00536,      0.587053,       1.26831,     -0.785832,     -0.318084,      0.696794,     -0.104066,      -0.42127,      -1.95235,      -0.36825,     0.0971004,      0.528329,     -0.590773,
			-0.233189,      0.252867,      0.899735,      0.254832,       1.22715,     -0.207686,     -0.993993,    -0.0412225,     0.0234352,     -0.235004,     -0.779418,      -1.16061,      -1.77459,      0.584868,      0.566674,
			-0.628159,     -0.280705,     0.0482022,      0.694287,        1.3084,      0.214685,      0.290038,      0.629987,     -0.626783,      0.577573,      -1.98008,     0.0367849,     -0.656095,      0.682744,       1.13123,
			-0.833953,      -0.22926,       1.38901,     -0.113628,        0.2097,     -0.426533,      -1.03602,        0.5377,      -1.33049,      -1.22339,     -0.730859,      -0.70146,      0.446892,       1.49525,       1.07706,
			-0.832054,      -1.69876,       0.80097,       0.39422,        2.0607,      -1.31296,      -1.16362,    -0.0492942,      -1.68716,     -0.343316,     -0.230858,       1.23637,      -1.52807,      0.997756,       1.47064,
			0.755912,     -0.464219,    -0.0500083,       1.86418,       1.04504,     -0.595939,      -1.21329,       2.22798,     -0.987053,     -0.399507,     -0.148837,     -0.499143,     -0.614258,     -0.199854,      0.604348,
			-1.95303,     -0.773308,      0.846723,      0.245268,       0.42804,     -0.795615,     -0.547756,      0.287163,       -1.7937,     0.0065707,     -0.199414,      0.554612,      -1.98376,       1.13973,       2.15689,
			-0.163181,      -1.06527,      0.430818,     -0.782545,       1.00586,      -1.06794,     0.0483124,      0.661506,      -1.10121,      -0.42002,      -1.00238,      -1.11814,     -0.333358,       1.06905,      0.458918,
			-0.674022,      -1.01993,       1.23149,     -0.195706,         1.108,    -0.0120518,      0.017637,    -0.0500292,    -0.0727943,     -0.490808,      -1.20616,      -1.30495,    -0.0153405,      0.155802,    0.00816301,
			-1.73752,      -1.21063,     0.0618489,      0.993754,       1.39787,     -0.513768,       -1.2895,       1.03407,      0.312695,      -1.69081,     -0.656142,     -0.241869,     -0.531094,       1.75026,     -0.825405,
			-1.03673,      -1.26294,      0.734005,       1.06365,       1.78842,      -1.00545,      -0.61338,      0.867974,     -0.775803,     -0.897859,      0.469798,    0.00165222,      0.456495,     -0.339453,      0.218942,
			-0.255064,       -1.7465,      0.185646,      0.811529,       1.73159,      -1.85932,     -0.190822,      0.946742,      0.364029,      0.424953,     -0.170258,      -1.34174,      -1.60874,       1.11818,       0.36795,
			-1.09839,      -1.70172,        1.2495,       1.44864,      0.887395,      -1.54713,    -0.0333158,       2.36962,      -1.04917,     -0.805941,      0.855852,     -0.736516,       0.36781,      0.824758,     0.0529315,
			-1.10158,     -0.161344,      -0.56902,       1.51973,     -0.149054,      0.276437,      0.341639,      0.650687,      -0.52976,      -1.12679,      -1.31878,      -2.16497,      -1.09103,      0.104465,      0.809803,
			-0.13279,     -0.562265,     -0.329816,      0.135253,        0.2714,      -1.49534,      -1.11122,       0.83262,      -1.13639,     -0.955266,      -0.52679,     -0.278197,       -1.0358,      0.696232,       1.45577,
			-1.40217,      -1.22676,        1.1768,     0.0895527,    -0.0757564,    -0.0209862,     -0.247286,     -0.771479,      -1.19495,     -0.976928,      -1.06962,      -1.01571,       -1.3153,       1.02821,      0.555102,
			-0.238579,      0.420028,      0.313941,     -0.157621,        1.3237,      -1.44502,     -0.212544,       1.00031,      -1.34941,     -0.840425,     -0.926003,      -1.73961,       -0.7186,       1.42669,      -1.00063,
			-0.992977,      -1.03395,       1.52689,      -0.23316,      0.887835,     -0.327608,     -0.270004,     -0.394426,       -1.1356,      -2.10929,     -0.848248,      -1.31152,     -0.896922,      0.839253,        0.3738,
			-1.37634,     -0.213578,       1.01778,     -0.130296,     0.0596708,      -1.86986,      -1.55827,    -0.0853474,      0.276352,      0.567553,     -0.628724,     -0.902994,     0.0740463,       1.74416,      0.164603,
			0.068018,      0.394182,       1.14236,       1.68651,        1.0609,      -1.01311,     -0.339776,      0.822365,     0.0346976,      -1.15441,     -0.490258,     -0.741471,     -0.111919,       1.02726,      0.370472,
			-1.91192,      -1.43143,      0.148173,      0.210507,      0.303723,      -1.45683,     -0.286254,      0.102864,      -1.12812,      -1.47937,      -1.22433,      -2.00056,      -1.13253,     -0.380317,    -0.0440959,
			0.263491,     -0.907501,      0.832449,     -0.239481,     -0.973634,      -2.30404,     -0.733899,      0.255921,      -1.76589,     -0.411146,     -0.369535,    -0.0791872,       -1.5972,       1.30885,       1.52615,
			-0.460479,     0.0684323,      0.248438,     -0.265957,      0.647176,     -0.710982,     -0.578308,       2.28503,     -0.838411,      0.392832,     -0.195676,      -2.00053,     -0.718944,      0.924912,     -0.363221,
			-0.544737,     -0.433804,       1.66727,       1.06468,      0.649823,       -1.5333,      0.500325,    -0.0513525,      0.715335,     -0.118435,     -0.893174,      -0.64194,     -0.854037,      0.523073,       0.87349,
			-0.226352,       -1.8369,       0.37364,      0.417293,     -0.443111,     -0.401369,       -1.6136,       1.50369,    -0.0739782,     -0.933577,      -1.14512,       0.23086,      -1.09717,       1.03276,        1.0545,
			-1.28922,      -1.54266,      0.838939,      0.298552,      0.588784,     -0.543503,     0.0458853,     -0.846414,     -0.324767,      -1.49853,      -1.98026,      -1.19444,      -1.31036,      0.620825,     0.0820856,
			0.623599,     -0.282106,       1.17326,       1.25971,      0.503276,     -0.563889,     -0.997871,       1.30923,     -0.264462,      -1.36163,      -1.56854,     -0.900793,     -0.418028,       1.03203,    -0.0867369,
			-0.570796,     -0.141643,       1.22659,      0.336032,      0.481824,     -0.512139,      0.566739,       1.70165,      -1.30142,     -0.266782,      -1.62234,      0.729603,      -0.52251,     0.0691787,       1.47571,
			-0.882789,      -2.08742,     -0.435896,        1.5188,     -0.557461,     -0.231917,     -0.803966,       1.09307,      -1.85939,      -1.12397,     -0.424673,     -0.805558,     -0.656904,     -0.270529,       1.11635,
			0.00652263,     -0.894258,       1.32286,       1.05623,      0.219566,     -0.211506,      -1.64664,      0.892856,      -0.09949,     -0.699463,      -1.19054,      -1.63279,       1.02761,      0.470151,     -0.132233,
			-0.77442,      -1.58572,      0.458388,     -0.490801,       1.16285,      -1.19974,     -0.122816,       0.83007,      -1.40522,      0.679483,      0.430082,      -0.73682,     -0.740769,      0.927776,      0.292883,
			1.52789,     -0.112337,      0.900093,       1.47804,      0.155643,      0.226591,       -0.3269,      0.938278,      -1.21444,     -0.290849,     -0.126655,      -1.07596,      -1.24523,       1.03085,       1.16218,
			-0.01019,     -0.425856,       1.22006,     -0.466433,      0.740738,      0.165583,     -0.272619,       1.38268,      0.224986,     -0.275837,      -1.20346,     -0.406407,      0.202881,      0.614757,       1.46759,
			-1.22612,     -0.457523,     -0.668532,     -0.352159,       1.62912,      -1.49844,      -1.11942,       1.26271,     -0.145515,       -1.4157,      0.405893,      -0.59079,     -0.671712,       0.93701,      0.648678,
			-1.3021,     -0.741108,      0.472653,      0.165534,      0.943572,      -1.31674,      0.650422,      0.204705,      -0.65831,     -0.706182,       -1.5168,      0.667835,      -1.35248,       1.27132,       1.07031,
			-0.352538,      -1.36082,       1.51413,      0.402837,      0.852836,     -0.116521,      -2.20653,     -0.557891,      -1.40805,      0.068418,    0.00713633,     -0.508583,     -0.845157,       1.93832,      0.703918,
			-1.46884,    -0.0266065,      0.659164,       1.90214,     -0.622615,      -1.16779,      -1.34884,     -0.160075,      -1.14993,     -0.893107,     -0.465096,      -0.30617,      0.200569,       0.52799,      0.895364,
			-1.01525,       0.45987,       0.81347,      0.995129,      0.115917,     -0.656611,       -2.2729,      0.132462,      -1.12966,     -0.831351,     -0.223842,      -1.21833,      -2.11808,     -0.228404,      0.954736,
			0.0581976,      -1.44033,       1.00589,      0.933725,      0.139758,    -0.0500724,      -0.31604,     -0.365484,      0.211458,     -0.275625,      -1.22439,     -0.955207,     -0.220702,       1.04409,     -0.341267,
			-1.035,     -0.757548,       1.14672,      0.925523,     -0.168942,     -0.407745,      -1.29054,      0.353171,     -0.199608,     -0.609547,      -1.30467,      0.426056,      -1.37928,     -0.864697,     -0.752834,
			-1.21331,     0.0793423,      0.935487,     0.0747527,      0.856164,      0.748673,      -0.83114,       1.51275,     -0.283681,     -0.325874,     -0.517548,      0.123658,    -0.0248652,      0.272223,       1.50971,
			0.117554,      -1.03785,       1.31358,       1.13644,      0.956212,     -0.707269,     -0.837737,      0.104583,      -1.47385,     -0.680494,     -0.667972,      0.454579,      -1.25187,      0.252649,      0.303929,
			-1.77841,      -1.12831,       0.10055,      0.116068,       1.78365,     -0.127406,      -2.01409,       1.35016,      0.280786,     -0.421889,     -0.374725,    -0.0720362,     -0.763039,        1.1631,      0.279605,
			0.215303,       -1.2856,      0.951406,       1.26427,      0.130711,     -0.540343,     -0.976638,      0.486956,      -1.63302,    -0.0205699,      -1.59162,      0.100799,      -1.62257,      0.470901,     -0.767669,
			-1.42068,     -0.863249,       1.31546,     -0.602406,      0.863305,      -0.73203,     -0.758344,     -0.186644,      -1.23531,     -0.931277,     -0.158445,     -0.222002,       0.41624,      0.170195,      0.465069,
			-1.91903,      -1.17177,      0.695103,     0.0913321,       1.28788,      -1.35467,      0.362741,       1.42276,     -0.980874,     -0.709023,      0.753466,      -1.15586,       0.30504,      0.510263,      0.329607,
			0.199025,      0.274936,      0.854088,       1.24832,       1.07876,      -1.67413,     -0.585958,      0.393471,     -0.534814,      -1.85915,     -0.676931,      -1.23907,      0.569124,      0.912498,      0.157726,
			0.373764,      -1.60811,       1.25988,       1.49587,      0.567885,     -0.988902,      -1.29215,     -0.197889,      -1.47572,      -1.21773,      -1.79678,     -0.625682,      -2.16722,     -0.782505,      0.736529,
			-0.767294,     -0.784507,      0.901798,     -0.160056,     0.0766524,      -1.96552,      -1.05539,      -0.29132,      -0.22625,    -0.0141122,     -0.588266,      -2.04506,      -2.28659,     -0.235002,     -0.292645,
			-1.48572,      -1.82935,      0.126413,       1.08966,       1.34212,      -1.19177,      -2.07893,       1.44883,      0.479165,      -0.39437,      -1.73052,     -0.885166,     -0.360581,      0.141173,     -0.154774,
			-1.15427,      0.431658,      0.257102,      0.212635,      0.290256,      -1.15571,     0.0461691,       1.04155,     -0.979678,      -1.63236,     -0.342336,     -0.766815,     -0.395737,      0.980067,     -0.430744,
			-1.15387,      0.163074,      0.610725,     -0.551969,      0.616653,      -1.94694,      0.466323,      0.736777,     -0.655594,      0.697225,      -1.36551,     -0.750635,      -1.49298,      0.627214,      0.868891,
			-1.34706,     -0.508659,       1.26876,      0.329634,      0.709447,      -1.78835,     -0.559478,      0.456309,     -0.445965,      -0.19641,      0.610976,       0.36673,       -1.1667,      0.956761,      0.893856,
			-0.505104,    -0.0174116,        1.6692,      0.226495,       1.90659,    -0.0717907,      -1.13951,     -0.274913,      0.255459,    -0.0995935,      -1.01645,     -0.917246,     -0.582391,      0.680968,       1.53788,
			-0.406351,      -1.64506,       1.29658,    -0.0530493,     -0.131537,      0.791374,     -0.668313,       1.03364,     -0.707411,      -1.24403,     -0.857292,      -1.32007,     -0.460183,      0.897468,     -0.255357,
			-0.0609312,      -0.98352,       1.55711,       1.15531,      0.114518,     -0.838221,      -1.54217,       1.22957,      0.516703,       -1.8205,     -0.262023,     -0.411007,     -0.175446,      0.811865,       1.09089,
			-0.824595,        0.1688,      0.892135,     -0.327565,     -0.218487,     -0.280269,      -1.58621,       1.14771,      -1.20043,       0.27398,     -0.721698,    -0.0441358,      -1.22837,     0.0925498,       1.55437,
			-0.91078,      -1.80863,       1.42962,       1.11755,      0.826765,      -1.54193,     -0.293815,       1.88145,      -2.05955,     -0.708345,     -0.515849,     -0.591009,      0.250594,      -1.10356,      0.397292,
			-0.35099,      0.331642,      0.592164,      0.307087,      0.757958,     -0.880327,      -1.93163,       2.24855,      -0.84826,      -1.75577,     -0.251643,     -0.599602,      -1.05353,      0.616071,     -0.927923,
			-1.88836,      -1.42704,       -0.6504,      0.740925,      0.421833,      0.175182,     -0.830791,       1.29936,     0.0355606,      -1.06206,      -2.19925,      -1.01218,      -0.36999,      0.811628,     -0.181078,
			-0.567672,      0.775319,       0.44045,      0.723107,       2.00661,      0.531847,     -0.678314,       1.17377,      -1.83299,     -0.625139,      -1.62232,     -0.732001,      0.264829,      0.677165,       1.34605,
			-0.949975,     -0.722189,      0.379069,      0.647796,      0.829271,        0.5678,      -1.18871,      0.953125,      0.361963,      -0.17264,      -1.00353,     -0.550722,      -1.19753,       1.09561,    -0.0870362,
			-0.686833,      0.639271,     -0.672129,       1.77611,      0.740196,     -0.681059,       0.13334,      0.996203,      -1.33449,     -0.500089,     -0.994468,      -1.34741,      -1.09886,     -0.210672,       1.24843,
			-1.38309,     -0.018849,      0.860961,      0.592258,      0.649527,      0.271069,      -0.98512,       1.07013,      -1.01737,     -0.261684,      -1.12215,      -0.15446,     -0.265893,     -0.259213,       1.20027,
			-0.454256,      0.296104,      0.362761,      0.934515,       1.53593,      -1.22518,      -1.17363,       1.37638,      -1.47093,      0.462206,     -0.890109,     -0.628413,      -1.21931,      0.346824,       1.67619,
			-0.253597,      -1.61601,       1.46306,       1.59797,      0.700271,      -1.73336,      -1.06062,     -0.135849,       -1.0876,     -0.496225,      -1.28306,       0.34707,     -0.918532,      0.221162,     0.0836565,
			-0.413916,     -0.196794,      0.989183,        1.1627,      0.135258,      -2.24402,       0.11753,      0.164685,      0.912739,      -1.07721,       -1.3366,      -2.01476,     -0.393515,      0.581557,      0.843724,
			-1.28495,       -1.8834,      0.551678,     0.0509215,     -0.218806,      -1.18239,     -0.415447,      0.152962,     -0.335199,       0.54733,     -0.867083,      -1.05084,     -0.163826,       1.94489,        1.0216,
			-0.792998,     -0.838705,     0.0932986,       0.95781,     0.0211142,      0.556256,      -1.03734,       1.32506,      0.412397,     -0.673552,     -0.730927,    0.00245951,       -1.4333,       2.04219,      0.294703,
			0.322991,      0.490005,        1.3598,       1.02225,      0.460116,     -0.637971,      0.194733,       1.19452,     -0.314086,     -0.833702,      -0.54427,       -0.5263,     -0.432385,       0.55459,     -0.469013,
			-0.99527,     0.0363857,      0.277427,       0.76314,       1.51241,     -0.103766,    0.00260088,     0.0631043,      -1.19585,     -0.816851,      0.354576,     -0.774207,     -0.727055,       1.51911,      0.820754,
			-1.83671,      0.609795,       1.87986,      0.154305,      0.208414,     -0.389892,      -1.39762,        1.6842,      -1.35301,      -2.53009,      -1.55933,      -0.89356,     -0.247596,      0.541494,     0.0622022,
			-1.48737,     -0.187033,     0.0167548,       0.65977,      0.769683,      0.184807,      0.684605,       1.70204,     -0.797554,     -0.627398,        -1.128,      -1.22559,     -0.758692,      0.026494,       1.48447,
			-2.06798,      0.118664,       1.38114,       1.60847,        1.4654,      0.371452,      -1.13374,      0.103742,     -0.241983,      0.269961,      -1.17588,      -2.18145,      -1.11246,      0.713739,    -0.0471614,
			-1.52802,     0.0681082,       1.31806,    -0.0321723,       1.88263,      -1.44326,      0.299602,     -0.104474,     -0.058716,    -0.0481299,      -1.22221,     -0.193541,     -0.923945,        1.7953,      0.151431,
			0.214386,      0.445939,      0.168239,       2.25878,    -0.0345405,      -1.08295,      -1.09029,       1.54748,     -0.933329,     -0.344862,      -1.14021,      -1.34605,     -0.823646,     0.0830204,       2.01756,
			-0.236354,     -0.509987,       1.04311,     -0.255529,       1.02584,      0.777809,      0.728534,      0.787091,      -1.07812,     -0.257517,      -1.80434,      -1.25194,     -0.894887,      0.742634,      0.767162,
			0.341759,     -0.685015,     -0.585945,      0.538017,      0.314887,     -0.743433,      -1.14851,      0.510574,     -0.267141,      -1.23717,     -0.670704,      0.225248,      -1.73185,       1.51997,      0.716708,
			-1.70553,     -0.824548,      0.747972,      0.715423,       1.05466,     -0.303242,      -1.76951,      0.811679,      -1.11093,      0.142829,      -1.12633,      0.103334,      0.710229,      0.624791,       1.02723,
			-1.28208,     -0.474614,       1.52947,       0.66345,       0.79952,      -1.51701,     -0.984487,      0.766629,      -1.21667,      -1.26414,      0.585733,      -1.39002,       1.38947,      0.152935,      0.758467,
			-1.07659,      -0.81202,     -0.187534,       0.32498,       1.71472,       -1.5363,      0.179692,      0.434999,     -0.127012,     -0.193508,      0.605466,      -1.28164,       -1.0232,     -0.216714,       1.28907,
			-1.18997,      -0.63035,      0.439014,     -0.190472,       1.65216,      -1.53756,      0.886761,      0.132058,     -0.100244,     -0.745047,      -1.71795,     -0.591945,      -1.28159,       1.98276,       1.53492,
			0.126847,     -0.977212,     -0.833793,     -0.594279,      0.991017,     -0.986959,      -1.08229,      0.482465,     -0.553993,      -1.11237,     -0.621512,      -1.56565,     -0.182919,      0.735694,       1.48251,
			-2.32206,       -1.7424,       1.18627,     -0.661431,      0.662608,      -1.17998,      -1.04533,     -0.434605,      -1.74989,      0.241728,     -0.135967,     -0.809401,      -1.17403,      0.605804,      0.633602,
			-0.656536,      -1.29593,      0.578357,         1.571,       1.58261,    -0.0570147,      0.504825,       1.52476,      -1.81431,    -0.0928947,      0.301903,      -2.10164,      0.402318,       0.11973,       0.67105,
			-0.796905,       -1.3379,      0.897095,       0.91732,       1.24711,      -1.01306,       0.27238,      0.906682,      -1.39327,     -0.429609,     -0.470585,     -0.994768,     -0.120255,    -0.0113366,     -0.595282,
			-1.13564,     -0.761977,       1.73408,      0.647782,       1.59035,      -1.28049,      -1.49452,      0.272459,     -0.111923,     -0.514415,      0.444515,       1.11919,      -1.22409,      0.940102,      0.832338,
			-1.45862,     -0.309813,      0.881624,      0.126106,      0.687739,     -0.799127,     -0.479175,      0.279209,      -1.30622,     -0.256836,    -0.0634967,      0.464278,      -0.85049,       1.06548,       1.83171,
			-1.02863,     -0.013288,     -0.108519,     0.0603884,       1.67076,     -0.852956,      -1.48297,      0.756806,     -0.230352,     -0.251648,      -0.33442,     -0.800843,      -1.36327,      0.885711,     -0.682597,
			-0.704691,     -0.582902,      0.779447,        1.8244,      0.505201,      -1.74191,       0.73605,       1.74205,      -1.77292,      -1.06977,    -0.0492514,    -0.0209397,     0.0534279,      0.733837,     0.0393118,
			0.213144,     -0.345257,      0.752701,       1.21998,       1.69863,      0.318319,      0.146464,      0.356776,      0.510041,      -1.49369,      -1.05214,     -0.900248,     -0.400047,      0.854042,     -0.237271,
			-0.480162,     -0.625512,      0.556941,    -0.0119859,      -1.29351,      0.108929,      -1.02792,      0.505206,      -1.28465,      -1.24341,      -0.54829,      -1.72832,     -0.436392,      0.164356,      0.669316,
			0.54828,     -0.141713,       0.63499,        1.1059,      0.638696,      -0.79344,     -0.567469,      0.570621,      -1.04329,     -0.699169,      -1.44552,     0.0755908,       -0.8981,      0.162011,     -0.342213,
			-0.255374,     -0.852614,     -0.427776,      0.881926,      0.994342,     -0.826339,       -1.2169,     -0.277988,      -1.16867,     -0.623082,     -0.426256,     -0.955069,     -0.912764,     -0.281818,     -0.713418,
			0.0740732,      0.429405,      0.805451,      0.581176,     -0.870163,       -0.7125,     -0.840246,       1.22369,     -0.751424,     -0.598438,      -1.11092,     -0.262799,     -0.578995,     -0.316979,       1.28846,
			-1.29783,     -0.192035,      0.947162,       1.85212,   0.000882778,     -0.577526,      0.120265,     0.0481956,      -1.73587,      -1.64855,      -1.01596,      -1.57215,      0.293682,       1.51702,    -0.0917284,
			0.583221,      -0.87464,      0.380953,       1.24451,      0.584561,     -0.411239,      -1.44761,       0.17075,      -1.70722,     -0.776644,      -1.01935,    -0.0957941,      -1.72571,     0.0918547,      0.296025,
			0.393677,     -0.232628,       1.71059,       1.92167,      0.814026,      -1.19037,      -1.66002,       1.03289,     0.0569119,     -0.990648,     -0.950339,     -0.780791,     -0.276213,        1.1063,     -0.304432,
			-0.626474,        1.0948,       1.10061,      0.150099,       1.13772,     -0.427555,      -1.45895,      0.538658,       0.16762,     -0.237629,     -0.958322,     -0.442995,     -0.646711,      0.452133,       1.45873,
			-0.679966,      0.472348,       1.89354,      0.447663,    -0.0144788,      0.170169,     -0.798249,       1.30077,      -1.15581,      -1.44343,      -2.04558,      -1.49691,     -0.396199,     0.0264414,        1.7039,
			-0.731073,     -0.101324,     -0.849117,     -0.434864,      0.978692,     -0.943568,      -0.25567,       1.36867,       -1.4061,      -1.07845,    -0.0491533,     -0.471553,      -1.31014,       1.43233,      0.637145,
			-0.317976,      -1.63019,   -0.00507553,       2.14814,      0.098676,      -1.88639,      -1.05994,      0.450247,     -0.747692,      -2.09669,     -0.111765,      0.480921,      -2.08167,      -0.77658,      0.479879,
			-0.923673,      -1.41988,     0.0939753,       1.16357,     -0.246299,     0.0749555,     0.0780948,      0.977762,     -0.430543,      -1.56793,      -1.47772,       -1.0155,     -0.278422,    -0.0033685,       1.71588,
			-0.949149,      0.475696,      0.996811,      0.484309,       1.18768,      0.031609,      -1.14787,     -0.240656,     -0.440743,     -0.269533,      -1.38845,       -1.7529,      -1.49427,       1.15216,      0.743257,
			-0.908532,      -0.67634,       1.47623,       1.50669,    0.00496364,     -0.607491,      -1.40406,      0.206293,     -0.371526,      -1.94054,      0.571153,     -0.439124,     0.0863337,      0.853688,       1.68602,
			0.669872,     0.0398194,       1.84495,       1.79957,       1.05393,      -1.20328,        -1.832,        1.2852,     -0.122546,     -0.545498,      -1.60793,      -0.19487,     -0.366894,       1.43604,      0.517859,
			-0.882397,      0.331625,     -0.440087,      0.781058,       1.35328,       0.58557,      0.359164,      0.825104,     -0.920952,     -0.910485,     -0.666563,      -1.51581,      -1.43896,       1.42141,        1.6102,
			-0.853375,      0.288032,      0.873327,      0.729786,      0.893387,      0.368446,      -1.26302,       0.58848,      0.107244,      0.150751,     -0.797214,     -0.336038,      -1.94721,       1.41582,       1.32385,
			-1.34785,      -1.09864,       0.06964,        1.9856,       -0.6539,    0.00170855,     -0.502003,       1.63154,     -0.822679,      -1.11919,      -1.27099,      -1.02176,     -0.794766,     -0.166202,       1.17557,
			0.993313,    -0.0181357,      0.583813,      0.555033,       1.63136,      -1.23922,    -0.0332807,     -0.205448,      -1.15243,      0.152966,      -1.46906,     -0.794717,     -0.483046,      0.039166,        1.4533,
			-0.664019,      -1.48148,    0.00179646,    -0.0834031,       1.27344,     -0.160773,       0.16056,     -0.311465,     -0.746131,      -1.39493,     -0.877595,      -2.02086,      -1.74683,       1.67838,      0.279335,
			-0.401585,      -1.79021,       0.96112,      0.305725,      0.925103,     -0.476181,       -0.5882,     -0.332849,      -2.42786,     -0.142659,      0.258258,     0.0512966,     -0.495886,       1.27325,      0.437794,
			-0.592501,     -0.736731,     -0.392304,      0.938686,       1.54991,     -0.190764,       -1.1292,       1.44477,     -0.283658,     -0.602738,     -0.723278,     -0.295857,       -1.2232,      0.823468,     -0.164363,
			0.263419,     -0.487573,      0.959039,      0.500768,       1.31939,     -0.820025,      0.915985,       1.92826,      -1.15742,      -0.20124,     -0.975602,     -0.635223,     0.0387943,       1.25113,       1.01436,
			-0.410897,     -0.464358,      0.277765,       1.28404,      0.346731,      0.909251,        1.1315,       1.41458,     -0.318815,      -1.06696,      -1.31306,      -2.09051,      -1.40421,       1.27975,       1.21372,
			-1.2687,      -1.39457,      0.442862,       1.13556,      0.777946,       0.27045,      0.375968,      0.637638,     -0.840895,     -0.473554,       -1.5655,     -0.771169,     -0.633736,     -0.684586,       0.60943,
			-0.905019,      0.518286,       1.15831,    -0.0815173,      0.307808,     -0.351352,      -1.17604,      0.767493,      -1.79143,     -0.292815,      -1.25072,      -1.00545,     -0.110585,      0.447356,       1.47628,
			-0.203743,     -0.680096,       1.34766,       1.12381,        1.6025,      -1.38056,     -0.971103,     0.0798288,      -1.36553,      -1.30607,     -0.680217,      0.253027,     -0.238609,     0.0872553,     -0.145986,
			-1.19152,     -0.488878,       0.76367,       1.34776,      0.136317,      0.107435,     -0.982275,      0.961231,       -1.2029,      -0.12372,     -0.490867,      0.140506,      0.182311,      0.134539,       1.19252,
			-1.28858,     -0.373943,       1.69469,    -0.0107414,     -0.882853,       -1.2378,       -1.7188,       1.53202,      -1.66412,     -0.240711,       -1.1523,      0.129253,      -1.12587,      0.581724,       1.15935,
			-1.14914,      -1.67328,       1.33587,        1.2337,    -0.0946285,     -0.491575,      -1.45512,      0.407904,      0.540011,    -0.0346314,     0.0161328,       -2.2772,     -0.775517,      0.511256,       1.69045,
			-0.319635,      0.219777,     -0.702188,       1.33626,       1.17218,      -1.18364,      -1.79473,      0.621986,     -0.653061,     -0.530498,     -0.562353,      -1.29324,     -0.846934,     0.0114817,      0.913073,
			-0.995062,      -1.07086,      0.510449,      0.191979,      0.680283,      -1.20919,      -1.10382,      0.364447,     -0.830626,      0.551459,      0.303792,      0.174746,     -0.960979,      0.864578,       1.75948,
			-0.549105,     -0.743039,       0.59787,      0.310984,      0.703185,     -0.614878,      -1.64432,      0.848169,      0.655222,     -0.459479,     -0.148104,      -1.08766,     -0.250702,      0.365924,     -0.557287,
			-0.120527,     -0.553215,       1.88168,     -0.222982,     -0.050275,     -0.392169,    -0.0326031,       1.45306,      -1.06094,     -0.125866,      -1.10811,      -1.07406,     -0.190305,       1.75298,     -0.593313,
			-0.815533,     0.0317984,       1.41163,       0.58226,      0.316482,      -1.12598,      -0.55269,      0.318874,      -1.50403,     -0.983122,      0.951912,     -0.624979,     -0.235169,     -0.408805,       1.55827,
			-1.14188,      -1.52933,      0.606316,      0.518752,       1.31393,      0.713301,      -1.23441,      0.531349,      0.910492,      -0.94534,      -1.45882,     -0.482351,      -1.40242,      0.146291,      0.520362,
			-1.49841,     -0.320041,       1.46019,      0.306468,      0.675463,     -0.605304,      0.804221,      0.285122,      0.149199,     -0.737339,      -1.20669,     0.0216587,      -1.26147,       1.74634,     -0.576526,
			-0.870479,       -1.2073,       1.12368,      0.423391,      0.126952,      -1.10026,      -1.19216,       1.34985,      0.428217,      -1.39924,      0.515917,      -0.19418,     -0.382095,     -0.533772,      0.222449,
			0.0416618,     -0.593253,      0.577359,       1.69934,       0.99478,      -1.38067,      -1.96976,     -0.725995,      -1.44763,      -1.30174,      -1.32543,     -0.666712,      -1.13098,    -0.0736963,      0.475792,
			0.245257,      -1.37802,      0.239739,    -0.0511969,       0.77022,     -0.317024,      -1.35766,      0.180833,      -0.99922,     -0.135993,      0.626751,      -1.10744,      -1.24589,      0.359323,      0.260165,
			-0.73421,      -1.80971,      0.909369,      0.363261,       0.96896,      -1.52014,     -0.341851,        0.7278,      -1.11513,     -0.928926,     0.0590561,      -1.83749,      0.342303,      0.113298,     -0.208289,
			-0.710357,      0.569466,       1.04794,     -0.208213,     -0.331586,    -0.0252024,     -0.720014,      0.594367,      -1.38202,    -0.0948741,     -0.774833,      -1.70446,     -0.819331,      -0.26376,       1.29007,
			0.433752,      -1.12718,    -0.0598945,      0.186788,      0.622876,      -1.49003,      -0.54398,       0.69441,     -0.322454,     -0.953578,    -0.0953122,      -1.26494,      0.854365,         1.543,       1.61725,
			-1.74233,       -1.5464,       1.46708,      0.663664,       1.04803,     -0.381855,       -1.2981,       1.27716,     -0.863967,       0.56312,      -1.87714,      0.122011,     -0.143465,       0.50056,       1.06693,
			-1.41133,      -1.78115,      0.605467,       2.03454,      0.226305,      -1.36241,     -0.808115,       1.66438,     -0.503142,     -0.429706,      0.625889,     -0.247667,     -0.341272,       1.43861,        1.5818,
			-0.935826,     -0.302708,       1.06359,      0.551794,     -0.531056,      -1.64242,      -1.09701,     -0.881222,      -2.02457,      -1.15756,     -0.410337,      -2.42918,    0.00516787,      0.451219,       1.67123,
			-1.45407,      -1.09073,      0.227369,       1.73243,     0.0775926,      0.315049,     -0.306553,      0.950548,      -1.71769,     -0.534154,      -1.53487,     -0.825237,      -1.40603,       0.20182,       2.07089,
			-0.0882364,    -0.0732137,      0.549485,       2.13273,      0.970781,      -0.89709,     -0.867001,        0.2688,     -0.261222,      -1.64734,      -1.16051,      -1.99249,      0.509469,      0.840659,      0.195183,
			-0.0933706,      -1.57444,       2.10167,       1.31573,      0.797545,      -1.51677,      -2.31029,      0.218861,      -1.20959,      -1.21056,      -0.13679,      0.104196,      -1.10624,      0.171975,     -0.226247,
			-0.146005,      -1.15066,      0.623144,        1.3853,       -0.7638,      -1.39765,     -0.642698,     0.0805745,     0.0187597,      -1.59021,     -0.303883,      -1.43747,      -1.58151,      0.632849,       1.15827,
			0.152601,      -2.36961,      0.374226,     -0.464747,      0.292434,    -0.0361984,     -0.138938,      0.790688,     -0.813702,     -0.935847,     -0.799513,     0.0730829,      -1.51939,       1.66739,        1.5026,
			-0.196404,     -0.596348,      0.787735,       1.60377,     -0.117836,     -0.964378,     -0.711294,      0.908551,     -0.113521,     -0.250558,      0.513477,     -0.293201,      -0.16963,      0.596376,        0.4238,
			1.03966,     -0.984911,       0.57323,     -0.135425,      0.634011,     -0.461738,      0.146308,       1.27429,    -0.0177474,     -0.588245,     -0.278282,      -0.95256,     -0.681392,       1.79739,       1.18626,
			-1.85166,     -0.175552,       1.61952,       0.17123,     -0.481499,     -0.261665,      -1.14453,      0.838567,      -2.21333,      -1.28847,       -1.1405,      -1.75302,      0.225129,      0.699835,       2.31372,
			-0.970602,      -1.22702,      0.950376,       1.38703,      0.849228,      -1.61612,     -0.158162,        2.2467,     -0.870499,     -0.389855,     -0.254549,     -0.930756,     -0.291389,     -0.425188,     0.0976011,
			-1.00953,      -1.56204,     -0.394762,     0.0705088,       1.29894,     -0.444062,     -0.174732,     -0.148291,     -0.874677,     -0.897738,      -1.26547,      -1.18979,      -1.26778,       1.28427,       -0.4478,
			-0.214089,      0.183904,      -1.06198,     -0.113447,       1.88464,      -0.18549,     -0.750936,       2.01831,      -1.44682,      0.332818,     -0.649082,      -1.11552,     -0.304511,      0.668242,     0.0283147,
			-0.253645,     -0.719517,       1.94309,     -0.233555,      0.713142,    -0.0404066,     0.0577595,    0.00575312,     -0.322046,      -1.14314,      -1.36072,      -1.15052,     -0.987971,       1.58908,      0.101073,
			-0.947954,     -0.727525,       0.86413,     -0.545745,      0.231764,      -1.35157,     -0.101407,      0.425669,     -0.549241,     0.0460293,      -2.18337,      -1.63353,      0.568033,      0.725305,       1.33905,
			-0.979134,    0.00847705,       0.93903,     -0.609326,     -0.225797,      0.168181,       -1.6411,      0.490987,     -0.881885,      -1.34857,      -1.77204,        -1.052,     -0.274494,       1.41125,    -0.0127812,
			-1.33258,       -1.3782,      0.587122,      0.307786,      0.101758,      0.687118,     -0.274779,       1.37856,      -1.36918,     -0.639784,      -1.32088,     -0.161441,     -0.674042,      0.148209,       1.17403,
			0.503905,      -0.39821,       2.01526,      0.967839,       1.05786,     -0.665397,     -0.255361,      0.825317,      -2.12449,     -0.394825,      0.785801,    -0.0687699,    -0.0489322,      0.806981,      0.877037,
			-1.15851,     -0.121347,      0.705347,       1.57298,      0.826242,      -0.58291,    -0.0860467,     -0.528146,     -0.406484,     -0.999565,      -1.14263,     -0.999693,    -0.0492986,      0.875733,       -0.7883,
			-1.43428,     -0.416832,       1.70272,       1.01377,     -0.262463,      -2.04936,      -1.65441,      0.267989,      -1.14351,      -1.11516,      0.435566,     -0.988638,     -0.274334,      0.260663,       1.19072,
			-1.44502,      0.553376,      0.866699,       2.06779,      -0.16132,     -0.600651,      -1.12018,      0.347053,      -2.04226,     0.0525032,     -0.214149,      -1.89772,     -0.924626,      0.482633,       2.02375,
			-0.80964,      -0.34614,       1.07237,     -0.478005,      0.395946,     0.0194471,      -1.52323,       1.02245,     -0.850444,     -0.756659,     -0.857424,      -1.49323,    -0.0237568,       1.30793,      0.874585,
			-1.30514,      -1.65185,       1.88523,      0.830975,      0.163901,     -0.149507,      0.526322,      0.988063,     -0.435984,       0.21022,      -2.06582,      -1.31954,      0.254651,       1.33882,      0.227837,
			0.00624556,      -1.26663,       0.39403,       2.12172,      0.360222,      -1.61011,      -1.85846,     -0.617583,       -1.5649,      -1.89179,     0.0205843,        0.2105,      -1.17101,     0.0385012,      0.597482,
			-0.742598,     -0.495527,      0.273267,     -0.764588,     0.0393151,     -0.314272,      0.263317,      0.684149,     -0.584878,      -1.93455,      -1.94736,      -1.16222,      0.100295,        1.2181,       1.35667,
			0.0110305,     -0.719706,     0.0573174,     -0.573576,      0.790239,      -1.24382,      -1.11728,       1.96167,       -0.4917,      -1.62359,      0.328664,     -0.849741,    -0.0668597,       1.97442,      0.175732,
			0.670529,     0.0280804,      0.666594,      0.917169,      0.627837,     -0.976665,      -1.45162,        1.2741,     -0.958636,      0.206955,      -1.07351,     -0.168041,      -1.40401,     -0.255482,      0.904402,
			-1.19025,      -0.72193,      0.350625,      0.949516,       1.07759,      -1.16838,     -0.438194,       2.02835,      0.469303,      0.429995,       -0.5761,       -1.5075,     -0.205443,      0.459108,     -0.199291,
			-0.634953,     -0.767737,       1.35592,       1.68668,       1.47753,     -0.372873,     -0.572099,      0.230653,      -1.06372,      0.521507,     -0.866835,     -0.520676,      0.164742,    -0.0251011,       1.68542,
			-1.94945,      0.565287,       1.56474,       1.07415,      0.346954,      0.129269,     -0.478744,       1.30032,      0.255413,      0.444856,     -0.579065,     -0.757352,     -0.767873,       1.26343,     -0.249456,
			-1.18875,      -1.33281,      0.646258,      0.306192,       1.39062,     -0.934848,       -1.2028,      0.143006,      -1.01226,     -0.151309,     0.0855735,      0.315944,       -1.1119,      0.654284,      0.562476,
			-0.870893,      -2.15874,     0.0240763,      0.693543,     -0.362682,      -0.21617,  -0.000904812,        1.2333,     -0.433097,      -1.42123,     -0.435275,       -1.4729,       -2.0833,     -0.271169,      0.640186,
			-1.74505,      0.263874,     -0.717968,    -0.0345249,      0.486287,      0.179901,     -0.681052,      0.724084,     -0.801788,      -0.78888,     -0.631377,      -1.73643,      -1.22233,       1.13112,       0.19853,
			-0.37185,      0.197179,      0.827896,       1.22602,       1.85958,      -1.54643,      0.628288,       1.57892,      -1.18835,     -0.820655,    -0.0626798,      -1.31984,      -0.32468,      0.249091,     -0.201489,
			-1.38558,      0.112545,      0.691983,        1.4146,     0.0214252,     -0.416255,      -0.87919,     -0.216624,      0.460944,      -0.44518,      -1.63653,      -1.18342,     -0.278035,       1.16057,       1.27254,
			-0.472955,     -0.775907,       1.24567,      0.806332,      0.967147,      -2.28854,      -1.70081,      0.505479,     -0.171172,      -1.42683,     -0.673271,    -0.0477808,      -1.03712,     -0.211031,      -1.03048,
			-0.39818,      0.460221,     0.0427841,       1.29998,         1.817,      -1.72896,      -1.41647,        1.5157,     -0.119612,     -0.776786,    -0.0296403,      -1.62324,      0.421559,      0.893425,      0.315745,
			-0.117662,     -0.641813,       0.48838,       1.06168,      0.477541,     -0.112084,     -0.569509,     0.0663246,      0.678679,      -1.06015,     -0.344857,       -1.5569,     -0.771945,    -0.0437308,        2.0684,
			-0.965894,     -0.901436,      0.989654,     -0.102293,       1.44062,     -0.907382,      0.940968,     -0.334211,     -0.139741,      -1.49768,      -1.76347,     0.0623678,      -1.38234,       0.78076,      0.252544,
			});
		*/

		Matrix<my_double> recv_matrix(simulation_times, n, {
			-1.71704,	-0.0355986,	0.953004,	1.03202,	1.27793,	-0.130623,	-0.81576,	1.57263,	-1.14507,	-0.151133,	-0.500191,	-0.54,	0.550899,	-0.165649,	1.15493,
		});

		// decoding parameter
		int order = 2;
		OSD_v2 osd(GM, true, d_min);		// have problem in high SNR

		// performance recording

		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> recv = recv_matrix.get_row(i);
			//Matrix<my_double> recv = AWGN::pass(c);
			Matrix<GF2> v_hat = osd.decode_v(recv, order);

			if (v == v_hat) {
				cout << "@@@ correct at i = " << i << endl;
			}
			else {
				cout << "error at i = " << i << endl;
			}
			//cout << "simulation_times = " << i + 1 << '\r';
		}

	}
	static void OSD_PSM_aided_GE_simulation(bool& is_first_time, my_double SNR_dB = 2.0) {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 5000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 5;
		const int t = 3;
		GF2e<m>::init();
		eBCH<m, t> bch;			// here we choose eBCH
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> PM = bch.get_parity_matrix();
		Matrix<GF2> GM = bch.get_generator_matrix();
		//bch.print_info();

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//cout << "v" << v;

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		int order = 2;
		int partition_num = 2;
		OSD_PSM_aided_GE osd(PM, false, partition_num, d_min);		// have problem in high SNR

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;

		if (is_first_time) {

			bch.print_info();
			cout << "max_error_frame          = " << max_error_frame << endl;
			cout << "simulation type          = OSD-" << order << endl;

			cout << "\n------------------------------------------------------------------------------------------------------------"
				<< "--------------------------------------------------------------------------------" << endl;


			cout << setw(7) << "SNR_dB" << setw(20) << "simulation_times" << setw(20) << "error_rate" \
				<< setw(20) << "ave_used_list" << setw(20) << "early_stop_rate";

#ifdef count_operation_number
#ifdef use_my_double

			cout << setw(20) << "double_ope_num";
#endif // use_my_double

			cout << setw(20) << "GF2_ope_num" << setw(20) << "GF2e_ope_num";
#endif // count_operation_number

			cout << setw(20) << "time_consume(s)" << endl;

			is_first_time = false;
		}

		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
		unsigned long long double_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;

		unsigned long long GF2e_ope_num_before = GF2e_auxiliary_storage::operation_number;
		unsigned long long GF2e_ope_num_after;
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> recv = AWGN::pass_standard(c);
			//Matrix<my_double> recv = AWGN::pass(c);			
			Matrix<GF2> v_hat = osd.decode_v(recv, order, PM);

			if (v == v_hat);
			else {
				// record the error frame
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
			//cout << "simulation_times = " << i + 1 << '\r';
		}
		end = clock();

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << setw(7) << fixed << setprecision(1) << SNR_dB;
		cout << setw(20) << simulation_times;
		cout << scientific << setprecision(2);
		cout << setw(20) << error_frame / (double)simulation_times;

		cout << setw(20) << osd.total_used_list_num / (double)simulation_times;
		cout << setw(20) << osd.num_early_stop / (double)simulation_times;

#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << setw(20) << (double_ope_num_after - double_ope_num_before) / (double)simulation_times;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << setw(20) << (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << setw(20) << (GF2e_ope_num_after - GF2e_ope_num_before) / (double)simulation_times;
#endif // count_operation_number

		cout << setw(20) << time_consume - pass_AWGN_time << endl;

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);

	}
	static void OSD_PSM_aided_GE_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "OSD_r_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(2, 0.5, 2, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			OSD_PSM_aided_GE_simulation(is_first_time, test_SNR(i));
		}
	}
	static void OSD_PSM_aided_GE_4_cyc_simulation(bool& is_first_time, my_double SNR_dB = 2.0) {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 5000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;			// here we choose BCH
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> PM = bch.get_parity_matrix();
		Matrix<GF2> GM = bch.get_generator_matrix();
		//bch.print_info();

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//cout << "v" << v;

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		int order = 1;
		int partition_num = 2;
		OSD_PSM_aided_GE_4_cyc osd(PM, false, partition_num, d_min);		// have problem in high SNR

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;

		if (is_first_time) {

			bch.print_info();
			cout << "max_error_frame          = " << max_error_frame << endl;
			cout << "simulation type          = OSD-" << order << endl;

			cout << "\n------------------------------------------------------------------------------------------------------------"
				<< "--------------------------------------------------------------------------------" << endl;


			cout << setw(7) << "SNR_dB" << setw(20) << "simulation_times" << setw(20) << "error_rate" \
				<< setw(20) << "ave_used_list" << setw(20) << "early_stop_rate";

#ifdef count_operation_number
#ifdef use_my_double

			cout << setw(20) << "double_ope_num";
#endif // use_my_double

			cout << setw(20) << "GF2_ope_num" << setw(20) << "GF2e_ope_num";
#endif // count_operation_number

			cout << setw(20) << "time_consume(s)" << endl;

			is_first_time = false;
		}

		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
		unsigned long long double_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;

		unsigned long long GF2e_ope_num_before = GF2e_auxiliary_storage::operation_number;
		unsigned long long GF2e_ope_num_after;
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> recv = AWGN::pass_standard(c);
			//Matrix<my_double> recv = AWGN::pass(c);
			Matrix<GF2> v_hat = osd.decode_v(recv, order, PM);

			if (v == v_hat);
			else {
				// record the error frame
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
			//cout << "simulation_times = " << i + 1 << '\r';
		}
		end = clock();

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << setw(7) << fixed << setprecision(1) << SNR_dB;
		cout << setw(20) << simulation_times;
		cout << scientific << setprecision(2);
		cout << setw(20) << error_frame / (double)simulation_times;

		cout << setw(20) << osd.total_used_list_num / (double)simulation_times;
		cout << setw(20) << osd.num_early_stop / (double)simulation_times;

#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << setw(20) << (double_ope_num_after - double_ope_num_before) / (double)simulation_times;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << setw(20) << (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << setw(20) << (GF2e_ope_num_after - GF2e_ope_num_before) / (double)simulation_times;
#endif // count_operation_number

		cout << setw(20) << time_consume - pass_AWGN_time << endl;

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);

	}
	static void OSD_PSM_aided_GE_4_cyc_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "OSD_r_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(2, 0.5, 6, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			OSD_PSM_aided_GE_4_cyc_simulation(is_first_time, test_SNR(i));
		}
	}
	static void OSD_PSM_aided_GE_4_ext_cyc_simulation(bool& is_first_time, my_double SNR_dB = 2.0) {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 100000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		eBCH<m, t> bch;			// here we choose eBCH
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> GM = bch.get_generator_matrix();
		Matrix<GF2> PM = bch.get_parity_matrix();
		//bch.print_info();

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//cout << "v" << v;

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		int order = 2;
		int partition_num = 2;
		OSD_PSM_aided_GE_4_ext_cyc osd(PM, false, partition_num, d_min);		// have problem in high SNR

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;

		if (is_first_time) {

			bch.print_info();
			cout << "max_error_frame          = " << max_error_frame << endl;
			cout << "simulation type          = OSD-" << order << endl;

			cout << "\n------------------------------------------------------------------------------------------------------------"
				<< "--------------------------------------------------------------------------------" << endl;


			cout << setw(7) << "SNR_dB" << setw(20) << "simulation_times" << setw(20) << "error_rate" \
				<< setw(20) << "ave_used_list" << setw(20) << "early_stop_rate";

#ifdef count_operation_number
#ifdef use_my_double

			cout << setw(20) << "double_ope_num";
#endif // use_my_double

			cout << setw(20) << "GF2_ope_num" << setw(20) << "GF2e_ope_num";
#endif // count_operation_number

			cout << setw(20) << "time_consume(s)" << endl;

			is_first_time = false;
		}

		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
		unsigned long long double_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;

		unsigned long long GF2e_ope_num_before = GF2e_auxiliary_storage::operation_number;
		unsigned long long GF2e_ope_num_after;
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> recv = AWGN::pass_standard(c);
			//Matrix<my_double> recv = AWGN::pass(c);
			Matrix<GF2> v_hat = osd.decode_v(recv, order, PM);

			if (v == v_hat);
			else {
				// record the error frame
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
			//cout << "simulation_times = " << i + 1 << '\r';
		}
		end = clock();

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << setw(7) << fixed << setprecision(1) << SNR_dB;
		cout << setw(20) << simulation_times;
		cout << scientific << setprecision(2);
		cout << setw(20) << error_frame / (double)simulation_times;

		cout << setw(20) << osd.total_used_list_num / (double)simulation_times;
		cout << setw(20) << osd.num_early_stop / (double)simulation_times;

#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << setw(20) << (double_ope_num_after - double_ope_num_before) / (double)simulation_times;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << setw(20) << (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << setw(20) << (GF2e_ope_num_after - GF2e_ope_num_before) / (double)simulation_times;
#endif // count_operation_number

		cout << setw(20) << time_consume - pass_AWGN_time << endl;

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);

	}
	static void OSD_PSM_aided_GE_4_ext_cyc_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "OSD_r_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(6, 0.5, 6, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			OSD_PSM_aided_GE_4_ext_cyc_simulation(is_first_time, test_SNR(i));
		}
	}

	static void OSD_v2_and_PSM_aided_double_ope_test() {
		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 1;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 6;
		const int t = 5;
		GF2e<m>::init();
		eBCH<m, t> bch;			// here we choose eBCH
		int k = bch.get_k();
		int n = bch.get_n();
		int d_min = bch.get_d();
		Matrix<GF2> PM = bch.get_parity_matrix();
		//bch.print_info();

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//cout << "v" << v;

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR_dB = 2;
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		int order = 2;
		OSD_v2 osd(PM, false, d_min);		// have problem in high SNR

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;


		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
		unsigned long long double_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;

		unsigned long long GF2e_ope_num_before = GF2e_auxiliary_storage::operation_number;
		unsigned long long GF2e_ope_num_after;
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> recv = AWGN::pass_standard(c);
			//Matrix<my_double> recv = AWGN::pass(c);
			//Matrix<GF2> v_hat = osd.decode_v(recv, order);
			Matrix<GF2> v_hat = (GF2)0;				// test on operation with out decoding

			if (v == v_hat);
			else {
				// record the error frame
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
			//cout << "simulation_times = " << i + 1 << '\r';
		}
		end = clock();

#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "double operation = " << (double_ope_num_after - double_ope_num_before) / (double)simulation_times << endl;
		// minus n is for operation number of passing AWGN channel
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "GF2 operation = " << (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times << endl;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << "GF2e operation = " << (GF2e_ope_num_after - GF2e_ope_num_before) / (double)simulation_times << endl;
#endif // count_operation_number

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void OSD_v2_simulation_complexity_distribution(my_double SNR_dB = 2.0) {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 100000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 7;
		const int t = 4;
		GF2e<m>::init();
		BCH<m, t> bch;			// here we choose BCH
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> PM = bch.get_parity_matrix();
		//bch.print_info();

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//cout << "v" << v;

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		int order = 2;
		OSD_v2 osd(PM, false, d_min);				// have problem in high SNR

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;			// no use 

		ofstream my_double_ope_file("my_double_ope.txt");
		ofstream GF2_ope_file("GF2_ope.txt");
		ofstream correct_flag_file("correct_flag_file.txt");

		if (my_double_ope_file.is_open() && GF2_ope_file.is_open()) { // check if file is opened
			// wirte the contents into file

			for (int i = 0; i < simulation_times; ++i) {
				Matrix<my_double> recv = AWGN::pass_standard(c);
				//Matrix<my_double> recv = AWGN::pass(c);


#ifdef count_operation_number
#ifdef use_my_double

				unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
				unsigned long long double_ope_num_after;
#endif // use_my_double

				unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
				unsigned long long GF2_ope_num_after;

#endif // count_operation_number


				Matrix<GF2> v_hat = osd.decode_v(recv, order);


#ifdef count_operation_number
#ifdef use_my_double

				double_ope_num_after = my_double_auxiliary_storage::operation_number;
				my_double_ope_file << (double_ope_num_after - double_ope_num_before) << " ";
#endif // use_my_double

				GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
				GF2_ope_file << (GF2_ope_num_after - GF2_ope_num_before) << " ";

#endif // count_operation_number

				if (v == v_hat) {
					correct_flag_file << 1 << " ";
				}
				else {
					// record the error frame
					correct_flag_file << 0 << " ";
				}
			}

			// close file
			my_double_ope_file.close();
			GF2_ope_file.close();
			correct_flag_file.close();
		}
		else {
			cout << "cannot open file(s)" << std::endl;
		}
	}
};

class test_IBU {
public:

	static void IBU_init_test() {
		test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_test_ini();
	}
	static void IBU_test() {
		test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_test_multi_trial();
	}

	static void IBUc_init_test() {
		test_NewGE_in_OSD::PreStored_Matrix_red_cycle_ini_test();
	}
	static void IBUc_test() {
		test_NewGE_in_OSD::PreStored_Matrix_red_cycle_test_multi_trial();
	}
};
