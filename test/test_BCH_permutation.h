/*****************************************************************//**
 * \file   test_BCH_permutation.h
 * \brief  BCH codeword permutation function
 * 
 * \author 26259
 * \date   March 2024
 *********************************************************************/

#pragma once

#include"test_common.h"

class test_BCH_permutation {
public:

	static int eBCH_affine_permutation() {
		// to be finished, index mapping: 
		//		X -> \alpha * X + \beta, \alpha,\beta\in GF(2^m)

		const int m = 3;
		const int t = 1;
		GF2e<m>::init();
		eBCH<m, t> ebch;
		Matrix<GF2> G = ebch.get_generator_matrix();
		Matrix<GF2> H = ebch.get_parity_matrix();
		cout << "G" << G;
		cout << "H" << H;
		cout << "G.multiply_transpose_of(H)" << G.multiply_transpose_of(H);

		//Matrix<GF2> H_sys = H;
		//Matrix<int> p = H_sys.GE_left_identity_4_GF2();

		//cout << "H_sys" << H_sys;
		//cout << "p" << p;
		//cout << "G.multiply_transpose_of(H_sys)" << G.multiply_transpose_of(H_sys);

		const int n = 1 << m;		// number of elements in field, also codeowrd length
		const int n_minus_1 = n - 1;
		int k = ebch.get_k();
		Matrix<int> mapping_permutation_all(n * (n - 1), n, '0');
		int insert_row_ind = 0;
		Matrix<int> mapping_permutation(1, n);

		// affine mapping for X -> t * X + r
		for (int wi = 1; wi < n; ++wi) {
			for (int ri = 0; ri < n; ++ri) {

				GF2e<m> w, r;
				w = wi;		// w can not set to 0, totally n-1 choices
				r = ri;		// n choices
				for (int i = 0; i < n_minus_1; ++i) {
					GF2e<m> orig;
					orig.set_by_alpha_power(i);
					GF2e<m> dist = w * orig + r;
					if (dist == 0) {
						mapping_permutation[i] = n_minus_1;
					}
					else {
						mapping_permutation[i] = GF2e_auxiliary_storage::alpha_table[(int)dist];
					}
				}
				if (r == 0) {
					mapping_permutation[n_minus_1] = n_minus_1;
				}
				else {
					mapping_permutation[n_minus_1] = GF2e_auxiliary_storage::alpha_table[(int)r];
				}
				//cout << "(" << wi << ", " << ri << "), " << "mapping_permutation" << mapping_permutation;
				// check duplicated permutation
				bool permutation_dup = false;
				for (int i = 0; i < insert_row_ind; ++i) {
					if (mapping_permutation == mapping_permutation_all.get_row(i)) {
						cout << "duplicated permutation found in (" << wi << ", " << ri << ")" << endl;
						cout << "------" << endl;
						permutation_dup = true;
						break;
					}
				}
				if (permutation_dup == false) {
					// insert the permutation
					mapping_permutation_all.set_row(insert_row_ind, mapping_permutation);
					insert_row_ind++;
				}

				Matrix<GF2> G_new = G;
				G_new.permute_col(mapping_permutation);
				//cout << "G_new" << G_new;
				if (!G_new.multiply_transpose_of(H).isZero()) {
					cout << "G_new.multiply_transpose_of(H)" << G_new.multiply_transpose_of(H);	
						// it validates that G_new is a generator matrix
				}
			}
		}

		mapping_permutation_all.resize(insert_row_ind, n);
		cout << "mapping_permutation_all" << mapping_permutation_all << endl;

		// now puncture position '7' and see what are the distinct permutation left
		Matrix<int> BCH_permutation(1, insert_row_ind * (n - 1), 'v');
		int push_num = mapping_permutation_all.size();
		for (int i = 0; i < push_num; ++i) {
			if (mapping_permutation_all[i] != n - 1) {
				BCH_permutation.push_back(mapping_permutation_all[i]);
			}
		}
		BCH_permutation.resize(insert_row_ind, n - 1);

		// recude the duplicated permutation in 'BCH_permutation'
		Matrix<int> duplicated_row_ind(1, insert_row_ind, 'v');
		for (int i = 0; i < insert_row_ind; ++i) {
			// check out the formmer duplicated permutation
			for (int j = 0; j < i; ++j) {
				if (BCH_permutation.get_row(i) == BCH_permutation.get_row(j)) {
					// detected duplication
					duplicated_row_ind.push_back(i);
					break;
				}
			}
		}

		cout << "BCH_permutation = " << BCH_permutation;
		cout << "duplicated_row_ind = " << duplicated_row_ind;

		BCH_permutation = BCH_permutation.erase_rows(duplicated_row_ind);
		cout << "BCH_permutation (erased dup rows) = " << BCH_permutation;

		// reduce the duplicated information set in 'BCH_permutation'
		duplicated_row_ind.resize(1, 0);
		cout << "k = " << k << endl;
		int BCH_permutation_rows = BCH_permutation.row();
		for (int i = 0; i < BCH_permutation_rows; ++i) {
			// check out the formmer duplicated permutation
			Matrix<int> ri = BCH_permutation.get_part(i, 0, i, k - 1);
			ri.sort('<');
			for (int j = 0; j < i; ++j) {

				Matrix<int> rj = BCH_permutation.get_part(j, 0, j, k - 1);
				rj.sort('<');

				if (ri == rj) {
					// detected duplication
					duplicated_row_ind.push_back(i);
					break;
				}
			}
		}

		cout << "duplicated_row_ind (info_Set) = " << duplicated_row_ind;
		BCH_permutation = BCH_permutation.erase_rows(duplicated_row_ind);
		cout << "BCH_permutation (erased dup rows and dup info_Set) = " << BCH_permutation;

		// verify that these permutation is really preserving BCH codes
		BCH<m, t> bch;
		Matrix<GF2> Gb = bch.get_generator_matrix();
		Matrix<GF2> Hb = bch.get_parity_matrix();

		cout << "Gb = " << Gb << endl;
		cout << "Hb = " << Hb << endl;

		cout << "Gb.multiply_transpose_of(Hb).isZero() = " << Gb.multiply_transpose_of(Hb).isZero() << endl;

		// check each permutation
		int BCH_permutation_size = BCH_permutation.row();
		for (int i = 0; i < BCH_permutation_size; ++i) {
			Matrix<GF2> Gb_permuted = Gb;
			Gb_permuted.permute_col(BCH_permutation.get_row(i));

			cout << "(i = " << i << ") Gb_permuted.multiply_transpose_of(Hb).isZero() = " \
				<< Gb_permuted.multiply_transpose_of(Hb).isZero() << endl;
		}

		/**
		 * .	You can see that not all permutation pass the parity-check of BCH code.
		 *		Because the row space of eBCH code after puncturing position 7, is not the row space of
		 *			the first 7 columns of the generator matrix of eBCH codes.
		 *		Therefore, puncturing position 7 in permutation that perserves eBCH code can not generate
		 *			the permutation that preserve BCH codes.
		 */

		// --------------------------------------------------------------------------------------------------------------

		 // reduce the duplicated information set in 'mapping_permutation_all'
		duplicated_row_ind.resize(1, 0);
		cout << "k = " << k << endl;
		int eBCH_permutation_rows = mapping_permutation_all.row();
		for (int i = 0; i < eBCH_permutation_rows; ++i) {
			// check out the formmer duplicated permutation
			Matrix<int> ri = mapping_permutation_all.get_part(i, 0, i, k - 1);
			ri.sort('<');
			for (int j = 0; j < i; ++j) {

				Matrix<int> rj = mapping_permutation_all.get_part(j, 0, j, k - 1);
				rj.sort('<');

				if (ri == rj) {
					// detected duplication
					duplicated_row_ind.push_back(i);
					break;
				}
			}
		}
		cout << "duplicated_row_ind (info_Set) = " << duplicated_row_ind;
		mapping_permutation_all = mapping_permutation_all.erase_rows(duplicated_row_ind);
		cout << "mapping_permutation_all (erased dup rows and dup info_Set) = " << mapping_permutation_all;

		/**
		 * . The generated eBCH preserving permutation has different information set
		 */

		return 0;
	}

	static int eBCH_affine_permutation2() {
		// to be finished, index mapping: 
		//		X -> \alpha * X + \beta, \alpha,\beta\in GF(2^m)

		const int m = 4;
		const int t = 2;
		GF2e<m>::init();
		eBCH<m, t> ebch;
		Matrix<GF2> G = ebch.get_generator_matrix();
		Matrix<GF2> H = ebch.get_parity_matrix();
		cout << "G" << G;
		cout << "H" << H;
		cout << "G.multiply_transpose_of(H)" << G.multiply_transpose_of(H);

		//Matrix<GF2> H_sys = H;
		//Matrix<int> p = H_sys.GE_left_identity_4_GF2();

		//cout << "H_sys" << H_sys;
		//cout << "p" << p;
		//cout << "G.multiply_transpose_of(H_sys)" << G.multiply_transpose_of(H_sys);

		const int n = 1 << m;		// number of elements in field, also codeowrd length
		const int n_minus_1 = n - 1;
		int k = ebch.get_k();
		Matrix<int> mapping_permutation_all(n * (n - 1), n, '0');
		int insert_row_ind = 0;
		Matrix<int> mapping_permutation(1, n);

		int wi = 1;
		// affine mapping for X -> X + r
		for (int ri = 0; ri < n; ++ri) {

			GF2e<m> w, r;
			w = wi;		// w can not set to 0, totally n-1 choices
			r = ri;		// n choices
			for (int i = 0; i < n_minus_1; ++i) {
				GF2e<m> orig;
				orig.set_by_alpha_power(i);
				GF2e<m> dist = w * orig + r;
				if (dist == 0) {
					mapping_permutation[i] = n_minus_1;
				}
				else {
					mapping_permutation[i] = GF2e_auxiliary_storage::alpha_table[(int)dist];
				}
			}
			if (r == 0) {
				mapping_permutation[n_minus_1] = n_minus_1;
			}
			else {
				mapping_permutation[n_minus_1] = GF2e_auxiliary_storage::alpha_table[(int)r];
			}
			//cout << "(" << wi << ", " << ri << "), " << "mapping_permutation" << mapping_permutation;
			
			// insert the permutation
			mapping_permutation_all.set_row(insert_row_ind, mapping_permutation);
			insert_row_ind++;

			Matrix<GF2> G_new = G;
			G_new.permute_col(mapping_permutation);
			//cout << "G_new" << G_new;
			if (!G_new.multiply_transpose_of(H).isZero()) {
				cout << "G_new.multiply_transpose_of(H)" << G_new.multiply_transpose_of(H);
				// it validates that G_new is a generator matrix
			}
		}

		mapping_permutation_all.resize(insert_row_ind, n);
		cout << "mapping_permutation_all" << mapping_permutation_all << endl;

		// extend the permutation by cyclic shift and *2 mapping, generating in total (3*7*8) = 168 permutations

		int s = (int)my::_log2(n);
		Matrix<int> ext_eBCH_permutations(1, s * n_minus_1 * n * n, 'v');

		Matrix<int> nat(1, n_minus_1, 'N');
		Matrix<int> U(s, n_minus_1);
		for (int j = 0; j < s; ++j) {
			// *2 permutations
			for (int w = 0; w < n_minus_1; ++w) {
				U(j, w) = (nat(w) << j) % (n_minus_1);
			}
		}
		cout << "U = " << U;

		int current_row_inserted = 0;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < s; ++j) {
				// *2 permutations
				Matrix<int> p_orig = mapping_permutation_all.get_part(i, 0, i, n_minus_1 - 1);
				Matrix<int> pU = p_orig;
				pU.permute(U.get_row(j));

				for (int w = 0; w < n_minus_1; ++w) {
					// shifting permutations
					Matrix<int> pUT(1, n_minus_1);
					pU.col_shift_right_cir(w, pUT);

					// insert permutations
					for (int q = 0; q < n_minus_1; ++q) {
						ext_eBCH_permutations.push_back(pUT(q));
					}
					ext_eBCH_permutations.push_back(mapping_permutation_all(i, n_minus_1));
				}
			}
		}

		ext_eBCH_permutations.resize(s * n_minus_1 * n, n);

		cout << "ext_eBCH_permutations = " << ext_eBCH_permutations;

		// verify that these permutation is really preserving BCH codes


		// check each permutation
		int eBCH_permutation_size = ext_eBCH_permutations.row();
		for (int i = 0; i < eBCH_permutation_size; ++i) {
			Matrix<GF2> G_permuted = G;
			G_permuted.permute_col(ext_eBCH_permutations.get_row(i));

			if (G_permuted.multiply_transpose_of(H).isZero() == false) {
				cout << "(i = " << i << ") Gb_permuted.multiply_transpose_of(Hb).isZero() = false" << endl;
			}
		}

		// recude the duplicated permutation in 'BCH_permutation'
		Matrix<int> duplicated_row_ind(1, eBCH_permutation_size, 'v');
		for (int i = 0; i < eBCH_permutation_size; ++i) {
			// check out the formmer duplicated permutation
			for (int j = 0; j < i; ++j) {
				if (ext_eBCH_permutations.get_row(i) == ext_eBCH_permutations.get_row(j)) {
					// detected duplication
					duplicated_row_ind.push_back(i);
					break;
				}
			}
		}

		cout << "duplicated_row_ind = " << duplicated_row_ind;
		ext_eBCH_permutations = ext_eBCH_permutations.erase_rows(duplicated_row_ind);

		// reduce the duplicated information set in 'ext_eBCH_permutations'
		duplicated_row_ind.resize(1, 0);
		cout << "k = " << k << endl;
		int eBCH_permutation_size_2 = ext_eBCH_permutations.row();
		for (int i = 0; i < eBCH_permutation_size_2; ++i) {
			// check out the formmer duplicated permutation
			Matrix<int> ri = ext_eBCH_permutations.get_part(i, 0, i, k - 1);
			ri.sort('<');
			for (int j = 0; j < i; ++j) {

				Matrix<int> rj = ext_eBCH_permutations.get_part(j, 0, j, k - 1);
				rj.sort('<');

				if (ri == rj) {
					// detected duplication
					duplicated_row_ind.push_back(i);

					cout << "duplicated info sets" << endl;
					cout << ext_eBCH_permutations.get_row(i);
					cout << ext_eBCH_permutations.get_row(j);
					cout << " ---------- " << endl;

					break;
				}
			}
		}

		cout << "duplicated_row_ind (info_Set) = " << duplicated_row_ind;
		ext_eBCH_permutations = ext_eBCH_permutations.erase_rows(duplicated_row_ind);
		cout << "ext_eBCH_permutations (erased dup rows and dup info_Set) = " << ext_eBCH_permutations;

		// re-initialize 'mapping_permutation_all'
		insert_row_ind = 0;
		for (int wi = 1; wi < n; ++wi) {
			// affine mapping for X -> X + r
			for (int ri = 0; ri < n; ++ri) {

				GF2e<m> w, r;
				w = wi;		// w can not set to 0, totally n-1 choices
				r = ri;		// n choices
				for (int i = 0; i < n_minus_1; ++i) {
					GF2e<m> orig;
					orig.set_by_alpha_power(i);
					GF2e<m> dist = w * orig + r;
					if (dist == 0) {
						mapping_permutation[i] = n_minus_1;
					}
					else {
						mapping_permutation[i] = GF2e_auxiliary_storage::alpha_table[(int)dist];
					}
				}
				if (r == 0) {
					mapping_permutation[n_minus_1] = n_minus_1;
				}
				else {
					mapping_permutation[n_minus_1] = GF2e_auxiliary_storage::alpha_table[(int)r];
				}
				//cout << "(" << wi << ", " << ri << "), " << "mapping_permutation" << mapping_permutation;

				// insert the permutation
				mapping_permutation_all.set_row(insert_row_ind, mapping_permutation);
				insert_row_ind++;

				Matrix<GF2> G_new = G;
				G_new.permute_col(mapping_permutation);
				//cout << "G_new" << G_new;
				if (!G_new.multiply_transpose_of(H).isZero()) {
					cout << "G_new.multiply_transpose_of(H)" << G_new.multiply_transpose_of(H);
					// it validates that G_new is a generator matrix
				}
			}
		}

		mapping_permutation_all.resize(insert_row_ind, n);
		cout << "mapping_permutation_all (new)" << mapping_permutation_all;
		
		// reduce the duplicated information set in 'mapping_permutation_all'

		duplicated_row_ind.resize(1, 0);
		cout << "k = " << k << endl;
		int mapping_permutation_all_size = mapping_permutation_all.row();
		for (int i = 0; i < mapping_permutation_all_size; ++i) {
			// check out the formmer duplicated permutation
			Matrix<int> ri = mapping_permutation_all.get_part(i, 0, i, k - 1);
			ri.sort('<');
			for (int j = 0; j < i; ++j) {

				Matrix<int> rj = mapping_permutation_all.get_part(j, 0, j, k - 1);
				rj.sort('<');

				if (ri == rj) {
					// detected duplication
					duplicated_row_ind.push_back(i);

					cout << "duplicated info sets" << endl;
					cout << mapping_permutation_all.get_row(i);
					cout << mapping_permutation_all.get_row(j);
					cout << " ---------- " << endl;

					break;
				}
			}
		}

		cout << "duplicated_row_ind (info_Set) = " << duplicated_row_ind;
		mapping_permutation_all = mapping_permutation_all.erase_rows(duplicated_row_ind);
		cout << "mapping_permutation_all (erased dup rows and dup info_Set) = " << mapping_permutation_all;

		return 0;
	}

	static int BCH_weight() {
		const int m = 5;
		const int t = 3;
		GF2e<m>::init();
		eBCH<m, t> bch;
		bch.print_info();

		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> G_sys = bch.get_generator_matrix();
		Matrix<int> G_sys_col_permute_record(1, n, 'N');
		G_sys.GJE_4_GF2_left(G_sys_col_permute_record);
		cout << "G_sys" << G_sys;

		Weight_Specturm ws(G_sys);

		int max_weight_allowed = n;
		bool codewords_needed = false;
		ws.solve(max_weight_allowed, codewords_needed);

		cout << "weight specturm = " << ws.ans;
		if (codewords_needed == true) {
			cout << "codewords = " << ws.codewords;
		}

		return 0;
	}

	static int BCH_dual_weight() {
		const int m = 7;
		const int t = 21;
		GF2e<m>::init();
		BCH<m, t> bch;
		bch.print_info();

		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> H_sys = bch.get_parity_matrix();
		Matrix<int> H_sys_col_permute_record(1, n, 'N');
		H_sys.GJE_4_GF2_left(H_sys_col_permute_record);
		cout << "H_sys" << H_sys;

		Weight_Specturm ws(H_sys);

		int max_weight_allowed = 6;
		bool codewords_needed = false;
		ws.solve(max_weight_allowed, codewords_needed);

		cout << "weight specturm = " << ws.ans;
		if (codewords_needed == true) {
			cout << "codewords = " << ws.codewords;
		}

		system("pause");
		return 0;
	}

	static int BCH_punctured_wieght() {
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		bch.print_info();

		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> G = bch.get_generator_matrix();

		// randomly punctured 2*t bits
		Matrix<int> nat(1, n, 'N');
		srand(64);
		Matrix<int> random_erased_col = nat.get_random_element(2 * t);
		Matrix<GF2> G_sys_punctured = G.erase_cols(random_erased_col, false);

		Matrix<int> G_sys_punctured_col_permute_record(1, n - 2 * t, 'N');
		G_sys_punctured.GJE_4_GF2_left(G_sys_punctured_col_permute_record);
		//cout << "G_sys_punctured" << G_sys_punctured;
		//cout << "G_sys_punctured_col_permute_record" << G_sys_punctured_col_permute_record;

		Weight_Specturm ws(G_sys_punctured);

		int max_weight_allowed = 6;
		bool codewords_needed = true;
		ws.solve(max_weight_allowed, codewords_needed);

		cout << "weight specturm = " << ws.ans;
		/*if (codewords_needed == true) {
			cout << "codewords = " << ws.codewords;
		}*/

		Matrix<int> flip_pos(1, 2);
		flip_pos(0) = 54;
		flip_pos(1) = 56;
		int max_weight_cnt = 3;
		Matrix<unsigned long long> bias_weight = Weight_Specturm_transform::solve(ws.codewords, flip_pos, max_weight_cnt);

		cout << "bias_weight" << bias_weight;

		return 0;
	}

	static int BCH_punctured_wieght_multi() {
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		bch.print_info();

		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<int> nat(1, n, 'N');

		int max_weight_allowed = 3;
		bool codewords_needed = true;

		Matrix<my_double> weight_ave(max_weight_allowed + 1, 2, '0');
		for (int i = 0; i <= max_weight_allowed; ++i) {
			weight_ave(i, 0) = i;
		}

		int max_test_num = 1000;
		Weight_Specturm ws(n - 2 * t, k);
		for (int i = 0; i < max_test_num; ++i) {

			// randomly punctured 2*t bits
			Matrix<int> random_erased_col = nat.get_random_element(2 * t);
			Matrix<GF2> G_sys_punctured = G.erase_cols(random_erased_col, false);

			Matrix<int> G_sys_punctured_col_permute_record(1, n - 2 * t, 'N');
			G_sys_punctured.GJE_4_GF2_left(G_sys_punctured_col_permute_record);
			//cout << "G_sys_punctured" << G_sys_punctured;
			//cout << "G_sys_punctured_col_permute_record" << G_sys_punctured_col_permute_record;

			ws.G_sys = G_sys_punctured;
			ws.solve(max_weight_allowed, codewords_needed);

			//cout << "weight specturm = " << ws.ans;
			for (int j = 0; j < ws.ans.row(); ++j) {
				weight_ave((unsigned)ws.ans(j, 0), 1) += ws.ans(j, 1);
			}
		}

		for (int i = 0; i <= max_weight_allowed; ++i) {
			weight_ave(i, 1) /= max_test_num;
		}
		cout << "weight_ave" << weight_ave;

		return 0;
	}

	static int BCH_punctured_wieght_flip_multi() {
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		bch.print_info();

		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<int> nat(1, n, 'N');

		int max_weight_cnt = 3;
		int flip_cnt = 1;
		int max_weight_allowed = max_weight_cnt + flip_cnt;
		bool codewords_needed = true;

		Matrix<my_double> weight_ave(max_weight_cnt + 1, 2, '0');
		for (int i = 0; i <= max_weight_cnt; ++i) {
			weight_ave(i, 0) = i;
		}

		Matrix<int> flip_pos(1, flip_cnt);
		for (int j = 0; j < flip_cnt; ++j) {
			flip_pos(j) = n - 2 * t - j - 1;
		}
		cout << "flip_pos" << flip_pos;

		int max_test_num = 10;
		Weight_Specturm ws(n - 2 * t, k);
		for (int i = 0; i < max_test_num; ++i) {

			// randomly punctured 2*t bits
			Matrix<int> random_erased_col = nat.get_random_element(2 * t);
			Matrix<GF2> G_sys_punctured = G.erase_cols(random_erased_col, false);

			Matrix<int> G_sys_punctured_col_permute_record(1, n - 2 * t, 'N');
			G_sys_punctured.GJE_4_GF2_left(G_sys_punctured_col_permute_record);
			//cout << "G_sys_punctured" << G_sys_punctured;
			//cout << "G_sys_punctured_col_permute_record" << G_sys_punctured_col_permute_record;

			ws.G_sys = G_sys_punctured;
			ws.solve(max_weight_allowed, codewords_needed);
			
			Matrix<unsigned long long> bias_weight = Weight_Specturm_transform::solve(ws.codewords, flip_pos, max_weight_cnt);

			//cout << "weight specturm = " << bias_weight;
			for (int j = 0; j < bias_weight.row(); ++j) {
				weight_ave((unsigned)bias_weight(j, 0), 1) += bias_weight(j, 1);
			}
		}

		for (int i = 0; i <= max_weight_cnt; ++i) {
			weight_ave(i, 1) /= max_test_num;
		}
		cout << "weight_ave" << weight_ave;

		system("pause");
		return 0;
	}

	static int nnsBCH_weight() {
		const int m = 5;
		const int t = 3;
		GF2e<m>::init();
		nnsBCH<m, t> bch;
		bch.print_info();

		int n = bch.get_n();
		int k = bch.get_k();

		Matrix<GF2> G_sys = bch.get_generator_matrix();
		Matrix<int> G_sys_col_permute_record(1, n, 'N');
		G_sys.GJE_4_GF2_left(G_sys_col_permute_record);
		cout << "G_sys" << G_sys;

		Weight_Specturm ws(G_sys);

		int max_weight_allowed = n;
		bool codewords_needed = false;
		ws.solve(max_weight_allowed, codewords_needed);

		cout << "weight specturm = " << ws.ans;
		if (codewords_needed == true) {
			cout << "codewords = " << ws.codewords;
		}

		return 0;
	}
};
