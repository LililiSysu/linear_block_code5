#pragma once
/*****************************************************************//**
 * \file   MRIP_Reprocessing.h
 * \brief  implementation of OSD decode algorithm
 * 
 * \author 26259
 * \date   October 2022
 *********************************************************************/
#include"reprocess_common.h"
#include"MRIP_PSM.h"

/**
 * .discarded, or not recommend to use
 */
class OSD{
public:
	static bool is_ML;

	static Matrix<int> sum_dependent_column;

	/**
	 * .get the most reliable independent positions
	 * 
	 * \param G: generator matrix from code instance, the process change it to systematic G1
	 * \param r_abs: reliability of each received bit
	 * \return the permutation result
	 */
	static Matrix<int> get_MRIP(Matrix<GF2>& G_or_H, Matrix<my_double> r_abs, bool is_by_G = true) {
		// include producing systematic G through parity check matrix


		// step 1: turn r to abs(r) and arrange r_abs from big to small
		//Matrix<my_double> r_abs_tmp = r_abs;

		/*cout << "G_or_H=" << G_or_H;*/
		Matrix<int> Pi = r_abs.sort_with_ind('>');		// r_abs_tmp are thrown away
		//r_abs.permute_back(Pi);		// this process return the original un-permuted r_abs

		if (is_by_G) {
			// G_or_H is G
			G_or_H.permute_col(Pi);
			/*cout << "G=" << G;
			cout << "Pi=" << Pi;*/

			// step2: Gaussian cancellation to find independent positions
			G_or_H.row_transformation_to_up_triangle();

			// compute number of dependent column, can be fix to delete
			/*for (int i = G_or_H.row() - 1; i < G_or_H.col(); ++i) {
				if (G_or_H(G_or_H.row() - 1, i) != 0) {
					sum_dependent_column(i - G_or_H.row() + 1)++;
					break;
				}
			}*/

			//cout << "G=" << G;
			//cout << "Pi=" << Pi;
			Matrix<int> nPi = G_or_H.col_permute_to_full_rank_on_left();
			Pi.permute(nPi);
			//G_or_H.col_permute_to_full_rank_on_left(Pi);

			//cout << "G=" << G;
			//cout << "Pi=" << Pi;
			G_or_H.row_transformation_left_up_triangle_to_identity();

			//cout << "G=" << G;
			//cout << "Pi=" << Pi;

			// by now the first k col of G are most reliable independent positions
			// that is the first k element as index in r, are most reliable independent positions
			// Pi is the record of perbutation

			return Pi;
		}
		else {
			// G_or_H is H

			G_or_H.permute_col(Pi);
			/*cout << "G=" << G;
			cout << "Pi=" << Pi;*/

			// step2: Gaussian cancellation to find independent positions
			//G_or_H.row_transformation_to_up_triangle();
			G_or_H.row_transformation_to_low_triangle();

			//cout << "G=" << G;
			//cout << "Pi=" << Pi;
			//G_or_H.col_permute_to_full_rank_on_left(Pi);
			Matrix<int> nPi = G_or_H.col_permute_to_full_rank_on_right();
			Pi.permute(nPi);

			//cout << "G=" << G;
			//cout << "Pi=" << Pi;
			//G_or_H.row_transformation_left_up_triangle_to_identity();
			G_or_H.row_transformation_right_low_triangle_to_identity();

			//cout << "G=" << G;
			//cout << "Pi=" << Pi;
			int r = G_or_H.row();
			int n = G_or_H.col();
			int k = n - r;

			r_abs.permute(nPi);
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
			G_or_H = G_or_H.get_part(0, 0, r - 1, k - 1);		// the parity check part
			G_or_H.permute_col(nPi);

			G_or_H = Matrix<GF2>(k, k, 'i').combine_right(G_or_H.Transpose());

			// by now the first k col of G are most reliable independent positions
			// that is the first k element as index in r, are most reliable independent positions
			// Pi is the record of perbutation

			return Pi;

		}
	}

	/**
	 * .a warp function of decode_v, see decode_v
	 */
	template<class code_type>
	static Matrix<GF2> decode(code_type& code_instance, Matrix<my_double> r, int order = -1) {
		return code_instance.v2u(decode_v(code_instance, r, order));
	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by flipping all possible <=order bits in k most reliable independent position
	 *
	 * \param code_instance: the coding and decoding scheme
	 * \param r: received codeword, after BPSK_demodulation
	 * \param order: order of OSD decoding, all possible bit combinations with number under order, will be test vectors
	 * \return decoded result, u
	 */
	template<class code_type>
	static Matrix<GF2> decode_v(code_type& code_instance, Matrix<my_double> r, int order = -1) {
		// a shell for code_type, such as BCH, return code_instance.decode(r);

		is_ML = false;
		int n = r.size();	// also known as 
		int k = code_instance.get_k();
		bool isby_G = (my_double)k / n <= 0.5;		// in general, false
		//bool isby_G = true;		// fix for compute dependent column

		Matrix<GF2> G = code_instance.get_generator_matrix();
		Matrix<GF2> H = code_instance.get_parity_matrix();
		Matrix<GF2> G_or_H = isby_G ? G : H;

		Matrix<my_double> r_abs(1, n);
		for (int i = 0; i < n; ++i) {
			r_abs(i) = my::abs(r(i));
		}
		Matrix<int> Pi = get_MRIP(G_or_H, r_abs, isby_G);	// finally, G_or_H will be G1
		Matrix<GF2> G1 = G_or_H;
		//cout << "Pi=" << Pi; cout << "r=" << r;

		// test code
		//for (int i = 0; i < n; ++i) {
		//	r_abs(i) = my::abs(r(i));
		//}
		//Matrix<int> Pit = get_MRIP(G, r_abs, true);	// finally, G_or_H will be G1
		//for (int i = 0; i < n; ++i) {
		//	r_abs(i) = my::abs(r(i));
		//}
		//Matrix<int> Pif = get_MRIP(H, r_abs, false);	// finally, G_or_H will be G1
		//if (G != H) {
		//	cout << "Pit-Pif" << Pit - Pif;
		//}

		r.permute(Pi);
		for (int i = 0; i < n; ++i) {
			r_abs(i) = my::abs(r(i));		// stay with r(i)
		}

		//cout << "r=" << r;
		Matrix<GF2> hdr = BPSK::demodulation(r);
		//cout << "hdr=" << hdr;
		int error_num;
		int d_min = code_instance.get_d();
		order = order == -1 ? (d_min/4) : order;		// if not given any parameter, set order=d_min/4 as default
		//order = 0;
		//cout << "order=" << order << endl;
		
		// packed G1 for paralelle computation
		int col_packed_len = (n - 1) / 32 + 1;												// n=G1.col()
		vector<Matrix<unsigned>> G1_row_packed(k, Matrix<unsigned>(1, col_packed_len));		// k=G1.row()
		
		//cout << "G1" << G1;
		data_ope::pack_row_data(G1, G1_row_packed);
		//cout << "G1_row_test";
		//for (auto iter = G1_row_packed.begin(); iter != G1_row_packed.end(); ++iter) {
		//	Matrix<GF2> G1_row_test(1, n);
		//	data_ope::unpack_data(*iter, G1_row_test);
		//	//cout << G1_row_test;
		//}

		// first not flip any position, for the first set we do demand in v space
		//cout << "hdr" << hdr;

		//cout << "col_packed_len=" << col_packed_len << endl;
		Matrix<unsigned> dv_packed(1, col_packed_len, '0');
		for (int rr = 0; rr < k; ++rr) {
			if (hdr(rr) != 0) {
				for (int i = 0; i < col_packed_len; ++i) {
					dv_packed(i) ^= G1_row_packed[rr](i);		// + in GF2 is ^ in unsigned
				}

#ifdef count_operation_number
				GF2_auxiliary_storage::operation_number += n;	// count the operation number correctly, add n bit in total
#endif // count_operation_number

				//Matrix<GF2> G1_row_test(1, n);
				//data_ope::unpack_data(dv_packed, G1_row_test);
				//cout << "G1_row_test" << G1_row_test;
			}
		}
		Matrix<GF2> best_dv(1, n);
		data_ope::unpack_data(dv_packed, best_dv);
		//cout << "best_dv: 0" << best_dv;

		//Matrix<GF2> origin_dv(1, n);
		//for (int i = 0; i < k; ++i) {
		//	origin_dv(i) = hdr(i);		// known that G1 is systemmatic
		//}
		//for (int i = k; i < n; ++i) {
		//	origin_dv(i) = 0;		// known that G1 is systemmatic
		//}
		//for (int j = 0; j < k; ++j) {
		//	for (int i = k; i < n; ++i) {
		//		origin_dv(i) += hdr(j) * G1(j, i);
		//	}
		//}
		//cout << "best_dv" << best_dv;
		//cout << "origin_dv" << origin_dv;
		//cout << "diff num=" << best_dv.Hamming_distance(hdr) << endl;
		//cout << "diff num=" << best_dv.Hamming_distance(origin_dv) << endl;

		//cout << "G1" << G1; cout << "origin_dv" << origin_dv;

		// best_dv after permute back must in v space
		/*if (best_dv.get_permute_back(Pi) != v) {
			cout << "detected error" << endl;
			cout << "e" << e;
			return code_instance.v2u(best_dv.get_permute_back(Pi));
		}*/
		my_double lambda_min = Measure::correlation_discrepancy_v(r_abs, hdr, best_dv);
		error_num = Measure::error_num;
		if (optimum_condition_OSD::is_optimum(d_min, k, best_dv, hdr, r_abs, -1, lambda_min, error_num)) {
			//cout << "early stop, o=0" << endl;
			is_ML = true;
			best_dv.permute_back(Pi);
			return best_dv;
		}

		int current_test_n = 0;
		Matrix<unsigned> iter_dv_packed(1, col_packed_len);
		Matrix<GF2> dv(1, n);
		for (int o = 1; o <= order; ++o) {
			Matrix<int> flipped_mat(1, o, 'N');
			//cout << "current_test_n=" << current_test_n << endl;
			//cout << "o=" << o << endl;
			//Matrix<GF2> dv = best_dv;		// this is alright
			iter_dv_packed = dv_packed;
			Matrix<int> last_flipped_mat(1, o);
			Matrix<int> fmd(1, 2 * o);
			int fmd_size = 0;				// actual useful size of fmd
			do{
				//cout << "flipped_mat" << flipped_mat;
				flip_mat_discrepancy(last_flipped_mat, flipped_mat, fmd, fmd_size);	// leave it, to be replaced.
				
				//cout << "fmd_size=" << fmd_size << "\t" << "fmd" << fmd;
				//cout << "fmd_size="<<fmd_size << endl;
				for (int kk = 0; kk < fmd_size; ++kk) {
					for (int i = 0; i < col_packed_len; ++i) {
						iter_dv_packed(i) ^= G1_row_packed[fmd(kk)](i);		// + in GF2 is ^ in unsigned
					}


#ifdef count_operation_number
					GF2_auxiliary_storage::operation_number += n;	// count the operation number correctly, add n bit in total
#endif // count_operation_number

				}
				data_ope::unpack_data(iter_dv_packed, dv);

				//for (int kk = 0; kk < fmd_size; ++kk) {
				//	int add_row = fmd(kk);
				//	// systematic part
				//	dv(add_row) += 1;
				//	for (int pp = k; pp < n; ++pp) {
				//		// this for loop cost 40% of total time in MRIP

				//		// parity check part
				//		dv(pp) += G1(add_row, pp);		// another way to speed up: pack the data as int having 32 GF2
				//	}
				//}
				//cout << "dv" << dv;

				//if (!code_instance.is_in_v_space(dv.get_permute_back(Pi)))	
					//cout << "error" << endl; // for later dv, we demand it in v space
				my_double lambda = Measure::correlation_discrepancy_v(r_abs, hdr, dv);
				error_num = Measure::error_num;
				if (lambda >= lambda_min);
				else{
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					//cout << "dv=" << dv; cout << "lambda_min=" << lambda << endl;
					if (optimum_condition_OSD::is_optimum(d_min, k, best_dv, hdr, r_abs, o - 1, lambda_min, error_num)) {
						//cout << "early stop, o=" << o << endl;
						is_ML = true;
						best_dv.permute_back(Pi);
						return best_dv;
					}
				}
			} while (flip_TEP::next(flipped_mat, k, o));
		}

		best_dv.permute_back(Pi);
		return best_dv;
	}

	/**
	 * .compute difference between 'last_flipped_mat' and 'flipped_mat'
	 * 
	 * \param last_flipped_mat: the last flipped pattern
	 * \param flipped_mat: current flip pattern
	 * \param fmd: this will be set as difference between 'last_flipped_mat' and 'flipped_mat'
	 * i.e., 'last_flipped_mat'={2,3}, and 'flipped_mat'={2,4}, then 'fmd'={3,4}
	 * \param fmd_size: valid size of fmd, in case above this will be set as 2
	 */
	static void flip_mat_discrepancy(Matrix<int>& last_flipped_mat,\
		const Matrix<int>& flipped_mat, Matrix<int>& fmd, int& fmd_size) {

		int len = flipped_mat.size();	// which is 'o'
		if (fmd_size != 0) {
			// flipped_mat sizes are equal
			fmd_size = 0;
			for (int i = 0; i < len; ++i) {
				if (last_flipped_mat(i) == flipped_mat(i));
				else {
					fmd(fmd_size) = last_flipped_mat(i);
					fmd_size++;
					fmd(fmd_size) = flipped_mat(i);
					fmd_size++;
					last_flipped_mat(i) = flipped_mat(i);		// set last = now
					//result.push_back(last_flipped_mat(i));
					//result.push_back(flipped_mat(i));
				}
			}
		}
		else {
			for (int i = 0; i < len; ++i) {
				fmd(i) = flipped_mat(i);
			}
			fmd_size = len; // set it equals to 'o'
			last_flipped_mat = flipped_mat;
		}
	}
};

bool OSD::is_ML = false;
Matrix<int> OSD::sum_dependent_column(1,43,'0'); //p<=N-K-d+1, case for eBCH(128,64,22)

class OSD_r {
public:
	
	Matrix<GF2> PM;		// we restrict ourself to consider PM only
	Matrix<GF2> right_systematic_PM;
	int d_min;
	int n;
	int k;
	int r;

	int num_early_stop;
	unsigned long long total_used_list_num;		// counting used list num
	bool is_ML;
	my_double lambda_min;
	
	OSD_r(const Matrix<GF2>& _PM, int _d_min) :PM(_PM), d_min(_d_min) {
		n = PM.col();
		k = n - PM.row();
		r = n - k;
		num_early_stop = 0;
		total_used_list_num = 0;
		is_ML = false;
		lambda_min = 0;		// doesn't matter
	}

	void init(const Matrix<GF2>& _PM, int _d_min){
		PM = _PM;
		d_min = _d_min;		// for ML stopping rule
		n = PM.col();
		k = n - PM.row();
		r = n - k;
		num_early_stop = 0;
		total_used_list_num = 0;
		is_ML = false;
	}

	Matrix<GF2> encode_by_right_systematic_PM(const Matrix<GF2>& information_bit) {
		Matrix<GF2> ans(1, n);		// unit test to be done

		for (int i = k; i < n; ++i) {	// parity check bit
			ans(i) = 0;
		}

		for (int i = 0; i < k; ++i) {	// information bit, systematically encoded
			ans(i) = information_bit(i);
			if (ans(i) != 0) {
				for (int j = 0; j < r; ++j) {	// parity check bit, add the column of Patrity check matrix
					ans(j + k) += right_systematic_PM(j, i);
				}
			}			
		}
		return ans;
	}

	// the input and output will be ref_codeword 
	void ref_encode_by_right_systematic_PM(Matrix<GF2>& ref_codeword, const Matrix<int>& flip_ind) {
		int flip_ind_size = flip_ind.size();
		//cout << "flip_ind_size = " << flip_ind_size << endl;
		for (int i = 0; i < flip_ind_size; ++i) {	// information bit, systematically encoded
			int actual_col = flip_ind(i);		// reverse the flip ind, starting from the least reliable opsition in MRIP
			ref_codeword(actual_col) += 1;
			for (int j = 0; j < r; ++j) {	// parity check bit, add the column of Patrity check matrix
				ref_codeword(j + k) += right_systematic_PM(j, actual_col);
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
		/*cout << "G=" << G;
		cout << "Pi=" << Pi;*/

		// step2: Gaussian cancellation to find independent positions
		//G_or_H.row_transformation_to_up_triangle();
		right_systematic_PM.row_transformation_to_low_triangle();

		//cout << "G=" << G;
		//cout << "Pi=" << Pi;
		//G_or_H.col_permute_to_full_rank_on_left(Pi);
		Matrix<int> nPi = right_systematic_PM.col_permute_to_full_rank_on_right();
		Pi.permute(nPi);
		r_abs.permute(nPi);

		//cout << "G=" << G;
		//cout << "Pi=" << Pi;
		//G_or_H.row_transformation_left_up_triangle_to_identity();
		right_systematic_PM.row_transformation_right_low_triangle_to_identity();

		//cout << "G=" << G;
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

		// by now the first k col of G are most reliable independent positions
		// that is the first k element as index in r, are most reliable independent positions
		// Pi is the record of perbutation

		return Pi;
	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by flipping all possible <=order bits in k most reliable independent position
	 *
	 * \param r: received codeword, after BPSK_demodulation
	 * \param order: order of OSD decoding, all possible bit combinations with number under order, will be test vectors
	 * \return decoded result, u
	 */
	Matrix<GF2> decode_v(Matrix<my_double> recv, int order) {

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
		Matrix<my_double> recv_abs = recv.get_abs();
		Matrix<int> Pi = get_MRIP(recv_abs);	// finally, G_or_H will be G1
		recv.permute(Pi);

		//cout << "Pi=" << Pi << "recv=" << recv;

		//cout << "r=" << r;
		Matrix<GF2> hdr = BPSK::demodulation(recv);
		//cout << "hdr=" << hdr;
		int error_num;

		// first not flip any position, for the first set we do demand in v space
		Matrix<GF2> best_dv = encode_by_right_systematic_PM(hdr);
		//cout << "best_dv: 0" << best_dv;

		// best_dv after permute back must in v space
		total_used_list_num++;
		lambda_min = Measure::correlation_discrepancy_v(recv_abs, hdr, best_dv);
		error_num = Measure::error_num;
		if (optimum_condition_OSD::is_optimum(d_min, k, best_dv, hdr, recv_abs, -1, lambda_min, error_num)) {
			//cout << "early stop, o=0" << endl;
			is_ML = true;
			num_early_stop++;
			
			best_dv.permute_back(Pi);
			return best_dv;
		}

		int current_test_n = 0;
		Matrix<GF2> dv_order_0 = best_dv;
		Matrix<GF2> dv(1, n);
		for (int o = 1; o <= order; ++o) {
			Matrix<int> flipped_mat(1, o, 'N');
			Matrix<int> real_flip_mat(1, o);

			if (o <= 2) {
				do {
					dv = dv_order_0;
					for (int i = 0, imax = flipped_mat.size(); i < imax; ++i) {
						real_flip_mat(i) = k - 1 - flipped_mat(i);
					}
					ref_encode_by_right_systematic_PM(dv, real_flip_mat);

					total_used_list_num++;
					my_double lambda = Measure::correlation_discrepancy_v(recv_abs, hdr, dv);
					error_num = Measure::error_num;
					if (lambda >= lambda_min);
					else {
						// update best_dv to this erasure decoding result
						best_dv = dv;
						lambda_min = lambda;
						//cout << "dv=" << dv; cout << "lambda_min=" << lambda << endl;
						if (optimum_condition_OSD::is_optimum(d_min, k, best_dv, hdr, recv_abs, o - 1, lambda_min, error_num)) {
							//cout << "early stop, o=" << o << endl;
							is_ML = true;
							num_early_stop++;
							best_dv.permute_back(Pi);
							return best_dv;
						}
					}
				} while (flip_TEP::next(flipped_mat, k, o));
			}
			else {
				dv = dv_order_0;
				Matrix<int> last_flipped_mat;
				Matrix<int> fmd(1, 2 * o);
				do {
					//cout << "last_flipped_mat" << last_flipped_mat;
					//cout << "flipped_mat" << flipped_mat;
					flip_TEP::diff(last_flipped_mat, flipped_mat, fmd);
					//cout << "fmd" << fmd;							// in case of order <= 2, to be optimized by storing ordered-zero vector
					
					int fmd_size = fmd.size();
					real_flip_mat.resize(1, fmd_size, false);		// be careful with size of fmd
					for (int i = 0; i < fmd_size; ++i) {
						real_flip_mat(i) = k - 1 - fmd(i);			// to be deleted and replaced by new TEP sequence
					}

					ref_encode_by_right_systematic_PM(dv, real_flip_mat);
					//cout << "dv" << dv;
					last_flipped_mat = flipped_mat;
					
					total_used_list_num++;
					my_double lambda = Measure::correlation_discrepancy_v(recv_abs, hdr, dv);
					error_num = Measure::error_num;
					if (lambda >= lambda_min);
					else {
						// update best_dv to this erasure decoding result
						best_dv = dv;
						lambda_min = lambda;
						//cout << "dv=" << dv; cout << "lambda_min=" << lambda << endl;
						if (optimum_condition_OSD::is_optimum(d_min, k, best_dv, hdr, recv_abs, o - 1, lambda_min, error_num)) {
							//cout << "early stop, o=" << o << endl;
							is_ML = true;
							num_early_stop++;
							best_dv.permute_back(Pi);
							return best_dv;
						}
					}
				} while (flip_TEP::next(flipped_mat, k, o));
			}			
		}

		best_dv.permute_back(Pi);


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

class OSD_v2 {
public:

	Matrix<GF2> A;					// A can be either a generator matrix or a parity check matrix
	bool is_generator;
	Matrix<GF2> systematic_A;		
		// left systematic, for A a generator (parity-check) matrix, column ordered with decreasing (increasing) reliability
	int d_min;
	int n;
	int k;
	int r;
	Matrix<int> Pi;

	int num_early_stop;
	unsigned long long total_used_list_num;		// counting used list num
	bool is_ML;
	my_double lambda_min;

	/**
	 * .constructor of OSD_v2
	 *
	 * \param _A:				the generator matrix or parity-check matrix of a code
	 * 
	 * \param is_generator:		ture -> _A is a generator matrix, 
	 *							false -> _A is a parity check matrix
	 * 
	 * \param _d_min:			the minimum Hamming distance of a code
	 */
	OSD_v2(const Matrix<GF2>& _A, bool _is_generator, int _d_min) :A(_A), is_generator(_is_generator), d_min(_d_min) {
		n = A.col();
		k = _is_generator ? A.row() : n - A.row();
		r = n - k;
		systematic_A.resize(A.row(), n, false);
		Pi.resize(1, n);

		num_early_stop = 0;
		total_used_list_num = 0;
		is_ML = false;
		lambda_min = 0;		// doesn't matter
	}

	// do the same work as the constructor
	void init(const Matrix<GF2>& _A, bool _is_generator, int _d_min) {
		A = _A;
		is_generator = _is_generator;
		d_min = _d_min;		// for ML stopping rule

		n = A.col();
		k = _is_generator ? A.row() : n - A.row();
		r = n - k;
		systematic_A.resize(A.row(), n, false);
		Pi.resize(1, n);

		num_early_stop = 0;
		total_used_list_num = 0;
		is_ML = false;
		lambda_min = 0;		// doesn't matter
	}

	Matrix<GF2> encode_by_systematic_A(const Matrix<GF2>& information_bits) {
		Matrix<GF2> ans(1, n);

		if (is_generator == false) {
			// this is the usual case, OSD deal with the code rate > 1/2, A is a parity check matrix

			for (int i = r; i < n; ++i) {
				ans(i) = information_bits(i);				// k information bits, systematically encoded
			}
			for (int i = 0; i < r; ++i) {
				ans(i) = 0;									// n-k parity check bits
			}

			for (int i = r; i < n; ++i) {
				if (information_bits(i) != 0) {
					for (int j = 0; j < r; ++j) {			// parity check bit, add the column of Patrity check matrix
						ans(j) += systematic_A(j, i);
					}
				}
			}
		}
		else {
			// A is a generator matrix
			
			for (int i = 0; i < k; ++i) {
				ans(i) = information_bits(i);				// k infomation bits, systematically encoded
			}
			for (int i = k; i < n; ++i) {
				ans(i) = 0;									// n-k parity check bits
			}

			for (int i = 0; i < k; ++i) {
				if (information_bits(i) != 0) {
					for (int j = k; j < n; ++j) {			// add the row i of systematic_A to the codeword
						ans(j) += systematic_A(i, j);
					}
				}
			}
		}
		return ans;
	}

	// the input and output will be ref_codeword 
	void ref_encode_by_systematic_A(Matrix<GF2>& ref_codeword, const Matrix<int>& flip_ind) {
		int flip_ind_size = flip_ind.size();
		//cout << "flip_ind_size = " << flip_ind_size << endl;
		if (is_generator == false) {
			for (int i = 0; i < flip_ind_size; ++i) {
				int actual_col = n - 1 - flip_ind(i);	// information bit, systematically encoded, to be check out
				ref_codeword(actual_col) += 1;
				for (int j = 0; j < r; ++j) {			// parity check bit, add the column 'actual_col' of systematic_A to codeword
					ref_codeword(j) += systematic_A(j, actual_col);
				}
			}
		}
		else {
			for (int i = 0; i < flip_ind_size; ++i) {
				int actual_col = flip_ind(i);
				ref_codeword(actual_col) += 1;
				for (int j = k; j < n; ++j) {			// add the row 'actual_col' of systematic_A to the codeword
					ref_codeword(j) += systematic_A(actual_col, j);
				}
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
	void get_MRIP(Matrix<my_double>& r_abs) {
		// include producing systematic G through parity check matrix

		if (is_generator == false) {
			// step 1: turn r to abs(r) and arrange r_abs in increasing order

			Pi = r_abs.sort_with_ind('<');
			//r_abs.permute_back(Pi);		// this process return the original un-permuted r_abs

			// A is a parity check matrix
			systematic_A = A;

			systematic_A.permute_col(Pi);

			Matrix<int> nPi = systematic_A.GE_left_identity_4_GF2();

			Pi.permute(nPi);
			r_abs.permute(nPi);

			// need to sort the latter k position in reliability increasing			
			for (int i = 0; i < n; ++i) {
				nPi(i) = i;										// nPi = Matrix<int>(1, n, 'N');
			}

			for (int i = r; i < n; ++i) {						// for k element in r_abs
				int j = i;
				while (j != r && r_abs(j - 1) > r_abs(j)) {		// bubble sort
					r_abs.switch_ele(j - 1, j);
					nPi.switch_ele(j - 1, j);
					j--;
				}
			}
			Pi.permute(nPi);									// partial permute
			systematic_A.permute_col(nPi);

			// by now the first r col of systematic_A are least reliable independent positions
			// and the last n-k positions are sorted in reliability increasing order
			// Pi is the record of permutation
		}
		else {
			// A is a genarator matrix
			// step 1: turn r to abs(r) and arrange r_abs in decreasing order

			Pi = r_abs.sort_with_ind('>');
			//r_abs.permute_back(Pi);		// this process return the original un-permuted r_abs

			// A is a generator matrix
			systematic_A = A;

			systematic_A.permute_col(Pi);

			Matrix<int> nPi = systematic_A.GE_left_identity_4_GF2();

			Pi.permute(nPi);
			r_abs.permute(nPi);

			// by now the first k col of systematic_A are most reliable independent positions
			// Pi is the record of permutation

			//cout << "systematic_A" << systematic_A;
			//cout << "Pi" << Pi;
		}
	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by flipping all possible <=order bits in k most reliable independent position
	 *
	 * \param r: received codeword, after BPSK_demodulation
	 * \param order: order of OSD decoding, all possible bit combinations with number under order, will be test vectors
	 * \return decoded result, u
	 */
	Matrix<GF2> decode_v(Matrix<my_double> recv, int order) {

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

		//cout << "recv" << recv;
		Matrix<GF2> hdr_orig = BPSK::demodulation(recv);
		//cout << "hdr_orig" << hdr_orig;

		is_ML = false;
		Matrix<my_double> recv_abs = recv.get_abs();
		get_MRIP(recv_abs);	// finally, G_or_H will be G1
		recv.permute(Pi);

		//cout << "Pi=" << Pi << "recv=" << recv;

		//cout << "recv=" << recv;
		Matrix<GF2> hdr = BPSK::demodulation(recv);
		//cout << "hdr=" << hdr;
		int error_num;

		// first not flip any position, for the first set we do demand in v space
		Matrix<GF2> best_dv = encode_by_systematic_A(hdr);
		//cout << "best_dv: 0" << best_dv;

		// best_dv after permute back must in v space
		total_used_list_num++;
		lambda_min = Measure::correlation_discrepancy_v(recv_abs, hdr, best_dv);
		//cout << "lambda_min = " << lambda_min << endl;
		error_num = Measure::error_num;
		if (optimum_condition_OSD::is_optimum_with_sort_way(d_min, k, best_dv, hdr, recv_abs, -1, lambda_min, error_num, is_generator)) {
			//cout << "early stop, o=0" << endl;
			is_ML = true;
			num_early_stop++;

			best_dv.permute_back(Pi);
			return best_dv;
		}

		//best_dv.permute_back(Pi);
		//cout << "best_dv" << best_dv;
		//return best_dv;

		int current_test_n = 0;
		Matrix<GF2> dv_order_0 = best_dv;
		Matrix<GF2> dv(1, n);
		for (int o = 1; o <= order; ++o) {
			Matrix<int> flipped_mat(1, o);
			for (int i = 0; i < o; ++i) {
				flipped_mat(i) = k - o + i;
			}

			while(flipped_mat(0) != -1) {
				dv = dv_order_0;				
				ref_encode_by_systematic_A(dv, flipped_mat);

				//cout << "flipped_mat" << flipped_mat;
				//cout << "dv" << dv;

				total_used_list_num++;
				my_double lambda = Measure::correlation_discrepancy_v(recv_abs, hdr, dv);
				//cout << "lambda = " << lambda << endl;
				error_num = Measure::error_num;
				if (lambda >= lambda_min);
				else {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					//cout << "dv=" << dv; cout << "lambda_min=" << lambda << endl;
					if (optimum_condition_OSD::is_optimum_with_sort_way(d_min, k, best_dv, hdr, \
						recv_abs, o - 1, lambda_min, error_num, is_generator)) {

						//cout << "early stop, o=" << o << endl;
						is_ML = true;
						num_early_stop++;
						//cout << "best_dv" << best_dv;

						best_dv.permute_back(Pi);

						return best_dv;
					}
				}
				flip_TEP::next_v2(flipped_mat, k);						// done
			}
			//cout << "flipped_mat" << flipped_mat;
		}

		//cout << "Pi" << Pi;
		best_dv.permute_back(Pi);

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

class OSD_PSM_aided_GE {
public:

	Matrix<GF2> A;					// A can be either a generator matrix or a parity check matrix
	bool is_generator;
	Matrix<GF2> systematic_A;
	// left systematic, for A a generator (parity-check) matrix, column ordered with decreasing (increasing) reliability
	int d_min;
	int n;
	int k;
	int r;
	PreStored_Matrix_red_equal_dist_fill psm_red_edf;
	Matrix<int> Pi;
	Matrix<int> post_permute;

	int num_early_stop;
	unsigned long long total_used_list_num;		// counting used list num
	bool is_ML;
	my_double lambda_min;

	/**
	 * .constructor of OSD_v2
	 *
	 * \param _A:				the generator matrix or parity-check matrix of a code
	 *
	 * \param is_generator:		ture -> _A is a generator matrix,
	 *							false -> _A is a parity check matrix
	 *
	 * \param _d_min:			the minimum Hamming distance of a code
	 */
	OSD_PSM_aided_GE(const Matrix<GF2>& _A, bool _is_generator, int _partition_num, int _d_min) : \
		A(_A), is_generator(_is_generator), d_min(_d_min), psm_red_edf(_A, _partition_num) {

		n = A.col();
		k = _is_generator ? A.row() : n - A.row();
		r = n - k;
		systematic_A.resize(A.row(), n, false);
		num_early_stop = 0;
		total_used_list_num = 0;
		is_ML = false;
		lambda_min = 0;		// doesn't matter
				 
		Pi.resize(1, n, false);
		post_permute.resize(1, n);
	}

	Matrix<GF2> encode_by_systematic_A(const Matrix<GF2>& information_bits) {
		Matrix<GF2> ans(1, n);

		if (is_generator == false) {
			// this is the usual case, OSD deal with the code rate > 1/2, A is a parity check matrix

			for (int i = r; i < n; ++i) {
				ans(i) = information_bits(i);				// k information bits, systematically encoded
			}
			for (int i = 0; i < r; ++i) {
				ans(i) = 0;									// n-k parity check bits
			}

			for (int i = r; i < n; ++i) {
				if (information_bits(i) != 0) {
					for (int j = 0; j < r; ++j) {			// parity check bit, add the column of Patrity check matrix
						ans(j) += systematic_A(j, i);
					}
				}
			}
		}
		else {
			// A is a generator matrix

			for (int i = 0; i < k; ++i) {
				ans(i) = information_bits(i);				// k infomation bits, systematically encoded
			}
			for (int i = k; i < n; ++i) {
				ans(i) = 0;									// n-k parity check bits
			}

			for (int i = 0; i < k; ++i) {
				if (information_bits(i) != 0) {
					for (int j = k; j < n; ++j) {			// add the row i of systematic_A to the codeword
						ans(j) += systematic_A(i, j);
					}
				}
			}
		}

		//cout << "information_bits" << information_bits;

		//cout << "correct ans" << information_bits.get_part(0, 0, 0, k - 1) * systematic_A;

		//cout << "ans" << ans;		// correct finally

		return ans;
	}

	// the input and output will be ref_codeword 
	void ref_encode_by_systematic_A(Matrix<GF2>& ref_codeword, const Matrix<int>& flip_ind) {
		int flip_ind_size = flip_ind.size();
		//cout << "flip_ind_size = " << flip_ind_size << endl;
		if (is_generator == false) {
			for (int i = 0; i < flip_ind_size; ++i) {
				int actual_col = post_permute(n - 1 - flip_ind(i));	// information bit, systematically encoded, to be check out
				ref_codeword(actual_col) += 1;
				for (int j = 0; j < r; ++j) {			// parity check bit, add the column 'actual_col' of systematic_A to codeword
					ref_codeword(j) += systematic_A(j, actual_col);
				}
			}
		}
		else {
			for (int i = 0; i < flip_ind_size; ++i) {
				int actual_col = post_permute(flip_ind(i));
				ref_codeword(actual_col) += 1;
				for (int j = k; j < n; ++j) {			// add the row 'actual_col' of systematic_A to the codeword
					ref_codeword(j) += systematic_A(actual_col, j);
				}
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
	void get_MRIP(Matrix<my_double>& r_abs) {
		// include producing systematic G through parity check matrix

		if (is_generator == false) {
			// step 1: turn r to abs(r) and arrange r_abs in increasing order

			Matrix<my_double> r_abs_tmp = r_abs;
			Matrix<int> Pi_tmp = r_abs_tmp.sort_with_ind('<');

			psm_red_edf.get_MRIP_sys_G(Pi_tmp, systematic_A, Pi);
				// the k first position of permute_target is not with reliability decreasing order, be careful
				// so as the last n - k position of permute_target

			post_permute.set_natural();		// the ith element stands for the ith most reliable position, in [0,n-k-1] and [n-k,n-1]
			psm_red_edf.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(0, n - k - 1, post_permute);
			psm_red_edf.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(n - k, n - 1, post_permute);

			r_abs.permute(Pi);

			// by now the first r col of systematic_A are least reliable independent positions
			// and the last n-k positions are sorted in reliability increasing order
			// Pi is the record of permutation
		}
		else {
			// A is a genarator matrix
			// step 1: turn r to abs(r) and arrange r_abs in decreasing order

			//cout << "r_abs" << r_abs;

			Matrix<my_double> r_abs_tmp = r_abs;
			Matrix<int> Pi_tmp = r_abs_tmp.sort_with_ind('>');
			//r_abs.permute_back(Pi);		// this process return the original un-permuted r_abs

			psm_red_edf.get_MRIP_sys_G(Pi_tmp, systematic_A, Pi);
				// the k first position of permute_target is not with reliability decreasing order, be careful
				// so as the last n - k position of permute_target

			//cout << "Pi" << Pi;
			//r_abs.permute(Pi);
			//cout << "r_abs (permute by Pi)" << r_abs;

			// by now the first k col of systematic_A are most reliable independent positions
			
			// we should permute systematic_A such that the first k and the last n-k column are in an reliability decreasing order
			
			// this change the elemnts in psm_red_edf, be careful not to rely on this variables any more

			//cout << "psm_red_edf.permute_target_col_ind_to_reliability" << psm_red_edf.permute_target_col_ind_to_reliability;

			post_permute.set_natural();		// the ith element stands for the ith most reliable position, in [0,k-1] and [k,n-1]
			psm_red_edf.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(0, k - 1, post_permute);
			psm_red_edf.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(k, n - 1, post_permute);
			//cout << "post_permute" << post_permute;
						
			// permute the systematic_A, Pi as well as r_abs
			//systematic_A.permute_col(post_permute);				// consider annexing this part to 'PreStored_Matrix_red_equal_dist_fill'

			//cout << "systematic_A" << systematic_A;
			r_abs.permute(Pi);
			//cout << "r_abs (final)" << r_abs;

		}
	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by flipping all possible <=order bits in k most reliable independent position
	 *
	 * \param r: received codeword, after BPSK_demodulation
	 * \param order: order of OSD decoding, all possible bit combinations with number under order, will be test vectors
	 * \return decoded result, u
	 */
	Matrix<GF2> decode_v(Matrix<my_double> recv, int order, Matrix<GF2> PM) {

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

		//cout << "recv" << recv;
		Matrix<GF2> hdr_orig = BPSK::demodulation(recv);
		//cout << "hdr_orig" << hdr_orig;

		is_ML = false;
		Matrix<my_double> recv_abs = recv.get_abs();
		get_MRIP(recv_abs);	// finally, G_or_H will be G1
		recv.permute(Pi);
		//cout << "Pi=" << Pi;

		// check out systematic_A
		//systematic_A.permute_col_back(Pi);
		//cout << "systematic_A.multiply_transpose_of(PM)" << systematic_A.multiply_transpose_of(PM);		// right

		//cout << "recv=" << recv;
		Matrix<GF2> hdr = BPSK::demodulation(recv);
		//cout << "hdr=" << hdr;
		int error_num;

		// first not flip any position, for the first set we do demand in v space
		Matrix<GF2> best_dv = encode_by_systematic_A(hdr);
		//cout << "best_dv: 0" << best_dv;

		// best_dv after permute back must in v space
		total_used_list_num++;
		lambda_min = Measure::correlation_discrepancy_v(recv_abs, hdr, best_dv);
		error_num = Measure::error_num;
		if (optimum_condition_OSD::is_optimum_with_sort_way_pp(d_min, k, best_dv, hdr, \
			recv_abs, -1, lambda_min, error_num, is_generator, post_permute)) {

			//cout << "early stop, o=0" << endl;
			is_ML = true;
			num_early_stop++;

			best_dv.permute_back(Pi);
			return best_dv;
		}

		//best_dv.permute_back(Pi);
		//cout << "best_dv (order 0, permute back)" << best_dv;
		//return best_dv;

		int current_test_n = 0;
		Matrix<GF2> dv_order_0 = best_dv;
		Matrix<GF2> dv(1, n);
		for (int o = 1; o <= order; ++o) {
			Matrix<int> flipped_mat(1, o);
			for (int i = 0; i < o; ++i) {
				flipped_mat(i) = k - o + i;
			}

			while (flipped_mat(0) != -1) {
				dv = dv_order_0;
				ref_encode_by_systematic_A(dv, flipped_mat);

				//cout << "flipped_mat" << flipped_mat;
				//cout << "dv" << dv;

				total_used_list_num++;
				my_double lambda = Measure::correlation_discrepancy_v(recv_abs, hdr, dv);		// why it costs more double operation
				error_num = Measure::error_num;
				if (lambda >= lambda_min);
				else {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					//cout << "dv=" << dv; cout << "lambda_min=" << lambda << endl;
					if (optimum_condition_OSD::is_optimum_with_sort_way_pp(d_min, k, best_dv, hdr, \
						recv_abs, o - 1, lambda_min, error_num, is_generator, post_permute)) {

						//cout << "early stop, o=" << o << endl;
						is_ML = true;
						num_early_stop++;
						best_dv.permute_back(Pi);
						return best_dv;
					}
				}
				flip_TEP::next_v2(flipped_mat, k);						// done
			}
			//cout << "flipped_mat" << flipped_mat;
		}

		//cout << "Pi" << Pi;
		best_dv.permute_back(Pi);

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

class OSD_PSM_aided_GE_4_cyc {
public:

	Matrix<GF2> A;					// A can be either a generator matrix or a parity check matrix
	bool is_generator;
	Matrix<GF2> systematic_A;
	// left systematic, for A a generator (parity-check) matrix, column ordered with decreasing (increasing) reliability
	int d_min;
	int n;
	int k;
	int r;
	PreStored_Matrix_red_equal_dist_cycle psm_red_ed_cyc;		// only for cycle code, like BCH code, eBCH code is not okay
	Matrix<int> Pi;
	Matrix<int> post_permute;

	int num_early_stop;
	unsigned long long total_used_list_num;		// counting used list num
	bool is_ML;
	my_double lambda_min;

	/**
	 * .constructor of OSD_v2
	 *
	 * \param _A:				the generator matrix or parity-check matrix of a code
	 *
	 * \param is_generator:		ture -> _A is a generator matrix,
	 *							false -> _A is a parity check matrix
	 *
	 * \param _d_min:			the minimum Hamming distance of a code
	 */
	OSD_PSM_aided_GE_4_cyc(const Matrix<GF2>& _A, bool _is_generator, int _partition_num, int _d_min) : \
		A(_A), is_generator(_is_generator), d_min(_d_min), psm_red_ed_cyc(_A,_partition_num) {

		n = A.col();
		k = _is_generator ? A.row() : n - A.row();
		r = n - k;
		systematic_A.resize(A.row(), n, false);
		num_early_stop = 0;
		total_used_list_num = 0;
		is_ML = false;
		lambda_min = 0;		// doesn't matter

		Pi.resize(1, n, false);
		post_permute.resize(1, n);
	}

	Matrix<GF2> encode_by_systematic_A(const Matrix<GF2>& information_bits) {
		Matrix<GF2> ans(1, n);

		if (is_generator == false) {
			// this is the usual case, OSD deal with the code rate > 1/2, A is a parity check matrix

			for (int i = r; i < n; ++i) {
				ans(i) = information_bits(i);				// k information bits, systematically encoded
			}
			for (int i = 0; i < r; ++i) {
				ans(i) = 0;									// n-k parity check bits
			}

			for (int i = r; i < n; ++i) {
				if (information_bits(i) != 0) {
					for (int j = 0; j < r; ++j) {			// parity check bit, add the column of Patrity check matrix
						ans(j) += systematic_A(j, i);
					}
				}
			}
		}
		else {
			// A is a generator matrix

			for (int i = 0; i < k; ++i) {
				ans(i) = information_bits(i);				// k infomation bits, systematically encoded
			}
			for (int i = k; i < n; ++i) {
				ans(i) = 0;									// n-k parity check bits
			}

			for (int i = 0; i < k; ++i) {
				if (information_bits(i) != 0) {
					for (int j = k; j < n; ++j) {			// add the row i of systematic_A to the codeword
						ans(j) += systematic_A(i, j);
					}
				}
			}
		}

		//cout << "information_bits" << information_bits;

		//cout << "correct ans" << information_bits.get_part(0, 0, 0, k - 1) * systematic_A;

		//cout << "ans" << ans;		// correct finally

		return ans;
	}

	// the input and output will be ref_codeword 
	void ref_encode_by_systematic_A(Matrix<GF2>& ref_codeword, const Matrix<int>& flip_ind) {
		int flip_ind_size = flip_ind.size();
		//cout << "flip_ind_size = " << flip_ind_size << endl;
		if (is_generator == false) {
			for (int i = 0; i < flip_ind_size; ++i) {
				int actual_col = post_permute(n - 1 - flip_ind(i));	// information bit, systematically encoded, to be check out
				ref_codeword(actual_col) += 1;
				for (int j = 0; j < r; ++j) {			// parity check bit, add the column 'actual_col' of systematic_A to codeword
					ref_codeword(j) += systematic_A(j, actual_col);
				}
			}
		}
		else {
			for (int i = 0; i < flip_ind_size; ++i) {
				int actual_col = post_permute(flip_ind(i));
				ref_codeword(actual_col) += 1;
				for (int j = k; j < n; ++j) {			// add the row 'actual_col' of systematic_A to the codeword
					ref_codeword(j) += systematic_A(actual_col, j);
				}
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
	void get_MRIP(Matrix<my_double>& r_abs) {
		// include producing systematic G through parity check matrix

		// A is a genarator matrix
		// step 1: turn r to abs(r) and arrange r_abs in decreasing order

		//cout << "r_abs" << r_abs;

		Matrix<my_double> r_abs_tmp = r_abs;
		Matrix<int> Pi_tmp = r_abs_tmp.sort_with_ind(is_generator ? '>' : '<');
		//r_abs.permute_back(Pi);		// this process return the original un-permuted r_abs

		psm_red_ed_cyc.get_MRIP_sys_G(Pi_tmp, systematic_A, Pi);
		// the k first position of permute_target is not with reliability decreasing order, be careful
		// so as the last n - k position of permute_target

		//cout << "Pi" << Pi;
		//r_abs.permute(Pi);
		//cout << "r_abs (permute by Pi)" << r_abs;

		// by now the first k col of systematic_A are most reliable independent positions

		// we should permute systematic_A such that the first k and the last n-k column are in an reliability decreasing order

		// this change the elemnts in psm_red_edf, be careful not to rely on this variables any more

		//cout << "psm_red_edf.permute_target_col_ind_to_reliability" << psm_red_edf.permute_target_col_ind_to_reliability;

		post_permute.set_natural();		// the ith element stands for the ith most reliable position, in [0,k-1] and [k,n-1]
		if (is_generator == false) {
			psm_red_ed_cyc.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(0, n - k - 1, post_permute);
			psm_red_ed_cyc.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(n - k, n - 1, post_permute);
		}
		else {
			psm_red_ed_cyc.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(0, k - 1, post_permute);
			psm_red_ed_cyc.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(k, n - 1, post_permute);
		}
		//cout << "post_permute" << post_permute;

		// permute the systematic_A, Pi as well as r_abs
		//systematic_A.permute_col(post_permute);				// consider annexing this part to 'PreStored_Matrix_red_equal_dist_fill'

		//cout << "systematic_A" << systematic_A;
		r_abs.permute(Pi);
		//cout << "r_abs (final)" << r_abs;

	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by flipping all possible <=order bits in k most reliable independent position
	 *
	 * \param r: received codeword, after BPSK_demodulation
	 * \param order: order of OSD decoding, all possible bit combinations with number under order, will be test vectors
	 * \return decoded result, u
	 */
	Matrix<GF2> decode_v(Matrix<my_double> recv, int order, Matrix<GF2> PM) {

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

		//cout << "recv" << recv;
		Matrix<GF2> hdr_orig = BPSK::demodulation(recv);
		//cout << "hdr_orig" << hdr_orig;

		is_ML = false;
		Matrix<my_double> recv_abs = recv.get_abs();
		get_MRIP(recv_abs);	// finally, G_or_H will be G1
		recv.permute(Pi);
		//cout << "Pi=" << Pi;

		// check out systematic_A
		//systematic_A.permute_col_back(Pi);
		//cout << "systematic_A.multiply_transpose_of(PM)" << systematic_A.multiply_transpose_of(PM);		// right

		//cout << "recv=" << recv;
		Matrix<GF2> hdr = BPSK::demodulation(recv);
		//cout << "hdr=" << hdr;
		int error_num;

		// first not flip any position, for the first set we do demand in v space
		Matrix<GF2> best_dv = encode_by_systematic_A(hdr);
		//cout << "best_dv: 0" << best_dv;

		// best_dv after permute back must in v space
		total_used_list_num++;
		lambda_min = Measure::correlation_discrepancy_v(recv_abs, hdr, best_dv);
		error_num = Measure::error_num;
		if (optimum_condition_OSD::is_optimum_with_sort_way_pp(d_min, k, best_dv, hdr, \
			recv_abs, -1, lambda_min, error_num, is_generator, post_permute)) {

			//cout << "early stop, o=0" << endl;
			is_ML = true;
			num_early_stop++;

			best_dv.permute_back(Pi);
			return best_dv;
		}

		//best_dv.permute_back(Pi);
		//cout << "best_dv (order 0, permute back)" << best_dv;		// problem to be solved, should be ended in order 0
		//return best_dv;

		int current_test_n = 0;
		Matrix<GF2> dv_order_0 = best_dv;
		Matrix<GF2> dv(1, n);
		for (int o = 1; o <= order; ++o) {
			Matrix<int> flipped_mat(1, o);
			for (int i = 0; i < o; ++i) {
				flipped_mat(i) = k - o + i;
			}

			while (flipped_mat(0) != -1) {
				dv = dv_order_0;
				ref_encode_by_systematic_A(dv, flipped_mat);

				//cout << "flipped_mat" << flipped_mat;
				//cout << "dv" << dv;

				total_used_list_num++;
				my_double lambda = Measure::correlation_discrepancy_v(recv_abs, hdr, dv);		// why it costs more double operation
				error_num = Measure::error_num;
				if (lambda >= lambda_min);
				else {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					//cout << "dv=" << dv; cout << "lambda_min=" << lambda << endl;
					if (optimum_condition_OSD::is_optimum_with_sort_way_pp(d_min, k, best_dv, hdr, \
						recv_abs, o - 1, lambda_min, error_num, is_generator, post_permute)) {

						//cout << "early stop, o=" << o << endl;
						is_ML = true;
						num_early_stop++;
						best_dv.permute_back(Pi);
						return best_dv;
					}
				}
				flip_TEP::next_v2(flipped_mat, k);						// done
			}
			//cout << "flipped_mat" << flipped_mat;
		}

		//cout << "Pi" << Pi;
		best_dv.permute_back(Pi);

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

		//cout << "best_dv" << best_dv;
		return best_dv;
	}
};

class OSD_PSM_aided_GE_4_ext_cyc {
public:

	Matrix<GF2> A;					// A can be either a generator matrix or a parity check matrix
	bool is_generator;
	Matrix<GF2> systematic_A;
	// left systematic, for A a generator (parity-check) matrix, column ordered with decreasing (increasing) reliability
	int d_min;
	int n;
	int k;
	int r;
	PreStored_Matrix_red_equal_dist_extend_cycle psm_red_ed_ext_cyc;		// only for cycle code, like BCH code, eBCH code is not okay
	Matrix<int> Pi;
	Matrix<int> post_permute;

	int num_early_stop;
	unsigned long long total_used_list_num;		// counting used list num
	bool is_ML;
	my_double lambda_min;

	/**
	 * .constructor of OSD_v2
	 *
	 * \param _A:				the generator matrix or parity-check matrix of a code
	 *
	 * \param is_generator:		ture -> _A is a generator matrix,
	 *							false -> _A is a parity check matrix
	 *
	 * \param _d_min:			the minimum Hamming distance of a code
	 */
	OSD_PSM_aided_GE_4_ext_cyc(const Matrix<GF2>& _A, bool _is_generator,int _partition_num, int _d_min) : \
		A(_A), is_generator(_is_generator), d_min(_d_min), psm_red_ed_ext_cyc(_A, _is_generator, _partition_num) {

		n = A.col();
		k = _is_generator ? A.row() : n - A.row();
		r = n - k;
		systematic_A.resize(A.row(), n, false);
		num_early_stop = 0;
		total_used_list_num = 0;
		is_ML = false;
		lambda_min = 0;		// doesn't matter

		Pi.resize(1, n, false);
		post_permute.resize(1, n);
	}

	Matrix<GF2> encode_by_systematic_A(const Matrix<GF2>& information_bits) {
		Matrix<GF2> ans(1, n);
		if (is_generator == false) {
			// this is the usual case, OSD deal with the code rate > 1/2, A is a parity check matrix

			for (int i = r; i < n; ++i) {
				ans(i) = information_bits(i);				// k information bits, systematically encoded
			}
			for (int i = 0; i < r; ++i) {
				ans(i) = 0;									// n-k parity check bits
			}

			for (int i = r; i < n; ++i) {
				if (information_bits(i) != 0) {
					for (int j = 0; j < r; ++j) {			// parity check bit, add the column of Patrity check matrix
						ans(j) += systematic_A(j, i);
					}
				}
			}
		}
		else {
			// A is a generator matrix

			for (int i = 0; i < k; ++i) {
				ans(i) = information_bits(i);				// k infomation bits, systematically encoded
			}
			for (int i = k; i < n; ++i) {
				ans(i) = 0;									// n-k parity check bits
			}

			for (int i = 0; i < k; ++i) {
				if (information_bits(i) != 0) {
					for (int j = k; j < n; ++j) {			// add the row i of systematic_A to the codeword
						ans(j) += systematic_A(i, j);
					}
				}
			}
		}

		//cout << "information_bits" << information_bits;

		//cout << "correct ans" << information_bits.get_part(0, 0, 0, k - 1) * systematic_A;

		//cout << "ans" << ans;		// correct finally
		return ans;
	}

	// the input and output will be ref_codeword 
	void ref_encode_by_systematic_A(Matrix<GF2>& ref_codeword, const Matrix<int>& flip_ind) {
		int flip_ind_size = flip_ind.size();
		//cout << "flip_ind_size = " << flip_ind_size << endl;
		if (is_generator == false) {
			for (int i = 0; i < flip_ind_size; ++i) {
				int actual_col = post_permute(n - 1 - flip_ind(i));	// information bit, systematically encoded, to be check out
				ref_codeword(actual_col) += 1;
				for (int j = 0; j < r; ++j) {			// parity check bit, add the column 'actual_col' of systematic_A to codeword
					ref_codeword(j) += systematic_A(j, actual_col);
				}
			}
		}
		else {
			for (int i = 0; i < flip_ind_size; ++i) {
				int actual_col = post_permute(flip_ind(i));
				ref_codeword(actual_col) += 1;
				for (int j = k; j < n; ++j) {			// add the row 'actual_col' of systematic_A to the codeword
					ref_codeword(j) += systematic_A(actual_col, j);
				}
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
	void get_MRIP(Matrix<my_double>& r_abs) {
		// include producing systematic G through parity check matrix

		
		// A is a genarator matrix
		// step 1: turn r to abs(r) and arrange r_abs in decreasing order

		//cout << "r_abs" << r_abs;

		Matrix<my_double> r_abs_tmp = r_abs;
		Matrix<int> Pi_tmp = r_abs_tmp.sort_with_ind(is_generator ? '>' : '<');
		//r_abs.permute_back(Pi);		// this process return the original un-permuted r_abs

		psm_red_ed_ext_cyc.get_MRIP_sys_G(Pi_tmp, systematic_A, Pi);
		// the k first position of permute_target is not with reliability decreasing order, be careful
		// so as the last n - k position of permute_target

		//cout << "Pi" << Pi;
		//r_abs.permute(Pi);
		//cout << "r_abs (permute by Pi)" << r_abs;

		// by now the first k col of systematic_A are most reliable independent positions

		// we should permute systematic_A such that the first k and the last n-k column are in an reliability decreasing order

		// this change the elemnts in psm_red_edf, be careful not to rely on this variables any more

		//cout << "psm_red_edf.permute_target_col_ind_to_reliability" << psm_red_edf.permute_target_col_ind_to_reliability;

		post_permute.set_natural();		// the ith element stands for the ith most reliable position, in [0,k-1] and [k,n-1]
		if (is_generator == false) {
			psm_red_ed_ext_cyc.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(0, n - k - 1, post_permute);
			psm_red_ed_ext_cyc.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(n - k, n - 1, post_permute);
		}
		else {
			psm_red_ed_ext_cyc.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(0, k - 1, post_permute);
			psm_red_ed_ext_cyc.BU.permute_target_col_ind_to_reliability.quick_sort_recur_lt_with_ind(k, n - 1, post_permute);
		}
		//cout << "post_permute" << post_permute;

		// permute the systematic_A, Pi as well as r_abs
		//systematic_A.permute_col(post_permute);				// consider annexing this part to 'PreStored_Matrix_red_equal_dist_fill'

		//cout << "systematic_A" << systematic_A;
		r_abs.permute(Pi);
		//cout << "r_abs (final)" << r_abs;
	}

	/**
	 * .apply decode method of 'code_instance' to 'r',
	 *		by flipping all possible <=order bits in k most reliable independent position
	 *
	 * \param r: received codeword, after BPSK_demodulation
	 * \param order: order of OSD decoding, all possible bit combinations with number under order, will be test vectors
	 * \return decoded result, u
	 */
	Matrix<GF2> decode_v(Matrix<my_double> recv, int order, Matrix<GF2> PM) {

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

		//cout << "recv" << recv;
		Matrix<GF2> hdr_orig = BPSK::demodulation(recv);
		//cout << "hdr_orig" << hdr_orig;

		is_ML = false;
		Matrix<my_double> recv_abs = recv.get_abs();
		get_MRIP(recv_abs);	// finally, G_or_H will be G1
		recv.permute(Pi);
		//cout << "Pi=" << Pi;

		// check out systematic_A
		//systematic_A.permute_col_back(Pi);
		//cout << "systematic_A.multiply_transpose_of(PM)" << systematic_A.multiply_transpose_of(PM);		// right

		//cout << "recv=" << recv;
		Matrix<GF2> hdr = BPSK::demodulation(recv);
		//cout << "hdr=" << hdr;
		int error_num;

		// first not flip any position, for the first set we do demand in v space
		Matrix<GF2> best_dv = encode_by_systematic_A(hdr);
		//cout << "best_dv: 0" << best_dv;

		// best_dv after permute back must in v space
		total_used_list_num++;
		lambda_min = Measure::correlation_discrepancy_v(recv_abs, hdr, best_dv);
		error_num = Measure::error_num;
		if (optimum_condition_OSD::is_optimum_with_sort_way_pp(d_min, k, best_dv, hdr, \
			recv_abs, -1, lambda_min, error_num, is_generator, post_permute)) {

			//cout << "early stop, o=0" << endl;
			is_ML = true;
			num_early_stop++;

			best_dv.permute_back(Pi);
			return best_dv;
		}

		//best_dv.permute_back(Pi);
		//cout << "best_dv (order 0, permute back)" << best_dv;		// problem to be solved, should be ended in order 0
		//return best_dv;

		int current_test_n = 0;
		Matrix<GF2> dv_order_0 = best_dv;
		Matrix<GF2> dv(1, n);
		for (int o = 1; o <= order; ++o) {
			Matrix<int> flipped_mat(1, o);
			for (int i = 0; i < o; ++i) {
				flipped_mat(i) = k - o + i;
			}

			while (flipped_mat(0) != -1) {
				dv = dv_order_0;
				ref_encode_by_systematic_A(dv, flipped_mat);

				//cout << "flipped_mat" << flipped_mat;
				//cout << "dv" << dv;

				total_used_list_num++;
				my_double lambda = Measure::correlation_discrepancy_v(recv_abs, hdr, dv);		// why it costs more double operation
				error_num = Measure::error_num;
				if (lambda >= lambda_min);
				else {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					//cout << "dv=" << dv; cout << "lambda_min=" << lambda << endl;
					if (optimum_condition_OSD::is_optimum_with_sort_way_pp(d_min, k, best_dv, hdr, \
						recv_abs, o - 1, lambda_min, error_num, is_generator, post_permute)) {

						//cout << "early stop, o=" << o << endl;
						is_ML = true;
						num_early_stop++;
						best_dv.permute_back(Pi);
						return best_dv;
					}
				}
				flip_TEP::next_v2(flipped_mat, k);						// done
			}
			//cout << "flipped_mat" << flipped_mat;
		}

		//cout << "Pi" << Pi;
		best_dv.permute_back(Pi);

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

		//cout << "best_dv" << best_dv;
		return best_dv;
	}
};