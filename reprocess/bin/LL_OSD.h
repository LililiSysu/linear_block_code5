#pragma once
/*****************************************************************//**
 * \file   Low_Latency_OSD.h
 * \brief  idea from Lijia's 2022 ITW paper "Low-lantency ordered statistics decoding of BCH codes"
 *			use the fact that BCH code is sub-field sub-code of RS code, and Generator Matrix of RS code
 *			can be construct by Lagrange polynomials, saving time for Gaussian elemination
 *
 * \author 26259
 * \date   March 2023
 *********************************************************************/
#include "../code/BCH.h"
#include "../code/RS.h"
#include"../GF/GF2.h"
#include"../my_lib/matrix.h"
#include"../channel/channel.h"

#ifdef use_my_double
#include"../my_lib/my_double.h"
#else
#define my_double double
#endif

/**
 * .special for BCH code, m is the super field's power, t is super error correcting capability of bch code
 * k_prime is 2^m - 1 - 2 * t, please compute it and set it a constant number before using the class
 */
template<int m,int t,int k_prime> class Low_Latency_OSD {
public:
	BCH<m, t> bch;
	RS<m, k_prime> super_code_rs;

	/**
	 * .a warp function of decode_v, see decode_v
	 */
	Matrix<GF2> decode(Matrix<my_double> r, int order) {
		return bch.v2u(decode_v(r, order));
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
	Matrix<GF2> decode_v(Matrix<my_double> r, int order) {
		// a shell for code_type, such as BCH, return code_instance.decode(r);

		int n = r.size();	// also known as 
		int k = bch.get_k();

		Matrix<my_double> r_abs(1, n);
		for (int i = 0; i < n; ++i) {
			r_abs(i) = my::abs(r(i));
		}

		Matrix<GF2> hdr = BPSK::demodulation<GF2>(r);
		//cout << "r=" << r;
		//cout << "hdr=" << hdr;
		Matrix<int> Pi = r_abs.sort_with_ind();

		int error_num;
		int d_min = bch.get_d();
		//order = 0;
		//cout << "order=" << order << endl;

		Matrix<int> information_set_ind(1, k_prime);
		for (int i = 0; i < k_prime; ++i) {
			information_set_ind(i) = Pi(i);
		}
		information_set_ind.sort();


		super_code_rs.generate_systematic_generator_any_pos(information_set_ind);

		Matrix<GF2e<m>> super_information_set(1, k_prime);
		for (int i = 0; i < k_prime; ++i) { 
			super_information_set(i) = hdr(information_set_ind(i));
		}
		//cout << "super_information_set" << super_information_set;

		Matrix<GF2e<m>> best_dv(1, n, '0');
		// check if encoded result is binary
		bool is_valid_codeword = true;
		for (int i = 0; i < n - k; ++i) {
			int result_ind = super_code_rs.redundancy_set_ind(i);
			for (int j = 0; j < k; ++j) {
				best_dv(result_ind) += super_information_set(j) * super_code_rs.generator_M_systematic_any_pos(j, result_ind);
			}
			if (best_dv(result_ind) != 0 && best_dv(result_ind) != 1) {
				is_valid_codeword = false;
				break;
			}
		}
		if (is_valid_codeword) {
			for (int i = 0; i < k; ++i) {
				best_dv(information_set_ind(i)) = super_information_set(i);
			}
		}

		//cout << "best_dv" << best_dv;		// same as BCH code

		// first not flip any position, for the first set we do demand in v space
		//cout << "hdr" << hdr;


		my_double lambda_min = Measure::correlation_discrepancy_BPSK_for_v(r_abs, hdr, best_dv, error_num);
		if (is_valid_codeword && optimum_condition_OSD::is_optimum(d_min, k, best_dv, hdr, r_abs, -1, lambda_min, error_num)) {
			//cout << "early stop, o=0" << endl;
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
			do {
				//cout << "flipped_mat" << flipped_mat;
				flip_mat_discrepancy(last_flipped_mat, flipped_mat, fmd, fmd_size);
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

				//if (!code_instance.is_in_v_space(dv.get_permute_back(Pi)))	cout << "error" << endl; // for later dv, we demand it in v space
				my_double lambda = Measure::correlation_discrepancy_BPSK_for_v(r_abs, hdr, dv, error_num);
				if (lambda >= lambda_min);
				else {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					//cout << "dv=" << dv; cout << "lambda_min=" << lambda << endl;
					if (optimum_condition_OSD::is_optimum(d_min, k, best_dv, hdr, r_abs, o - 1, lambda_min, error_num)) {
						//cout << "early stop, o=" << o << endl;
						best_dv.permute_back(Pi);
						return best_dv;
					}
				}
			} while (next_flip(flipped_mat, k, o));
		}

		best_dv.permute_back(Pi);
		return best_dv;
	}

	/**
	 * .find the next bit flip pattern in OSD
	 *
	 * \param is_flipped_mat: current flip pattern
	 * \param k: information bit numbers
	 * \param order: order in OSD, aka. number of flip pattern
	 * \return whether travelled all flipp pattern
	 */
	bool next_flip(Matrix<int>& is_flipped_mat, int k, int order) {
		bool need_update_former;
		int former_ind = 0;
		do {
			former_ind++;
			is_flipped_mat(order - former_ind)++;
			need_update_former = (is_flipped_mat(order - former_ind) == (k - former_ind + 1));
		} while (need_update_former && former_ind != order);
		if (!(need_update_former && former_ind == order)) {
			for (int i = former_ind - 1; i > 0; --i) {
				is_flipped_mat(order - i) = is_flipped_mat(order - i - 1) + 1;
			}
			return true;
		}
		else
			return false;
	}

	/**
	 * .compute difference between 'last_flipped_mat' and 'flipped_mat'
	 *
	 * \param last_flipped_mat: the last flipped pattern
	 * \param flipped_mat: current flip pattern
	 * \param fmd: this will be set as difference between 'last_flipped_mat' and 'flipped_mat'
	 * ie. 'last_flipped_mat'={2,3}, and 'flipped_mat'={2,4}, then 'fmd'={3,4}
	 * \param fmd_size: valid size of fmd, in case above this will be set as 2
	 */
	void flip_mat_discrepancy(Matrix<int>& last_flipped_mat, \
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

