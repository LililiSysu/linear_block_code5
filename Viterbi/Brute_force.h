/*****************************************************************//**
 * \file   Brute_force.h
 * \brief  brute force search to decode any linear block code
 * 
 * \author 26259
 * \date   April 2023
 *********************************************************************/

#pragma once

#include "Viterbi_common.h"

class Brute_force {		// this is also very important, for short linear block code, only valid for k < 31
private:

	int n;
	int k;
	Matrix<GF2> GM;
	int pow_2_k;
	Matrix<my_double> dist;// initialize correlation distance, considering both case for soft metric and hard metric

	Matrix<GF2> codeword_space;				// for brutle force

public:

	Brute_force(const Matrix<GF2>& _GM) {
		// get the code parameters, for a linear block code of (n,k) and r = n - k, 
		// with parity check matrix PM, generator check matrix GM

		GM = _GM;
		n = GM.col();
		k = GM.row();

		pow_2_k = 1 << k;
		dist = Matrix<my_double>(2, n);
	}

	/**
	 * .brutle force decoding
	 */
	template<class T>
	Matrix<GF2> decode_v(const Matrix<T>& r_or_hdr, int list_size = 1) {

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

		// to merge the method of hard and soft decoding, if T is ordered structure, then use soft decoidng, else use hard deocding
		bool is_soft = !(T() < T());
		//cout << "is_soft = " << is_soft << endl;

		if (codeword_space.row() != pow_2_k) {

			codeword_space = Matrix<GF2>(pow_2_k, n);
			for (int j = 0; j < n; ++j) {
				codeword_space(0, j) = 0;
			}

			for (int i = 1; i < pow_2_k; ++i) {
				// decompose i into two small weight numbers, by taking out the least significant digit 1

				int least_1_pos = 0;
				for (int store_i = i; (store_i & 1) == 0; store_i >>= 1) {
					least_1_pos++;
				}
				int seperate_num = 1 << least_1_pos;
				int seperate_num2 = i - seperate_num;
				if (seperate_num2 == 0) {
					for (int j = 0; j < n; ++j) {
						codeword_space(i, j) = GM(least_1_pos, j);
					}
				}
				else {
					// linear combination of 2 row in codeword space
					for (int j = 0; j < n; ++j) {
						codeword_space(i, j) = codeword_space(seperate_num, j) + codeword_space(seperate_num2, j);
					}
				}
			}
		}

		//cout << "codeword_space" << codeword_space;
		// compute distance for each symbol to 0 and 1
		for (int i = 0; i < n; ++i) {
			// not key complexity
			T ri = r_or_hdr(i);

			my_double d0 = (double)(is_soft ? (ri - 1) : (ri != 0));
			my_double d1 = (double)(is_soft ? (ri + 1) : (ri == 0));

			dist(0, i) = d0 * d0;
			dist(1, i) = d1 * d1;
		}

		//ope_num_after = my_double_auxiliary_storage::operation_number;
		//cout << "ope_num = " << ope_num_after - ope_num_before << endl;

		// compute each codewords metric
		Matrix<my_double> metric(pow_2_k, 1, '0');
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < pow_2_k; ++i) {
				metric(i) += codeword_space(i, j) == 1 ? dist(1, j) : dist(0, j);
			}
		}
		//cout << "codeword_space" << codeword_space;
		//cout << "metric" << metric;

		// all sort  --- (down)
		Matrix<int> permute_ind = metric.sort_with_ind();
		codeword_space.permute_row(permute_ind);

		// all sort  --- (up)


#ifdef RUN_MSG
		cout << "codeword_list" << codeword_space.get_part(0, 0, list_size - 1, n - 1);
		cout << "metric" << metric.get_part(0, 0, list_size - 1, 0);	// right sorting

#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(Brute_force) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Brute_force) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return codeword_space.get_part(0, 0, list_size - 1, n - 1);
	}
};
