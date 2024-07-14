#pragma once
/*****************************************************************//**
 * \file   reprocess_common.h
 * \brief  common class for reprocess function
 * 
 * \author lilili
 * \date   April 2023
 *********************************************************************/
#include"../GF/GF2.h"
#include"../my_lib/Matrix.h"
#include"../my_lib/Matrix_common.h"
#include"../my_lib/Vector_static.h"
#include"../channel/channel.h"

 // but we use my_double, not my_float in these class, to be fixed, but not important
#ifdef use_my_double
#include"../my_lib/my_double.h"
#include"../my_lib/my_float.h"
#else
#define my_double double
#define my_float float
#endif

class optimum_condition
{
public:
	/**
	 * .judge whether current codeword 'dv' is optimum
	 *
	 * \param d_min: minimum hamming distance of code
	 * \param dv: codeword to judge, length of n
	 * \param hdr: hard decoding result of received word
	 * \param r_abs: absolote value of received word, indicating the reialbility
	 * \param r_ind: sorted index indicating r_abs from small to big, i.e., r_abs(r_ind(0)) has smallest value in r_abs
	 * \param lambda: correlation discripency between 'dv' and received word
	 * \param error_num: number of difference between 'dv' and 'hdr'
	 * \return  whether 'dv' is an optimum codeword
	 */
	static bool is_optimum(int d_min, const Matrix<GF2>& dv, const Matrix<GF2>& hdr, \
		const Matrix<my_double>& r_abs, const Matrix<int>& r_ind, my_double lambda, int error_num) {

		int n = dv.size();
		my_double G_vw1 = 0;
		// consider j=1 only, this judge if is optimum
		int q1 = d_min - error_num;		// d_min
		int cnt = 0;
		for (int i = 0; i < n && cnt < q1; ++i) {
			if (dv(r_ind(i)) == hdr(r_ind(i))) {
				G_vw1 += r_abs(i);
				cnt++;
			}
		}

		//cout << "G_vw1=" << G_vw1 << endl;
		return lambda <= G_vw1;
	}
};

class optimum_condition_OSD
{
public:
	/**
	 * .judge whether current codeword 'dv' is optimum, special for OSD use only
	 *
	 * \param d_min: minimum hamming distance of code
	 * \param k: information bit number
	 * \param dv: codeword to judge, length of n
	 * \param hdr: hard decoding result of received word
	 * \param r_abs: absolote value of received word, indicating the reialbility
	 * \param finished_order: the order OSD that has been finished, set it be current order - 1
	 * \param lambda: correlation discripency between 'dv' and received word
	 * \param error_num: number of difference between 'dv' and 'hdr'
	 * \return whether 'dv' is an optimum codeword
	 */
	static bool is_optimum(int d_min, int k, const Matrix<GF2>& dv, const Matrix<GF2>& hdr, \
		const Matrix<my_double>& r_abs, int finished_order, my_double lambda, int error_num) {	// this is designed for OSD

		//return false;			// for worse case complexity testing

		int n = dv.size();
		my_double G_vw1 = 0;
		// consider j=1 only, this judge if is optimum
		int q1 = d_min - error_num - (finished_order + 1);		// d_min
		int cnt = 0;
		for (int i = n - 1; i >= 0 && cnt < q1; --i) {
			if (dv(i) == hdr(i)) {
				G_vw1 += r_abs(i);		// it seems the optimum condition has no relation to the reliability ordering shortly after k
				cnt++;
			}
		}
		for (int i = 1; i <= finished_order + 1; ++i) {
			G_vw1 += r_abs(k - i);		// this has relation to the reliability ordering shortly before k
		}

		//cout << "G_vw1=" << G_vw1 << endl;
		return lambda <= G_vw1;
	}

	/**
	 * .judge whether current codeword 'dv' is optimum, special for OSD use only.
	 * takes the variable indicating whether reliability sort in decrease or increas order,
	 * for using generator matrix or parity-check matrix to re-encode, a more general function.
	 *
	 * \param d_min:						minimum hamming distance of code
	 * \param k:							information bit number
	 * \param dv:							codeword to judge, length of n
	 * \param hdr:							hard decoding result of received word
	 * \param r_abs:						absolote value of received word, indicating the reialbility
	 * \param finished_order:				the order OSD that has been finished, set it be current order - 1
	 * \param lambda:						correlation discripency between 'dv' and received word
	 * \param error_num:					number of difference between 'dv' and 'hdr'
	 * \param is_reliability_decreasing:	is_generator == true, this is ture, else this is false
	 * \return whether 'dv' is an optimum codeword
	 */
	static bool is_optimum_with_sort_way(int d_min, int k, const Matrix<GF2>& dv, const Matrix<GF2>& hdr, \
		const Matrix<my_double>& r_abs, int finished_order, my_double lambda, int error_num, bool is_reliability_decreasing) {	
		// this is designed for OSD

		//return false;			// for worse case complexity testing

		if (is_reliability_decreasing == false) {
			int n = dv.size();
			my_double G_vw1 = 0;
			// consider j=1 only, this judge if is optimum
			int q1 = d_min - error_num - (finished_order + 1);		// d_min
			int cnt = 0;
			for (int i = 0; i < n && cnt < q1; ++i) {
				// this can be optimized, test for correctness first
				if (dv(i) == hdr(i)) {
					G_vw1 += r_abs(i);		// it seems the optimum condition has no relation to the reliability ordering shortly after k
					cnt++;
				}
			}
			for (int i = n - k, imax = n - k + finished_order; i <= imax; ++i) {
				G_vw1 += r_abs(i);		// this has relation to the reliability ordering shortly before k
			}

			//cout << "G_vw1=" << G_vw1 << endl;
			return lambda <= G_vw1;
		}
		else {
			int n = dv.size();
			my_double G_vw1 = 0;
			// consider j=1 only, this judge if is optimum
			int q1 = d_min - error_num - (finished_order + 1);		// d_min
			int cnt = 0;
			for (int i = n - 1; i >= 0 && cnt < q1; --i) {
				if (dv(i) == hdr(i)) {
					G_vw1 += r_abs(i);		// it seems the optimum condition has no relation to the reliability ordering shortly after k
					cnt++;
				}
			}
			for (int i = 1; i <= finished_order + 1; ++i) {
				G_vw1 += r_abs(k - i);		// this has relation to the reliability ordering shortly before k
			}

			//cout << "G_vw1=" << G_vw1 << endl;
			return lambda <= G_vw1;
		}
		
	}

	/**
	 * .judge whether current codeword 'dv' is optimum, special for OSD use only.
	 * takes the variable indicating whether reliability sort in decrease or increas order,
	 * for using generator matrix or parity-check matrix to re-encode, a more general function.
	 * This is designed for iterative basis updated (IBU) algorithm used in OSD, 
	 * where the generator (parity-check) matrix's first k (n-k) columns are not sorted in reliability decreasing (increasing) order.
	 * 
	 * \param d_min:						minimum hamming distance of code
	 * \param k:							information bit number
	 * \param dv:							codeword to judge, length of n
	 * \param hdr:							hard decoding result of received word
	 * \param r_abs:						absolote value of received word, indicating the reialbility
	 * \param finished_order:				the order OSD that has been finished, set it be current order - 1
	 * \param lambda:						correlation discripency between 'dv' and received word
	 * \param error_num:					number of difference between 'dv' and 'hdr'
	 * \param is_reliability_decreasing:	is_generator == true, this is ture, else this is false
	 * \param post_permute:					if r_abs permute after that, the info and redundancy will be sorted by reliability 
	 * \return whether 'dv' is an optimum codeword
	 */
	static bool is_optimum_with_sort_way_pp(int d_min, int k, const Matrix<GF2>& dv, const Matrix<GF2>& hdr, 
		const Matrix<my_double>& r_abs, int finished_order, my_double lambda, int error_num, \
		bool is_reliability_decreasing, const Matrix<int>& post_permute) {

		// this is designed for OSD

		//return false;			// for worse case complexity testing

		if (is_reliability_decreasing == false) {
			int n = dv.size();
			my_double G_vw1 = 0;
			// consider j=1 only, this judge if is optimum
			int q1 = d_min - error_num - (finished_order + 1);		// d_min
			int cnt = 0;
			for (int i = 0; i < n && cnt < q1; ++i) {
				int i_post = post_permute(i);

				// this can be optimized, test for correctness first
				if (dv(i_post) == hdr(i_post)) {
					G_vw1 += r_abs(i_post);		// it seems the optimum condition has no relation to the reliability ordering shortly after k
					cnt++;
				}
			}
			for (int i = n - k, imax = n - k + finished_order; i <= imax; ++i) {

				int i_post = post_permute(i);
				G_vw1 += r_abs(i_post);		// this has relation to the reliability ordering shortly before k
			}

			//cout << "G_vw1=" << G_vw1 << endl;
			return lambda <= G_vw1;
		}
		else {
			int n = dv.size();
			my_double G_vw1 = 0;
			// consider j=1 only, this judge if is optimum
			int q1 = d_min - error_num - (finished_order + 1);		// d_min
			int cnt = 0;
			for (int i = n - 1; i >= 0 && cnt < q1; --i) {
				int i_post = post_permute(i);

				if (dv(i_post) == hdr(i_post)) {
					G_vw1 += r_abs(i_post);		// it seems the optimum condition has no relation to the reliability ordering shortly after k
					cnt++;
				}
			}
			for (int i = k - 1, imin = k - finished_order - 1; i >= imin; --i) {
				int i_post = post_permute(i);
				G_vw1 += r_abs(i_post);		// this has relation to the reliability ordering shortly before k
			}

			//cout << "G_vw1=" << G_vw1 << endl;
			return lambda <= G_vw1;
		}

	}
};

class data_ope {
public:
	/**
	 * .packed GF2 data into unsigned data, to compress or accelerate computing
	 *
	 * \param unpacked: a GF2 Matrix to be packed
	 * \param packed: will set to packed result, should have enough size, >= (unpacked.size()-1)/32+1
	 */
	static void pack_data(const Matrix<GF2>& unpacked, Matrix<unsigned>& packed) {
		// Matrix packed should have enough size, >= (unpacked.size()-1)/32+1
		int len = unpacked.size();
		int i = 0;
		int packed_ind = 0;
		while (i < len - 31) {
			int packed_num = 0;
			for (int j = 0; j < 32; ++j) {
				packed_num <<= 1;
				packed_num |= (int)unpacked(i + j);		// 0 or 1
			}
			packed(packed_ind) = packed_num;
			packed_ind++;
			i += 32;
		}

		// the last packed data
		if (i < len) {
			int packed_num = 0;
			for (int j = i; j < len; ++j) {
				packed_num <<= 1;
				packed_num |= (int)unpacked(j);		// 0 or 1
			}
			packed(packed_ind) = packed_num;
		}
	}

	/**
	 * .unpack unsigned data into GF2 data
	 *
	 * \param packed: data to be unpacked
	 * \param unpacked: will be set to GF2 data, as result, should have size as desired
	 */
	static void unpack_data(const Matrix<unsigned>& packed, Matrix<GF2>& unpacked) {
		// Matrix unpacked should have size as desired
		int packed_len = unpacked.size();
		int packed_ind = 0;
		int unpacked_block_ind = 31;
		while (unpacked_block_ind < packed_len) {
			int num = packed(packed_ind);
			for (int j = 0; j < 32; ++j) {
				unpacked(unpacked_block_ind - j) = num & 1;
				num >>= 1;
			}
			packed_ind++;
			unpacked_block_ind += 32;
		}

		// deal with the last block
		int left_len = 32 - (unpacked_block_ind - packed_len + 1);
		if (left_len != 0) {
			int num = packed(packed_ind);
			for (int j = 1; j <= left_len; ++j) {
				unpacked(packed_len - j) = num & 1;
				num >>= 1;
			}
		}
	}

	/**
	 * .designed for matrix data packed by rows, like generator matrix of a code
	 *
	 * \param unpacked: matrix to be packed by rows
	 * \param packed: should have size unpacked.row(), and its elements should have
	 * enough size, >= (unpacked.col()-1)/32+1, will be set as result, with element of row packed data
	 */
	static void pack_row_data(const Matrix<GF2>& unpacked, vector<Matrix<unsigned>>& packed) {
		int r = unpacked.row();
		for (int i = 0; i < r; ++i) {
			pack_data(unpacked.get_row(i), packed[i]);
		}
	}

	/**
	 * .turn a 'num' to Matrix of GF2, the size of the matrix is fix to 'total_bit', natual order
	 */
	static Matrix<GF2> i2GF2_natual(unsigned long long num, int total_bit) {
		Matrix<GF2> ans(1, total_bit, '0');
		int i = total_bit - 1;
		while (num != 0) {
			ans(i) = num & 1;
			num >>= 1;
			i--;
		}
		return ans;
	}

	/**
	 * .turn a 'num' to Matrix of GF2, the size of the matrix is fix to 'total_bit', OSD flip order
	 */
	static Matrix<GF2> i2GF2_osd(unsigned long long index, int total_bit) {
		Matrix<GF2> ans(1, total_bit, '0');
		if (index == 0) {
			return ans;
		}

		// decide which order currently, the wirting is confusing, should be simplified, if have time
		if (total_bit_computed != total_bit) {
			total_bit_combinatorial_accumulation.resize(total_bit + 1, 1, true);
			total_bit_combinatorial_accumulation(0) = 0;
			for (int i = 1; i <= total_bit; ++i) {
				total_bit_combinatorial_accumulation(i) = total_bit_combinatorial_accumulation(i - 1) + my::n_choose_k_long(total_bit, i);
			}
			//cout << "total_bit_combinatorial_accumulation" << total_bit_combinatorial_accumulation;
			total_bit_computed = total_bit;
		}

		// use the binary search method
		int start = 0;
		int end = total_bit;

		// binary search over v, O(log(end)) comparison and program jump
		while (start < end) {		// the standard binary search template
			int mid = (start + end) >> 1;
			if (total_bit_combinatorial_accumulation(mid) < index) {
				start = mid + 1;
			}
			else {
				end = mid;
			}
		}
		int outer_order = total_bit_combinatorial_accumulation(start) >= index ? start : start - 1;
		//cout << "outer_order = " << outer_order << endl;

		// we return here as all bit at last position
		//cout << "outer_order = " << outer_order << endl;

		// decide each bit position
		// recursive computation, find the last bit position
		unsigned long long inner_index = index - total_bit_combinatorial_accumulation(outer_order - 1);
		//cout << "inner_index = " << inner_index << endl;

		//if (inner_index > 0.5 * index) {
		//	// we can use one more outer_order by experiance
		//	outer_order++;		// causing performance degradation
		//}
		//cout << "outer_order new = " << outer_order << endl;

		ans.set_part(0, total_bit - outer_order, -1, -1, 1);
		return ans;

		// find the exact flip pattern in OSD order, abandoned

		int inner_total_bit = total_bit;
		int inner_order = outer_order;
		while (inner_order != 0) {
			/*cout << "-------" << endl;
			cout << "inner_index = " << inner_index << endl;
			cout << "inner_total_bit = " << inner_total_bit << endl;
			cout << "inner_order = " << inner_order << endl;*/

			ans(last_bit_ind(inner_index, inner_total_bit, inner_order)) = 1;
			//cout << "ans" << ans;
		}
		return ans;
	}

	static int last_bit_ind(unsigned long long& inner_index, int& inner_total_bit, int& inner_order) {
		int max_group_num = inner_total_bit - inner_order + 1;
		unsigned long long last_accumulation = 0;
		unsigned long long accumulation = last_accumulation;
		int i = 1;
		while (i <= max_group_num) {
			last_accumulation = accumulation;
			accumulation += my::n_choose_k_long(inner_total_bit - i, inner_order - 1);
			//cout << "accumulation = " << accumulation << endl;
			if (accumulation >= inner_index) {
				// the last bit at index: total_bit - i
				break;
			}
			++i;
		}
		inner_index -= last_accumulation;
		inner_total_bit = inner_total_bit - i;
		inner_order--;
		return inner_total_bit;			// return the index of last bit
	}

	static Matrix<unsigned long long> total_bit_combinatorial_accumulation;
	static int total_bit_computed;
};
Matrix<unsigned long long> data_ope::total_bit_combinatorial_accumulation;
int data_ope::total_bit_computed = -1;

/**
 * .the sequence of TEP can be adjusted, with better early stoping property
 */
class flip_TEP {
public:

	/**
	* .find the next bit flip pattern in OSD
	*
	* \param is_flipped_mat: current flip pattern
	* \param k: information bit numbers
	* \param order: order in OSD, aka. number of flip pattern
	* \return whether travelled all flipp pattern
	*/
	static bool next(Matrix<int>& is_flipped_mat, int k, int order) {		// this should be replaced by Matrix_common::OSD_next_TEP()
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

	static void next_v2(Matrix<int>& TEP_now, int n) {
		Matrix_common::OSD_next_TEP(TEP_now, n);
	}

	/**
	 * .compute difference between 'last_flipped_mat' and 'flipped_mat'
	 *
	 * \param last_flipped_mat: the last flipped pattern, same size as the last fipped pattern or empty, sorted
	 * \param flipped_mat: current flip pattern, sorted
	 * \param fmd: this will be set as difference between 'last_flipped_mat' and 'flipped_mat'
	 * i.e., 'last_flipped_mat'={2,3}, and 'flipped_mat'={2,4}, then 'fmd'={3,4}
	 */
	static void diff(const Matrix<int>& last_flipped_mat, const Matrix<int>& flipped_mat, Matrix<int>& fmd) {
		if (last_flipped_mat.size() != 0) {
			int total_size = last_flipped_mat.size();
			int ind_last = 0;
			int ind_flipped = 0;
			fmd.resize(1, 0, false);
			while (ind_last < total_size && ind_flipped < total_size) {
				if (last_flipped_mat(ind_last) == flipped_mat(ind_flipped)) {
					ind_last++;
					ind_flipped++;
				}
				else if (last_flipped_mat(ind_last) < flipped_mat(ind_flipped)) {
					fmd.push_back(last_flipped_mat(ind_last));
					ind_last++;
				}
				else {
					fmd.push_back(flipped_mat(ind_flipped));
					ind_flipped++;
				}
			}

			if (ind_flipped == total_size) {
				while (ind_last < total_size) {
					fmd.push_back(last_flipped_mat(ind_last));
					ind_last++;
				}
			}
			else {
				while (ind_flipped < total_size) {
					fmd.push_back(flipped_mat(ind_flipped));
					ind_flipped++;
				}
			}
		}
		else {
			fmd = flipped_mat;
		}
	}

	/**
	 * .counting the cdf according to osd order 
	 */
	static Matrix<GF2> get_estimated_flip_pattern(int len, my_double Pe_ML) {
		// counting the flip index
		//Pe_ML /= 10;	// if the probability of best codeword occuring is less than one tenth of Pe_ML, we stop finding next codeword
		my_double Pr_ZRi = pow(2, -len);
		//cout << "len = " << len << endl;
		//cout << "Pr_ZRi = " << Pr_ZRi << endl;

		unsigned long long flip_index = unsigned long long(Pe_ML / Pr_ZRi);		// do not need to be accurate
		//cout << "flip_index = " << flip_index << endl;

		//return data_ope::i2GF2_natual(flip_index, len);
		return data_ope::i2GF2_osd(flip_index, len);
		
			// actually it is not sorted flip list, for simplicity we use direct <int to GF2> transformation
			// there are some uncertainty for performance degradation, to be tested
	}
};

class jiabao_TEP{
public:
	static int n_choose_k(int _n, int _k) {
		if (_k <= _n - _k) {
			// cancell (n-k)!
			int result = 1;
			for (int i = 0; i < _k; ++i) {
				result = result * (_n - i) / (i + 1);
			}
			return result;
		}
		else {
			// cancell k!, for (n,k) = (n,n-k)			
			return n_choose_k(_n, _n - _k);
		}
	}
	static int** generate(int& list, int K = 10, int order = 3) {
		// parameter input
		list = 0;
		for (int i = 1; i <= order; ++i) {
			list += my::n_choose_k(K, i);
		}
		int** TEP = new int* [list];
		for (int i = 0; i < list; ++i) {
			TEP[i] = new int[order];
			memset(TEP[i], 0, order * sizeof(int));
		}

		int c = 0;
		int temp = 1;
		int* e = new int[K];
		while (temp <= order) {
			for (int i = 0; i < K; i++) {
				e[i] = i >= K - temp;
			}
			do {
				int ii = 0;
				for (int i = 0; i < K; i++) {
					if (e[i] == 1) {
						TEP[c][ii] = i + 1;
						ii++;
					}
				}
				c++;
			} while (next_permutation(e, e + K));
			temp++;
		}
		delete[]e;
		return TEP;
	}
	static int** combine(int** TEP_first, int list_first, int order_first, int** TEP_second, int list_second, int order_second) {
		int first_K = TEP_first[0][0];

		int total_list = list_first * list_second + list_first + list_second;
		int total_order = order_first + order_second;
		int** TEP = new int* [total_list];
		for (int i = 0; i < total_list; ++i) {
			TEP[i] = new int[total_order];
			memset(TEP[i], 0, total_order * sizeof(int));
		}

		for (int j = 0; j < list_second; ++j) {
			for (int k = 0; k < order_second && TEP_second[j][k] != 0; ++k) {
				TEP[j][k] = TEP_second[j][k] + first_K;
			}
		}
		for (int i = 0; i < list_first; ++i) {
			for (int k = 0; k < order_first && TEP_first[i][k] != 0; ++k) {
				TEP[list_second + i][k] = TEP_first[i][k];
			}
		}

		int list_shift = list_first + list_second;

		for (int i = 0; i < list_first; ++i) {
			for (int j = 0; j < list_second; ++j) {
				int kk = 0;
				for (int k = 0; k < order_first && TEP_first[i][k] != 0; ++k) {
					TEP[list_shift + i * list_second + j][kk] = TEP_first[i][k];
					kk++;
				}
				for (int k = 0; k < order_second && TEP_second[j][k] != 0; ++k) {
					TEP[list_shift + i * list_second + j][kk] = TEP_second[j][k] + first_K;
					kk++;
				}
			}
		}
		return TEP;
	}
	static void delete_matrix(int** pointer, int row_num) {
		for (int i = 0; i < row_num; ++i) {
			delete[] pointer[i];
		}
		delete[] pointer;
	}
	static int** generate_segmentation(int _k, int order_first, int _K, int order_second, int& _list_total) {
		_list_total = 0;
		for (int i = 0; i <= order_first; ++i) {
			_list_total += n_choose_k(_k, i);
		}

		int _list_total_mul = 0;
		int k_diff = _K - _k;
		for (int i = 0; i <= order_second; ++i) {
			_list_total_mul += n_choose_k(k_diff, i);
		}

		_list_total *= _list_total_mul;
		_list_total--;
		//cout << "_list_total = " << _list_total << endl;

		int order_total = order_first + order_second;
		int** TEP = new int* [_list_total];
		for (int i = 0; i < _list_total; ++i) {
			TEP[i] = new int[order_total];
			memset(TEP[i], 0, order_total * sizeof(int));
		}

		int c = 0;
		int temp = 1;
		int* e = new int[_K];
		while (temp <= order_total) {
			for (int i = 0; i < _K; i++) {
				e[i] = i >= _K - temp;
			}
			do {
				int ii = 0;
				for (int i = 0; i < _K; i++) {
					if (e[i] == 1) {
						TEP[c][ii] = i + 1;
						ii++;
					}
				}
				/*for (int j = 0; j < temp; j++) {
					cout << TEP[c][j] << '\t';
				}
				cout << endl;*/

				// select the TEP
				int first_bin_order = 0;
				int second_bin_order = 0;
				for (int i = 0; i < temp; ++i) {
					if (TEP[c][i] > _k) {
						second_bin_order++;
					}
					else {
						first_bin_order++;
					}
				}

				if (first_bin_order <= order_first && second_bin_order <= order_second){
					c++;
				}
			} while (next_permutation(e, e + _K) && c != _list_total);
			temp++;
		}
		delete[]e;
		return TEP;
	}
};

struct info_bin_continous_record {
	int start_ind;
	int end_ind;
	info_bin_continous_record(int _start_ind = 0, int _end_ind = 0) :start_ind(_start_ind), end_ind(_end_ind) {}
	bool operator < (const info_bin_continous_record& ibcr) const {
		return true;				// we use true to indicate it has no ordered structure

		//return start_ind < ibcr.start_ind;
	}
	friend ostream& operator << (ostream& out, const info_bin_continous_record& ibcr) {
		out << "(" << ibcr.start_ind << ", " << ibcr.end_ind << ")";
		return out;
	}
};

// a function suitable for all info set distribution scheme
class Basis_update {
public:

	int n;
	int k;
	int red;											// these 3 int is duplicated for the class only

	// the variable during get_MRIP_sys_G
	Matrix<int> col_ind_to_reliability;					// the smaller the reliability number, indicating that its order is more forward
	Matrix<int> permute_target_col_ind_to_reliability;
	Matrix<int> permute_target_reliability_to_col_ind;
	Matrix<int> redundancy_set_reliability;
	sList<int> info_set_reliability_sList;

	Basis_update(int _n, int _k) :n(_n), k(_k), red(_n - _k) {

		// variable to store mapping from index to col_ind and its reverse mapping
		col_ind_to_reliability.resize(1, n, false);
		permute_target_col_ind_to_reliability.resize(1, n, false);
		permute_target_reliability_to_col_ind.resize(1, n, false);
		redundancy_set_reliability.resize(1, red, false);
		info_set_reliability_sList.resize(1, k, false);
	}

	int solve(const Matrix<int>& reliability_to_col_ind, Matrix<GF2>& G_target, Matrix<int>& permute_target) {

		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		for (int i = 0; i < n; ++i) {
			col_ind_to_reliability(reliability_to_col_ind(i)) = i;
		}
		//cout << "col_ind_to_reliability" << col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_col_ind_to_reliability(i) = col_ind_to_reliability(permute_target(i));
		}
		//cout << "permute_target_col_ind_to_reliability" << permute_target_col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(i)) = i;
		}
		//cout << "permute_target_reliability_to_col_ind" << permute_target_reliability_to_col_ind;

		// how to mantain the max reliability number of info set? may be using priority queue, or heap of myself
		info_set_reliability_sList.clear();
		for (int i = 0; i < k; ++i) {
			info_set_reliability_sList.push_back(permute_target_col_ind_to_reliability(i));
		}
		info_set_reliability_sList.sort('>');
		info_set_reliability_sList.build();

		for (int i = k; i < n; ++i) {
			redundancy_set_reliability(i - k) = permute_target_col_ind_to_reliability(i);
		}
		redundancy_set_reliability.sort('<');
		//cout << "redundancy_set_reliability" << redundancy_set_reliability;		// has problem here, since we donot delet the element

		// note that this is not the MRB yet

		// what if some MRP are out? we can continue permuting the redundancy position over the info position

		// check from the MRP of redundancy set, and try to insert it over the info set

		// some new idea should be presented here.	

		// we only consider redundancy_set with reliability <= k, since info_set with reliability > k has finished is decision
		int q = 0;
		for (;q < red && info_set_reliability_sList.size() != 0 && redundancy_set_reliability(q) < info_set_reliability_sList.start_val();
			++q) {

			// count the 1 position of column w.r.t. iter_red, find the max reliability
			int accept_col_ind_can = permute_target_reliability_to_col_ind(redundancy_set_reliability(q));
			bool if_update_basis = false;

			// test if can switch the iter_red to the info_set
			//cout << "info_set_reliability_sList" << info_set_reliability_sList;
			//cout << "info_set_reliability_sList.start_val() = " << info_set_reliability_sList.start_val() << endl;

			int throw_col_ind = permute_target_reliability_to_col_ind(info_set_reliability_sList.start_val());
			if (G_target(throw_col_ind, accept_col_ind_can) == 1) {
				// we have half of the chance getting inside this if clause

				// the first index should be removed
				info_set_reliability_sList.remove_start();	// this should be fix and you should erase the if clause
				if_update_basis = true;						// it is sure to update basis since we are in the for loop
			}
			else {
				// and half of the chance getting inside this else clause

				info_set_reliability_sList.move_start();								// move current index to start
				int sList_next_ind_tmp = info_set_reliability_sList.get_next_ind();		// to get next index to test if it is valid

				while (sList_next_ind_tmp != -1) {		// the list is not traversed to the end

					int switch_reliability_tmp = info_set_reliability_sList(sList_next_ind_tmp);
						// get information set element reliability order

					if (redundancy_set_reliability(q) < switch_reliability_tmp) {	
						// the redundancy set element has small reliability order, i.e., more reliable

						throw_col_ind = permute_target_reliability_to_col_ind(switch_reliability_tmp);
						if (G_target(throw_col_ind, accept_col_ind_can) == 1) {
							// this entry equals to 1 make sures redundancy set element 
							// can be switched into the information set, forming a new information set

							info_set_reliability_sList.remove_next();
								// we have to set the current index former to the remove node, hence wirte 2 case, wheter removing head
							if_update_basis = true;
							break;
						}
						else {
							// this entry equals to 0, not forming a new information set, 
							// continue to the next less reliable element in information set
							// the number of times of this while loop will decrease exponentially, half by half. 
							// so the expectation of times of this while loop is only 2.

							info_set_reliability_sList.move_next();								// find the next index
							sList_next_ind_tmp = info_set_reliability_sList.get_next_ind();		// update next index to storage
						}
					}
					else {
						// the redundancy set element is not reliable compared to all information set to enter the information set

						break;
					}
				}
			}

			if (if_update_basis) {

				// we can switch the redundant set, to get a more reliable base

				/*. the following operations are perfromed on G_target, 
					same as Algorithm 2: A Basis Update on $\check{\mathbf G}$ in the paper

					parameter transfer:

					lambda = throw_col_ind
					mu = accept_col_ind_can
				*/

				/*
				// switch columns of throw_col_ind and accept_col_ind_can
				G_target.switch_col(throw_col_ind, accept_col_ind_can);			// O(k)


				// do row transformation on G_target to keep identity on left
				for (int i = 0; i < k; ++i) {
					if (G_target(i, throw_col_ind) == 1 && i != throw_col_ind) {// if the left-part is not diagonal, this will be a problem

						// let the hardware implementation people to consider how to optimized, 
						// we permute parity part to ensure that matrix has identity part on the left hand side.
						// may be this problem is in the range of column permutation and should be solved at the beginning
						G_target(i, throw_col_ind) = 0;
						for (int j = k; j < n; ++j) {							// O(k*(n-k))
							G_target(i, j) += G_target(throw_col_ind, j);		// this can be done in parallel and save computation time
						}
					}
				}
				*/

				// a new and simpler equivalent operations on G_target, doing the same things as above, using the structure of G_target
				// sturcture: G_target = [I^(k)    P]. no column switch are needed, row operation is all right for column throw_col_ind and 
				// column throw_col_ind are unchanged after all the operations.

				int row_ind_tmp_for_G_target = 0;
				int row_ind_for_G_target_throw_col_ind = throw_col_ind * n;		// just record for acceleration, optimized for vital codes

				G_target[row_ind_for_G_target_throw_col_ind + accept_col_ind_can] = 0;
				for (int i = 0; i < k; ++i) {
					if (G_target[row_ind_tmp_for_G_target +  accept_col_ind_can] == 1) {
						for (int j = k; j < n; ++j) {
							G_target[row_ind_tmp_for_G_target + j] += G_target[row_ind_for_G_target_throw_col_ind + j];
								// this will be more meaningful in hardware implementation
						}
					}
					row_ind_tmp_for_G_target += n;
				}
				G_target[row_ind_for_G_target_throw_col_ind + accept_col_ind_can] = 1;		// it's verified that this perform rightly


				//cout << "G_target (after)" << G_target;
				//cout << "G_target (2)" << G_target;

				// mantian the info_set_reliability_heap

				// update permute_target, permute_target_col_ind_to_reliability, permute_target_reliability_to_col_ind
				permute_target.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_col_ind_to_reliability.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(throw_col_ind)) = throw_col_ind;
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(accept_col_ind_can)) = accept_col_ind_can;

				//cout << "permute_target" << permute_target;
				//cout << "G_target" << G_target;
				//cout << "accept_col_ind_can = " << accept_col_ind_can << endl;
				//cout << "throw_col_ind = " << throw_col_ind << endl;
			}
		}
		//cout << "info_set_reliability_sList" << info_set_reliability_sList;

		//cout << "G_target (last)" << G_target;

		return q;
	}
};

/**
 * .record for target index and shift value, for info set that has the most common elements with MRPs
 */
struct target_ind_and_shift {
	int target_ind;
	int shift;
	target_ind_and_shift(int _target_ind = 0, int _shift = 0) :target_ind(_target_ind), shift(_shift) {}
	bool operator < (const target_ind_and_shift& tias) const {
		return true;				// we use true to indicate it has no ordered structure

		//return start_ind < ibcr.start_ind;
	}
	friend ostream& operator << (ostream& out, const target_ind_and_shift& tias) {
		out << "(" << tias.target_ind << ", " << tias.shift << ")";
		return out;
	}
};
