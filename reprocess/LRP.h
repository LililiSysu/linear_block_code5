#pragma once
/*****************************************************************//**
 * \file   LRP_Reprocessing.h
 * \brief  Least Reliable Position Reprocessing algorithm, for at most 3dB gain versus hard decoding
 * \param	'u' --[channel_coding]--> 'v' --[BPSK_modulation]--> 'c' --[AWGN_channel]--> 'r'('recv')
 *		   --[BPSK_demodulation]--> 'hdr'('z') --[decode_v]--> 'dv'('v') --[v2u]--> 'du'
 * 
 * \author lilili
 * \date   October 2022
 *********************************************************************/

#include"reprocess_common.h"

class GMD {			// only valid for GF2 field coding and decoding scheme
public:
	static int num_early_stop;
	static unsigned long long total_used_list_num;		// counting used list num
	static bool is_ML;
	static my_double lambda_min;

	/**
	 * .apply decode method of 'code_instance' to 'r', by erasing bit to improve
	 *
	 * \param code_instance: the coding and decoding scheme
	 * \param r: received codeword
	 * \return decoded result, u
	 */
	template<class code_type>
	static Matrix<GF2> decode(code_type& code_instance, const Matrix<my_double>& r)
	{
		return code_instance.v2u(decode_v(code_instance, r));
	}

	/**
	 * .apply decode method of 'code_instance' to 'r', by erasing bit to improve
	 *
	 * \param code_instance: the coding and decoding scheme
	 * \param r: received codeword
	 * \return decoded result, u
	 */
	template<class code_type>
	static Matrix<GF2> decode_v(code_type& code_instance, const Matrix<my_double>& r)
	{
		// a shell for code_type, such as BCH
		//return code_instance.decode(r);

		is_ML = false;
		int n = r.size();
		Matrix<my_double> r_abs(1, n);
		for (int i = 0; i < n; ++i) {
			r_abs(i) = my::abs(r(i));
		}
		Matrix<int> r_ind = r_abs.sort_with_ind();		// note that 'r_abs' is sorted
		Matrix<GF2> hdr = BPSK::demodulation(r);

		int error_num;
		int d_min = code_instance.get_d();
		if (d_min % 2 == 1) {
			// d_min is odd, then erase 0,2,...,d_min-1 least reliable positions
			// number of them is (d_min+1)/2

			// first not erase any position
			lambda_min = r_abs(n - 1) * n;	// the upper bound for lambda_min
			Matrix<GF2> best_dv = code_instance.decode_v(hdr);	// for the first set is we do demand in v space
			if (code_instance.is_in_v_space(best_dv)) {
				lambda_min = Measure::correlation_discrepancy_v(r, best_dv);
				error_num = Measure::error_num;
				if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda_min, error_num)) {
					//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
					is_ML = true;
					num_early_stop++;
					total_used_list_num++;
					return best_dv;
				}
			}

			Matrix<bool> is_erased(1, n, '0');
			for (int erased_num = 0; erased_num < d_min - 1; erased_num += 2) {
				is_erased(r_ind(erased_num)) = 1;
				is_erased(r_ind(erased_num + 1)) = 1;
				Matrix<GF2> dv = decode_erasure_v(code_instance, hdr, is_erased);
				if (!code_instance.is_in_v_space(dv))	continue; // for later dv, we demand it in v space
				my_double lambda = Measure::correlation_discrepancy_v(r, dv);
				error_num = Measure::error_num;
				if (lambda < lambda_min) {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;/*
					cout << "dv=" << dv;
					cout << "lambda_min=" << lambda << endl;*/
					if (optimum_condition::is_optimum(d_min, dv, hdr, r_abs, r_ind, lambda, error_num)) {
						//cout << "optimum condition satisfied, quit testing at iter " << erased_num / 2 + 2 << endl;
						is_ML = true;
						num_early_stop++;
						total_used_list_num += (unsigned long long) erased_num + 3;
						return best_dv;
					}
				}
			}

			//cout << "best_dv=" << best_dv;
			total_used_list_num += d_min;
			return best_dv;
		}
		else {
			// d_min is even, then erase 1,3,...,d_min-1 least reliable positions
			// number of them is (d_min+1)/2

			// first erase only one position
			Matrix<bool> is_erased(1, n, '0');
			is_erased(r_ind(0)) = 1;

			lambda_min = r_abs(n - 1) * n;	// the upper bound for lambda_min
			Matrix<GF2> best_dv = decode_erasure_v(code_instance, hdr, is_erased);	// for the first set is we do demand in v space
			if (code_instance.is_in_v_space(best_dv)) {
				lambda_min = Measure::correlation_discrepancy_v(r, best_dv);
				error_num = Measure::error_num;
				if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda_min, error_num)) {
					//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
					is_ML = true;
					num_early_stop++;
					total_used_list_num++;
					return best_dv;
				}
			}

			for (int erased_num = 1; erased_num < d_min - 1; erased_num += 2) {
				is_erased(r_ind(erased_num)) = 1;
				is_erased(r_ind(erased_num + 1)) = 1;
				Matrix<GF2> dv = decode_erasure_v(code_instance, hdr, is_erased);
				if (!code_instance.is_in_v_space(dv))	continue; // for later dv, we demand it in v space
				my_double lambda = Measure::correlation_discrepancy_v(r, dv);
				error_num = Measure::error_num;
				if (lambda < lambda_min) {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda, error_num)) {
						//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
						is_ML = true;
						num_early_stop++;
						total_used_list_num += (unsigned long long) erased_num + 3;
						return best_dv;
					}
				}
			}
			total_used_list_num += d_min;
			return best_dv;
		}
	}

	/**
	 * .apply decode method of 'code_instance' to 'r' with erasure position, return v
	 *  if 2v+e+1<=d_{min}, decode success, where v is number of errors
	 *  out of erasure and e is number of erasure bits
	 *
	 * \param code_instance: the coding and decoding scheme
	 * \param r: received codeword, after BPSK_demodulation
	 * \param is_erased: a matrix of size (1,n), with 1 means an erase position and 0 means versa
	 * \return decoded result, u
	 */
	template<class code_type>
	static Matrix<GF2> decode_erasure_v(code_type& code_instance, const Matrix<GF2>& hdr, \
		const Matrix<bool>& is_erased /* 1 means an erase position and 0 means versa */)
	{
		/* choose the decode codeword that has minimum errors out of erase position */

		/* set all erase positions to be 0 or 1 */
		Matrix<GF2> r0 = hdr;
		Matrix<GF2> r1 = hdr;
		int rlen = hdr.size();
		for (int i = 0; i < rlen; ++i) {
			if (is_erased(i) == false);
			else {
				r0(i)= 0;
				r1(i)= 1;
			}
		}
		Matrix<GF2> v0 = code_instance.decode_v(r0);
		Matrix<GF2> v1 = code_instance.decode_v(r1);
		if (!code_instance.is_in_v_space(v0)) {
			return v1;
		}
		if (!code_instance.is_in_v_space(v1)) {
			return v0;
		}

		/* decide which one to stay */
		/* choose the decode codeword that has minimum errors out of erase position */
		/* NOTE: different with GMD statement that keep the codeword with least correlation discrepancy*/
		int error_cnt0 = 0;
		int error_cnt1 = 0;
		for (int i = 0; i < rlen; ++i) {
			if (is_erased(i) == false) {
				error_cnt0 += (v0(i) != r0(i));
				error_cnt1 += (v1(i) != r1(i));
			}
		}

		if (error_cnt0 <= error_cnt1) {	// keep the erase result with 0 codeword
			return v0;
		}
		else {							// keep the erase result with 1 codeword
			return v1;
		}
	}

	/**
	 * .apply decode method of 'code_instance' to 'r' with erasure position, return u
	 *  if 2v+e+1<=d_{min}, decode success, where v is number of errors
	 *  out of erasure and e is number of erasure bits
	 *
	 * \param code_instance: the coding and decoding scheme
	 * \param r: received codeword, after BPSK_demodulation
	 * \param is_erased: a matrix of size (1,n), with 1 means an erase position and 0 means versa
	 * \return decoded result, u
	 */
	template<class code_type>
	static Matrix<GF2> decode_erasure(code_type& code_instance, const Matrix<GF2>& hdr, \
		const Matrix<bool>& is_erased /* 1 means an erase position and 0 means versa */)
	{
		/* choose the decode codeword that has minimum errors out of erase position */

		return code_instance.v2u(decode_erasure_v(code_instance, hdr, is_erased));

	}
};

bool GMD::is_ML = false;
int GMD::num_early_stop = 0;
unsigned long long GMD::total_used_list_num = 0;
my_double GMD::lambda_min = 0;

// write the method as static class function, to be fixed
class Chase{
public:

	static int num_early_stop;
	static unsigned long long total_used_list_num;		// counting used list num
	static bool is_ML;
	static my_double lambda_min;

	/**
	 * .apply decode method of 'code_instance' to 'r', by flipping bit to improve
	 *
	 * \param code_instance: the coding and decoding scheme
	 * \param r: received codeword, after BPSK_demodulation
	 * \return decoded result, u
	 */
	template<class code_type>
	static Matrix<GF2> decode3(code_type& code_instance, const Matrix<my_double>& r)
	{
		// a shell for code_type, such as BCH
		//return code_instance.decode(r);

		is_ML = false;
		int n = r.size();
		Matrix<my_double> r_abs(1, n);
		for (int i = 0; i < n; ++i) {
			r_abs(i)= r(i) >= 0 ? r(i) : -r(i);
		}
		Matrix<int> r_ind = r_abs.sort_with_ind();		// note that 'r_abs' is sorted
		Matrix<GF2> hdr = BPSK::demodulation(r);

		int error_num;
		int d_min = code_instance.get_d();
		if (d_min % 2 == 1) {
			// d_min is odd, then flip 0,2,...,d_min-1 least reliable positions
			// number of them is (d_min+1)/2

			// first not flip any position
			lambda_min = r_abs(n - 1) * n;	// the upper bound for lambda_min
			Matrix<GF2> best_dv = code_instance.decode_v(hdr);	// for the first set is we do demand in v space
			if (code_instance.is_in_v_space(best_dv)) {
				lambda_min = Measure::correlation_discrepancy_v(r, best_dv);
				error_num = Measure::error_num;
				if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda_min, error_num)) {
					//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
					is_ML = true;
					return code_instance.v2u(best_dv);
				}
			}

			Matrix<bool> is_flipped(1, n, '0');
			for (int flipped_num = 0; flipped_num < d_min - 1; flipped_num += 2) {
				is_flipped(r_ind(flipped_num))= 1;
				is_flipped(r_ind(flipped_num + 1))= 1;
				Matrix<GF2> dv = decode_flip_v(code_instance, hdr, is_flipped);
				if (!code_instance.is_in_v_space(dv))	continue; // for later dv, we demand it in v space
				my_double lambda = Measure::correlation_discrepancy_v(r, dv);
				error_num = Measure::error_num;
				if (lambda < lambda_min) {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;/*
					cout << "dv=" << dv;
					cout << "lambda_min=" << lambda << endl;*/
					if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda, error_num)) {
						//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
						is_ML = true;
						return code_instance.v2u(best_dv);
					}
				}
			}

			//cout << "best_dv=" << best_dv;
			return code_instance.v2u(best_dv);
		}
		else {
			// d_min is even, then flip 1,3,...,d_min-1 least reliable positions
			// number of them is (d_min+1)/2

			// first not flip only one position
			Matrix<bool> is_flipped(1, n, '0');
			is_flipped(r_ind(0))= 1;

			lambda_min = r_abs(n - 1) * n;	// the upper bound for lambda_min
			Matrix<GF2> best_dv = decode_flip_v(code_instance, hdr, is_flipped);	// for the first set is we do demand in v space
			if (code_instance.is_in_v_space(best_dv)) {
				lambda_min = Measure::correlation_discrepancy_v(r, best_dv);
				error_num = Measure::error_num;
				if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda_min, error_num)) {
					//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
					is_ML = true;
					return code_instance.v2u(best_dv);
				}
			}

			for (int flipped_num = 1; flipped_num < d_min - 1; flipped_num += 2) {
				is_flipped(r_ind(flipped_num))= 1;
				is_flipped(r_ind(flipped_num + 1))= 1;
				Matrix<GF2> dv = decode_flip_v(code_instance, hdr, is_flipped);
				if (!code_instance.is_in_v_space(dv))	continue; // for later dv, we demand it in v space
				my_double lambda = Measure::correlation_discrepancy_v(r, dv);
				error_num = Measure::error_num;
				if (lambda < lambda_min) {
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;
					if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda, error_num)) {
						//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
						is_ML = true;
						return code_instance.v2u(best_dv);
					}
				}
			}
			return code_instance.v2u(best_dv);
		}
	}

	/**
	 * .apply decode method of 'code_instance' to 'r' with erasure position, return v
	 *  if 2v+e+1<=d_{min}, decode success, where v is number of errors
	 *  out of erasure and e is number of erasure bits
	 *
	 * \param code_instance: the coding and decoding scheme
	 * \param r: received codeword, after BPSK_demodulation
	 * \param is_flipped: a matrix of size (1,n), with 1 means an flipped position and 0 means versa
	 * \return decoded result, u
	 */
	template<class code_type>
	static Matrix<GF2> decode_flip_v(code_type& code_instance, const Matrix<GF2>& hdr, \
			const Matrix<bool>& is_flipped /* 1 means an flip position and 0 means versa */)
	{
		/* set all flip positions flipped */
		Matrix<GF2> r1 = hdr;
		int rlen = hdr.size();
		for (int i = 0; i < rlen; ++i) {
			if (is_flipped(i) == false);
			else {
				r1(i) += 1;/* flip this bit */
			}
		}
		//Matrix<GF2> v1 = code_instance.decode_v(r1);
		return code_instance.decode_v(r1);		// v1 must stay
	}

	/**
	 * .apply decode method of 'code_instance' to 'r' with erasure position, return u
	 *  if 2v+e+1<=d_{min}, decode success, where v is number of errors
	 *  out of erasure and e is number of erasure bits
	 *
	 * \param code_instance: the coding and decoding scheme
	 * \param r: received codeword, after BPSK_demodulation
	 * \param is_flipped: a matrix of size (1,n), with 1 means an flipped position and 0 means versa
	 * \return decoded result, u
	 */
	template<class code_type>
	static Matrix<GF2> decode_flip(code_type& code_instance, const Matrix<GF2>& hdr, \
			const Matrix<bool>& is_flipped /* 1 means an flip position and 0 means versa */)
	{
		return code_instance.v2u(decode_flip_v(code_instance, hdr, is_flipped));
	}

	// return decoded result, u
	template<class code_type>
	static Matrix<GF2> decode2(code_type& code_instance, const Matrix<my_double>& r, int eta = -1) {
		return code_instance.v2u(decode2_v(code_instance, r, eta));
	}

	/**
	 * .apply decode method of 'code_instance' to 'r', 
	 * by flipping all possible d_min/2 bits in d_min/2 least reliable position
	 *
	 * \param code_instance: the coding and decoding scheme
	 * \param r: received codeword, after BPSK_demodulation
	 * \param eta: number of bit to flip on LRP, default with d_min/2
	 * \return decoded result, v
	 */
	template<class code_type>
	static Matrix<GF2> decode2_v(code_type& code_instance, const Matrix<my_double>& r, int eta = -1)
	{

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

		is_ML = false;
		// a shell for code_type, such as BCH, return code_instance.decode(r);
		int n = r.size();
		Matrix<my_double> r_abs(1, n);
		for (int i = 0; i < n; ++i) {
			r_abs(i)= r(i) >= 0 ? r(i) : -r(i);
		}
		//cout << "r_abs" << r_abs;
		Matrix<int> r_ind = r_abs.sort_with_ind();		// note that 'r_abs' is sorted
		Matrix<GF2> hdr = BPSK::demodulation(r);

		int error_num;
		int d_min = code_instance.get_d();
		eta = eta == -1 ? d_min / 2 : eta;	// set eta default as d_min / 2 

		// first not flip any position
		Matrix<GF2> best_dv;
		lambda_min = r_abs(n - 1) * n;	// the upper bound for lambda_min
		best_dv = code_instance.decode_v(hdr);			// for the first set is we do demand in v space
		if (code_instance.is_in_v_space(best_dv)) {
			lambda_min = Measure::correlation_discrepancy_v(r, best_dv);
			error_num = Measure::error_num;
			if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda_min, error_num)) {
				//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
				is_ML = true;
				num_early_stop++;
				total_used_list_num++;
				return best_dv;
			}
		}

		int total_test = 1 << eta;
		//cout << "total_test=" << total_test << endl;
		// generate test vectors
		//Matrix<bool> is_flipped_mat(total_test, eta, '0');
		Matrix<GF2> flipped_hdr(hdr);
		Matrix<GF2> dv(1, n, '0');
		//Matrix<my_double> cd(1, total_test, '0');
		//cout << "r_ind" << r_ind;
		//cout << "hdr" << hdr;
		for (int i = 1; i < total_test; ++i) {
			int pos = 0;
			int flip_num = i;
			while (flip_num != 0) {
				flipped_hdr(r_ind(pos)) += flip_num & 1;
				//is_flipped_mat(i, pos) = flip_num & 1;
				//cd(i) += ((flip_num & 1) ? r_abs(pos) : 0);

				flip_num >>= 1;
				++pos;
			}
			//cout << "flipped_hdr" << flipped_hdr;
			dv = code_instance.decode_v(flipped_hdr);
			if (!code_instance.is_in_v_space(dv));
			else{	 // for later dv, we demand it in v space
				my_double lambda = Measure::correlation_discrepancy_v(r, dv);
				error_num = Measure::error_num;
				if (lambda >= lambda_min);
				else{
					// update best_dv to this erasure decoding result
					best_dv = dv;
					lambda_min = lambda;/*
					cout << "dv=" << dv;
					cout << "lambda_min=" << lambda << endl;*/
					if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda, error_num)){
						//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
						is_ML = true;
						num_early_stop++;
						total_used_list_num += (unsigned long long) i + 1;
						return best_dv;
					}
				}
			}
			flipped_hdr = hdr;
		}

		// if function like sort used, note that rand function is called and noise and bit data generated is different
		// the error rate may variate and be different from decoding method not using it (although they are same method)
		
		//Matrix<int> cd_ind = cd.sort_with_ind();
		//Matrix<bool> is_flipped(1, n, '0');
		//for (int i = 1; i < total_test; i++) {
		//	for (int j = 0; j < eta; ++j) {
		//		is_flipped(r_ind(j))= is_flipped_mat(i, j);
		//	}
		//	Matrix<GF2> dv = decode_flip_v(code_instance, hdr, is_flipped);
		//	if (!code_instance.is_in_v_space(dv))	continue; // for later dv, we demand it in v space
		//	my_double lambda = Measure::correlation_discrepancy_v(r, dv, error_num);
		//	if (lambda < lambda_min) {
		//		// update best_dv to this erasure decoding result
		//		best_dv = dv;
		//		lambda_min = lambda;/*
		//		cout << "dv=" << dv;
		//		cout << "lambda_min=" << lambda << endl;*/
		//		if (optimum_condition::is_optimum(d_min, best_dv, hdr, r_abs, r_ind, lambda, error_num)) {
		//			//cout << "optimum condition satisfied, quit testing at iter " << 1 << endl;
		//			return best_dv;
		//		}
		//	}
		//}

#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(Chase2) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Chase2) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << "(Chase2) GF2e_ope_num = " << GF2e_ope_num_after - GF2e_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		//cout << "best_dv=" << best_dv;
		total_used_list_num += total_test;
		return best_dv;
	}

	/* we do not implement algorithm 1 since it is too complex in computation*/
};

bool Chase::is_ML = false;
int Chase::num_early_stop = 0;
unsigned long long Chase::total_used_list_num = 0;
my_double Chase::lambda_min = 0;
