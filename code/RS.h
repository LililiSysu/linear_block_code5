/*****************************************************************//**
 * \file   RS.h
 * \brief  class of RS code
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"../GF/GF2e.h"
#include"../GF/polynomial.h"
#include"../GF/polynomial_2v.h"

template<class T> class Kotter_interpolation;
template<class T> class RR_root_finding;
class Lagrange_interpolation;

template<int m, int k> class RS {
private:
	int r;		// number of redundancy, not rate
	int n;		// length of code, n=2^m-1
	int d;		// minimum Hamming distance

	polynomial<GF2e<m>> gX;	// generator polynomial, length of r+1, that is max order r
	polynomial<GF2e<m>> pX;	// parity polynomial, length of k+1, that is max order k
	Matrix<GF2e<m>> generator_M;// generator matrix, size of k*n
	Matrix<GF2e<m>> parity_M;	// parity matrix, size of r*n

	Matrix<GF2e<m>> alpha;		// code locator, can be in any sequence
	Matrix<GF2e<m>> generator_M_by_code_locators;	// generator matrix, size of k*n
	Matrix<GF2e<m>> decoder_M_by_code_locators;		// decoder matrix, size of k*k

	Matrix<GF2e<m>> generator_M_systematic;	// systematic generator matrix, size of k*n, we generate it by Lagrange interpolation polynomial

	/**
	 * . for detail, SEE "Shift-Register Synthesis and BCH Decoding"
	 *
	 * \param Sn: the target to generate
	 * \return the coefficients of LFSR that generate Sn
	 */
	polynomial<GF2e<m>> LFSR_Synthesis_Algorithm(const polynomial<GF2e<m>>& Sn) {
		Matrix<GF2e<m>> S = Sn.get_coeff();
		polynomial<GF2e<m>> C(1, '1');
		polynomial<GF2e<m>> B(1, '1');
		int x = 1;		// n-m, interval between new change and the last change,
		int L = 0;		// minimum length
		GF2e<m> b = 1;	// discrepancy since last change
		int N = 0;		// the number of 'Sn' processed
		GF2e<m> d;

		int nn = Sn.size();
		while (N != nn) {
			d = S(N);
			for (int i = 1; i <= L; ++i) {
				if (i < C.size()) {
					d = d + C(i) * S(N - i);
				}
			}
			if (d == 0) {
				x++;
			}
			else if (2 * L > N) {
				C = C - d / b * B.get_shift_right(x);
				x++;
			}
			else {
				polynomial<GF2e<m>> T(C);
				C = C - d / b * B.get_shift_right(x);
				L = N + 1 - L;
				B = T;
				b = d;
				x = 1;
			}
			N++;
		}
		C.simplify();

		return C;
	}

public:
	Matrix<int> information_set_ind;
	Matrix<int> redundancy_set_ind;
	Matrix<GF2e<m>> generator_M_systematic_any_pos;	// systematic generator matrix with any k position be information set
	Matrix<GF2e<m>> generator_compute_up_store_fix;
	Matrix<GF2e<m>> generator_compute_down_store;

	// systematic generator matrix with any n-k position be information set, for dual RS code
	Matrix<GF2e<m>> generator_M_systematic_any_pos_dual;

	RS() {
		n = (1 << m) - 1;
		r = n - k;
		d = n - k + 1;	// MDS property of RS code
		Matrix<GF2e<m>> tmp(1, 2, { 1,0 });
		gX = tmp;
		//cout << "gX=" << gX;
		tmp(1)= 1;
		for (int i = 1; i <= n - k; ++i) {
			tmp(0).set_by_alpha_power(i);
			//tmp(0)= -tmp(0);
			gX = gX * tmp;
			//cout << "gX=" << gX;
		}
		polynomial<GF2e<m>> xup(n + 1, '0');
		xup(0) = 1;
		xup(n) = 1;

		pX = xup / gX;
		Matrix<GF2e<m>> revp = pX.get_coeff();
		revp.rev();
		pX = revp;		// reverse

		int glen = gX.size();
		generator_M = Matrix<GF2e<m>>(k, n, '0');
		for (int i = 0; i < k; ++i) 
			for (int j = 0; j < glen; ++j)
				generator_M(i * n + i + j) = gX(j);

		int plen = pX.size();
		parity_M = Matrix<GF2e<m>>(r, n, '0');
		for (int i = 0; i < r; ++i)
			for (int j = 0; j < plen; ++j)
				parity_M(i * n + i + j) = pX(j);

		// the code locator of how we encode an RS code
		alpha = Matrix<GF2e<m>>(1, n);
		for (int i = 0; i < n; ++i) {
			alpha(i).set_by_alpha_power(i);		// start from 1, and go to a^{n}=a^{q-1}
		}

		// the corresponding generator matrix. Note that the encoding mapping 
		// by this generator matrix is different from the one by cyclic generator polynomial
		generator_M_by_code_locators = Matrix<GF2e<m>>(k, n, '0');
		GF2e<m> sigma_0 = 2;
		for (int i = 0; i < k; ++i)
			for (int j = 0; j < n; ++j)
				generator_M_by_code_locators(i, j) = pow(sigma_0, i * j);

		decoder_M_by_code_locators = generator_M_by_code_locators.get_part(0, 0, k - 1, k - 1).inv();

		// generator matrix of a different RS code, by Lagrange interpolation polynomial
		// the advantage of this is Gaussian elemination free
		generator_M_systematic = Matrix<GF2e<m>>(k, n, 'i');
		Matrix<GF2e<m>> alpha(1, k);
		for (int i = 0; i < k; ++i) {
			alpha(i).set_by_alpha_power(i);
		}
		for (int i = 0; i < k; ++i) {
			polynomial<GF2e<m>> p = Lagrange_interpolation::generate(alpha, i);

			GF2e<m> tmp_for_new_GF2e;
			for (int j = k; j < n; ++j) {
				tmp_for_new_GF2e.set_by_alpha_power(j);
				generator_M_systematic(i, j) = p.evaluate(tmp_for_new_GF2e);
			}
		}
		//cout << "generator_M_systematic" << generator_M_systematic;
	}

	~RS() {}


	inline int get_k() const {
		return k;
	}

	// number of redundancy
	inline int get_r() const {
		return r;
	}

	inline int get_n() const {
		return n;
	}

	inline int get_d() const {
		return d;
	}

	Matrix<GF2e<m>> get_generator_matrix() const {
		// size of k*n matrix
		return generator_M;
	}

	Matrix<GF2e<m>> get_parity_matrix() const {
		// size of r*n matrix
		return parity_M;
	}

	Matrix<GF2e<m>> get_generator_matrix_by_code_locators() const {
		// size of k*n matrix
		return generator_M_by_code_locators;
	}

	void print_info() const {
		cout << "RS (" << n << "," << k << "," << d << ") code, rate = " << k / (double)n << endl;
	}

	// check if in v space
	bool is_in_v_space(const Matrix<GF2e<m>>& v) const {
		polynomial<GF2e<m>> result = v % gX;
		if (result.size() == 1 && result(0) == 0) {
			return true;
		}
		else {
			return false;
		}
	}

	/* encode scheme by generator matrix, where v2u and u2v is fast */

	// u2c
	Matrix<GF2e<m>> encode(const Matrix<GF2e<m>>& u) {

		/* encode to codeword by generator polynomial */
		polynomial<GF2e<m>> pc = u * gX;
		pc.format_len(n);
		return pc.get_coeff();
		/*for (int i = 0; i < n; ++i) {
			c(i)= pu.evaluate(i + 1);
		}

		return c;*/
	}

	// decode by PGZ method
	Matrix<GF2e<m>> decode_naive(const Matrix<GF2e<m>>& hdr) {
		return v2u(decode_naive_v(hdr));
	}

	// PGZ method by solving linear equations, decode to u
	Matrix<GF2e<m>> decode_naive_v(const Matrix<GF2e<m>>& r) {
		polynomial<GF2e<m>> pr(r);			// store 'r' as polynomial

		/* syndrome based decode */

		/* compute syndrome S1, ..., S2t by c(\sigma), ..., c(\sigma^ {n - k}) */
		Matrix<GF2e<m>> S(1, n - k);
		for (int i = 0; i < n - k; ++i) {
			S(i)= pr.evaluate(i + 2);
		}
		//cout << "S=" << S;
		bool error_free = true;
		for (int i = 0; i < n - k; ++i) {
			if (S(i) != 0) {
				error_free = false;
				break;
			}
		}
		if (error_free) {
			//polynomial<GF2e<m>> pu = pr / gX;
			//cout << "pu=" << pu;
			//cout << "mod=" << pr % gX;		// no problem
			//pu.format_len(k);
			return r;
		}

		/* solve equatinos to find error locations */
		int t = (n - k) / 2;	// 2t=n-k£¬we need to determine t, number of errors

		Matrix<GF2e<m>> Lambda_left_tmp(t, t);		// left matrix of equations to find Lambda, the error locations
		for (int i = 0; i < t; ++i) {
			for (int j = 0; j < t; ++j) {
				Lambda_left_tmp(i, j)= S(i + j);
			}
		}
		t = Lambda_left_tmp.rank();
		//Lambda_left_tmp.clear();					// naive, doesn't matter

		Matrix<GF2e<m>> Lambda_left(t, t);
		for (int i = 0; i < t; ++i) {
			for (int j = 0; j < t; ++j) {
				Lambda_left(i, j)= S(i + j);
			}
		}

		//cout << "Lambda_left_tmp=" << Lambda_left_tmp;
		//Lambda_left = Lambda_left_tmp;
		// if Lambda_left is singular, which means y=0 and the second equation is not needed
		// only one place is error, the Lambda_left need to shrink down to the left-up part as
		// well as Lambda_right, until the matrix Lambda_left is non-singular and equation can be solved


		// Note !!!!!!!!
		// we have to determine t in this while loop which is time wasting

		Matrix<GF2e<m>> Lambda_right(t, 1);		// right matrix of equations to find Lambda, the error locations
		for (int i = 0; i < t; ++i) {
			Lambda_right(i)= -S(i + t);
		}
		//cout << "Lambda_right=" << Lambda_right;

		Matrix<GF2e<m>> Lambda = Lambda_left.inv() * Lambda_right;
		//cout << "Lambda=" << Lambda;

		/* form the locator polynomial */
		polynomial<GF2e<m>> p_locator(t + 1, '0');
		p_locator(0)= 1;					// Lambda 0 is sure to be 1
		for (int i = 1; i <= t; ++i) {
			p_locator(i)= Lambda(t - i);
		}
		//cout << "p_locator=" << p_locator;

		/* find roots in locator polynomial */
		Matrix<GF2e<m>> root_pl = p_locator.find_roots(1 << m);
		//cout << "root_pl=" << root_pl;
		if (root_pl.size()==0) {
			//cout << "encounter uncorrectable errors" << endl;
			return Matrix<GF2e<m>>(1, n, '0');
		}
		/* compute root^{-1}, then these are the locator!!! */
		int root_len = root_pl.size();
		for (int i = 0; i < root_len; ++i) {
			root_pl(i)= 1 / root_pl(i);
		}

		/* solve equatinos to find error magnitude */
		Matrix<GF2e<m>> Magnitude_left(t, t);		// left matrix of equations to find Magnitude, the error locations
		for (int j = 0; j < t; ++j) {
			Magnitude_left(j)= root_pl(j);
		}
		for (int i = 1; i < t; ++i) {
			for (int j = 0; j < t; ++j) {
				Magnitude_left(i, j)= Magnitude_left(i - 1, j) * root_pl(j);
			}
		}
		//cout << "Magnitude_left=" << Magnitude_left;

		Matrix<GF2e<m>> Magnitude_right(t, 1);		// right matrix of equations to find Magnitude, the error locations
		for (int i = 0; i < t; ++i) {
			Magnitude_right(i)= -S(i);
		}
		//cout << "Magnitude_right=" << Magnitude_right;

		Matrix<GF2e<m>> Magnitude = Magnitude_left.inv() * Magnitude_right;
		//cout << "Magnitude=" << Magnitude;

		/* correct errors in r */
		Matrix<GF2e<m>> corrected_r(r);
		int error_pos_len = root_pl.size();
		for (int i = 0; i < error_pos_len; ++i) {
			int error_pos = (int)root_pl(i) - 1;
			corrected_r(error_pos)= corrected_r(error_pos) + Magnitude(i);
		}

		//polynomial<GF2e<m>> pu = corrected_r / gX;
		//cout << "pu=" << pu;
		//cout << "corrected_r=" << corrected_r;
		//pu.format_len(k);
		//return pu.get_coeff();
		return corrected_r;
	}

	// BM decoding as a main decode method, decode to u
	Matrix<GF2e<m>> decode_BM(const Matrix<GF2e<m>>& hdr) {
		return v2u(decode_BM_v(hdr));
	}

	// decode to v
	Matrix<GF2e<m>> decode_BM_v(const Matrix<GF2e<m>>& r) {
		//int n_start = GF2e_auxiliary_storage::num_operations;
		polynomial<GF2e<m>> pr(r);			// store 'r' as polynomial

		/* syndrome based decode */

		/* compute syndrome S1, ..., S2t by c(\sigma), ..., c(\sigma^ {n - k}) */
		Matrix<GF2e<m>> S(1, n - k);
		GF2e<m> tmp_for_new_GF2e;
		for (int i = 0; i < n - k; ++i) {
			tmp_for_new_GF2e.set_by_alpha_power(i + 1);
			S(i) = pr.evaluate(tmp_for_new_GF2e);
		}
		//cout << "S" << S;
		//int n_syndrome = GF2e_auxiliary_storage::num_operations;
		//cout << "syndrome_operations=" << n_syndrome - n_start << endl;
		bool error_free = true;
		for (int i = 0; i < n - k; ++i) {
			if (S(i) != 0) {
				error_free = false;
				break;
			}
		}
		if (error_free) {
			//n_start = GF2e_auxiliary_storage::num_operations;
			//polynomial<GF2e<m>> pu = pr / gX;
			//cout << "pu=" << pu;
			//cout << "mod=" << pr % gX;		// no problem

			//int n_divide_gX = GF2e_auxiliary_storage::num_operations;
			//cout << "divide_gX_operations=" << n_divide_gX - n_start << endl;
			//pu.format_len(k);
			return r;
		}

		polynomial<GF2e<m>> pS(S);			// store 'S' as polynomial

		//cout << "pS=" << pS;
		polynomial<GF2e<m>> p_locator = LFSR_Synthesis_Algorithm(pS);
		//cout << "p_locator=" << p_locator;
		Matrix<GF2e<m>> root_pl = p_locator.find_roots(1 << m);				// inverse of error locator

		//cout << "root_pl=" << root_pl;
		if (root_pl.size()==0) {
			//cout << "encounter uncorrectable errors" << endl;
			return Matrix<GF2e<m>>(1, n, '0');
		}
		int t = root_pl.size();

		/* using Forney's alogrithm to compute error magnitude */
		polynomial<GF2e<m>> X_2t(2 * t + 1, '0');
		X_2t(2 * t)= 1;
		polynomial<GF2e<m>> Omega_x = pS * p_locator % X_2t;	// the key equation

		//cout << "Omega_x=" << Omega_x;
		polynomial<GF2e<m>> d_p_locator = p_locator.get_derivative();
		//cout << "d_p_locator=" << d_p_locator;

		Matrix<GF2e<m>> Magnitude(1, t);
		for (int i = 0; i < t; ++i) {
			//cout << "Omega_x.evaluate(root_pl(i))=" << Omega_x.evaluate(root_pl(i)) << endl;
			//cout << "d_p_locator.evaluate(root_pl(i))=" << d_p_locator.evaluate(root_pl(i)) << endl;
			if (d_p_locator.evaluate(root_pl(i)) == 0) {
				// decode error
				return Matrix<GF2e<m>>(1, n, '0');
			}

			Magnitude(i) = -Omega_x.evaluate(root_pl(i)) / d_p_locator.evaluate(root_pl(i));
		}
		//cout << "Magnitude=" << Magnitude;

		/* correct errors in r */
		Matrix<GF2e<m>> corrected_r(r);
		int error_pos_len = root_pl.size();
		tmp_for_new_GF2e = 1;
		for (int i = 0; i < error_pos_len; ++i) {
			int error_pos = GF2e_auxiliary_storage::alpha_table[(int)(tmp_for_new_GF2e / root_pl(i))];
			//cout << "error_pos = " << error_pos << endl;
			corrected_r(error_pos)= corrected_r(error_pos) + Magnitude(i);
		}
		//cout << "corrected_r" << corrected_r;
		return corrected_r;
	}

	// decode to u
	Matrix<GF2e<m>> v2u(const Matrix<GF2e<m>>& v) {
		polynomial<GF2e<m>> result = v / gX;
		result.format_len(k);
		return result.get_coeff();
	}

	friend class GMD;
	friend class Chase;
	friend class OSD;

	// function with dual code
	Matrix<GF2e<m>> encode_dual(const Matrix<GF2e<m>>& u_perp) const {

		/* encode to codeword by parity polynomial */
		polynomial<GF2e<m>> pc = u_perp * pX;
		pc.format_len(n);
		return pc.get_coeff();
	}

	bool is_in_dual_space(const Matrix<GF2e<m>>& v_perp) const {

		polynomial<GF2e<m>> result = v_perp % pX;
		return result.size() == 1 && result(0) == 0;
	}

	Matrix<GF2e<m>> v2u_dual(const Matrix<GF2e<m>>& v_perp) const {

		polynomial<GF2e<m>> result = v_perp / pX;
		result.format_len(r);
		return result.get_coeff();
	}


	// GS algorithm

	/* encode scheme by code locators, where v2u and u2v is slow */

	// well we need a different encoding, maybe can merge into the generator polynomial encoding in the future
	Matrix<GF2e<m>> encode_by_code_locators(const Matrix<GF2e<m>>& u) {
		Matrix<GF2e<m>> result(1, n, '0');
		polynomial<GF2e<m>> pu(u);
		GF2e<m> tmp_for_new_GF2e;
		for (int i = 0; i < n; ++i) {
			tmp_for_new_GF2e.set_by_alpha_power(i);
			result(i) = pu.evaluate(tmp_for_new_GF2e);
		}
		return result;
	}

	Matrix<GF2e<m>> v2u_by_code_locators(const Matrix<GF2e<m>>& r) {
		return r.get_part(0, 0, 0, k - 1) * decoder_M_by_code_locators;		// this method is naive and cannot correct any error
	}

	Matrix<GF2e<m>> decode_naive_by_code_locators(const Matrix<GF2e<m>>& hdr) {
		return v2u_by_code_locators(decode_naive_v(hdr));
	}

	Matrix<GF2e<m>> decode_BM_by_code_locators(const Matrix<GF2e<m>>& hdr) {
		//cout << "decode_v(hdr)" << decode_v(hdr);
		return v2u_by_code_locators(decode_BM_v(hdr));
	}

	int compute_lm(int multiplicity) {
		int C = n * my::n_choose_k(multiplicity + 1, 2);
		int result = 1;
		while (ipair::B(result, k - 1) <= C) {
			result++;		// less than 100, very fast
		}
		return result - 1;
	}

	int compute_tau(int multiplicity) {
		int C = n * my::n_choose_k(multiplicity + 1, 2);
		int result = 1;
		while (ipair::A(result, k - 1) <= C) {
			result++;		// may be several thousands, but less than 1 second to finish
		}
		return n - 1 - (result - 1) / multiplicity;
	}

	/**
	 * .GS algorithm with multiplicity of "multiplicity"
	 */
	Matrix<GF2e<m>> decode_GS_by_code_locators(const Matrix<GF2e<m>>& hdr, int multiplicity) {
		int L = compute_lm(multiplicity);
		//cout << "L=" << L << endl;

		//clock_t start, end;
		//start = clock();	//start timing

		polynomial_2v<GF2e<m>> Q0 = Kotter_interpolation<GF2e<m>>::solve_for_GF2e(alpha, hdr, multiplicity, L, k - 1);
		//cout << "Q0 = " << Q0;
		
		// Debug mode: about 5 times slower than polynomial_2v, Release mode: 7 times slower
		//polynomial_2v_map<GF2e<m>> Q1 = Kotter_interpolation<GF2e<m>>::solve_for_GF2e_map(alpha, hdr, multiplicity, L, k - 1);
		//cout << "Q1 = " << Q1;

		// Debug mode: about 3.3 times slower than polynomial_2v, Release mode: 1.7 times slower, 
		// for m=16, Release mode: 3 times slower, much more memory (4 times at least)
		
		//polynomial_2v_umap<GF2e<m>> Q2 = Kotter_interpolation<GF2e<m>>::solve_for_GF2e_umap(alpha, hdr, multiplicity, L, k - 1);
		//cout << "Q2 = " << Q2;

		//end = clock();		//end timing
		//cout.setf(ios::scientific);
		//cout << scientific << "kotter time consume: " << ((double)end - start) / CLOCKS_PER_SEC << "s\n\n";
		//cout.unsetf(ios::floatfield);

		if (Q0.size() == 0) {	// no Q0 foundd
			return Matrix<GF2e<m>>(1, k, '0');
		}

		RR_root_finding<GF2e<m>> problem_RR;
		problem_RR.set_input(Q0, k, 1 << m);

		//clock_t start, end;
		//start = clock();	//start timing

		problem_RR.solve_for_GF2e();

		//end = clock();		//end timing
		//cout.setf(ios::scientific);
		//cout << scientific << "RR time consume: " << ((double)end - start) / CLOCKS_PER_SEC << "s\n\n";
		//cout.unsetf(ios::floatfield);

		//cout << "problem_RR.y_roots" << problem_RR.y_roots;

		int root_size = (int)problem_RR.y_roots.size();
		if (root_size == 0) {
			// fail to decode
			return Matrix<GF2e<m>>(1, k, '0');
		}
		else if (root_size == 1) {
			int sss = problem_RR.y_roots(0).size();
			return problem_RR.y_roots(0).combine_right(Matrix<GF2e<m>>(1, k - sss, '0'));
		}
		else {
			// if root_size > 1, re-encode the polynomial, choose the closest one with received vector
			Matrix<GF2e<m>> encoded_back = encode_by_code_locators(problem_RR.y_roots(0));
			int min_distance = encoded_back.Hamming_distance(hdr);
			int min_ind = 0;
			for (int i = 1; i < root_size; ++i) {
				encoded_back = encode_by_code_locators(problem_RR.y_roots(i));
				int can_distance = encoded_back.Hamming_distance(hdr);
				/*cout << "min_distance: " << min_distance << endl;b
				cout << "can_distance: " << can_distance << endl;*/
				if (can_distance < min_distance) {
					min_distance = can_distance;
					min_ind = i;
				}
			}
			int sss = problem_RR.y_roots(min_ind).size();
			return problem_RR.y_roots(min_ind).combine_right(Matrix<GF2e<m>>(1, k - sss, '0'));
		}
	}

	/* systematic encoding by Lagrange interpolation polynomials */

	Matrix<GF2e<m>> encode_by_Lagrange_interpolation(const Matrix<GF2e<m>>& u) {
		Matrix<GF2e<m>> result(1, n, '0');
		for (int i = 0; i < k; ++i) {
			result(i) = u(i);
		}
		for (int i = k; i < n; ++i) {
			for (int j = 0; j < k; ++j) {
				result(i) += u(j) * generator_M_systematic(j, i);
			}
		}
		return result;
	}

	Matrix<GF2e<m>> v2u_by_Lagrange_interpolation(const Matrix<GF2e<m>>& r) {
		return r.get_part(0, 0, 0, k - 1);		// this method is naive and cannot correct any error
	}

	/**
	 * .information set should be size of 1*k, and with each different element in {0,1,...,n-1}
	 */
	void generate_systematic_generator_any_pos_naive(const Matrix<int>& _information_set_ind) {

		information_set_ind = _information_set_ind;		// must be size of 1*k
		information_set_ind.sort();
		
		//cout << "information_set_ind_naive" << information_set_ind;
		redundancy_set_ind = Matrix<int>(1, n, 'N');
		redundancy_set_ind = redundancy_set_ind.erase_cols(information_set_ind);		// this should be ordered with '<'
		//cout << "redundancy_set_ind" << redundancy_set_ind;

		// generator matrix of a different RS code, by Lagrange interpolation polynomial
		// the advantage of this is Gaussian elemination free
		generator_M_systematic_any_pos = Matrix<GF2e<m>>(k, n, '0');
		Matrix<GF2e<m>> alpha(1, k);
		for (int i = 0; i < k; ++i) {
			alpha(i).set_by_alpha_power(information_set_ind(i));
		}

		// further simplified, knowing the position of 1 and 0

		/*cout << "GF2e_auxiliary_storage::operation_number: before polynomial generation = " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/

		vector<polynomial<GF2e<m>>> p(k);
		for (int i = 0; i < k; ++i) {
			p[i] = Lagrange_interpolation::generate(alpha, i);
		}
		/*cout << "GF2e_auxiliary_storage::operation_number: after polynomial generation = " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/

		int info_ind = 0;
		GF2e<m> tmp_for_new_GF2e;
		for (int j = 0; j < n; ++j) {
			if (info_ind < k && j == information_set_ind(info_ind)) {
				for (int i = 0; i < k; ++i) {
					//cout << "ahpha" << alpha;
					generator_M_systematic_any_pos(i, j) = i == info_ind;
				}
				info_ind++;
			}
			else {
				for (int i = 0; i < k; ++i) {
					//cout << "ahpha" << alpha;
					tmp_for_new_GF2e.set_by_alpha_power(j);
					generator_M_systematic_any_pos(i, j) = p[i].evaluate(tmp_for_new_GF2e);
				}
			}
		}
		//cout << "generator_M_systematic_any_pos" << generator_M_systematic_any_pos;
	}

	/**
	 * .information set should be size of 1*k, and with each different element in {0,1,...,n-1}
	 * if k < r, this cause more computation than naive method
	 */
	void generate_systematic_generator_any_pos(const Matrix<int>& _information_set_ind, bool is_sorted = false) {

		information_set_ind = _information_set_ind;		// must be size of 1*k
		if (is_sorted);
		else {
			information_set_ind.sort();
		}
		
		//cout << "information_set_ind" << information_set_ind;
		redundancy_set_ind = Matrix<int>(1, n, 'N');
		redundancy_set_ind = redundancy_set_ind.erase_cols(information_set_ind, is_sorted);		// this should be ordered with '<'
		
		//cout << "redundancy_set_ind" << redundancy_set_ind;

		// generator matrix of a different RS code, by Lagrange interpolation polynomial
		// the advantage of this is Gaussian elemination free
		generator_M_systematic_any_pos.resize(k, n, false);

		// further simplified, knowing the position of 1 and 0

		/*cout << "GF2e_auxiliary_storage::operation_number: before polynomial generation = " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/

		/*cout << "GF2e_auxiliary_storage::operation_number: after polynomial generation = " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/

		int info_ind = 0;
		for (int j = 0; j < n; ++j) {
			if (info_ind < k && j == information_set_ind(info_ind)) {
				for (int i = 0; i < k; ++i) {
					//cout << "ahpha" << alpha;
					generator_M_systematic_any_pos(i, j) = i == info_ind;
				}
				info_ind++;
			}
			else {
				for (int i = 0; i < k; ++i) {
					//cout << "ahpha" << alpha;
					generator_M_systematic_any_pos(i, j) = H_u(i, j);
				}
			}
		}
		//cout << "generator_M_systematic_any_pos" << generator_M_systematic_any_pos;
	}

	GF2e<m> H_u(int i, int j) {
		int i_prime = information_set_ind(i);
		GF2e<m> result(1);
		GF2e<m> alpha_i_prime;
		alpha_i_prime.set_by_alpha_power(i_prime);
		GF2e<m> alpha_j;
		alpha_j.set_by_alpha_power(j);

		if (k < r) {		// less computation, can be further optimized via storage
			int information_set_ind_size = information_set_ind.size();
			for (int w = 0; w < information_set_ind_size; ++w) {
				int j_prime = information_set_ind(w);
				if (j_prime != i_prime) {
					GF2e<m> alpha_j_prime;
					alpha_j_prime.set_by_alpha_power(j_prime);
					result *= (alpha_j - alpha_j_prime) / (alpha_i_prime - alpha_j_prime);
				}
			}
		}
		else {
			int redundancy_set_ind_size = redundancy_set_ind.size();
			for (int w = 0; w < redundancy_set_ind_size; ++w) {
				int j_prime = redundancy_set_ind(w);
				if (j_prime != j) {		// i_prime will never equals to j_prime
					/*if (i_prime == j_prime) {
						cout << "(i_prime == j_prime) out of expectation" << endl;
					}*/

					GF2e<m> alpha_j_prime;
					alpha_j_prime.set_by_alpha_power(j_prime);
					result *= (alpha_i_prime - alpha_j_prime) / (alpha_j - alpha_j_prime);
				}
			}
			result *= alpha_i_prime / alpha_j;
		}

		return result;
	}

	/**
	 * .information_set_ind must be set, efficient for k>r, i.e., rate > 1/2
	 * 
	 */
	void compute_generator_store() {

		GF2e<m> alpha_i_prime;
		GF2e<m> alpha_j;
		GF2e<m> alpha_j_prime;

		if (k > r) {
			generator_compute_up_store_fix.resize(1, k, false);
			generator_compute_down_store.resize(1, r, false);

			// for k>r, i.e., rate>1/2
			for (int i = 0; i < k; ++i) {
				alpha_i_prime.set_by_alpha_power(information_set_ind(i));
				generator_compute_up_store_fix(i) = alpha_i_prime;
				for (int w = 0; w < r; ++w) {		
					// can be further optimized by switching the for loop, accelerate computation speed but not saving operation times
					alpha_j_prime.set_by_alpha_power(redundancy_set_ind(w));
					generator_compute_up_store_fix(i) *= alpha_i_prime - alpha_j_prime;
				}
			}

			for (int j = 0; j < r; ++j) {
				alpha_j.set_by_alpha_power(redundancy_set_ind(j));
				generator_compute_down_store(j) = alpha_j;
				for (int w = 0; w < r; ++w) {
					if (j != w) {
						alpha_j_prime.set_by_alpha_power(redundancy_set_ind(w));		// this seems ugly, but fixing it cost a lot of force
						generator_compute_down_store(j) *= alpha_j - alpha_j_prime;
					}
				}
			}
		}
		else {
			generator_compute_up_store_fix.resize(1, r, false);
			generator_compute_down_store.resize(1, k, false);

			// for k<=r, i.e., rate<=1/2
			for (int i = 0; i < r; ++i) {
				alpha_i_prime.set_by_alpha_power(redundancy_set_ind(i));
				alpha_j_prime.set_by_alpha_power(information_set_ind(0));
				generator_compute_up_store_fix(i) = alpha_i_prime - alpha_j_prime;
				for (int w = 1; w < k; ++w) {
					alpha_j_prime.set_by_alpha_power(information_set_ind(w));
					generator_compute_up_store_fix(i) *= alpha_i_prime - alpha_j_prime; 
				}
			}

			// we should have k > 1
			alpha_j.set_by_alpha_power(information_set_ind(0));
			alpha_j_prime.set_by_alpha_power(information_set_ind(1));
			generator_compute_down_store(0) = alpha_j - alpha_j_prime;
			for (int w = 2; w < k; ++w) {
				alpha_j_prime.set_by_alpha_power(information_set_ind(w));
				generator_compute_down_store(0) *= alpha_j - alpha_j_prime;
			}

			for (int j = 1; j < k; ++j) {
				alpha_j.set_by_alpha_power(information_set_ind(j));
				alpha_j_prime.set_by_alpha_power(information_set_ind(0));
				generator_compute_down_store(j) = alpha_j - alpha_j_prime;
				for (int w = 1; w < k; ++w) {
					if (w != j) {
						alpha_j_prime.set_by_alpha_power(information_set_ind(w));
						generator_compute_down_store(j) *= alpha_j - alpha_j_prime;
					}
				}
			}
		}		
	}

	/**
	 * .information set should be size of 1*k, and with each different element in {0,1,...,n-1}
	 */
	void generate_systematic_generator_any_pos_best(const Matrix<int>& _information_set_ind, bool is_sorted = false) {

		information_set_ind = _information_set_ind;		// must be size of 1*k
		if (is_sorted);
		else {
			information_set_ind.sort();
		}

		//cout << "information_set_ind" << information_set_ind;
		redundancy_set_ind = Matrix<int>(1, n, 'N');
		redundancy_set_ind = redundancy_set_ind.erase_cols(information_set_ind, is_sorted);		// this should be ordered with '<'

		//cout << "redundancy_set_ind" << redundancy_set_ind;

		compute_generator_store();

		// generator matrix of a different RS code, by Lagrange interpolation polynomial
		// the advantage of this is Gaussian elemination free
		generator_M_systematic_any_pos.resize(k, n, false);

		// further simplified, knowing the position of 1 and 0

		/*cout << "GF2e_auxiliary_storage::operation_number: before polynomial generation = " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/

		/*cout << "GF2e_auxiliary_storage::operation_number: after polynomial generation = " \
			<< GF2e_auxiliary_storage::operation_number << endl;*/

		int info_ind = 0;
		GF2e<m> alpha_j;
		GF2e<m> alpha_i;

		if (k > r) {
			for (int j = 0; j < n; ++j) {
				if (info_ind < k && j == information_set_ind(info_ind)) {
					for (int i = 0; i < k; ++i) {
						//cout << "ahpha" << alpha;
						generator_M_systematic_any_pos(i, j) = i == info_ind;
					}
					info_ind++;
				}
				else {
					alpha_j.set_by_alpha_power(redundancy_set_ind(j - info_ind));
					for (int i = 0; i < k; ++i) {
						alpha_i.set_by_alpha_power(information_set_ind(i));
						//cout << "ahpha" << alpha;
						generator_M_systematic_any_pos(i, j) = generator_compute_up_store_fix(i) / \
							generator_compute_down_store(j - info_ind) / (alpha_j - alpha_i);
					}
				}
			}
		}
		else {
			for (int j = 0; j < n; ++j) {
				if (info_ind < k && j == information_set_ind(info_ind)) {
					for (int i = 0; i < k; ++i) {
						//cout << "ahpha" << alpha;
						generator_M_systematic_any_pos(i, j) = i == info_ind;
					}
					info_ind++;
				}
				else {
					alpha_j.set_by_alpha_power(redundancy_set_ind(j - info_ind));
					for (int i = 0; i < k; ++i) {
						alpha_i.set_by_alpha_power(information_set_ind(i));
						//cout << "ahpha" << alpha;
						generator_M_systematic_any_pos(i, j) = generator_compute_up_store_fix(j - info_ind) / \
							generator_compute_down_store(i) / (alpha_j - alpha_i);
					}
				}
			}
		}
		//cout << "generator_M_systematic_any_pos" << generator_M_systematic_any_pos;
	}

	Matrix<GF2e<m>> encode_by_Lagrange_interpolation_any_pos(const Matrix<GF2e<m>>& u) {
		Matrix<GF2e<m>> result(1, n, '0');
		for (int i = 0; i < k; ++i) {
			result(information_set_ind(i)) = u(i);
		}
		for (int i = 0; i < n - k; ++i) {
			int result_ind = redundancy_set_ind(i);
			for (int j = 0; j < k; ++j) {
				result(result_ind) += u(j) * generator_M_systematic_any_pos(j, result_ind);
			}
		}
		return result;
		//return u * generator_M_systematic_any_pos;
	}

	Matrix<GF2e<m>> v2u_by_Lagrange_interpolation_any_pos(const Matrix<GF2e<m>>& r) {
		return r.get_cols(information_set_ind);
	}

	/**
	 * .(*this) should be (n,n-k) RS code, we generate the Generator matrix of dual code of (n,k) RS code
	 * 
	 */
	void generate_systematic_generator_any_pos_dual(const Matrix<int>& _information_set_ind, bool is_sorted = false) {

		generate_systematic_generator_any_pos_best(_information_set_ind, is_sorted);
		generator_M_systematic_any_pos_dual = generator_M_systematic_any_pos;

		int j_i;
		int p;
		GF2e<m> tmp;
		for (int i = 0; i < k; ++i) {
			j_i = information_set_ind(i);
			// r si the size of theta^c
			for (int pp = 0; pp < r; ++pp) {
				p = redundancy_set_ind(pp);
				tmp.set_by_alpha_power(p - j_i);
				generator_M_systematic_any_pos_dual(i, p) *= tmp;
			}
		}
	}
};

// well, use a new class

template<class T> class Kotter_interpolation {
public:

	/**
	 * .Kotter's interpolation algorithm, can set each point's multiplicity differently, general use
	 *
	 * \param alpha: vector of points' x-value
	 * \param beta: vector of points' y-value
	 * \param ms: vector of multiplicity at each point
	 * \param L: the max degree of y, the output list will contain at most L polynomials
	 * \param v: will find the polynomials of minimun (1,v)-revlex order
	 * 
	 * \return the interpolation polynomial
	 */
	static polynomial_2v<T> solve(const Matrix<T>& alpha, const Matrix<T>& beta, const Matrix<int>& ms, int L, int v) {		
		int n = alpha.size();
		Matrix<polynomial_2v<T>> g(1, L + 1);	// use matrix of matrix
		for (int j = 0; j <= L; ++j) {
			g(j) = polynomial_2v<T>(1, j + 1, '0');
			g(j)(0, j) = 1;
		}

		//cout << "g" << g;		// the initial polynomail, by iteration, the polymonial will satisfy all the multiplicity condition

		Matrix<T> Delta(1, L + 1, '0');	// discrepancy

		for (int i = 0; i < n; ++i) {

			// (ms-1,1)-lex order, if ms(i) is invariant of i, the following lines can be removed out of the loop
			// if we enforce the conditions D_{r,s}(alpha,beta)=0 for r+s<ms(i) in 
			// an order in which (r-1,s) always precedes (r,s), the cumulative kernels 
			// will be F(x)-modules, interpolation admits of a less complex solution. (McEliece,p23)
			ipair::weighted_i = (ms(i) - 1) == 0 ? 1 : (ms(i) - 1);		// prevent ms(i)==1, which induce weighted_i be 0
			ipair::weighted_j = 1;
			ipair::is_lex = true;

			ipair top_p(ms(i) - 1, 0);
			Matrix<ipair> sorted_ipair = ipair::under_top(top_p, ms(i));
				// find the ipair from (0,0) to (ms(i) - 1, 0), under constrain r+s<ms(i)

			int sorted_ipair_size = sorted_ipair.size();			// sort them under (ms(i)-1,1)-lex order

			//cout << "sorted_ipair_size=" << sorted_ipair_size << endl;

			for (int p = 0; p < sorted_ipair_size; ++p) {
				int r = sorted_ipair(p).i;
				int s = sorted_ipair(p).j;
				Matrix<int> J(1, 0, 'v');
				for (int j = 0; j <= L; ++j) {
					Delta(j) = g(j).Hasse_D(r, s, alpha(i), beta(i));	// jth discrepancy
					if (Delta(j) != 0) {
						J.push_back(j);		// store the index where discrepency is not 0, as J
					}
				}
				/*cout << "---------------- i = " << i << " ---------------" << endl;
				cout << "g" << g;
				cout << "Delta" << Delta;*/

				if (J.size() != 0) {

					// change it back to (1,k-1) wdeg monomial order, revlex, we compare polynomials under this order
					ipair::weighted_i = 1;
					ipair::weighted_j = v;
					ipair::is_lex = false;

					int Js = J.size();
					ipair min_monomial = g(J(0)).leading_monomial();
					int j_star = 0;
					for (int j = 1; j < Js; ++j) {
						ipair can_monomial = g(J(j)).leading_monomial();
						if (min_monomial < can_monomial);
						else {
							min_monomial = can_monomial;
							j_star = j;		// j_star=arg min{g_j:j\in J}
						}
					}
					polynomial_2v<T> f = g(J(j_star));		// as the book says

					// update g(j), make it in pass the point (alpha(i),beta(i))
					for (int j = 0; j < Js; ++j) {
						if (j != j_star) {
							g(J(j)) = Delta(J(j_star)) * g(J(j)) - Delta(J(j)) * f;	// no change in wdeg
						}
						else {
							polynomial_2v<T> tmp(2, 1);
							tmp(0, 0) = -alpha(i);
							tmp(1, 0) = 1;		// tmp=x - alpha(i)
							g(J(j)) = Delta(J(j_star)) * tmp * f;					// wdeg increases by 1
						}
					}
				}
			}
		}

		// change it back to (1,k-1) wdeg monomial order, revlex, we compare polynomials under this order
		ipair::weighted_i = 1;
		ipair::weighted_j = v;
		ipair::is_lex = false;
		//cout << "g" << g;

		// find the polynomial with minimun rank under (1,k-1)-revlex order
		ipair min_monomial = g(0).leading_monomial();
		int min_ind = 0;
		for (int j = 1; j <= L; ++j) {
			ipair can_monomial = g(j).leading_monomial();
			if (min_monomial < can_monomial);
			else {
				min_monomial = can_monomial;
				min_ind = j;
			}
		}
		polynomial_2v<T> Q0 = g(min_ind);

		//cout << "Q0" << Q0;
		return Q0;
	}

	/**
	 * .Kotter's interpolation algorithm, each point's multiplicity is same
	 *
	 * \param alpha: vector of points' x-value
	 * \param beta: vector of points' y-value
	 * \param multiplicity: target multiplicity at each point
	 * \param L: the max degree of y, the output list will contain at most L polynomials
	 * \param v: will find the polynomials of minimun (1,v)-revlex order
	 *
	 * \return the interpolation polynomial
	 */
	static polynomial_2v<T> solve_for_GF2e(const Matrix<T>& alpha, const Matrix<T>& beta, int multiplicity, int L, int v) {
		int n = alpha.size();
		Matrix<polynomial_2v<T>> g(1, L + 1);	// use matrix of matrix
		for (int j = 0; j <= L; ++j) {
			g(j) = polynomial_2v<T>(1, j + 1, '0');
			g(j)(0, j) = 1;
		}

		//cout << "g" << g;		// the initial polynomail, by iteration, the polymonial will satisfy all the multiplicity condition

		Matrix<T> Delta(1, L + 1, '0');	// discrepancy


		// if we enforce the conditions D_{r,s}(alpha,beta)=0 for r+s<multiplicity in 
		// an order in which (r-1,s) always precedes (r,s), the cumulative kernels 
		// will be F(x)-modules, interpolation admits of a less complex solution. (McEliece,p23)
		ipair::weighted_i = (multiplicity - 1) == 0 ? 1 : (multiplicity - 1);	// prevent multiplicity==1, which induce weighted_i be 0
		ipair::weighted_j = 1;
		ipair::is_lex = true;

		ipair top_p(multiplicity - 1, 0);
		Matrix<ipair> sorted_ipair = ipair::under_top(top_p, multiplicity);	
			// find the ipair from (0,0) to (multiplicity - 1, 0), under constrain r+s<multiplicity

		int sorted_ipair_size = sorted_ipair.size();			// sort them under (multiplicity-1,1)-lex order

		//cout << "sorted_ipair_size=" << sorted_ipair_size << endl;
		//cout << "sorted_ipair" << sorted_ipair;

		// change it back to (1,k-1) wdeg monomial order, revlex, we compare polynomials under this order
		ipair::weighted_i = 1;
		ipair::weighted_j = v;
		ipair::is_lex = false;

		for (int i = 0; i < n; ++i) {
			//cout << "g" << g;
			for (int p = 0; p < sorted_ipair_size; ++p) {
				int r = sorted_ipair(p).i;
				int s = sorted_ipair(p).j;
				Matrix<int> J(1, 0, 'v');
				for (int j = 0; j <= L; ++j) {
					Delta(j) = g(j).Hasse_D_for_GF2e(r, s, alpha(i), beta(i));	// jth discrepancy

					//Delta(j) = g(j).Hasse_D(r, s, alpha(i), beta(i));	// jth discrepancy
					if (Delta(j) != 0) {
						J.push_back(j);		// store the index where discrepency is not 0, as J
					}
				}
				

				if (J.size() != 0) {
					int Js = J.size();
					ipair min_monomial = g(J(0)).leading_monomial();
					int j_star = 0;
					for (int j = 1; j < Js; ++j) {
						ipair can_monomial = g(J(j)).leading_monomial();
						if (min_monomial < can_monomial);
						else {
							min_monomial = can_monomial;
							j_star = j;		// j_star=arg min{g_j:j\in J}
						}
					}
					polynomial_2v<T> f = g(J(j_star));		// as the book says

					polynomial_2v<T> tmp(2, 1);
					tmp(0, 0) = -alpha(i);
					tmp(1, 0) = 1;		// tmp = x - alpha(i)

					// update g(j), make it pass the point (alpha(i),beta(i))
					for (int j = 0; j < Js; ++j) {				// do not need to optimize it, a waste of time
						g(J(j)) *= Delta(J(j_star));
						if (j != j_star) {
							g(J(j)) -= Delta(J(j)) * f;			// no change in wdeg
						}
						else {
							g(J(j)) = tmp * g(J(j));			// times x - alpha(i), wdeg increases by 1
						}
					}

				}
			}
		}

		//cout << "g" << g;

		// find the polynomial with minimun rank under (1,k-1)-revlex order
		ipair min_monomial = g(0).leading_monomial();
		int min_ind = 0;
		//cout << "can_monomial = " << min_monomial << ", can_monomial.weighted_degree() = " << min_monomial.weighted_degree << endl;
		for (int j = 1; j <= L; ++j) {
			ipair can_monomial = g(j).leading_monomial();
			//cout << "j = " << j << endl;
			//cout << "can_monomial = " << can_monomial << ", can_monomial.weighted_degree() = " << can_monomial.weighted_degree << endl;
			if (min_monomial < can_monomial);
			else {
				min_monomial = can_monomial;
				min_ind = j;
			}
		}
		//cout << "min_ind = " << min_ind << endl;
		polynomial_2v<T> Q0 = g(min_ind);

		return Q0;
	}

	/**
	 * .Kotter's interpolation algorithm, each point's multiplicity is same
	 *
	 * \param alpha: vector of points' x-value
	 * \param beta: vector of points' y-value
	 * \param multiplicity: target multiplicity at each point
	 * \param L: the max degree of y, the output list will contain at most L polynomials
	 * \param v: will find the polynomials of minimun (1,v)-revlex order
	 *
	 * \return the interpolation polynomial
	 */
	static polynomial_2v_map<T> solve_for_GF2e_map(const Matrix<T>& alpha, const Matrix<T>& beta, int multiplicity, int L, int v) {

		// (1,k-1) wdeg monomial order, revlex, we compare polynomials under this order
		monomial_order::weighted_i = 1;
		monomial_order::weighted_j = v;
		monomial_order::is_lex = false;

		int n = alpha.size();
		Matrix<polynomial_2v_map<T>> g(1, L + 1);	// use matrix of map
		for (int j = 0; j <= L; ++j) {
			g(j).assign(0, j, 1);		// this is only okay with 1 being the * identity element for type T
		}

		//cout << "g" << g;		// the initial polynomail, by iteration, the polymonial will satisfy all the multiplicity condition

		Matrix<T> Delta(1, L + 1, '0');	// discrepancy


		// if we enforce the conditions D_{r,s}(alpha,beta)=0 for r+s<multiplicity in 
		// an order in which (r-1,s) always precedes (r,s), the cumulative kernels 
		// will be F(x)-modules, interpolation admits of a less complex solution. (McEliece,p23)
		ipair::weighted_i = (multiplicity - 1) == 0 ? 1 : (multiplicity - 1);	// prevent multiplicity==1, which induce weighted_i be 0
		ipair::weighted_j = 1;
		ipair::is_lex = true;

		ipair top_p(multiplicity - 1, 0);
		Matrix<ipair> sorted_ipair = ipair::under_top(top_p, multiplicity);	
			// find the ipair from (0,0) to (multiplicity - 1, 0), under constrain r+s<multiplicity

		int sorted_ipair_size = sorted_ipair.size();			// sort them under (multiplicity-1,1)-lex order

		//cout << "sorted_ipair_size=" << sorted_ipair_size << endl;
		//cout << "sorted_ipair" << sorted_ipair;


		for (int i = 0; i < n; ++i) {
			//cout << "g" << g;
			for (int p = 0; p < sorted_ipair_size; ++p) {
				int r = sorted_ipair(p).i;
				int s = sorted_ipair(p).j;
				Matrix<int> J(1, 0, 'v');
				for (int j = 0; j <= L; ++j) {
					Delta(j) = g(j).Hasse_D_for_GF2e(r, s, alpha(i), beta(i));	// jth discrepancy

					//Delta(j) = g(j).Hasse_D(r, s, alpha(i), beta(i));	// jth discrepancy
					if (Delta(j) != 0) {
						J.push_back(j);		// store the index where discrepency is not 0, as J
					}
				}


				if (J.size() != 0) {
					int Js = J.size();
					pair<int,int> min_monomial = g(J(0)).leading_monomial();
					int j_star = 0;
					for (int j = 1; j < Js; ++j) {
						pair<int, int> can_monomial = g(J(j)).leading_monomial();
						if (monomial_order::lt(min_monomial, can_monomial));
						else {
							min_monomial = can_monomial;
							j_star = j;		// j_star=arg min{g_j:j\in J}
						}
					}
					//polynomial_2v<T> tmp(2, 1);
					//tmp(0, 0) = -alpha(i);
					//tmp(1, 0) = 1;		// tmp = x - alpha(i)

					// update g(j), make it pass the point (alpha(i),beta(i))
					for (int j = 0; j < Js; ++j) {
						if (j != j_star) {
							g(J(j_star)) *= Delta(J(j));		// prepare for updating g(J(j))

							g(J(j)) *= Delta(J(j_star));
							g(J(j)) -= g(J(j_star));			// no change in wdeg

							g(J(j_star)) /= Delta(J(j));		// recover g(J(star))
						}
					}
					// finally update g(J(j_star))
					g(J(j_star)) *= Delta(J(j_star));
					g(J(j_star)).times_x_plus_alpha(-alpha(i));			// wdeg increases by 1
				}
			}
		}

		//cout << "g" << g;

		// find the polynomial with minimun rank under (1,k-1)-revlex order
		pair<int, int> min_monomial = g(0).leading_monomial();
		int min_ind = 0;
		//cout << "can_monomial = " << min_monomial << ", can_monomial.weighted_degree() = " << min_monomial.weighted_degree << endl;
		for (int j = 1; j <= L; ++j) {
			pair<int, int> can_monomial = g(j).leading_monomial();
			//cout << "j = " << j << endl;
			//cout << "can_monomial = " << can_monomial << ", can_monomial.weighted_degree() = " << can_monomial.weighted_degree << endl;
			if (monomial_order::lt(min_monomial, can_monomial));
			else {
				min_monomial = can_monomial;
				min_ind = j;
			}
		}
		//cout << "min_ind = " << min_ind << endl;

		//polynomial_2v<T> Q0 = g(min_ind);
		return g(min_ind);

		//cout << "Q0" << Q0;
		//return Q0;
	}
	
	/**
	 * .Kotter's interpolation algorithm, each point's multiplicity is same
	 *
	 * \param alpha: vector of points' x-value
	 * \param beta: vector of points' y-value
	 * \param multiplicity: target multiplicity at each point
	 * \param L: the max degree of y, the output list will contain at most L polynomials
	 * \param v: will find the polynomials of minimun (1,v)-revlex order
	 *
	 * \return the interpolation polynomial
	 */
	static polynomial_2v_umap<T> solve_for_GF2e_umap(const Matrix<T>& alpha, const Matrix<T>& beta, int multiplicity, int L, int v) {

		// (1,k-1) wdeg monomial order, revlex, we compare polynomials under this order
		monomial_order::weighted_i = 1;
		monomial_order::weighted_j = v;
		monomial_order::is_lex = false;

		int n = alpha.size();
		vector<polynomial_2v_umap<T>> g(L + 1);	// use matrix of unordered map
		//Matrix<polynomial_2v<T>> g2(1, L + 1);
		for (int j = 0; j <= L; ++j) {
			g[j].set_max_y_pow(j);
			g[j].assign(0, j, 1);		// this is only okay with 1 being the * identity element for type T
			
			//g2[j] = polynomial_2v<T>(1, j + 1, '0');
			//g2[j](0, j) = 1;
		}

		//cout << "g" << g;		// the initial polynomail, by iteration, the polymonial will satisfy all the multiplicity condition

		Matrix<T> Delta(1, L + 1, '0');	// discrepancy


		// if we enforce the conditions D_{r,s}(alpha,beta)=0 for r+s<multiplicity in 
		// an order in which (r-1,s) always precedes (r,s), the cumulative kernels 
		// will be F(x)-modules, interpolation admits of a less complex solution. (McEliece,p23)
		ipair::weighted_i = (multiplicity - 1) == 0 ? 1 : (multiplicity - 1);	// prevent multiplicity==1, which induce weighted_i be 0
		ipair::weighted_j = 1;
		ipair::is_lex = true;

		ipair top_p(multiplicity - 1, 0);
		Matrix<ipair> sorted_ipair = ipair::under_top(top_p, multiplicity);	
		// find the ipair from (0,0) to (multiplicity - 1, 0), under constrain r+s<multiplicity
		int sorted_ipair_size = sorted_ipair.size();			// sort them under (multiplicity-1,1)-lex order

		//cout << "sorted_ipair_size=" << sorted_ipair_size << endl;
		//cout << "sorted_ipair" << sorted_ipair;

		// back to compare monomial
		/*ipair::weighted_i = 1;
		ipair::weighted_j = v;
		ipair::is_lex = false;*/

		for (int i = 0; i < n; ++i) {
			/*for (int b = 0; b < L + 1; ++b)
				cout << "g[" << b << "]" << g[b];*/

			for (int p = 0; p < sorted_ipair_size; ++p) {
				int r = sorted_ipair(p).i;
				int s = sorted_ipair(p).j;
				Matrix<int> J(1, 0, 'v');
				for (int j = 0; j <= L; ++j) {
					Delta(j) = g[j].Hasse_D_for_GF2e(r, s, alpha(i), beta(i));	// jth discrepancy
					
					//T tmp = g2[j].Hasse_D_for_GF2e(r, s, alpha(i), beta(i));	// jth discrepancy
					/*if (Delta(j) != tmp) {
						cout << "error" << endl;
					}*/

					//Delta(j) = g(j).Hasse_D(r, s, alpha(i), beta(i));	// jth discrepancy
					if (Delta(j) != 0) {
						J.push_back(j);		// store the index where discrepency is not 0, as J
					}
				}


				if (J.size() != 0) {
					int Js = J.size();
					pair<int, int> min_monomial = g[J(0)].leading_monomial();
					//pair<int, int> min_monomial2 = make_pair(g2[J(0)].leading_monomial().i, g2[J(0)].leading_monomial().j);
					/*if (min_monomial != min_monomial2) {
						cout << "error" << endl;
					}*/

					int j_star = 0;
					for (int j = 1; j < Js; ++j) {
						pair<int, int> can_monomial = g[J(j)].leading_monomial();
						/*pair<int, int> min_monomial2 = make_pair(g2[J(j)].leading_monomial().i, g2[J(j)].leading_monomial().j);

						if (can_monomial != min_monomial2) {
							cout << "error" << endl;
						}*/

						if (monomial_order::lt(min_monomial, can_monomial));
						else {
							min_monomial = can_monomial;
							j_star = j;		// j_star=arg min{g_j:j\in J}
						}
					}
					//polynomial_2v<T> tmp(2, 1);
					//tmp(0, 0) = -alpha(i);
					//tmp(1, 0) = 1;		// tmp = x - alpha(i)

					// update g(j), make it pass the point (alpha(i),beta(i))
					for (int j = 0; j < Js; ++j) {
						if (j != j_star) {						// this is all good
							g[J(j_star)] *= Delta(J(j));		// prepare for updating g(J(j))
							
							//g2[J(j_star)] *= Delta(J(j));		// prepare for updating g(J(j))
							//if (!g[J(j_star)].contain_Matrix(g2[J(j_star)].get_coeff())) {
							//	cout << "error" << endl;
							//}

							g[J(j)] *= Delta(J(j_star));
							/*g2[J(j)] *= Delta(J(j_star));
							if (!g[J(j_star)].contain_Matrix(g2[J(j_star)].get_coeff())) {
								cout << "error" << endl;
							}*/

							g[J(j)] -= g[J(j_star)];			// no change in wdeg
							
							//g2[J(j)] -= g2[J(j_star)];			// no change in wdeg
							//if (!g[J(j_star)].contain_Matrix(g2[J(j_star)].get_coeff())) {
							//	cout << "error" << endl;
							//}

							g[J(j_star)] /= Delta(J(j));		// recover g(J(star))
							
							//g2[J(j_star)] /= Delta(J(j));		// recover g(J(star))
							//if (!g[J(j_star)].contain_Matrix(g2[J(j_star)].get_coeff())) {
							//	cout << "error" << endl;
							//}
						}
					}
					// finally update g(J(j_star))
					g[J(j_star)] *= Delta(J(j_star)); 
					
					//g2[J(j_star)] *= Delta(J(j_star));		// recover g(J(star))
					//if (!g[J(j_star)].contain_Matrix(g2[J(j_star)].get_coeff())) {
					//	cout << "error" << endl;
					//}

					g[J(j_star)].times_x_plus_alpha(-alpha(i));			// wdeg increases by 1
					
					//g2[J(j_star)].times_x_plus_alpha(-alpha(i));			// wdeg increases by 1
					//if (!g[J(j_star)].contain_Matrix(g2[J(j_star)].get_coeff())) {
					//	cout << "g[J(j_star)]" << g[J(j_star)] << endl;
					//	cout << "g2[J(j_star)]" << g2[J(j_star)] << endl;
					//	cout << "error" << endl;
					//}
					
				}
			}
		}

		//cout << "g" << g;

		// find the polynomial with minimun rank under (1,k-1)-revlex order
		pair<int, int> min_monomial = g[0].leading_monomial();
		int min_ind = 0;
		//cout << "can_monomial = " << min_monomial << ", can_monomial.weighted_degree() = " << min_monomial.weighted_degree << endl;
		for (int j = 1; j <= L; ++j) {
			pair<int, int> can_monomial = g[j].leading_monomial();
			//cout << "j = " << j << endl;
			//cout << "can_monomial = " << can_monomial << ", can_monomial.weighted_degree() = " << can_monomial.weighted_degree << endl;
			if (monomial_order::lt(min_monomial, can_monomial));
			else {
				min_monomial = can_monomial;
				min_ind = j;
			}
		}
		//cout << "min_ind = " << min_ind << endl;

		//polynomial_2v<T> Q0 = g(min_ind);
		return g[min_ind];

		//cout << "Q0" << Q0;
		//return Q0;
	}

};

template<class T> class RR_root_finding {
public:

	/* input */
	polynomial_2v<T> Q;
	int D;
	int field_cardinality;		// make sure the field region be {0,1,...,field_cardinality-1}

	/* variables during process */
	Matrix<polynomial_2v<T>> Q_j;
	Matrix<T> u;

	/* output */
	Matrix<Matrix<T>> y_roots;

	RR_root_finding():D(0), field_cardinality(0) {}

	void set_input(const polynomial_2v<T>& _Q, int _D, int _field_cardinality) {
		Q = _Q;
		D = _D;
		field_cardinality = _field_cardinality;
		//cout << "field_cardinality="<<field_cardinality << endl;
	}

	/**
	 * .compute Q_next(x,y)=<<Q(x,xy+a)>>
	 */
	static polynomial_2v<T> compute_Q_next(const polynomial_2v<T>& _Q,T a){
		int r = _Q.row();
		int c = _Q.col();

		int Q_next_r = r + c - 1;
		int Q_next_c = c;

		polynomial_2v<T> Q_next(Q_next_r, Q_next_c, '0');

		// Q_next(x,y)=Q(x,xy+a)
		for (int j = 0; j < c; ++j) {
			for (int i = j; i - j < r; ++i) {
				Q_next(i, j) = _Q.Hasse_D(i - j, j, 0, a);
			}
		}
		// Q_next(x,y)=<<Q_next(x,y)>>
		Q_next.simplify(); 
		Q_next_r = Q_next.row();
		Q_next_c = Q_next.col();

		int Qs = Q_next_r * Q_next_c;
		int i = 0;
		for (; i < Qs; ++i) {
			if (Q_next(i) == 0);
			else
				break;
		}
		int start_row = i / Q_next_c;
		if (start_row == 0)
			return Q_next;
		else 
			return Q_next.get_part(start_row, 0, Q_next_r - 1, Q_next_c - 1);
	}

	/**
	 * .RR root finding, find f(x) such that (y-f(x))|Q(x,y), and deg f(x)<=D
	 *
	 * \param Q: 2 polynomial of 2 variables
	 * \param D: max degree of f(x)
	 */
	void solve() {

		// remove x of Q if x|Q(x,y)
		int Qr = Q.row();
		int Qc = Q.col();
		int sss = Qr * Qc;
		int i = 0;
		for (; i < sss; ++i) {
			if (Q(i) == 0);
			else {
				break;
			}
		}
		int start_row = i / Qc;

		Q_j.push_back(polynomial_2v<T>(Q.get_part(start_row, 0, Qr - 1, Qc - 1)));
		RR_DFS();
		Q_j.pop_back();
	}
	void RR_DFS() {
		/* deep first search beginning at u */
		if (Q_j.back().is_0_when_y_eq_0()) {
			// add u into the result list
			y_roots.push_back(u);
			//cout << "y_roots" << y_roots;
		}
		else if (u.size() < D) {
			Matrix<T> roots = Q_j.back().root_when_x_eq_0(field_cardinality);
			//cout << "roots" << roots;
			int roots_num = roots.size();
			for (int i = 0; i < roots_num; ++i) {
				u.push_back(roots(i));
				Q_j.push_back(compute_Q_next(Q_j.back(), roots(i)));
				//cout << "Q_j.back()" << Q_j.back();
				RR_DFS();		// DFS
				u.pop_back();
				Q_j.pop_back();
			}
		}
	}

	/**
	 * .compute Q_next(x,y)=<<Q(x,xy+a)>>
	 */
	static polynomial_2v<T> compute_Q_next_for_GF2e(const polynomial_2v<T>& _Q, T a) {
		int r = _Q.row();
		int c = _Q.col();

		int Q_next_r = r + c - 1;
		int Q_next_c = c;

		polynomial_2v<T> Q_next(Q_next_r, Q_next_c, '0');

		// Q_next(x,y)=Q(x,xy+a)
		for (int j = 0; j < c; ++j) {
			for (int i = j; i - j < r; ++i) {
				Q_next(i, j) = _Q.Hasse_D_for_GF2e(i - j, j, 0, a);
			}
		}
		// Q_next(x,y)=<<Q_next(x,y)>>
		Q_next.simplify();
		Q_next_r = Q_next.row();
		Q_next_c = Q_next.col();

		int Qs = Q_next_r * Q_next_c;
		int i = 0;
		for (; i < Qs; ++i) {
			if (Q_next(i) == 0);
			else
				break;
		}
		int start_row = i / Q_next_c;
		if (start_row == 0)
			return Q_next;
		else
			return Q_next.get_part(start_row, 0, Q_next_r - 1, Q_next_c - 1);
	}

	/**
	 * .RR root finding, find f(x) such that (y-f(x))|Q(x,y), and deg f(x)<=D
	 *
	 * \param Q: 2 polynomial of 2 variables
	 * \param D: max degree of f(x)
	 */
	void solve_for_GF2e() {

		// remove x of Q if x|Q(x,y)
		int Qr = Q.row();
		int Qc = Q.col();
		int sss = Qr * Qc;
		int i = 0;
		for (; i < sss; ++i) {
			if (Q(i) == 0);
			else {
				break;
			}
		}
		int start_row = i / Qc;

		Q_j.push_back(polynomial_2v<T>(Q.get_part(start_row, 0, Qr - 1, Qc - 1)));
		RR_DFS_for_GF2e();
		Q_j.pop_back();
	}
	void RR_DFS_for_GF2e() {
		/* deep first search beginning at u */
		if (Q_j.back().is_0_when_y_eq_0()) {
			// add u into the result list
			y_roots.push_back(u);
			//cout << "y_roots" << y_roots;
		}
		else if (u.size() < D) {
			Matrix<T> roots = Q_j.back().root_when_x_eq_0(field_cardinality);
			//cout << "roots" << roots;
			int roots_num = roots.size();
			for (int i = 0; i < roots_num; ++i) {
				u.push_back(roots(i));
				Q_j.push_back(compute_Q_next_for_GF2e(Q_j.back(), roots(i)));
				//cout << "Q_j.back()" << Q_j.back();
				RR_DFS_for_GF2e();		// DFS
				u.pop_back();
				Q_j.pop_back();
			}
		}
	}
};

class Lagrange_interpolation{
public:
	/**
	 * .the returned polynomial will have 0 at point alpha expect for the point alpha(alpha_ind), at which the evaluation is 1
	 */
	template<class T> static polynomial<T> generate(const Matrix<T>& alpha, int alpha_ind) {
		T beta = alpha(alpha_ind);
		int alpha_size = alpha.size();
		
		polynomial<T> result(Matrix<T>(1, alpha_size, '0'));	// allocate enough size, we will not need to declear new array
		result(0) = 1;
		result.simplify();								// result = 1
		polynomial<T> term(2, '0');
		term(1) = 1;
		for (int i = 0; i < alpha_size; ++i) {
			if (i != alpha_ind) {
				term(0) = -alpha(i);					// term = -alpha(i) + x
				result *= term / (alpha(i) - beta);
			}
		}
		//cout << "result" << result;
		return result;
	}
};
