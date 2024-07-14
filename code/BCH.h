/*****************************************************************//**
 * \file   BCH.h
 * \brief  BCH coding and decoding
 * \param	'u' --[channel_coding]--> 'v' --[BPSK_modulation]--> 'c' --[AWGN_channel]--> 'r'
 *		   --[BPSK_demodulation]--> 'hdr'('z') --[decode_v]--> 'dv'('v') --[v2u]--> 'du'
 * 
 * \author lilili
 * \date   October 2022
 *********************************************************************/

#pragma once

#include"../GF/GF2e.h"		//already include "GF2.h"
#include"../GF/polynomial.h"
#include"../reprocess/LRP.h"
#include"../reprocess/MRIP.h"

/**
 * assume beta is a primitive nth root of unity in GF2e<m>
 * 
 * .BCH code C(q,n,d,h) has generator polynomial
 *		lcm{M_{beta^h}(x),M_{beta^{h+1}}(x), . . . ,M_{beta^{h+d-2}}(x)},
 * 
 *	where M_{beta^j}(x) denotes the minimal polynomial of alpha^j over
 *	GF(q), and lcm denotes the least common multiple of a set
 *	of polynomials. For the different choices of beta, the BCH codes
 *	C(q,n,d,h) are equal up to a fixed permutation on the codeword coordinates.
 *	
 * but if beta is not a proimtive nth root of unity in GF2e<m>, BCH code becomes non-primitive one
 */

// the narrow sense primative BCH code, represented as C(q,n,d,h) => C(2,2^m-1,2t+1,1)
template<int m, int t> class BCH {
private:
	// constrain: m>3 and t<2^{m-1}
	// n=2^m-1=(1>>m)-1
	int n;		// length of code
	int k;		// number of information bits
	int d;		// minimum Hamming distance
	int r;		// number of redundancy, not rate
	int d_dual;	// minimum Hamming distance of dual code
	polynomial<GF2> gX2;	// generator polynomial, length of r+1, that is max order r
	polynomial<GF2> pX2;	// parity polynomial, length of k+1, that is max order k
	Matrix<GF2> generator_M;// generator matrix, size of k*n
	Matrix<GF2> parity_M;	// parity matrix, size of r*n

	/**
	 * . for detail, SEE "Shift-Register Synthesis and BCH Decoding"
	 *
	 * \param Sn: the target to generate
	 * \return the coefficients of LFSR that generate Sn
	 */
	polynomial<GF2e<m>> LFSR_Synthesis_Algorithm_BCH(const polynomial<GF2e<m>>& Sn) const {

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
			N += 2;		// acceletate for BCH code for d=0 if N is odd
			x++;
		}
		C.simplify();

		return C;
	}

public:
	BCH() {
		Matrix<GF2e<m>> tmp(1, 2, { 1,0 });
		n = (1 << m) - 1;
		Matrix<bool> alpha_used(1, n + 1, '0');
		polynomial<GF2e<m>> gX(tmp);		// gX=1
		
		//cout << "gX=" << gX;
		tmp(1) =1;				// tmp={1,1}
		for (int i = 1; i <= 2 * t; ++i) {
			GF2e<m> alpha;
			alpha.set_by_alpha_power(i);
			while (!alpha_used((int)alpha)) {
				//cout << "alpha=" << alpha << endl;
				tmp(0)= alpha;		//tmp={alpha,1}

				gX = gX * tmp;					//gX=gX*(alpha+X)
				
				//cout << "gX=" << gX;

				alpha_used((int)alpha)= true;
				alpha = alpha * alpha;
			}
			//tmp(0)= -tmp(0);		// set tmp={-alpha,1}, for GF2e, this is not need
			//cout << "###" << endl;
			//cout << "gX=" << gX;
		}
		//cout << "-----" << endl;
		//cout << "gX=" << gX;

		r = gX.size() - 1;		// number of redundancy, gX coeff len is 1+r
		Matrix<GF2> gX2_coeff(1, r + 1);
		for (int i = 0; i < r + 1; ++i) {
			gX2_coeff(i)= (int)(gX(i));
		}
		gX2 = gX2_coeff;
		//cout << "gX2=" << gX2;
		k = n - r;
		//cout << "n=" << n << "\tk=" << k << "\tr=" << r << endl;
		d = 2 * t + 1;			// the designed distance, lower bound of the true distance, be careful
		polynomial<GF2> xup(n + 1, '0');
		xup(0) = 1;
		xup(n) = 1;

		pX2 = xup / gX2;
		pX2.rev();
		//Matrix<GF2> revp = pX2.get_coeff();
		//revp.rev();
		//pX2 = revp;		// reverse

		int glen = gX2.size();
		generator_M = Matrix<GF2>(k, n, '0');
		for (int i = 0; i < k; ++i) {
			for (int j = 0; j < glen; ++j)
				generator_M(i * n + i + j) = gX2(j);
		}

		int plen = pX2.size();
		parity_M = Matrix<GF2>(r, n, '0');
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < plen; ++j)
				parity_M(i * n + i + j) = pX2(j);
		}

		//cout << "parity_M" << parity_M;
		//cout << "generator_M.multiply_transpose_of(parity_M).isZero() = " << generator_M.multiply_transpose_of(parity_M).isZero() << endl;
	}

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

	Matrix<GF2> get_generator_matrix() const {
		// size of k*n matrix
		return generator_M;
	}

	Matrix<GF2> get_parity_matrix() const {
		// size of r*n matrix
		return parity_M;
	}

	polynomial<GF2> get_gX2() const {
		return gX2;
	}

	polynomial<GF2> get_pX2() const {
		return pX2;
	}

	void print_info() const {
		cout << "BCH (" << n << "," << k << "," << d << ") code, rate = " << k / (double)n << endl;
	}

	// check if in v space
	bool is_in_v_space(const Matrix<GF2>& v) const {


#ifdef RUN_MSG
#ifdef count_operation_number
		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		polynomial<GF2> result = v % gX2;

#ifdef RUN_MSG
#ifdef count_operation_number
		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(is_in_v_space) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		if (result.size() == 1 && result(0) == 0) {
			return true;
		}
		else {
			return false;
		}
	}

	// also u2c
	Matrix<GF2> encode(const Matrix<GF2>& u) const {		// encode like cyclic code
		//polynomial<GF2> pu(u);			// store 'u' as polynomial
		//cout << "u=" << u;

		/* encode to codeword by generator polynomial */
		polynomial<GF2> pc = u * gX2;
		pc.format_len(n);
		return pc.get_coeff();
		/*for (int i = 0; i < n; ++i) {
			c(i)= pu.evaluate(i + 1);
		}

		return c;*/
	}

	// BM decoding as a main decode method, decode to u
	Matrix<GF2> decode(const Matrix<GF2>& hdr) const {
		//cout << "decode_v(hdr)" << decode_v(hdr);
		return v2u(decode_v(hdr));
	}

	// decode to v
	Matrix<GF2> decode_v(const Matrix<GF2>& hdr) const {		// decode result may not in v space

#ifdef RUN_MSG
#ifdef count_operation_number
		int GF2e_ope_num_before = GF2e_auxiliary_storage::operation_number;
		int GF2e_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		polynomial<GF2e<m>> pr(n, '0');
		for (int i = 0; i < n; ++i) {
			pr(i) = (int)hdr(i);	// store r as polynomial
		}

		/* compute syndrome */
		bool flag_no_error = 1;
		polynomial<GF2e<m>> pS(2 * t, '0');
		int tt_1 = 2 * t - 1;
		GF2e<m> tmp_for_new_GF2e;
		for (int i = 0; i < tt_1; i += 2) {
			tmp_for_new_GF2e.set_by_alpha_power(i + 1);
			pS(i) = pr.evaluate(tmp_for_new_GF2e);		// pr(alpha),pr(alpha^1), ... ,pr(alpha^{2t})
			pS(i + 1) = pS(i / 2) * pS(i / 2);	// for S_{2i}=S_i^2 in BCH code
			if (pS(i) != 0)flag_no_error = 0;	// has some error
		}
		if (flag_no_error) {
			// no error
			return hdr;
		}

		/* call 'LFSR_Synthesis_Algorithm' to generate location polynomial Lambda(X)*/
		//polynomial<GF2e<m>> pLambda = LFSR_Synthesis_Algorithm_BCH(pS);
		polynomial<GF2e<m>> pLambda = LFSR_Synthesis_Algorithm_BCH(pS);

		/* find root's reciprocal value, then it is error location*/
		Matrix<GF2e<m>> error_location = pLambda.find_roots(1 << m);

		/* correct error in r */
		Matrix<GF2> v(hdr);
		int t_real = error_location.size();

		tmp_for_new_GF2e = 1;
		for (int i = 0; i < t_real; ++i) {
			int error_pos = GF2e_auxiliary_storage::alpha_table[(int)(tmp_for_new_GF2e / error_location(i))];
			//int error_pos = (int)(1 / error_location(i)) - 1;
			v(error_pos) = v(error_pos) == 0/* flip this bit */;
		}

#ifdef RUN_MSG
#ifdef count_operation_number
		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << "(Algebra) GF2e_ope_num = " << GF2e_ope_num_after - GF2e_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return v;
	}

	// decode to u
	Matrix<GF2> v2u(const Matrix<GF2>& v) const {
#ifdef RUN_MSG
#ifdef count_operation_number
		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG


		polynomial<GF2> result = v / gX2;
		result.format_len(k);

#ifdef RUN_MSG
#ifdef count_operation_number
		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(v2u) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return result.get_coeff();
	}

	friend class GMD;
	friend class Chase;
	friend class OSD;

	// function with dual code
	Matrix<GF2> encode_dual(const Matrix<GF2>& u_perp) const {

		/* encode to codeword by parity polynomial */
		polynomial<GF2> pc = u_perp * pX2;
		pc.format_len(n);
		return pc.get_coeff();
	}

	bool is_in_dual_space(const Matrix<GF2>& v_perp) const {

		polynomial<GF2> result = v_perp % pX2;
		return result.size() == 1 && result(0) == 0;
	}

	Matrix<GF2> v2u_dual(const Matrix<GF2>& v_perp) const {

		polynomial<GF2> result = v_perp / pX2;
		result.format_len(r);
		return result.get_coeff();
	}

};

// the non-narrow sense primative BCH code, represented as C(q,n,d,h) => C(2,2^m-1,2t+2,0), same as shortening one bit of eBCH code
template<int m, int t> class nnsBCH {
private:
	// constrain: m>3 and t<2^{m-1}
	// n=2^m-1=(1>>m)-1
	int n;		// length of code
	int k;		// number of information bits
	int d;		// minimum Hamming distance
	int r;		// number of redundancy, not rate
	int d_dual;	// minimum Hamming distance of dual code
	polynomial<GF2> gX2;	// generator polynomial, length of r+1, that is max order r
	polynomial<GF2> pX2;	// parity polynomial, length of k+1, that is max order k
	Matrix<GF2> generator_M;// generator matrix, size of k*n
	Matrix<GF2> parity_M;	// parity matrix, size of r*n

public:
	nnsBCH() {
		Matrix<GF2e<m>> tmp(1, 2, { 1,0 });
		n = (1 << m) - 1;
		Matrix<bool> alpha_used(1, n + 1, '0');
		polynomial<GF2e<m>> gX(tmp);		// gX=1

		//cout << "gX=" << gX;
		tmp(1) = 1;				// tmp={1,1}
		for (int i = 0; i <= 2 * t; ++i) {		// ******** i start from 0 ********
			GF2e<m> alpha;
			alpha.set_by_alpha_power(i);
			//cout << "alpha = " << alpha << endl;
			while (!alpha_used((int)alpha)) {
				//cout << "alpha=" << alpha << endl;
				tmp(0) = alpha;		//tmp={alpha,1}

				gX = gX * tmp;					//gX=gX*(alpha+X)

				//cout << "gX=" << gX;

				alpha_used((int)alpha) = true;
				alpha = alpha * alpha;			// traverse the 2-cyclotomic coset
			}
			//tmp(0)= -tmp(0);		// set tmp={-alpha,1}, for GF2e, this is not need
			//cout << "###" << endl;
			//cout << "gX=" << gX;
		}
		//cout << "-----" << endl;
		//cout << "gX=" << gX;

		r = gX.size() - 1;		// number of redundancy, gX coeff len is 1+r
		Matrix<GF2> gX2_coeff(1, r + 1);
		for (int i = 0; i < r + 1; ++i) {
			gX2_coeff(i) = (int)(gX(i));
		}
		gX2 = gX2_coeff;
		//cout << "gX2=" << gX2;
		k = n - r;
		//cout << "n=" << n << "\tk=" << k << "\tr=" << r << endl;
		d = 2 * t + 2;			// the designed distance, lower bound of the true distance, be careful
		polynomial<GF2> xup(n + 1, '0');
		xup(0) = 1;
		xup(n) = 1;

		pX2 = xup / gX2;
		pX2.rev();
		//Matrix<GF2> revp = pX2.get_coeff();
		//revp.rev();
		//pX2 = revp;		// reverse

		int glen = gX2.size();
		generator_M = Matrix<GF2>(k, n, '0');
		for (int i = 0; i < k; ++i) {
			for (int j = 0; j < glen; ++j)
				generator_M(i * n + i + j) = gX2(j);
		}

		int plen = pX2.size();
		parity_M = Matrix<GF2>(r, n, '0');
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < plen; ++j)
				parity_M(i * n + i + j) = pX2(j);
		}

		//cout << "parity_M" << parity_M;
		//cout << "generator_M.multiply_transpose_of(parity_M).isZero() = " << generator_M.multiply_transpose_of(parity_M).isZero() << endl;
	}

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

	Matrix<GF2> get_generator_matrix() const {
		// size of k*n matrix
		return generator_M;
	}

	Matrix<GF2> get_parity_matrix() const {
		// size of r*n matrix
		return parity_M;
	}

	polynomial<GF2> get_gX2() const {
		return gX2;
	}

	polynomial<GF2> get_pX2() const {
		return pX2;
	}

	void print_info() const {
		cout << "BCH (" << n << "," << k << "," << d << ") code, rate = " << k / (double)n << endl;
	}

	// check if in v space
	bool is_in_v_space(const Matrix<GF2>& v) const {


#ifdef RUN_MSG
#ifdef count_operation_number
		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		polynomial<GF2> result = v % gX2;

#ifdef RUN_MSG
#ifdef count_operation_number
		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(is_in_v_space) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		if (result.size() == 1 && result(0) == 0) {
			return true;
		}
		else {
			return false;
		}
	}

	// also u2c
	Matrix<GF2> encode(const Matrix<GF2>& u) const {		// encode like cyclic code
		//polynomial<GF2> pu(u);			// store 'u' as polynomial
		//cout << "u=" << u;

		/* encode to codeword by generator polynomial */
		polynomial<GF2> pc = u * gX2;
		pc.format_len(n);
		return pc.get_coeff();
		/*for (int i = 0; i < n; ++i) {
			c(i)= pu.evaluate(i + 1);
		}

		return c;*/
	}

	// decode to u
	Matrix<GF2> v2u(const Matrix<GF2>& v) const {
#ifdef RUN_MSG
#ifdef count_operation_number
		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG


		polynomial<GF2> result = v / gX2;
		result.format_len(k);

#ifdef RUN_MSG
#ifdef count_operation_number
		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(v2u) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return result.get_coeff();
	}

	friend class GMD;
	friend class Chase;
	friend class OSD;

	// function with dual code
	Matrix<GF2> encode_dual(const Matrix<GF2>& u_perp) const {

		/* encode to codeword by parity polynomial */
		polynomial<GF2> pc = u_perp * pX2;
		pc.format_len(n);
		return pc.get_coeff();
	}

	bool is_in_dual_space(const Matrix<GF2>& v_perp) const {

		polynomial<GF2> result = v_perp % pX2;
		return result.size() == 1 && result(0) == 0;
	}

	Matrix<GF2> v2u_dual(const Matrix<GF2>& v_perp) const {

		polynomial<GF2> result = v_perp / pX2;
		result.format_len(r);
		return result.get_coeff();
	}

};

// the non-primative BCH code, with beta = alpha^b, represented as C(q,n,d,h) => C(2,2^m-1,2t+1,1), when b=1, become primative bch code
template<int m, int t, int b> class npBCH {
private:
	// constrain: m>3 and t<2^{m-1}
	// n=2^m-1=(1>>m)-1
	int n;		// length of code
	int k;		// number of information bits
	int d;		// minimum Hamming distance
	int r;		// number of redundancy, not rate
	int d_dual;	// minimum Hamming distance of dual code
	polynomial<GF2> gX2;	// generator polynomial, length of r+1, that is max order r
	polynomial<GF2> pX2;	// parity polynomial, length of k+1, that is max order k
	Matrix<GF2> generator_M;// generator matrix, size of k*n
	Matrix<GF2> parity_M;	// parity matrix, size of r*n

public:
	npBCH() {
		Matrix<GF2e<m>> tmp(1, 2, { 1,0 });
		// decide the order alpha^b, which gives the size of non-primative BCH code
		GF2e<m> beta;
		beta.set_by_alpha_power(b);
		GF2e<m> alpha = beta;
		GF2e<m> identity_ele;
		identity_ele.set_by_alpha_power(0);
		for (n = 1; alpha != identity_ele; ++n) {
			alpha *= beta;
		}

		// store the alpha that in the least common multiple of a set of polinomials (x-alpha)
		Matrix<bool> alpha_used(1, n + 1, '0');
		int alpha_store_divide = ((1 << m) - 1) / n;
		polynomial<GF2e<m>> gX(tmp);		// gX=1

		//cout << "gX=" << gX;
		tmp(1) = 1;				// tmp={1,1}
		int beta_power = b;
		for (int i = 1; i <= 2 * t; ++i) {		// i start from 1
			GF2e<m> alpha;
			beta_power %= GF2e_auxiliary_storage::q_minus_1;
			alpha.set_by_alpha_power(beta_power);
			//cout << "alpha = " << alpha << endl;

			//cout << "(outer) test 0 = " << alpha.alpha_power() % alpha_store_divide << endl;		
				// must be a multiple of alpha_store_divide

			while (!alpha_used(alpha.alpha_power() / alpha_store_divide)) {
				//cout << "alpha=" << alpha << endl;
				tmp(0) = alpha;		//tmp={alpha,1}

				gX = gX * tmp;					//gX=gX*(alpha+X)

				//cout << "gX=" << gX;

				//cout << "test 0 = " << alpha.alpha_power() % alpha_store_divide << endl;

				alpha_used(alpha.alpha_power() / alpha_store_divide) = true;
				alpha = alpha * alpha;			// traverse the 2-cyclotomic coset
			}
			beta_power += b;
			//tmp(0)= -tmp(0);		// set tmp={-alpha,1}, for GF2e, this is not need
			//cout << "###" << endl;
			//cout << "gX=" << gX;
		}
		//cout << "-----" << endl;
		//cout << "gX=" << gX;

		r = gX.size() - 1;		// number of redundancy, gX coeff len is 1+r
		Matrix<GF2> gX2_coeff(1, r + 1);	// it is stupid to do that
		for (int i = 0; i < r + 1; ++i) {
			gX2_coeff(i) = (int)(gX(i));	// in raw form, elements in GF2e<m> is same as element in sub field GF2, hence using (int)
		}
		gX2 = gX2_coeff;
		//cout << "gX2=" << gX2;
		k = n - r;
		//cout << "n=" << n << "\tk=" << k << "\tr=" << r << endl;
		d = 2 * t + 1;			// the designed distance, lower bound of the true distance, be careful
		polynomial<GF2> xup(n + 1, '0');
		xup(0) = 1;
		xup(n) = 1;

		pX2 = xup / gX2;
		pX2.rev();
		//Matrix<GF2> revp = pX2.get_coeff();
		//revp.rev();
		//pX2 = revp;		// reverse

		int glen = gX2.size();
		generator_M = Matrix<GF2>(k, n, '0');
		for (int i = 0; i < k; ++i) {
			for (int j = 0; j < glen; ++j)
				generator_M(i * n + i + j) = gX2(j);
		}

		int plen = pX2.size();
		parity_M = Matrix<GF2>(r, n, '0');
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < plen; ++j)
				parity_M(i * n + i + j) = pX2(j);
		}

		//cout << "parity_M" << parity_M;
		//cout << "generator_M.multiply_transpose_of(parity_M).isZero() = " << generator_M.multiply_transpose_of(parity_M).isZero() << endl;
	}

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

	Matrix<GF2> get_generator_matrix() const {
		// size of k*n matrix
		return generator_M;
	}

	Matrix<GF2> get_parity_matrix() const {
		// size of r*n matrix
		return parity_M;
	}

	polynomial<GF2> get_gX2() const {
		return gX2;
	}

	polynomial<GF2> get_pX2() const {
		return pX2;
	}

	void print_info() const {
		cout << "BCH (" << n << "," << k << "," << d << ") code, rate = " << k / (double)n << endl;
	}

	// check if in v space
	bool is_in_v_space(const Matrix<GF2>& v) const {


#ifdef RUN_MSG
#ifdef count_operation_number
		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		polynomial<GF2> result = v % gX2;

#ifdef RUN_MSG
#ifdef count_operation_number
		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(is_in_v_space) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		if (result.size() == 1 && result(0) == 0) {
			return true;
		}
		else {
			return false;
		}
	}

	// also u2c
	Matrix<GF2> encode(const Matrix<GF2>& u) const {		// encode like cyclic code
		//polynomial<GF2> pu(u);			// store 'u' as polynomial
		//cout << "u=" << u;

		/* encode to codeword by generator polynomial */
		polynomial<GF2> pc = u * gX2;
		pc.format_len(n);
		return pc.get_coeff();
		/*for (int i = 0; i < n; ++i) {
			c(i)= pu.evaluate(i + 1);
		}

		return c;*/
	}

	// decode to u
	Matrix<GF2> v2u(const Matrix<GF2>& v) const {
#ifdef RUN_MSG
#ifdef count_operation_number
		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG


		polynomial<GF2> result = v / gX2;
		result.format_len(k);

#ifdef RUN_MSG
#ifdef count_operation_number
		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(v2u) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return result.get_coeff();
	}

	friend class GMD;
	friend class Chase;
	friend class OSD;

	// function with dual code
	Matrix<GF2> encode_dual(const Matrix<GF2>& u_perp) const {

		/* encode to codeword by parity polynomial */
		polynomial<GF2> pc = u_perp * pX2;
		pc.format_len(n);
		return pc.get_coeff();
	}

	bool is_in_dual_space(const Matrix<GF2>& v_perp) const {

		polynomial<GF2> result = v_perp % pX2;
		return result.size() == 1 && result(0) == 0;
	}

	Matrix<GF2> v2u_dual(const Matrix<GF2>& v_perp) const {

		polynomial<GF2> result = v_perp / pX2;
		result.format_len(r);
		return result.get_coeff();
	}

};