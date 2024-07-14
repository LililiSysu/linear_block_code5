/*****************************************************************//**
 * \file   eBCH.h
 * \brief  eBCH coding and decoding
 * 
 * \author lilili
 * \date   January 2023
 *********************************************************************/

#pragma once

#include"../GF/GF2e.h"		//already include "GF2.h"
#include"../GF/polynomial.h"
#include"../reprocess/LRP.h"
#include"../reprocess/MRIP.h"

template<int m, int t> class eBCH {		// noly defined for primitive BCH code
private:
	// constrain: m>3 and t<2^{m-1}
	// n=2^m-1=(1>>m)-1
	int r;		// number of redundancy, not rate
	int k;		// number of information bits
	int n;		// length of code
	int d;		// minimum Hamming distance
	polynomial<GF2> gX2;	// generator polynomial, length of r+1, that is max order r
	polynomial<GF2> pX2;	// parity polynomial, length of k+1, that is max order k
	Matrix<GF2> generator_M;// generator matrix, size of k*n
	Matrix<GF2> parity_M;	// parity matrix, size of r*n

public:
	eBCH() {
		Matrix<GF2e<m>> tmp(1, 2, { 1,0 });
		n = (1 << m) - 1;
		Matrix<bool> alpha_used(1, n + 1, '0');
		polynomial<GF2e<m>> gX(tmp);		// gX=1

		//cout << "gX=" << gX;
		tmp(1) = 1;				// tmp={1,1}
		for (int i = 1; i <= 2 * t; ++i) {
			GF2e<m> alpha;
			alpha.set_by_alpha_power(i);
			while (!alpha_used((int)alpha)) {
				//cout << "alpha=" << alpha << endl;
				tmp(0) = alpha;		//tmp={alpha,1}

				gX = gX * tmp;					//gX=gX*(alpha+X)

				//cout << "gX=" << gX;

				alpha_used((int)alpha) = true;
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
			gX2_coeff(i) = (int)(gX(i));
		}
		gX2 = gX2_coeff;
		//cout << "gX2=" << gX2;
		k = n - r;
		//cout << "n=" << n << "\tk=" << k << "\tr=" << r << endl;
		d = 2 * t + 1;
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

		// extended BCH code generator matrix and parity check matrix
		Matrix<GF2> k1_colv = Matrix<GF2>(k, 1, '1');

		// be a row vector
		Matrix<GF2> q = (generator_M.get_part(0, 0, k - 1, k - 1).inv() * k1_colv).Transpose().combine_right(Matrix<GF2>(1, n - k, '0'));
			
		//cout << "q" << q;
		generator_M = generator_M.combine_right(k1_colv);
		parity_M = parity_M.combine(Matrix<GF2>(n - k, 1, '0'), q, GF2(1));
		n++;
		r++;
		d++;

		//cout << "generator_M" << generator_M;
		//cout << "parity_M" << parity_M;
		//cout << "zeros" << generator_M * parity_M.Transpose();
	}

	~eBCH() {}

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

	void print_info() const {
		cout << "eBCH (" << n << "," << k << "," << d << ") code, rate = " << k / (double)n << endl;
	}

	// check if in v space
	bool is_in_v_space(const Matrix<GF2>& v) {
		polynomial<GF2> result = v.get_part(0, 0, 0, n - 2) % gX2;
		if (result.size() == 1 && result(0) == 0 && \
			parity_M.get_row(n - k - 1).dot_product(v) == 0) {
			return true;
		}
		else {
			return false;
		}
	}

	// also u2c
	Matrix<GF2> encode(const Matrix<GF2>& u) {		// encode like cyclic code
		//polynomial<GF2> pu(u);			// store 'u' as polynomial
		//cout << "u=" << u;

		/* encode to codeword by generator polynomial */
		polynomial<GF2> pc = u * gX2;
		pc.format_len(n);
		pc(n - 1) = u.dot_product(Matrix<GF2>(k, 1, '1'));
		return pc.get_coeff();
		/*for (int i = 0; i < n; ++i) {
			c(i)= pu.evaluate(i + 1);
		}

		return c;*/
	}
	
	// decode to u
	Matrix<GF2> v2u(const Matrix<GF2>& v) {
		polynomial<GF2> result = v.get_part(0, 0, 0, n - 2) / gX2;
		result.format_len(k);
		return result.get_coeff();
	}

	// the false decode
	Matrix<GF2> decode_v(const Matrix<GF2>& v) {
		return v;
	}
	Matrix<GF2> decode(const Matrix<GF2>& v) {
		return v2u(v);
	}

	friend class GMD;
	friend class Chase;
	friend class OSD;
};
