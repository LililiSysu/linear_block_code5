/*****************************************************************//**
 * \file   GS_Viterbi_newOSD.h
 * \brief  Algorithm of state-of-art research.
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"test_common.h"

class test_GS_Viterbi_newOSD {
public:

	static void plain_old_data() {

		cout << boolalpha;
		cout << "\n-------is_trivial-------\n";
		cout << "int               " << is_trivial<int>::value << endl;
		cout << "my_float          " << is_trivial<my_float>::value << endl;
		cout << "my_double         " << is_trivial<my_double>::value << endl;
		cout << "Complex           " << is_trivial<Complex>::value << endl;
		cout << "GF2               " << is_trivial<GF2>::value << endl;
		cout << "GF<7>             " << is_trivial<GF<7>>::value << endl;
		cout << "GF2e<5>           " << is_trivial<GF2e<5>>::value << endl;
		cout << "Matrix<GF2e<5>>   " << is_trivial<Matrix<GF2e<5>>>::value << endl;
		cout << "vector<int>       " << is_trivial<vector<int>>::value << endl;
		cout << "polynomial_2v<int>" << is_trivial<polynomial_2v<int>>::value << endl;

		cout << "\n-------is_trivially_destructible-------\n";
		cout << "int               " << is_trivially_destructible<int>::value << endl;
		cout << "my_float          " << is_trivially_destructible<my_float>::value << endl;
		cout << "my_double         " << is_trivially_destructible<my_double>::value << endl;
		cout << "Complex           " << is_trivially_destructible<Complex>::value << endl;
		cout << "GF2               " << is_trivially_destructible<GF2>::value << endl;
		cout << "GF<7>             " << is_trivially_destructible<GF<7>>::value << endl;
		cout << "GF2e<5>           " << is_trivially_destructible<GF2e<5>>::value << endl;
		cout << "Matrix<GF2e<5>>   " << is_trivially_destructible<Matrix<GF2e<5>>>::value << endl;
		cout << "vector<int>       " << is_trivially_destructible<vector<int>>::value << endl;
		cout << "polynomial_2v<int>" << is_trivially_destructible<polynomial_2v<int>>::value << endl;

		cout << "\n-------is_standard_layout-------\n";
		cout << "int               " << is_standard_layout<int>::value << endl;
		cout << "my_float          " << is_standard_layout<my_float>::value << endl;
		cout << "my_double         " << is_standard_layout<my_double>::value << endl;
		cout << "Complex           " << is_standard_layout<Complex>::value << endl;
		cout << "GF2               " << is_standard_layout<GF2>::value << endl;
		cout << "GF<7>             " << is_standard_layout<GF<7>>::value << endl;
		cout << "GF2e<5>           " << is_standard_layout<GF2e<5>>::value << endl;
		cout << "Matrix<GF2e<5>>   " << is_standard_layout<Matrix<GF2e<5>>>::value << endl;
		cout << "vector<int>       " << is_standard_layout<vector<int>>::value << endl;
		cout << "polynomial_2v<int>" << is_standard_layout<polynomial_2v<int>>::value << endl;
	}
	static void int_overflow() {

		int a = -2147483648;
		int b = 2147483647;

		cout << "a=" << a << endl;
		cout << "b=" << b << endl;
		cout << "a-b=" << a - b << endl;
	}
	static void rank_test() {
		Matrix<GF2> A(5, 3, {
			1,0,1,0,0,
			0,0,0,0,1,
			1,1,0,1,0
			});
		cout << "A" << A;
		cout << "A.rank()=" << A.rank() << endl;

		Matrix<GF2> B(3, 5, {
			0,0,1,0,0,
			0,0,0,0,0,
			1,1,0,1,1
			});
		cout << "B" << B;
		cout << "B.rank()=" << B.rank() << endl;

		Matrix<double> C(6, 6, {
			35,     1,     6,    26,    19,    24,
			 3,    32,     7,    21,    23,    25,
			31,     9,     2,    22,    27,    20,
			 8,    28,    33,    17,    10,    15,
			30,     5,    34,    12,    14,    16,
			 4,    36,    29,    13,    18,    11
			});
		cout << "C" << C;
		cout << "C.rank()=" << C.rank() << endl;
	}
	static void polynomial_test() {
		typedef GF2e<4> ty;
		ty::init();

		polynomial<ty> p1(Matrix<ty>(1, 6, { 4,2,1,6,7,8 }));
		cout << "p1" << p1;
		polynomial<ty> p2(Matrix<ty>(1, 6, { 0,0,1,0,1,1 }));
		cout << "p2" << p2;
		polynomial<ty> p3 = p1 + p2;
		cout << "p3" << p3;
		p3.format_len(5);
		cout << "p3" << p3;
		p3 = p3.get_derivative();
		cout << "p3" << p3;
		p1 = p1.get_derivative();
		cout << "p1" << p1;

		polynomial<ty> px(Matrix<ty>(1, 2, { 4,2 }));
		polynomial<ty> py(Matrix<ty>(1, 3, { 6,7,8 }));
		polynomial<ty> pz = px * py;
		cout << "pz" << pz;		// {9,6,7,9}

		px *= py;
		cout << "px" << px;		// {9.6.7.9}

		px = px / py;
		cout << "px" << px;		// {4,2}

		pz %= px;
		cout << "pz" << pz;		// {6,7,8}
	}
	static void GF_test() {
		GF2e<3>::init();
		GF2e<3> a = 6;
		GF2e<3> b = 2;
		int c = 3;
		GF2e<3> d = c * a;		// this means 3 'a's add together
		cout << d << endl;
	}
	static void BCH_test() {
		typedef GF2 ty;
		GF2e<3>::init();

		BCH<3, 1> bch;
		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();
		int d = bch.get_d();

		/*Matrix<ty> u(1, k, { 1,0,1,0,1,0,1 });
		cout << "u" << u;
		Matrix<ty> v = bch.encode(u);
		cout << "v" << v;
		Matrix<ty> e(1, n, { 0,0,0,0,0,0,0,1,0,0,0,0,0,0,1 });
		Matrix<ty> r = e + v;
		cout << "r" << r;

		Matrix<ty> u_hat = bch.decode(r);
		cout << "u_hat" << u_hat;*/

		cout << "bch.get_generator_matrix() = " << endl;
		bch.get_generator_matrix().print();
		cout << "bch.get_parity_matrix() = " << endl;
		bch.get_parity_matrix().print();
		cout << "bch.get_generator_matrix() * bch.get_parity_matrix() = " << endl;
		(bch.get_generator_matrix() * bch.get_parity_matrix().Transpose()).print();

		Matrix<ty> Gs = bch.get_generator_matrix();
		int cp;
		Gs.GJE_left_identity_4_GF2(cp);
		cout << "Gs = " << endl;
		Gs.print();
		cout << "Gs * bch.get_parity_matrix() = " << endl;
		(Gs * bch.get_parity_matrix().Transpose()).print();
	}
	static void BCH_info() {
		const int m = 10;	// GF field of 2^m
		const int t = 171;	// error that can correct
		GF2e<m>::init();
		eBCH<m, t> ebch;
		ebch.print_info();
	}
	static void BCH_info2() {
		const int m = 4;	// GF field of 2^m
		const int t = 2;	// error that can correct
		GF2e<m>::init();

		BCH<m, t> bch;
		bch.print_info();
		cout << "bch.get_gX2() = " << bch.get_gX2();

		nnsBCH<m, t> nnsbch;
		nnsbch.print_info();
		cout << "nnsbch.get_gX2() = " << nnsbch.get_gX2();

		const int b = 3;
		npBCH<m, t, b> npbch;
		npbch.print_info();
		cout << "npbch.get_gX2() = " << npbch.get_gX2();
	}
	static void BCH_info3() {
		const int m = 5;	// GF field of 2^m
		const int t = 3;	// error that can correct
		GF2e<m>::init();
		BCH<m, t> bch;
		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();

		// get standard form of generator matrix for matlab
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<int> nat(1, n, 'N');
		G.GJE_4_GF2_right(nat);
		cout << "G = " << G;
	}
	static void RS_test() {
		typedef GF2e<3> ty;
		ty::init();

		RS<3, 3> rs;
		rs.print_info();
		int n = rs.get_n();
		int k = rs.get_k();
		int d = rs.get_d();

		Matrix<ty> u(1, k, { 5,1,6 });
		cout << "u" << u;
		Matrix<ty> v = rs.encode(u);
		cout << "v" << v;
		Matrix<ty> e(1, n, { 0,0,0,0,7,1,0 });
		Matrix<ty> r = e + v;
		cout << "r" << r;

		Matrix<ty> u_hat = rs.decode_BM(r);
		cout << "u_hat" << u_hat;

		cout << "rs.get_generator_matrix()" << rs.get_generator_matrix();
		cout << "rs.get_parity_matrix()" << rs.get_parity_matrix();
		cout << "rs.get_generator_matrix() * rs.get_parity_matrix()" << rs.get_generator_matrix() * rs.get_parity_matrix().Transpose();
	}
	static void RS_test2() {
		typedef GF2e<3> ty;
		ty::init();

		RS<3, 3> rs;
		rs.print_info();
		int n = rs.get_n();
		int k = rs.get_k();
		int d = rs.get_d();

		Matrix<ty> u(1, k, { 0,2,7 });
		cout << "u" << u;
		Matrix<ty> v = rs.encode_by_code_locators(u);
		cout << "v" << v;
		Matrix<ty> u_back = rs.v2u_by_code_locators(v);
		cout << "u_back" << u_back;

		Matrix<ty> e(1, 7, { 0,0,4,0,0,7,0 });
		Matrix<ty> r = v + e;
		cout << "r" << r;
		Matrix<ty> u_hat = rs.decode_BM_by_code_locators(r);
		cout << "u_hat" << u_hat;
	}
	static void RS_test3() {
		typedef GF2e<3> ty;
		ty::init();

		RS<3, 3> rs;
		rs.print_info();
		int n = rs.get_n();
		int k = rs.get_k();
		int d = rs.get_d();

		Matrix<ty> u(1, k, { 0,2,7 });
		cout << "u" << u;
		Matrix<ty> v = rs.encode_by_code_locators(u);
		cout << "v" << v;

		// error pattern
		Matrix<ty> e(1, n, { 0,0,4,0,0,7,0 });
		cout << "e" << e;

		Matrix<ty> r = v + e;
		cout << "r" << r;
		Matrix<ty> u_hat = rs.decode_GS_by_code_locators(r, 1);
		cout << "u_hat" << u_hat;

		Matrix<ty> v_hat = rs.encode_by_code_locators(u_hat);
		cout << "v_hat" << v_hat;
		cout << "v_hat - v" << v_hat - v;
	}
	static void RS_test4() {
		const int m = 6;
		typedef GF2e<m> ty;
		ty::init();

		RS<m, 21> rs;
		rs.print_info();
		int n = rs.get_n();		// n=(1<<m)-1
		int k = rs.get_k();
		int d = rs.get_d();

		Matrix<ty> u(1, k);
		for (int i = 0; i < k; ++i) {
			u(i) = my::rand_int(0, (1 << m) - 1);
		}

		cout << "u" << u;
		Matrix<ty> v = rs.encode_by_code_locators(u);
		cout << "v" << v;

		// error pattern
		int error_num = 26;
		Matrix<ty> e(1, n, '0');
		for (int i = 0; i < error_num; ++i) {
			e(i) = my::rand_int(1, (1 << m) - 1);	// random error
		}
		cout << "e" << e;

		Matrix<ty> r = v + e;
		cout << "r" << r;
		const int multiplicity = 5;
		int tau = rs.compute_tau(multiplicity);
		cout << "tau = " << tau << endl;

		clock_t start, end;
		start = clock();	//start timing

		Matrix<ty> u_hat = rs.decode_GS_by_code_locators(r, multiplicity);

		end = clock();		//end timing
		cout.setf(ios::scientific);
		cout << scientific << "time consume: " << ((double)end - start) / CLOCKS_PER_SEC << "s\n\n";
		cout.unsetf(ios::floatfield);

		cout << "mod2_special::max_factorial = " << mod2_special::max_factorial << endl;

		//for (int i = 0; i < 10000; ++i)
		//	rs.decode_GS_by_code_locators(r, 5);		// testing memory leakage
		cout << "u_hat" << u_hat;
		cout << "u_hat - u" << u_hat - u;

		Matrix<ty> v_hat = rs.encode_by_code_locators(u_hat);
		cout << "v_hat" << v_hat;
		cout << "v_hat - r" << v_hat - r;
	}
	static void RS_test5() {
		typedef GF2e<6> ty;
		ty::init();

		RS<6, 21> rs;

		Matrix<int> multiplicity(1, 16);
		Matrix<int> tau_m(1, 16);
		Matrix<int> l_m(1, 16);
		for (int i = 0; i < 16; ++i) {
			multiplicity(i) = i + 1;
			tau_m(i) = rs.compute_tau(i + 1);
			l_m(i) = rs.compute_lm(i + 1);
		}
		cout << "multiplicity" << multiplicity;
		cout << "tau_m" << tau_m;
		cout << "l_m" << l_m;
	}
	static void RS_test6() {
		typedef GF2e<6> ty;
		ty::init();

		RS<6, 21> rs;
		// 22 error, BM decoding make errors
		Matrix<ty> hdr(1, 63, { 16,4,0,0,0,0,0,0,0,0,0,0,16,16,1,0,0,0,0,16,2,0,32,0,2,0,16,2,4,0,0,40,0,0,0,2,32,9,0,0,16,0,0,0,0,0,0,0,0,0,0,16,0,8,4,0,0,0,0,0,4,0,16 });
		cout << "hdr" << hdr;
		cout << "hdr.size()=" << hdr.size() << endl;
		Matrix<ty> u_hat_2 = rs.decode_naive_by_code_locators(hdr);
		cout << "u_hat_2" << u_hat_2;		// decoding fail will return all 0 vector
		Matrix<ty> u_hat = rs.decode_BM_by_code_locators(hdr);
		cout << "u_hat" << u_hat;
	}
	static void RS_test7() {
		const int m = 4;
		const int k = 9;
		typedef GF2e<m> ty;
		ty::init();

		RS<m, k> rs;
		rs.print_info();
		int n = rs.get_n();
		int d = rs.get_d();

		Matrix<int> info_set(1, 9, { 2,1,4,6,10,9,8,13,7 });
		info_set.sort();
		rs.generate_systematic_generator_any_pos(info_set);

		Matrix<ty> u(1, k);
		for (int i = 0; i < k; ++i) {
			u(i) = my::rand_int(0, (1 << m) - 1);
		}
		cout << "u" << u;
		Matrix<ty> v = rs.encode_by_Lagrange_interpolation_any_pos(u);
		cout << "v" << v;

		// error pattern
		Matrix<ty> e(1, n, '0');
		for (int i = 0; i < d / 2; ++i) {
			e(i) = my::rand_int(1, (1 << m) - 1);
		}
		cout << "e" << e;

		Matrix<ty> r = v + e;
		cout << "r" << r;
		Matrix<ty> v_hat = rs.decode_BM_v(r);
		cout << "v_hat" << v_hat;

		Matrix<ty> u_hat = rs.v2u_by_Lagrange_interpolation_any_pos(v_hat);
		cout << "u_hat" << u_hat;
		cout << "u_hat - u" << u_hat - u;

	}
	static void RS_code_Matrix() {

		typedef GF2e<3> ty;
		ty::init();

		Matrix<ty> G1(3, 7, {
			1,1,1,1,1,1,1,
			1,2,3,4,5,6,7,
			1,3,5,7,2,4,6
			});
		cout << "G1" << G1;
		cout << "G1*G1.Transpose()" << G1 * G1.Transpose();

		// extended RS code, to be finished in a new class, important extension
		Matrix<ty> G(3, 8, {
			1,1,1,1,1,1,1,1,
			1,2,3,4,5,6,7,0,
			1,3,5,7,2,4,6,0
			});
		cout << "G" << G;
		cout << "G*G.Transpose()" << G * G.Transpose();

		// parity check matrix of RS code
		Matrix<ty> H(5, 8, {
			1,2,3,4,5,6,7,0,
			1,3,5,7,2,4,6,0,
			1,4,7,3,6,2,5,0,
			1,5,2,6,3,7,4,0,
			1,1,1,1,1,1,1,1
			});
		cout << "H" << H;
		cout << "G * H.Transpose()" << G * H.Transpose();
	}
	static void polynomial_v2_test() {
		typedef GF2e<4> ty;
		ty::init();

		polynomial_2v<ty> p1(Matrix<ty>(2, 3, {
			4,2,1,
			3,7,2
			}));
		cout << "p1" << p1;
		polynomial_2v<ty> p2(Matrix<ty>(1, 3, {
			3,0,1
			}));
		cout << "p2" << p2;
		polynomial_2v<ty> p3 = p1 * p2;
		cout << "p3" << p3;
		cout << "p3/p2" << p3 / p2;
		cout << "p3%p2" << p3 % p2;

		int x = 4, y = 3;
		cout << "evaluate=" << p1.evaluate(4, 3) << endl;

		ipair::weighted_i = 1;
		ipair::weighted_j = 3;
		ipair::is_lex = true;
		polynomial_2v<ty> p4(Matrix<ty>(5, 2, {
			4,2,
			3,1,
			1,1,
			1,0,
			1,0
			}));
		cout << "p4" << p4;
		cout << "p4.leading_monomial()=" << p4.leading_monomial() << endl;

		polynomial_2v<ty> p5;
		cout << "p5" << p5;

		cout << "p5*p4" << p5 * p4;
	}
	static void ipair_test() {
		ipair::weighted_i = 1;
		ipair::weighted_j = 3;
		ipair::is_lex = true;
		ipair p1(2, 0);
		ipair p2(0, 1);
		cout << "p1 < p2 = " << (p1 < p2) << endl;

		ipair top_p(5, 1);
		Matrix<ipair> sorted_ipair = ipair::under_top(top_p, 1);
		cout << "sorted_ipair" << sorted_ipair;
		ipair move = top_p;
		cout << "move" << move << endl;
	}
	static void koetter_test() {
		typedef GF2e<4> ty;
		ty::init();

		// begin (given L, (alpha_i,beta_i)_{i=1:n}, (m_i)_{i=1:n}, (1, k-1) wdeg monomial order
		int L = 5;
		// (m-1,1)-lex order
		int n = 15;
		int k = 7;
		Matrix<ty> alpha(1, n);
		for (int i = 0; i < n; ++i) {
			alpha(i) = 6;
		}
		Matrix<int> beta(1, n);
		for (int i = 0; i < n; ++i) {
			beta(i) = 12;
		}
		Matrix<int> m(1, n);
		for (int i = 0; i < n; ++i) {
			m(i) = 2;
		}

		Matrix<polynomial_2v<ty>> g(1, L + 1);
		for (int j = 0; j <= L; ++j) {
			g(j) = polynomial_2v<ty>(1, j + 1, '0');
			g(j)(0, j) = 1;
		}
		cout << "g" << g;

		Matrix<ty> Delta(1, L + 1, '0');	// discrepancy
		for (int i = 0; i < n; ++i) {
			ipair::weighted_i = (m(i) - 1) == 0 ? 1 : (m(i) - 1);		// prevent ms(i)==1, which induce weighted_i be 0
			ipair::weighted_j = 1;
			ipair::is_lex = true;

			ipair top_p(m(i) - 1, 1);
			Matrix<ipair> sorted_ipair = ipair::under_top(top_p, m(i));	// by (m_i-1,1) lex order
			int sorted_ipair_size = sorted_ipair.size();

			for (int p = 0; p < sorted_ipair_size; ++p) {
				int r = sorted_ipair(p).i;
				int s = sorted_ipair(p).j;
				Matrix<int> J(1, 0, 'v');
				for (int j = 0; j <= L; ++j) {
					Delta(j) = g(j).Hasse_D(r, s, alpha(i), beta(i));	// jth discrepancy
					if (Delta(j) != 0) {
						J.push_back(j);
					}
				}

				if (J.size() != 0) {

					ipair::weighted_i = 1;
					ipair::weighted_j = k - 1;
					ipair::is_lex = false;			// change it back to (1,k-1) wdeg monomial order, revlex

					int Js = J.size();
					ipair min_monomial = g(J(0)).leading_monomial();
					int j_star = 0;
					for (int j = 1; j < Js; ++j) {
						ipair can_monomial = g(J(j)).leading_monomial();
						if (min_monomial < can_monomial);
						else {
							min_monomial = can_monomial;
							j_star = j;
						}
					}
					polynomial_2v<ty> f = g(J(j_star));

					for (int j = 0; j < Js; ++j) {
						if (j != j_star) {
							g(j) = Delta(J(j_star)) * g(j) - Delta(J(j)) * f;	// no change in wdeg
						}
						else if (j == j_star) {
							polynomial_2v<ty> tmp(2, 1, '0');
							tmp(0, 0) = -alpha(i);
							tmp(1, 0) = 1;		// tmp=x - alpha_i
							g(j) = Delta(J(j_star)) * tmp * f;					// wdeg increases by 1
						}
					}
				}
			}
		}
		ipair::weighted_i = 1;
		ipair::weighted_j = k - 1;
		ipair::is_lex = false;			// change it back to (1,k-1) wdeg monomial order, revlex

		ipair min_monomial = g(0).leading_monomial();
		int min_ind = 0;
		for (int j = 1; j <= L; ++j) {
			ipair can_monomial = g(j).leading_monomial();
			if (min_monomial < can_monomial);
			else {
				min_monomial = can_monomial;
				min_ind = 0;
			}
		}
		polynomial_2v<ty> Q0 = g(min_ind);
		cout << "Q0" << Q0;

		cout << "Q0.evaluate(6, 12)=" << Q0.evaluate(6, 12) << endl;
	}
	static void Q_next_test() {
		typedef GF<19> ty;
		ty::init();

		polynomial_2v<ty> Q(Matrix<ty>(6, 5, {
			4,14,14,2,17,
			12,14,13,11,0,
			5,9,1,1,0,
			11,16,0,0,0,
			8,8,0,0,0,
			13,0,0,0,0
			}));
		polynomial_2v<ty> Q_next = RR_root_finding<ty>::compute_Q_next(Q, 18);
		cout << "Q_next" << Q_next;
		cout << "check" << endl;
	}
	static void RR_root_finding_test() {
		typedef GF<19> ty;
		ty::init();
		polynomial_2v<ty> Q(Matrix<ty>(6, 5, {
			4,14,14,2,17,
			12,14,13,11,0,
			5,9,1,1,0,
			11,16,0,0,0,
			8,8,0,0,0,
			13,0,0,0,0
			}));
		for (int qq = 0; qq < 1000000000; ++qq) 		// memory loss, donnot know why
		{
			cout << "qq=" << qq << endl;
			RR_root_finding<ty> problem_RR;
			problem_RR.set_input(Q, 2, 19);
			problem_RR.solve();

		}

		//cout << "problem_RR.y_roots" << problem_RR.y_roots;
	}
	static void n_choose_k_test() {
		int n = 127;
		int k = 9;
		cout << my::n_choose_k(n, k) << endl;
	}
	static void num_test() {
		int n = 7;
		int k = 3;
		int m = 1;
		int v = k - 1;
		double L = (sqrt(n / (double)v * m * (m + 1) + pow(((v + 2) / 2.0 / v), 2)) - (v + 2.0) / 2 / v);

		cout << "L=" << L << endl;
	}
	static void matrix_push_back_test() {
		for (int j = 0; j < 10000; ++j) {
			Matrix<polynomial_2v<int>> v;
			for (int i = 0; i < 40; ++i) {
				cout << "---------- i = " << i << " ------------\n";
				v.push_back(Matrix<int>(2, 5, '1'));
				//cout << v << endl;
			}
		}

		/*for (int j = 0; j < 10000000; ++j) {
			Matrix<int> v(1, 0, '0');
			for (int i = 0; i < 20; ++i) {
				v.push_back(4);
			}
			cout << "j=" << j << endl;
		}*/
	}
	static void GF2e_to_GF_test() {
		typedef GF2e<3> ty;
		ty::init();
		GF_trans<ty, 1 << 3> gt;
		Matrix<ty> v(1, 3, { 1,4,3 });
		Matrix<GF2> x = gt.to_bits(v);
		cout << "v" << v;
		cout << "x" << x;
		Matrix<ty> v_hat = gt.to_symbol(x);
		cout << "v_hat" << v_hat;
	}
	static void RS_simulate_test() {
		const int m = 6;
		const int k = 21;
		typedef GF2e<m> ty;
		ty::init();
		Matrix<my_double> SNR_db(5.0, 0.5, 7.0, 'd');
		cout << "SNR_db" << SNR_db;
		simulation::_error_frame_upper_bound = 100;
		simulation::_iteration = 5000;
		int SNR_db_size = SNR_db.size();
		Matrix<my_double> FER_v(1, SNR_db_size);
		Matrix<my_double> frame_proceed_v(1, SNR_db_size);

		clock_t start, end;
		start = clock();	//start timing

		for (int i = 0; i < SNR_db_size; ++i) {
			simulation::RS_code<m, k>(SNR_db(i), RS_decode_type::GS, 5);		// basically all-right
			FER_v(i) = simulation::FER;
			frame_proceed_v(i) = simulation::frame_proceed;
		}

		cout << "FER_v" << FER_v;
		cout << "frame_proceed_v" << frame_proceed_v;


		end = clock();		//end timing
		cout.setf(ios::scientific);
		cout << scientific << "time consume: " << ((double)end - start) / CLOCKS_PER_SEC << "s\n\n";
		cout.unsetf(ios::floatfield);
	}
	static void mod_2_special_test() {
		cout << "mod2_special::n_choose_k(x, y) = " << mod2_special::n_choose_k(0, 0) << endl;
		cout << "mod2_special::n_choose_k(x, y) = " << mod2_special::n_choose_k(3, 1) << endl;
		cout << mod2_special::max_factorial << endl;
	}

	static void polynomial_map_test() {
		typedef GF2e<6> ty;
		ty::init();

		monomial_order::weighted_i = 1;
		monomial_order::weighted_j = 1;
		monomial_order::is_lex = true;

		polynomial_2v_map<ty> p(Matrix<ty>(5, 5, {
			0,5,1,5,7,
			7,0,3,0,0,
			0,10,9,0,0,
			7,0,3,2,0,
			8,4,21,0,0
			}));
		cout << "p=" << p;
		cout << "a<b? " << monomial_order::lt(make_pair(1, 0), make_pair(4, 0)) << endl;

		cout << "p.leading_monomial() = " << p.leading_monomial().first << ", " << p.leading_monomial().second << endl;
		cout << "p=" << p;
		cout << "p(4, 0) = " << p(4, 0) << endl;
		cout << "p(1, 1) = " << p(1, 1) << endl;
		cout << "p=" << p;

		cout << "p.max_x_pow=" << p.max_x_pow << ", p.max_y_pow=" << p.max_y_pow << endl;

		polynomial_2v<ty> p2(Matrix<ty>(5, 5, {
			0,5,1,5,7,
			7,0,3,0,0,
			0,10,9,0,0,
			7,0,3,2,0,
			8,4,21,0,0
			}));
		cout << "p2" << p2;

		ty alpha = 11;
		ty beta = 37;
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 4; ++j) {

				cout << "(" << i << "," << j << ") " << "(p compute) - (p2 compute)=" << \
					p.Hasse_D_for_GF2e(i, j, alpha, beta) - p2.Hasse_D_for_GF2e(i, j, alpha, beta) << endl;
			}
		}

		cout << "p=" << p;
		cout << "p.max_x_pow=" << p.max_x_pow << ", p.max_y_pow=" << p.max_y_pow << endl;

		polynomial_2v_map<ty> p3_empty;
		cout << "p3=" << p3_empty;
		pair<int, int> p3_lead = p3_empty.leading_monomial();
		cout << "p3_lead: (" << p3_lead.first << ", " << p3_lead.second << ")" << endl;
	}
	static void polynomial_umap_test() {
		typedef GF2e<6> ty;
		ty::init();

		monomial_order::weighted_i = 1;
		monomial_order::weighted_j = 1;
		monomial_order::is_lex = true;

		polynomial_2v_umap<ty> p(Matrix<ty>(5, 5, {
			0,5,1,5,7,
			7,0,3,0,0,
			0,10,9,0,0,
			7,0,3,2,0,
			8,4,21,0,0
			}));
		cout << "p=" << p;
		cout << "a<b? " << monomial_order::lt(make_pair(1, 0), make_pair(4, 0)) << endl;

		cout << "p.leading_monomial() = " << p.leading_monomial().first << ", " << p.leading_monomial().second << endl;
		cout << "p=" << p;
		cout << "p(4, 0) = " << p(4, 0) << endl;
		cout << "p(1, 1) = " << p(1, 1) << endl;
		cout << "p=" << p;

		polynomial_2v<ty> p2(Matrix<ty>(5, 5, {
			0,5,1,5,7,
			7,0,3,0,0,
			0,10,9,0,0,
			7,0,3,2,0,
			8,4,21,0,0
			}));
		cout << "p2" << p2;

		ty alpha = 11;
		ty beta = 37;
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 4; ++j) {

				cout << "(" << i << "," << j << ") " << "(p compute) - (p2 compute)=" << \
					p.Hasse_D_for_GF2e(i, j, alpha, beta) - p2.Hasse_D_for_GF2e(i, j, alpha, beta) << endl;
			}
		}

		cout << "p=" << p;

		polynomial_2v_umap<ty> p3_empty;
		cout << "p3=" << p3_empty;
		pair<int, int> p3_lead = p3_empty.leading_monomial();
		cout << "p3_lead: (" << p3_lead.first << ", " << p3_lead.second << ")" << endl;

		p3_empty = p;
		p *= 2;
		cout << "p" << p;
		p /= 2;
		cout << "p" << p;

		p3_empty -= p;		// -= can not use to itself
		cout << "p3_empty" << p3_empty;

		polynomial_2v_umap<ty> up(Matrix<ty>(5, 5, {
			0,5,1,5,7,
			7,0,3,0,0,
			0,10,9,0,0,
			7,0,3,2,0,
			8,4,21,0,0
			}));

		polynomial_2v_map<ty> mp(Matrix<ty>(5, 5, {
			0,5,1,5,7,
			7,0,3,0,0,
			0,10,9,0,0,
			7,0,3,2,0,
			8,4,21,0,0
			}));
		cout << "up" << up;
		cout << "mp" << mp;

		bool test_eq = up.contain_Matrix(Matrix<ty>(5, 5, {
			0,5,1,5,7,
			7,0,3,0,0,
			0,10,9,0,0,
			7,0,3,2,0,
			8,4,21,0,0
			}));
		cout << "test_eq = " << test_eq << endl;

		up.times_x_plus_alpha(2);
		mp.times_x_plus_alpha(2);		// same !

		cout << "up" << up;
		cout << "mp" << mp;


	}
	static void Lagrange_test() {
		typedef GF2e<4> ty;
		ty::init();
		Matrix<ty> alpha(1, 8, { 3,2,6,8,15,13,10,4 });
		int alpha_size = alpha.size();
		polynomial<ty> p = Lagrange_interpolation::generate(alpha, 7);		// p(4)=1, index is 7
		cout << "p" << p;
		for (int i = 0; i < alpha_size; ++i) {
			cout << "p(" << alpha(i) << ")=" << p.evaluate(alpha(i)) << endl;
		}

		ty beta = 9;
		cout << "p(" << beta << ")=" << p.evaluate(beta) << endl;

		// my evaluation
		ty result = 1;
		ty alpha_i;
		alpha_i.set_by_alpha_power(3); // alpha^3
		ty alpha_j = beta;

		for (int w = 0; w < alpha_size; ++w) {
			ty alpha_j_prime = alpha(w);
			if (alpha_j_prime != alpha_i) {
				result *= (alpha_j - alpha_j_prime) / (alpha_i - alpha_j_prime);
			}
		}
		cout << "result = " << result << endl;

		//RS<4, 10> rs;
	}
	static void BCH_test2() {
		typedef GF2 ty;
		GF2e<4>::init();

		BCH<4, 2> bch;
		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();
		int d = bch.get_d();

		Matrix<ty> u(1, k, { 1,0,1,0,1,0,1 });
		cout << "u" << u;
		Matrix<ty> v = bch.encode(u);
		cout << "v" << v;
		Matrix<ty> e(1, n, { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 });
		Matrix<ty> r = e + v;
		cout << "r" << r;

		Matrix<ty> u_hat = bch.decode(r);
		cout << "u_hat" << u_hat;
		cout << "u_hat - u" << u_hat - u;

		RS<4, 11> rs;
		rs.print_info();
		int k_prime = rs.get_k();
		Matrix<int> information_set_ind(1, k_prime, { 3,5,11,6,9,2,14,0,1,10,8 });
		information_set_ind.sort();
		//cout << "information_set_ind" << information_set_ind;
		rs.generate_systematic_generator_any_pos(information_set_ind);
		Matrix<ty> information_set = r.get_cols(information_set_ind);
		//cout << "information_set" << information_set;

		Matrix<GF2e<4>> super_information_set(1, k_prime);
		for (int i = 0; i < k_prime; ++i) {
			super_information_set(i) = r(information_set_ind(i));
		}
		//cout << "super_information_set" << super_information_set;
		Matrix<GF2e<4>> rs_encoded_v = rs.encode_by_Lagrange_interpolation_any_pos(super_information_set);
		cout << "rs_encoded_v" << rs_encoded_v;		// same as BCH code
	}
	static void BCH_test3() {
		const int m = 4;
		typedef GF2e<m> ty;
		ty::init();

		BCH<m, 1> bch;		// which is (7,4) Hamming code

		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();
		int d = bch.get_d();
		int r = bch.get_r();

		Matrix<GF2> GM = bch.get_generator_matrix();
		//cout << "GM" << GM;
		Matrix<GF2> PM = bch.get_parity_matrix();
		cout << "PM" << PM;

		Matrix<GF2> v(1, n, 'b');
		cout << "v" << v;

		Matrix<GF2> v_hat = bch.decode_v(v);
		cout << "v_hat" << v_hat;
		cout << "bch.is_in_v_space(v_hat) = " << bch.is_in_v_space(v_hat) << endl;
#ifdef use_my_double
		unsigned long long ope_num_before = my_double_auxiliary_storage::operation_number;
#endif
		// viterbi algorithm
		int pow_2_r = 1 << r;
		cout << "pow_2_r = " << pow_2_r << endl;
		Matrix<int> past_state(pow_2_r, n);
		past_state.reset(-3000);		// negtive 1000 means invalid past state
		Matrix<my_double> correlation_distance(pow_2_r, n);
		correlation_distance.reset(3000);		// 3000 means invalid correlation distance

		// for i==0, the first received symbol
		int PM_col_repsent = 0;
		for (int k = 0; k < r; ++k) {
			PM_col_repsent <<= 1;
			PM_col_repsent += (int)PM(k, 0);
		}
		//cout << "PM_col_repsent = " << PM_col_repsent << endl;


		// initialize correlation distance
		correlation_distance(0, 0) = v_hat(0) != 0;		// using hamming distance first
		correlation_distance(PM_col_repsent, 0) = v_hat(0) == 0;

		// initialize past state
		past_state(0, 0) = 0;
		past_state(PM_col_repsent, 0) = 0;

		// starting trellises, special for bch code (cyclical code)

		for (int i = 1; i < r; ++i) {
			int save_ope_param = r - i - 1;
			PM_col_repsent = 0;
			for (int k = 0; k <= i; ++k) {
				PM_col_repsent <<= 1;
				PM_col_repsent += (int)PM(k, i);
			}
			PM_col_repsent <<= save_ope_param;

			int jump_state_num = 1 << save_ope_param;
			for (int j = 0; j < pow_2_r; j += jump_state_num) {
				// update correlation distance and past state
				my_double can_1 = correlation_distance(j, i - 1) + (v(i) != 0);
				my_double can_2 = correlation_distance(j ^ PM_col_repsent, i - 1) + (v(i) == 0);
				if (can_1 < can_2) {
					correlation_distance(j, i) = can_1;
					past_state(j, i) = j;
				}
				else {
					correlation_distance(j, i) = can_2;
					past_state(j, i) = j ^ PM_col_repsent;
				}

			}
		}

		// trellises with full state

		for (int i = r; i < n - r; ++i) {
			PM_col_repsent = 0;
			for (int k = 0; k < r; ++k) {
				PM_col_repsent <<= 1;
				PM_col_repsent += (int)PM(k, i);
			}

			// can be shrinked when consider both the initial state and final state are 0
			for (int j = 0; j < pow_2_r; ++j) {
				// update correlation distance and past state
				my_double can_1 = correlation_distance(j, i - 1) + (v(i) != 0);
				my_double can_2 = correlation_distance(j ^ PM_col_repsent, i - 1) + (v(i) == 0);
				if (can_1 < can_2) {
					correlation_distance(j, i) = can_1;
					past_state(j, i) = j;
				}
				else {
					correlation_distance(j, i) = can_2;
					past_state(j, i) = j ^ PM_col_repsent;
				}

			}
		}

		// ending trellises, special for bch code (cyclical code)

		int new_start = n - r > r ? n - r : r;
		for (int i = new_start; i < n; ++i) {
			int save_ope_param = i - n + r;
			PM_col_repsent = 0;
			for (int k = save_ope_param; k < r; ++k) {
				PM_col_repsent <<= 1;
				PM_col_repsent += (int)PM(k, i);
			}

			// can be shrinked when consider both the initial state and final state are 0
			int range = pow_2_r >> (save_ope_param + 1);
			for (int j = 0; j < range; ++j) {
				// update correlation distance and past state
				my_double can_1 = correlation_distance(j, i - 1) + (v(i) != 0);
				my_double can_2 = correlation_distance(j ^ PM_col_repsent, i - 1) + (v(i) == 0);
				if (can_1 < can_2) {
					correlation_distance(j, i) = can_1;
					past_state(j, i) = j;
				}
				else {
					correlation_distance(j, i) = can_2;
					past_state(j, i) = j ^ PM_col_repsent;
				}

			}
		}

		//cout << "correlation_distance" << correlation_distance;
		//cout << "past_state" << past_state;

		// get the decoded codeword, final state is 0
		Matrix<GF2> v_hat_viterbi(1, n);
		int tract_state = 0;
		for (int i = n - 1; i >= 0; --i) {
			v_hat_viterbi(i) = past_state(tract_state, i) != tract_state;
			tract_state = past_state(tract_state, i);
		}

#ifdef use_my_double
		unsigned long long ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "ope_num = " << ope_num_after - ope_num_before << endl;
#endif

		cout << "v_hat_viterbi" << v_hat_viterbi;
		cout << "v_hat - v_hat_viterbi" << v_hat - v_hat_viterbi;

		cout << "bch.is_in_v_space(v_hat_viterbi) = " << bch.is_in_v_space(v_hat_viterbi) << endl;
	}
	static void BCH_test4() {
		const int m = 4;
		const int t = 3;
		typedef GF2e<m> ty;
		ty::init();

		BCH<m, t> bch;

		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();
		int d = bch.get_d();
		int r = bch.get_r();

		Matrix<GF2> GM = bch.get_generator_matrix();
		cout << "GM" << GM;
		Matrix<GF2> PM = bch.get_parity_matrix();
		cout << "PM" << PM;

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);

		cout << "v" << v;
		Matrix<my_double> c = BPSK::modulation(v);
		cout << "c" << c;
		Matrix<my_double> e(1, n);
		for (int i = 0; i < n; ++i) {
			e(i) = my::rand_ga() * 0.0;		// in this case v_hat is error but v_hat_viterbi is correct
		}

		Matrix<my_double> recv = c + e;
		cout << "recv" << recv;
		Matrix<GF2> hdr = BPSK::demodulation(recv);
		cout << "hdr" << hdr;

		Matrix<GF2> v_hat = bch.decode_v(hdr);
		cout << "v_hat" << v_hat;
		cout << "bch.is_in_v_space(v_hat) = " << bch.is_in_v_space(v_hat) << endl;

		Viterbi_cyclic_code_r vit_bch(bch.get_parity_matrix());
		Matrix<GF2> v_hat_viterbi = vit_bch.decode_v(recv);
		//Matrix<GF2> v_hat_viterbi = Viterbi_cyclic_code_k::decode_v_brutle_forcce(bch, recv);
		cout << "v_hat_viterbi" << v_hat_viterbi;
		cout << "bch.is_in_v_space(v_hat) = " << bch.is_in_v_space(v_hat_viterbi) << endl;

		cout << "v_hat -v" << v_hat - v;
		cout << "v_hat_viterbi - v" << v_hat_viterbi - v;
	}
	static void BCH_test5() {
		cout << boolalpha;

		const int m = 6;
		const int t = 2;
		typedef GF2e<m> ty;
		ty::init();

		BCH<m, t> bch;

		bch.print_info();
		int n = bch.get_n();
		int k = bch.get_k();
		int d = bch.get_d();
		int r = bch.get_r();

		Matrix<GF2> GM = bch.get_generator_matrix();
		cout << "GM" << GM;
		Matrix<GF2> PM = bch.get_parity_matrix();
		cout << "PM" << PM;

		// for cyclic code, permutation may make decoding more complex
		/*Matrix<int> permute_test = Matrix<int>(1, n, {
			0,	1,	2,	3,	4,	5,	6,	7,	8,	9,	10,	11,	12,	13,	14,	15,	24,	16,	25,	17,	26,	18,	20,	27,	29,	19,	21,	22,	23,	28,	30
			});*/
			//PM.permute_col(permute_test);
			//cout << "PM" << PM;

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//v.permute(permute_test);

		cout << "v" << v;
		Matrix<my_double> c = BPSK::modulation(v);
		cout << "c" << c;
		Matrix<my_double> e(1, n);
		for (int i = 0; i < n; ++i) {
			e(i) = my::rand_ga() * 0.6;		// in this case v_hat is error but v_hat_viterbi is correct
		}

		Matrix<my_double> recv = c + e;
		cout << "recv" << recv;
		Matrix<GF2> hdr = BPSK::demodulation(recv);
		cout << "hdr" << hdr;

		Matrix<GF2> v_hat = bch.decode_v(hdr);			// using BM decoding
		cout << "v_hat" << v_hat;
		cout << "bch.is_in_v_space(v_hat) = " << bch.is_in_v_space(v_hat) << endl;

		int list_size = 4096;
		/*
		Brute_force brute_bch(GM);
		Matrix<GF2> list_v_hat_brute = brute_bch.decode_v(recv, list_size);
		*/

		Viterbi_optimized vit_bch_r(bch.get_parity_matrix());
		Matrix<GF2> list_v_hat_viterbi_new = vit_bch_r.decode_v(recv, list_size, list_size);

		Viterbi_unordered_map vit_bch3(PM);
		vit_bch3.initialize_state();		// this must be called before decode_v_repeatedly
		Matrix<GF2> list_v_hat_viterbi_test = vit_bch3.decode_v_once(recv, list_size);
		//Matrix<GF2> list_v_hat_viterbi_test = vit_bch3.decode_v_repeatedly(recv, 16);

		//Viterbi_Sorted_vector vit_bch4(PM);		// this is not valid for k<n-k, not fixing, discarding this class
		//vit_bch4.initialize_state();
		//Matrix<GF2> list_v_hat_viterbi_test2 = vit_bch4.decode_v_repeatedly_space_saving(recv, list_size);

		//cout << "(list_v_hat_viterbi_bf - list_v_hat_viterbi_new)"
		//	<< (list_v_hat_brute - list_v_hat_viterbi_new);
		//cout << "(list_v_hat_viterbi_bf - list_v_hat_viterbi_test)"
		//	<< (list_v_hat_brute - list_v_hat_viterbi_test);
		//cout << "(list_v_hat_viterbi_test - list_v_hat_viterbi_new)"
		//	<< (list_v_hat_viterbi_test - list_v_hat_viterbi_new);

		/*cout << "(list_v_hat_viterbi_bf - list_v_hat_viterbi_new).isZero() = "
			<< (list_v_hat_brute - list_v_hat_viterbi_new).isZero() << endl;*/
			/*cout << "(list_v_hat_viterbi_bf - list_v_hat_viterbi_test).isZero() = "
				<< (list_v_hat_brute - list_v_hat_viterbi_test).isZero() << endl;*/
		cout << "(list_v_hat_viterbi_test - list_v_hat_viterbi_new).isZero() = "
			<< (list_v_hat_viterbi_test - list_v_hat_viterbi_new).isZero() << endl;

		// actually okay, difference is induced by two identical float metric
		/*cout << "(list_v_hat_viterbi_test - list_v_hat_viterbi_new)"
			<< (list_v_hat_viterbi_test - list_v_hat_viterbi_new) << endl;*/


			/*cout << "(list_v_hat_viterbi_bf - list_v_hat_viterbi_test2)"
				<< (list_v_hat_brute - list_v_hat_viterbi_test2);*/

		Matrix<GF2> v_hat_viterbi = list_v_hat_viterbi_new.get_row(0);	// the best decoded codeword
		cout << "v_hat_viterbi" << v_hat_viterbi;
		cout << "bch.is_in_v_space(v_hat) = " << bch.is_in_v_space(v_hat_viterbi) << endl;

		cout << "v_hat_viterbi - v_hat" << v_hat_viterbi - v_hat;
		cout << "v_hat_viterbi - v" << v_hat_viterbi - v;
	}
	static void unordered_map_testing() {
		const int m = 6;
		typedef GF2e<6> ty;
		ty::init();

		unordered_map<int, ipair> um;
		cout << "um.size() = " << um.size() << endl;
		// well, visit a not found key value will call value's constructor with no parameter and create it in unordered_map!

		// for interger, parameterless constructor will assign it to zero, the c++ rule

		cout << "um[0] = " << um[0] << endl;
		cout << "um[0] = " << um[345] << endl;
		cout << "um[0] = " << um[2] << endl;
		cout << "um[0] = " << um[98] << endl;
		cout << "um.size() = " << um.size() << endl;
	}
	static void BCH_non_Gaussian_elimination() {
		const int m = 6;
		typedef GF2e<m> ty;
		ty::init();

		const int t = 3;
		BCH<m, t> bch;
		int n = bch.get_n();
		int k = bch.get_k();
		int d = bch.get_d();

		bch.print_info();
		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		cout << "v" << v;
		Matrix<my_double> c = BPSK::modulation(v);
		AWGN::sigma = 0.6;
		Matrix<my_double> r = AWGN::pass(c);
		cout << "r" << r;
		Matrix<my_double> r_abs = r.get_abs();

		Matrix<int> permute_ind = r_abs.sort_with_ind('>');		// from big to small
		cout << "permute_ind" << permute_ind;
		const int k_prime = (1 << m) - 1 - 2 * t;
		RS<m, k_prime> rs;
		int r_prime = rs.get_r();		// equals to 2*t

		// rs generator matrix from most reliable information set of bch received codeword
		rs.generate_systematic_generator_any_pos(permute_ind.get_part(0, 0, 0, k_prime - 1));
		//cout << "rs.generator_M_systematic_any_pos" << rs.generator_M_systematic_any_pos;
		//cout << "rs.information_set_ind" << rs.information_set_ind;
		// the non-systematic part of rs generator matrix, and transposed, n * r_prime matrix
		Matrix<ty> parity_part = rs.generator_M_systematic_any_pos.get_cols(rs.redundancy_set_ind).Transpose();
		//cout << "parity_part" << parity_part;
		// extend GF2e<m> to GF2 and apply to each element of parity_part, generate a (m*r_prime) * n matrix
		Matrix<GF2> extend_parity_part = ty::to_bits_row_extention(parity_part);
		//cout << "extend_parity_part" << extend_parity_part;
		// throw away the rows having ending 1 at last r_prime columns
		Matrix<int> erase_row_ind(1, r_prime);
		for (int i = 0; i < r_prime; ++i) {
			erase_row_ind(i) = i * m;
		}
		//cout << "erase_row_ind" << erase_row_ind;
		Matrix<GF2> sub_parity_matrix = extend_parity_part.erase_rows(erase_row_ind);
		//cout << "sub_parity_matrix" << sub_parity_matrix;
		// randomly choose (k_prime - k) line, we directly choose the first (k_prime - k)
		Matrix<int> row_ind(1, sub_parity_matrix.row(), 'N');
		Matrix<int> chosen_row_ind = row_ind.get_random_element(k_prime - k);
		//cout << "chosen_row_ind" << chosen_row_ind;

		// by size of (k_prime - k) * k_prime
		Matrix<GF2> chosen_parity_check_matrix = sub_parity_matrix.get_rows(chosen_row_ind);
		//cout << "choosen_parity_check_matrix" << chosen_parity_check_matrix;
		Viterbi_unordered_map vit(chosen_parity_check_matrix);
		Matrix<my_double> corresponding_r = r.get_cols(rs.information_set_ind);
		//cout << "rs.information_set_ind" << rs.information_set_ind;
		//cout << "corresponding_r" << corresponding_r;

		int list_size = 6;
		Matrix<GF2> viterbi_list = vit.decode_v_once(corresponding_r, list_size);		// 25M memory, 337k operation
		//cout << "viterbi_list" << viterbi_list;
		viterbi_list = viterbi_list.get_row(0);

		//cout << "extend_parity_part.get_rows(erase_row_ind).Transpose()" << extend_parity_part.get_rows(erase_row_ind).Transpose();
		Matrix<GF2> Zr = viterbi_list * extend_parity_part.get_rows(erase_row_ind).Transpose();
		//cout << "Zr" << Zr;
		//cout << "rs.redundancy_set_ind" << rs.redundancy_set_ind;
		Matrix<GF2> v_hat = viterbi_list.insert_cols(Zr, rs.redundancy_set_ind);
		cout << "v_hat" << v_hat;

		// problematic
		//Matrix<int> permuted_r_ind = rs.information_set_ind.combine_right(rs.redundancy_set_ind);
		//v_hat.permute_back(permuted_r_ind);
		//cout << "v_hat" << v_hat;

		cout << "v_hat-v" << v_hat - v;

	}
	static void repeatedly_test_BCH_non_Gaussian_elimination() {
		for (int i = 0; i < 10; ++i) {
			BCH_non_Gaussian_elimination();
		}
	}
	static void BCH_non_Gaussian_elimination_simulation(my_double SNR_dB = 2.0) {

		const int m = 6;
		const int t = 3;
		const int n = (1 << m) - 1;
		const int d = 2 * t + 1;
		const int k_prime = n - d + 1;
		const int r_prime = d - 1;
		const int simulation_times = 10;

		/* decoding parameter */
		const int selected_row_num = 7;
		const int max_list_num = 128;

		typedef GF2e<m> ty;
		ty::init();

		BCH<m, t> bch;
		int k = bch.get_k();

		bch.print_info();
		cout << "selected_row_num = " << selected_row_num << endl;
		cout << "max_list_num = " << max_list_num << endl;

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//cout << "v" << v;
		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoder
		RS<m, k_prime> rs;
		int error_frame = 0;
		int num_invalid_list = 0;
		double viterbi_time_ave = 0;

		Viterbi_unordered_map vit(k_prime, selected_row_num);
		clock_t start, end;
		//clock_t	viterbi_start, viterbi_end;
		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> r = AWGN::pass_standard(c);

			// following states the procees of decoding

			//cout << "r" << r;
			Matrix<my_double> r_abs = r.get_abs();

			Matrix<int> permute_ind = r_abs.sort_with_ind('>');		// from big to small
			//cout << "permute_ind" << permute_ind;

			// rs generator matrix from most reliable information set of bch received codeword
			rs.generate_systematic_generator_any_pos(permute_ind.get_part(0, 0, 0, k_prime - 1));
			// this is not optimal, change to 'generate_systematic_generator_any_pos_best_r' if you want

			//cout << "rs.generator_M_systematic_any_pos" << rs.generator_M_systematic_any_pos;
			//cout << "rs.information_set_ind" << rs.information_set_ind;
			// the non-systematic part of rs generator matrix, and transposed, n * r_prime matrix
			Matrix<ty> parity_part = rs.generator_M_systematic_any_pos.get_cols(rs.redundancy_set_ind).Transpose();
			//cout << "parity_part" << parity_part;
			// extend GF2e<m> to GF2 and apply to each element of parity_part, generate a (m*r_prime) * n matrix
			Matrix<GF2> extend_parity_part = ty::to_bits_row_extention(parity_part);
			//cout << "extend_parity_part" << extend_parity_part;
			// throw away the rows having ending 1 at last r_prime columns
			Matrix<int> erase_row_ind(1, r_prime);
			for (int i = 0; i < r_prime; ++i) {
				erase_row_ind(i) = i * m;
			}
			//cout << "erase_row_ind" << erase_row_ind;
			Matrix<GF2> sub_parity_matrix = extend_parity_part.erase_rows(erase_row_ind);
			//cout << "sub_parity_matrix" << sub_parity_matrix;
			// randomly choose (k_prime - k) line
			Matrix<int> row_ind(1, sub_parity_matrix.row(), 'N');
			Matrix<int> chosen_row_ind = row_ind.get_random_element(selected_row_num);
			//cout << "chosen_row_ind" << chosen_row_ind;

			// by size of (k_prime - k) * k_prime
			Matrix<GF2> chosen_parity_check_matrix = sub_parity_matrix.get_rows(chosen_row_ind);
			//cout << "choosen_parity_check_matrix" << chosen_parity_check_matrix;

			vit.change_PM(chosen_parity_check_matrix);

			// select the valid codeword from list

			//Matrix<GF2> backup_not_processed_parity_check_matrix = sub_parity_matrix.erase_rows(chosen_row_ind);
			Matrix<GF2> not_processed_parity_check_matrix = sub_parity_matrix.erase_rows(chosen_row_ind);

			/*if (backup_not_processed_parity_check_matrix != not_processed_parity_check_matrix) {
				cout << "sub_parity_matrix" << sub_parity_matrix;

				cout << "backup_not_processed_parity_check_matrix" << backup_not_processed_parity_check_matrix;
				cout << "not_processed_parity_check_matrix" << not_processed_parity_check_matrix;

				cout << "error: backup_not_processed_parity_check_matrix - not_processed_parity_check_matrix"
					<< backup_not_processed_parity_check_matrix - not_processed_parity_check_matrix;

				Matrix<GF2> a_revisit = sub_parity_matrix.erase_rows(chosen_row_ind);

				cout << "a_revisit" << a_revisit;
				cout << "a_revisit - backup_not_processed_parity_check_matrix" << a_revisit - backup_not_processed_parity_check_matrix;
				cout << "a_revisit - not_processed_parity_check_matrix" << a_revisit - not_processed_parity_check_matrix;
			}*/

			// compute the inner product of vector in viterbi_list and each row of not_processed_parity_check_matrix
			// if one result is non zero, then the codeword is invalid and should be discarded

			vit.change_unused_PM(not_processed_parity_check_matrix);

			Matrix<my_double> corresponding_r = r.get_cols(rs.information_set_ind);
			//cout << "rs.information_set_ind" << rs.information_set_ind;
			//cout << "corresponding_r" << corresponding_r;

			//viterbi_start = clock();
			Matrix<GF2> viterbi_list = vit.decode_v_once(corresponding_r, max_list_num);		// 25M memory, 337k operation			
			//viterbi_end = clock();
			//viterbi_time_ave += (double)viterbi_end - viterbi_start;

			//cout << "viterbi_list" << viterbi_list;
			Matrix<GF2> v_hat;
			Matrix<GF2> extend_generator = extend_parity_part.get_rows(erase_row_ind).Transpose();

			int best_ind_list = -1;
			my_double min_soft_distance = 3000;

			if (viterbi_list.size() != 0) {
				// we should not directly choose the fist row
				Matrix<GF2> can_v_hat;

				int list_size = viterbi_list.row();
				for (int i = 0; i < list_size; ++i) {
					//cout << "extend_parity_part.get_rows(erase_row_ind).Transpose()" 
					//	<< extend_parity_part.get_rows(erase_row_ind).Transpose();
					can_v_hat = viterbi_list.get_row(i);
					//cout << "can_v_hat(1)" << can_v_hat << endl;
					Matrix<GF2> Zr = can_v_hat * extend_generator;
					//cout << "Zr" << Zr;
					//cout << "rs.redundancy_set_ind" << rs.redundancy_set_ind;
					can_v_hat = can_v_hat.insert_cols(Zr, rs.redundancy_set_ind);
					//cout << "can_v_hat(2)" << can_v_hat << endl;
					//cout << "v_hat" << v_hat;

					my_double can_soft_distance = Measure::correlation_discrepancy_v(r, can_v_hat);
					bool take_can = can_soft_distance < min_soft_distance;
					min_soft_distance = take_can ? can_soft_distance : min_soft_distance;
					best_ind_list = take_can ? i : best_ind_list;
					v_hat = take_can ? can_v_hat : v_hat;		//we can also compute the list error probability, which is smaller than ML
				}

				/*if (best_ind_list != 0) {
					cout << "(not 0) best_ind_list / selected_ind_list = " << best_ind_list << " / " << viterbi_list.row() << endl;
				}*/
			}
			else {
				// there is no valid list, decode fail
				//cout << "there is no valid list, i=" << i << endl;
				v_hat = Matrix<GF2>(1, n, '0');
				num_invalid_list++;
			}

			//cout << "v_hat" << v_hat;

			//cout << "v_hat-v" << v_hat - v;

			if (v == v_hat);
			else {
				/*bool test = bch.is_in_v_space(v_hat);
				cout << "error frame: i=" << i << "\tlog2_max_state_num = " << vit.log2_max_state_num << "\t"
					<< "best_ind_list / selected_ind_list = " << best_ind_list << " / " << viterbi_list.row() - 1
					<< "\nbch.is_in_v_space(v_hat) = " << test << "\tmin_soft_distance = "<< min_soft_distance << endl;*/

					//// problem!
					//cout << "sub_parity_matrix" << sub_parity_matrix;
					//cout << "chosen_row_ind" << chosen_row_ind;
					//cout << "chosen_parity_check_matrix" << chosen_parity_check_matrix;
					//cout << "not_processed_parity_check_matrix" << not_processed_parity_check_matrix;
					//cout << "sub_parity_matrix.rank() = " << sub_parity_matrix.rank() << endl;
					//cout << "chosen_parity_check_matrix.rank() = " << chosen_parity_check_matrix.rank() << endl;
					//cout << "not_processed_parity_check_matrix.rank() = " << not_processed_parity_check_matrix.rank() << endl;
					//cout << "backup_not_processed_parity_check_matrix.rank() = " << backup_not_processed_parity_check_matrix.rank() << endl;
					//cout << "backup_not_processed_parity_check_matrix - not_processed_parity_check_matrix" 
					//	<< backup_not_processed_parity_check_matrix - not_processed_parity_check_matrix;
					//cout << "\n-------\n";
					//cout << "test_zero" << vit.test_zero;
				error_frame++;
			}

			/*if (i % 100 == 0) {
				cout << "\n-------------\n";
				cout << "\n-----------i=" << i << "-----------\n";
				cout << "\n-------------\n";
			}*/
		}

		end = clock();
		vit.print_rank_distribution();
		cout << "------" << endl;


		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		//viterbi_time_ave /= CLOCKS_PER_SEC * (double)simulation_times;
		cout << fixed << setprecision(2);
		//cout << "viterbi_time_ave = " << viterbi_time_ave * 1000 << "\tms/iteration" << endl;
		//cout << "time_consume = " << time_consume * 1000 << "\tms/iteration;" << endl;

		cout.unsetf(ios_base::fixed);
		cout << setprecision(6);

		cout << "num_invalid_list = " << num_invalid_list << endl;
		cout << "error_frame = " << error_frame << endl;

		cout << scientific << setprecision(2);
		cout << "error_rate = " << error_frame / (my_double)simulation_times << endl;
		cout.unsetf(ios::floatfield);

	}
	static void BCH_non_Gaussian_elimination_simulation2(my_double SNR_dB = 2.0) {

		cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;

		//cout << "GF2e_auxiliary_storage::operation_number = " << GF2e_auxiliary_storage::operation_number << endl;

		const int m = 3;
		const int t = 1;
		const int n = (1 << m) - 1;
		const int d = 2 * t + 1;
		const int k_prime = n - d + 1;
		const int r_prime = d - 1;
		int max_simulation_times = 1000;			// 50000000
		const int max_error_num = 200;

		/* decoding parameter */
		const int selected_row_num = 2;
		const int max_list_num = 4;
		//const my_double beta = 3;			// for stopping rule

		typedef GF2e<m> ty;
		ty::init();

		//cout << "GF2e_auxiliary_storage::operation_number = " << GF2e_auxiliary_storage::operation_number << endl;
		BCH<m, t> bch;
		int k = bch.get_k();

		//cout << "GF2e_auxiliary_storage::operation_number = " << GF2e_auxiliary_storage::operation_number << endl;

		bch.print_info();
		cout << "selected_row_num = " << selected_row_num << endl;
		cout << "max_list_num = " << max_list_num << endl;
		//cout << "beta = " << beta << endl;

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//cout << "v" << v;
		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoder
		RS<m, k_prime> rs;		// mother code
		LL_OSD_Viterbi<m, t, k_prime, selected_row_num> LL_vit;
		int error_frame = 0;
		double viterbi_time_ave = 0;

		//cout << "GF2e_auxiliary_storage::operation_number = " << GF2e_auxiliary_storage::operation_number << endl;

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

		clock_t start, end;
		//clock_t	viterbi_start, viterbi_end;
		start = clock();

		for (int i = 0; i < max_simulation_times; ++i) {
			Matrix<my_double> r = AWGN::pass_standard(c);

			// following states the procees of decoding

			//cout << "GF2e_auxiliary_storage::operation_number: before = " << GF2e_auxiliary_storage::operation_number << endl;
			Matrix<GF2> v_hat = LL_vit.decode_v(r, max_list_num);
			//cout << "GF2e_auxiliary_storage::operation_number: end = " << GF2e_auxiliary_storage::operation_number << endl;

			if (v == v_hat);
			else {
				error_frame++;
			}

			if (error_frame == max_error_num) {
				max_simulation_times = i + 1;
				break;
			}

		}
		end = clock();

#ifdef count_operation_number
#ifdef use_my_double

		cout << "------" << endl;
		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(LL_OSD_Viterbi) double_ope_num = "  \
			<< (double_ope_num_after - double_ope_num_before) / (double)max_simulation_times << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(LL_OSD_Viterbi) GF2_ope_num = "  \
			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)max_simulation_times << endl;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << "(LL_OSD_Viterbi) GF2e_ope_num = "  \
			<< (GF2e_ope_num_after - GF2e_ope_num_before) / (double)max_simulation_times << endl;
#endif // count_operation_number

		cout << "------" << endl;
		cout << "max_simulation_times = " << max_simulation_times << endl;
		cout << "list_ave = " << LL_vit.total_used_list_num / (double)max_simulation_times << endl;

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)max_simulation_times;
		cout << scientific << setprecision(2);
		cout << "rate of num_invalid_list = " << (double)LL_vit.num_invalid_list / max_simulation_times << endl;
		cout << "error_rate = " << error_frame / (double)max_simulation_times << endl;
		cout << "time_consume = " << time_consume << "\ts/iteration;" << endl;
		cout.unsetf(ios::floatfield);
		cout << setprecision(6);

	}
	static void BCH_non_Gaussian_elimination_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[155] = { 0 };
			sprintf_s(file_name, 155, "BCH_non_Gaussian_elimination_simulation2.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(0, 0.5, 0, 'd');
		cout << "test_SNR" << test_SNR;
		int len = test_SNR.size();
		for (int i = 0; i < len; ++i) {
			BCH_non_Gaussian_elimination_simulation2(test_SNR(i));
		}
	}
	static void viterbi3_test_any_PM() {
		const int n = 16;
		const int r = 10;		// good

		Viterbi_unordered_map vit(n, r);
		for (int w = 0; w < 100; ++w) {
			// recursively test

			Matrix<GF2> PM(r, n, 'b');
			//cout << "PM" << PM;
			int rank_PM = PM.rank();
			//cout << "rank_PM = " << rank_PM << endl;

			Matrix<double> recv(1, n, my::rand_ga);
			//cout << "recv" << recv;

			vit.change_PM(PM);

			Matrix<GF2> decode_list = vit.decode_v_once(recv, 1 << (n - rank_PM));
			//cout << "decode_list" << decode_list;
			cout << "------------(w=" << w << ")ending-------------" << endl;
		}

		vit.print_rank_distribution();
	}
	static void erase_rows_test() {
		for (int i = 0; i < 10000; ++i) {
			Matrix<GF2> big(30, 57, 'b');
			Matrix<int> row_ind(1, 30, 'N');
			row_ind = row_ind.get_random_element(12);

			//cout << "big" << big;
			//cout << "row_ind" << row_ind;

			Matrix<GF2> ebig = big.erase_rows(row_ind);
			//cout << "ebig" << ebig;

			Matrix<GF2> e2big = big.erase_rows(row_ind);

			if (ebig != e2big) {
				cout << "big.rank()=" << big.rank() << endl;
				cout << "ebig.rank()=" << ebig.rank() << endl;
				cout << "error_occur" << endl;
				cout << "ebig" << ebig;
				cout << "e2big" << e2big;
				cout << "ebig - e2big" << ebig - e2big;
			}
			else {
				cout << "-------" << i << "-------" << endl;
			}
		}

	}
	static void BCH_test_viterbi3_repeatedly_time() {
		cout << boolalpha;

		const int m = 7;
		const int t = 2;
		typedef GF2e<m> ty;
		ty::init();

		BCH<m, t> bch;

		bch.print_info();		// BCH(127,113,5)
		int n = bch.get_n();
		int k = bch.get_k();
		int d = bch.get_d();
		int r = bch.get_r();

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//v.permute(permute_test);

		Matrix<my_double> c = BPSK::modulation(v);
		Matrix<my_double> e(1, n);
		for (int i = 0; i < n; ++i) {
			e(i) = my::rand_ga() * 0.6;		// in this case v_hat is error but v_hat_viterbi is correct
		}

		Matrix<my_double> recv = c + e;
		Matrix<GF2> hdr = BPSK::demodulation(recv);

		Matrix<GF2> v_hat = bch.decode_v(hdr);

		clock_t start, end;
		int simulation_times = 400;

		//Viterbi_cyclic_code_r vit_bch(bch.get_parity_matrix());		// (0.039s, 31MB)
		//start = clock();
		//for (int i = 0; i < simulation_times; ++i) {
		//	Matrix<GF2> list_v_hat_viterbi = vit_bch.decode_v(recv, 4096);
		//}
		//end = clock();
		//double time_consume_Viterbi = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		//cout << "time_consume_Viterbi = " << time_consume_Viterbi << endl;

		Viterbi_optimized vit_bch(bch.get_parity_matrix());	// the fastest and least memory
		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			//Matrix<GF2> list_v_hat_viterbi = vit_bch.decode_v(recv, 4096, 1);	// (0.011s, 17MB)
			//Matrix<GF2> list_v_hat_viterbi = vit_bch.decode_v(recv, 4096, 4096);	// (0.038s, 30MB)
			//Matrix<GF2> list_v_hat_viterbi = vit_bch.decode_v2(recv, 4096, 4096);	// (0.035s, 19MB)
			//Matrix<GF2> list_v_hat_viterbi = vit_bch.decode_v3(recv, 4096, 4096);	// (0.031s, 17MB)
			//Matrix<GF2> list_v_hat_viterbi = vit_bch.decode_v4(recv, 4096, 4096);	// (0.041s, 19MB)
			//Matrix<GF2> list_v_hat_viterbi = vit_bch.decode_v5(recv, 4096, 4096);	// (0.041s, 19MB)
			Matrix<GF2> list_v_hat_viterbi = vit_bch.decode_v(recv, 4096, 4096);	// (0.027s, 17MB), the best
		}
		end = clock();
		double time_consume_Viterbi = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		cout << "time_consume_Viterbi = " << time_consume_Viterbi << endl;

		//Viterbi_unordered_map vit_bch3(bch.get_parity_matrix());	// (0.071s, 120MB)
		//vit_bch3.initialize_state();
		//start = clock();
		//for (int i = 0; i < simulation_times; ++i) {
		//	//Matrix<GF2> list_v_hat_viterbi_test = vit_bch3.decode_v_once(recv, 16);
		//	Matrix<GF2> list_v_hat_viterbi_tmp = vit_bch3.decode_v_repeatedly(recv, 16);
		//	//Matrix<GF2> list_v_hat_viterbi = vit_bch3.decode_v_repeatedly(recv, 4);		// the second time will be faster
		//}
		//end = clock();
		//double time_consume_Viterbi3 = ((double)end - start) / CLOCKS_PER_SEC / (my_double)simulation_times;
		//cout << "time_consume_Viterbi3 = " << time_consume_Viterbi3 << endl;

		//Viterbi_Sorted_vector vit_bch4(bch.get_parity_matrix());			// (0.123s, 40MB)
		//vit_bch4.initialize_state();
		//start = clock();
		//for (int i = 0; i < simulation_times; ++i) {
		//	Matrix<GF2> list_v_hat_viterbi_tmp = vit_bch4.decode_v_repeatedly_space_saving(recv, 16);
		//}
		//end = clock();
		//double time_consume_Viterbi4 = ((double)end - start) / CLOCKS_PER_SEC / (my_double)simulation_times;
		//cout << "time_consume_Viterbi4 = " << time_consume_Viterbi4 << endl;
	}
	static void Sorted_vector_find_test() {
		// it also provides an example to use sorted vector
		Sorted_vector<double> sv(7);
		sv.assign_ind(0, 2);
		sv.val(0) = 0.2;
		sv.assign_ind(1, 3);
		sv.val(1) = 0.3;
		sv.assign_ind(2, 4);
		sv.val(2) = 0.4;
		sv.assign_ind(3, 6);
		sv.val(3) = 0.6;
		sv.assign_ind(4, 8);
		sv.val(4) = 0.8;
		sv.assign_ind(5, 11);
		sv.val(5) = 0.11;
		sv.assign_ind(6, 15);
		sv.val(6) = 0.15;

		cout << "sv" << sv;
		cout << sv.find(2) << endl;
		cout << sv.find(3) << endl;
		cout << sv.find(4) << endl;
		cout << sv.find(5) << endl;
		cout << sv.find(6) << endl;
		cout << sv.find(8) << endl;
		cout << sv.find(11) << endl;
		cout << sv.find(15) << endl;

		cout << sv[2] << endl;
		cout << sv[3] << endl;
		cout << sv[4] << endl;
		//cout << sv[5] << endl;		// non-exist, output 0.6
		cout << sv[6] << endl;
		cout << sv[8] << endl;
		cout << sv[11] << endl;
		cout << sv[15] << endl;
	}
	static void OSD_new_test() {
		const int m = 6;
		const int t = 3;
		typedef GF2e<m> ty;
		ty::init();
		BCH<m, t> bch;
		bch.print_info();
		int k = bch.get_k();
		int n = bch.get_n();

		Matrix<GF2> u(1, k, 'b');
		cout << "u" << u;
		Matrix<GF2> v = bch.encode(u);
		cout << "v" << v;
		Matrix<my_double> c = BPSK::modulation(v);
		cout << "c" << c;
		Matrix<my_double> r = AWGN::pass(c, 0.7);
		cout << "r" << r;
		Matrix<GF2> hdr = BPSK::demodulation(r);
		cout << "hdr" << hdr;
		cout << "e" << hdr - v;

		Matrix<GF2> du = bch.decode(hdr);
		cout << "du" << du;
		cout << "du - u" << du - u;

		OSD_r osd(bch.get_parity_matrix(), bch.get_d());
		Matrix<GF2> dv = osd.decode_v(r, 1);
		cout << "dv - v" << dv - v;

		Matrix<GF2> dv2 = OSD::decode_v(bch, r, 1);
		cout << "dv2 - v" << dv2 - v;

		cout << "dv2 - dv" << dv2 - dv;
	}
	static void OSD_encode_test() {
		const int m = 6;
		const int t = 2;
		typedef GF2e<m> ty;
		ty::init();
		BCH<m, t> bch;
		bch.print_info();
		int k = bch.get_k();
		int n = bch.get_n();

		Matrix<GF2> u(1, k, 'b');		// note that the random seed change with k
		cout << "u" << u;

		Matrix<GF2> v = bch.encode(u);
		cout << "v" << v;
		Matrix<my_double> c = BPSK::modulation(v);
		cout << "c" << c;
		Matrix<my_double> r = AWGN::pass(c, 1);
		cout << "r" << r;

		Matrix<GF2> hdr = BPSK::demodulation(r);
		cout << "hdr" << hdr;

		cout << "----------- osd_new -------" << endl;
		OSD_r osd(bch.get_parity_matrix(), bch.get_d());
		Matrix<GF2> dv = osd.decode_v(r, 1);

		cout << "----------- osd_old -------" << endl;
		Matrix<GF2> dv2 = OSD::decode_v(bch, r, 1);

		my_double dist_v = Measure::correlation_discrepancy_v(r, v);
		my_double dist_dv = Measure::correlation_discrepancy_v(r, dv);
		my_double dist_dv2 = Measure::correlation_discrepancy_v(r, dv2);
		cout << "dist_v = " << dist_v << endl;
		cout << "dist_dv = " << dist_dv << endl;		// if this is less, ML decoding error occur
		cout << "dist_dv2 = " << dist_dv2 << endl;
	}
	static void OSD_new_simulation() {
		const int m = 8;
		const int t = 11;
		typedef GF2e<m> ty;
		ty::init();
		BCH<m, t> bch;
		bch.print_info();
		int k = bch.get_k();
		int n = bch.get_n();

		Matrix<GF2> u(1, k, 'b');		// note that the random seed change with k
		Matrix<GF2> v = bch.encode(u);
		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR_dB = 2.0;
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		const int simulation_times = 300;
		int error_frame = 0;
		OSD_r osd(bch.get_parity_matrix(), bch.get_d());

		clock_t start, end;
		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> r = AWGN::pass(c);

			cout << "i=" << i << endl;
			Matrix<GF2> dv = osd.decode_v(r, 3);
			if (v == dv);
			else {
				error_frame++;
			}
		}

		end = clock();
		cout << "------" << endl;
		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		cout << fixed << setprecision(2);
		cout << "time_consume = " << time_consume * 1000 << "\tms/iteration;" << endl;
		cout << "error_rate = " << error_frame / (my_double)simulation_times << endl;
		cout.unsetf(ios_base::fixed);
		cout << setprecision(6);
		cout << "error_frame = " << error_frame << endl;
	}
	static void viterbi_preprocessing_test() {
		const int m = 5;
		const int t = 2;
		typedef GF2e<m> ty;
		ty::init();
		BCH<m, t> bch;
		Matrix<GF2> PM = bch.get_parity_matrix();
		Matrix<int> natual(1, (1 << m) - 1, 'N');
		Matrix<int> rand_permutation = natual.get_random_element(31);
		cout << "rand_permutation" << rand_permutation;
		PM.permute_col(rand_permutation);

		//Matrix<GF2> PM(6, 17, 'b');
		int r = PM.row();
		int n = PM.col();

		cout << "PM" << PM;
		PM.row_transformation_to_low_triangle();
		PM.col_permute_to_full_rank_on_right();
		cout << "PM: full rank on right" << PM;		// if the right part of the Matrix is rank defficient, error occur
		PM.row_transformation_right_low_triangle_to_identity();
		cout << "PM: low triangle" << PM;

		PM.row_transformation_to_up_triangle();
		PM.col_permute_to_full_rank_on_left();
		cout << "PM: up triangle" << PM;

		Matrix<int> right_part_permute_ind(1, n, 'N');
		// arange the column of the right part according to the first 1 of the column, 
		// i.e., {0,0,1,0,0,1} will be arranged to the 3rd column of the right part
		for (int j = n - r; j < n; ++j) {
			// find the starting 1 of column j
			int i = 0;		// index of the starting 1
			for (; i < r; ++i) {
				if (PM(i, j) == 0);
				else {
					break;		// it must break at some point
				}
			}

			// set permute ind
			right_part_permute_ind(n - r + i) = j;
			/*for (int k = 0; k < r; ++k) {
				right_part(k, i) = PM(k, j);
			}*/
		}

		// equal to permute column at right part of PM
		//PM.set_part(0, n - r, right_part);
		PM.permute_col(right_part_permute_ind);
		cout << "PM: switching right part" << PM;

		// continuing permute the columns of Parity check Matrix, it seems meaningless, find a book to read
		Matrix<int> middle_part_permute_ind(1, n, 'N');
		// permute the middle n-2r column
		for (int j = 0; j < n; ++j) {
			// get the starting 1 and ending 1 of the column

			// find the starting 1 of column j
			int i = 0;			// index of the starting 1
			for (; i < r; ++i) {
				if (PM(i, j) == 0);
				else {
					break;		// it must break at some point
				}
			}
			// find the ending 1 of column j
			int k = r - 1;		// index of the ending 1
			for (; k >= 0; --k) {
				if (PM(k, j) == 0);
				else {
					break;		// it must break at some point
				}
			}

			if (i < r - 1 - k) {
				// permute the column to the front, deciding parameter r-1-k
				middle_part_permute_ind(j) = k;
			}
			else {
				// permute the column to the back, deciding parameter i
				middle_part_permute_ind(j) = n - r + i;
			}
		}
		cout << "middle_part_permute_ind" << middle_part_permute_ind;
		Matrix<int> sort_permute = middle_part_permute_ind.sort_with_ind('<');
		PM.permute_col(sort_permute);
		cout << "PM: switching middle part" << PM;
	}
	static bool viterbi_preprocessing_test2(int& min_computation) {
		const int m = 5;
		const int t = 3;
		typedef GF2e<m> ty;
		ty::init();
		BCH<m, t> bch;
		Matrix<GF2> PM = bch.get_parity_matrix();
		Matrix<int> natual(1, (1 << m) - 1, 'N');
		Matrix<int> rand_permutation = natual.get_random_element((1 << m) - 1);
		//cout << "rand_permutation" << rand_permutation;
		PM.permute_col(rand_permutation);

		int r = PM.row();
		int n = PM.col();

		//cout << "PM" << PM;
		PM.row_transformation_to_low_triangle();
		Matrix<int> permutation_record_2 = PM.col_permute_to_full_rank_on_right();
		//cout << "PM: full rank on right" << PM;		// if the right part of the Matrix is rank defficient, error occur
		PM.row_transformation_right_low_triangle_to_identity();
		//cout << "PM: low triangle" << PM;

		PM.row_transformation_to_up_triangle();

		//PM.col_permute_to_full_rank_on_left();
		//cout << "PM: up triangle" << PM;

		//// continuing permute the columns of Parity check Matrix, it seems meaningless, find a book to read
		//Matrix<int> permute_ind(1, n, 'N');
		//// permute the middle n-2r column
		//for (int j = 0; j < n; ++j) {
		//	// get the starting 1 and ending 1 of the column
		//	// find the starting 1 of column j
		//	int i = 0;			// index of the starting 1
		//	for (; i < r; ++i) {
		//		if (PM(i, j) == 0);
		//		else {
		//			break;		// it must break at some point
		//		}
		//	}
		//	// find the ending 1 of column j
		//	int k = r - 1;		// index of the ending 1
		//	for (; k >= 0; --k) {
		//		if (PM(k, j) == 0);
		//		else {
		//			break;		// it must break at some point
		//		}
		//	}
		//	if (i < r - 1 - k) {
		//		// permute the column to the front, deciding parameter r-1-k
		//		permute_ind(j) = k;
		//	}
		//	else {
		//		// permute the column to the back, deciding parameter i
		//		permute_ind(j) = n - r + i;
		//	}
		//}
		////cout << "middle_part_permute_ind" << permute_ind;
		//Matrix<int> sort_permute = permute_ind.sort_with_ind('<');
		//PM.permute_col(sort_permute);
		////cout << "PM: switching middle part" << PM;

		// continuing permute the columns of Parity check Matrix, it seems meaningless, find a book to read
		Matrix<int> permute_ind(1, n, 'N');
		// permute the middle n-2r column
		for (int j = 0; j < n; ++j) {
			// get the starting 1 and ending 1 of the column
			// find the starting 1 of column j
			int i = 0;			// index of the starting 1
			for (; i < r; ++i) {
				if (PM(i, j) == 0);
				else {
					break;		// it must break at some point
				}
			}
			// find the ending 1 of column j
			int k = r - 1;		// index of the ending 1
			for (; k >= 0; --k) {
				if (PM(k, j) == 0);
				else {
					break;		// it must break at some point
				}
			}
			if (i < r - 1 - k) {
				// permute the column to the front, deciding parameter r-1-k
				permute_ind(j) = k * r + i;
			}
			else {
				// permute the column to the back, deciding parameter i
				permute_ind(j) = (n - r + i) * r + k;
			}
		}
		//cout << "middle_part_permute_ind" << permute_ind;
		Matrix<int> sort_permute = permute_ind.sort_with_ind('<');
		PM.permute_col(sort_permute);

		// record all permutation
		natual.permute(rand_permutation);
		natual.permute(permutation_record_2);
		natual.permute(sort_permute);


		//cout << "PM" << PM;

		// count saving states, note that 1 may not continue at border
		Matrix<int> states_num(1, n);

		int ending_1 = 0;
		int time_instance = 0;
		for (; time_instance < n && ending_1 != r - 1; ++time_instance) {
			ending_1 = r - 1;		// index of the ending 1 of a column
			for (; ending_1 >= 0; --ending_1) {
				if (PM(ending_1, time_instance) == 0);
				else {
					break;
				}
			}
			states_num(time_instance) = 1 << (ending_1 + 1);
		}
		for (; time_instance < n; ++time_instance) {
			states_num(time_instance) = 1 << r;
		}

		states_num(n - 1) = 1;
		int starting_1 = r - 1;
		time_instance--;
		for (; time_instance >= 0 && starting_1 != 0; --time_instance) {
			starting_1 = 0;			// index of the starting 1 of a column
			for (; starting_1 < r; ++starting_1) {
				if (PM(starting_1, time_instance) == 0);
				else {
					break;
				}
			}
			states_num(time_instance - 1) >>= starting_1;
		}
		//cout << "states_num" << states_num;

		//bool starting = true;
		//bool ending = false;
		//for (int j = 0; j < n; ++j) {
		//	// find the starting 1 of column j, but starting and ending may overlap
		//	if (starting) {
		//		// find the ending 1 of column j
		//		int k = r - 1;		// index of the ending 1
		//		for (; k >= 0; --k) {
		//			if (PM(k, j) == 0);
		//			else {
		//				break;		// it must break at some point
		//			}
		//		}
		//		states_num(j) = 1 << (1 + k);
		//		starting = k != r - 1;
		//	}	
		//	if (starting == false && ending == false) {
		//		states_num(j) = 1 << r;
		//		int i = 0;			// index of the starting 1
		//		for (; i < r; ++i) {
		//			if (PM(i, j) == 0);
		//			else {
		//				break;		// it must break at some point
		//			}
		//		}
		//		ending = i != 0;
		//	}
		//	if (ending) {
		//		int i = 0;			// index of the starting 1
		//		for (; i < r; ++i) {
		//			if (PM(i, j) == 0);
		//			else {
		//				break;		// it must break at some point
		//			}
		//		}
		//		states_num(j - 1) = 1 << (r - i);
		//	}
		//	
		//	// find the ending 1 of column j			
		//}
		//states_num(n - 1) = 1;

		// compute total computation
		int total_compare_computation = 0;
		int total_add_computation = 0;
		total_add_computation = 2;
		for (int j = 1; j < n; ++j) {
			if (states_num(j) == 2 * states_num(j - 1)) {
				total_add_computation += states_num(j - 1) * 2;		// one trellis for each state
			}
			else if (states_num(j) == states_num(j - 1) || states_num(j) == states_num(j - 1) / 2) {
				total_add_computation += states_num(j) * 2;			// 2 trellis for each state
				total_compare_computation += states_num(j);
			}
			else if (states_num(j) > states_num(j - 1)) {
				total_add_computation += states_num(j - 1) * 2;
				//cout << "unexpected state, over expending. j = " << j << endl;
			}
			else {
				total_add_computation += states_num(j) * 2;			// 2 trellis for each state
				total_compare_computation += states_num(j);
				//cout << "unexpected state, over shrinking. j = " << j << endl;
			}
			// for optimized network, there should be no expending and shrinking
		}
		if (total_add_computation < min_computation) {
			cout << "rand_permutation" << rand_permutation;
			cout << "PM" << PM;
			cout << "states_num" << states_num;
			cout << "total_add_computation = " << total_add_computation << endl;
			min_computation = total_add_computation;

			// validation of PM
			cout << "the best permutation:" << natual;
			// permuting back
			PM.permute_col_back(natual);
			cout << "(bch.get_generator_matrix() * PM.Transpose()).isZero() = "
				<< (bch.get_generator_matrix() * PM.Transpose()).isZero() << endl;
			return true;
		}
		else {
			return false;
		}
		//cout << "total_compare_computation = " << total_compare_computation << endl;
		//int total_computation = total_add_computation + total_compare_computation;
		//cout << "total_computation = " << total_computation << endl;
	}
	static void viterbi_preprocessing_test_repeatedly() {
		// random test with viterbi_preprocessing_test2
		int min_computation = INT_MAX;
		for (int i = 0; i < 100000; ++i) {
			//cout << "i=" << i << endl;
			if (viterbi_preprocessing_test2(min_computation)) {
				cout << "i = " << i << endl;
			}
		}
	}
	static void viterbi_preprocessing_test_orig() {
		const int m = 5;
		const int t = 3;
		typedef GF2e<m> ty;
		ty::init();
		BCH<m, t> bch;
		bch.print_info();
		Matrix<GF2> PM = bch.get_parity_matrix();
		cout << "PM" << PM;

		int r = PM.row();
		int n = PM.col();
		Matrix<int> states_num(1, n);
		for (int i = 0; i < r; ++i) {
			states_num(i) = 1 << (i + 1);
		}
		for (int i = r; i < n - r; ++i) {
			states_num(i) = 1 << r;
		}
		for (int i = n - r; i < n; ++i) {
			states_num(i) = 1 << (n - 1 - i);
		}
		cout << "states_num" << states_num;

		// compute total computation
		int total_compare_computation = 0;
		int total_add_computation = 0;
		total_add_computation = 2;
		for (int j = 1; j < n; ++j) {
			if (states_num(j) == 2 * states_num(j - 1)) {
				total_add_computation += states_num(j - 1) * 2;		// one trellis for each state
			}
			else if (states_num(j) == states_num(j - 1) || states_num(j) == states_num(j - 1) / 2) {
				total_add_computation += states_num(j) * 2;			// 2 trellis for each state
				total_compare_computation += states_num(j);
			}
			else if (states_num(j) > states_num(j - 1)) {
				total_add_computation += states_num(j - 1) * 2;
				cout << "unexpected state, over expending. j = " << j << endl;
			}
			else {
				total_add_computation += states_num(j) * 2;			// 2 trellis for each state
				total_compare_computation += states_num(j);
				cout << "unexpected state, over shrinking. j = " << j << endl;
			}
			// for optimized network, there should be no expending and shrinking
		}
		cout << "total_add_computation = " << total_add_computation << endl;
		cout << "total_compare_computation = " << total_compare_computation << endl;
		int total_computation = total_add_computation + total_compare_computation;
		cout << "total_computation = " << total_computation << endl;
	}
	static void BCH_generate_test() {
		const int m = 10;
		const int t = 170;
		GF2e<m>::init();
		BCH<m, t> bch;
		bch.print_info();
	}
	static void viterbi_preprocessing_row_transform_test() {
		const int m = 6;
		const int t = 3;
		typedef GF2e<m> ty;
		ty::init();
		BCH<m, t> bch;
		int n = bch.get_n();
		int r = bch.get_r();
		Matrix<GF2> PM(18, 63, {
			1,1,0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,1,0,0,0,1,0,1,0,1,1,0,0,1,0,0,1,1,1,0,1,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,1,0,0,1,0,0,1,0,0,1,0,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,0,0,1,1,1,0,1,0,1,0,1,1,1,0,1,0,1,1,0,0,1,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,1,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,0,1,1,1,1,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,1,1,1,0,1,1,0,1,1,1,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,1,1,1,0,1,1,1,1,0,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,1,1,1,1,0,1,1,0,1,0,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,0,1,0,1,1,0,1,0,1,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,1,0,1,0,1,0,0,0,1,0,0,1,0,0,1,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0,0,1,0,1,1,1,1,0,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,0,1,0,0,1,1,1,0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1,0,1,1,0,0,1,1,0,0,1,1,1,1,1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,0,1,0,1,1,0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,0,1,1,1,0,0,1,0,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,1,1,1,1,0,1,1,1,1,1,0,0,0,1,0,0,1,1,1,1,1,0,1,1,0,0,1,0,0,0,1,1,1,0,1,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,1,1,1,0,1,0,1,0,1,0,1,1,1,1,1,0,0,1,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,1,0,1,0,1,0,0,1,0,0,0,0,1,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,0,1,0,1,0,1,0,0,1,0,0,0,0,1,1,1,1,1,0,1,1
			});

		// continuing permute the columns of Parity check Matrix, it seems meaningless, find a book to read
		Matrix<int> permute_ind(1, n, 'N');
		// permute the middle n-2r column
		for (int j = 0; j < n; ++j) {
			// get the starting 1 and ending 1 of the column
			// find the starting 1 of column j
			int i = 0;			// index of the starting 1
			for (; i < r; ++i) {
				if (PM(i, j) == 0);
				else {
					break;		// it must break at some point
				}
			}
			// find the ending 1 of column j
			int k = r - 1;		// index of the ending 1
			for (; k >= 0; --k) {
				if (PM(k, j) == 0);
				else {
					break;		// it must break at some point
				}
			}
			if (i < r - 1 - k) {
				// permute the column to the front, deciding parameter r-1-k
				permute_ind(j) = k * r + i;
			}
			else {
				// permute the column to the back, deciding parameter i
				permute_ind(j) = (n - r + i) * r + k;
			}
		}
		//cout << "middle_part_permute_ind" << permute_ind;
		Matrix<int> sort_permute = permute_ind.sort_with_ind('<');
		PM.permute_col(sort_permute);
		cout << "PM" << PM;

		// transform again, will be better
		int former_ending_1 = 0;
		for (int j = 1; j < n; ++j) {
			int ending_1 = r - 1;		// index of the ending 1 of a column
			for (; ending_1 >= 0; --ending_1) {
				if (PM(ending_1, j) == 0);
				else {
					break;
				}
			}

			// for normal case, ending_1 - former_ending_1==0 || ending_1 - former_ending_1== 1
			if (ending_1 - former_ending_1 > 1) {
				// find the main row to do transformation
				for (int i = ending_1; i > former_ending_1 + 1; --i) {
					GF2 ratio = PM(i, j) / PM(former_ending_1 + 1, j);		// Gaussian elemination to eliminate the row below
					for (int q = j; q < n; ++q) {
						PM(i, q) -= ratio * PM(former_ending_1 + 1, q);
					}
				}
				ending_1 = former_ending_1 + 1;
			}
			if (ending_1 == r - 1) {
				break;
			}

			former_ending_1 = ending_1;
		}

		cout << "PM" << PM;

		int later_starting_1 = r - 1;
		for (int j = n - 2; j >= 0; --j) {
			int starting_1 = 0;		// index of the starting 1 of a column
			for (; starting_1 < r; ++starting_1) {
				if (PM(starting_1, j) == 0);
				else {
					break;
				}
			}
			// for normal case, starting_1 == later_starting_1 || later_starting_1 - starting_1== 1
			if (later_starting_1 - starting_1 > 1) {
				// find the main row to do transformation
				for (int i = starting_1; i < later_starting_1 - 1; ++i) {
					GF2 ratio = PM(i, j) / PM(later_starting_1 - 1, j);		// Gaussian elemination to eliminate the row upon
					for (int q = j; q >= 0; --q) {
						PM(i, q) -= ratio * PM(later_starting_1 - 1, q);
					}
				}
				starting_1 = later_starting_1 - 1;
			}

			if (starting_1 == 0) {
				break;
			}
			later_starting_1 = starting_1;
		}

		//PM.row_transformation_to_up_triangle();
		//PM.row_transformation_to_low_triangle();

		// count saving states
		Matrix<int> states_num(1, n);

		int ending_1 = 0;
		int time_instance = 0;
		for (; time_instance < n && ending_1 != r - 1; ++time_instance) {
			ending_1 = r - 1;		// index of the ending 1 of a column
			for (; ending_1 >= 0; --ending_1) {
				if (PM(ending_1, time_instance) == 0);
				else {
					break;
				}
			}
			states_num(time_instance) = 1 << (ending_1 + 1);
		}
		for (; time_instance < n; ++time_instance) {
			states_num(time_instance) = 1 << r;
		}

		states_num(n - 1) = 1;
		int starting_1 = r - 1;
		time_instance--;
		for (; time_instance >= 0 && starting_1 != 0; --time_instance) {
			starting_1 = 0;			// index of the starting 1 of a column
			for (; starting_1 < r; ++starting_1) {
				if (PM(starting_1, time_instance) == 0);
				else {
					break;
				}
			}
			states_num(time_instance - 1) >>= starting_1;
		}


		// compute total computation
		int total_compare_computation = 0;
		int total_add_computation = 0;
		total_add_computation = 2;
		for (int j = 1; j < n; ++j) {
			if (states_num(j) == 2 * states_num(j - 1)) {
				total_add_computation += states_num(j - 1) * 2;		// one trellis for each state
			}
			else if (states_num(j) == states_num(j - 1) || states_num(j) == states_num(j - 1) / 2) {
				total_add_computation += states_num(j) * 2;			// 2 trellis for each state
				total_compare_computation += states_num(j);
			}
			else if (states_num(j) > states_num(j - 1)) {
				total_add_computation += states_num(j - 1) * 2;
				//cout << "unexpected state, over expending. j = " << j << endl;
			}
			else {
				total_add_computation += states_num(j) * 2;			// 2 trellis for each state
				total_compare_computation += states_num(j);
				//cout << "unexpected state, over shrinking. j = " << j << endl;
			}
			// for optimized network, there should be no expending and shrinking
		}
		cout << "PM" << PM;
		cout << "states_num" << states_num;
		cout << "total_add_computation = " << total_add_computation << endl;		// to be finished


		//cout << "total_compare_computation = " << total_compare_computation << endl;
		//int total_computation = total_add_computation + total_compare_computation;
		//cout << "total_computation = " << total_computation << endl;
	}
	static void viterbi_preprocessing_row_transform_test2() {
		Matrix<GF2> PM(15, 31, {
			1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,1,1,0,1,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,0,1,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,1,1,1,0,1,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0,1,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,1,0,1,0,0,0,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,1,0,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,1,0,0,1,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,1,1,0,1,0,0,1
			});
		cout << "PM" << PM;
		cout << "find_opt_PM::is_final_form(PM) = " << find_opt_PM::is_final_form(PM) << endl;

		PM.row_transformation_to_up_triangle();
		cout << "PM" << PM;
		cout << "find_opt_PM::is_final_form(PM) = " << find_opt_PM::is_final_form(PM) << endl;

		Matrix<int> permute_tmp = find_opt_PM::column_permute(PM);
		cout << "PM" << PM;
		cout << "permute_tmp" << permute_tmp;
		cout << "find_opt_PM::is_final_form(PM) = " << find_opt_PM::is_final_form(PM) << endl;


		Matrix<int> state_num = find_opt_PM::counting_states(PM);
		cout << "state_num" << state_num;

		long long ope_num = find_opt_PM::ope_num_estimation(state_num);
		cout << "ope_num = " << ope_num << endl;
	}
	static void viterbi_optimized_PM_test() {

		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		bch.print_info();

		Matrix<GF2> PM(18, 63, {
			1,0,0,0,0,1,0,1,0,1,1,1,0,1,1,1,0,0,1,0,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,0,1,1,0,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,1,0,0,1,0,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,1,1,0,1,1,1,1,0,1,0,1,1,1,1,1,1,0,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,1,0,1,1,0,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,1,1,0,0,0,1,1,1,1,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,0,1,1,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,1,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,0,0,0,0,1,1,0,0,1,1,0,1,1,0,1,1,0,1,1,0,1,0,1,0,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,0,1,1,0,1,1,1,0,0,1,0,1,1,1,0,1,1,1,0,1,0,1,1,0,1,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,0,1,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,1,0,1,1,1,0,1,0,0,1,0,1,1,1,0,1,1,1,0,1,1,1,1,0,0,1,0,0,0,0,1,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,1,0,0,1,1,0,0,1,1,1,0,1,1,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,0,0,1,0,1,1,0,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,1,1,1,0,0,0,1,1,0,0,0,1,0,1,1,0,0,0,1,0,0,0,1,1
			});
		cout << "PM" << PM;
		//vit.initialize_state();


		Matrix<GF2> u(1, k, 'b');
		cout << "u" << u;
		Matrix<GF2> v = bch.encode(u);
		cout << "v" << v;
		Matrix<int> best_permutation(1, 63, {
			37,  4, 34, 50, 21, 42,
			60, 27, 46, 12, 55,  9,
			 7, 16, 61, 47, 48, 38,
			29, 44, 30, 31, 10, 26,
			 8, 24, 49, 45,  2, 23,
			36, 35, 57, 11,  0, 28,
			39, 51, 15, 43, 53, 22,
			32, 25, 56, 19,  3,  5,
			58, 33,  6,  1, 62, 13,
			18, 59, 40, 20, 14, 17,
			41, 52, 54
			});
		v.permute(best_permutation);
		cout << "v: permute" << v;
		Matrix<my_double> c = BPSK::modulation(v);
		// no error test

		Matrix<my_double> r = c;

		int list_size = 128;
		Viterbi_cyclic_code_r vit(PM);
		Matrix<GF2> dvr = vit.decode_v(r, list_size);

		Viterbi_unordered_map vit3(PM);
		Matrix<GF2> dv3 = vit3.decode_v_once(r, list_size);
		cout << "dvr - dv3" << dvr - dv3;

		Matrix<GF2> dv = dvr.get_row(0);
		cout << "dv" << dv;
		dv.permute_back(best_permutation);
		Matrix<GF2> du = bch.v2u(dv);
		cout << "du" << du;
		cout << "du-u" << du - u;

	}
	static void viterbi_optimized_PM_test2() {
		Matrix<GF2> PM(4, 7, {
			1,1,1,1,0,0,0,
			0,1,1,0,1,0,0,
			0,0,0,1,0,1,0,
			0,0,0,0,1,1,1
			});
		cout << "PM" << PM;
		Viterbi_cyclic_code_r vit(PM);
	}
	static void Matrix_flex_col_test() {
		int r = 5;
		Matrix<int> col_size(1, 5, { 4,2,6,5,8 });
		Matrix_flex_col<double> A(r, col_size);
		A.reset(21);

		cout << "A" << A;
		cout << "A.size() = " << A.size() << endl;
		for (int i = 0; i < r; ++i) {
			int Acol = A.col(i);
			for (int j = 0; j < Acol; ++j) {
				cout << "A(" << i << ", " << j << ") = " << A(i, j) << endl;
			}
		}

		Matrix_flex_col<bool> B;
		cout << "B" << B;
		B.resize(col_size);
		cout << "B" << B;
	}
	static void find_opt_PM_test() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "LL_OSD_hybrid_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		unsigned rand_seed = (unsigned)time(0);
		cout << "rand_seed = " << rand_seed << endl;
		my::set_seed_adv(rand_seed);
		const int m = 6;
		const int t = 5;
		typedef GF2e<m> ty;
		ty::init();
		eBCH<m, t> ebch;
		ebch.print_info();
		Matrix<GF2> PM = ebch.get_parity_matrix();
		Matrix<GF2> GM = ebch.get_generator_matrix();

		find_opt_PM::solve(PM, 1000, GM);

		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;
	}
	static void find_opt_PM_test2() {

		const int row = 14;
		const int col = 50;

		// set cout into file
		char file_name[55] = { 0 };
		sprintf_s(file_name, 55, "eBCH(64,36,12)-sub-PM(%d,%d).txt", row, col);
		FILE* stream1;
		freopen_s(&stream1, file_name, "w", stdout);

		unsigned rand_seed = (unsigned)time(0);
		cout << "rand_seed = " << rand_seed << endl;
		my::set_seed_adv(rand_seed);

		// the sub Matrix of optimum parity Matrix of eBCH(64,36)
		Matrix<GF2> PM(row, col, {
			1,0,0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,0,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,1,1,0,0,1,1,1,0,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,1,0,1,1,0,1
			});

		find_opt_PM::solve(PM, 100000000);

		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;

		fclose(stdout);
	}
	static void find_opt_PM_test3() {

		const int row = 14;
		const int col = 52;

		// set cout into file
		char file_name[55] = { 0 };
		sprintf_s(file_name, 55, "eBCH(64,36,12)-sub-PM(%d,%d).txt", row, col);
		FILE* stream1;
		freopen_s(&stream1, file_name, "w", stdout);

		unsigned rand_seed = (unsigned)time(0);
		cout << "rand_seed = " << rand_seed << endl;
		my::set_seed_adv(rand_seed);

		// the sub Matrix of optimum parity Matrix of eBCH(64,36)
		Matrix<GF2> PM(row, col, {
			1,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,1,0,0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,1,0,1,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,1,1,0,1,1,1,0,0,1
			});

		find_opt_PM::solve(PM, 100000000);

		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;

		fclose(stdout);
	}
	static void find_opt_PM_test4() {

		const int row = 18;
		const int col = 56;
		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "eBCH(64,36,12)-sub-PM(%d,%d).txt", row, col);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		unsigned rand_seed = (unsigned)time(0);
		cout << "rand_seed = " << rand_seed << endl;
		my::set_seed_adv(rand_seed);

		// the sub Matrix of optimum parity Matrix of eBCH(64,36)
		Matrix<GF2> PM(row, col, {
			1,0,1,0,1,0,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,1,0,0,1,1,0,0,1,1,1,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,1,0,0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,1,0,1,1,0,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,1,1,0,1,1,1,0,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,0,1,0,1,1,1,1,1,1,0,0,1,1,0,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,1
			});

		find_opt_PM::solve(PM, 100000000);

		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;

		fclose(stdout);
	}
	static void find_opt_PM_test5() {

		const int row = 14;
		const int col = 50;
		bool if_output_to_file = true;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "eBCH(64,36,12)-sub-PM(%d,%d)-D5.txt", row, col);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		unsigned rand_seed = (unsigned)time(0);
		cout << "rand_seed = " << rand_seed << endl;
		my::set_seed_adv(rand_seed);

		// the sub Matrix of optimum parity Matrix of eBCH(64,36)
		Matrix<GF2> PM(row, col, {
			1,0,1,0,1,0,0,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,0,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,1,1,0,0,1,1,1,0,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,0,1,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,1,0,1,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,1,1,0,1,1,1,0,0,1
			});

		find_opt_PM::solve(PM, 100000000);

		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;

		fclose(stdout);
	}
	static void find_opt_PM_test6() {

		const int row = 12;
		const int col = 48;
		bool if_output_to_file = true;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "eBCH(64,36,12)-sub-PM(%d,%d)-D6.txt", row, col);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		unsigned rand_seed = (unsigned)time(0);
		cout << "rand_seed = " << rand_seed << endl;
		my::set_seed_adv(rand_seed);

		// the sub Matrix of optimum parity Matrix of eBCH(64,36)
		Matrix<GF2> PM(row, col, {
			1,0,1,0,1,0,0,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,
			0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,
			0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,0,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,1,1,0,0,1,1,1,0,1,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,0,1,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1,0,0,0,0,
			0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,0,1
			});

		find_opt_PM::solve(PM, 100000000);

		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;

		fclose(stdout);
	}
	static void find_opt_PM_test7() {

		const int row = 10;
		const int col = 46;
		bool if_output_to_file = true;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "eBCH(64,36,12)-sub-PM(%d,%d)-D7.txt", row, col);		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		unsigned rand_seed = (unsigned)time(0);
		cout << "rand_seed = " << rand_seed << endl;
		my::set_seed_adv(rand_seed);

		// the sub Matrix of optimum parity Matrix of eBCH(64,36)
		Matrix<GF2> PM(row, col, {
			1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,
			0,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,
			0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,
			0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,
			0,0,0,0,1,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,
			0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,0,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,1,1,0,0,1,1,1,0,1,0,0,0,0,
			0,0,0,0,0,0,1,1,0,1,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1,0,0,0,
			0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,
			0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1
			});

		find_opt_PM::solve(PM, 100000000);

		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;

		fclose(stdout);
	}
	static void find_opt_PM_test8() {

		const int row = 12;
		const int col = 52;
		bool if_output_to_file = true;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "eBCH(64,36,12)-sub-PM(%d,%d)-D8.txt", row, col);		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		unsigned rand_seed = (unsigned)time(0);
		cout << "rand_seed = " << rand_seed << endl;
		my::set_seed_adv(rand_seed);

		// the sub Matrix of optimum parity Matrix of eBCH(64,36)
		Matrix<GF2> PM(row, col, {
			1,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,1,0,0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,1,0,1,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,1,1,0,1,1,1,0,0,1
			});

		find_opt_PM::solve(PM, 100000000);

		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;

		fclose(stdout);
	}
	static void find_opt_PM_test9() {

		const int row = 10;
		const int col = 52;
		bool if_output_to_file = true;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "eBCH(64,36,12)-sub-PM(%d,%d)-D9.txt", row, col);		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		unsigned rand_seed = (unsigned)time(0);
		cout << "rand_seed = " << rand_seed << endl;
		my::set_seed_adv(rand_seed);

		// the sub Matrix of optimum parity Matrix of eBCH(64,36)
		Matrix<GF2> PM(row, col, {
			1,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,1,0,0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,1,0,1,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,1,1,0,1,1,1,0,0,1
			});

		find_opt_PM::solve(PM, 100000000);

		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;

		fclose(stdout);
	}
	static void viterbi_optimized_PM_test3() {

		// this is a unsolvable Matrix
		Matrix<GF2> PM(15, 31, {
			1,0,1,1,0,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,1,0,1,0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,1,0,1,0,0,1,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,1,1,1,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,1,1,1,0,0,1,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,0,1,1,1,1,1,1,0,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,1,0,1,1,1,1,1,0,0,1,1,0,1,0,1,1,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,0,1,1,0,0,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,1,1,0,1,0,1,0,1,1,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,1,1,0,0,0,0,0,0,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,0,1,0,0,1,1,0,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,1,0,0,0,0,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1
			});

		cout << "PM" << PM;
		int PM_form = find_opt_PM::is_final_form(PM);
		cout << "PM_form = " << PM_form << endl;

		PM.row_transformation_to_up_triangle();
		cout << "PM" << PM;
		PM_form = find_opt_PM::is_final_form(PM);
		cout << "PM_form = " << PM_form << endl;

		find_opt_PM::column_permute(PM);
		cout << "PM" << PM;
		PM_form = find_opt_PM::is_final_form(PM);
		cout << "PM_form = " << PM_form << endl;

		PM.row_transformation_to_low_triangle();
		cout << "PM" << PM;
		PM_form = find_opt_PM::is_final_form(PM);
		cout << "PM_form = " << PM_form << endl;

		Matrix<int> state_num = find_opt_PM::counting_states(PM);
		long long can_ope_num = find_opt_PM::ope_num_estimation(state_num);

		cout << "state_num" << state_num;
		cout << "can_ope_num = " << can_ope_num << endl;
	}
	static void viterbi_optimized_PM_test4() {
		// redo the case to check that our method is of no problem

		const int m = 5;
		const int t = 3;
		typedef GF2e<m> ty;
		ty::init();
		BCH<m, t> bch;
		Matrix<GF2> PM = bch.get_parity_matrix();
		Matrix<GF2> GM = bch.get_generator_matrix();

		Matrix<int> rand_permutation(1, 31, {
			16,28,22, 2,15,13,
			11, 8,21,30,17,24,
			18,27,26,25, 0,29,
			20,10,12,14,19, 6,
			5, 4, 9,23, 3, 7 ,
			1
			});

		int max_column_iteration = 6;		// test whether it affect the result
		find_opt_PM::further_column_permutation_iteration.resize(1, max_column_iteration + 1, false);
		find_opt_PM::further_column_permutation_iteration.reset(0);
		cout << "rand_permutation" << rand_permutation;
		Matrix<GF2> PM_store = PM;
		PM_store.permute_col(rand_permutation);
		int PM_form;
		int column_iteration;
		Matrix<int> optimze_permutation = find_opt_PM::optimize_PM(PM_store, PM_form, column_iteration, max_column_iteration);
		Matrix<int> state_num = find_opt_PM::counting_states(PM_store);
		long long can_ope_num = find_opt_PM::ope_num_estimation(state_num);

		cout << "PM_form = " << PM_form << endl;
		cout << "column_iteration = " << column_iteration << endl;
		cout << "PM_store" << PM_store;
		cout << "rand_permutation" << rand_permutation;
		cout << "optimze_permutation" << optimze_permutation;

		Matrix<int> total_permutation = rand_permutation;
		total_permutation.permute(optimze_permutation);
		cout << "total_permutation" << total_permutation;
		cout << "state_num" << state_num;
		cout << "can_ope_num = " << can_ope_num << endl;

		// validation
		PM_store.permute_col_back(total_permutation);
		cout << "(GM * PM_store.Transpose()).isZero() = " << (GM * PM_store.Transpose()).isZero() << endl;
		cout << "find_opt_PM::further_column_permutation_iteration" << find_opt_PM::further_column_permutation_iteration;
	}
	static void viterbi_optimized_PM_test5() {

		const int m = 5;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		bch.print_info();

		Matrix<GF2> PM(15, 31, {
			1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,0,1,1,1,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1
			});
		cout << "PM" << PM;
		//vit.initialize_state();


		Matrix<GF2> u(1, k, 'b');
		cout << "u" << u;
		Matrix<GF2> v = bch.encode(u);
		cout << "v" << v;
		Matrix<int> best_permutation(1, 31, {
			25,   30,    4,   26,   14,   22,
			24,   13,   16,   15,    9,   21,
			10,    7,   17,   28,    3,    5,
			0,   20,    8,   29,    1,    6 ,
			12,   11,   27,   23,   18,    2,
			19
			});
		v.permute(best_permutation);
		cout << "v: permute" << v;
		Matrix<my_double> c = BPSK::modulation(v);
		// no error test

		Matrix<my_double> r = c;

		int list_size = 2;
		Viterbi_cyclic_code_r vit(PM);
		Matrix<GF2> dvr = vit.decode_v(r, list_size);

		Viterbi_unordered_map vit3(PM);
		Matrix<GF2> dv3 = vit3.decode_v_once(r, list_size);
		cout << "dvr - dv3" << dvr - dv3;

		Matrix<GF2> dv = dvr.get_row(0);
		cout << "dv" << dv;
		dv.permute_back(best_permutation);
		Matrix<GF2> du = bch.v2u(dv);
		cout << "du" << du;
		cout << "du-u" << du - u;

	}
	static void viterbi_optimized_PM_test6() {
		// stability test

		const int m = 5;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		bch.print_info();
		Matrix<GF2> GM = bch.get_generator_matrix();

		// this PM's border is not totally warped around by 1
		Matrix<GF2> PM(15, 31, {
			1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,0,1,0,0,0,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,1,0,1,0,0,0,0,0,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,0,1,1,0,0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,0,1,0,1,0,1,1,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,1,0,0,1,0,0,0,1,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,1,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,1,1,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0,0,1,0,1,0,1,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,0,0,0,1,1,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,1,0,1,0,0,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,0,1,1,1
			});
		cout << "PM" << PM;
		//vit.initialize_state();


		Matrix<GF2> u(1, k, 'b');
		cout << "u" << u;
		Matrix<GF2> v = bch.encode(u);
		cout << "v" << v;
		Matrix<int> best_permutation(1, 31, {
			10,17,12,28,26,20,
			3,22, 6,23, 0,14,
			5,21,13,15,11, 7 ,
			30,16, 8,27,19, 4 ,
			9, 1,25, 2,18,29 ,
			24
			});

		v.permute(best_permutation);
		cout << "v: permute" << v;
		Matrix<my_double> c = BPSK::modulation(v);
		// no error test

		Matrix<my_double> r = c;

		int list_size = 256;
		Viterbi_cyclic_code_r vit(PM);
		Matrix<GF2> dvr = vit.decode_v(r, list_size);

		Viterbi_unordered_map vit3(PM);
		Matrix<GF2> dv3 = vit3.decode_v_once(r, list_size);
		cout << "dvr - dv3" << dvr - dv3;

		Matrix<GF2> dv = dvr.get_row(0);
		cout << "dv" << dv;
		dv.permute_back(best_permutation);
		Matrix<GF2> du = bch.v2u(dv);
		cout << "du" << du;
		cout << "du-u" << du - u;

		// test if permutation is valid
		PM.permute_col_back(best_permutation);
		cout << "(GM * PM_store.Transpose()).isZero() = " << (GM * PM.Transpose()).isZero() << endl;
	}
	static void viterbi_optimized_PM_test7() {

		const int m = 5;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		bch.print_info();

		Matrix<GF2> PM(15, 31, {
			1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,0,1,1,1,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1
			});
		cout << "PM" << PM;
		//vit.initialize_state();


		Matrix<GF2> u(1, k, 'b');
		cout << "u" << u;
		Matrix<GF2> v = bch.encode(u);
		cout << "v" << v;
		Matrix<int> best_permutation(1, 31, {
			25,   30,    4,   26,   14,   22,
			24,   13,   16,   15,    9,   21,
			10,    7,   17,   28,    3,    5,
			0,   20,    8,   29,    1,    6 ,
			12,   11,   27,   23,   18,    2,
			19
			});
		v.permute(best_permutation);
		cout << "v: permute" << v;
		Matrix<my_double> c = BPSK::modulation(v);
		// no error test

		Matrix<my_double> r = AWGN::pass(c, 0.9);
		cout << "r" << r;

		Matrix<GF2> hdr = BPSK::demodulation(r);
		cout << "hdr" << hdr;

		cout << "error bit" << hdr - v;

		int max_list_size = 32;
		int valid_list_size = 1;

		// now fetch the columns to form PM and unused PM, that is great, i believe

		Matrix<GF2> used_PM = PM.get_rows(Matrix<int>(1, 10, { 0,1,2,3,4,10,11,12,13,14 }));
		Matrix<GF2> unused_PM = PM.get_rows(Matrix<int>(1, 5, { 5,6,7,8,9 }));		// original: 6,7,8
		cout << "used_PM" << used_PM;
		cout << "unused_PM" << unused_PM;
		Matrix<int> states_num = find_opt_PM::counting_states(used_PM);
		cout << "states_num" << states_num;
		long long ope_num_estimate = find_opt_PM::ope_num_estimation(states_num);
		cout << "ope_num_estimate = " << ope_num_estimate << endl;

		Viterbi_optimized vit(used_PM);
		vit.change_unused_PM(unused_PM);
		Matrix<GF2> dvr = vit.decode_v(r, max_list_size, valid_list_size);
		cout << "dvr" << dvr;

		Viterbi_unordered_map vit3(used_PM);
		vit3.change_unused_PM(unused_PM);
		Matrix<GF2> dv3 = vit3.decode_v_once(r, max_list_size);
		cout << "dv3" << dv3;

		Matrix<GF2> dv = dvr.get_row(0);
		cout << "dv" << dv;
		dv.permute_back(best_permutation);
		cout << "dv: permute back" << dv;
		Matrix<GF2> du = bch.v2u(dv);
		cout << "du" << du;
		cout << "du-u" << du - u;


		// testing OSD decoding
		OSD_r osd(bch.get_parity_matrix(), bch.get_d());
		Matrix<my_double> r_permute_back = r;
		r_permute_back.permute_back(best_permutation);
		Matrix<GF2> dv_osd = osd.decode_v(r_permute_back, 1);
		cout << "dv_osd" << dv_osd;
		Matrix<GF2> du_osd = bch.v2u(dv_osd);
		cout << "du-du_osd" << du - du_osd;

		// testing Chase2 decoding
		Matrix<GF2> dv_Chase2 = Chase::decode2_v(bch, r_permute_back, 1);
		cout << "dv_Chase2" << dv_Chase2;
		Matrix<GF2> du_Chase2 = bch.v2u(dv_Chase2);
		cout << "du - du_Chase2" << du - du_Chase2;

		//// testing algebra decoding
		//Matrix<GF2> dv_Algebra = bch.decode_v(hdr);
		//cout << "dv_Algebra" << dv_Algebra;
		//Matrix<GF2> du_Algebra = bch.v2u(dv_Algebra);
		//cout << "du - du_Algebra" << du - du_Algebra;

#ifdef RUN_MSG
		cout << "RUN_MSG" << endl;
#endif // RUN_MSG

	}
	static void viterbi_optimized_PM_simulation() {

		const int simulation_times = 100000;
		const my_double SNR_dB = 2.0;

		// encoding parameter
		const int m = 5;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		int n = bch.get_n();
		bch.print_info();			// BCH(31,16,7)

		Matrix<GF2> PM(15, 31, {
			1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,1,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,1,0,1,1,1,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,1,0,1,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1
			});
		cout << "PM" << PM;

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v_orig = bch.encode(u);
		Matrix<int> best_permutation(1, 31, {
			25,   30,    4,   26,   14,   22,
			24,   13,   16,   15,    9,   21,
			10,    7,   17,   28,    3,    5,
			0,   20,    8,   29,    1,    6 ,
			12,   11,   27,   23,   18,    2,
			19
			});
		cout << "best permutation" << best_permutation;

		Matrix<GF2> v = v_orig;
		v.permute(best_permutation);
		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter

		// fetch the columns to form PM and unused PM, that is great, i believe


		Matrix<GF2> used_PM;
		Matrix<GF2> unused_PM;

		int selected_PM_row = 10;
		switch (selected_PM_row) {
		case (15):
			used_PM = PM.get_rows(Matrix<int>(1, 15, 'N'));
			break;
		case(14):
			used_PM = PM.get_rows(Matrix<int>(1, 14, { 0,1,2,3,4,5,6,8,9,10,11,12,13,14 }));
			unused_PM = PM.get_rows(7);
			break;
		case(12):
			used_PM = PM.get_rows(Matrix<int>(1, 12, { 0,1,2,3,4,5,9,10,11,12,13,14 }));
			unused_PM = PM.get_rows(Matrix<int>(1, 3, { 6,7,8 }));
			break;
		case(10):
			used_PM = PM.get_rows(Matrix<int>(1, 10, { 0,1,2,3,4,10,11,12,13,14 }));
			unused_PM = PM.get_rows(Matrix<int>(1, 5, { 5,6,7,8,9 }));
			break;
		case(8):
			// NOTE: we cannot deal with the case of used_PM having zero column now
			used_PM = PM.get_rows(Matrix<int>(1, 8, { 0,1,2,3,11,12,13,14 }));
			unused_PM = PM.get_rows(Matrix<int>(1, 7, { 4,5,6,7,8,9,10 }));
			break;
		default:
			used_PM = PM.get_rows(Matrix<int>(1, 15, 'N'));
		}

		cout << "used_PM" << used_PM;
		cout << "unused_PM" << unused_PM;

		int max_list_size = 1 << (unused_PM.row() + 1);		// 2*2^{number of unused constrain}
		int valid_list_size = 1;

		Viterbi_optimized vit(used_PM);
		vit.change_unused_PM(unused_PM);

		// performance recording
		int error_frame = 0;
		int num_invalid_list = 0;
		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long float_ope_num_before = my_float_auxiliary_storage::operation_number;
		unsigned long long float_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> r = AWGN::pass_standard(c);

			Matrix<GF2> v_hat = vit.decode_v(r, max_list_size, valid_list_size);

			if (v_hat.size() == 0) {
				num_invalid_list++;
				error_frame++;
			}
			else {
				v_hat.permute_back(best_permutation);
				if (v_orig == v_hat);
				else {
					error_frame++;
				}
			}
		}
		end = clock();

#ifdef count_operation_number
#ifdef use_my_double

		float_ope_num_after = my_float_auxiliary_storage::operation_number;
		cout << "(Viterbi_optimized) float_ope_num = "  \
			<< (float_ope_num_after - float_ope_num_before) / (double)simulation_times << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Viterbi_optimized) GF2_ope_num = "  \
			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times << endl;
#endif // count_operation_number

		cout << "------" << endl;
		cout << "(Viterbi_optimized) ratio_invalid_list = " << num_invalid_list / (double)simulation_times << endl;
		//cout << "(Viterbi_optimized) error_frame = " << error_frame << endl;

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6;
		cout << scientific << setprecision(2);
		cout << "(Viterbi_optimized) time_consume = " << time_consume - pass_AWGN_time << "\ts/iteration;" << endl;
		cout << "(Viterbi_optimized) error_rate = " << error_frame / (my_double)simulation_times << endl;
		cout.unsetf(ios::floatfield);
	}
	static void OSD_r_simulation_old(my_double SNR_dB = 2.0) {

		cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 5000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		int n = bch.get_n();
		int d_min = bch.get_d();
		Matrix<GF2> PM = bch.get_parity_matrix();
		bch.print_info();			// BCH(31,16,7)

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		OSD_r osd(PM, d_min);
		int order = 2;

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
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> recv = AWGN::pass_standard(c);

			Matrix<GF2> v_hat = osd.decode_v(recv, order);
			//Matrix<GF2> v_hat = OSD::decode_v(bch, recv, order);

			if (v == v_hat);
			else {
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();


#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(OSD_r) double_ope_num = "  \
			<< (double_ope_num_after - double_ope_num_before) / (double)simulation_times << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(OSD_r) GF2_ope_num = "  \
			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times << endl;
#endif // count_operation_number

		cout << "------" << endl;

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6;
		cout << "(OSD_r) max_simulation_times = " << simulation_times << endl;
		cout << scientific << setprecision(2);
		cout << "(OSD_r) time_consume = " << time_consume - pass_AWGN_time << "\ts/iteration;" << endl;
		cout << "(OSD_r) error_rate = " << error_frame / (my_double)simulation_times << endl;
		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void OSD_r_simulation_multi_SNR_old() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "OSD_r_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(2, 1, 2, 'd');
		cout << "test_SNR" << test_SNR;
		int len = test_SNR.size();
		for (int i = 0; i < len; ++i) {
			OSD_r_simulation_old(test_SNR(i));
		}
	}
	static void LC_OSD_r_simulation_eBCH_256_131(my_double SNR_dB = 1.0) {

		cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int max_simulation_times = 50000000;		// 1
		const int max_error_frame = 200;

		const int delta = 12;
		//const my_double alpha = 3.2;
		const int Lmax = 1 << (12 * 2);

		// encoding parameter
		const int m = 8;
		const int t = 18;
		GF2e<m>::init();
		eBCH<m, t> ebch;
		int k = ebch.get_k();
		int n = ebch.get_n();
		int d_min = ebch.get_d();
		Matrix<GF2> PM = ebch.get_parity_matrix();
		ebch.print_info();			// eBCH(128,64,22)

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = ebch.encode(u);

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		LC_OSD_r lc_osd(PM, d_min, delta);
		//int len = n - k - delta;
		//my_double Pe_ML = 1.03e-4;		// for SNR = 3 dB
		//Matrix<GF2> p = flip_TEP::get_estimated_flip_pattern(len, Pe_ML);
		//cout << "p" << p;		// have problem, not as good as DAI
		//lc_osd.set_R_pattern(p);
		lc_osd.set_sigma(sigma);

		// performance recording
		int error_frame = 0;
		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
		unsigned long long double_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < max_simulation_times; ++i) {
			Matrix<my_double> r = AWGN::pass_standard(c);

			Matrix<GF2> v_hat = lc_osd.decode_v(r, Lmax);

			//cout << "v-v_hat" << v - v_hat;
			if (v == v_hat);
			else {
				error_frame++;
				if (error_frame == max_error_frame) {			// 200
					max_simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

#ifdef count_operation_number
#ifdef use_my_double

		cout << "------" << endl;
		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) double_ope_num = "  \
			<< (double_ope_num_after - double_ope_num_before) / (double)max_simulation_times << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) GF2_ope_num = "  \
			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)max_simulation_times << endl;
#endif // count_operation_number

		cout << "------" << endl;

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)max_simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << "(LC_OSD_r) max_simulation_times = " << max_simulation_times << endl;
		cout << "(LC_OSD_r) ave_used_list_num = " << lc_osd.total_used_list_num / (double)max_simulation_times << endl;
		cout << "(LC_OSD_r) ML_used_list_num = " << lc_osd.ML_used_list_num / (double)max_simulation_times << endl;
		cout << "(LC_OSD_r) reach_max_list_num_rate= " << lc_osd.times_reach_max_list_num / (double)max_simulation_times << endl;
		cout << scientific << setprecision(2);
		cout << "(LC_OSD_r) time_consume = " << time_consume - pass_AWGN_time << "\ts/iteration" << endl;
		cout << "(LC_OSD_r) error = " << error_frame / (my_double)max_simulation_times << endl;
		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void LC_OSD_r_simulation_eBCH_128_64(my_double SNR_dB = 1.0) {

		cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int max_simulation_times = 50000000;		// 1
		const int max_error_frame = 200;

		const int delta = 8;
		//const my_double alpha = 3.2;
		const int Lmax = 16384;

		// encoding parameter
		const int m = 7;
		const int t = 10;
		GF2e<m>::init();
		eBCH<m, t> ebch;
		int k = ebch.get_k();
		int n = ebch.get_n();
		int d_min = ebch.get_d();
		Matrix<GF2> PM = ebch.get_parity_matrix();
		ebch.print_info();			// eBCH(128,64,22)

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = ebch.encode(u);

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		LC_OSD_r lc_osd(PM, d_min, delta);
		//int len = n - k - delta;
		//my_double Pe_ML = 1.03e-4;		// for SNR = 3 dB
		//Matrix<GF2> p = flip_TEP::get_estimated_flip_pattern(len, Pe_ML);
		//cout << "p" << p;		// have problem, not as good as DAI
		//lc_osd.set_R_pattern(p);
		lc_osd.set_sigma(sigma);

		// performance recording
		int error_frame = 0;
		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
		unsigned long long double_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < max_simulation_times; ++i) {
			Matrix<my_double> r = AWGN::pass_standard(c);

			Matrix<GF2> v_hat = lc_osd.decode_v(r, Lmax);

			//cout << "v-v_hat" << v - v_hat;
			if (v == v_hat);
			else {
				error_frame++;
				if (error_frame == max_error_frame) {			// 200
					max_simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

#ifdef count_operation_number
#ifdef use_my_double

		cout << "------" << endl;
		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) double_ope_num = "  \
			<< (double_ope_num_after - double_ope_num_before) / (double)max_simulation_times << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) GF2_ope_num = "  \
			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)max_simulation_times << endl;
#endif // count_operation_number

		cout << "------" << endl;

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)max_simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << "(LC_OSD_r) max_simulation_times = " << max_simulation_times << endl;
		cout << "(LC_OSD_r) ave_used_list_num = " << lc_osd.total_used_list_num / (double)max_simulation_times << endl;
		cout << "(LC_OSD_r) ML_used_list_num = " << lc_osd.ML_used_list_num / (double)max_simulation_times << endl;
		cout << "(LC_OSD_r) reach_max_list_num_rate= " << lc_osd.times_reach_max_list_num / (double)max_simulation_times << endl;
		cout << scientific << setprecision(2);
		cout << "(LC_OSD_r) time_consume = " << time_consume - pass_AWGN_time << "\ts/iteration" << endl;
		cout << "(LC_OSD_r) error_rate = " << error_frame / (my_double)max_simulation_times << endl;
		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void LC_OSD_r_simulation_BCH_63_45(my_double SNR_dB = 1.0) {

		cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int max_simulation_times = 5000000;		// 1

		// decoding parameter
		const int delta = 6;
		//const my_double beta = 4;
		const int Lmax = 1024;

		// encoding parameter
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		int n = bch.get_n();
		int d_min = bch.get_d();
		Matrix<GF2> PM = bch.get_parity_matrix();
		bch.print_info();			// eBCH(128,64,22)

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		LC_OSD_r lc_osd(PM, d_min, delta);
		lc_osd.set_sigma(sigma);		// this is important

		// performance recording
		int error_frame = 0;
		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
		unsigned long long double_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;
#endif // count_operation_number

		//srand(10);
		start = clock();
		for (int i = 0; i < max_simulation_times; ++i) {
			Matrix<my_double> r = AWGN::pass_standard(c);

			Matrix<GF2> v_hat = lc_osd.decode_v(r, Lmax);

			//cout << "v-v_hat" << v - v_hat;
			if (v == v_hat);
			else {
				error_frame++;
				if (error_frame == 200) {
					max_simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

#ifdef count_operation_number
#ifdef use_my_double

		cout << "------" << endl;
		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) double_ope_num = "  \
			<< (double_ope_num_after - double_ope_num_before) / (double)max_simulation_times << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) GF2_ope_num = "  \
			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)max_simulation_times << endl;
#endif // count_operation_number

		cout << "------" << endl;

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)max_simulation_times;
		double pass_AWGN_time = 6e-6;
		cout << "(LC_OSD_r) max_simulation_times = " << max_simulation_times << endl;
		cout << "(LC_OSD_r) ave_used_list_num = " << lc_osd.total_used_list_num / (double)max_simulation_times << endl;
		cout << "(LC_OSD_r) ave_ML_used_list_num = " << lc_osd.ML_used_list_num / (double)max_simulation_times << endl;
		cout << "(LC_OSD_r) reach_max_list_num_rate= " << lc_osd.times_reach_max_list_num / (double)max_simulation_times << endl;
		cout << scientific << setprecision(2);
		cout << "(LC_OSD_r) time_consume = " << time_consume - pass_AWGN_time << "\ts/iteration" << endl;
		cout << "(LC_OSD_r) error_rate = " << error_frame / (my_double)max_simulation_times << endl;
		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void LC_OSD_r_simulation_BCH_63_30(my_double SNR_dB = 1.0) {

		cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int max_simulation_times = 50000000;		// 1
		const int max_error_frame = 200;

		const int delta = 3;
		//const my_double alpha = 3.2;
		const int Lmax = 128;

		// encoding parameter
		const int m = 6;
		const int t = 6;
		GF2e<m>::init();
		BCH<m, t> ebch;
		int k = ebch.get_k();
		int n = ebch.get_n();
		int d_min = ebch.get_d();
		Matrix<GF2> PM = ebch.get_parity_matrix();
		ebch.print_info();			// eBCH(128,64,22)

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = ebch.encode(u);

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		LC_OSD_r lc_osd(PM, d_min, delta);
		//int len = n - k - delta;
		//my_double Pe_ML = 1.03e-4;		// for SNR = 3 dB
		//Matrix<GF2> p = flip_TEP::get_estimated_flip_pattern(len, Pe_ML);
		//cout << "p" << p;		// have problem, not as good as DAI
		//lc_osd.set_R_pattern(p);
		lc_osd.set_sigma(sigma);

		// performance recording
		int error_frame = 0;
		clock_t start, end;

#ifdef count_operation_number
#ifdef use_my_double

		unsigned long long double_ope_num_before = my_double_auxiliary_storage::operation_number;
		unsigned long long double_ope_num_after;
#endif // use_my_double

		unsigned long long GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		unsigned long long GF2_ope_num_after;
#endif // count_operation_number

		start = clock();
		for (int i = 0; i < max_simulation_times; ++i) {
			Matrix<my_double> r = AWGN::pass_standard(c);

			Matrix<GF2> v_hat = lc_osd.decode_v(r, Lmax);

			//cout << "v-v_hat" << v - v_hat;
			if (v == v_hat);
			else {
				error_frame++;
				if (error_frame == max_error_frame) {			// 200
					max_simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

#ifdef count_operation_number
#ifdef use_my_double

		cout << "------" << endl;
		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) double_ope_num = "  \
			<< (double_ope_num_after - double_ope_num_before) / (double)max_simulation_times << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(LC_OSD_r) GF2_ope_num = "  \
			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)max_simulation_times << endl;
#endif // count_operation_number

		cout << "------" << endl;

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)max_simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << "(LC_OSD_r) max_simulation_times = " << max_simulation_times << endl;
		cout << "(LC_OSD_r) ave_used_list_num = " << lc_osd.total_used_list_num / (double)max_simulation_times << endl;
		cout << "(LC_OSD_r) ML_used_list_num = " << lc_osd.ML_used_list_num / (double)max_simulation_times << endl;
		cout << "(LC_OSD_r) reach_max_list_num_rate= " << lc_osd.times_reach_max_list_num / (double)max_simulation_times << endl;
		cout << scientific << setprecision(2);
		cout << "(LC_OSD_r) time_consume = " << time_consume - pass_AWGN_time << "\ts/iteration" << endl;
		cout << "(LC_OSD_r) error_rate = " << error_frame / (my_double)max_simulation_times << endl;
		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void LC_OSD_r_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "LC_OSD_r_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(1, 0.5, 4, 'd');
		//Matrix<my_double> test_SNR(0, 0.5, 6, 'd');
		cout << "test_SNR" << test_SNR;
		int len = test_SNR.size();
		for (int i = 0; i < len; ++i) {
			//LC_OSD_r_simulation_eBCH_256_131(test_SNR(i));			// fail to decode such a long code
			LC_OSD_r_simulation_eBCH_128_64(test_SNR(i));
			//LC_OSD_r_simulation_BCH_63_45(test_SNR(i));
			//LC_OSD_r_simulation_BCH_63_30(test_SNR(i));
		}

		if (if_output_to_file == false) {
			system("pause");
		}

		/*Matrix<my_double> test_SNR(2, 1, 5, 'd');
		cout << "test_SNR" << test_SNR;
		int len = test_SNR.size();
		for (int i = 0; i < len; ++i) {
			LC_OSD_r_simulation_BCH_63_45(test_SNR(i));
		}*/
	}
	static void Chase2_simulatio_old() {

		const int simulation_times = 100000;
		const my_double SNR_dB = 2.0;

		// encoding parameter
		const int m = 5;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		int n = bch.get_n();
		int d_min = bch.get_d();
		Matrix<GF2> PM = bch.get_parity_matrix();
		bch.print_info();			// BCH(31,16,7)

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		int eta = 4;

		// performance recording
		int error_frame = 0;
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
			Matrix<my_double> r = AWGN::pass_standard(c);

			Matrix<GF2> v_hat = Chase::decode2_v(bch, r, eta);

			if (v == v_hat);
			else {
				error_frame++;
			}
		}
		end = clock();

#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(Chase2) double_ope_num = " \
			<< (double_ope_num_after - double_ope_num_before) / (double)simulation_times << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Chase2) GF2_ope_num = " \
			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times << endl;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << "(Chase2) GF2e_ope_num = " \
			<< (GF2e_ope_num_after - GF2e_ope_num_before) / (double)simulation_times << endl;
#endif // count_operation_number

		cout << "------" << endl;

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6;
		cout << scientific << setprecision(2);
		cout << "(Chase2) time_consume = " << time_consume - pass_AWGN_time << "\ts/iteration;" << endl;
		cout << "(Chase2) error_rate = " << error_frame / (my_double)simulation_times << endl;
		cout.unsetf(ios::floatfield);
	}
	static void my_float_test() {
		my_float mf = 12;
		my_double md = 15;
		double x = 3;
		double y = 1;

		cout << "mf=" << mf << endl;
		cout << "md=" << md << endl;

		mf = float(md);		// you have to convert my_double to float explicitly, then use my_float
		cout << "mf: new = " << mf << endl;

		mf = 5;
		md = mf;
		cout << "md: new = " << md << endl;

		mf = float(mf + md);	// looking good and nromal
		cout << "mf: new2 = " << mf << endl;

		md = mf + md;
		cout << "md: new2 = " << md << endl;

		md = -md;
		cout << "md: negtive = " << md << endl;
	}
	static void rand_test() {
		//srand((unsigned)time(0));		// you can set the random seed, to generate different random number for each run
		my::set_seed_adv((unsigned)time(0));
		int n = 64;
		Matrix<int> natual(1, n, 'N');
		Matrix<int> rand_permutation = natual.get_random_element(n);
		cout << "rand_permutation" << rand_permutation;
	}
	static void sort_test() {
		Matrix<int> A(1, 15, my::rand_int);
		cout << "A" << A << endl;

		A.sort('>');
		cout << "A" << A;

		A.sort();
		cout << "A" << A;


		Matrix<int> sort_ind = A.sort_with_ind('>');
		cout << "A" << A;
		cout << "sort_ind" << sort_ind;

		sort_ind = A.sort_with_ind();
		cout << "A" << A;
		cout << "sort_ind" << sort_ind;

		Matrix<int> B(1, 8, my::rand_int);
		cout << "B" << B;
		B.sort();
		cout << "B" << B;

		Matrix<int> C;
		C.merge_sort_lt(A, B, 12);
		cout << "C" << C;

		Matrix<int> D(1, 7, { 1,2,4,6,6,6,8 });
		cout << "D" << D;
		for (int i = 0; i < 10; ++i) {
			int lb = D.binary_search(i);
			cout << "lb(" << i << ") = " << lb << endl;
		}
	}
	static void generate_systematic_generator_any_pos_test() {
		const int m = 5;
		const int k = 27;
		typedef GF2e<m> ty;
		ty::init();
		cout << "GF2e_auxiliary_storage::operation_number: = " << GF2e_auxiliary_storage::operation_number << endl;

		RS<m, k> rs;
		rs.print_info();
		int n = rs.get_n();
		int d = rs.get_d();

		Matrix<int> info_choose(1, (1 << m) - 1, 'N');
		Matrix<int> info_set = info_choose.get_random_element(k);
		//Matrix<int> info_set(1, 7, { 2,1,4,14,10,9,13 });
		//Matrix<int> info_set(1, 5, { 2,1,4,14,10 });
		//info_set.sort();
		cout << "info_set" << info_set;

		unsigned long long before, after, before_mul, after_mul, before_add, after_add;
		cout << "GF2e_auxiliary_storage::operation_number: before polynomial generation = " \
			<< GF2e_auxiliary_storage::operation_number << endl;

		before = GF2e_auxiliary_storage::operation_number;
		before_mul = GF2e_auxiliary_storage::mul_number;
		before_add = GF2e_auxiliary_storage::add_number;
		rs.generate_systematic_generator_any_pos_naive(info_set);
		cout << "rs" << rs.generator_M_systematic_any_pos << endl;

		after = GF2e_auxiliary_storage::operation_number;
		after_mul = GF2e_auxiliary_storage::mul_number;
		after_add = GF2e_auxiliary_storage::add_number;
		cout << "computation cost (add & mul) = " << after - before << endl;
		cout << "computation cost (mul) = " << after_mul - before_mul << endl;
		cout << "computation cost (add) = " << after_add - before_add << endl;


		before = GF2e_auxiliary_storage::operation_number;
		before_mul = GF2e_auxiliary_storage::mul_number;
		before_add = GF2e_auxiliary_storage::add_number;
		rs.generate_systematic_generator_any_pos(info_set);
		Matrix<ty> store = rs.generator_M_systematic_any_pos;
		cout << "rs" << rs.generator_M_systematic_any_pos << endl;

		after = GF2e_auxiliary_storage::operation_number;
		after_mul = GF2e_auxiliary_storage::mul_number;
		after_add = GF2e_auxiliary_storage::add_number;
		cout << "computation cost (add & mul) = " << after - before << endl;
		cout << "computation cost (mul) = " << after_mul - before_mul << endl;
		cout << "computation cost (add) = " << after_add - before_add << endl;


		before = GF2e_auxiliary_storage::operation_number;
		before_mul = GF2e_auxiliary_storage::mul_number;
		before_add = GF2e_auxiliary_storage::add_number;
		rs.generate_systematic_generator_any_pos_best(info_set);
		cout << "rs" << rs.generator_M_systematic_any_pos << endl;

		after = GF2e_auxiliary_storage::operation_number;
		after_mul = GF2e_auxiliary_storage::mul_number;
		after_add = GF2e_auxiliary_storage::add_number;
		cout << "computation cost (add & mul) = " << after - before << endl;
		cout << "computation cost (mul) = " << after_mul - before_mul << endl;
		cout << "computation cost (add) = " << after_add - before_add << endl;

		// checking
		cout << "(store - rs.generator_M_systematic_any_pos).isZero() = " << (store - rs.generator_M_systematic_any_pos).isZero() << endl;
	}
	static void viterbi_demonstrate() {
		GF2e<3>::init();
		BCH<3, 1> bch;
		Matrix<GF2> PM = bch.get_parity_matrix();
		cout << "PM" << PM;

		Viterbi_unordered_map vit(PM);
		vit.initialize_state();
		vit.print_state();
	}
	static void non_0_pos_test() {
		const int m = 6;
		const int t = 2;

		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		int n = bch.get_n();

		Matrix<GF2> u(1, k, '0');
		u(0) = 1;
		Matrix<GF2> v = bch.encode(u);
		Matrix<GF2> PM = bch.get_parity_matrix();
		cout << "PM" << PM;

		cout << "v" << v;

		cout << "GF2_auxiliary_storage::operation_number = " << GF2_auxiliary_storage::operation_number << endl;
		cout << "PM.check_inner_product(v) = " << PM.check_inner_product(v) << endl;
		cout << "GF2_auxiliary_storage::operation_number = " << GF2_auxiliary_storage::operation_number << endl;

		Matrix<int> non_0_ind_aux(1, v.size());
		cout << "GF2_auxiliary_storage::operation_number = " << GF2_auxiliary_storage::operation_number << endl;
		cout << "PM.check_inner_product_4_GF2(v, Matrix<GF2>(), non_0_ind_aux) = " \
			<< PM.check_inner_product_4_GF2(v, Matrix<GF2>(), non_0_ind_aux) << endl;
		cout << "GF2_auxiliary_storage::operation_number = " << GF2_auxiliary_storage::operation_number << endl;
	}
	static void LL_OSD_simulation(my_double SNR_dB = 2.0) {

		cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 5000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> PM = bch.get_parity_matrix();
		bch.print_info();			// BCH(31,16,7)

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
		LL_OSD<m, t, k_prime> llosd;
		int order = 3;

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;
		//Matrix<int> error_frame_ind(1, 200, 'v');

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

			/*int early_stop_before = llosd.num_early_stop;
			unsigned long long total_used_list_num_before = llosd.total_used_list_num;*/

			Matrix<GF2> v_hat = llosd.decode_v(recv, order);

			/*unsigned long long total_used_list_num_after = llosd.total_used_list_num;
			int early_stop_after = llosd.num_early_stop;*/

			//Matrix<GF2> v_hat = OSD::decode_v(bch, recv, order);

			if (v == v_hat);
			else {
				// record the error frame

				//error_frame_ind.push_back(i);
				/*cout << "--------- error frame index i = " << i << "--------" << endl;
				if (early_stop_before != early_stop_after) {
					cout << "early stop using" << endl;
					cout << "recv" << recv;
					cout << "used_list = " << total_used_list_num_after - total_used_list_num_before << endl;
				}*/

				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

		/*cout << "------" << endl;
		cout << "error_frame_ind" << error_frame_ind;
		cout << "------" << endl;*/

#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(LLOSD) double_ope_num = "  \
			<< (double_ope_num_after - double_ope_num_before) / (double)simulation_times << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(LLOSD) GF2_ope_num = "  \
			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times << endl;

		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		cout << "(LLOSD) GF2e_ope_num = "  \
			<< (GF2e_ope_num_after - GF2e_ope_num_before) / (double)simulation_times << endl;
#endif // count_operation_number

		cout << "------" << endl;

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << "(LLOSD) max_simulation_times = " << simulation_times << endl;
		cout << scientific << setprecision(2);
		cout << "(LLOSD) time_consume = " << time_consume - pass_AWGN_time << "\ts/iteration;" << endl;
		cout << "(LLOSD) error_rate = " << error_frame / (my_double)simulation_times << endl;

		cout << "------" << endl;

		cout << "(LLOSD) invalid_list_rate = " << llosd.num_invalid_list / (my_double)simulation_times << endl;
		cout << "(LLOSD) early_stop_rate = " << llosd.num_early_stop / (my_double)simulation_times << endl;
		cout << "(LLOSD) ave_used_list_num = " << llosd.total_used_list_num / (my_double)simulation_times << endl;

		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void LL_OSD_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "LL_OSD_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(2, 0.5, 5, 'd');
		cout << "test_SNR" << test_SNR;
		int len = test_SNR.size();
		for (int i = 0; i < len; ++i) {
			LL_OSD_simulation(test_SNR(i));
		}
	}
	static void LL_OSD_hybrid_simulation(bool& is_first_time, my_double SNR_dB = 2.0) {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 5000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> PM = bch.get_parity_matrix();
		//bch.print_info();			// BCH(31,16,7)

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
		const int select_row_num = 6;
		//const my_double beta = 3;
		const int Lmax = 128;
		const int OSD_order = 1;
		LL_OSD_Hybrid_flip_Viterbi<m, t, k_prime, select_row_num> llosd_hybrid;
		llosd_hybrid.set_sigma(sigma);

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;

		if (is_first_time) {

			bch.print_info();
			cout << "max_error_frame          = " << max_error_frame << endl;
			cout << "simulation type          = LL_OSD_Hybrid_flip_Viterbi" << endl;
			cout << "select_row_num           = " << select_row_num << endl;
			cout << "Lmax                     = " << Lmax << endl;

			cout << "\n---------------------------------------------------------------------------------------------------------------------"
				<< "---------------------------------------------------------------------------------------------------------------" << endl;


			/*	cout << setw(7) << "SNR_dB" << setw(20) << "simulation_times" << setw(20) << \
					"error_rate" << setw(20) << "" << setw(20) << "invalid_list_rate" << setw(20) << "early_stop_rate" \
					<< setw(20) << "ML_used_list_num" << setw(20) << "max_list_rate";*/


			cout << setw(7) << "SNR_dB" << setw(20) << "simulation_times" << setw(20) << \
				"error_rate" << setw(20) << "invalid_list_rate" << setw(20) << "ave_used_list" << setw(20) << "ML_used_list_num" \
				<< setw(20) << "early_stop_rate" << setw(20) << "max_list_rate";

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

			Matrix<GF2> v_hat = llosd_hybrid.decode_v(recv, Lmax, OSD_order);
			//Matrix<GF2> v_hat = OSD::decode_v(bch, recv, order);

			if (v == v_hat);
			else {
				// record the error frame
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

		//#ifdef count_operation_number
		//#ifdef use_my_double
		//		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		//		cout << "(LL_OSD_Hybrid_flip_Viterbi) double_ope_num = "  \
		//			<< (double_ope_num_after - double_ope_num_before) / (double)simulation_times << endl;
		//#endif // use_my_double
		//		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		//		cout << "(LL_OSD_Hybrid_flip_Viterbi) GF2_ope_num = "  \
		//			<< (GF2_ope_num_after - GF2_ope_num_before) / (double)simulation_times << endl;
		//		GF2e_ope_num_after = GF2e_auxiliary_storage::operation_number;
		//		cout << "(LL_OSD_Hybrid_flip_Viterbi) GF2e_ope_num = "  \
		//			<< (GF2e_ope_num_after - GF2e_ope_num_before) / (double)simulation_times << endl;
		//#endif // count_operation_number
		//		cout << "------" << endl;
		//		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		//		double pass_AWGN_time = 2.85e-6;
		//		cout << "(LL_OSD_Hybrid_flip_Viterbi) max_simulation_times = " << simulation_times << endl;
		//		cout << scientific << setprecision(2);
		//		cout << "(LL_OSD_Hybrid_flip_Viterbi) time_consume = " << time_consume - pass_AWGN_time << "\ts/iteration;" << endl;
		//		cout << "(LL_OSD_Hybrid_flip_Viterbi) error_rate = " << error_frame / (my_double)simulation_times << endl;
		//		cout << "------" << endl;
		//		cout << "(LLOSD) invalid_list_rate = " << llosd_hybrid.num_invalid_list / (my_double)simulation_times << endl;
		//		cout << "(LLOSD) early_stop_rate = " << llosd_hybrid.num_early_stop / (my_double)simulation_times << endl;
		//		cout << "(LLOSD) ave_used_list = " << llosd_hybrid.total_used_list_num / (my_double)simulation_times << endl;
		//		cout.unsetf(ios::floatfield);
		//		cout << setprecision(6);

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << setw(7) << fixed << setprecision(1) << SNR_dB;
		cout << setw(20) << simulation_times;
		cout << scientific << setprecision(2);
		cout << setw(20) << error_frame / (double)simulation_times;

		cout << setw(20) << llosd_hybrid.num_invalid_list / (double)simulation_times;
		cout << setw(20) << llosd_hybrid.total_used_list_num / (double)simulation_times;
		cout << setw(20) << llosd_hybrid.ML_used_list_num / (double)simulation_times;
		cout << setw(20) << llosd_hybrid.num_early_stop / (double)simulation_times;
		cout << setw(20) << llosd_hybrid.times_reach_max_list_num / (double)simulation_times;

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
	static void LL_OSD_hybrid_simulation_multi_SNR() {

		bool if_output_to_file = true;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "LL_OSD_hybrid_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(2, 0.5, 5, 'd');
		//Matrix<my_double> test_SNR(2, 0.5, 5, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			LL_OSD_hybrid_simulation(is_first_time, test_SNR(i));
		}
	}

	static void OSD_r_simulation(bool& is_first_time, my_double SNR_dB = 2.0) {

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
		OSD_r osd(PM, d_min);

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
	static void OSD_r_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "OSD_r_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(1, 1, 6, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			OSD_r_simulation(is_first_time, test_SNR(i));
		}
	}
	static void Chase2_simulation(bool& is_first_time, my_double SNR_dB = 2.0) {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 5000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		//bch.print_info();			// BCH(31,16,7)

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
		int eta = 6;
		Chase::num_early_stop = 0;
		Chase::total_used_list_num = 0;

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;

		if (is_first_time) {

			bch.print_info();
			cout << "max_error_frame          = " << max_error_frame << endl;
			cout << "simulation type          = Chase2-" << eta << endl;

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
			Matrix<GF2> v_hat = Chase::decode2_v(bch, recv, eta);

			if (v == v_hat);
			else {
				// record the error frame
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << setw(7) << fixed << setprecision(1) << SNR_dB;
		cout << setw(20) << simulation_times;
		cout << scientific << setprecision(2);
		cout << setw(20) << error_frame / (double)simulation_times;

		cout << setw(20) << Chase::total_used_list_num / (double)simulation_times;
		cout << setw(20) << Chase::num_early_stop / (double)simulation_times;

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
	static void Chase2_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "Chase2_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(2, 0.5, 5, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			Chase2_simulation(is_first_time, test_SNR(i));
		}
	}
	static void GMD_simulation(bool& is_first_time, my_double SNR_dB = 2.0) {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 5000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		//bch.print_info();			// BCH(31,16,7)

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
		GMD::num_early_stop = 0;
		GMD::total_used_list_num = 0;

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;

		if (is_first_time) {

			bch.print_info();
			cout << "max_error_frame          = " << max_error_frame << endl;
			cout << "simulation type          = GMD" << endl;

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
			Matrix<GF2> v_hat = GMD::decode_v(bch, recv);

			if (v == v_hat);
			else {
				// record the error frame
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << setw(7) << fixed << setprecision(1) << SNR_dB;
		cout << setw(20) << simulation_times;
		cout << scientific << setprecision(2);
		cout << setw(20) << error_frame / (double)simulation_times;

		cout << setw(20) << GMD::total_used_list_num / (double)simulation_times;
		cout << setw(20) << GMD::num_early_stop / (double)simulation_times;

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
	static void GMD_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "GMD_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(2, 0.5, 5, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			GMD_simulation(is_first_time, test_SNR(i));
		}
	}
	static void Hybrid_Chase2_OSD_simulation(bool& is_first_time, my_double SNR_dB = 2.0) {

		//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
		int simulation_times = 5000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 6;
		const int t = 3;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> PM = bch.get_parity_matrix();
		//bch.print_info();			// BCH(31,16,7)

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
		int eta = 1;
		int order = 1;
		Hybrid_Chase2_OSD::num_Chase2_stop = 0;
		Hybrid_Chase2_OSD::num_early_stop = 0;
		Hybrid_Chase2_OSD::total_used_list_num = 0;

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;

		if (is_first_time) {

			bch.print_info();
			cout << "max_error_frame          = " << max_error_frame << endl;
			cout << "simulation type          = Hybrid_Chase2_OSD: eta = " << eta << ", tau = " << order << endl;

			cout << "\n------------------------------------------------------------------------------------------------------------"
				<< "--------------------------------------------------------------------------------" << endl;


			cout << setw(7) << "SNR_dB" << setw(20) << "simulation_times" << setw(20) << "error_rate" \
				<< setw(20) << "ave_used_list" << setw(20) << "early_stop_rate" << setw(20) << "Chase2_stop_rate";

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
			Matrix<GF2> v_hat = Hybrid_Chase2_OSD::decode_v(bch, recv, eta, order);

			if (v == v_hat);
			else {
				// record the error frame
				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << setw(7) << fixed << setprecision(1) << SNR_dB;
		cout << setw(20) << simulation_times;
		cout << scientific << setprecision(2);
		cout << setw(20) << error_frame / (double)simulation_times;

		cout << setw(20) << Hybrid_Chase2_OSD::total_used_list_num / (double)simulation_times;
		cout << setw(20) << Hybrid_Chase2_OSD::num_early_stop / (double)simulation_times;
		cout << setw(20) << Hybrid_Chase2_OSD::num_Chase2_stop / (double)simulation_times;


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
	static void Hybrid_Chase2_OSD_simulation_multi_SNR() {

		bool if_output_to_file = true;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "Hybrid_Chase2_OSD_simulation_multi_SNR.txt");		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(2, 0.5, 5, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			Hybrid_Chase2_OSD_simulation(is_first_time, test_SNR(i));
		}
	}
	static void Hamming_code_example() {
		GF2e<3>::init();
		BCH<3, 1> bch;
		bch.print_info();


		int k = bch.get_k();
		const int n = bch.get_n();
		const int d_min = bch.get_d();

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		cout << "v" << v;

		Matrix<my_double> c = BPSK::modulation(v);
		cout << "c" << c;
		//Matrix<my_double> channel_error = Matrix<my_double>(1, n, { 1.4,0.1,0.5,-0.2,-0.3,0,-1.04 });
		//Matrix<my_double> channel_error = Matrix<my_double>(1, n, { 0.6,0.1,0.5,1.6,-0.3,0,-1.04 });
		Matrix<my_double> channel_error = Matrix<my_double>(1, n, { 0.6,0.1,0.5,2,-0.3,-0.1,-1.04 });
		Matrix<my_double> recv = c + channel_error;
		cout << "recv" << recv;
		Matrix<GF2> z = BPSK::demodulation(recv);
		cout << "z" << z;
		cout << "e" << z - v;

		Matrix<GF2> dv_algebra = bch.decode_v(z);
		cout << "dv_algebra" << dv_algebra;
		cout << "dv_algebra == v : " << (dv_algebra == v) << endl;

		Matrix<GF2> dv_GMD = GMD::decode_v(bch, recv);
		cout << "dv_GMD" << dv_GMD;
		cout << "dv_GMD == v : " << (dv_GMD == v) << endl;

		Matrix<GF2> dv_Chase2 = Chase::decode2_v(bch, recv, 2);
		cout << "dv_Chase2" << dv_Chase2;
		cout << "dv_Chase2 == v : " << (dv_Chase2 == v) << endl;

		Matrix<GF2> tz = Matrix<GF2>(1, n, { 0,0,1,0,1,0,0 });
		dv_algebra = bch.decode_v(tz);
		cout << "dv_algebra_tz" << dv_algebra;
		cout << "dv_algebra_tz == v : " << (dv_algebra == v) << endl;

		Matrix<GF2> tz2 = Matrix<GF2>(1, n, { 1,0,1,0,1,0,1 });
		dv_algebra = bch.decode_v(tz2);
		cout << "dv_algebra_tz2" << dv_algebra;
		cout << "dv_algebra_tz2 == v : " << (dv_algebra == v) << endl;

		Matrix<GF2> tz3 = Matrix<GF2>(1, n, { 0,0,1,0,1,0,1 });
		dv_algebra = bch.decode_v(tz3);
		cout << "dv_algebra_tz3" << dv_algebra;
		cout << "dv_algebra_tz3 == v : " << (dv_algebra == v) << endl;

		OSD_r osd(bch.get_parity_matrix(), bch.get_d());
		Matrix<GF2> dv_OSD_1 = osd.decode_v(recv, 1);
		cout << "dv_OSD_1" << dv_OSD_1;
		cout << "dv_OSD_1 == v : " << (dv_OSD_1 == v) << endl;
	}
	static void list_viterbi_zero_columns() {
		const int n = 51;
		const int r = 6;		// good

		Viterbi_unordered_map vit(n, r);

		Matrix<GF2> PM(r, n, {
			0,1,1,1,0,0,1,1,1,1,1,0,1,1,0,1,0,0,1,1,1,1,0,0,0,1,1,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,
			1,1,1,1,1,0,0,1,1,0,0,1,0,0,0,1,0,1,1,0,0,0,1,0,1,1,0,0,1,0,1,1,1,0,0,0,1,0,1,0,1,1,1,0,0,0,1,0,0,0,0,
			0,0,1,1,1,1,0,1,1,1,0,0,1,0,1,1,1,0,1,0,1,0,0,1,0,0,0,1,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,
			1,0,1,0,0,1,0,0,1,0,0,0,0,1,1,0,1,1,1,1,0,0,0,0,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,
			0,1,0,1,1,0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,1,0,
			1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1
			});

		// this PM without zero columns
		/*Matrix<GF2> PM(r, n, {
			0,1,1,1,0,0,1,1,1,1,1,0,1,1,0,1,0,0,1,1,1,1,0,0,0,1,1,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,
			1,1,1,1,1,0,0,1,1,0,0,1,0,0,0,1,0,1,1,0,0,0,1,0,1,1,0,0,1,0,1,1,1,0,0,0,1,0,1,0,1,1,1,0,0,0,1,0,0,0,0,
			0,0,1,1,1,1,0,1,1,1,0,0,1,0,1,1,1,0,1,0,1,0,0,1,0,0,0,1,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,0,1,0,0,0,
			1,0,1,0,0,1,0,0,1,0,0,0,0,1,1,0,1,1,1,1,0,0,0,0,0,1,1,0,0,1,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,
			0,1,0,1,1,0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,0,0,0,1,0,1,1,0,0,1,0,1,0,0,0,0,1,0,
			1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1
			});*/

		cout << "PM" << PM;
		Matrix<GF2> PM_preocessed = PM;
		Matrix<int> col_permute_record = find_opt_PM::column_permute(PM_preocessed);

		cout << "PM_preocessed" << PM_preocessed;
		cout << "col_permute_record = " << col_permute_record;

		// try to permute back
		PM_preocessed = PM_preocessed.combine_right(Matrix<GF2>(r, n - PM_preocessed.col(), '0'));
		PM_preocessed.permute_col_back(col_permute_record);
		cout << "PM_preocessed" << PM_preocessed;
		cout << "PM - PM_preocessed" << PM - PM_preocessed;

		Matrix<my_double> recv(1, n, {
			2.369,-2.195    ,2.115,-1.774,-1.681    ,1.672     ,1.56    ,1.539    ,1.527    ,1.506,
			-1.479,-1.465,-1.431    ,1.431,-1.416,-1.359    ,1.347,-1.327,-1.326,-1.32    ,1.276  ,
			1.273    ,1.215,-1.213    ,1.196    ,1.194,-1.177    ,1.157,-1.137,-1.074    ,1.058,-1.045,
			0.9864,-0.9861,-0.9528    ,0.9297,-0.9147,-0.9084,-0.8679    ,0.8561,-0.8509,-0.8268,-0.7669,
			0.6123,-0.5966,-0.7832    ,0.709    ,0.7086    ,0.7037    ,0.5261,-0.4992
			});
		cout << "recv" << recv;

		vit.change_PM(PM);

		Matrix<GF2> decode_list = vit.decode_v_once(recv, 64);
		for (int i = 0, imax = decode_list.row(); i < imax; ++i) {
			my_double d_CL = Measure::Euclidean_distance(recv, BPSK::modulation(decode_list.get_row(i)));
			cout << "d_CL = " << d_CL << endl;
		}
	}

	static void BCH_BM_simulation_multi_SNR() {

		bool if_output_to_file = false;

		// encoding parameter
		const int m = 5;
		const int t = 1;			// set BCH parameter
		GF2e<m>::init();
		BCH<m, t> bch;				// here we choose BCH
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		Matrix<GF2> PM = bch.get_parity_matrix();
		//bch.print_info();

		my_double SNR_start = 3;
		my_double SNR_gap = 1;
		my_double SNR_end = 8;		// set SNR range

		Matrix<my_double> test_SNR(SNR_start, SNR_gap, SNR_end, 'd');
		int len = test_SNR.size();

		// set cout into file
		if (if_output_to_file) {
			char file_name[555] = { 0 };
			sprintf_s(file_name, 555, "BCH(%d,%d)_count_BM_simulation_SNR(%.2f-%.2f)dB.txt", n, k, (double)SNR_start, (double)SNR_end);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}


		bool is_first_time = true;
		for (int ii = 0; ii < len; ++ii) {

			my_double SNR_dB = test_SNR(ii);

			//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
			int simulation_times = 50000000;	// it is not until 10000 times did the operation conunt become stable, using pass_standard


			Matrix<GF2> u(1, k, 'b');
			Matrix<GF2> v = bch.encode(u);
			//cout << "v" << v;

			Matrix<my_double> c = BPSK::modulation(v);

			// channel
			my_double SNR = pow(10, SNR_dB / 10);
			//cout << "SNR=" << SNR << endl;
			my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
			AWGN::sigma = sigma;

			// performance recording
			int error_frame = 0;
			int error_frame_half_distance = 0;
			const int max_error_frame = 300;

			if (is_first_time) {

				bch.print_info();
				cout << "max_error_frame          =  " << max_error_frame << endl;
				cout << "simulation type          =  count_BM" << endl;

				cout << "\n------------------------------------------------------------------------------------------------------------"
					<< "--------------------------------------------------------------------------------" << endl;


				cout << setw(20) << setiosflags(ios::left) << "SNR_dB" << setw(20) << "simulation_times" 
					<< setw(20) << "error_rate" << setw(20) << "error_rate_hd";

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
				Matrix<GF2> hdr = BPSK::demodulation(recv);
				Matrix<GF2> v_hat = bch.decode_v(hdr);


				if (v != v_hat) {
					// record the error frame
					error_frame++;
					if (error_frame == max_error_frame) {
						simulation_times = i + 1;
						break;
					}
				}

				// errors over half distaance, check BM correctness
				//if ((hdr - v).Hamming_weight() > (d_min - 1) / 2) {		// Remind that this is (d-1)/2
				//	error_frame_half_distance++; 
				//	if (error_frame_half_distance == max_error_frame) {
				//		simulation_times = i + 1;
				//		break;
				//	}
				//}


				//cout << "simulation_times = " << i + 1 << '\r';
			}
			end = clock();

			double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
			double pass_AWGN_time = 2.85e-6 / 31 * n;
			cout << setw(20) << fixed << setprecision(1) << SNR_dB;
			cout << setw(20) << simulation_times;
			cout << scientific << setprecision(2);
			cout << setw(20) << error_frame / (double)simulation_times;
			cout << setw(20) << error_frame_half_distance / (double)simulation_times;


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
	}
	static void RS_BM_simulation_multi_SNR() {

		bool if_output_to_file = true;

		// encoding parameter
		const int m = 8;
		const int k = 128;			// set RS parameter
		GF2e<m>::init();
		RS<m, k> rs;				// here we choose RS
		const int n = (1 << m) - 1;
		const int d_min = n - k + 1;// MDS property
		Matrix<GF2e<m>> PM = rs.get_parity_matrix();
		//bch.print_info();

		my_double SNR_start = 6.3;
		my_double SNR_gap = 0.5;
		my_double SNR_end = 6.3;		// set SNR range

		Matrix<my_double> test_SNR(SNR_start, SNR_gap, SNR_end, 'd');
		int len = test_SNR.size();

		// set cout into file
		if (if_output_to_file) {
			char file_name[555] = { 0 };
			sprintf_s(file_name, 555, "RS(%d,%d)_count_BM_simulation_SNR(%0.2f-%0.2f)dB.txt", n, k, (double)SNR_start, (double)SNR_end);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}


		bool is_first_time = true;
		GF_trans<GF2e<m>, m> gft;		// For RS code transmitted to the AWGN channel
		for (int ii = 0; ii < len; ++ii) {

			my_double SNR_dB = test_SNR(ii);

			//cout << endl << "------------- SNR_dB = " << SNR_dB << " ------------" << endl;
			int simulation_times = 30000000;	// it is not until 10000 times did the operation conunt become stable, using pass_standard


			Matrix<GF2e<m>> u(1, k);
			for (int i = 0; i < k; ++i) {
				u(i) = my::rand_int(0, (1 << m) - 1);
			}
			//cout << "u" << u;

			Matrix<GF2e<m>> v_GF = rs.encode(u);
			Matrix<GF2> v = gft.to_bits(v_GF);	// turn GF2e<m> vector into binary vector

			//cout << "v" << v;

			Matrix<my_double> c = BPSK::modulation(v);

			// channel
			my_double SNR = pow(10, SNR_dB / 10);
			//cout << "SNR=" << SNR << endl;
			my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
			AWGN::sigma = sigma;

			// performance recording
			int error_frame = 0;
			int error_frame_half_distance = 0;
			const int max_error_frame = 300;

			if (is_first_time) {

				rs.print_info();
				cout << "max_error_frame          =  " << max_error_frame << endl;
				cout << "simulation type          =  BM" << endl;

				cout << "\n------------------------------------------------------------------------------------------------------------"
					<< "--------------------------------------------------------------------------------" << endl;

				// for output in excel
				cout << setw(20) << setiosflags(ios::left) << "SNR_dB" << setw(20) << "simulation_times" 
					<< setw(20) << "error_rate" << setw(20) << "error_rate_hd";

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
				Matrix<GF2> hdr = BPSK::demodulation(recv);
				Matrix<GF2e<m>> hdr_GF = gft.to_symbol(hdr);
				//cout << "hdr_GF" << hdr_GF;
				//cout << "v_GF - hdr_GF" << v_GF - hdr_GF;
				//cout << "(v_GF - hdr_GF).Hamming_weight() = " << (v_GF - hdr_GF).Hamming_weight() << endl;

				Matrix<GF2e<m>> v_hat_GF = rs.decode_BM_v(hdr_GF);

				if (v_GF != v_hat_GF) {
					// record the error frame
					error_frame++;
				}

				// errors over half distaance
				if ((v_GF - hdr_GF).Hamming_weight() > (d_min - 1) / 2) {	// Remind that this is (d-1)/2
					error_frame_half_distance++;

					if (error_frame_half_distance == max_error_frame) {
						simulation_times = i + 1;
						break;
					}
				}

				//cout << "simulation_times = " << i + 1 << '\r';
			}
			end = clock();

			double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
			double pass_AWGN_time = 2.85e-6 / 31 * n;
			
			// for convenience in output to Excel
			cout << setw(20) << fixed << setprecision(1) << SNR_dB;
			cout << setw(20) << simulation_times;
			cout << scientific << setprecision(2);
			cout << setw(20) << error_frame / (double)simulation_times;
			cout << setw(20) << error_frame_half_distance / (double)simulation_times;


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
	}
	static void RS_BM_test() {
		// encoding parameter
		const int m = 4;
		const int k = 5;			// set RS parameter
		GF2e<m>::init();
		RS<m, k> rs;				// here we choose RS
		const int n = (1 << m) - 1;
		const int d_min = n - k + 1;	// MDS property
		rs.print_info();

		Matrix<GF2e<m>> hdr(1, 15);
		hdr(0).set_by_alpha_power(3);
		hdr(1).set_by_alpha_power(9);
		hdr(2).set_by_alpha_power(13);
		hdr(3) = 0;
		hdr(4).set_by_alpha_power(6);
		hdr(5).set_by_alpha_power(0);
		hdr(6).set_by_alpha_power(5);
		hdr(7).set_by_alpha_power(2);
		hdr(8).set_by_alpha_power(7);
		hdr(9).set_by_alpha_power(4);
		hdr(10).set_by_alpha_power(13);
		hdr(11).set_by_alpha_power(1);
		hdr(12).set_by_alpha_power(8);
		hdr(13).set_by_alpha_power(2);
		hdr(14) = 0;

		cout << "hdr" << hdr;

		Matrix<GF2e<m>> v = rs.decode_BM_v(hdr);
		cout << "v" << v;

		cout << "rs.is_in_v_space(v) = " << rs.is_in_v_space(v) << endl;
	}

	static void RS_7_3_code_test() {
		typedef GF2e<3> ty;
		ty::init();

		RS<3, 3> rs;
		rs.print_info();
		int n = rs.get_n();
		int k = rs.get_k();
		int d = rs.get_d();

		Matrix<ty> u(1, k);
		u(0).set_by_alpha_power(3);
		u(1).set_by_alpha_power(5);
		u(2).set_by_alpha_power(0);

		cout << "u" << u;
		Matrix<ty> v = rs.encode(u);

		cout << "rs.is_in_v_space(v) = " << rs.is_in_v_space(v) << endl;
		/*Matrix<ty> v(1, n);
		v(0).set_by_alpha_power(5);
		v(1).set_by_alpha_power(4);
		v(2) = 0;
		v(3).set_by_alpha_power(0);
		v(4).set_by_alpha_power(4);
		v(5).set_by_alpha_power(0);
		v(6).set_by_alpha_power(5);*/

		cout << "v" << v;
		Matrix<ty> e(1, n, '0');
		e(3).set_by_alpha_power(1);
		e(4).set_by_alpha_power(6);
		Matrix<ty> r = e + v;
		/*Matrix<ty> r(1, n);
		r(0).set_by_alpha_power(5);
		r(1).set_by_alpha_power(4);
		r(2) = 0;
		r(3).set_by_alpha_power(0);
		r(4).set_by_alpha_power(4);
		r(5).set_by_alpha_power(2);
		r(6).set_by_alpha_power(1);*/

		cout << "r" << r;

		Matrix<ty> u_hat = rs.decode_BM(r);
		cout << "u_hat" << u_hat;
		cout << "u - u_hat" << u - u_hat;

		cout << "rs.get_generator_matrix()" << rs.get_generator_matrix();
		cout << "rs.get_parity_matrix()" << rs.get_parity_matrix();
		cout << "rs.get_generator_matrix() * rs.get_parity_matrix()" << rs.get_generator_matrix() * rs.get_parity_matrix().Transpose();
	}
	static void RS_7_5_code_test() {
		typedef GF2e<3> ty;
		ty::init();

		RS<3, 5> rs;
		rs.print_info();
		int n = rs.get_n();
		int k = rs.get_k();
		int d = rs.get_d();

		Matrix<ty> u(1, k);
		u(0).set_by_alpha_power(4);
		u(1).set_by_alpha_power(0);
		u(2).set_by_alpha_power(5);

		cout << "u" << u;
		Matrix<ty> v = rs.encode(u);

		cout << "rs.is_in_v_space(v) = " << rs.is_in_v_space(v) << endl;
		/*Matrix<ty> v(1, n);
		v(0).set_by_alpha_power(5);
		v(1).set_by_alpha_power(4);
		v(2) = 0;
		v(3).set_by_alpha_power(0);
		v(4).set_by_alpha_power(4);
		v(5).set_by_alpha_power(0);
		v(6).set_by_alpha_power(5);*/

		cout << "v" << v;
		Matrix<ty> e(1, n, '0');
		e(3).set_by_alpha_power(1);
		e(4).set_by_alpha_power(6);
		Matrix<ty> r = e + v;
		/*Matrix<ty> r(1, n);
		r(0).set_by_alpha_power(5);
		r(1).set_by_alpha_power(4);
		r(2) = 0;
		r(3).set_by_alpha_power(0);
		r(4).set_by_alpha_power(4);
		r(5).set_by_alpha_power(2);
		r(6).set_by_alpha_power(1);*/

		cout << "r" << r;

		Matrix<ty> u_hat = rs.decode_BM(r);
		cout << "u_hat" << u_hat;
		cout << "u - u_hat" << u - u_hat;

		cout << "rs.get_generator_matrix()" << rs.get_generator_matrix();
		cout << "rs.get_parity_matrix()" << rs.get_parity_matrix();
		cout << "rs.get_generator_matrix() * rs.get_parity_matrix()" << rs.get_generator_matrix() * rs.get_parity_matrix().Transpose();
	}
	static void RS_7_5_dual_code_test() {
		const int m = 3;
		const int k = 5;
		const int field_size = 1 << m;
		const int n = field_size - 1;
		const int k_perp = n - k;

		typedef GF2e<m> ty;
		ty::init();
		RS<m, k> rs;
		Matrix<ty> G = rs.get_generator_matrix();
		Matrix<ty> H = rs.get_parity_matrix();

		//Matrix<ty> u(1, k);
		//for (int i = 0; i < k; ++i) {
		//	u(i) = my::rand_int_adv(0, field_size - 1);
		//}
		//cout << "u" << u;
		//
		//Matrix<ty> c = rs.encode(u);
		//cout << "c" << c;

		RS<m, k_perp> rs_perp;
		Matrix<ty> G_perp = rs_perp.get_generator_matrix();
		Matrix<ty> H_perp = rs_perp.get_parity_matrix();

		//cout <<"c.multiply_transpose_of(G_perp)" <<  c.multiply_transpose_of(G_perp);

		cout << "G" << G;
		cout << "H" << H;
		cout << "G_perp" << G_perp;
		cout << "H_perp" << H_perp;

		cout << "G.multiply_transpose_of(H)" << G.multiply_transpose_of(H);
		cout << "G_perp.multiply_transpose_of(H_perp)" << G_perp.multiply_transpose_of(H_perp);
		cout << "G.multiply_transpose_of(G_perp)" << G.multiply_transpose_of(G_perp);
		cout << "H.multiply_transpose_of(H_perp)" << H.multiply_transpose_of(H_perp);

		// enumerate all codeword that generated by H, i.e., enumerate the dual code of rs
		Matrix<ty> u(1, k_perp, '0');

		int total_codewords = (int)round(pow(field_size, k_perp));
		int ii = 0;
		int min_non_zero_cnt = n;
		while (true) {

			if (ii == total_codewords - 1)
				break;

			// generate all possible vector in GF2e<m> of length k_perp
			for (int j = 0; j < k_perp; ++j) {
				u(j) = (int)u(j) + 1;
				if (((int)u(j)) < field_size) {
					break;
				}
				else {
					u(j) = 0;
				}
			}
			ii++;

			cout << "u" << u;

			// encode with u
			Matrix<ty> v = u * H;
			cout << "v" << v;

			// counting number of zeros in v
			int non_zero_cnt = 0;
			for (int i = 0; i < n; ++i) {
				if (v(i) != 0) {
					non_zero_cnt++;
				}
			}
			min_non_zero_cnt = my::min(min_non_zero_cnt, non_zero_cnt);
		}

		cout << "min_non_zero_cnt = " << min_non_zero_cnt << endl;
	}
	static void linear_block_code_6_3_3_test() {
		Matrix<GF2> G(3, 6, {
			1,0,0,0,1,1,
			0,1,0,1,0,1,
			0,0,1,1,1,0
			});
		Matrix<GF2> H(3, 6, {
			0,1,1,1,0,0,
			1,0,1,0,1,0,
			1,1,0,0,0,1
			});

		cout << "G" << G;
		cout << "H" << H;
		cout << "G.multiply_transpose_of(H)" << G.multiply_transpose_of(H);

		G.row_transformation_to_low_triangle();
		H.row_transformation_to_up_triangle();


		cout << "G" << G;
		cout << "H" << H;
		cout << "G.multiply_transpose_of(H)" << G.multiply_transpose_of(H);
	}
	static void Hamming_code_7_4_3_test() {
		GF2e<3>::init();
		BCH<3, 1> bch;
		int k = bch.get_k();
		int n = bch.get_n();
		int r = bch.get_r();
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();

		cout << "G" << G;
		cout << "H" << H;
		cout << "G.multiply_transpose_of(H)" << G.multiply_transpose_of(H);

		Matrix<GF2> G_new = G.get_part(0, 0, k - 2, n - 2);
		Matrix<GF2> H_new = H.get_part(0, 0, -1, n - 2);
		for (int i = 0; i < n - 1; ++i) {
			H_new(1, i) += H_new(2, i);
		}

		cout << "G_new" << G_new;
		cout << "H_new" << H_new;
		cout << " G_new.multiply_transpose_of(H_new)" << G_new.multiply_transpose_of(H_new);
	}
	static void GF2t_test() {
		const int m = 4;
		GF2e<4>::init();
		int r = 4;
		int c = 4;
		Matrix<int> A(r, c);
		Matrix<GF2e<m>> B(r, c);		// in final, we only keep GF2e !!
		//Matrix<GF2t<m>> C(r, c);

		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
				A(i, j) = my::rand_int_adv(0, (1 << m) - 1);

				// set matrix of GF2e<m>
				B(i, j) = A(i, j);

				// set matrix of GF2t<m>
				//C(i, j).set_by_alpha_power(A(i, j) - 1);
			}
		}

		cout << "A" << A;
		cout << "B" << B;
		//cout << "C" << C;

		cout << "B*B" << B * B;
		//cout << "C*C" << C * C;

		/*int var_1 = 13, var_2 = 15;
		GF2e<m> x(var_1);
		GF2e<m> y(var_2);
		cout << "x-y = " << x - y << endl;

		GF2t<m> w;
		w.set_by_alpha_power(var_1 - 1);
		GF2t<m> z;
		z.set_by_alpha_power(var_2 - 1);
		cout << "w-z = " << w - z << endl;*/
		/*
		for (int i = 0; i < 16; ++i) {
			cout << "GF2e_auxiliary_storage::alpha_table[" << i << "]: " << GF2e_auxiliary_storage::alpha_table[i] << endl;
		}

		cout << "-------------" << endl;

		for (int i = 0; i < 15; ++i) {
			cout << "GF2e_auxiliary_storage::polynomial_table[" << i << "]: " << GF2e_auxiliary_storage::polynomial_table[i] << endl;
		}*/

		/*for (int i = 0; i < 16; ++i) {
			cout << "i=" << i << ",\trev(i)=" << my::rev_bits(i, 4) << endl;
		}*/

		cout << "B.inv()" << B.inv();
		//cout << "C.inv()" << C.inv();

		cout << "B * B.inv()" << B * B.inv();
		//cout << "C * C.inv()" << C * C.inv();
	}
	static void GF2e_to_bits_row_extention() {
		const int m = 4;
		const int pow_m = 1 << m;
		GF2e<m>::init();
		Matrix<GF2e<m>> orig(1, pow_m, 'N');
		cout << "orig" << orig;

		Matrix<GF2> M;
		M = GF2e<m>::to_bits_row_extention(orig);
		cout << "M" << M;
	}
	static void GF2e_trace_test() {
		const int m = 3;
		const int k = 5;		// RS(7,5) code d = 3
		const int t = 1;		// BCH(7,4) code, error correcting capability = 1
		const int field_size = 1 << m;
		const int n = field_size - 1;

		typedef GF2e<m> ty;
		ty::init();

		for (int i = 0; i < field_size; ++i) {
			cout << "----------------" << endl;
			ty a = i;
			GF2 b = a.trace();
			cout << "a = " << a << endl;
			cout << "b = " << b << endl;
		}

		RS<m, k> rs;
		Matrix<ty> G = rs.get_generator_matrix();
		Matrix<ty> H = rs.get_parity_matrix();

		cout << "G" << G;
		cout << "H" << H;
		cout << "G.multiply_transpose_of(H)" << G.multiply_transpose_of(H);

		// extend G to GF2 matrix and do Gaussian elemination to get GF2 subfield subcode Generator

		Matrix<GF2> H2 = ty::to_bits_row_extention(H);
		cout << "H2" << H2;

		H2.row_transformation_to_up_triangle();
		cout << "H2" << H2;

		BCH<m, t> bch;
		int k_bch = bch.get_k();
		Matrix<GF2> G_bch = bch.get_generator_matrix();
		Matrix<GF2> H_bch = bch.get_parity_matrix();
		cout << "G_bch" << G_bch;
		// prove G_bch is in the subspace of G, BCH code is subfield subcode of RS code
		Matrix<ty> G_bch_ext(k_bch, n);		// GF2e<m> as element
		for (int i = 0, imax = G_bch_ext.size(); i < imax; ++i) {
			G_bch_ext(i) = G_bch(i);
		}
		cout << "G_bch_ext" << G_bch_ext;
		cout << "G_bch_ext.multiply_transpose_of(H)" << G_bch_ext.multiply_transpose_of(H);
		// zero matrix indicates that BCH(7,4) is a subfield subcode of RS(7,5)

		cout << "H_bch" << H_bch;

		// prove BCH dual code is subfield subcode of RS dual code
		int k_bch_dual = n - k_bch;
		H2 = H2.get_part(0, 0, k_bch_dual - 1, -1);		// erase zero rows in H2
		cout << "H2 without zeros" << H2;

		cout << "H2.multiply_transpose_of(G_bch)" << H2.multiply_transpose_of(G_bch);

		// to see if we can turn an RS codeword into its subfield subcode, not ture
		Matrix<ty> u(1, k);
		for (int i = 0; i < k; ++i) {
			u(i) = my::rand_int(0, field_size - 1);
		}
		cout << "u" << u;
		Matrix<ty> v = rs.encode(u);
		cout << "v" << v;

		Matrix<GF2> v_tr(1, n);
		for (int i = 0; i < n; ++i) {
			v_tr(i) = v(i).trace();
		}
		cout << "v_tr" << v_tr;
		cout << "v_tr.multiply_transpose_of(H_bch)" << v_tr.multiply_transpose_of(H_bch);		// not ture

		// to see if we can turn an RS dual codeword into its subfield subcode, i.e. dual BCH code, true
		int k_perp = n - k;
		Matrix<ty> u_perp(1, k_perp);
		for (int i = 0; i < k_perp; ++i) {
			u_perp(i) = my::rand_int(0, field_size - 1);
		}
		cout << "u_perp" << u_perp;
		Matrix<ty> v_perp = rs.encode_dual(u_perp);
		cout << "v_perp" << v_perp;

		Matrix<GF2> v_perp_tr(1, n);
		for (int i = 0; i < n; ++i) {
			v_perp_tr(i) = v_perp(i).trace();
		}
		cout << "v_perp_tr" << v_perp_tr;
		cout << "v_perp_tr.multiply_transpose_of(G_bch)" << v_perp_tr.multiply_transpose_of(G_bch);

		// generate dual BCH code from all dual RS code
		u_perp.reset(0);

		int total_codewords = (int)round(pow(field_size, k_perp));
		int ii = 0;
		int min_non_zero_cnt = n;
		Matrix<GF2> dual_BCH_G_from_codeword(total_codewords, n, '0');
		while (true) {

			if (ii == total_codewords - 1)
				break;

			// generate all possible vector in GF2e<m> of length k_perp
			for (int j = 0; j < k_perp; ++j) {
				u_perp(j) = (int)u_perp(j) + 1;
				if (((int)u_perp(j)) < field_size) {
					break;
				}
				else {
					u_perp(j) = 0;
				}
			}
			ii++;

			//cout << "u_perp" << u_perp;

			// encode with u
			Matrix<ty> v_perp = u_perp * H;
			//cout << "v_perp" << v_perp;

			// get subfield subcode by trace function
			for (int i = 0; i < n; ++i) {
				v_perp_tr(i) = v_perp(i).trace();
			}
			cout << "ii = " << ii << endl;
			cout << "v_perp_tr" << v_perp_tr;
			//cout << "v_perp_tr.multiply_transpose_of(G_bch)" << v_perp_tr.multiply_transpose_of(G_bch);	// test if it is a valid dual bch code

			dual_BCH_G_from_codeword.set_row(ii, v_perp_tr);

			// counting number of zeros in v_perp_tr
			int non_zero_cnt = 0;
			for (int i = 0; i < n; ++i) {
				if (v_perp_tr(i) != 0) {
					non_zero_cnt++;
				}
			}
			if (non_zero_cnt == 0) {
				non_zero_cnt = n;
			}

			min_non_zero_cnt = my::min(min_non_zero_cnt, non_zero_cnt);
		}

		// turn the set into a dual BCH generator

		cout << "min_non_zero_cnt (v_perp_tr) = " << min_non_zero_cnt << endl;
		cout << "dual_BCH_G_from_codeword" << dual_BCH_G_from_codeword;
		dual_BCH_G_from_codeword.row_transformation_to_up_triangle();
		cout << "dual_BCH_G_from_codeword" << dual_BCH_G_from_codeword;
		cout << "dual_BCH_G_from_codeword.multiply_transpose_of(G_bch)" << dual_BCH_G_from_codeword.multiply_transpose_of(G_bch);
		dual_BCH_G_from_codeword.resize(k_bch_dual, n, true);
		dual_BCH_G_from_codeword.col_permute_to_full_rank_on_left();
		dual_BCH_G_from_codeword.row_transformation_left_up_triangle_to_identity();
		cout << "dual_BCH_G_from_codeword" << dual_BCH_G_from_codeword;
		cout << "dual_BCH_G_from_codeword.multiply_transpose_of(G_bch)" << dual_BCH_G_from_codeword.multiply_transpose_of(G_bch);
	}
	static void BCH_dual_test() {
		const int m = 3;
		const int t = 1;		// BCH(15,11) code, error correcting capability = 1
		const int field_size = 1 << m;
		const int n = field_size - 1;

		typedef GF2e<m> ty;
		ty::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		int k_perp = bch.get_r();

		polynomial<GF2> gX = bch.get_gX2();
		polynomial<GF2> pX = bch.get_pX2();

		cout << "gX" << gX;
		cout << "pX" << pX;

		// the distance of dual code is unknown, we search the codeword to get the weight distribution

		// the dual code dimention is k_perp
		Matrix<GF2> u_perp(1, k_perp, '0');

		const int GF2_size = 2;
		int total_codewords = (int)round(pow(GF2_size, k_perp));
		int ii = 0;
		int min_non_zero_cnt = n;
		while (true) {

			if (ii == total_codewords - 1)
				break;

			// generate all possible vector in GF2e<m> of length k_perp
			for (int j = 0; j < k_perp; ++j) {
				u_perp(j) = (int)u_perp(j) + 1;
				if (((int)u_perp(j)) < GF2_size) {
					break;
				}
				else {
					u_perp(j) = 0;
				}
			}
			ii++;

			cout << "u_perp" << u_perp;

			// encode with u
			Matrix<GF2> v = bch.encode_dual(u_perp);
			cout << "v" << v;

			// counting number of zeros in v
			int non_zero_cnt = 0;
			for (int i = 0; i < n; ++i) {
				if (v(i) != 0) {
					non_zero_cnt++;
				}
			}
			cout << "non_zero_cnt = " << non_zero_cnt << endl;
			min_non_zero_cnt = my::min(min_non_zero_cnt, non_zero_cnt);
		}
		cout << "min_non_zero_cnt = " << min_non_zero_cnt << endl;
	}
	static void GF2e_trace_test_only() {
		const int m = 9;
		const int field_size = 1 << m;
		typedef GF2e<m> ty;
		ty::init();

		for (int i = 0; i < field_size; ++i) {
			cout << "----------------" << endl;
			ty a = i;
			GF2 b = a.trace();
			cout << "a = " << a << endl;
			cout << "b = " << b << endl;
		}
	}
};
