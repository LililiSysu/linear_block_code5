#pragma once
/*****************************************************************//**
 * \file   Matrix_analysis.h
 * \brief  for matrix decompsition algorithm testing mainly
 *
 * \author 26259
 * \date   March 2023
 *********************************************************************/

#include"../GF/polynomial.h"
using namespace std;

class Eigen_decomposition_of_block_tridiagonal_Matrices {
public:

	/**
	 * . A should be square tridiagonal matrix, k is the block size, eigen_value and eigen_vector will be returned
	 */
	template <class T>
	static void solve(const Matrix<T>& A, const int k, Matrix<T>& eigen_value, Matrix<T>& eigen_vector) {
		const int n = A.col();
		const int L = n / k;

		Matrix<polynomial<T>> A_to_polynomials(n, n);
		Matrix_polynomials::to(A, A_to_polynomials);

		// get k*k block D[0,L-2], B[0,L-1], C[0,L-2]
		Matrix<Matrix<T>> D(1, L - 1);
		Matrix<polynomial<T>> Di(k, k);
		Matrix<Matrix<polynomial<T>>> B(1, L);
		Matrix<Matrix<polynomial<T>>> C(1, L - 1);

		// get the block over the Matrix
		B(0) = A_to_polynomials.get_part(0, 0, k - 1, k - 1);
		int ik = 0;
		for (int i = 0; i < L - 1; ++i) {		
			int i1_k = ik + k;
			int i2_k = i1_k + k;
			D(i) = A.get_part(ik, i1_k, i1_k - 1, i2_k - 1);
			B(i + 1) = A_to_polynomials.get_part(i1_k, i1_k, i2_k - 1, i2_k - 1);
			C(i) = A_to_polynomials.get_part(i1_k, ik, i2_k - 1, i1_k - 1);
			ik += k;
		}

		// polynomial matrix index from in [0,L]
		Matrix<Matrix<polynomial<T>>> P(1, L + 1);
		// P(-1)=0, doesnot need to store
		P(0) = Matrix<polynomial<T>>(k, k, 'i');
		polynomial<T> x(Matrix<T>(1, 2, { 0,1 }));

		Matrix_polynomials::to(D(0).inv(), Di);

		// starting recurrance
		P(1) = Di * (x * P(0) - B(0) * P(0));

		// recurrance
		for (int i = 1; i < L - 1; ++i) {


			Matrix_polynomials::to(D(i).inv(), Di);

			P(i + 1) = Di * (x * P(i) - B(i) * P(i) - C(i - 1) * P(i - 1));

		}

		// ending recurrance

		// D(L-1)=I, the identity matrix		
		P(L) = x * P(L - 1) - B(L - 1) * P(L - 1) - C(L - 2) * P(L - 2);

		// compute det P(L), using Gaussian elimination
		polynomial<T> det_PL = P(L).det_4_matrix_polynomials();

		// find the roots
		Matrix<T> R(n, n, '0');
		int n_minus_1 = n - 1;
		for (int i = 0; i < n_minus_1; ++i) {
			R(i + 1, i) = 1;
		}
		for (int i = 0; i < n; ++i) {
			R(0, i) = -det_PL(n - 1 - i) / det_PL(n);
		}
		Matrix<T> eval;
		R.ritz_Hessenberg(eval);		// you can adjust the precision, to get accurate result

		// now get the eigen vectors
		Matrix<T> PL_lambda(k, k);
		Matrix<T> orth;
		Matrix<T> u;
		eigen_vector.resize(n, n, false);
		int eval_size = eval.size();	// consider the case of non-deficient
		for (int eval_ind = 0; eval_ind < eval_size; ) {

			// plug in eigen value into P(L), to get a k by k defficient Matrix, where its null space can generate eigen vectors
			Matrix_polynomials::from(P(L), eval(eval_ind), PL_lambda);

			// find the null space of PL_lambda
			orth = PL_lambda.Schimidt_orthogonalization(false);
			orth = orth.erase_rows(orth.find_zero_row());
			int orth_row_dim = orth.row();
			orth = orth.add_orthogonal_complementary();		// the null space is added to the bottom of the orth

			// get the null space
			u = orth.get_part(orth_row_dim, 0, -1, -1).Hermit();		// use Hermit for general case!
			int u_col = u.col();							// there may be more than one eigen vector corresponding to one eigen values

			for (int i = 0; i <= L - 1; ++i) {				// compute eigen vectors part by part and set them into result evec
				Matrix_polynomials::from(P(i), eval(eval_ind), PL_lambda);
				eigen_vector.set_part(i * k, eval_ind, PL_lambda * u);
			}

			eval_ind += u_col;								// skip the identical eigen values
		}

		eigen_vector.unitize_col();
		eigen_value = eval.place_diag();
	}
};

class Matrix_analysis {
public:

	static void polybomial_determinant_test() {
		typedef Complex ty;

		Matrix<polynomial<ty>> A(2, 2);
		A(0, 0) = Matrix<ty>(1, 2, { 1,1 });
		A(0, 1) = Matrix<ty>(1, 2, { -5,1 });
		A(1, 0) = Matrix<ty>(1, 2, { -2,3 });
		A(1, 1) = Matrix<ty>(1, 2, { 0,6 });

		cout << "A" << A;
		polynomial<ty> det_A = A.det_4_matrix_polynomials();
		cout << "det_A" << det_A;		// desired: {-10,23,3}
	}
	static void polynomial_matrix_test() {
		typedef Complex ty;
		const int k = 2;
		const int n = 6;
		const int L = 3;

		Matrix<ty> A(n, n, {
			3,	7,	3,	8,	0,	0,
			0,	4,	5,	0,	0,	0,
			1,	4,	2,	5,	-5,	-7,
			9,	-1,	4,	4,	2,	7,
			0,	0,	-2,	6,	2,	10,
			0,	0,	4,	1,	-3,	6
			});
		/*Matrix<ty> A(n, n, {
			1,	0,	1,	0,	0,	0,
			0,	1,	0,	1,	0,	0,
			1,	0,	1,	0,	1,	0,
			0,	1,	0,	1,	0,	1,
			0,	0,	1,	0,	1,	0,
			0,	0,	0,	1,	0,	1
			});*/

		cout << "A" << A;
		cout << "A.rank() = " << A.rank() << endl;
		Matrix<polynomial<ty>> A_to_polynomials(n, n);
		Matrix_polynomials::to(A, A_to_polynomials);

		// get k*k block D[0,L-2], B[0,L-1], C[0,L-2]
		Matrix < Matrix<ty>> D(1, L - 1);
		Matrix<polynomial<ty>> Di(k, k);
		Matrix < Matrix<polynomial<ty>>> B(1, L);
		Matrix < Matrix<polynomial<ty>>> C(1, L - 1);
		B(0) = A_to_polynomials.get_part(0, 0, k - 1, k - 1);
		int ik = 0;
		for (int i = 0; i < L - 1; ++i) {
			int i1_k = ik + k;
			int i2_k = i1_k + k;
			D(i) = A.get_part(ik, i1_k, i1_k - 1, i2_k - 1);
			B(i + 1) = A_to_polynomials.get_part(i1_k, i1_k, i2_k - 1, i2_k - 1);
			C(i) = A_to_polynomials.get_part(i1_k, ik, i2_k - 1, i1_k - 1);
			ik += k;
		}

		/*cout << "B" << B;
		cout << "C" << C;
		cout << "D" << D;*/

		// polynomial matrix index from in [0,L]
		Matrix<Matrix<polynomial<ty>>> P(1, L + 1);
		// P(-1)=0, doesnot need to store
		P(0) = Matrix<polynomial<ty>>(k, k, 'i');
		polynomial<ty> x(Matrix<ty>(1, 2, { 0,1 }));

		/*cout << "x" << x;
		cout << "P(0)" << P(0);
		cout << "B(0)" << B(0);
		cout << "D(0)" << D(0);*/

		//cout << "x * P(0)" << x * P(0);
		//cout << "B(0) * P(0)" << B(0) * P(0);
		//cout << "(x * P(0) - B(0) * P(0))" << (x * P(0) - B(0) * P(0));
		//cout << "D(0).inv()" << D(0).inv();
		Matrix_polynomials::to(D(0).inv(), Di);
		//cout << "Di" << Di;

		// starting recurrance
		P(1) = Di * (x * P(0) - B(0) * P(0));

		//cout << "(x * P(0) - B(0) * P(0))" << (x * P(0) - B(0) * P(0));
		//cout << "P(1)" << P(1);


		for (int i = 1; i < L - 1; ++i) {
			/*cout << "P(" << i << ")" << P(i);
			cout << "B(" << i << ") * P(" << i << ")" << B(i) * P(i);
			cout << "C(" << i - 1 << ") * P(" << i - 1 << ")" << C(i - 1) * P(i - 1);
			cout << "B(" << i << ") * P(" << i << ") + C(" << i - 1 << ") * P(" << i - 1 << ")" << B(i) * P(i)+ C(i - 1) * P(i - 1);
			cout << "x*P(" << i << ")" << x * P(i);
			cout << "x * P(" << i << ") - B(" << i << ") * P(" << i << ") - C(" << i - 1 << ") * P(" << i - 1 << ")"
				<< x * P(i) -  B(i) * P(i) - C(i - 1) * P(i - 1);*/


			Matrix_polynomials::to(D(i).inv(), Di);

			/*cout << "Di" << Di;
			cout << "Di * (x * P(i) - B(i) * P(i) - C(i - 1) * P(i - 1))"
				<< Di * (x * P(i) - B(i) * P(i) - C(i - 1) * P(i - 1));*/

			P(i + 1) = Di * (x * P(i) - B(i) * P(i) - C(i - 1) * P(i - 1));

			//cout << "P(" << i + 1 << ")" << P(i + 1);
		}

		// ending recurrance

		//cout << "P(" << L-1 << ")" << P(L-1);
		//cout << "B(" << L-1 << ") * P(" << L-1 << ")" << B(L-1) * P(L-1);
		//cout << "C(" << L-1 - 1 << ") * P(" << L-1 - 1 << ")" << C(L-1 - 1) * P(L-1 - 1);
		//cout << "B(" << L-1 << ") * P(" << L-1 << ") + C(" << L-1 - 1 << ") * P(" << L-1 - 1 << ")" << B(L-1) * P(L-1) + C(L-1 - 1) * P(L-1 - 1);
		//cout << "x*P(" << L-1 << ")" << x * P(L-1);
		//cout << "x * P(" << L - 1 << ") - B(" << L - 1 << ") * P(" << L - 1 << ") - C(" << L - 1 - 1 << ") * P(" << L - 1 - 1 << ")"
		//	<< x * P(L - 1) - B(L - 1) * P(L - 1) - C(L - 1 - 1) * P(L - 1 - 1);

		// D(L-1)=I, the identity matrix		
		P(L) = x * P(L - 1) - B(L - 1) * P(L - 1) - C(L - 2) * P(L - 2);

		/*for (int i = 0; i <= L; ++i) {
			cout << "P(" << i << ")" << P(i);
		}*/

		polynomial<ty> real_det = P(L)(0, 0) * P(L)(1, 1) - P(L)(0, 1) * P(L)(1, 0);
		//cout << "real_det" << real_det;					// match with det_PL

		// compute det P(L), using Gaussian elimination
		polynomial<ty> det_PL = P(L).det_4_matrix_polynomials();
		//cout << "det_PL" << det_PL;

		// verify
		polynomial<ty> desired_det(7, '0');
		desired_det(0) = 1;
		desired_det.simplify();		// make sure it is simplified during computation
		desired_det *= Matrix<ty>(1, 2, { Complex(8.0942, +0.0000) ,1 });		// transform Matrix into polynomial implicitly
		desired_det *= Matrix<ty>(1, 2, { Complex(-14.2751, +0.0000) ,1 });
		desired_det *= Matrix<ty>(1, 2, { Complex(-6.0613, +6.6387) ,1 });
		desired_det *= Matrix<ty>(1, 2, { Complex(-6.0613, -6.6387) ,1 });
		desired_det *= Matrix<ty>(1, 2, { Complex(-1.3482, -0.7949) ,1 });
		desired_det *= Matrix<ty>(1, 2, { Complex(-1.3482, +0.7949) ,1 });

		desired_det *= det_PL(0) / desired_det(0);		// normalize

		//cout << "desired_det" << desired_det;			// this matches with det_PL

		// find the roots
		Matrix<ty> R(n, n, '0');
		int n_minus_1 = n - 1;
		for (int i = 0; i < n_minus_1; ++i) {
			R(i + 1, i) = 1;
		}
		for (int i = 0; i < n; ++i) {
			R(0, i) = -det_PL(n - 1 - i) / det_PL(n);
		}
		cout << "R" << R;
		Matrix<ty> eval;
		R.ritz(eval, false);		// this is costly, we lower the accuracy. How to converge faster?
		cout << "eval" << eval;

		// now get the eigen vectors
		Matrix<ty> PL_lambda(k, k);
		Matrix<ty> orth;
		Matrix<ty> u;
		Matrix<ty> evec(n, n, '0');
		int col_set = 0;
		int eval_size = eval.size();
		Matrix<ty> real_eig_val(1, n, 'v');		// consider the case of non-deficient
		for (int eval_ind = 0; eval_ind < eval_size; ++eval_ind) {

			cout << "---------------------------------------------------" << endl;
			cout << "eval_ind = " << eval_ind << ", eigen value = " << eval(eval_ind) << endl;
			cout << "---------------------------------------------------" << endl;
			Matrix_polynomials::from(P(L), eval(eval_ind), PL_lambda);
			//cout << "PL_lambda" << PL_lambda;

			orth = PL_lambda.Schimidt_orthogonalization(false, 1e-4);
			//cout << "orth: Schimidt_orthogonalization" << orth;
			orth = orth.erase_rows(orth.find_zero_row());
			//cout << "orth: erase_rows" << orth;
			int orth_row_dim = orth.row();
			orth = orth.add_orthogonal_complementary();
			/*cout << "roth: add_orthogonal_complementary" << orth;
			cout << "orth * orth.Transpose()" << orth * orth.Transpose();
			cout << "orth * orth.Hermit()" << orth * orth.Hermit();*/

			u = orth.get_part(orth_row_dim, 0, orth.row() - 1, orth.col() - 1).Hermit();		// use Hermit for general case!
			int u_col = u.col();							// there may be more than one eigen vector corresponding to one eigen values
			/*cout << "u" << u;
			cout << "PL_lambda * u" << PL_lambda * u;*/

			for (int i = 0; i <= L - 1; ++i) {
				Matrix_polynomials::from(P(i), eval(eval_ind), PL_lambda);
				/*cout << "PL_lambda" << PL_lambda;
				cout << "PL_lambda * u" << PL_lambda * u;*/
				evec.set_part(i * k, col_set, PL_lambda * u);
				//cout << "v" << v;
			}
			for (int i = 0; i < u_col; ++i) {
				real_eig_val.push_back(eval(eval_ind));		// for multiple identical eigen values
				evec.unitize_col(col_set + i);
			}
			// column normalization, unnecessary

			col_set += u_col;

			cout << "evec" << evec;
		}
		//cout << "v * v.inv()" << v * v.inv();
		cout << "eval.place_diag()" << eval.place_diag();
		cout << "v * eval.place_diag() * v.inv()" << evec * eval.place_diag() * evec.inv();		// finishing... 

		Matrix<ty> diff = evec * eval.place_diag() * evec.inv() - A;
		cout << "diff" << diff;

		cout << "diff.norm_F() / A.norm_F() =" << diff.norm_F() / A.norm_F() << endl;


		Matrix<ty> rval;
		Matrix<ty> rvec;
		A.ritz(rval, rvec, false);		// unfinished
		cout << "rval" << rval;
		cout << "rvec" << rvec;
		cout << "v * eval.place_diag() * v.inv()" << rvec * rval.place_diag() * rvec.inv();

		Matrix<ty> diff2 = rvec * rval.place_diag() * rvec.inv() - A;
		cout << "diff2" << diff2;

		cout << "diff2.norm_F() / A.norm_F() =" << diff2.norm_F() / A.norm_F() << endl;
	}
	static void polynomial_matrix_test_simple() {
		typedef Complex ty;
		const int k = 1;
		const int n = 4;
		const int L = 4;

		Matrix<ty> A(n, n, {
			4,22,0,0,
			-2,9,-6,0,
			0,5,1,-2,
			0,0,5,13
			});
		cout << "A" << A;
		cout << "A.rank() = " << A.rank() << endl;
		Matrix<polynomial<ty>> A_to_polynomials(n, n);
		Matrix_polynomials::to(A, A_to_polynomials);

		// get k*k block D[0,L-2], B[0,L-1], C[0,L-2]
		Matrix < Matrix<ty>> D(1, L - 1);
		Matrix<polynomial<ty>> Di_inv(k, k);
		Matrix < Matrix<polynomial<ty>>> B(1, L);
		Matrix < Matrix<polynomial<ty>>> C(1, L - 1);
		B(0) = A_to_polynomials.get_part(0, 0, k - 1, k - 1);
		int ik = 0;
		for (int i = 0; i < L - 1; ++i) {
			int i1_k = ik + k;
			int i2_k = i1_k + k;
			D(i) = A.get_part(ik, i1_k, i1_k - 1, i2_k - 1);
			B(i + 1) = A_to_polynomials.get_part(i1_k, i1_k, i2_k - 1, i2_k - 1);
			C(i) = A_to_polynomials.get_part(i1_k, ik, i2_k - 1, i1_k - 1);
			ik += k;
		}

		cout << "B" << B;
		cout << "C" << C;
		cout << "D" << D;

		// polynomial matrix index from in [0,L]
		Matrix<Matrix<polynomial<ty>>> P(1, L + 1);
		// P(-1)=0, doesnot need to store
		P(0) = Matrix<polynomial<ty>>(k, k, 'i');
		polynomial<ty> x(Matrix<ty>(1, 2, { 0,1 }));

		cout << "x" << x;
		cout << "P(0)" << P(0);
		cout << "B(0)" << B(0);
		cout << "D(0)" << D(0);

		//cout << "x * P(0)" << x * P(0);
		//cout << "B(0) * P(0)" << B(0) * P(0);
		//cout << "(x * P(0) - B(0) * P(0))" << (x * P(0) - B(0) * P(0));
		//cout << "D(0).inv()" << D(0).inv();
		Matrix_polynomials::to(D(0).inv(), Di_inv);
		cout << "Di" << Di_inv;

		// starting recurrance
		P(1) = Di_inv * (x * P(0) - B(0) * P(0));

		cout << "(x * P(0) - B(0) * P(0))" << (x * P(0) - B(0) * P(0));
		cout << "P(1)" << P(1);

		for (int i = 1; i < L - 1; ++i) {
			cout << "P(" << i << ")" << P(i);
			cout << "x*P(" << i << ")" << x * P(i);
			cout << "B(" << i << ") * P(" << i << ")" << B(i) * P(i);
			cout << "C(" << i - 1 << ") * P(" << i - 1 << ")" << C(i - 1) * P(i - 1);


			Matrix_polynomials::to(D(i).inv(), Di_inv);
			P(i + 1) = Di_inv * (x * P(i) - B(i) * P(i) - C(i - 1) * P(i - 1));

			cout << "P(" << i + 1 << ")" << P(i + 1);
		}

		// ending recurrance
		P(L) = x * P(L - 1) - B(L - 1) * P(L - 1) - C(L - 2) * P(L - 2);
		cout << "P(" << L << ")" << P(L);

		// compute det P(L), use Gaussian elimination
		polynomial<ty> det_PL = P(L).det_4_matrix_polynomials();
		cout << "det_PL" << det_PL;
		polynomial<ty> real_det = P(L)(0, 0)/* * P(L)(1, 1) - P(L)(0, 1) * P(L)(1, 0)*/;		// P is a 2 by 2 matrix
		cout << "real_det" << real_det;

		// verify
		polynomial<ty> desired_det(7, '0');
		desired_det(0) = 1;
		desired_det.simplify();		// make sure it is simplified during computation
		desired_det *= Matrix<ty>(1, 2, { -Complex(5.9810, +8.1181) ,1 });		// transform Matrix into polynomial implicitly
		desired_det *= Matrix<ty>(1, 2, { -Complex(5.9810, -8.1181) ,1 });
		desired_det *= Matrix<ty>(1, 2, { -Complex(2.7132, +0.0000) ,1 });
		desired_det *= Matrix<ty>(1, 2, { -Complex(12.3247, +0.0000) ,1 });

		desired_det *= det_PL(0) / desired_det(0);		// normalize
		cout << "desired_det" << desired_det;
	}
	static void add_orthogonal_complementary_test() {
		typedef Complex ty;
		Matrix<ty> v1(1, 2, { Complex(0.8103),Complex(0.3905,-0.4368) });
		Matrix<ty> Q = v1.add_orthogonal_complementary();
		cout << "Q" << Q;
		cout << "Q * Q.Transpose()" << Q * Q.Transpose();
		cout << "Q * Q.Hermit()" << Q * Q.Hermit();
		cout << "Q.get_row(0) * Q.get_row(1).Transpose()" << Q.get_row(0) * Q.get_row(1).Transpose();
		cout << "Q.get_row(0) * Q.get_row(1).Hermit()" << Q.get_row(0) * Q.get_row(1).Hermit();
	}
	static void ritz_test() {
		typedef Complex ty;
		const int n = 6;
		// this matrix has two identical eigen values. 4
		Matrix<ty> A(n, n, {
			{1.5868 } ,{   1.1197 } ,{   0.6456 } ,{   0.4152 } ,{   2.4565 } ,{   1.5618 } ,
			{ -1.9699  },{   4.6478 } ,{   3.5547 } ,{   1.1028 } ,{ -0.8435 } ,{ -0.5575 } ,
			{   3.0900 } ,{ -0.1346 } ,{   2.8177 } ,{ -0.7208 } ,{ -2.4416 } ,{ -5.3426 } ,
			{   1.1126 } ,{ -2.2170 } ,{ -2.5785 } ,{   3.3851 } ,{   0.4612 } ,{   5.9459 } ,
			{   0.5622 } ,{ -4.6417 } ,{ -2.5129 } ,{ -0.1981 } ,{   4.3912 } ,{  11.5370 } ,
			{  -0.4229 } ,{ -0.2169 } ,{   3.6158 } ,{   0.7668 } ,{ -3.0093 } ,{   5.9904} ,
			});
		cout << "A" << A;

		Matrix<ty> rval;
		Matrix<ty> rvec;
		A.ritz(rval, rvec, false);
		cout << "rval" << rval;
		cout << "rvec" << rvec;
		cout << "v * eval.place_diag() * v.inv()" << rvec * rval.place_diag() * rvec.inv();

		Matrix<ty> diff2 = rvec * rval.place_diag() * rvec.inv() - A;
		cout << "diff2" << diff2;

		cout << "diff2.norm_F() / A.norm_F() =" << diff2.norm_F() / A.norm_F() << endl;
	}
	static void ritz_inner_test() {
		typedef Complex ty;
		const int n = 4;
		// verify
		polynomial<ty> desired_det(7, '0');
		desired_det(0) = 1;
		desired_det.simplify();		// make sure it is simplified during computation
		desired_det *= Matrix<ty>(1, 2, { Complex(2) ,1 });		// transform Matrix into polynomial implicitly
		desired_det *= Matrix<ty>(1, 2, { Complex(4) ,1 });		// this will be undecomposible
		desired_det *= Matrix<ty>(1, 2, { Complex(5) ,1 });
		desired_det *= Matrix<ty>(1, 2, { Complex(6) ,1 });

		//cout << "desired_det" << desired_det;			// this matches with det_PL

		// find the roots
		Matrix<ty> R(n, n, '0');
		int n_minus_1 = n - 1;
		for (int i = 0; i < n_minus_1; ++i) {
			R(i + 1, i) = 1;
		}
		for (int i = 0; i < n; ++i) {
			R(0, i) = -desired_det(n - 1 - i) / desired_det(n);
		}
		cout << "R" << R;

		Matrix<ty> rval;
		Matrix<ty> rvec;
		//R.ritz_inner(rval, false, false, 500, 1e-8);		// eig is almost the same
		//cout << "rval" << rval;

		R.ritz_Hessenberg(rval, rvec, false);
		cout << "rval" << rval;
		cout << "rvec" << rvec;
		// the Matrix R may be undiagonalizable, hence the inverse of revc may be invalid
		cout << "rvec * rval.place_diag() * rvec.inv()" << rvec * rval.place_diag() * rvec.inv();
		cout << "R - rvec * rval.place_diag() * rvec.inv()" << R - rvec * rval.place_diag() * rvec.inv();
	}
	static void edobtm_test() {
		typedef Complex ty;
		const int n = 6;
		const int k = 2;
		Matrix<ty> A(n, n, {
			3,	7,	3,	8,	0,	0,
			0,	4,	5,	0,	0,	0,
			1,	4,	2,	5,	-5,	-7,
			9,	-1,	4,	4,	2,	7,
			0,	0,	-2,	6,	2,	10,
			0,	0,	4,	1,	-3,	6
			});
		Matrix<ty> eigval;
		Matrix<ty> eigvec;
		Eigen_decomposition_of_block_tridiagonal_Matrices::solve(A, k, eigval, eigvec);
		cout << "eigval" << eigval;
		cout << "eigvec" << eigvec;

		// not very accuracy...
		cout << "A - eigvec * eigval * eigvec.inv()" << A - eigvec * eigval * eigvec.inv();		// finishing... 
	}

	static void Hermitian_test() {
		typedef Complex ty;
		Matrix<ty> A(2, 2, {
			{1.0, 1.0}, { 0.4,-5.0 },
			{ 4.0,-0.6 }, { 1.2,0.9 },
			});
		cout << "A" << A;

		cout << "A.Hermit()" << A.Hermit();
	}
	static void givens_Hessenberg_test() {
		typedef Complex ty;

		int r = 5;
		//Matrix<ty> A(r, r, my::rand_u);
		//for (int j = 0; j < r; ++j) {
		//	for (int i = j + 2; i < r; ++i) {
		//		A(i, j) = 0;		// make it be Hassenberg
		//	}
		//}

		Matrix<ty> A(5, 5, {
			{0.8147,0.9},0.0975,0.1576,0.1419,0.6557,
			0.9058,0.2785,0.9706,0.4218,0.0357,
			0,0.5469,0.9572,0.9157,0.8491,
			0, 0,0.4854,0.7922,0.9340,
			0, 0, 0,0.9595,0.6787
			});
		/*Matrix<ty> A(3, 3, {
			{1,0.5},{4,-0.7},{0.5,-3.7},
			{0.7,-1},{0.33,-5}, {2.6,1.7},
			{0,0},{2,-0.7},{3.21,-1.7}
			});*/
		cout << "A" << A;

		Matrix<ty> Q, R;
		A.QR(Q, R);
		cout << "Q1" << Q;
		cout << "R1" << R;
		cout << "Q1 * Q1.Hermit()" << Q * Q.Hermit();
		cout << "A-Q1 * R1 " << A - Q * R;

		A.givens_Hessenberg(Q, R);
		cout << "Q2" << Q;
		cout << "R2" << R;
		cout << "Q2 * Q2.Hermit()" << Q * Q.Hermit();
		cout << "A-Q2.multiply_result_Hessenberg(R2) " << A - Q.multiply_result_Hessenberg(R);
		cout << "R * Q - R.multiply_result_Hessenberg(Q)" << R * Q - R.multiply_result_Hessenberg(Q);

	}
	static void givens_Tridiagonal_test() {
		typedef Complex ty;

		int r = 5;
		/* for double */
		Matrix<ty> A(r, r, '0');
		A(0, 0) = ty(my::rand_u());
		for (int i = 1; i < r; ++i) {
			// further more force A to be hermitian, and tridiagonalize
			A(i, i) = ty(my::rand_u());
			A(i, i - 1) = ty(my::rand_u(), my::rand_u());
			A(i - 1, i) = my::conj(A(i, i - 1));
		}

		/*Matrix<ty> A(5, 5, {
			{0.8147,0.9},0.0975,0.1576,0.1419,0.6557,
			0.9058,0.2785,0.9706,0.4218,0.0357,
			0,0.5469,0.9572,0.9157,0.8491,
			0, 0,0.4854,0.7922,0.9340,
			0, 0, 0,0.9595,0.6787
			});*/

		/*Matrix<ty> A(3, 3, {
			{1,0.5},{4,-0.7},{0,0},
			{0.7,-1},{0.33,-5}, {2.6,1.7},
			{0,0},{2,-0.7},{3.21,-1.7}
			});*/
		cout << "A" << A;

		Matrix<ty> Q, R;
		A.QR(Q, R);
		cout << "Q1" << Q;
		cout << "R1" << R;
		cout << "Q1 * Q1.Hermit()" << Q * Q.Hermit();
		cout << "A-Q1 * R1 " << A - Q * R;
		cout << "R * Q" << R * Q;

		A.givens_Tridiagonal(Q, R);
		cout << "Q2" << Q;
		cout << "R2" << R;
		cout << "Q2 * Q2.Hermit()" << Q * Q.Hermit();
		cout << "A-Q2.multiply_result_Tridiagonal(R2) " << A - Q.multiply_result_Tridiagonal(R);

		cout << "R * Q" << R * Q;
		cout << "R.multiply_result_Tridiagonal(Q)" << R.multiply_result_Tridiagonal(Q);
		cout << "R * Q - R.multiply_result_Tridiagonal(Q)" << R * Q - R.multiply_result_Tridiagonal(Q);
	}
	static void ritz_inner_test2() {
		typedef Complex ty;

		int r = 10;											// fail at 15
		Matrix<ty> A(r, r, '0');
		for (int j = 0; j < r; ++j) {
			for (int i = 0; i < j + 2; ++i) {				// make it be Hassenberg
				A(i, j) = ty(my::rand_u(), my::rand_u());
			}
		}

		/*Matrix<ty> A(5, 5, {
			{0.8147,+0.9},0.0975,0.1576,0.1419,0.6557,
			0.9058,0.2785,0.9706,0.4218,0.0357,
			0,0.5469,0.9572,0.9157,0.8491,
			0, 0,0.4854,0.7922,0.9340,
			0, 0, 0,0.9595,0.6787
			});*/

		/*Matrix<ty> A(3, 3, {
			{1,0.5},{4,-0.7},{0.5,-3.7},
			{0.7,-1},{0.33,-5}, {2.6,1.7},
			{0,0},{2,-0.7},{3.21,-1.7}
			});*/

		/*Matrix<ty> A(4, 4, {
			8, 24, 32, 16,
			0, 4, 6, 1,
			0, 0, 3, 0,
			0, 0, 0, 3
			});*/
		/*Matrix<ty> A(4, 4, {
			8, 0, 0, 0,
			0, 4, 0, 1,
			0, 0, 3, 0,
			0, 0, 0, 3
			});*/
		cout << "A" << A;
		Matrix<ty> A_ritz_val, A_ritz_vec;
		A.ritz_Hessenberg(A_ritz_val, A_ritz_vec, false);
		cout << "A_ritz_val" << A_ritz_val;
		cout << "A_ritz_vec" << A_ritz_vec;
		cout << "A - A_ritz_vec * A_ritz_val.place_diag() * A_ritz_vec.inv()" << A - A_ritz_vec * A_ritz_val.place_diag() * A_ritz_vec.inv();
	}
	static void set_cols_test() {
		Matrix<int> A(5, 5, 'i');
		cout << "A" << A;

		Matrix<int> B(3, 3, 'N');
		B(0, 0) = 9;
		cout << "B" << B;

		Matrix<int> col_ind(1, 2, { 3,4 });
		A.set_cols(col_ind, B, true);
		cout << "A" << A;
	}
	static void ritz_outer_test() {
		typedef Complex ty;

		int r = 6;
		Matrix<ty> A(r, r, '0');
		for (int j = 0; j < r; ++j) {
			for (int i = 0; i < r; ++i) {					// make it be Hassenberg
				A(i, j) = ty(my::rand_u(), my::rand_u());
			}
		}
		A = A * A.Hermit();
		cout << "A" << A;
		Matrix<ty> A_ritz_val, A_ritz_vec;
		A.ritz_non_defficient(A_ritz_val, A_ritz_vec, true);		// using the default precision
		cout << "A_ritz_val" << A_ritz_val;
		cout << "A_ritz_vec" << A_ritz_vec;
		cout << "A - A_ritz_vec * A_ritz_val.place_diag() * A_ritz_vec.inv()" << A - A_ritz_vec * A_ritz_val.place_diag() * A_ritz_vec.inv();
	}
	static void ritz_test2() {
		typedef Complex ty;

		int r = 4;
		Matrix<ty> A(r + 1, r + 1, '0');						// this will any Matrix
		for (int j = 0; j < r; ++j) {
			for (int i = 0; i < r; ++i) {			// make it be Hassenberg
				A(i, j) = ty(my::rand_u()/*, my::rand_u()*/);
				//A(i, j) = 0;
			}
		}
		A = A * A.Hermit();
		cout << "A" << A;
		Matrix<ty> A_ritz_val, A_ritz_vec;
		A.ritz(A_ritz_val, A_ritz_vec, true, -1);		// using the default precision
		cout << "A_ritz_val" << A_ritz_val;
		cout << "A_ritz_vec" << A_ritz_vec;
		// if Matrix is not squared, inv() return pseudo inverse
		cout << "A - A_ritz_vec * A_ritz_val.place_diag() * A_ritz_vec.inv()" << A - A_ritz_vec * A_ritz_val.place_diag() * A_ritz_vec.inv();
	}
	static void svd_val_test() {
		typedef Complex ty;

		int r = 7;
		int c = 4;
		Matrix<ty> A(r, c);			// this will any Matrix
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {			// make it be Hassenberg
				A(i, j) = ty(my::rand_u(), my::rand_u());
				//A(i, j) = 0;
			}
		}
		cout << "A" << A;
		// interface to be designed.
		Matrix<ty> Sigma;
		A.svd(Sigma);
		cout << "Sigma" << Sigma;
	}
	static void svd_test() {		// all okay, done
		typedef Complex ty;

		int r = 7;
		int c = 4;
		Matrix<ty> A(r, c);			// this will any Matrix
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {		// make it be Hassenberg
				A(i, j) = ty(my::rand_u(), my::rand_u());
				//A(i, j) = 0;
			}
		}
		cout << "A" << A;
		// interface to be designed.
		Matrix<ty> U, Sigma, VH;
		A.svd(U, Sigma, VH, -1);
		cout << "U" << U;
		cout << "Sigma" << Sigma;
		cout << "VH" << VH;

		/* validation */
		cout << "U.Hermit() * U" << U.Hermit() * U;
		cout << "VH * VH.Hermit()" << VH * VH.Hermit();
		cout << "A - U * Sigma * VH" << A - U * Sigma * VH;
	}
};
