/*
***********************************************************************
This file is program of MSSA.
***********************************************************************

Oct 13, 2023
Copyright 2023

Zuwei Huang
hzw1498218560@tongji.edu.cn
School of Ocean and Earth Science, Tongji University
Integrated Geophysics Group

version 1.1.0

* ***********************************************************************
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/Core>


using namespace std;
using namespace Eigen;

void read_grd_data(string fin_1, MatrixXd& T, int& nx, int& ny, double& xmin, double& xmax, double& ymin, double& ymax);
void write_grd_data(string fou_1, MatrixXd& T, int& nx, int& ny, double& xmin, double& xmax, double& ymin, double& ymax);
void construct_Hankel(MatrixXd& H, MatrixXd D, int Lr, int Kr, int i);
double diagave(MatrixXd A);
void reverse_diag_aver(MatrixXd& H, VectorXd& d);
MatrixXd diagave_block(MatrixXd A, int row, int col);
void reverse_diag_aver_block(MatrixXd& W, MatrixXd& H, MatrixXd& D);

int main()
{
	std::cout << std::endl;
	std::cout << "***************************************************************************" << std::endl << std::endl;
	std::cout << "                   SSA2D                                  " << std::endl << std::endl;
	std::cout << "  Copyright 2023, Zuwei Huang         " << std::endl << std::endl;

	std::cout << "  SSA2D is free software: you can redistribute it and/or modify it under the " << std::endl <<
		"  terms of the GNU Lesser General Public License as published by the Free " << std::endl <<
		"  Software Foundation, version 3 of the License. " << std::endl << std::endl;

	std::cout << "  SSA2D is distributed in the hope that it will be useful, but WITHOUT ANY " << std::endl <<
		"  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS " << std::endl <<
		"  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for " << std::endl <<
		"  more details. " << std::endl << std::endl;
	std::cout << "***************************************************************************" << std::endl << std::endl;

	int Nr, Nc, Lr, Lc, Kr, Kc, Nrm, Ncm, *k, reconstruct,num_k;
	double xmin, xmax, ymin, ymax;
	string tmp, input_filename, output_regional, output_eigenvalue, output_local;
	MatrixXd D, W, H, W1, H_tmp,D1,D2,S1;

	//read parameter.dat
	ifstream fpar("par.dat");
	//grd filename
	fpar >> tmp >> input_filename;
	fpar >> tmp >> output_regional;
	fpar >> tmp >> output_local;
	fpar >> tmp >> output_eigenvalue;
	//MSSA parameter
	fpar >> tmp >> Lr;
	fpar >> tmp >> Lc;
	fpar >> tmp >> num_k;
	k = new int[num_k];
	fpar >> tmp;
	for (int i = 0; i < num_k; i++)fpar >> k[i];
	fpar >> tmp >> reconstruct;
	fpar.close();
	//read grd files
	read_grd_data(input_filename, D, Nr, Nc, xmin, xmax, ymin, ymax);
	Kr = Nr - Lr + 1;
	Kc = Nc - Lc + 1;
	Ncm = Kr * Kc;
	Nrm = Lr * Lc;

	//constructing block Hankel matrix//
	cerr << "Inputfile Name: " << input_filename << endl;
	cerr << "Window Size: " << Lr<<" "<<Lc << endl;
	cerr << "Data Matrix[Nr(nx) Nc(ny)]: " << Nr << " " << Nc << endl;
	cerr << "Constructed Block Hankel Matrix[Kr*Kc Nr*Nc]: " << Nrm << " " << Ncm << endl;
	

	W.resize(Nrm, Ncm);
	H.resize(Lr, Kr);
	D1.resize(Nr, Nc);
	D2.resize(Nr, Nc);

	for (int i = 0; i < Lc; i++)
	{
		for (int j = 0; j < Kc; j++)
		{
			construct_Hankel(H, D, Lr, Kr, j + i);
			W.block(i * H.rows(), j * H.cols(), H.rows(), H.cols()) = H;
		}
	}

	//SVD decompsition
	cerr << endl;
	cerr << "--------------------------------------" << endl;
	cerr << "SVD Decompsition Start!" << endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	JacobiSVD<MatrixXd> svd(W, Eigen::ComputeThinU | Eigen::ComputeThinV);
	MatrixXd U = svd.matrixU();
	MatrixXd S = svd.singularValues().asDiagonal();
	MatrixXd V = svd.matrixV();
	cerr << "SVD Decompsition Complete!" << endl;
	cerr << "--------------------------------------" << endl;
	std::cerr << endl;
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	std::cerr << "SVD Decompsition Cost Time: " << duration.count() << " s" << endl;
	cerr << endl;
	ofstream ofeigen(output_eigenvalue);
	for (int ii = 0; ii < S.cols(); ii++)
	{
		ofeigen << ii << " " << S(ii, ii) << endl;
	}
	//SVD decopmsition down!

	if (reconstruct == 1)
	{
		//reconstruct the signal
		for (int i = 0; i < num_k; i++)
		{
			cerr << "--------------------------------------" << endl;
			cerr << "Reconstructing Order: "<<k[i] << endl;
			cerr << "signal reconstructing start!" << endl;
			S1 = S;
			S1.block(k[i], k[i], S.rows() - k[i], S.cols() - k[i]).setZero();
			W1 = U * S1 * V.transpose();
			//reconstruct(colum)
			reverse_diag_aver_block(W1, H, D1);
			D2 = D - D1;
			string output = output_regional + to_string(k[i]) + ".grd";
			string output2 = output_local + to_string(k[i]) + ".grd";
			write_grd_data(output, D1, Nr, Nc, xmin, xmax, ymin, ymax);
			write_grd_data(output2, D2, Nr, Nc, xmin, xmax, ymin, ymax);
			cerr << "signal reconstructing complete!" << endl;
			cerr << "--------------------------------------" << endl;
			cerr << endl;
		}
	}
	//reconstruct the signal down!

}


void read_grd_data(string fin_1, MatrixXd& T, int& nx, int& ny, double& xmin, double& xmax, double& ymin, double& ymax)
{
	ifstream ifin_1(fin_1);
	if (!ifin_1)
	{
		cout << "can't open the file: " << fin_1 << endl;
		exit(0);
	}
	char tmpchar[256];
	int nx_tmp, ny_tmp, ix, iy;
	double tmp1, tmp2;
	ifin_1 >> tmpchar;
	ifin_1 >> nx >> ny;
	ifin_1 >> xmin >> xmax ;
	ifin_1 >> ymin >> ymax;
	ifin_1 >> tmp1 >> tmp2;
	T.resize(nx, ny);
	for (iy = 0; iy < ny; iy++)
		for (ix = 0; ix < nx; ix++)
			ifin_1 >> T.coeffRef(ix, iy);
	ifin_1.close();
	ifin_1.clear();
}
void write_grd_data(string fou_1, MatrixXd& T, int& nx, int& ny, double& xmin, double& xmax, double& ymin, double& ymax)
{
	ofstream ofou_1(fou_1);
	ofou_1 << "DSAA" << endl;
	ofou_1 << nx << " " << ny << endl;
	ofou_1 << xmin << " " << xmax << endl;
	ofou_1 << ymin << " " << ymax << endl;
	double fmin, fmax;
	fmin = T.minCoeff();
	fmax = T.maxCoeff();
	ofou_1 << fmin << " " << fmax << endl;
	for (int iy = 0; iy < ny; iy++)
	{
		for (int ix = 0; ix < nx; ix++)
			ofou_1 << T.coeff(ix,iy) << " ";
		ofou_1 << endl;
		ofou_1 << endl;
	}
	ofou_1.close();
	ofou_1.clear();
}
void construct_Hankel(MatrixXd& H, MatrixXd D, int Lr, int Kr, int i)
{
	for (int ii = 0; ii < Lr; ii++)
	{
		for (int jj = 0; jj < Kr; jj++)
		{
			H.coeffRef(ii, jj) = D.coeff(jj + ii, i);
		}
	}
}
double diagave(MatrixXd A)
{
	double sum = 0.0;
	for (int i = 0; i < A.rows(); i++) {
		sum += A(i, A.rows() - i - 1);
	}
	sum /= A.rows();
	return sum;
}
void reverse_diag_aver(MatrixXd& H, VectorXd& d)
{
	MatrixXd tmp;
	d.resize(H.cols() + H.rows() - 1);

	for (int i = 0; i < H.rows(); i++)
	{
		tmp = H.block(0, 0, i + 1, i + 1);
		d(i) = diagave(tmp);
	}

	for (int i = H.rows(); i < H.cols(); i++)
	{
		tmp = H.block(0, i - (H.rows()) + 1, H.rows(), H.rows());
		d(i) = diagave(tmp);
	}

	for (int i = 1; i < H.rows(); i++)
	{
		tmp = H.block(i, i + H.cols() - H.rows(), H.rows() - i, H.rows() - i);
		d(H.cols() - 1 + i) = diagave(tmp);
	}

}
MatrixXd diagave_block(MatrixXd A, int row, int col)
{
	MatrixXd sum;
	sum.resize(row, col);
	sum.setZero();

	for (int i = 0; i < A.cols() / col; i++) {
		sum += A.block(i * row, (A.cols() / col - i - 1) * col, row, col);
		/*cerr << endl;
		cerr<< A.block(i * row, (A.cols() / col - i - 1) * col, row, col)<<endl;*/
	}
	sum /= A.rows() / row;
	/*cerr << endl;
	cerr << sum<< endl;
	cerr << endl;*/
	return sum;

}
void reverse_diag_aver_block(MatrixXd& W, MatrixXd& H, MatrixXd& D)
{
	int Kc = W.cols() / H.cols();
	int Lc = W.rows() / H.rows();

	MatrixXd  tmp;
	VectorXd dd;
	for (int i = 0; i < Lc; i++)
	{
		tmp = W.block(0, 0, (i + 1) * H.rows(), (i + 1) * H.cols());
		H = diagave_block(tmp, H.rows(), H.cols());
		reverse_diag_aver(H, dd);
		D.col(i) = dd;
	}

	for (int i = Lc; i < Kc; i++)
	{
		tmp = W.block(0, (i - Lc + 1) * H.cols(), Lc * H.rows(), Lc * H.cols());
		H = diagave_block(tmp, H.rows(), H.cols());
		reverse_diag_aver(H, dd);
		D.col(i) = dd;
	}

	for (int i = 1; i < Lc; i++)
	{
		tmp = W.block(i * H.rows(), (i + Kc - Lc) * H.cols(), (Lc - i) * H.rows(), (Lc - i) * H.cols());
		H = diagave_block(tmp, H.rows(), H.cols());
		reverse_diag_aver(H, dd);
		D.col(Kc - 1 + i) = dd;
	}

}

