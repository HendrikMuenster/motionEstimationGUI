/*
% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.1
% Date: 2015-09-26

% All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author and University of Muenster not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
*/
#include "mex.h"
#include "math.h"
#include <omp.h>
#include "tools.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//L2L2OpticalFlow(u1,u2,tol,lambda,ux,uy,ut)
{
	double *u1 = mxGetPr(prhs[0]);
	double *u2 = mxGetPr(prhs[1]);

	float tol = (float)mxGetScalar(prhs[2]);
	float lambda = (float)mxGetScalar(prhs[3]);

	int maxIterations = (int)mxGetScalar(prhs[4]);

	const size_t *sizeImage = mxGetDimensions(prhs[0]);

	double *inputV = 0;
	double *inputY = 0;

	int typeNorm = 2;

	//if (nrhs > 5)
	//{
	//	typeNorm = (int)mxGetScalar(prhs[5]);
	//}
	
	if (nrhs > 6)
	{
		inputV = mxGetPr(prhs[6]);
	}

	if (nrhs > 7)
	{
		inputY = mxGetPr(prhs[7]);
	}

	int nPx = (int)(sizeImage[0] * sizeImage[1]);

	const size_t sizeY[2] = { 4 * nPx, 1 };

	// Output v1
	plhs[0] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
	double *Outv1 = mxGetPr(plhs[0]);

	// Output v2
	plhs[1] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
	double *Outv2 = mxGetPr(plhs[1]);

	// Output  Y
	plhs[2] = mxCreateNumericArray(2, sizeY, mxDOUBLE_CLASS, mxREAL);
	double *YOut = mxGetPr(plhs[2]);

	float* v1 = new float[nPx];
	float* v2 = new float[nPx];

	float* ux = new float[nPx];
	float* uy = new float[nPx];
	float* ut = new float[nPx];

	float* uxut = new float[nPx];
	float* uyut = new float[nPx];

	float* c1 = new float[nPx];
	float* c2 = new float[nPx];
	float* c3 = new float[nPx];
	float* teiler = new float[nPx];

	float* y11 = new float[nPx];
	float* y12 = new float[nPx];
	float* y21 = new float[nPx];
	float* y22 = new float[nPx];

	float* Kty1 = new float[nPx];
	float* Kty2 = new float[nPx];

	float* Kx11 = new float[nPx];
	float* Kx12 = new float[nPx];
	float* Kx21 = new float[nPx];
	float* Kx22 = new float[nPx];

	float tau = 0.25;
	float sigma = 0.5;
	float dTau = 1.0f / tau;
	float dSigma = 1.0f / sigma;

	//residuals
	float p = 0;
	float d = 0;
	float err = 1.0;

	float ssl = 1.0f - sigma / (sigma + lambda);

	#pragma omp parallel for
	for (int j = 0; j < sizeImage[1]; ++j)
	{
		for (int i = 0; i < sizeImage[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeImage, i, j);

			//Index for gradients
			ut[tmpIndex] = (float)(u2[tmpIndex] - u1[tmpIndex]);

			if (i>0 && i < sizeImage[0] - 1)
			{
				uy[tmpIndex] = (float)(0.5f * (u1[index2DtoLinear(sizeImage, i + 1, j)] - u1[index2DtoLinear(sizeImage, i - 1, j)]));
			}
			else
			{
				uy[tmpIndex] = 0.0;
			}

			if (j>0 && j < sizeImage[1] - 1)
			{
				ux[tmpIndex] = (float)(0.5 * (u1[index2DtoLinear(sizeImage, i, j + 1)] - u1[index2DtoLinear(sizeImage, i, j - 1)]));
			}
			else
			{
				ux[tmpIndex] = 0.0f;
			}

			uxut[tmpIndex] = ux[tmpIndex] * ut[tmpIndex];
			uyut[tmpIndex] = uy[tmpIndex] * ut[tmpIndex];

			c1[tmpIndex] = 1.0f + tau * ux[tmpIndex] * ux[tmpIndex];
			c2[tmpIndex] = tau * ux[tmpIndex] * uy[tmpIndex];
			c3[tmpIndex] = 1.0f + tau * uy[tmpIndex] * uy[tmpIndex];

			teiler[tmpIndex] = 1.0f / (c1[tmpIndex] * c3[tmpIndex] - c2[tmpIndex] * c2[tmpIndex]);


			if (nrhs > 6)
			{
				v1[tmpIndex] = (float)inputV[tmpIndex];
				v2[tmpIndex] = (float)inputV[nPx + tmpIndex];
			}
			else
			{
				v1[tmpIndex] = 0;
				v2[tmpIndex] = 0;
			}

			Kty1[tmpIndex] = 0;
			Kty2[tmpIndex] = 0;

			if (nrhs > 7)
			{
				y11[tmpIndex] = (float)inputY[tmpIndex];
				y12[tmpIndex] = (float)inputY[nPx + tmpIndex];
				y21[tmpIndex] = (float)inputY[2 * nPx + tmpIndex];
				y22[tmpIndex] = (float)inputY[3 * nPx + tmpIndex];
			}
			else
			{
				y11[tmpIndex] = 0;
				y12[tmpIndex] = 0;
				y21[tmpIndex] = 0;
				y22[tmpIndex] = 0;
			}

			Kx11[tmpIndex] = 0;
			Kx12[tmpIndex] = 0;
			Kx21[tmpIndex] = 0;
			Kx22[tmpIndex] = 0;
		}
	}

	int iterations = 0;

	while (err > tol && iterations <= maxIterations)
	{
		++iterations;

		if (iterations % 50 == 0)
		{
			p = 0;
			d = 0;
		}

		//primal step
		#pragma omp parallel for reduction(+:p)
		for (int j = 0; j < sizeImage[1]; ++j)
		{
			for (int i = 0; i < sizeImage[0]; ++i)
			{
				int tmpIndex = index2DtoLinear(sizeImage, i, j);

				float Kty1Old = Kty1[tmpIndex];
				float Kty2Old = Kty2[tmpIndex];

				//transpose equals -div
				Kty1[tmpIndex] = -(dxm(y11, sizeImage, i, j) + dym(y12, sizeImage, i, j));
				Kty2[tmpIndex] = -(dxm(y21, sizeImage, i, j) + dym(y22, sizeImage, i, j));

				float b1 = v1[tmpIndex] - tau*(Kty1[tmpIndex] + uxut[tmpIndex]);
				float b2 = v2[tmpIndex] - tau*(Kty2[tmpIndex] + uyut[tmpIndex]);

				float v1Old = v1[tmpIndex];
				float v2Old = v2[tmpIndex];

				v1[tmpIndex] = (b1 * c3[tmpIndex] - c2[tmpIndex] * b2) * teiler[tmpIndex];
				v2[tmpIndex] = (b2 * c1[tmpIndex] - c2[tmpIndex] * b1) * teiler[tmpIndex];

				if (iterations % 50 == 0)
				{
					//residuals
					p += myAbs((v1Old - v1[tmpIndex]) * dTau - Kty1Old + Kty1[tmpIndex])
						+ myAbs((v2Old - v2[tmpIndex]) * dTau - Kty2Old + Kty2[tmpIndex]);
				}
			}
		}

		//dual step
		#pragma omp parallel for reduction(+:d)
		for (int j = 0; j < sizeImage[1]; ++j)
		{
			for (int i = 0; i < sizeImage[0]; ++i)
			{
				int tmpIndex = index2DtoLinear(sizeImage, i, j);

				float Kx11Old = Kx11[tmpIndex];
				float Kx12Old = Kx12[tmpIndex];
				float Kx21Old = Kx21[tmpIndex];
				float Kx22Old = Kx22[tmpIndex];

				Kx11[tmpIndex] = dxp(v1, sizeImage, i, j);
				Kx12[tmpIndex] = dyp(v1, sizeImage, i, j);
				Kx21[tmpIndex] = dxp(v2, sizeImage, i, j);
				Kx22[tmpIndex] = dyp(v2, sizeImage, i, j);

				float y11Old = y11[tmpIndex];
				float y12Old = y12[tmpIndex];
				float y21Old = y21[tmpIndex];
				float y22Old = y22[tmpIndex];

				float y11Tilde = y11[tmpIndex] + sigma*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old);
				float y12Tilde = y12[tmpIndex] + sigma*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old);
				float y21Tilde = y21[tmpIndex] + sigma*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old);
				float y22Tilde = y22[tmpIndex] + sigma*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old);

				if (typeNorm == 1)
				{
					y11[tmpIndex] = myMax(-lambda, myMin(lambda, y11Tilde));
					y12[tmpIndex] = myMax(-lambda, myMin(lambda, y12Tilde));
					y21[tmpIndex] = myMax(-lambda, myMin(lambda, y21Tilde));
					y22[tmpIndex] = myMax(-lambda, myMin(lambda, y22Tilde));
				}
				else if (typeNorm == 2)
				{
					float divisor1 = myMax(1.0f, sqrtf(y11Tilde*y11Tilde + y12Tilde*y12Tilde) / lambda);
					float divisor2 = myMax(1.0f, sqrtf(y21Tilde*y21Tilde + y22Tilde*y22Tilde) / lambda);

					y11[tmpIndex] = y11Tilde / divisor1;
					y12[tmpIndex] = y12Tilde / divisor1;
					y21[tmpIndex] = y21Tilde / divisor2;
					y22[tmpIndex] = y22Tilde / divisor2;
				}
				else if (typeNorm == 3)
				{
					float norm1 = sqrtf(y11Tilde*y11Tilde + y12Tilde*y12Tilde + y21Tilde*y21Tilde + y22Tilde*y22Tilde) / lambda;

					if (norm1 > 1.0f)
					{
						y11[tmpIndex] = y11Tilde / norm1;
						y12[tmpIndex] = y12Tilde / norm1;
						y21[tmpIndex] = y21Tilde / norm1;
						y22[tmpIndex] = y22Tilde / norm1;
					}
					else
					{
						y11[tmpIndex] = y11Tilde;
						y12[tmpIndex] = y12Tilde;
						y21[tmpIndex] = y21Tilde;
						y22[tmpIndex] = y22Tilde;
					}
				}

				if (iterations % 50 == 0)
				{
					d += myAbs((y11Old - y11[tmpIndex]) * dSigma - Kx11Old + Kx11[tmpIndex]) +
						myAbs((y12Old - y12[tmpIndex]) * dSigma - Kx12Old + Kx12[tmpIndex]) +
						myAbs((y21Old - y21[tmpIndex]) * dSigma - Kx21Old + Kx21[tmpIndex]) +
						myAbs((y22Old - y22[tmpIndex]) * dSigma - Kx22Old + Kx22[tmpIndex]);
				}
			}
		}

		if (iterations % 50 == 0)
		{
			err = (d*d + p*p) / nPx;
		}

		if (iterations % 1000 == 0)
		{
			mexPrintf("Iteration %d,Residual %e\n", iterations, err);
			mexEvalString("drawnow;");
		}
	}

	//write output
	#pragma omp parallel for
	for (int j = 0; j < sizeImage[1]; ++j)
	{
		for (int i = 0; i < sizeImage[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeImage, i, j);

			YOut[tmpIndex] = (double)y11[tmpIndex];
			YOut[tmpIndex + nPx] = (double)y12[tmpIndex];
			YOut[tmpIndex + 2 * nPx] = (double)y21[tmpIndex];
			YOut[tmpIndex + 3 * nPx] = (double)y22[tmpIndex];

			Outv1[tmpIndex] = (double)v1[tmpIndex];
			Outv2[tmpIndex] = (double)v2[tmpIndex];
		}
	}

	delete[] v1;
	delete[] v2;

	delete[] ux;
	delete[] uy;
	delete[] ut;

	delete[] uxut;
	delete[] uyut;

	delete[] c1;
	delete[] c2;
	delete[] c3;
	delete[] teiler;

	delete[] y11;
	delete[] y12;
	delete[] y21;
	delete[] y22;

	delete[] Kty1;
	delete[] Kty2;

	delete[] Kx11;
	delete[] Kx12;
	delete[] Kx21;
	delete[] Kx22;
}

