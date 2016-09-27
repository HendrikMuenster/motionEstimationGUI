/*
% Author Information:
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.0
% Date: 2015-10-03
%
% History:

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
//L1TVTVOpticalFlowCPP(u1,u2,tol,lambda0,lambda1,maxIterations)
{
	double *u1 = mxGetPr(prhs[0]);
	double *u2 = mxGetPr(prhs[1]);

	float tol = (float)mxGetScalar(prhs[2]);
	float lambda0 = (float)mxGetScalar(prhs[3]);
	float lambda1 = (float)mxGetScalar(prhs[4]);

	int maxIterations = (int)mxGetScalar(prhs[5]);

	const size_t *sizeImage = mxGetDimensions(prhs[0]);

	int typeNorm = 3;
	if (nrhs > 6)
	{
		typeNorm = (int)mxGetScalar(prhs[6]);
	}

	double *inputV = 0;
	double *inputY = 0;

	if (nrhs > 7)
	{
		inputV = mxGetPr(prhs[7]);
	}

	if (nrhs > 8)
	{
		inputY = mxGetPr(prhs[8]);
	}

	int nPx = (int)(sizeImage[0] * sizeImage[1]);

	const size_t sizeV[2] = { 6 * nPx, 1 };
	const size_t sizeY[2] = { 12 * nPx, 1 };

	// Output v
	plhs[0] = mxCreateNumericArray(2, sizeV, mxDOUBLE_CLASS, mxREAL);
	double *VOut = mxGetPr(plhs[0]);

	// Output  Y
	plhs[1] = mxCreateNumericArray(2, sizeY, mxDOUBLE_CLASS, mxREAL);
	double *YOut = mxGetPr(plhs[1]);

	float* v1 = new float[nPx];
	float* v2 = new float[nPx];
	float* v3 = new float[nPx];
	float* v4 = new float[nPx];
	float* v5 = new float[nPx];
	float* v6 = new float[nPx];

	float* tauAquad = new float[nPx];
	float* teiler1 = new float[nPx];
	float* teiler2 = new float[nPx];

	float* ux = new float[nPx];
	float* uy = new float[nPx];
	float* ut = new float[nPx];

	float* tauUx = new float[nPx];
	float* tauUy = new float[nPx];

	float* y1 = new float[nPx];
	float* y2 = new float[nPx];
	float* y3 = new float[nPx];
	float* y4 = new float[nPx];
	float* y5 = new float[nPx];
	float* y6 = new float[nPx];
	float* y7 = new float[nPx];
	float* y8 = new float[nPx];
	float* y9 = new float[nPx];
	float* y10 = new float[nPx];
	float* y11 = new float[nPx];
	float* y12 = new float[nPx];

	float* Kty1 = new float[nPx];
	float* Kty2 = new float[nPx];
	float* Kty3 = new float[nPx];
	float* Kty4 = new float[nPx];
	float* Kty5 = new float[nPx];
	float* Kty6 = new float[nPx];

	float* Kx1 = new float[nPx];
	float* Kx2 = new float[nPx];
	float* Kx3 = new float[nPx];
	float* Kx4 = new float[nPx];
	float* Kx5 = new float[nPx];
	float* Kx6 = new float[nPx];
	float* Kx7 = new float[nPx];
	float* Kx8 = new float[nPx];
	float* Kx9 = new float[nPx];
	float* Kx10 = new float[nPx];
	float* Kx11 = new float[nPx];
	float* Kx12 = new float[nPx];


	float tau1 = 0.25f;
	float tau2 = 0.2f;

	float sigma1 = 0.33f;
	float sigma2 = 0.5f;


	//residuals
	float p = 0;
	float d = 0;
	float err = 1.0;

	#pragma omp parallel for
	for (int j = 0; j < sizeImage[1]; ++j)
	{
		for (int i = 0; i < sizeImage[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeImage, i, j);

			ut[tmpIndex] = (float)(u2[tmpIndex] - u1[tmpIndex]);

			if (i>0 && i < sizeImage[0] - 1)
			{
				uy[tmpIndex] = (float)(0.5f * (u1[index2DtoLinear(sizeImage, i + 1, j)] - u1[index2DtoLinear(sizeImage, i - 1, j)]));
			}
			else
			{
				uy[tmpIndex] = 0.0f;
			}

			if (j>0 && j < sizeImage[1] - 1)
			{
				ux[tmpIndex] = (float)(0.5 * (u1[index2DtoLinear(sizeImage, i, j + 1)] - u1[index2DtoLinear(sizeImage, i, j - 1)]));
			}
			else
			{
				ux[tmpIndex] = 0.0f;
			}

			tauAquad[tmpIndex] = (ux[tmpIndex] * ux[tmpIndex] + uy[tmpIndex] * uy[tmpIndex]);


			teiler1[tmpIndex] = ux[tmpIndex] / (tauAquad[tmpIndex]);
			teiler2[tmpIndex] = uy[tmpIndex] / (tauAquad[tmpIndex]);

			tauAquad[tmpIndex] = tau1*tauAquad[tmpIndex];

			tauUx[tmpIndex] = ux[tmpIndex] * tau1;
			tauUy[tmpIndex] = uy[tmpIndex] * tau1;

			if (nrhs > 7)
			{
				v1[tmpIndex] = (float)inputV[0 * nPx + tmpIndex];
				v2[tmpIndex] = (float)inputV[1 * nPx + tmpIndex];
				v3[tmpIndex] = (float)inputV[2 * nPx + tmpIndex];
				v4[tmpIndex] = (float)inputV[3 * nPx + tmpIndex];
				v5[tmpIndex] = (float)inputV[4 * nPx + tmpIndex];
				v6[tmpIndex] = (float)inputV[5 * nPx + tmpIndex];
			}
			else
			{
				v1[tmpIndex] = 0.0;
				v2[tmpIndex] = 0.0;
				v3[tmpIndex] = 0.0;
				v4[tmpIndex] = 0.0;
				v5[tmpIndex] = 0.0;
				v6[tmpIndex] = 0.0;
			}

			Kty1[tmpIndex] = 0.0f;
			Kty2[tmpIndex] = 0.0f;
			Kty3[tmpIndex] = 0.0f;
			Kty4[tmpIndex] = 0.0f;
			Kty5[tmpIndex] = 0.0f;
			Kty6[tmpIndex] = 0.0f;

			if (nrhs > 8)
			{
				y1[tmpIndex] = (float)inputY[0 * nPx + tmpIndex];
				y2[tmpIndex] = (float)inputY[1 * nPx + tmpIndex];
				y3[tmpIndex] = (float)inputY[2 * nPx + tmpIndex];
				y4[tmpIndex] = (float)inputY[3 * nPx + tmpIndex];
				y5[tmpIndex] = (float)inputY[4 * nPx + tmpIndex];
				y6[tmpIndex] = (float)inputY[5 * nPx + tmpIndex];
				y7[tmpIndex] = (float)inputY[6 * nPx + tmpIndex];
				y8[tmpIndex] = (float)inputY[7 * nPx + tmpIndex];
				y9[tmpIndex] = (float)inputY[8 * nPx + tmpIndex];
				y10[tmpIndex] = (float)inputY[9 * nPx + tmpIndex];
				y11[tmpIndex] = (float)inputY[10 * nPx + tmpIndex];
				y12[tmpIndex] = (float)inputY[11 * nPx + tmpIndex];
			}
			else
			{
				y1[tmpIndex] = 0.0f;
				y2[tmpIndex] = 0.0f;
				y3[tmpIndex] = 0.0f;
				y4[tmpIndex] = 0.0f;
				y5[tmpIndex] = 0.0f;
				y6[tmpIndex] = 0.0f;
				y7[tmpIndex] = 0.0f;
				y8[tmpIndex] = 0.0f;
				y9[tmpIndex] = 0.0f;
				y10[tmpIndex] = 0.0f;
				y11[tmpIndex] = 0.0f;
				y12[tmpIndex] = 0.0f;
			}

			Kx1[tmpIndex] = 0.0f;
			Kx2[tmpIndex] = 0.0f;
			Kx3[tmpIndex] = 0.0f;
			Kx4[tmpIndex] = 0.0f;
			Kx5[tmpIndex] = 0.0f;
			Kx6[tmpIndex] = 0.0f;
			Kx7[tmpIndex] = 0.0f;
			Kx8[tmpIndex] = 0.0f;
			Kx9[tmpIndex] = 0.0f;
			Kx10[tmpIndex] = 0.0f;
			Kx11[tmpIndex] = 0.0f;
			Kx12[tmpIndex] = 0.0f;
		}
	}

	int iterations = 0;

	printf("lambda0=%f,lambda1=%f,maxIt=%d,tol=%f\n", lambda0, lambda1, maxIterations, tol);

	while (err > tol && iterations <= maxIterations)
	{
		++iterations;

		if (iterations % 50 == 0)
		{
			p = 0.0f;
			d = 0.0f;
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
				float Kty3Old = Kty3[tmpIndex];
				float Kty4Old = Kty4[tmpIndex];
				float Kty5Old = Kty5[tmpIndex];
				float Kty6Old = Kty6[tmpIndex];

				//transpose equals -div
				Kty1[tmpIndex] = -dxm(y1, sizeImage, i, j) - dym(y2, sizeImage, i, j);
				Kty2[tmpIndex] = -dxm(y3, sizeImage, i, j) - dym(y4, sizeImage, i, j);
				Kty3[tmpIndex] = -y1[tmpIndex] - dxm(y5, sizeImage, i, j) - dym(y6, sizeImage, i, j);
				Kty4[tmpIndex] = -y2[tmpIndex] - dxm(y7, sizeImage, i, j) - dym(y8, sizeImage, i, j);
				Kty5[tmpIndex] = -y3[tmpIndex] - dxm(y9, sizeImage, i, j) - dym(y10, sizeImage, i, j);
				Kty6[tmpIndex] = -y4[tmpIndex] - dxm(y11, sizeImage, i, j) - dym(y12, sizeImage, i, j);

				float v1Old = v1[tmpIndex];
				float v2Old = v2[tmpIndex];
				float v3Old = v3[tmpIndex];
				float v4Old = v4[tmpIndex];
				float v5Old = v5[tmpIndex];
				float v6Old = v6[tmpIndex];

				float v1Hat = v1Old - tau1*Kty1[tmpIndex];
				float v2Hat = v2Old - tau1*Kty2[tmpIndex];

				float rho = ut[tmpIndex] + ux[tmpIndex] * v1Hat + uy[tmpIndex] * v2Hat;

				if (rho <= -tauAquad[tmpIndex])
				{
					v1[tmpIndex] = v1Hat + tauUx[tmpIndex];
					v2[tmpIndex] = v2Hat + tauUy[tmpIndex];
				}
				else if (rho >= tauAquad[tmpIndex])
				{
					v1[tmpIndex] = v1Hat - tauUx[tmpIndex];
					v2[tmpIndex] = v2Hat - tauUy[tmpIndex];
				}
				else
				{
					v1[tmpIndex] = v1Hat - teiler1[tmpIndex] * rho;
					v2[tmpIndex] = v2Hat - teiler2[tmpIndex] * rho;
				}

				v3[tmpIndex] = v3Old - tau2*Kty3[tmpIndex];
				v4[tmpIndex] = v4Old - tau2*Kty4[tmpIndex];
				v5[tmpIndex] = v5Old - tau2*Kty5[tmpIndex];
				v6[tmpIndex] = v6Old - tau2*Kty6[tmpIndex];

				if (iterations % 50 == 0)
				{
					//residuals
					p += myAbs((v1Old - v1[tmpIndex]) / tau1 - Kty1Old + Kty1[tmpIndex])
						+ myAbs((v2Old - v2[tmpIndex]) / tau1 - Kty2Old + Kty2[tmpIndex])
						+ myAbs((v3Old - v3[tmpIndex]) / tau2 - Kty3Old + Kty3[tmpIndex])
						+ myAbs((v4Old - v4[tmpIndex]) / tau2 - Kty4Old + Kty4[tmpIndex])
						+ myAbs((v5Old - v5[tmpIndex]) / tau2 - Kty5Old + Kty5[tmpIndex])
						+ myAbs((v6Old - v6[tmpIndex]) / tau2 - Kty6Old + Kty6[tmpIndex]);
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

				float Kx1Old = Kx1[tmpIndex];
				float Kx2Old = Kx2[tmpIndex];
				float Kx3Old = Kx3[tmpIndex];
				float Kx4Old = Kx4[tmpIndex];
				float Kx5Old = Kx5[tmpIndex];
				float Kx6Old = Kx6[tmpIndex];
				float Kx7Old = Kx7[tmpIndex];
				float Kx8Old = Kx8[tmpIndex];
				float Kx9Old = Kx9[tmpIndex];
				float Kx10Old = Kx10[tmpIndex];
				float Kx11Old = Kx11[tmpIndex];
				float Kx12Old = Kx12[tmpIndex];

				float y1Old = y1[tmpIndex];
				float y2Old = y2[tmpIndex];
				float y3Old = y3[tmpIndex];
				float y4Old = y4[tmpIndex];
				float y5Old = y5[tmpIndex];
				float y6Old = y6[tmpIndex];
				float y7Old = y7[tmpIndex];
				float y8Old = y8[tmpIndex];
				float y9Old = y9[tmpIndex];
				float y10Old = y10[tmpIndex];
				float y11Old = y11[tmpIndex];
				float y12Old = y12[tmpIndex];

				Kx1[tmpIndex] = dxp(v1, sizeImage, i, j) - v3[tmpIndex];
				Kx2[tmpIndex] = dyp(v1, sizeImage, i, j) - v4[tmpIndex];
				Kx3[tmpIndex] = dxp(v2, sizeImage, i, j) - v5[tmpIndex];
				Kx4[tmpIndex] = dyp(v2, sizeImage, i, j) - v6[tmpIndex];

				// we need second order derivatives
				Kx5[tmpIndex] = dxp(v3, sizeImage, i, j);
				Kx6[tmpIndex] = dyp(v3, sizeImage, i, j);
				Kx7[tmpIndex] = dxp(v4, sizeImage, i, j);
				Kx8[tmpIndex] = dyp(v4, sizeImage, i, j);
				Kx9[tmpIndex] =  dxp(v5, sizeImage, i, j);
				Kx10[tmpIndex] = dyp(v5, sizeImage, i, j);
				Kx11[tmpIndex] = dxp(v6, sizeImage, i, j);
				Kx12[tmpIndex] = dyp(v6, sizeImage, i, j);

				float y1Tilde = y1[tmpIndex] + sigma1*(2 * Kx1[tmpIndex] - Kx1Old);
				float y2Tilde = y2[tmpIndex] + sigma1*(2 * Kx2[tmpIndex] - Kx2Old);
				float y3Tilde = y3[tmpIndex] + sigma1*(2 * Kx3[tmpIndex] - Kx3Old);
				float y4Tilde = y4[tmpIndex] + sigma1*(2 * Kx4[tmpIndex] - Kx4Old);


				/*
				float y5Tilde = y5[tmpIndex] + sigma2*(2 * Kx5[tmpIndex] - Kx5Old);
				float y6Tilde = y6[tmpIndex] + sigma2*(2 * Kx6[tmpIndex] - Kx6Old);
				float y7Tilde = y7[tmpIndex] + sigma2*(2 * Kx7[tmpIndex] - Kx7Old);
				float y8Tilde = y8[tmpIndex] + sigma2*(2 * Kx8[tmpIndex] - Kx8Old);
				float y9Tilde = y9[tmpIndex] + sigma2*(2 * Kx9[tmpIndex] - Kx9Old);
				float y10Tilde = y10[tmpIndex] + sigma2*(2 * Kx10[tmpIndex] - Kx10Old);
				float y11Tilde = y11[tmpIndex] + sigma2*(2 * Kx11[tmpIndex] - Kx11Old);
				float y12Tilde = y12[tmpIndex] + sigma2*(2 * Kx12[tmpIndex] - Kx12Old);
				*/

				//float norm1 = sqrtf(y1Tilde*y1Tilde + y2Tilde*y2Tilde + y3Tilde*y3Tilde + y4Tilde*y4Tilde);

				//float divisor1 = myMax(1.0f, norm1 / lambda0);
				
				//float norm2 = sqrtf(y5Tilde*y5Tilde + y6Tilde*y6Tilde + y7Tilde*y7Tilde + y8Tilde*y8Tilde + y9Tilde*y9Tilde + y10Tilde*y10Tilde + y11Tilde*y11Tilde + y12Tilde*y12Tilde);

				float divisor1 = myMax(1.0f, sqrtf(y1Tilde*y1Tilde + y2Tilde*y2Tilde) / lambda0);
				float divisor2 = myMax(1.0f, sqrtf(y3Tilde*y3Tilde + y4Tilde*y4Tilde) / lambda0);

				y1[tmpIndex] = y1Tilde / divisor1;
				y2[tmpIndex] = y2Tilde / divisor1;
				y3[tmpIndex] = y3Tilde / divisor2;
				y4[tmpIndex] = y4Tilde / divisor2;

				y5[tmpIndex] = myMax(-lambda1, myMin(lambda1, y5[tmpIndex] + sigma2*(2 * Kx5[tmpIndex] - Kx5Old)));
				y6[tmpIndex] = myMax(-lambda1, myMin(lambda1, y6[tmpIndex] + sigma2*(2 * Kx6[tmpIndex] - Kx6Old)));
				y7[tmpIndex] = myMax(-lambda1, myMin(lambda1, y7[tmpIndex] + sigma2*(2 * Kx7[tmpIndex] - Kx7Old)));
				y8[tmpIndex] = myMax(-lambda1, myMin(lambda1, y8[tmpIndex] + sigma2*(2 * Kx8[tmpIndex] - Kx8Old)));
				y9[tmpIndex] = myMax(-lambda1, myMin(lambda1, y9[tmpIndex] + sigma2*(2 * Kx9[tmpIndex] - Kx9Old)));
				y10[tmpIndex] = myMax(-lambda1, myMin(lambda1, y10[tmpIndex] + sigma2*(2 * Kx10[tmpIndex] - Kx10Old)));
				y11[tmpIndex] = myMax(-lambda1, myMin(lambda1, y11[tmpIndex] + sigma2*(2 * Kx11[tmpIndex] - Kx11Old)));
				y12[tmpIndex] = myMax(-lambda1, myMin(lambda1, y12[tmpIndex] + sigma2*(2 * Kx12[tmpIndex] - Kx12Old)));

				/*
				y5[tmpIndex] = y5Tilde / divisor2;
				y6[tmpIndex] = y6Tilde / divisor2;
				y7[tmpIndex] = y7Tilde / divisor2;
				y8[tmpIndex] = y8Tilde / divisor2;
				y9[tmpIndex] = y9Tilde / divisor2;
				y10[tmpIndex] = y10Tilde / divisor2;
				y11[tmpIndex] = y11Tilde / divisor2;
				y12[tmpIndex] = y12Tilde / divisor2;*/
		

				if (iterations % 50 == 0)
				{
					d += myAbs((y1Old - y1[tmpIndex]) / sigma1 - Kx1Old + Kx1[tmpIndex]) +
						myAbs((y2Old - y2[tmpIndex]) / sigma1 - Kx2Old + Kx2[tmpIndex]) +
						myAbs((y3Old - y3[tmpIndex]) / sigma1 - Kx3Old + Kx3[tmpIndex]) +
						myAbs((y4Old - y4[tmpIndex]) / sigma1 - Kx4Old + Kx4[tmpIndex]) +
						myAbs((y5Old - y5[tmpIndex]) / sigma2 - Kx5Old + Kx5[tmpIndex]) +
						myAbs((y6Old - y6[tmpIndex]) / sigma2 - Kx6Old + Kx6[tmpIndex]) +
						myAbs((y7Old - y7[tmpIndex]) / sigma2 - Kx7Old + Kx7[tmpIndex]) +
						myAbs((y8Old - y8[tmpIndex]) / sigma2 - Kx8Old + Kx8[tmpIndex]) +
						myAbs((y9Old - y9[tmpIndex]) / sigma2 - Kx9Old + Kx9[tmpIndex]) +
						myAbs((y10Old - y10[tmpIndex]) / sigma2 - Kx10Old + Kx10[tmpIndex]) +
						myAbs((y11Old - y11[tmpIndex]) / sigma2 - Kx11Old + Kx11[tmpIndex]) +
						myAbs((y12Old - y12[tmpIndex]) / sigma2 - Kx12Old + Kx12[tmpIndex]);
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

			YOut[tmpIndex] = (double)y1[tmpIndex];
			YOut[tmpIndex + nPx] = (double)y2[tmpIndex];
			YOut[tmpIndex + 2 * nPx] = (double)y3[tmpIndex];
			YOut[tmpIndex + 3 * nPx] = (double)y4[tmpIndex];
			YOut[tmpIndex + 4 * nPx] = (double)y5[tmpIndex];
			YOut[tmpIndex + 5 * nPx] = (double)y6[tmpIndex];
			YOut[tmpIndex + 6 * nPx] = (double)y7[tmpIndex];
			YOut[tmpIndex + 7 * nPx] = (double)y8[tmpIndex];
			YOut[tmpIndex + 8 * nPx] = (double)y9[tmpIndex];
			YOut[tmpIndex + 9 * nPx] = (double)y10[tmpIndex];
			YOut[tmpIndex + 10 * nPx] = (double)y11[tmpIndex];
			YOut[tmpIndex + 11 * nPx] = (double)y12[tmpIndex];

			VOut[tmpIndex + 0 * nPx] = (double)v1[tmpIndex];
			VOut[tmpIndex + 1 * nPx] = (double)v2[tmpIndex];
			VOut[tmpIndex + 2 * nPx] = (double)v3[tmpIndex];
			VOut[tmpIndex + 3 * nPx] = (double)v4[tmpIndex];
			VOut[tmpIndex + 4 * nPx] = (double)v5[tmpIndex];
			VOut[tmpIndex + 5 * nPx] = (double)v6[tmpIndex];
		}
	}


	delete[] v1;
	delete[] v2;
	delete[] v3;
	delete[] v4;
	delete[] v5;
	delete[] v6;

	delete[] tauAquad;
	delete[] teiler1;
	delete[] teiler2;

	delete[] ux;
	delete[] uy;
	delete[] ut;

	delete[] tauUx;
	delete[] tauUy;

	delete[] y1;
	delete[] y2;
	delete[] y3;
	delete[] y4;
	delete[] y5;
	delete[] y6;
	delete[] y7;
	delete[] y8;

	delete[] Kx1;
	delete[] Kx2;
	delete[] Kx3;
	delete[] Kx4;
	delete[] Kx5;
	delete[] Kx6;
	delete[] Kx7;
	delete[] Kx8;

	delete[] Kty1;
	delete[] Kty2;
	delete[] Kty3;
	delete[] Kty4;
	delete[] Kty5;
	delete[] Kty6;
}