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
% Date: 2015-06-17

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
//L1L2OpticalFlow(u1,u2,tol,lambda,ux,uy,ut)
{
	double *u1 = mxGetPr(prhs[0]);
	double *u2 = mxGetPr(prhs[1]);

	float tol = (float)mxGetScalar(prhs[2]);
	float lambda = (float)mxGetScalar(prhs[3]);

	mexPrintf("tol: %f, lambda: %f; ", tol, lambda);

	int maxIterations = (int)mxGetScalar(prhs[4]);

	const size_t *sizeImage = mxGetDimensions(prhs[0]);

	double *inputV = 0;
	double *inputY = 0;

	int typeNorm = 1;
	if (nrhs > 5)
	{
		typeNorm = (int)mxGetScalar(prhs[5]);
	}

	if (nrhs > 6)
	{
		inputV = mxGetPr(prhs[6]);
	}

	if (nrhs > 7)
	{
		inputY = mxGetPr(prhs[7]);
	}

	float stepsize[3] = { 1.0f, 1.0f, 1.0f };
	if (nrhs > 8)
	{
		double *tmpStepsize = mxGetPr(prhs[8]);

		stepsize[0] = (float)tmpStepsize[0];
		stepsize[1] = (float)tmpStepsize[1];
		stepsize[2] = (float)tmpStepsize[2];
	}
	float stepsizeD[3] = { 1.0f / stepsize[0], 1.0f / stepsize[1], 1.0f / stepsize[2] };

	mexPrintf("Stepsize: dt=%f, dx=%f, dy=%f; ", stepsize[0], stepsize[1], stepsize[2]);

	int discretizationScheme = 1; //forward time, central space
	if (nrhs > 9)
	{
		// 2: forward time, upwind for space (only in combination with given v)
		discretizationScheme = (int)mxGetScalar(prhs[9]);
	}

	if (discretizationScheme == 2)
	{
		mexPrintf("Upwind Discretization chosen");
	}
	else
	{
		mexPrintf("Regular centered Discretization chosen");
	}

	mexPrintf("\n");

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

	float *tauBetaQuad = new float[nPx];

	float* tauUx = new float[nPx];
	float* tauUy = new float[nPx];

	float* teiler1 = new float[nPx];
	float* teiler2 = new float[nPx];

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

	float tau = 1.0f / sqrt(8.0f)*stepsize[1];
	float sigma = tau;
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

			if (discretizationScheme == 2)
			{
				ut[tmpIndex] = (float)(stepsizeD[0] * (u2[tmpIndex] - u1[tmpIndex]));

				if (v1[tmpIndex] > 0)
				{
					if (j > 0)
					{
						ux[tmpIndex] = (float)(stepsizeD[1] * (u1[index2DtoLinear(sizeImage, i, j)] - u1[index2DtoLinear(sizeImage, i, j - 1)]));
					}
					else
					{
						ux[tmpIndex] = 0.0f;
					}
				}
				else
				{
					if (j < sizeImage[1] - 1)
					{
						ux[tmpIndex] = (float)(stepsizeD[1] * (u1[index2DtoLinear(sizeImage, i, j + 1)] - u1[index2DtoLinear(sizeImage, i, j)]));
					}
					else
					{
						ux[tmpIndex] = 0.0f;
					}
				}

				if (v2[tmpIndex] > 0)
				{
					if (i > 0)
					{
						uy[tmpIndex] = (float)(stepsizeD[2] * (u1[index2DtoLinear(sizeImage, i, j)] - u1[index2DtoLinear(sizeImage, i - 1, j)]));
					}
					else
					{
						uy[tmpIndex] = 0.0f;
					}
				}
				else
				{
					if (i < sizeImage[0] - 1)
					{
						uy[tmpIndex] = (float)(stepsizeD[2] * (u1[index2DtoLinear(sizeImage, i + 1, j)] - u1[index2DtoLinear(sizeImage, i, j)]));
					}
					else
					{
						uy[tmpIndex] = 0.0f;
					}
				}
			}
			else
			{
				//forward time, central space
				ut[tmpIndex] = (float)(stepsizeD[0] * (u2[tmpIndex] - u1[tmpIndex]));

				if (i>0 && i < sizeImage[0] - 1)
				{
					uy[tmpIndex] = (float)(stepsizeD[2] * 0.5f * (u1[index2DtoLinear(sizeImage, i + 1, j)] - u1[index2DtoLinear(sizeImage, i - 1, j)]));
				}
				else
				{
					uy[tmpIndex] = 0.0;
				}
				
				if (j>0 && j < sizeImage[1] - 1)
				{
					ux[tmpIndex] = (float)(stepsizeD[2] * 0.5f * (u1[index2DtoLinear(sizeImage, i, j + 1)] - u1[index2DtoLinear(sizeImage, i, j - 1)]));
				}
				else
				{
					ux[tmpIndex] = 0.0f;
				}
			}

			float betaQuad = myMax(0.0f, (float)(ux[tmpIndex] * ux[tmpIndex]) + (float)(uy[tmpIndex] * uy[tmpIndex]));
			tauBetaQuad[tmpIndex] = tau*betaQuad;

			teiler1[tmpIndex] = ux[tmpIndex] / betaQuad;
			teiler2[tmpIndex] = uy[tmpIndex] / betaQuad;

			tauUx[tmpIndex] = tau*ux[tmpIndex];
			tauUy[tmpIndex] = tau*uy[tmpIndex];
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
				Kty1[tmpIndex] = -stepsizeD[1] * (dxm(y11, sizeImage, i, j) + dym(y12, sizeImage, i, j));
				Kty2[tmpIndex] = -stepsizeD[1] * (dxm(y21, sizeImage, i, j) + dym(y22, sizeImage, i, j));

				float v1Old = v1[tmpIndex];
				float v2Old = v2[tmpIndex];
				
				float v1Hat = v1Old - tau*Kty1[tmpIndex];
				float v2Hat = v2Old - tau*Kty2[tmpIndex];

				float rho = ut[tmpIndex] + ux[tmpIndex] * v1Hat + uy[tmpIndex] * v2Hat;

				if (rho <= -tauBetaQuad[tmpIndex])
				{
					v1[tmpIndex] = v1Hat + tauUx[tmpIndex];
					v2[tmpIndex] = v2Hat + tauUy[tmpIndex];
				}
				else if (rho >= tauBetaQuad[tmpIndex])
				{
					v1[tmpIndex] = v1Hat - tauUx[tmpIndex];
					v2[tmpIndex] = v2Hat - tauUy[tmpIndex];
				}
				else
				{
					v1[tmpIndex] = v1Hat - teiler1[tmpIndex] * rho;
					v2[tmpIndex] = v2Hat - teiler2[tmpIndex] * rho;
				}

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

				y11[tmpIndex] = ssl * (y11[tmpIndex] + sigma*(2 * Kx11[tmpIndex] - Kx11Old));
				y12[tmpIndex] = ssl * (y12[tmpIndex] + sigma*(2 * Kx12[tmpIndex] - Kx12Old));
				y21[tmpIndex] = ssl * (y21[tmpIndex] + sigma*(2 * Kx21[tmpIndex] - Kx21Old));
				y22[tmpIndex] = ssl * (y22[tmpIndex] + sigma*(2 * Kx22[tmpIndex] - Kx22Old));

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

	delete v1;
	delete v2;

	delete ux;
	delete uy;
	delete ut;

	delete tauBetaQuad;

	delete tauUx;
	delete tauUy;

	delete teiler1;
	delete teiler2;

	delete y11;
	delete y12;
	delete y21;
	delete y22;

	delete Kty1;
	delete Kty2;

	delete Kx11;
	delete Kx12;
	delete Kx21;
	delete Kx22;
}

