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
#include "linearInterpolation.h"
#include "cubicInterpolation.h"
#include <stdio.h>
#include <sys/types.h>

#include <cstddef>
using namespace std;

void doWarp(const float *image1, const float *image2, float *v1, float *v2, float *ut, float *ux, float *uy, const size_t *sizeMat);
//float interpolateLinear(const double *image, float x, float y, const size_t *sizeMat);
//float interpolateLinear2(const float *image, float x, float y, const size_t *sizeMat);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//L1TVOpticalFlow(image1,image2,tol,lambda,maxIterations,norm,inputV,inputY,stepsize,discretization,numberOfWarps)
{
	double *u1 = mxGetPr(prhs[0]);
	double *u2 = mxGetPr(prhs[1]);

	float tol = (float)mxGetScalar(prhs[2]);
	float lambda = (float)mxGetScalar(prhs[3]);

	//mexPrintf("tol: %f, lambda: %f; ", tol, lambda);

	int maxIterations = (int)mxGetScalar(prhs[4]);

	const size_t *sizeImage = mxGetDimensions(prhs[0]);

	double *inputV = 0;
	double *inputY = 0;

	int typeNorm = 3;
	if (nrhs > 5)
	{
		typeNorm = (int)mxGetScalar(prhs[5]);
	}
	//mexPrintf("TypeNorm is %d\n", typeNorm);

	if (nrhs > 6)
	{
		inputV = mxGetPr(prhs[6]);
		//mexPrintf("Given input vector field\n");
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

	int numberOfWarps = 5;
	if (nrhs > 10)
	{
		
		numberOfWarps = (int)mxGetScalar(prhs[10]);
	}

	float huberEpsilon = 0.01f;
	if (nrhs > 11)
	{
		huberEpsilon = (float)mxGetScalar(prhs[11]);
	}
	//mexPrintf("huberEpsilon is %f\n", huberEpsilon);

	

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

	float* image1f = new float[nPx];
	float* image2f = new float[nPx];

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

	float tau = 0.25f;
	float sigma = 0.5f;
	float dTau = 1.0f / tau;
	float dSigma = 1.0f / sigma;

	//Huber Factor
	
	float huberFactor = 1.0f / (1.0f + sigma* huberEpsilon / lambda);

	//residuals
	float p = 0.0f;
	float d = 0.0f;
	
	#pragma omp parallel for
	for (int j = 0; j < sizeImage[1]; ++j)
	{
		for (int i = 0; i < sizeImage[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeImage, i, j);

			image1f[tmpIndex] = (float)u1[tmpIndex];
			image2f[tmpIndex] = (float)u2[tmpIndex];

			if (nrhs > 6)
			{
				v1[tmpIndex] = (float)inputV[tmpIndex];
				v2[tmpIndex] = (float)inputV[nPx + tmpIndex];
			}
			else
			{
				v1[tmpIndex] = 0.0f;
				v2[tmpIndex] = 0.0f;
			}

			Kty1[tmpIndex] = 0.0f;
			Kty2[tmpIndex] = 0.0f;

			if (nrhs > 7)
			{
				y11[tmpIndex] = (float)inputY[0 * nPx + tmpIndex];
				y12[tmpIndex] = (float)inputY[1 * nPx + tmpIndex];
				y21[tmpIndex] = (float)inputY[2 * nPx + tmpIndex];
				y22[tmpIndex] = (float)inputY[3 * nPx + tmpIndex];
			}
			else
			{
				y11[tmpIndex] = 0.0f;
				y12[tmpIndex] = 0.0f;
				y21[tmpIndex] = 0.0f;
				y22[tmpIndex] = 0.0f;
			}

			Kx11[tmpIndex] = 0.0f;
			Kx12[tmpIndex] = 0.0f;
			Kx21[tmpIndex] = 0.0f;
			Kx22[tmpIndex] = 0.0f;
		}
	}

	//do k warpings
	for (int k = 0; k < numberOfWarps; ++k)
	{
		doWarp(image1f, image2f, v1, v2, ut, ux, uy, sizeImage);

		for (int i = 0; i < nPx; ++i)
		{
			float betaQuad = myMax(0.0f, (float)(ux[i] * ux[i]) + (float)(uy[i] * uy[i]));
			tauBetaQuad[i] = tau*betaQuad;

			teiler1[i] = ux[i] / betaQuad;
			teiler2[i] = uy[i] / betaQuad;

			tauUx[i] = tau*ux[i];
			tauUy[i] = tau*uy[i];

			ut[i] = ut[i] - ux[i] * v1[i] - uy[i] * v2[i];
		}


		int iterations = 0;
		float err = 1.0f;

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

					float y11Tilde, y12Tilde, y21Tilde, y22Tilde;

					float Kx11Old = Kx11[tmpIndex];
					float Kx12Old = Kx12[tmpIndex];
					float Kx21Old = Kx21[tmpIndex];
					float Kx22Old = Kx22[tmpIndex];

					//need multiplication with step length!!!
					Kx11[tmpIndex] = stepsizeD[1] * dxp(v1, sizeImage, i, j);
					Kx12[tmpIndex] = stepsizeD[1] * dyp(v1, sizeImage, i, j);
					Kx21[tmpIndex] = stepsizeD[1] * dxp(v2, sizeImage, i, j);
					Kx22[tmpIndex] = stepsizeD[1] * dyp(v2, sizeImage, i, j);

					float y11Old = y11[tmpIndex];
					float y12Old = y12[tmpIndex];
					float y21Old = y21[tmpIndex];
					float y22Old = y22[tmpIndex];

					if (typeNorm == 4) // Huber
					{
						y11Tilde = (y11[tmpIndex] + sigma*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old)) * huberFactor;
						y12Tilde = (y12[tmpIndex] + sigma*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old)) * huberFactor;
						y21Tilde = (y21[tmpIndex] + sigma*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old)) * huberFactor;
						y22Tilde = (y22[tmpIndex] + sigma*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old)) * huberFactor;
					}
					else
					{
						y11Tilde = (y11[tmpIndex] + sigma*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old));
						y12Tilde = (y12[tmpIndex] + sigma*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old));
						y21Tilde = (y21[tmpIndex] + sigma*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old));
						y22Tilde = (y22[tmpIndex] + sigma*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old));
					}
					
					//Huber


					float norm1 = sqrtf(y11Tilde*y11Tilde + y12Tilde*y12Tilde);
					float norm2 = sqrtf(y21Tilde*y21Tilde + y22Tilde*y22Tilde);

					float divisor1 = myMax(1.0f, norm1 / lambda);
					float divisor2 = myMax(1.0f, norm2 / lambda);

					y11[tmpIndex] = y11Tilde / divisor1;
					y12[tmpIndex] = y12Tilde / divisor1;
					y21[tmpIndex] = y21Tilde / divisor2;
					y22[tmpIndex] = y22Tilde / divisor2;


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
				err = (d*d + p*p) / (float)nPx;
			}

			if (iterations % 1000 == 0)
			{
				mexPrintf("Iteration %d,Residual %e\n", iterations, err);
				mexEvalString("drawnow;");
			}
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

	delete[] tauBetaQuad;

	delete[] tauUx;
	delete[] tauUy;

	delete[] teiler1;
	delete[] teiler2;

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

void doWarp(const float *image1, const float *image2, float *v1, float *v2, float *ut, float *ux, float *uy, const size_t *sizeMat)
{
	int nPx = (int)(sizeMat[0] * sizeMat[1]);

	#pragma omp parallel for
	for (int j = 0; j < sizeMat[1]; ++j)
	{
		for (int i = 0; i < sizeMat[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeMat, i, j);

			float xp = cubicInterpolation(image2, (float)i + v2[tmpIndex], (float)j + v1[tmpIndex] + 0.5f, sizeMat);
			float xm = cubicInterpolation(image2, (float)i + v2[tmpIndex], (float)j + v1[tmpIndex] - 0.5f, sizeMat);

			float yp = cubicInterpolation(image2, (float)i + v2[tmpIndex] + 0.5f, (float)j + v1[tmpIndex], sizeMat);
			float ym = cubicInterpolation(image2, (float)i + v2[tmpIndex] - 0.5f, (float)j + v1[tmpIndex], sizeMat);

			float warp = cubicInterpolation(image2, (float)i + v2[tmpIndex], (float)j + v1[tmpIndex], sizeMat);

			uy[tmpIndex] = yp - ym;
			ux[tmpIndex] = xp - xm;
			
			ut[tmpIndex] = warp - image1[tmpIndex];
		}
	}

	//mexPrintf("Image 2 is %f; %f\n", grad1_2Pr[tmpIndex], interpolateLinear2(image2, g1 + v2[tmpIndex], g2 + v1[tmpIndex] - 0.5, sizeMat));
}