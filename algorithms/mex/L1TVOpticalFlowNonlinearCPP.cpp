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

#include <algorithm>

using namespace std;

void doWarp(const float *image1, const float *image2, float *v1, float *v2, float *ut, float *ux, float *uy, float *uxt, float *uyt, float *uxx, float *uxy, float *uyx, float *uyy, const size_t *sizeMat);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//L1TVOpticalFlow(image1,image2,tol,lambda,maxIterations,norm,inputV,inputY,stepsize,discretization,numberOfWarps)
{
	double *u1 = mxGetPr(prhs[0]);
	double *u2 = mxGetPr(prhs[1]);

	const float tol = (float)mxGetScalar(prhs[2]);
	const float lambda = (float)mxGetScalar(prhs[3]);

	//mexPrintf("tol: %f, lambda: %f; ", tol, lambda);

	const int maxIterations = (int)mxGetScalar(prhs[4]);

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
    
    int gradientConstancy = 0;
	if (nrhs > 12)
	{
		gradientConstancy = (int)mxGetScalar(prhs[12]);
	}
	
	const int nPx = (int)(sizeImage[0] * sizeImage[1]);

	const size_t sizeY[2] = { 7 * nPx, 1 };

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
    
	float* v1Old = new float[nPx];
	float* v2Old = new float[nPx];

	float* image1f = new float[nPx];
	float* image2f = new float[nPx];

	float* ux = new float[nPx];
	float* uy = new float[nPx];
	float* ut = new float[nPx];
    
    float* uxx = new float[nPx];
    float* uxy = new float[nPx];
    float* uyx = new float[nPx];
    float* uyy = new float[nPx];
    
    float* uxt = new float[nPx];
    float* uyt = new float[nPx];
    
	float* y11 = new float[nPx];
	float* y12 = new float[nPx];
	float* y21 = new float[nPx];
	float* y22 = new float[nPx];
    float* y3 = new float[nPx];
    
    float* y4 = new float[nPx];
    float* y5 = new float[nPx];
    
	float* y11Old = new float[nPx];
	float* y12Old = new float[nPx];
	float* y21Old = new float[nPx];
	float* y22Old = new float[nPx];
    float* y3Old = new float[nPx];
    
    float* y4Old = new float[nPx];
    float* y5Old = new float[nPx];

	float* Kty1 = new float[nPx];
	float* Kty2 = new float[nPx];
    
	float* Kty1Old = new float[nPx];
	float* Kty2Old = new float[nPx];

	float* Kx11 = new float[nPx];
	float* Kx12 = new float[nPx];
	float* Kx21 = new float[nPx];
	float* Kx22 = new float[nPx];
    float* Kx3 = new float[nPx];
    float* Kx4 = new float[nPx];
    float* Kx5 = new float[nPx];
    
	float* Kx11Old = new float[nPx];
	float* Kx12Old = new float[nPx];
	float* Kx21Old = new float[nPx];
	float* Kx22Old = new float[nPx];
    float* Kx3Old = new float[nPx];
    float* Kx4Old = new float[nPx];
    float* Kx5Old = new float[nPx];
    
	float sigma1 = myMin(stepsize[1]/2.0f,stepsize[2]/2.0f);
    float* sigma2 = new float[nPx];
    
    //for gradient constancy
    float* sigma3 = new float[nPx];
    float* sigma4 = new float[nPx];
    
    float* tau1 = new float[nPx];
	float* tau2 = new float[nPx];
    
    int * tableI = new int[nPx];
    int * tableJ = new int[nPx];

    //Huber Factor
	
	const float huberFactor = 1.0f / (1.0f + sigma1* huberEpsilon / lambda);

	//residuals
	float p = 0.0f;
	float d = 0.0f;
	
	#pragma omp parallel for
	for (int j = 0; j < sizeImage[1]; ++j)
	{
		for (int i = 0; i < sizeImage[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeImage, i, j);
            
            tableI[tmpIndex] = i;
            tableJ[tmpIndex] = j;

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
                y3[tmpIndex] = (float)inputY[4 * nPx + tmpIndex];
                
                //for gradient constancy
                y4[tmpIndex] = (float)inputY[5 * nPx + tmpIndex];
                y5[tmpIndex] = (float)inputY[6 * nPx + tmpIndex];
			}
			else
			{
				y11[tmpIndex] = 0.0f;
				y12[tmpIndex] = 0.0f;
				y21[tmpIndex] = 0.0f;
				y22[tmpIndex] = 0.0f;
                y3[tmpIndex] = 0.0f;
                
                //for gradient constancy
                y4[tmpIndex] = 0.0f;
                y5[tmpIndex] = 0.0f;
			}

			Kx11[tmpIndex] = 0.0f;
			Kx12[tmpIndex] = 0.0f;
			Kx21[tmpIndex] = 0.0f;
			Kx22[tmpIndex] = 0.0f;
            Kx3[tmpIndex] = 0.0f;
            
            //for gradient constancy
            Kx4[tmpIndex] = 0.0f;
            Kx5[tmpIndex] = 0.0f;
		}
	}

	//do k warpings
	for (int k = 0; k < numberOfWarps; ++k)
	{
		doWarp(image1f, image2f, v1, v2, ut, ux, uy,uxt,uyt,uxx,uxy,uyx,uyy, sizeImage);
        
        #pragma omp parallel for
		for (int i = 0; i < nPx; ++i)
		{
            //adjust step sizes
            
            tau1[i] = 4.0f/std::min(stepsize[1],stepsize[2]) + myAbs(ux[i]);
            tau2[i] = 4.0f/std::min(stepsize[1],stepsize[2]) + myAbs(uy[i]);
            
            if (gradientConstancy>0)
            {
                tau1[i] += std::abs(uxx[i]) + std::abs(uyx[i]);
                tau2[i] += std::abs(uxy[i]) + std::abs(uyy[i]);
            }
            
            sigma2[i] = std::abs(ux[i]) + std::abs(uy[i]);
            
            //for gradient constancy
            sigma3[i] = std::abs(uxx[i]) + std::abs(uxy[i]);
            sigma4[i] = std::abs(uyx[i]) + std::abs(uyy[i]);
            
            tau1[i] = 1.0f / tau1[i];
            tau2[i] = 1.0f / tau2[i];
            sigma2[i] = 1.0f / sigma2[i];
            
            //for gradient constancy
            sigma3[i] = 1.0f / sigma3[i];
            sigma4[i] = 1.0f / sigma4[i];
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
            
            #pragma omp parallel for reduction(+:p)
            for (int i = 0; i < nPx; ++i)
            {
                if (iterations % 50 == 0)
                {
                    v1Old[i] = v1[i];
                    v2Old[i] = v2[i];

                    Kty1Old[i] = Kty1[i];
                    Kty2Old[i] = Kty2[i];
                }

                //transpose equals -div
                Kty1[i] = -stepsizeD[1] * dxm(y11, sizeImage, tableI[i], tableJ[i]) - stepsizeD[2] * dym(y12, sizeImage, tableI[i], tableJ[i]) + ux[i] * y3[i];
                Kty2[i] = -stepsizeD[1] * dxm(y21, sizeImage, tableI[i], tableJ[i]) - stepsizeD[2] * dym(y22, sizeImage, tableI[i], tableJ[i]) + uy[i] * y3[i];

                if (gradientConstancy>0)
                {
                    Kty1[i] += uxx[i] * y4[i] + uyx[i] * y5[i];
                    Kty2[i] += uxy[i] * y4[i] + uyy[i] * y5[i];
                }

                v1[i] -= tau1[i]*Kty1[i];
                v2[i] -= tau2[i]*Kty2[i];

                if (iterations % 50 == 0)
                {
                    //residuals
                    p += std::abs((v1Old[i] - v1[i]) / tau1[i] - Kty1Old[i] + Kty1[i])
                       + std::abs((v2Old[i] - v2[i]) / tau2[i] - Kty2Old[i] + Kty2[i]);
                }
			}

			//dual step
			#pragma omp parallel for reduction(+:d)
			for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
			{
                float y11Tilde, y12Tilde, y21Tilde, y22Tilde;

                if (iterations % 50 == 0)
                {
                    y11Old[tmpIndex] = y11[tmpIndex];
                    y12Old[tmpIndex] = y12[tmpIndex];
                    y21Old[tmpIndex] = y21[tmpIndex];
                    y22Old[tmpIndex] = y22[tmpIndex];

                    y3Old[tmpIndex] = y3[tmpIndex];
                }

                Kx11Old[tmpIndex] = Kx11[tmpIndex];
                Kx12Old[tmpIndex] = Kx12[tmpIndex];
                Kx21Old[tmpIndex] = Kx21[tmpIndex];
                Kx22Old[tmpIndex] = Kx22[tmpIndex];

                Kx3Old[tmpIndex] = Kx3[tmpIndex];

                Kx11[tmpIndex] = stepsizeD[1] * dxp(v1, sizeImage, tableI[tmpIndex], tableJ[tmpIndex]);
                Kx12[tmpIndex] = stepsizeD[2] * dyp(v1, sizeImage, tableI[tmpIndex], tableJ[tmpIndex]);
                Kx21[tmpIndex] = stepsizeD[1] * dxp(v2, sizeImage, tableI[tmpIndex], tableJ[tmpIndex]);
                Kx22[tmpIndex] = stepsizeD[2] * dyp(v2, sizeImage, tableI[tmpIndex], tableJ[tmpIndex]);

                Kx3[tmpIndex] = ux[tmpIndex] * v1[tmpIndex] + uy[tmpIndex] * v2[tmpIndex];

                if (gradientConstancy>0)
                {
                    Kx4Old[tmpIndex] = Kx4[tmpIndex];
                    Kx5Old[tmpIndex] = Kx5[tmpIndex];

                    Kx4[tmpIndex] = uxx[tmpIndex] * v1[tmpIndex] + uxy[tmpIndex] * v2[tmpIndex];
                    Kx5[tmpIndex] = uyx[tmpIndex] * v1[tmpIndex] + uyy[tmpIndex] * v2[tmpIndex];

                    y4Old[tmpIndex] = y4[tmpIndex];
                    y5Old[tmpIndex] = y5[tmpIndex];

                    y4[tmpIndex] = std::max(-1.0f,std::min(1.0f, y4[tmpIndex] + sigma3[tmpIndex]*(Kx4[tmpIndex] + Kx4[tmpIndex] - Kx4Old[tmpIndex] + uxt[tmpIndex])));
                    y5[tmpIndex] = std::max(-1.0f,std::min(1.0f, y5[tmpIndex] + sigma4[tmpIndex]*(Kx5[tmpIndex] + Kx5[tmpIndex] - Kx5Old[tmpIndex] + uyt[tmpIndex])));
                }

                if (typeNorm == 4) // Huber
                {
                    y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old[tmpIndex])) * huberFactor;
                    y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old[tmpIndex])) * huberFactor;
                    y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old[tmpIndex])) * huberFactor;
                    y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old[tmpIndex])) * huberFactor;
                }
                else
                {
                    y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old[tmpIndex]));
                    y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old[tmpIndex]));
                    y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old[tmpIndex]));
                    y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old[tmpIndex]));
                }

                float divisor1 = std::max(1.0f, std::sqrt(y11Tilde*y11Tilde + y12Tilde*y12Tilde) / lambda);
                float divisor2 = std::max(1.0f, std::sqrt(y21Tilde*y21Tilde + y22Tilde*y22Tilde) / lambda);

                y11[tmpIndex] = y11Tilde / divisor1;
                y12[tmpIndex] = y12Tilde / divisor1;
                y21[tmpIndex] = y21Tilde / divisor2;
                y22[tmpIndex] = y22Tilde / divisor2;

                y3[tmpIndex] = std::max(-1.0f,std::min(1.0f, y3[tmpIndex] + sigma2[tmpIndex]*(Kx3[tmpIndex] + Kx3[tmpIndex] - Kx3Old[tmpIndex] + ut[tmpIndex])));

                if (iterations % 50 == 0)
                {
                    d += std::abs((y11Old[tmpIndex] - y11[tmpIndex]) / sigma1 - Kx11Old[tmpIndex] + Kx11[tmpIndex]) +
                        std::abs((y12Old[tmpIndex] - y12[tmpIndex]) / sigma1 - Kx12Old[tmpIndex] + Kx12[tmpIndex]) +
                        std::abs((y21Old[tmpIndex] - y21[tmpIndex]) / sigma1 - Kx21Old[tmpIndex] + Kx21[tmpIndex]) +
                        std::abs((y22Old[tmpIndex] - y22[tmpIndex]) / sigma1 - Kx22Old[tmpIndex] + Kx22[tmpIndex]) + 
                        std::abs((y3Old[tmpIndex] - y3[tmpIndex]) / sigma2[tmpIndex] - Kx3Old[tmpIndex] + Kx3[tmpIndex]);

                        if (gradientConstancy>0)
                        {
                            d += std::abs((y4Old[tmpIndex] - y4[tmpIndex]) / sigma3[tmpIndex] - Kx4Old[tmpIndex] + Kx4[tmpIndex]) + 
                                 std::abs((y5Old[tmpIndex] - y5[tmpIndex]) / sigma4[tmpIndex] - Kx5Old[tmpIndex] + Kx5[tmpIndex]);
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
	for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
	{
        YOut[tmpIndex + 0 * nPx] = (double)y11[tmpIndex];
        YOut[tmpIndex + 1 * nPx] = (double)y12[tmpIndex];
        YOut[tmpIndex + 2 * nPx] = (double)y21[tmpIndex];
        YOut[tmpIndex + 3 * nPx] = (double)y22[tmpIndex];
        YOut[tmpIndex + 4 * nPx] = (double)y3[tmpIndex];

        YOut[tmpIndex + 5 * nPx] = (double)y4[tmpIndex];
        YOut[tmpIndex + 6 * nPx] = (double)y5[tmpIndex];

        Outv1[tmpIndex] = (double)v1[tmpIndex];
        Outv2[tmpIndex] = (double)v2[tmpIndex];
	}
    
    delete[] tableI;
    delete[] tableJ;
    
    delete[] image1f;
    delete[] image2f;
    
    delete[] sigma2;
    delete[] sigma3;
    delete[] sigma4;
    delete[] tau1;
    delete[] tau2;

	delete[] v1;
	delete[] v2;
	delete[] v1Old;
	delete[] v2Old;
    
	delete[] ux;
	delete[] uy;
	delete[] ut;

	delete[] uxx;
    delete[] uxy;
	delete[] uyx;
    delete[] uyy;
    
    delete[] uxt;
    delete[] uyt;
    
	delete[] y11;
	delete[] y12;
	delete[] y21;
	delete[] y22;
    delete[] y3;
    delete[] y4;
    delete[] y5;
    
	delete[] y11Old;
	delete[] y12Old;
	delete[] y21Old;
	delete[] y22Old;
    delete[] y3Old;
    delete[] y4Old;
    delete[] y5Old;
    
	delete[] Kty1;
	delete[] Kty2;
    
	delete[] Kty1Old;
	delete[] Kty2Old;

	delete[] Kx11;
	delete[] Kx12;
	delete[] Kx21;
	delete[] Kx22;
    delete[] Kx3;
    delete[] Kx4;
    delete[] Kx5;
    
	delete[] Kx11Old;
	delete[] Kx12Old;
	delete[] Kx21Old;
	delete[] Kx22Old;
    delete[] Kx3Old;
    delete[] Kx4Old;
    delete[] Kx5Old;
}

void doWarp(const float *image1, const float *image2, float *v1, float *v2, float *ut, float *ux, float *uy, float *uxt, float *uyt, float *uxx, float *uxy, float *uyx, float *uyy, const size_t *sizeMat)
{
	int nPx = (int)(sizeMat[0] * sizeMat[1]);
    
    float* origX = new float[nPx];
    float* origY = new float[nPx];
    float* shiftX = new float[nPx];
    float* shiftY = new float[nPx];
    
    float* shiftXp = new float[nPx];
    float* shiftXm = new float[nPx];
    float* shiftYp = new float[nPx];
    float* shiftYm = new float[nPx];
    
    float* origXp = new float[nPx];
    float* origXm = new float[nPx];
    float* origYp = new float[nPx];
    float* origYm = new float[nPx];
    
    
    //create indicator grids
	#pragma omp parallel for
	for (int j = 0; j < sizeMat[1]; ++j)
	{
		for (int i = 0; i < sizeMat[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeMat, i, j);
            
            origX[tmpIndex] = float(i);
            origY[tmpIndex] = float(j);
            shiftX[tmpIndex] = float(i) + v2[tmpIndex];
            shiftY[tmpIndex] = float(j) + v1[tmpIndex];
            
            shiftXp[tmpIndex] = myMin(sizeMat[0],myMax(0,shiftX[tmpIndex] + 0.5f));
            shiftXm[tmpIndex] = myMin(sizeMat[0],myMax(0,shiftX[tmpIndex] - 0.5f));
            shiftYp[tmpIndex] = myMin(sizeMat[1],myMax(0,shiftY[tmpIndex] + 0.5f));
            shiftYm[tmpIndex] = myMin(sizeMat[1],myMax(0,shiftY[tmpIndex] - 0.5f));
            
            origXp[tmpIndex] = myMin(sizeMat[0],myMax(0,origX[tmpIndex] + 0.5f));
            origXm[tmpIndex] = myMin(sizeMat[0],myMax(0,origX[tmpIndex] - 0.5f));
            origYp[tmpIndex] = myMin(sizeMat[1],myMax(0,origY[tmpIndex] + 0.5f));
            origYm[tmpIndex] = myMin(sizeMat[1],myMax(0,origY[tmpIndex] - 0.5f));
        }
    }

    #pragma omp parallel for
	for (int j = 0; j < sizeMat[1]; ++j)
	{
		for (int i = 0; i < sizeMat[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeMat, i, j);
            
			uy[tmpIndex] = cubicInterpolation(image2,shiftXp[tmpIndex],shiftY[tmpIndex], sizeMat) - cubicInterpolation(image2,shiftXm[tmpIndex],shiftY[tmpIndex], sizeMat);
			ux[tmpIndex] = cubicInterpolation(image2,shiftX[tmpIndex],shiftYp[tmpIndex], sizeMat) - cubicInterpolation(image2,shiftX[tmpIndex],shiftYm[tmpIndex], sizeMat);
			ut[tmpIndex] = cubicInterpolation(image2,shiftX[tmpIndex],shiftY[tmpIndex], sizeMat) - image1[tmpIndex] - ux[tmpIndex] * v1[tmpIndex] - uy[tmpIndex] * v2[tmpIndex];
		}
	}
    
    #pragma omp parallel for
	for (int j = 0; j < sizeMat[1]; ++j)
	{
		for (int i = 0; i < sizeMat[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeMat, i, j);

            float image1x = cubicInterpolation(image1,origX[tmpIndex],origYp[tmpIndex], sizeMat) - cubicInterpolation(image1,origX[tmpIndex],origYm[tmpIndex], sizeMat);
            float image1y = cubicInterpolation(image1,origXp[tmpIndex],origY[tmpIndex], sizeMat) - cubicInterpolation(image1,origXm[tmpIndex],origY[tmpIndex], sizeMat);
            
			uyy[tmpIndex] = cubicInterpolation(uy,origXp[tmpIndex],origY[tmpIndex], sizeMat) - cubicInterpolation(uy,origXm[tmpIndex],origY[tmpIndex], sizeMat);
            uyx[tmpIndex] = cubicInterpolation(uy,origX[tmpIndex],origYp[tmpIndex], sizeMat) - cubicInterpolation(uy,origX[tmpIndex],origYm[tmpIndex], sizeMat);
            uxy[tmpIndex] = cubicInterpolation(ux,origXp[tmpIndex],origY[tmpIndex], sizeMat) - cubicInterpolation(ux,origXm[tmpIndex],origY[tmpIndex], sizeMat);
            uxx[tmpIndex] = cubicInterpolation(ux,origX[tmpIndex],origYp[tmpIndex], sizeMat) - cubicInterpolation(ux,origX[tmpIndex],origYm[tmpIndex], sizeMat);
            
            uxt[tmpIndex] = ux[tmpIndex]-image1x - uxx[tmpIndex] * v1[tmpIndex] - uxy[tmpIndex] * v2[tmpIndex];
            uyt[tmpIndex] = uy[tmpIndex]-image1y - uyx[tmpIndex] * v1[tmpIndex] - uyy[tmpIndex] * v2[tmpIndex];
            
            //check out of image:
            //due to derivative calculation, go one pixel inside
            if (shiftX[tmpIndex]<1.0f || shiftY[tmpIndex]<1.0f || shiftX[tmpIndex] > (sizeMat[0]-2.0f) || shiftY[tmpIndex] > (sizeMat[1]-2.0f)
            || origX[tmpIndex]<1.0f || origY[tmpIndex]<1.0f || origX[tmpIndex] > (sizeMat[0]-2.0f) || origY[tmpIndex] > (sizeMat[1]-2.0f))
            {
                uy[tmpIndex] = 0.0f;
                ux[tmpIndex] = 0.0f;
                ut[tmpIndex] = 0.0f;
                
                uxx[tmpIndex] = 0.0f;
                uxy[tmpIndex] = 0.0f;
                uyx[tmpIndex] = 0.0f;
                uyy[tmpIndex] = 0.0f;
                
                uxt[tmpIndex] = 0.0f;
                uyt[tmpIndex] = 0.0f;
            }
		}
	}
    
    delete[] origX;
    delete[] origY;
    delete[] shiftX;
    delete[] shiftY;

    delete[] shiftXp;
    delete[] shiftXm;
    delete[] shiftYp;
    delete[] shiftYm;

    delete[] origXp;
    delete[] origXm;
    delete[] origYp;
    delete[] origYm;
}