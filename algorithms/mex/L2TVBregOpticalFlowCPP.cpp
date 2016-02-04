#include "mex.h"
#include "math.h"
#include <omp.h>
#include "tools.h"

#include <vector> 
using namespace std;

//#define INCREMENT(x) x++

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//L2TVBregOpticalFlow(u1,u2,tol,lambda,maxIterations,numBreg)
{
	double *u1 = mxGetPr(prhs[0]);
	double *u2 = mxGetPr(prhs[1]);

	float tol = (float)mxGetScalar(prhs[2]);
	float lambda = (float)mxGetScalar(prhs[3]);

	int maxIterations = (int)mxGetScalar(prhs[4]);

	int numBreg = (int)mxGetScalar(prhs[5]);

	const mwSize *sizeImage = mxGetDimensions(prhs[0]);

	int typeNorm = 2;

	double *inputV = 0;
	double *inputY = 0;
	double *inputB = 0;

	if (nrhs > 6)
	{
		inputV = mxGetPr(prhs[6]);
	}

	if (nrhs > 7)
	{
		inputY = mxGetPr(prhs[7]);
	}
	if (nrhs > 8)
	{
		inputB = mxGetPr(prhs[8]);
	}

	int nPx = (int)(sizeImage[0] * sizeImage[1]);

	const mwSize sizeY[2] = {4 * nPx, 1 };
	const mwSize sizeB[2] = {2 * nPx, 1 };

	// Output v1
	plhs[0] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
	double *Outv1 = mxGetPr(plhs[0]);

	// Output v2
	plhs[1] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
	double *Outv2 = mxGetPr(plhs[1]);

	// Output  Y
	plhs[2] = mxCreateNumericArray(2, sizeY, mxDOUBLE_CLASS, mxREAL);
	double *YOut = mxGetPr(plhs[2]);

	plhs[3] = mxCreateNumericArray(2, sizeB, mxDOUBLE_CLASS, mxREAL);
	double *BOut = mxGetPr(plhs[3]);

	float* v1 = new float[nPx];
	float* v2 = new float[nPx];

	float* uxut = new float[nPx];
	float* uyut = new float[nPx];

	float* c1 = new float[nPx];
	float* c2 = new float[nPx];
	float* c3 = new float[nPx];
	float* teiler = new float[nPx];

	float* ux = new float[nPx];
	float* uy = new float[nPx];
	float* ut = new float[nPx];

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

	float* breg1 = new float[nPx];
	float* breg2 = new float[nPx];

	float* q = new float[nPx];


	float tau = 0.25f;
	float sigma = 0.5f;

	int i = 0; 
	int j = 0;
	
	int tmpIndex = 0;

	//residuals
	float p = 0.0f;
	float d = 0.0f;


	#pragma omp parallel for private(i,tmpIndex)
	for (j = 0; j < sizeImage[1]; ++j)
	{
		for (i = 0; i < sizeImage[0]; ++i)
		{
			tmpIndex = index2DtoLinear(sizeImage, i, j);

			q[tmpIndex] = 0.0f;


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
				v1[tmpIndex] = 0.0f;
				v2[tmpIndex] = 0.0f;
			}

			Kty1[tmpIndex] = 0.0f;
			Kty2[tmpIndex] = 0.0f;

			if (nrhs > 7)
			{
				y11[tmpIndex] = (float)inputY[tmpIndex];
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

			if (nrhs > 8)
			{
				breg1[tmpIndex] = (float)inputB[0 * nPx + tmpIndex];
				breg2[tmpIndex] = (float)inputB[1 * nPx + tmpIndex];
			}
			else
			{
				breg1[tmpIndex] = 0.0f;
				breg2[tmpIndex] = 0.0f;
			}

			Kx11[tmpIndex] = 0.0f;
			Kx12[tmpIndex] = 0.0f;
			Kx21[tmpIndex] = 0.0f;
			Kx22[tmpIndex] = 0.0f;
		}

	}


	for (int breg = 0; breg <= numBreg; ++breg)
	{
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
					Kty1[tmpIndex] = -dxm(y11, sizeImage, i, j) - dym(y12, sizeImage, i, j);
					Kty2[tmpIndex] = -dxm(y21, sizeImage, i, j) - dym(y22, sizeImage, i, j);

					float b1 = v1[tmpIndex] - tau*(Kty1[tmpIndex] + uxut[tmpIndex] - lambda * breg1[tmpIndex]);
					float b2 = v2[tmpIndex] - tau*(Kty2[tmpIndex] + uyut[tmpIndex] - lambda * breg2[tmpIndex]);

					float v1Old = v1[tmpIndex];
					float v2Old = v2[tmpIndex];

					v1[tmpIndex] = (b1 * c3[tmpIndex] - c2[tmpIndex] * b2) * teiler[tmpIndex];
					v2[tmpIndex] = (b2 * c1[tmpIndex] - c2[tmpIndex] * b1) * teiler[tmpIndex];

					if (iterations % 50 == 0)
					{
						//residuals
						p += myAbs((v1Old - v1[tmpIndex]) / tau - Kty1Old + Kty1[tmpIndex])
							+ myAbs((v2Old - v2[tmpIndex]) / tau - Kty2Old + Kty2[tmpIndex]);
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

					float y11Old = y11[tmpIndex];
					float y12Old = y12[tmpIndex];
					float y21Old = y21[tmpIndex];
					float y22Old = y22[tmpIndex];

					Kx11[tmpIndex] = dxp(v1, sizeImage, i, j);
					Kx12[tmpIndex] = dyp(v1, sizeImage, i, j);
					Kx21[tmpIndex] = dxp(v2, sizeImage, i, j);
					Kx22[tmpIndex] = dyp(v2, sizeImage, i, j);

					float y11Tilde = y11[tmpIndex] + sigma*(2 * Kx11[tmpIndex] - Kx11Old);
					float y12Tilde = y12[tmpIndex] + sigma*(2 * Kx12[tmpIndex] - Kx12Old);
					float y21Tilde = y21[tmpIndex] + sigma*(2 * Kx21[tmpIndex] - Kx21Old);
					float y22Tilde = y22[tmpIndex] + sigma*(2 * Kx22[tmpIndex] - Kx22Old);

					float divisor1 = myMax(1.0f, sqrtf(y11Tilde*y11Tilde + y12Tilde*y12Tilde) / lambda);
					float divisor2 = myMax(1.0f, sqrtf(y21Tilde*y21Tilde + y22Tilde*y22Tilde) / lambda);

					y11[tmpIndex] = y11Tilde / divisor1;
					y12[tmpIndex] = y12Tilde / divisor1;
					y21[tmpIndex] = y21Tilde / divisor2;
					y22[tmpIndex] = y22Tilde / divisor2;

					if (iterations % 50 == 0)
					{
						d += myAbs((y11Old - y11[tmpIndex]) / sigma - Kx11Old + Kx11[tmpIndex]) +
							myAbs((y12Old - y12[tmpIndex]) / sigma - Kx12Old + Kx12[tmpIndex]) +
							myAbs((y21Old - y21[tmpIndex]) / sigma - Kx21Old + Kx21[tmpIndex]) +
							myAbs((y22Old - y22[tmpIndex]) / sigma - Kx22Old + Kx22[tmpIndex]);
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

		if (breg < numBreg)
		{
			for (int j = 0; j < sizeImage[1]; ++j)
			{
				for (int i = 0; i < sizeImage[0]; ++i)
				{
					tmpIndex = index2DtoLinear(sizeImage, i, j);

					breg1[tmpIndex] = breg1[tmpIndex] - (1.0f / lambda) * (ut[tmpIndex] + ux[tmpIndex] * v1[tmpIndex] + uy[tmpIndex] * v2[tmpIndex]) * ux[tmpIndex];
					breg2[tmpIndex] = breg2[tmpIndex] - (1.0f / lambda) * (ut[tmpIndex] + ux[tmpIndex] * v1[tmpIndex] + uy[tmpIndex] * v2[tmpIndex]) * uy[tmpIndex];

					//q[tmpIndex] = (b1[tmpIndex] + b2[tmpIndex]) / sqrtf(ux[tmpIndex] * ux[tmpIndex] + uy[tmpIndex] * uy[tmpIndex] + 1e-6);
				}
			}
		}
		mexEvalString("drawnow;");

	}

	//write output
	#pragma omp parallel for private(i,tmpIndex)
	for (j = 0; j < sizeImage[1]; ++j)
	{
		for (i = 0; i < sizeImage[0]; ++i)
		{
			tmpIndex = index2DtoLinear(sizeImage, i, j);

			YOut[tmpIndex] = (double)y11[tmpIndex];
			YOut[tmpIndex + nPx] = (double)y12[tmpIndex];
			YOut[tmpIndex + 2 * nPx] = (double)y21[tmpIndex];
			YOut[tmpIndex + 3 * nPx] = (double)y22[tmpIndex];

			BOut[tmpIndex + 0 * nPx] = (double)breg1[tmpIndex];
			BOut[tmpIndex + 1 * nPx] = (double)breg2[tmpIndex];

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

	delete[] breg1;
	delete[] breg2;
	delete[] q;
}

