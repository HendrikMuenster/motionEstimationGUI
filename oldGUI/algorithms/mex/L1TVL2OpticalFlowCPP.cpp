 #include "mex.h"
#include "math.h"
#include <omp.h>
#include "tools.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//L1TVL2OpticalFlowCPP(u1,u2,tol,lambda0,lambda1,maxIterations)
{
	double *u1 = mxGetPr(prhs[0]);
	double *u2 = mxGetPr(prhs[1]);

	float tol = (float)mxGetScalar(prhs[2]);

	float lambda0 = (float)(mxGetScalar(prhs[3]));
	float lambda1 = (float)(mxGetScalar(prhs[4]));

	printf("lambda0=%f,lambda1=%f\n", lambda0, lambda1);

	int maxIterations = (int)mxGetScalar(prhs[5]);

	const mwSize *sizeImage = mxGetDimensions(prhs[0]);

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
	const size_t sizeY[2] = { 4 * nPx, 1 };

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

	float* Kx1 = new float[nPx];
	float* Kx2 = new float[nPx];
	float* Kx3 = new float[nPx];
	float* Kx4 = new float[nPx];

	float* Kty1 = new float[nPx];
	float* Kty2 = new float[nPx];
	float* Kty3 = new float[nPx];
	float* Kty4 = new float[nPx];
	float* Kty5 = new float[nPx];
	float* Kty6 = new float[nPx];

	float tau = 0.25f;
	float tau2 = 1.0f;

	float sigma = 0.33f;



	//residuals
	float p = 0;
	float d = 0;
	float err = 1.0;

	float teilerW = 1.0f / (1.0f + tau2*lambda1);

	#pragma omp parallel
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

			tauAquad[tmpIndex] = ux[tmpIndex] * ux[tmpIndex] + uy[tmpIndex] * uy[tmpIndex];

			teiler1[tmpIndex] = ux[tmpIndex] / tauAquad[tmpIndex];
			teiler2[tmpIndex] = uy[tmpIndex] / tauAquad[tmpIndex];

			tauAquad[tmpIndex] = tau*tauAquad[tmpIndex];

			tauUx[tmpIndex] = ux[tmpIndex] * tau;
			tauUy[tmpIndex] = uy[tmpIndex] * tau;

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
				v1[tmpIndex] = 0.0f;
				v2[tmpIndex] = 0.0f;
				v3[tmpIndex] = 0.0f;
				v4[tmpIndex] = 0.0f;
				v5[tmpIndex] = 0.0f;
				v6[tmpIndex] = 0.0f;
			}

			Kty1[tmpIndex] = 0.0f;
			Kty2[tmpIndex] = 0.0f;
			Kty3[tmpIndex] = 0.0f;
			Kty4[tmpIndex] = 0.0f;
			Kty5[tmpIndex] = 0.0f;
			Kty6[tmpIndex] = 0.0f;

			if (nrhs > 8)
			{
				y1[tmpIndex] = (float)inputY[tmpIndex];
				y2[tmpIndex] = (float)inputY[nPx + tmpIndex];
				y3[tmpIndex] = (float)inputY[2 * nPx + tmpIndex];
				y4[tmpIndex] = (float)inputY[3 * nPx + tmpIndex];
			}
			else
			{
				y1[tmpIndex] = 0.0f;
				y2[tmpIndex] = 0.0f;
				y3[tmpIndex] = 0.0f;
				y4[tmpIndex] = 0.0f;
			}

			Kx1[tmpIndex] = 0.0f;
			Kx2[tmpIndex] = 0.0f;
			Kx3[tmpIndex] = 0.0f;
			Kx4[tmpIndex] = 0.0f;
		}
	}

	int iterations = 0;
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

				float v1Old = v1[tmpIndex];
				float v2Old = v2[tmpIndex];
				float v3Old = v3[tmpIndex];
				float v4Old = v4[tmpIndex];
				float v5Old = v5[tmpIndex];
				float v6Old = v6[tmpIndex];

				//transpose equals -div
				Kty1[tmpIndex] = -(dxm(y1, sizeImage, i, j) + dym(y2, sizeImage, i, j));
				Kty2[tmpIndex] = -(dxm(y3, sizeImage, i, j) + dym(y4, sizeImage, i, j));
				Kty3[tmpIndex] = y1[tmpIndex];
				Kty4[tmpIndex] = y2[tmpIndex];
				Kty5[tmpIndex] = y3[tmpIndex];
				Kty6[tmpIndex] = y4[tmpIndex];

				float v1Hat = v1[tmpIndex] - tau*Kty1[tmpIndex];
				float v2Hat = v2[tmpIndex] - tau*Kty2[tmpIndex];

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

				v3[tmpIndex] = (v3[tmpIndex] + tau2*Kty3[tmpIndex]) * teilerW;
				v4[tmpIndex] = (v4[tmpIndex] + tau2*Kty4[tmpIndex]) * teilerW;
				v5[tmpIndex] = (v5[tmpIndex] + tau2*Kty5[tmpIndex]) * teilerW;
				v6[tmpIndex] = (v6[tmpIndex] + tau2*Kty6[tmpIndex]) * teilerW;

				if (iterations % 50 == 0)
				{
					//residuals
					p += myAbs((v1Old - v1[tmpIndex]) / tau - Kty1Old + Kty1[tmpIndex])
						+ myAbs((v2Old - v2[tmpIndex]) / tau - Kty2Old + Kty2[tmpIndex])
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

				float y1Old = y1[tmpIndex];
				float y2Old = y2[tmpIndex];
				float y3Old = y3[tmpIndex];
				float y4Old = y4[tmpIndex];

				Kx1[tmpIndex] = dxp(v1, sizeImage, i, j) - v3[tmpIndex];
				Kx2[tmpIndex] = dyp(v1, sizeImage, i, j) - v4[tmpIndex];
				Kx3[tmpIndex] = dxp(v2, sizeImage, i, j) - v5[tmpIndex];
				Kx4[tmpIndex] = dyp(v2, sizeImage, i, j) - v6[tmpIndex];

				float y1Tilde = y1[tmpIndex] + sigma*(2 * Kx1[tmpIndex] - Kx1Old);
				float y2Tilde = y2[tmpIndex] + sigma*(2 * Kx2[tmpIndex] - Kx2Old);
				float y3Tilde = y3[tmpIndex] + sigma*(2 * Kx3[tmpIndex] - Kx3Old);
				float y4Tilde = y4[tmpIndex] + sigma*(2 * Kx4[tmpIndex] - Kx4Old);

				float divisor1 = myMax(1.0f, sqrtf(y1Tilde*y1Tilde + y2Tilde*y2Tilde) / lambda0);
				float divisor2 = myMax(1.0f, sqrtf(y3Tilde*y3Tilde + y4Tilde*y4Tilde) / lambda0);

				y1[tmpIndex] = y1Tilde / divisor1;
				y2[tmpIndex] = y2Tilde / divisor1;
				y3[tmpIndex] = y3Tilde / divisor2;
				y4[tmpIndex] = y4Tilde / divisor2;

				if (iterations % 50 == 0)
				{
					d += myAbs((y1Old - y1[tmpIndex]) / sigma - Kx1Old + Kx1[tmpIndex]) +
						myAbs((y2Old - y2[tmpIndex]) / sigma - Kx2Old + Kx2[tmpIndex]) +
						myAbs((y3Old - y3[tmpIndex]) / sigma - Kx3Old + Kx3[tmpIndex]) +
						myAbs((y4Old - y4[tmpIndex]) / sigma - Kx4Old + Kx4[tmpIndex]);
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

	//write output
	#pragma omp parallel
	for (int j = 0; j < sizeImage[1]; ++j)
	{
		for (int i = 0; i < sizeImage[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeImage, i, j);

			YOut[tmpIndex + 0 * nPx] = (double)y1[tmpIndex];
			YOut[tmpIndex + 1 * nPx] = (double)y2[tmpIndex];
			YOut[tmpIndex + 2 * nPx] = (double)y3[tmpIndex];
			YOut[tmpIndex + 3 * nPx] = (double)y4[tmpIndex];

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

	delete[] Kx1;
	delete[] Kx2;
	delete[] Kx3;
	delete[] Kx4;

	delete[] Kty1;
	delete[] Kty2;
	delete[] Kty3;
	delete[] Kty4;
	delete[] Kty5;
	delete[] Kty6;
}