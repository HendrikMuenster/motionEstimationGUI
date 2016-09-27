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
//#include "tools.h"

#include <vector>

int index2DtoLinear(const mwSize *sizeMat, int i, int j);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//L1L2OpticalFlow(u1,u2,tol,lambda,ux,uy,ut)
{
	double *dimsImagePr = mxGetPr(prhs[0]);
    int dimsImage[2] = { (int)dimsImagePr[0],(int)dimsImagePr[1]};
    
    int filterSize = (int)mxGetScalar(prhs[1]);
    
    int filterSizeH = (filterSize-1)/2;
    
    std::vector<double> listI,listJ,listVal;
    
    int rowNumber = 1;
    
    mexPrintf("Dims: [%d,%d]\n",dimsImage[0],dimsImage[1]);
    
    for (int j = 1; j <= dimsImage[1]; ++j)
	{
		for (int i = 1; i <= dimsImage[0]; ++i)
		{
            int currentPoint = index2DtoLinear(dimsImage,i,j);
            //go through all neighbours
            for (int nx = -filterSizeH; nx<=filterSizeH; ++nx)
            {
                for (int ny = -filterSizeH; ny<=filterSizeH; ++ny)
                {
                    if ((i+nx) >= 1 && (j+ny) >= 1 && (i+nx) <= dimsImage[0] && (j+ny) <= dimsImage[1])
                    {
                        // no difference of point with himself
                        if (nx != 0 || ny != 0)
                        {
                            listI.push_back((double)rowNumber);
                            listJ.push_back((double)currentPoint);
                            listVal.push_back(1.0);

                            listI.push_back((double)rowNumber);
                            listJ.push_back((double)index2DtoLinear(dimsImage,i+nx,j+ny));
                            listVal.push_back(-1.0);

                            rowNumber += 1;
                        }
                    }
                }
            }
            
        }
    }
    
    const int result[2] = { listI.size(), 1 };
    
    plhs[0] = mxCreateNumericArray(1,result , mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(1,result , mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(1,result , mxDOUBLE_CLASS, mxREAL);
	
    double *outListI = mxGetPr(plhs[0]);
    double *outListJ = mxGetPr(plhs[1]);
    double *outListVal = mxGetPr(plhs[2]);
    
    
    for (int i = 0; i < listI.size(); ++i)
	{
        outListI[i] = listI[i];
        outListJ[i] = listJ[i];
        outListVal[i] = listVal[i];
	}
    

}

int index2DtoLinear(const mwSize *sizeMat, int i, int j)
{
	return (int)(i + (j-1)*sizeMat[0]);
}

