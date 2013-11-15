#include "mex.h"

/*
 * gridMetrics.c
 * Input: Two matrices describing the [x] and [y] coordinates of a grid
 * Outout: Five matrices--In Order Below
 * 		dx_dxsi and dy_dxsi on G-Flux faces
 * 		dx_deta and dy_deta on F-Flux faces
 *		Jacobians at cell centers 
 *
 * This is a MEX-file for MATLAB
*/

double Jacobian(double x_xsi_L, double x_xsi_R, double y_xsi_L, double y_xsi_R, double x_eta_T, double x_eta_B, double y_eta_T, double y_eta_B)
{
	double x_xsi = 0.5*(x_xsi_L + x_xsi_R);
	double y_xsi = 0.5*(y_xsi_L + y_xsi_R);
	double x_eta = 0.5*(x_eta_T + x_eta_B);
	double y_eta = 0.5*(y_eta_T + y_eta_B);

	double J = x_xsi*y_eta - x_eta*y_xsi;	/*Currently computes the INVERSE Jacobian*/
	return J;
} 

void metrics(double *xMat, double *yMat, double *dx_dxsi, double *dy_dxsi, double *dx_deta, double *dy_deta, double *J, mwSize m, mwSize n)
{
	mwSize i;   /*i indexes the row--xsi direction*/
	mwSize j;	/*j indexes the column--eta direction*/
    
	double xTopLeft,xTopRight,xBottomLeft,xBottomRight;		/*These variables hold the x & y locations of the four*/
	double yTopLeft,yTopRight,yBottomLeft,yBottomRight;		/*corners (nodes) of each cell under computation*/

	/*Complete first row of dx_deta and dy_deta*/
	xTopLeft = xMat[0];
	yTopLeft = yMat[0];

	for (j=0; j<(n-1); j++) {
	
		xTopRight = xMat[m*(j+1)];	/*Get next element in*/		
		yTopRight = yMat[m*(j+1)];	/*top row of grid vals*/

		dx_deta[(m)*j] = xTopRight - xTopLeft;
		dy_deta[(m)*j] = yTopRight - yTopLeft;	

		xTopLeft = xTopRight;		/*Set left node value to right node value*/	
		yTopLeft = yTopRight;		/*In preparation for next cell to the right*/	
	}

	/*Complete first column of dx_dxsi and dy_dxsi*/
	xTopLeft = xMat[0];
	yTopLeft = yMat[0];

	for (i=0; i<(m-1); i++) {
	
		xBottomLeft = xMat[i+1];	/*Get next element in*/		
		yBottomLeft = yMat[i+1];	/*left col of grid vals*/

		dx_dxsi[i] = xBottomLeft - xTopLeft;
		dy_dxsi[i] = yBottomLeft - yTopLeft;

		xTopLeft = xBottomLeft;		/*Set top node value to bottom node value*/
		yTopLeft = yBottomLeft;		/*in preparation for next cell below the current*/
	}


	/*Loop through cells, calculate right xsi metrics and bottom eta metrics*/
	/*and then compute Jacobians at cell centers--Note that there are (M-1)x(N-1) cells*/
	double x_xsi_left, x_xsi_right, y_xsi_left, y_xsi_right;		/*These variables hold the xsi metrics on the left & right faces and*/
	double x_eta_top, x_eta_bottom, y_eta_top, y_eta_bottom;		/*the eta metrics on the top and bottom faces of each cell under computation*/

	for (j=0; j<(n-1); j++) {

		x_eta_top = dx_deta[j*(m)];		/*For each new column, set top metrics equal to*/
		y_eta_top = dy_deta[j*(m)];		/*the first element of the d?_deta[] array column*/
	
		for (i=0; i<(m-1); i++) {

			x_xsi_left = dx_dxsi[i + j*(m-1)];	/*Set xsi left metric for left face of cell*/
			y_xsi_left = dy_dxsi[i + j*(m-1)];

			xTopRight = xMat[i + (j+1)*m];				/*Top right node of cell (i,j)*/
			yTopRight = yMat[i + (j+1)*m];				/*Top right node of cell (i,j)*/
 
			xBottomRight = xMat[(i+1) + (j+1)*m];		/*Bottom right node of cell (i,j)*/
			yBottomRight = yMat[(i+1) + (j+1)*m];		/*Bottom right node of cell (i,j)*/

			xBottomLeft = xMat[(i+1) + j*m];			/*Bottom left node of cell (i,j)*/
			yBottomLeft = yMat[(i+1) + j*m];			/*Bottom left node of cell (i,j)*/
 
 			x_xsi_right = xBottomRight - xTopRight;		/*Compute new metric at right face*/
 			y_xsi_right = yBottomRight - yTopRight;		/*Compute new metric at right face*/

			x_eta_bottom = xBottomRight - xBottomLeft;	/*Compute new metric at bottom face*/
			y_eta_bottom = yBottomRight - yBottomLeft;	/*Compute new metric at bottom face*/

			J[i + j*(m-1)] = Jacobian(x_xsi_left, x_xsi_right, y_xsi_left, y_xsi_right, x_eta_top, x_eta_bottom, y_eta_top, y_eta_bottom);	/*Compute Jacobian at cell center*/

			dx_dxsi[i + (j+1)*(m-1)] = x_xsi_right;		/*Write right face metrics to xsi array*/	
			dy_dxsi[i + (j+1)*(m-1)] = y_xsi_right;		/*Write right face metrics to xsi array*/

			dx_deta[(i+1) + j*(m)] = x_eta_bottom;		/*Write bottom face metrics to eta array*/	
			dy_deta[(i+1) + j*(m)] = y_eta_bottom;		/*Write bottom face metrics to eta array*/

			x_eta_top = x_eta_bottom;					/*Set top metrics to bottom metrics to prepare for next lower cell ((i+1),j)*/	
			y_eta_top = y_eta_bottom;					/*Set top metrics to bottom metrics to prepare for next lower cell ((i+1),j)*/	

		}
	}   
}

/* Gateway Function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* variable declarations here */
	double *xMatrix;		/* MxN input matrix of x-values */
	double *yMatrix;		/* MxN input matrix of y-values */
	mwSize ncols;			/* size of grid */
	mwSize nrows;			/* size of grid */

	/* create pointers to the real data in the input matrices */
	xMatrix = mxGetPr(prhs[0]);
	yMatrix = mxGetPr(prhs[1]);

	/* get dimensions of the input matrix */
	ncols = mxGetN(prhs[0]);
	nrows = mxGetM(prhs[0]);

	double *dx_dxsi;		/* (M-1)xN output matrix  */
	double *dy_dxsi;		/* (M-1)xN output matrix  */
	double *dx_deta;		/* Mx(N-1) output matrix  */
	double *dy_deta;		/* Mx(N-1) output matrix  */
	double *J;				/* (M-1)x(N-1) output matrix */

	/* create the output matrices */
	plhs[0] = mxCreateDoubleMatrix((nrows-1),ncols,mxREAL);	/*dx_dxsi: (M-1)xN */
	plhs[1] = mxCreateDoubleMatrix((nrows-1),ncols,mxREAL);	/*dy_dxsi: (M-1)xN */
	plhs[2] = mxCreateDoubleMatrix(nrows,(ncols-1),mxREAL);	/*dx_deta: Mx(N-1) */
	plhs[3] = mxCreateDoubleMatrix(nrows,(ncols-1),mxREAL);	/*dy_deta: Mx(N-1) */
	plhs[4] = mxCreateDoubleMatrix((nrows-1),(ncols-1),mxREAL);		/*J: (M-1)x(N-1) */

	/* get a pointer to the real data in the output matrix */
	dx_dxsi = mxGetPr(plhs[0]);
	dy_dxsi = mxGetPr(plhs[1]);
	dx_deta = mxGetPr(plhs[2]);
	dy_deta = mxGetPr(plhs[3]);
	J = mxGetPr(plhs[4]);
	
	/* code here */
	/* check for proper number of arguments */
	if (nlhs!=5) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "Five outputs required.");
	}
	if (nrhs!=2) {	
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Two inputs required.");
	}

	/* call computational routine */
	metrics(xMatrix, yMatrix, dx_dxsi, dy_dxsi, dx_deta, dy_deta, J, nrows, ncols);


}
