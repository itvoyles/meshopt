#include "mex.h"
#include <math.h>

/*
 * truncErrorU1.c
 * This function will return the total residual (TE) on the given mesh for the U1 function (continuity equation)
 * Input: Two matrices describing the [x] and [y] coordinates of a grid plus a "parameters" vector holding solution details
 * Output: Nine matrices--In Order Below
 * 		dx_dxsi and dy_dxsi on G-Flux faces
 * 		dx_deta and dy_deta on F-Flux faces
 *		Jacobians at cell centers
 *		exact1,exact2,exact3,exact4--holding the exact solutions to 4 equations
 *		
 *		All parameters are actually passed through RHS
 *		exactSolver(x,y,x_xsi,y_xsi,x_eta,y_eta,J,exact1,exact2,exact3,exact4,parameters); 
 *
 * This is a MEX-file for MATLAB
*/

__inline double Jacobian(double x_xsi_L, double x_xsi_R, double y_xsi_L, double y_xsi_R, double x_eta_T, double x_eta_B, double y_eta_T, double y_eta_B)
{
	double x_xsi = 0.5*(x_xsi_L + x_xsi_R);
	double y_xsi = 0.5*(y_xsi_L + y_xsi_R);
	double x_eta = 0.5*(x_eta_T + x_eta_B);
	double y_eta = 0.5*(y_eta_T + y_eta_B);

	double J = 1.0 / (x_xsi*y_eta - x_eta*y_xsi);	/*Jacobian*/
	return J;
} 

__inline double calcExact(double coeff, double param1, double param2)
{
	double sol = (param2 - param1)*coeff + param1;

	return sol;
}

__inline double primToU2(double rho, double uvel)
{
	double u2 = rho*uvel;
	return u2;
}

__inline double primToU3(double rho, double vvel)
{
	double u3 = rho*vvel;
	return u3;
}

__inline double primToU4(double rho, double uvel, double vvel, double Prs)
{
	/*Caclulate rho*et*/
	double u4 = 1.0/(1.4-1.0)*Prs + 0.5*rho*(uvel*uvel+vvel*vvel);
	return u4;
}

__inline double primToht(double rho, double uvel, double vvel, double Prs)
{
	/*Caclulate ht*/
	double ht = 1.4/(1.4-1.0)*(Prs/rho) + 0.5*(uvel*uvel+vvel*vvel);
	return ht;
}

__inline double modEval(double lambda, double e, double Roe_a)
{
	/*Keeps eigenvalues from going to zero and avoids an expansion shock in Roe method*/
	double lam_check = fabs(lambda) - 2.0*e*Roe_a;
	double lam_limit = lambda*lambda / (4.0*e*Roe_a) + e*Roe_a;
	double lam_sign = copysign(lambda, lam_check);
	double lam_smoothP = 0.5*(lam_sign + fabs(lam_sign));
	lam_sign = copysign(lam_limit, lam_check);
	double lam_smoothN = 0.5*(lam_sign - fabs(lam_sign));
	lambda = lam_smoothP - lam_smoothN;

	return lambda;
}

void PM_function(double M1, double V1, double angle_diff, double M2, double *local_M2, double *local_turn)
{
	/*Find local M2 and turn angle corresponding to a fan with the width of angle_diff.   This requires newton iteration*/

	double V2_1, theta2_1, u2_1, f_1, dfdx;

	double epsilon = 3.0*pow(10.0,-15.0);
	double u1 = asin(1.0/M1);		//calculate upstream mach angle u1

	//Guess second mach number M2--assuming M1 < (local)M2 < M2
	//initial guesses
	double M2_0 = 0.5*(M2+M1);
	double M2_1 = 0.5*(M2+M1)+fabs(M2-M1)*0.01;

	//get v(M2)
	double V2_0 = sqrt( (1.4+1.0)/(1.4-1.0) )*atan( sqrt( (1.4-1.0)*(M2_0*M2_0-1.0)/(1.4+1.0) ) ) - atan( sqrt( M2_0*M2_0-1.0 ) );

	//get theta2=v2-v1
	double theta2_0 = V2_0 - V1;

	//get second mach angle
	double u2_0 = asin(1.0/M2_0);

	//(theta2 - u2) - (fan_width - u1) = 0 {function to be solved}
	double f_0 = (theta2_0 - u2_0) - (angle_diff - u1);

	while( fabs(M2_1-M2_0) > epsilon ){

		V2_1 = sqrt( (1.4+1.0)/(1.4-1.0) )*atan( sqrt( (1.4-1.0)*(M2_1*M2_1-1.0)/(1.4+1.0) ) ) - atan( sqrt( M2_1*M2_1-1.0 ) );
		theta2_1 = V2_1 - V1;
		u2_1 = asin(1.0/M2_1);

		f_1 = (theta2_1 - u2_1) - (angle_diff - u1);

		dfdx = (f_1 - f_0)/(M2_1 - M2_0);
    
    	f_0 = f_1;
    	M2_0 = M2_1;
    
    	M2_1 = M2_1 - f_1/dfdx;
	}
	
	V2_1 = sqrt( (1.4+1.0)/(1.4-1.0) )*atan( sqrt( (1.4-1.0)*(M2_1*M2_1-1.0)/(1.4+1.0) ) ) - atan( sqrt( M2_1*M2_1-1.0 ) );
	theta2_1 = V2_1 - V1;

	*local_M2 = M2_1;
	*local_turn = theta2_1;
}



void metrics(double *xMat, double *yMat, double *dx_dxsi, double *dy_dxsi, double *dx_deta, double *dy_deta, double *J, mwSize m, mwSize n)
{
	mwSize i;   /*i indexes the row--xsi direction*/
	mwSize j;	/*j indexes the column--eta direction*/

	mwSize index;  /*simplified index for optimizing--will be set to i+j*m for 2D array*/
    
	double xTopLeft,xTopRight,xBottomLeft,xBottomRight;		/*These variables hold the x & y locations of the four*/
	double yTopLeft,yTopRight,yBottomLeft,yBottomRight;		/*corners (nodes) of each cell under computation*/

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
		
		xTopLeft = xMat[m*(j)];	/*Get next element in*/		
		yTopLeft = yMat[m*(j)];	/*top row of grid vals*/

		xTopRight = xMat[m*(j+1)];	/*Get next element in*/		
		yTopRight = yMat[m*(j+1)];	/*top row of grid vals*/

		x_eta_top = xTopRight - xTopLeft;		/*For each new column, find top metrics*/
		y_eta_top = yTopRight - yTopLeft;

		dx_deta[(m)*j] = x_eta_top;
		dy_deta[(m)*j] = y_eta_top;	
	
		for (i=0; i<(m-1); i++) {

			index = i + j*m;

			x_xsi_left = dx_dxsi[index - j];	/*Set xsi left metric for left face of cell--[index - j]=[i + j(m-1)]*/
			y_xsi_left = dy_dxsi[index - j];

			xTopRight = xMat[index + m];				/*Top right node of cell (i,j)--[index + m]=[i + (j+1)m]*/
			yTopRight = yMat[index + m];				/*Top right node of cell (i,j)*/
 
			xBottomRight = xMat[index + m + 1];		/*Bottom right node of cell (i,j)--[index + m + 1]=[(i+1) + (j+1)m]*/
			yBottomRight = yMat[index + m + 1];		/*Bottom right node of cell (i,j)*/

			xBottomLeft = xMat[index + 1];			/*Bottom left node of cell (i,j)--[index + 1]=[(i+1) + j*m]*/
			yBottomLeft = yMat[index + 1];			/*Bottom left node of cell (i,j)*/

			x_xsi_right = xBottomRight - xTopRight;		/*Compute new metric at right face*/
 			y_xsi_right = yBottomRight - yTopRight;		/*Compute new metric at right face*/

			x_eta_bottom = xBottomRight - xBottomLeft;	/*Compute new metric at bottom face*/
			y_eta_bottom = yBottomRight - yBottomLeft;	/*Compute new metric at bottom face*/

			J[index - j] = Jacobian(x_xsi_left, x_xsi_right, y_xsi_left, y_xsi_right, x_eta_top, x_eta_bottom, y_eta_top, y_eta_bottom);	/*Compute Jacobian at cell center--[index - j]=[i + j(m-1)]*/

			dx_dxsi[index - j + m - 1] = x_xsi_right;		/*Write right face metrics to xsi array--[index - j + m - 1]=[i + (j+1)(m-1)]*/	
			dy_dxsi[index - j + m - 1] = y_xsi_right;		/*Write right face metrics to xsi array*/

			dx_deta[index + 1] = x_eta_bottom;		/*Write bottom face metrics to eta array--[index + 1]=[(i+1) + j*m]*/	
			dy_deta[index + 1] = y_eta_bottom;		/*Write bottom face metrics to eta array*/

			x_eta_top = x_eta_bottom;					/*Set top metrics to bottom metrics to prepare for next lower cell ((i+1),j)*/	
			y_eta_top = y_eta_bottom;					/*Set top metrics to bottom metrics to prepare for next lower cell ((i+1),j)*/	

		}
	}   
}

void setExact(double *xMat, double *yMat, double *exact1, double *exact2, double *exact3, double *exact4, mwSize m, mwSize n, double *param)
{
	mwSize i;   /*i indexes the row--xsi direction*/
	mwSize j;	/*j indexes the column--eta direction*/

	mwSize index;  /*simplified index for optimizing--will be set to i+j*m for 2D array*/

	double x_node, y_node;  /*The x,y position of the node at which the exact solution is calculated*/
	double solution1, solution2, solution3, solution4;  /*These temporarily will hold the four exact solution values at each node*/

	double yDomainSize = fabs( yMat[0 + (n-1)*m] - yMat[0] );   /*Get length of domain in y direction*/	
	double xDomainSize = fabs( xMat[m-1] - xMat[0] );			/*Get length of domain in x direction*/

	/*Get the parameters needed for exact solution from the param vector*/
	double rho1 = param[0];  /*Upstream conditions*/
	double uvel1 = param[1];
	double vvel1 = param[2];
	double P1 = param[3];

	double rho2 = param[4]; /*Downstream condtions*/
	double uvel2 = param[5];
	double vvel2 = param[6];
	double P2 = param[7];

	double u1 = param[8];   /*Forward mach line angle*/
	double u2 = param[9];	/*Rear angle*/
	double v1 = param[12];	/*Prandtl meyer function v(M1) for the upstream conditions*/

	double M1 = ( sqrt(uvel1*uvel1+vvel1*vvel1) ) / sqrt(1.4*P1/rho1);	/*Upstream mach number*/
	double M2 = ( sqrt(uvel2*uvel2+vvel2*vvel2) ) / sqrt(1.4*P2/rho2);	/*Downstream mach number*/
	double T1 = P1 / (rho1*287.04);		/*Upstream temperature*/

	double local_M2;	/*The mach number within the fan */
	double local_turn;	/*Turn angle for use inside the fan, as if the fan terminated at the local angle/current M2 mach number*/
	double local_T2;	/*Local temp in fan*/
	double local_a;		/*Local sound speed*/
	double local_speed;	/*Local flow speed (magnitude of veloc)*/

	double angle;           /*Will be used to hold current node angle*/
	double angle_diff;	/*Angle between u1 and the current node*/
	
	double global_turn = atan(vvel1/uvel1);		/*true turn angle of entire flow*/

	bool inFan;		/*Indicates whether current node is inside the fan*/

	/*Note that the exact solution matrices include a row of ghost cells around the entire domain.*/
	/*This means that the exact solution size is (M+1)*(N+1) for (M-1)*(N-1) total cells.   There are*/
	/*There are M*N nodes (x,y points) on the grid, or (M+2)*(N+2) including "ghost" nodes.*/
	
	/*Compute interior columns of exact solution.   Add 1/4*ExactSolution to each of the surrounding cells */	
   	for (j=0; j<n; j++) {
		for (i=0; i<m; i++) {
			
			index = i + j*m;

			inFan = false;

			x_node = xMat[index];   /*Current node--[index]=[i + j*m]*/
			y_node = yMat[index];

			/*Calculate current angle of vector from origin to node, relative to bottom of domain--note atan2(0,0)=0 at origin*/
			angle = atan2(y_node, x_node);

			/*If at origin, angle is half way between the two mach lines u1 and u2*/
			if( (fabs(y_node)+fabs(x_node))==0.0)
				angle = 0.5*(u1+u2);
			
			if(angle>=u1){			//if ahead of the fan
				solution1 = rho1;
				solution2 = uvel1;
				solution3 = vvel1;
				solution4 = P1;
			}
			else if(angle<=u2){		//if behind the fan
				solution1 = rho2;
				solution2 = uvel2;
				solution3 = vvel2;
				solution4 = P2;
			}
			else{				//if within the fan
				inFan = true;

				/*Find the angle between the lead mach line u1 and the current angle*/
				angle_diff = u1 - angle;

				/*Find local M2 and turn angle corresponding to a fan with the width of angle_diff.   This requires iteration*/
				PM_function(M1, v1, angle_diff, M2, &local_M2, &local_turn);	

				/*Calculate pressure (solution4)*/
				solution4 = P1/pow(( (1.0+(1.4-1.0)*local_M2*local_M2/2.0)/(1.0+(1.4-1.0)*M1*M1/2.0) ),(1.4/(1.4-1.0)));
			}
			
			solution4 = 0.25*solution4;	/*Pressure*/

			/*Write 1/4 of exact solution #1 to four surrounding cells*/
			exact4[index + j] += solution4;
			exact4[index + j + 1] += solution4;
			exact4[index + j + m + 1] += solution4;
			exact4[index + j + m + 2] += solution4;

			solution4 = 4.0*solution4;

			if(inFan){
				/*Compute local temperature*/
				local_T2 = T1*( (1.0+(1.4-1.0)*M1*M1/2.0)/(1.0+(1.4-1.0)*local_M2*local_M2/2.0) );
				/*Compute local density (solution1)*/
				solution1 = solution4/(287.04*local_T2);
			}

			solution1 = 0.25*solution1;	/*Density*/

			/*Write 1/4 of exact solution #1 to four surrounding cells*/
			exact1[index + j] += solution1;             /*[index + j]=[i + j*(m+1)]*/
			exact1[index + j + 1] += solution1;			/*[index + j + 1]=[(i+1) + j*(m+1)]*/
			exact1[index + j + m + 1] += solution1;		/*[index + j + m + 1]=[i + (j+1)*(m+1)]*/
			exact1[index + j + m + 2] += solution1;		/*[index + j + m + 2]=[(i+1) + (j+1)*(m+1)]*/

			solution1 = 4.0*solution1; 

			if(inFan){
				//use local_turn angle, local_M2, and local_a to get the velocity components
				local_a = sqrt(1.4*solution4/solution1);
				local_speed = local_a*local_M2;
				
				//Note that local_turn is now being used to store the difference between the local flow direction and the 
				//bottom of the domain.   Thus local_turn is now the angle between the flow and the x-axis of domain
				local_turn = fabs(global_turn - local_turn);
	
				//veloc components
				solution2 = cos(local_turn)*local_speed;
				solution3 = sin(local_turn)*local_speed;
			}
			
			solution2 = 0.25*solution2;

			/*Write 1/4 of exact solution #2 to four surrounding cells*/
			exact2[index + j] += solution2;
			exact2[index + j + 1] += solution2;
			exact2[index + j + m + 1] += solution2;
			exact2[index + j + m + 2] += solution2; 


			solution3 = 0.25*solution3;

			/*Write 1/4 of exact solution #3 to four surrounding cells*/
			exact3[index + j] += solution3;
			exact3[index + j + 1] += solution3;
			exact3[index + j + m + 1] += solution3;
			exact3[index + j + m + 2] += solution3;

		}
	}
	
	/*Compute first and last columns of exact solution*/
	for (i=0; i<m; i++) {
		
		/*First Column of exact solution--j=0*/

		exact1[i] = rho2;
		exact1[i+1] = rho2;

		exact2[i] = uvel2;
		exact2[i+1] = uvel2;

		exact3[i] = vvel2;
		exact3[i+1] = vvel2;

		exact4[i] = P2;
		exact4[i+1] = P2;

		/*Last Column of exact solution--j=n < n+1*/

		inFan = false;		

		x_node = xMat[i + (n-1)*m];
		y_node = yMat[i + (n-1)*m] + yDomainSize/(double)n;

		/*Calculate current angle of vector from origin to node, relative to bottom of domain*/
		angle = atan2(y_node, x_node);

		if(angle>=u1){			//if ahead of the fan
			solution1 = rho1;
			solution2 = uvel1;
			solution3 = vvel1;
			solution4 = P1;
		}
		else if(angle<=u2){		//if behind the fan
			solution1 = rho2;
			solution2 = uvel2;
			solution3 = vvel2;
			solution4 = P2;
		}
		else{				//if within the fan
			inFan = true;

			/*Find the angle between the lead mach line u1 and the current angle*/
			angle_diff = u1 - angle;

			/*Find local M2 and turn angle corresponding to a fan with the width of angle_diff.   This requires iteration*/
			PM_function(M1, v1, angle_diff, M2, &local_M2, &local_turn);	

			/*Calculate pressure (solution4)*/
			solution4 = P1*pow(( (1.0+(1.4-1.0)*local_M2*local_M2/2.0)/(1.0+(1.4-1.0)*M1*M1/2.0) ),(-1.4/(1.4-1.0)));
		}

		solution4 = 0.25*solution4;

		exact4[i + (n)*(m+1)] += solution4;
		exact4[i+1 + (n)*(m+1)] += solution4;

		solution4 = 4.0*solution4;

		if(inFan){
			/*Compute local temperature*/
			local_T2 = T1*( (1.0+(1.4-1.0)*M1*M1/2.0)/(1.0+(1.4-1.0)*local_M2*local_M2/2.0) );
			/*Compute local density (solution1)*/
			solution1 = solution4/(287.04*local_T2);
		}

		solution1 = 0.25*solution1;

		exact1[i + (n)*(m+1)] += solution1;
		exact1[i+1 + (n)*(m+1)] += solution1;

		solution1 = 4.0*solution1;

		if(inFan){
			//use local_turn angle, local_M2, and local_a to get the velocity components
			local_a = sqrt(1.4*solution4/solution1);
			local_speed = local_a*local_M2;
			
			//Note that local_turn is now being used to store the difference between the local flow direction and the 
			//bottom of the domain.   Thus local_turn is now the angle between the flow and the x-axis of domain
			local_turn = fabs(global_turn - local_turn);
	
			//veloc components
			solution2 = cos(local_turn)*local_speed;
			solution3 = sin(local_turn)*local_speed;
		}

		solution2 = 0.25*solution2;

		exact2[i + (n)*(m+1)] += solution2;
		exact2[i+1 + (n)*(m+1)] += solution2;

		solution3 = 0.25*solution3;

		exact3[i + (n)*(m+1)] += solution3;
		exact3[i+1 + (n)*(m+1)] += solution3;


		

	}

	/*Calculate top and bottom rows of exact solution*/
	for (j=0; j<n; j++) {
		
		/*Top row of exact solution--i=0*/

		exact1[j*(m+1)] = rho1;
		exact1[(j+1)*(m+1)] = rho1;

		exact2[j*(m+1)] = uvel1;
		exact2[(j+1)*(m+1)] = uvel1;

		exact3[j*(m+1)] = vvel1;
		exact3[(j+1)*(m+1)] = vvel1;

		exact4[j*(m+1)] = P1;
		exact4[(j+1)*(m+1)] = P1;

		/*Bottom Row of exact sol--i=m < m+1*/
	
		inFan = false;

		x_node = xMat[(m-1) + j*m] + xDomainSize/(double)m;
		y_node = yMat[(m-1) + j*m];

		/*Calculate current angle of vector from origin to node, relative to bottom of domain*/
		angle = atan2(y_node, x_node);

		if(angle>=u1){			//if ahead of the fan
			solution1 = rho1;
			solution2 = uvel1;
			solution3 = vvel1;
			solution4 = P1;
		}
		else if(angle<=u2){		//if behind the fan
			solution1 = rho2;
			solution2 = uvel2;
			solution3 = vvel2;
			solution4 = P2;
		}
		else{				//if within the fan
			inFan = true;

			/*Find the angle between the lead mach line u1 and the current angle*/
			angle_diff = u1 - angle;

			/*Find local M2 and turn angle corresponding to a fan with the width of angle_diff.   This requires iteration*/
			PM_function(M1, v1, angle_diff, M2, &local_M2, &local_turn);	

			/*Calculate pressure (solution4)*/
			solution4 = P1*pow(( (1.0+(1.4-1.0)*local_M2*local_M2/2.0)/(1.0+(1.4-1.0)*M1*M1/2.0) ),(-1.4/(1.4-1.0)));
		}

		solution4 = 0.25*solution4;

		exact4[m + j*(m+1)] += solution4;
		exact4[m + (j+1)*(m+1)] += solution4;

		solution4 = 4.0*solution4;
		
		if(inFan){
			/*Compute local temperature*/
			local_T2 = T1*( (1.0+(1.4-1.0)*M1*M1/2.0)/(1.0+(1.4-1.0)*local_M2*local_M2/2.0) );
			/*Compute local density (solution1)*/
			solution1 = solution4/(287.04*local_T2);
		}

		solution1 = 0.25*solution1;

		exact1[m + j*(m+1)] += solution1;
		exact1[m + (j+1)*(m+1)] += solution1;

		solution1 = 4.0*solution1;

		if(inFan){
			//use local_turn angle, local_M2, and local_a to get the velocity components
			local_a = sqrt(1.4*solution4/solution1);
			local_speed = local_a*local_M2;
			
			//Note that local_turn is now being used to store the difference between the local flow direction and the 
			//bottom of the domain.   Thus local_turn is now the angle between the flow and the x-axis of domain
			local_turn = fabs(global_turn - local_turn);
	
			//veloc components
			solution2 = cos(local_turn)*local_speed;
			solution3 = sin(local_turn)*local_speed;
		}

		solution2 = 0.25*solution2;

		exact2[m + j*(m+1)] += solution2;
		exact2[m + (j+1)*(m+1)] += solution2;

		solution3 = 0.25*solution3;

		exact3[m + j*(m+1)] += solution3;
		exact3[m + (j+1)*(m+1)] += solution3;

	}

}



void interiorFlux(double *fluxF, double *fluxG, double *dx_dxsi, double *dy_dxsi, double *dx_deta, double *dy_deta, double *J, double *ex1, double *ex2, double *ex3, double *ex4, mwSize m, mwSize n, double *param)
{
	mwSize i;
	mwSize j;
	int index;

	double epsilon = param[10];   /*Epsilon and kappa--determines type of MUSCL extrapolation*/
	double kappa = param[11];
	double gamma = 1.4;

	/*L&R primitive state variables*/
	double rhoL, rhoR, uvelL, uvelR, vvelL, vvelR, PrsL, PrsR, htR, htL;
	/*L&R conserved variables*/
	double U1_L, U2_L, U3_L, U4_L, U1_R, U2_R, U3_R, U4_R;
	/*Roe avg variables*/
	double R, Roe_rho, Roe_u, Roe_v, Roe_ht, Roe_a;
	/*Grid metrics*/
	double x_eta, y_eta, x_xsi, y_xsi, J_F, J_G;
	double xsi_x, xsi_y, eta_x, eta_y;
	double xsi_x_hat, xsi_y_hat, eta_x_hat, eta_y_hat;
	/*Contravarient velocity*/
	double contra_UL, contra_UR, contra_VL, contra_VR;
	/*Eigenvalues*/
	double lambda1, lambda2, lambda3, lambda4;
	/*Eigenvectors (only the first element of each of 4 vectors since this function calculates TE for eq 1)*/
	double r1, r3, r4;
	/*Variables for wave amplitudes*/
	double dU1, dU2, dU3, dU4;
	double dw1, dw2, dw3, dw4;
	/*Left and right fluxes--only need first flux equation*/
	double F1L, G1L;
	double F1R, G1R;
	/*Interface flux*/
	double F1, G1;
	/*Small parameter for avoiding expansion shock*/
	double e = 0.001;

	/*Calculate first column of F-Fluxes*/
	j = 0;
	for(i=1; i<(m-1); i++){
		
		/*MUSCL extrapolate left and right states for F-flux*/
		/*The current cell index (the cell to the left of the interface being calculated) is calculated below*/
		index = i + (j+1)*(m+1);  /*Note: there are m+1 rows of cells (including ghost).   First column is a ghost col. so it is skipped*/

		rhoL = ex1[index] + 0.25*epsilon*( (1.0-kappa)*(ex1[index]-ex1[index-1]) + (1.0+kappa)*(ex1[index+1]-ex1[index]) );
		uvelL = ex2[index] + 0.25*epsilon*( (1.0-kappa)*(ex2[index]-ex2[index-1]) + (1.0+kappa)*(ex2[index+1]-ex2[index]) );
		vvelL = ex3[index] + 0.25*epsilon*( (1.0-kappa)*(ex3[index]-ex3[index-1]) + (1.0+kappa)*(ex3[index+1]-ex3[index]) );
		PrsL = ex4[index] + 0.25*epsilon*( (1.0-kappa)*(ex4[index]-ex4[index-1]) + (1.0+kappa)*(ex4[index+1]-ex4[index]) );
 	
		rhoR = ex1[index+1] - 0.25*epsilon*( (1.0+kappa)*(ex1[index+1]-ex1[index]) + (1.0-kappa)*(ex1[index+2]-ex1[index+1]) );
		uvelR = ex2[index+1] - 0.25*epsilon*( (1.0+kappa)*(ex2[index+1]-ex2[index]) + (1.0-kappa)*(ex2[index+2]-ex2[index+1]) );
		vvelR = ex3[index+1] - 0.25*epsilon*( (1.0+kappa)*(ex3[index+1]-ex3[index]) + (1.0-kappa)*(ex3[index+2]-ex3[index+1]) );
		PrsR = ex4[index+1] - 0.25*epsilon*( (1.0+kappa)*(ex4[index+1]-ex4[index]) + (1.0-kappa)*(ex4[index+2]-ex4[index+1]) );

		/*Calculate the left and right total enthalpy and conserved variables*/
		U1_L = rhoL;
		U2_L = primToU2(rhoL, uvelL);
		U3_L = primToU3(rhoL, vvelL);
		U4_L = primToU4(rhoL, uvelL, vvelL, PrsL);
		htL = primToht(rhoL, uvelL, vvelL, PrsL);

		U1_R = rhoR;
		U2_R = primToU2(rhoR, uvelR);
		U3_R = primToU3(rhoR, vvelR);
		U4_R = primToU4(rhoR, uvelR, vvelR, PrsR);
		htR = primToht(rhoR, uvelR, vvelR, PrsR);

		/*Roe avg variables*/
		R = sqrt(rhoR/rhoL);
		Roe_rho = R*rhoL;
		Roe_u = (R*uvelR + uvelL)/(R + 1.0);
		Roe_v = (R*vvelR + vvelL)/(R + 1.0);
		Roe_ht = (R*htR + htL)/(R + 1.0);
		Roe_a = sqrt( (gamma-1.0)*fabs(Roe_ht - 0.5*(Roe_u*Roe_u + Roe_v*Roe_v) ) );

		/*Get gridmetrics for F-flux interface*/
		x_eta = dx_deta[i];
		y_eta = dy_deta[i];

		J_F = 0.5*(J[i-1] + J[i]);   /*Interpolate J to current F interface*/

		xsi_x = y_eta*J_F;			/*Calculate inverse transform metrics at current F interface*/
		xsi_y = -1.0*x_eta*J_F;

		/*Eigenvalues, F flux*/
		lambda1 = xsi_x*Roe_u + xsi_y*Roe_v;
			/*lambda 2, same as 1 but un-needed for this function*/
		lambda3 = lambda1 + Roe_a*sqrt(xsi_x*xsi_x + xsi_y*xsi_y);
		lambda4 = lambda1 - Roe_a*sqrt(xsi_x*xsi_x + xsi_y*xsi_y);

		/*Right eigenvectors, F-flux*/
		r1 = 1.0;
		r3 = Roe_rho / (2.0*Roe_a);
		r4 = r3;

		/*Wave amplitudes*/
		dU1 = U1_R - U1_L;
		dU2 = U2_R - U2_L;
		dU3 = U3_R - U3_L;
		dU4 = U4_R - U4_L;

		xsi_x_hat = xsi_x / sqrt(xsi_x*xsi_x + xsi_y*xsi_y); 
		xsi_y_hat = xsi_y / sqrt(xsi_x*xsi_x + xsi_y*xsi_y);

		dw1 = 2.0*(dU1 + (dU4 + dU1*Roe_ht - dU2*Roe_u - dU3*Roe_v)/(-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v) );
		/*dw2 = (-1.0*(dU3*xsi_x_hat) + dU1*Roe_v*xsi_x_hat + dU2*xsi_y_hat - dU1*Roe_u*xsi_y_hat) / (Roe_rho*(xsi_x_hat*xsi_x_hat + xsi_y_hat*xsi_y_hat  ) );*/
		dw3 = (-1.0*((Roe_a*(2.0*dU4 - 2.0*dU2*Roe_u - 2.0*dU3*Roe_v + dU1*(Roe_u*Roe_u + Roe_v*Roe_v)))/(-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v)) + (dU2*xsi_x_hat + dU3*xsi_y_hat - dU1*(Roe_u*xsi_x_hat + Roe_v*xsi_y_hat)) / (xsi_x_hat + xsi_y_hat)) / Roe_rho;
		dw4 = (-1.0*((Roe_a*(2.0*dU4 - 2.0*dU2*Roe_u - 2.0*dU3*Roe_v + dU1*(Roe_u*Roe_u + Roe_v*Roe_v))) / (-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v)) + (-1.0*(dU2*xsi_x_hat) + dU1*Roe_u*xsi_x_hat - dU3*xsi_y_hat + dU1*Roe_v*xsi_y_hat) / (xsi_x_hat*xsi_x_hat + xsi_y_hat*xsi_y_hat)) / Roe_rho;

		/*Get contravar velocity*/
		contra_UL = xsi_x*uvelL + xsi_y*vvelL;
		contra_UR = xsi_x*uvelR + xsi_y*vvelR;

		/*Left and right flux*/
		F1L = rhoL*contra_UL / J_F;
		F1R = rhoR*contra_UR / J_F;

		/*avoid expansion shock by modifying eigenvalues*/
		lambda1 = modEval(lambda1, e, Roe_a);
		lambda3 = modEval(lambda3, e, Roe_a);
		lambda4 = modEval(lambda4, e, Roe_a);

		/*Interface flux*/
		F1 = 0.5*(F1L + F1R) - 0.5*( lambda1*dw1*r1 + lambda3*dw3*r3 + lambda4*dw4*r4 )/J_F;

		fluxF[i] = F1;
	}

	/*Calculate other columns of F and G fluxes*/
	for(j=1; j<(n-1); j++){
		for(i=1; i<(m-1); i++){
			/*F flux*/
			/*MUSCL extrapolate left and right states for F-flux*/
			/*The current cell index (the cell to the left of the interface being calculated) is calculated below*/
			index = i + (j+1)*(m+1);  /*Note: there are m+1 rows of cells (including ghost).   First column is a ghost col. so it is skipped*/

			rhoL = ex1[index] + 0.25*epsilon*( (1.0-kappa)*(ex1[index]-ex1[index-1]) + (1.0+kappa)*(ex1[index+1]-ex1[index]) );
			uvelL = ex2[index] + 0.25*epsilon*( (1.0-kappa)*(ex2[index]-ex2[index-1]) + (1.0+kappa)*(ex2[index+1]-ex2[index]) );
			vvelL = ex3[index] + 0.25*epsilon*( (1.0-kappa)*(ex3[index]-ex3[index-1]) + (1.0+kappa)*(ex3[index+1]-ex3[index]) );
			PrsL = ex4[index] + 0.25*epsilon*( (1.0-kappa)*(ex4[index]-ex4[index-1]) + (1.0+kappa)*(ex4[index+1]-ex4[index]) );
 	
			rhoR = ex1[index+1] - 0.25*epsilon*( (1.0+kappa)*(ex1[index+1]-ex1[index]) + (1.0-kappa)*(ex1[index+2]-ex1[index+1]) );
			uvelR = ex2[index+1] - 0.25*epsilon*( (1.0+kappa)*(ex2[index+1]-ex2[index]) + (1.0-kappa)*(ex2[index+2]-ex2[index+1]) );
			vvelR = ex3[index+1] - 0.25*epsilon*( (1.0+kappa)*(ex3[index+1]-ex3[index]) + (1.0-kappa)*(ex3[index+2]-ex3[index+1]) );
			PrsR = ex4[index+1] - 0.25*epsilon*( (1.0+kappa)*(ex4[index+1]-ex4[index]) + (1.0-kappa)*(ex4[index+2]-ex4[index+1]) );

			/*Temp Test Area*/			
			/*rhoL = 1.1;
			uvelL = 400.0;
			vvelL = 25.5;
			PrsL = 100000.0;

			rhoR = 1.11;
			uvelR = 405.5;
			vvelR = 25.6;
			PrsR = 110000.0;*/


			/*Calculate the left and right total enthalpy and conserved variables*/
			U1_L = rhoL;
			U2_L = primToU2(rhoL, uvelL);
			U3_L = primToU3(rhoL, vvelL);
			U4_L = primToU4(rhoL, uvelL, vvelL, PrsL);
			htL = primToht(rhoL, uvelL, vvelL, PrsL);

			U1_R = rhoR;
			U2_R = primToU2(rhoR, uvelR);
			U3_R = primToU3(rhoR, vvelR);
			U4_R = primToU4(rhoR, uvelR, vvelR, PrsR);
			htR = primToht(rhoR, uvelR, vvelR, PrsR);
	
			/*Roe avg variables*/
			R = sqrt(rhoR/rhoL);
			Roe_rho = R*rhoL;
			Roe_u = (R*uvelR + uvelL)/(R + 1.0);
			Roe_v = (R*vvelR + vvelL)/(R + 1.0);
			Roe_ht = (R*htR + htL)/(R + 1.0);
			Roe_a = sqrt( (gamma-1.0)*fabs(Roe_ht - 0.5*(Roe_u*Roe_u + Roe_v*Roe_v) ) );

			/*Get gridmetrics for F-flux interface*/
			x_eta = dx_deta[i + j*m];
			y_eta = dy_deta[i + j*m];
	
			J_F = 0.5*(J[i-1 + j*(m-1)] + J[i + j*(m-1)]);   /*Interpolate J to current F interface*/
	
			xsi_x = y_eta*J_F;			/*Calculate inverse transform metrics at current F interface*/
			xsi_y = -1.0*x_eta*J_F;

			/*Eigenvalues, F flux*/
			lambda1 = xsi_x*Roe_u + xsi_y*Roe_v;
				/*lambda 2, same as 1 but un-needed for this function*/
			lambda3 = lambda1 + Roe_a*sqrt(xsi_x*xsi_x + xsi_y*xsi_y);
			lambda4 = lambda1 - Roe_a*sqrt(xsi_x*xsi_x + xsi_y*xsi_y);
	
			/*Right eigenvectors, F-flux*/
			r1 = 1.0;
			r3 = Roe_rho / (2.0*Roe_a);
			r4 = r3;
	
			/*Wave amplitudes*/
			dU1 = U1_R - U1_L;
			dU2 = U2_R - U2_L;
			dU3 = U3_R - U3_L;
			dU4 = U4_R - U4_L;
	
			xsi_x_hat = xsi_x / sqrt(xsi_x*xsi_x + xsi_y*xsi_y); 
			xsi_y_hat = xsi_y / sqrt(xsi_x*xsi_x + xsi_y*xsi_y);

			dw1 = 2.0*(dU1 + (dU4 + dU1*Roe_ht - dU2*Roe_u - dU3*Roe_v)/(-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v) );
			/*dw2 = (-1.0*(dU3*xsi_x_hat) + dU1*Roe_v*xsi_x_hat + dU2*xsi_y_hat - dU1*Roe_u*xsi_y_hat) / (Roe_rho*(xsi_x_hat*xsi_x_hat + xsi_y_hat*xsi_y_hat  ) );*/
			dw3 = (-1.0*((Roe_a*(2.0*dU4 - 2.0*dU2*Roe_u - 2.0*dU3*Roe_v + dU1*(Roe_u*Roe_u + Roe_v*Roe_v)))/(-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v)) + (dU2*xsi_x_hat + dU3*xsi_y_hat - dU1*(Roe_u*xsi_x_hat + Roe_v*xsi_y_hat)) / (xsi_x_hat + xsi_y_hat)) / Roe_rho;
			dw4 = (-1.0*((Roe_a*(2.0*dU4 - 2.0*dU2*Roe_u - 2.0*dU3*Roe_v + dU1*(Roe_u*Roe_u + Roe_v*Roe_v))) / (-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v)) + (-1.0*(dU2*xsi_x_hat) + dU1*Roe_u*xsi_x_hat - dU3*xsi_y_hat + dU1*Roe_v*xsi_y_hat) / (xsi_x_hat*xsi_x_hat + xsi_y_hat*xsi_y_hat)) / Roe_rho;

			/*Get contravar velocity*/
			contra_UL = xsi_x*uvelL + xsi_y*vvelL;
			contra_UR = xsi_x*uvelR + xsi_y*vvelR;
	
			/*Left and right flux*/
			F1L = rhoL*contra_UL / J_F;
			F1R = rhoR*contra_UR / J_F;

			/*avoid expansion shock by modifying eigenvalues*/
			lambda1 = modEval(lambda1, e, Roe_a);
			lambda3 = modEval(lambda3, e, Roe_a);
			lambda4 = modEval(lambda4, e, Roe_a);

			/*Interface flux*/
			F1 = 0.5*(F1L + F1R) - 0.5*( lambda1*dw1*r1 + lambda3*dw3*r3 + lambda4*dw4*r4 )/J_F;
	
			fluxF[i+j*m] = F1;


			/*G flux*/
			/*MUSCL extrapolate left and right states for G-flux*/
			/*The current cell index (the cell to the left of the interface being calculated) is calculated below*/
			index = i + (j)*(m+1);  /*Note: there are m+1 rows of cells (including ghost).   First column is a ghost col. so it is skipped*/

			rhoL = ex1[index] + 0.25*epsilon*( (1.0-kappa)*(ex1[index]-ex1[index-(m+1)]) + (1.0+kappa)*(ex1[index+(m+1)]-ex1[index]) );
			uvelL = ex2[index] + 0.25*epsilon*( (1.0-kappa)*(ex2[index]-ex2[index-(m+1)]) + (1.0+kappa)*(ex2[index+(m+1)]-ex2[index]) );
			vvelL = ex3[index] + 0.25*epsilon*( (1.0-kappa)*(ex3[index]-ex3[index-(m+1)]) + (1.0+kappa)*(ex3[index+(m+1)]-ex3[index]) );
			PrsL = ex4[index] + 0.25*epsilon*( (1.0-kappa)*(ex4[index]-ex4[index-(m+1)]) + (1.0+kappa)*(ex4[index+(m+1)]-ex4[index]) );
 	
			rhoR = ex1[index+(m+1)] - 0.25*epsilon*( (1.0+kappa)*(ex1[index+(m+1)]-ex1[index]) + (1.0-kappa)*(ex1[index+2*(m+1)]-ex1[index+(m+1)]) );
			uvelR = ex2[index+(m+1)] - 0.25*epsilon*( (1.0+kappa)*(ex2[index+(m+1)]-ex2[index]) + (1.0-kappa)*(ex2[index+2*(m+1)]-ex2[index+(m+1)]) );
			vvelR = ex3[index+(m+1)] - 0.25*epsilon*( (1.0+kappa)*(ex3[index+(m+1)]-ex3[index]) + (1.0-kappa)*(ex3[index+2*(m+1)]-ex3[index+(m+1)]) );
			PrsR = ex4[index+(m+1)] - 0.25*epsilon*( (1.0+kappa)*(ex4[index+(m+1)]-ex4[index]) + (1.0-kappa)*(ex4[index+2*(m+1)]-ex4[index+(m+1)]) );

			/*Calculate the left and right total enthalpy and conserved variables*/
			U1_L = rhoL;
			U2_L = primToU2(rhoL, uvelL);
			U3_L = primToU3(rhoL, vvelL);
			U4_L = primToU4(rhoL, uvelL, vvelL, PrsL);
			htL = primToht(rhoL, uvelL, vvelL, PrsL);

			U1_R = rhoR;
			U2_R = primToU2(rhoR, uvelR);
			U3_R = primToU3(rhoR, vvelR);
			U4_R = primToU4(rhoR, uvelR, vvelR, PrsR);
			htR = primToht(rhoR, uvelR, vvelR, PrsR);
	
			/*Roe avg variables*/
			R = sqrt(rhoR/rhoL);
			Roe_rho = R*rhoL;
			Roe_u = (R*uvelR + uvelL)/(R + 1.0);
			Roe_v = (R*vvelR + vvelL)/(R + 1.0);
			Roe_ht = (R*htR + htL)/(R + 1.0);
			Roe_a = sqrt( (gamma-1.0)*fabs(Roe_ht - 0.5*(Roe_u*Roe_u + Roe_v*Roe_v) ) );

			/*Get gridmetrics for G-flux interface*/
			x_xsi = dx_dxsi[(i-1) + j*(m-1)];
			y_xsi = dy_dxsi[(i-1) + j*(m-1)];
	
			J_G = 0.5*(J[(i-1) + j*(m-1)] + J[(i-1) + (j-1)*(m-1)]);   /*Interpolate J to current G interface*/
	
			eta_x = -1.0*y_xsi*J_G;			/*Calculate inverse transform metrics at current G interface*/
			eta_y = x_xsi*J_G;

			/*Eigenvalues, F flux*/
			lambda1 = eta_x*Roe_u + eta_y*Roe_v;
				/*lambda 2, same as 1 but un-needed for this function*/
			lambda3 = lambda1 + Roe_a*sqrt(eta_x*eta_x + eta_y*eta_y);
			lambda4 = lambda1 - Roe_a*sqrt(eta_x*eta_x + eta_y*eta_y);
	
			/*Right eigenvectors, F-flux*/
			r1 = 1.0;
			r3 = Roe_rho / (2.0*Roe_a);
			r4 = r3;
	
			/*Wave amplitudes*/
			dU1 = U1_R - U1_L;
			dU2 = U2_R - U2_L;
			dU3 = U3_R - U3_L;
			dU4 = U4_R - U4_L;

			eta_x_hat = eta_x / sqrt(eta_x*eta_x + eta_y*eta_y); 
			eta_y_hat = eta_y / sqrt(eta_x*eta_x + eta_y*eta_y);

			dw1 = 2.0*(dU1 + (dU4 + dU1*Roe_ht - dU2*Roe_u - dU3*Roe_v)/(-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v));
			/*dw2 = (-1.0*(dU3*eta_x_hat) + dU1*Roe_v*eta_x_hat + dU2*eta_y_hat - dU1*Roe_u*eta_y_hat)/(Roe_rho*(eta_x_hat*eta_x_hat + eta_y_hat*eta_y_hat)); */
			dw3 = (-1.0*((Roe_a*(2.0*dU4 - 2.0*dU2*Roe_u - 2.0*dU3*Roe_v + dU1*(Roe_u*Roe_u + Roe_v*Roe_v))) / (-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v)) + (dU2*eta_x_hat + dU3*eta_y_hat - dU1*(Roe_u*eta_x_hat + Roe_v*eta_y_hat)) / (eta_x_hat*eta_x_hat + eta_y_hat*eta_y_hat))/Roe_rho;
			dw4 = (-1.0*((Roe_a*(2.0*dU4 - 2.0*dU2*Roe_u - 2.0*dU3*Roe_v + dU1*(Roe_u*Roe_u + Roe_v*Roe_v))) / (-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v)) + (-1.0*(dU2*eta_x_hat) + dU1*Roe_u*eta_x_hat - dU3*eta_y_hat + dU1*Roe_v*eta_y_hat) / (eta_x_hat*eta_x_hat + eta_y_hat*eta_y_hat))/Roe_rho;

			
			/*Get contravar velocity*/
			contra_VL = eta_x*uvelL + eta_y*vvelL;
			contra_VR = eta_x*uvelR + eta_y*vvelR;
	
			/*Left and right flux*/
			G1L = rhoL*contra_VL / J_G;
			G1R = rhoR*contra_VR / J_G;

			/*avoid expansion shock by modifying eigenvalues*/
			lambda1 = modEval(lambda1, e, Roe_a);
			lambda3 = modEval(lambda3, e, Roe_a);
			lambda4 = modEval(lambda4, e, Roe_a);

			/*Interface flux*/
			G1 = 0.5*(G1L + G1R) - 0.5*( lambda1*dw1*r1 + lambda3*dw3*r3 + lambda4*dw4*r4 )/J_G;
	
			fluxG[(i-1)+j*(m-1)] = G1;
		}
		/*G-flux, bottom row (cells at i=m-1 ; the above loop ends with i=m-2)*/
		/*MUSCL extrapolate left and right states for G-flux*/
		/*The current cell index (the cell to the left of the interface being calculated) is calculated below*/
		index = (m-1) + (j)*(m+1);  /*Note: there are m+1 rows of cells (including ghost).   First column is a ghost col. so it is skipped*/

		rhoL = ex1[index] + 0.25*epsilon*( (1.0-kappa)*(ex1[index]-ex1[index-(m+1)]) + (1.0+kappa)*(ex1[index+(m+1)]-ex1[index]) );
		uvelL = ex2[index] + 0.25*epsilon*( (1.0-kappa)*(ex2[index]-ex2[index-(m+1)]) + (1.0+kappa)*(ex2[index+(m+1)]-ex2[index]) );
		vvelL = ex3[index] + 0.25*epsilon*( (1.0-kappa)*(ex3[index]-ex3[index-(m+1)]) + (1.0+kappa)*(ex3[index+(m+1)]-ex3[index]) );
		PrsL = ex4[index] + 0.25*epsilon*( (1.0-kappa)*(ex4[index]-ex4[index-(m+1)]) + (1.0+kappa)*(ex4[index+(m+1)]-ex4[index]) );
 	
		rhoR = ex1[index+(m+1)] - 0.25*epsilon*( (1.0+kappa)*(ex1[index+(m+1)]-ex1[index]) + (1.0-kappa)*(ex1[index+2*(m+1)]-ex1[index+(m+1)]) );
		uvelR = ex2[index+(m+1)] - 0.25*epsilon*( (1.0+kappa)*(ex2[index+(m+1)]-ex2[index]) + (1.0-kappa)*(ex2[index+2*(m+1)]-ex2[index+(m+1)]) );
		vvelR = ex3[index+(m+1)] - 0.25*epsilon*( (1.0+kappa)*(ex3[index+(m+1)]-ex3[index]) + (1.0-kappa)*(ex3[index+2*(m+1)]-ex3[index+(m+1)]) );
		PrsR = ex4[index+(m+1)] - 0.25*epsilon*( (1.0+kappa)*(ex4[index+(m+1)]-ex4[index]) + (1.0-kappa)*(ex4[index+2*(m+1)]-ex4[index+(m+1)]) );

		/*Calculate the left and right total enthalpy and conserved variables*/
		U1_L = rhoL;
		U2_L = primToU2(rhoL, uvelL);
		U3_L = primToU3(rhoL, vvelL);
		U4_L = primToU4(rhoL, uvelL, vvelL, PrsL);
		htL = primToht(rhoL, uvelL, vvelL, PrsL);

		U1_R = rhoR;
		U2_R = primToU2(rhoR, uvelR);
		U3_R = primToU3(rhoR, vvelR);
		U4_R = primToU4(rhoR, uvelR, vvelR, PrsR);
		htR = primToht(rhoR, uvelR, vvelR, PrsR);
	
		/*Roe avg variables*/
		R = sqrt(rhoR/rhoL);
		Roe_rho = R*rhoL;
		Roe_u = (R*uvelR + uvelL)/(R + 1.0);
		Roe_v = (R*vvelR + vvelL)/(R + 1.0);
		Roe_ht = (R*htR + htL)/(R + 1.0);
		Roe_a = sqrt( (gamma-1.0)*fabs(Roe_ht - 0.5*(Roe_u*Roe_u + Roe_v*Roe_v) ) );

		/*Get gridmetrics for G-flux interface*/
		x_xsi = dx_dxsi[(m-1-1) + j*(m-1)];
		y_xsi = dy_dxsi[(m-1-1) + j*(m-1)];
	
		J_G = 0.5*(J[(m-1-1) + j*(m-1)] + J[(m-1-1) + (j-1)*(m-1)]);   /*Interpolate J to current G interface*/
	
		eta_x = -1.0*y_xsi*J_G;			/*Calculate inverse transform metrics at current G interface*/
		eta_y = x_xsi*J_G;

		/*Eigenvalues, F flux*/
		lambda1 = eta_x*Roe_u + eta_y*Roe_v;
			/*lambda 2, same as 1 but un-needed for this function*/
		lambda3 = lambda1 + Roe_a*sqrt(eta_x*eta_x + eta_y*eta_y);
		lambda4 = lambda1 - Roe_a*sqrt(eta_x*eta_x + eta_y*eta_y);
	
		/*Right eigenvectors, F-flux*/
		r1 = 1.0;
		r3 = Roe_rho / (2.0*Roe_a);
		r4 = r3;
	
		/*Wave amplitudes*/
		dU1 = U1_R - U1_L;
		dU2 = U2_R - U2_L;
		dU3 = U3_R - U3_L;
		dU4 = U4_R - U4_L;

		eta_x_hat = eta_x / sqrt(eta_x*eta_x + eta_y*eta_y); 
		eta_y_hat = eta_y / sqrt(eta_x*eta_x + eta_y*eta_y);

		dw1 = 2.0*(dU1 + (dU4 + dU1*Roe_ht - dU2*Roe_u - dU3*Roe_v)/(-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v));
		/*dw2 = (-1.0*(dU3*eta_x_hat) + dU1*Roe_v*eta_x_hat + dU2*eta_y_hat - dU1*Roe_u*eta_y_hat)/(Roe_rho*(eta_x_hat*eta_x_hat + eta_y_hat*eta_y_hat)); */
		dw3 = (-1.0*((Roe_a*(2.0*dU4 - 2.0*dU2*Roe_u - 2.0*dU3*Roe_v + dU1*(Roe_u*Roe_u + Roe_v*Roe_v))) / (-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v)) + (dU2*eta_x_hat + dU3*eta_y_hat - dU1*(Roe_u*eta_x_hat + Roe_v*eta_y_hat)) / (eta_x_hat*eta_x_hat + eta_y_hat*eta_y_hat))/Roe_rho;
		dw4 = (-1.0*((Roe_a*(2.0*dU4 - 2.0*dU2*Roe_u - 2.0*dU3*Roe_v + dU1*(Roe_u*Roe_u + Roe_v*Roe_v))) / (-2.0*Roe_ht + Roe_u*Roe_u + Roe_v*Roe_v)) + (-1.0*(dU2*eta_x_hat) + dU1*Roe_u*eta_x_hat - dU3*eta_y_hat + dU1*Roe_v*eta_y_hat) / (eta_x_hat*eta_x_hat + eta_y_hat*eta_y_hat))/Roe_rho;

			
		/*Get contravar velocity*/
		contra_VL = eta_x*uvelL + eta_y*vvelL;
		contra_VR = eta_x*uvelR + eta_y*vvelR;
	
		/*Left and right flux*/
		G1L = rhoL*contra_VL / J_G;
		G1R = rhoR*contra_VR / J_G;

		/*avoid expansion shock by modifying eigenvalues*/
		lambda1 = modEval(lambda1, e, Roe_a);
		lambda3 = modEval(lambda3, e, Roe_a);
		lambda4 = modEval(lambda4, e, Roe_a);

		/*Interface flux*/
		G1 = 0.5*(G1L + G1R) - 0.5*( lambda1*dw1*r1 + lambda3*dw3*r3 + lambda4*dw4*r4 )/J_G;
	
		fluxG[(m-2)+j*(m-1)] = G1;
	}	

}

void boundaryFlux(double *fluxF, double *fluxG, mwSize m, mwSize n, double *xMatrix, double *yMatrix, double *param)
{
	mwSize i, j;
	int index;
	double x,y;

	/*Get the parameters needed for exact solution from the param vector*/
	double rho1 = param[0];  /*Upstream conditions*/
	double uvel1 = param[1];
	double vvel1 = param[2];
	double P1 = param[3];

	double rho2 = param[4]; /*Downstream condtions*/
	double uvel2 = param[5];
	double vvel2 = param[6];
	double P2 = param[7];

	double u1 = param[8];   /*Forward angle*/
	double u2 = param[9];	/*Rear angle*/

	double angle, ACF;
	double solution1, solution2;
	double deltax, deltay;

	/*Set inflow boundary (first row of F flux--i=0, j=0 to (n-2) )*/
	for(j=0; j<(n-1); j++){
		deltay = yMatrix[(j+1)*m] - yMatrix[j*m];
		fluxF[j*m] = rho1*uvel1*deltay;  /*The inflow boundary is always upstream of the expansion fan*/
	}

	/*Set lower wall boundary (first col of G flux--i=0 to (m-2), j=0 )*/
	for(i=0; i<(m-1); i++){
		fluxG[i] = 0.0;		/*There is no mass flux through the wall*/
	}

	/*Set outflow boundary (last row of F flux--i=m-1, j=0 to (n-2) )*/
	x = xMatrix[(m-1)];   /*x is constant along outflow boundary*/
	for(j=0; j<(n-1); j++){
		y = 0.5*(yMatrix[(m-1) + j*m] + yMatrix[(m-1) + (j+1)*m]);  /*Get center of flux face*/
		deltay = yMatrix[(m-1) + (j+1)*m] - yMatrix[(m-1) + j*m];

		/*Calculate current angle of vector from origin to node, relative to bottom of domain--note atan2(0,0)=0 at origin*/
		angle = atan2(y, x);
		/*Calculate averaging coefficient--this will be used in the following equation to linearly interpolate between upstream*/
		/*and downstream conditions within the fan*/
		ACF = fmax( (angle - u1)/(u2 - u1) , 0.0 );
		ACF = fmin( ACF, 1.0 ); 

		/*Exact solution for rho and uvel*/
		solution1 = calcExact(ACF, rho1, rho2);
		solution2 = calcExact(ACF, uvel1, uvel2);

		/*F-flux = rho*uvel*/
		fluxF[(m-1) + j*m] = solution1*solution2*deltay;
	}

	/*Set 'upper' outflow boundary (last col of G flux--i=0 to (m-2), j=n-1 )*/
	y = yMatrix[(n-1)*m];   /*y is constant along outflow boundary*/
	for(i=0; i<(m-1); i++){
		x = 0.5*(xMatrix[i + (n-1)*m] + xMatrix[(i+1) + (n-1)*m]);  /*Get center of flux face*/
		deltax = xMatrix[(i+1) + (n-1)*m] - xMatrix[i + (n-1)*m];

		/*Calculate current angle of vector from origin to node, relative to bottom of domain--note atan2(0,0)=0 at origin*/
		angle = atan2(y, x);
		/*Calculate averaging coefficient--this will be used in the following equation to linearly interpolate between upstream*/
		/*and downstream conditions within the fan*/
		ACF = fmax( (angle - u1)/(u2 - u1) , 0.0 );
		ACF = fmin( ACF, 1.0 ); 

		/*Exact solution for rho and vvel*/
		solution1 = calcExact(ACF, rho1, rho2);
		solution2 = calcExact(ACF, vvel1, vvel2);

		/*G-flux = rho*vvel*/
		fluxG[i + (n-1)*(m-1)] = solution1*solution2*deltax;
	}

}

double TruncError(double *fluxF, double *fluxG, double *J, mwSize m, mwSize n, double *TEMat)
{
	int i,j;
	double TE = 0.0;	/*total TE for the entire grid*/
	double c = 0.0;		/*Temp variable for Kahan Algorithm--holds error in sum*/
	double y,t; 		/*Temp vars, Kahan algorithm*/

	double FL, FR, GB, GT;  /*Left, right, top, bottom fluxes for given cell*/
	double localTE; 	/*Truncation error of current cell*/
	
	for(j=0; j<(n-1); j++){
		for(i=0; i<(m-1); i++){
			
			FL = fluxF[i + j*(m)];
			FR = fluxF[i+1 + j*(m)];
			GB = fluxG[i + j*(m-1)];
			GT = fluxG[i + (j+1)*(m-1)];

			localTE = FL - FR + GB - GT;
			TEMat[i + j*(m-1)] = localTE;	/*Store the TE of this cell in the TE matrix*/
			//localTE = localTE/J[i + j*(m-1)];  /*Multiply by grid metric of current cell*/
			localTE = localTE*localTE;	/*Square it, to take L2 norm*/

			/*Sum the local TE's over all cells*/
			/*Uses Kahan's algorithm for better accuracy*/
			y = localTE - c;
			t = TE + y;
			c = (t - TE) - y;
			TE = t;

		}
	}
	
	double cells = (double)((m-1)*(n-1));  /*size of domain*/

	TE = sqrt(TE / cells );		/*sqrt(sum(TE_i^2)/cells) = L2 norm*/
	return TE;
}



/* Gateway Function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* variable declarations here */
	double * __restrict__ xMatrix;		/* MxN input matrix of x-values */
	double * __restrict__ yMatrix;		/* MxN input matrix of y-values */
	mwSize ncols;			/* size of grid -- x,y nodes*/
	mwSize nrows;			/* size of grid -- x,y nodes*/

	/* create pointers to the real data in the input matrices */
	xMatrix = mxGetPr(prhs[0]);
	yMatrix = mxGetPr(prhs[1]);

	/* get dimensions of the input matrix */
	ncols = mxGetN(prhs[0]);
	nrows = mxGetM(prhs[0]);

	double * __restrict__ dx_dxsi;		/* (M-1)xN output matrix  */
	double * __restrict__ dy_dxsi;		/* (M-1)xN output matrix  */
	double * __restrict__ dx_deta;		/* Mx(N-1) output matrix  */
	double * __restrict__ dy_deta;		/* Mx(N-1) output matrix  */
	double * __restrict__ J;				/* (M-1)x(N-1) output matrix */
	double * __restrict__ exact1;			/* (M+1)x(N+1) matrix*/
	double * __restrict__ exact2;			/* (M+1)x(N+1) matrix*/
	double * __restrict__ exact3;			/* (M+1)x(N+1) matrix*/
	double * __restrict__ exact4;			/* (M+1)x(N+1) matrix*/
	double * __restrict__ parameters;			/*Parameters vector*/
	double * __restrict__ fluxF;
	double * __restrict__ fluxG;
	double * __restrict__ TEMatrix;

	/* get a pointer to the real data in the output matrices--!!note they are passed on RHS!! */
	dx_dxsi = mxGetPr(prhs[2]);
	dy_dxsi = mxGetPr(prhs[3]);
	dx_deta = mxGetPr(prhs[4]);
	dy_deta = mxGetPr(prhs[5]);
	J = mxGetPr(prhs[6]);
	exact1 = mxGetPr(prhs[7]);
	exact2 = mxGetPr(prhs[8]);
	exact3 = mxGetPr(prhs[9]);
	exact4 = mxGetPr(prhs[10]);
	parameters = mxGetPr(prhs[11]);
	fluxF = mxGetPr(prhs[12]);
	fluxG = mxGetPr(prhs[13]);
	TEMatrix = mxGetPr(prhs[14]);
	
	/* computational subroutines */
	/* calculate grid metrics */
	metrics(xMatrix, yMatrix, dx_dxsi, dy_dxsi, dx_deta, dy_deta, J, nrows, ncols);

	/* calculate exact solution*/
	setExact(xMatrix, yMatrix, exact1, exact2, exact3, exact4, nrows, ncols, parameters);

	/* compute interior fluxes*/
	interiorFlux(fluxF, fluxG, dx_dxsi, dy_dxsi, dx_deta, dy_deta, J, exact1, exact2, exact3, exact4, nrows, ncols, parameters);

	/* set boundary fluxes*/
	boundaryFlux(fluxF, fluxG, nrows, ncols, xMatrix, yMatrix, parameters);	

	/* compute L2 of residuals to estimate TE*/
	double TruncErr =  TruncError(fluxF, fluxG, J, nrows, ncols, TEMatrix);

	/* return TE*/
	plhs[0] = mxCreateDoubleScalar(TruncErr);
	
}
