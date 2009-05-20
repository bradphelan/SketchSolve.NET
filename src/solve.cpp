 /*
 * solve.cpp
 *
 *  Created on: May 4, 2009
 *      Author: Jonathan George
 *      This
 *      Copyright (c) 2009, Jonathan George
 *      This program is released under the BSD license. See the file COPYING for details.
 *
 */

#include "solve.h"
#include <math.h>
#include <stdlib.h>

using namespace std;

int solve(double  x[],int xLength, constraint * cons, int consLength, int isFine)
{
	//integer to keep track of how many times calc is called
	int ftimes=0;
	//Calculate Function at the starting point:
	double f0;
	f0 = calc(cons,consLength);
	ftimes++;
	//Calculate the gradient at the starting point:

	//Calculate the gradient
	//gradF=x;
	double grad[xLength]; //The gradient vector (1xn)
	double norm; //The norm of the gradient vector
	double f1,f2,f3,alpha1,alpha2,alpha3,alphaStar;
	norm = 0;
	for(int j=0;j<xLength;j++)
	{
		x[j]= x[j]+pert;
		grad[j]=(calc(cons,consLength)-f0)/pert;
		ftimes++;
		//cout<<"gradient: "<<grad[j]<<endl;
		x[j]-=pert;
		norm = norm+(grad[j]*grad[j]);
	}
	norm = sqrt(norm);
	//Estimate the norm of N

	//Initialize N and calculate s
	double s [xLength]; //The current search direction
	double N [xLength][xLength]; //The estimate of the Hessian inverse
	for(int i=0;i<xLength;i++)
	{
		for(int j=0;j<xLength;j++)
		{
			if(i==j)
			{
				//N[i][j]=norm; //Calculate a scaled identity matrix as a Hessian inverse estimate
				N[i][j]=grad[i]/norm;
				if(N[i][j]<0) N[i][j]=-N[i][j];
				s[i]=-grad[i]/norm; //Calculate the initial search vector

			}
			else N[i][j]=0;
		}
	}
	double fnew;
	fnew=f0+1; 	//make fnew greater than fold
	double alpha=1; //Initial search vector multiplier

	double xold[xLength]; //Storage for the previous design variables
	double fold;
	for(int i=0;i<xLength;i++)
	{
		xold[i]=x[i];//Copy last values to xold
	}

	///////////////////////////////////////////////////////
	/// Start of line search
	///////////////////////////////////////////////////////

	//Make the initial position alpha1
	alpha1=0;
	f1 = f0;

	//Take a step of alpha=1 as alpha2
	alpha2=1;
	for(int i=0;i<xLength;i++)
	{
		x[i]=xold[i]+alpha2*s[i];//calculate the new x
	}
	f2 = calc(cons,consLength);
	ftimes++;

	//Take a step of alpha 3 that is 2*alpha2
	alpha3 = alpha*2;
	for(int i=0;i<xLength;i++)
	{
		x[i]=xold[i]+alpha3*s[i];//calculate the new x
	}
	f3=calc(cons,consLength);
	ftimes++;

	//Now reduce or lengthen alpha2 and alpha3 until the minimum is
	//Bracketed by the triplet f1>f2<f3
	while(f2>f1 || f2>f3)
	{
		if(f2>f1)
		{
			//If f2 is greater than f1 then we shorten alpha2 and alpha3 closer to f1
			//Effectively both are shortened by a factor of two.
			alpha3=alpha2;
			f3=f2;
			alpha2=alpha2/2;
			for(int i=0;i<xLength;i++)
			{
				x[i]=xold[i]+alpha2*s[i];//calculate the new x
			}
			f2=calc(cons,consLength);
			ftimes++;
		}

		else if(f2>f3)
		{
			//If f2 is greater than f3 then we length alpah2 and alpha3 closer to f1
			//Effectively both are lengthened by a factor of two.
			alpha2=alpha3;
			f2=f3;
			alpha3=alpha3*2;
			for(int i=0;i<xLength;i++)
			{
				x[i]=xold[i]+alpha3*s[i];//calculate the new x
			}
			f3=calc(cons,consLength);
			ftimes++;

		}
	}
	// get the alpha for the minimum f of the quadratic approximation
	alphaStar= alpha2+((alpha2-alpha1)*(f1-f3))/(3*(f1-2*f2+f3));

	//Guarantee that the new alphaStar is within the bracket
	if(alphaStar>alpha3 || alphaStar<alpha1) alphaStar=alpha2;

	/// Set the values to alphaStar
	for(int i=0;i<xLength;i++)
	{
		x[i]=xold[i]+alphaStar*s[i];//calculate the new x
	}
	fnew=calc(cons,consLength);
	ftimes++;
	fold=fnew;
	/*
	cout<<"F at alphaStar: "<<fnew<<endl;
	cout<<"alphaStar: "<<alphaStar<<endl;
	cout<<"F0: "<<f0<<endl;
	cout<<"F1: "<<f1<<endl;
	cout<<"F2: "<<f2<<endl;
	cout<<"F3: "<<f3<<endl;
	cout<<"Alpha1: "<<alpha1<<endl;
	cout<<"Alpha2: "<<alpha2<<endl;
	cout<<"Alpha3: "<<alpha3<<endl;
	*/

	/////////////////////////////////////
	///end of line search
	/////////////////////////////////////






	double deltaX[xLength];
	double gradnew[xLength];
	double gamma[xLength];
	double bottom=0;
	double deltaXtDotGamma;
	double gammatDotN[1][xLength];
	double gammatDotNDotGamma=0;
	double firstTerm=0;
	double FirstSecond[xLength][xLength];
	double deltaXDotGammatDotN[xLength][xLength];
	double gammatDotDeltaXt[xLength][xLength];
	double NDotGammaDotDeltaXt[xLength][xLength];
	double deltaXnorm=1;

	int iterations=1;
	int steps;

	///Calculate deltaX
	for(int i=0;i<xLength;i++)
	{
		deltaX[i]=x[i]-xold[i];//Calculate the difference in x for the Hessian update
	}
	double convergence;
	if(isFine>0) convergence = XconvergenceFine;
	else convergence = XconvergenceRough;
	while(deltaXnorm>convergence && fnew>smallF)
	{
	//////////////////////////////////////////////////////////////////////
	///Start of main loop!!!!
	//////////////////////////////////////////////////////////////////////
	bottom=0;
	deltaXtDotGamma = 0;

	for(int i=0;i<xLength;i++)
	{
		//Calculate the new gradient vector
		x[i]=x[i]+pert;
		gradnew[i]=(calc(cons,consLength)-fnew)/pert;
		ftimes++;
		x[i]=x[i]-pert;
		//Calculate the change in the gradient
		gamma[i]=gradnew[i]-grad[i];
		bottom+=deltaX[i]*gamma[i];

		deltaXtDotGamma += deltaX[i]*gamma[i];

	}

	//make sure that bottom is never 0
	if (bottom==0) bottom=.0000000001;

	//calculate all (1xn).(nxn)

	for(int i=0;i<xLength;i++)
	{
		gammatDotN[0][i]=0;
		for(int j=0;j<xLength;j++)
		{
			gammatDotN[0][i]+=gamma[j]*N[i][j];//This is gammatDotN transpose
		}

	}
	//calculate all (1xn).(nx1)

	gammatDotNDotGamma=0;
	for(int i=0;i<xLength;i++)
	{
		gammatDotNDotGamma+=gammatDotN[0][i]*gamma[i];
	}

	//Calculate the first term

	firstTerm=0;
	firstTerm=1+gammatDotNDotGamma/bottom;

	//Calculate all (nx1).(1xn) matrices
	for(int i=0;i<xLength;i++)
	{
		for(int j=0;j<xLength;j++)
		{
			FirstSecond[i][j]=((deltaX[j]*deltaX[i])/bottom)*firstTerm;
			deltaXDotGammatDotN[i][j]=deltaX[i]*gammatDotN[0][j];
			gammatDotDeltaXt[i][j]=gamma[i]*deltaX[j];
		}
	}

	//Calculate all (nxn).(nxn) matrices

	for(int i=0;i<xLength;i++)
	{
		for(int j=0;j<xLength;j++)
		{
			NDotGammaDotDeltaXt[i][j]=0;
			for(int k=0;k<xLength;k++)
			{
				NDotGammaDotDeltaXt[i][j]+=N[i][k]*gammatDotDeltaXt[k][j];
			}
		}
	}
	//Now calculate the BFGS update on N
	//cout<<"N:"<<endl;
	for(int i=0;i<xLength;i++)
	{

		for(int j=0;j<xLength;j++)
		{
			N[i][j]=N[i][j]+FirstSecond[i][j]-(deltaXDotGammatDotN[i][j]+NDotGammaDotDeltaXt[i][j])/bottom;
			//cout<<" "<<N[i][j]<<" ";
		}
		//cout<<endl;
	}

	//Calculate s
	for(int i=0;i<xLength;i++)
	{
		s[i]=0;
		for(int j=0;j<xLength;j++)
		{
			s[i]+=-N[i][j]*gradnew[j];
		}
	}

	alpha=1; //Initial search vector multiplier


	//copy newest values to the xold
	for(int i=0;i<xLength;i++)
	{
		xold[i]=x[i];//Copy last values to xold
	}
	steps=0;

	///////////////////////////////////////////////////////
	/// Start of line search
	///////////////////////////////////////////////////////

	//Make the initial position alpha1
	alpha1=0;
	f1 = fnew;

	//Take a step of alpha=1 as alpha2
	alpha2=1;
	for(int i=0;i<xLength;i++)
	{
		x[i]=xold[i]+alpha2*s[i];//calculate the new x
	}
	f2 = calc(cons,consLength);
	ftimes++;

	//Take a step of alpha 3 that is 2*alpha2
	alpha3 = alpha2*2;
	for(int i=0;i<xLength;i++)
	{
		x[i]=xold[i]+alpha3*s[i];//calculate the new x
	}
	f3=calc(cons,consLength);
	ftimes++;

	//Now reduce or lengthen alpha2 and alpha3 until the minimum is
	//Bracketed by the triplet f1>f2<f3
	steps=0;
	while(f2>f1 || f2>f3)
	{
		if(f2>f1)
		{
			//If f2 is greater than f1 then we shorten alpha2 and alpha3 closer to f1
			//Effectively both are shortened by a factor of two.
			alpha3=alpha2;
			f3=f2;
			alpha2=alpha2/2;
			for(int i=0;i<xLength;i++)
			{
				x[i]=xold[i]+alpha2*s[i];//calculate the new x
			}
			f2=calc(cons,consLength);
			ftimes++;
		}

		else if(f2>f3)
		{
			//If f2 is greater than f3 then we length alpah2 and alpha3 closer to f1
			//Effectively both are lengthened by a factor of two.
			alpha2=alpha3;
			f2=f3;
			alpha3=alpha3*2;
			for(int i=0;i<xLength;i++)
			{
				x[i]=xold[i]+alpha3*s[i];//calculate the new x
			}
			f3=calc(cons,consLength);
			ftimes++;
		}
		if(steps==4)
			{
				alpha2=1;
				alpha3=2;

				for(int i=0;i<xLength;i++)
				{
					for(int j=0;j<xLength;j++)
					{
						if(i==j)
						{
							N[i][j]=1;
							s[i]=-gradnew[i]; //Calculate the initial search vector
						}
						else N[i][j]=0;
					}
				}
			}
		steps=steps+1;
	}

	// get the alpha for the minimum f of the quadratic approximation
	alphaStar= alpha2+((alpha2-alpha1)*(f1-f3))/(3*(f1-2*f2+f3));

	//Guarantee that the new alphaStar is within the bracket
	if(alphaStar>alpha3 || alphaStar<alpha1)
	{
		alphaStar=alpha2;
	}

	/// Set the values to alphaStar
	for(int i=0;i<xLength;i++)
	{
		x[i]=xold[i]+alphaStar*s[i];//calculate the new x
	}
	fnew=calc(cons,consLength);
	ftimes++;

	/*
	cout<<"F at alphaStar: "<<fnew<<endl;
	cout<<"alphaStar: "<<alphaStar<<endl;
	cout<<"F1: "<<f1<<endl;
	cout<<"F2: "<<f2<<endl;
	cout<<"F3: "<<f3<<endl;
	cout<<"Alpha1: "<<alpha1<<endl;
	cout<<"Alpha2: "<<alpha2<<endl;
	cout<<"Alpha3: "<<alpha3<<endl;
	*/

	/////////////////////////////////////
	///end of line search
	////////////////////////////////////

	deltaXnorm=0;
	for(int i=0;i<xLength;i++)
	{
		deltaX[i]=x[i]-xold[i];//Calculate the difference in x for the hessian update
		deltaXnorm+=deltaX[i]*deltaX[i];
		grad[i]=gradnew[i];
	}
	deltaXnorm=sqrt(deltaXnorm);
	iterations++;
	/////////////////////////////////////////////////////////////
	///End of Main loop
	/////////////////////////////////////////////////////////////
	}
	////Debug

	for(int i=0;i<xLength;i++)
	{
		cout<<"Parameter("<<i<<"): "<<x[i]<<endl;
		//cout<<xold[i]<<endl;
	}
	cout<<"Fnew: "<<fnew<<endl;
	cout<<"Number of Iterations: "<<iterations<<endl;
	cout<<"Number of function calls: "<<ftimes<<endl;

	///End of function
	if(fnew<validSolution) return succsess;
	else return noSolution;
}



double calc(constraint * cons, int consLength)
{
	double error=0;
	double temp,dx,dy,m,n,Ex,Ey,rad1,rad2,t,Xint,Yint,dx2,dy2,hyp1,hyp2,temp2;
	for(int i=0;i<consLength;i++)
	{
		if((cons[i]).type==pointOnPoint)
		{
			//Hopefully avoid this constraint, make coincident points use the same parameters
			error += (P1_x - P2_x) * (P1_x - P2_x) + (P1_y - P2_y) * (P1_y - P2_y);
		}

		if(cons[i].type==P2PDistance)
		{
			error+= (P1_x - P2_x) * (P1_x - P2_x) + (P1_y - P2_y) * (P1_y - P2_y) - distance * distance;

		}

		if(cons[i].type==P2PDistanceVert)
		{
			error+= (P1_y - P2_y) * (P1_y - P2_y) - distance * distance;
		}

		if(cons[i].type==P2PDistanceHorz)
		{
			error+= (P1_x - P2_x) * (P1_x - P2_x) - distance * distance;
		}


		if((cons[i]).type==pointOnCurve)
		{
			dx = L1_P2_x - L1_P1_x;
			dy = L1_P2_y - L1_P1_y;

			m=dy/dx; //Slope
			n=dx/dy; //1/Slope

			if(m<=1)
			{
				//Calculate the expected y point given the x coordinate of the point
				Ey=L1_P1_y+m*(P1_x-L1_P1_x);
				error+=(Ey-P1_y)*(Ey-P1_y);
			}
			else
			{
				//Calculate the expected x point given the y coordinate of the point
				Ex=L1_P1_x+n*(P1_y-L1_P1_y);
				error+=(Ex-P1_x)*(Ex-P1_x);
			}
		}

		if((cons[i]).type==P2LDistance)
		{
			dx = L1_P2_x - L1_P1_x;
			dy = L1_P2_y - L1_P1_y;

			double Xint,Yint;
			t=-(L1_P1_x*dx-P1_x*dx+L1_P1_y*dy-P1_y*dy)/(dx*dx+dy*dy);
			Xint=L1_P1_x+dx*t;
			Yint=L1_P1_y+dy*t;
			temp= hypot((P1_x - Xint),(P1_y - Yint)) - distance;
			error += temp*temp/10;

		}
		if((cons[i]).type==P2LDistanceVert)
		{
			dx = L1_P2_x - L1_P1_x;
			dy = L1_P2_y - L1_P1_y;

			t=(P1_x- L1_P1_x)/dx;
			Yint=L1_P1_y+dy*t;
			temp= fabs((P1_y - Yint)) - distance;
			error += temp*temp;

		}
		if((cons[i]).type==P2LDistanceHorz)
		{
			dx = L1_P2_x - L1_P1_x;
			dy = L1_P2_y - L1_P1_y;

			t=(P1_y- L1_P1_y)/dy;
			Xint=L1_P1_x+dx*t;
			temp= fabs((P1_x - Xint)) - distance;
			error += temp*temp/10;

		}


		if(cons[i].type==horizontal)
		{
			temp= (L1_P1_y-L1_P2_y);
			error+=temp*temp;
		}

		if(cons[i].type==vertical)
		{
			temp=(L1_P1_x-L1_P2_x);
			error+=temp*temp;
		}

		if(cons[i].type==tangentToCircle)
		{
			double dx,dy,Rpx,Rpy,RpxN,RpyN,hyp,error1,error2;
			dx = L1_P2_x-L1_P1_x;
			dy = L1_P2_y-L1_P1_y;
			hyp=hypot(dx,dy);
			//Calculate the expected tangent intersection points
			Rpx =C1_Center_x - dy / hyp * C1_rad;
			Rpy =C1_Center_y + dx / hyp * C1_rad;
			RpxN=C1_Center_x + dy / hyp * C1_rad;
			RpyN=C1_Center_y - dx / hyp * C1_rad;

			error1=(-dy * Rpx + dx * Rpy   + (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y))/hyp;
			error2=(-dy * RpxN + dx * RpyN + (L1_P1_x * L1_P2_y - L1_P2_x * L1_P1_y))/hyp;
			error1=error1 * error1;
			error2=error2 * error2;
			if(error1<error2) error+=error1;
			else error+=error2;

		}

		if(cons[i].type==tangentToArc)
		{
			double dx,dy,Rpx,Rpy,RpxN,RpyN,hyp,error1,error2,rad;
			dx = L1_P2_x - L1_P1_x;
			dy = L1_P2_y - L1_P1_y;
			hyp=hypot(dx,dy);

			rad=hypot(A1_Center_x - A1_Start_x , A1_Center_y - A1_Start_y);
			Rpx=A1_Center_x - dy / hyp * rad;
			Rpy=A1_Center_y + dx / hyp * rad;
			RpxN=A1_Center_x + dy / hyp * rad;
			RpyN=A1_Center_y - dx / hyp * rad;

			m=dy / dx;
			n=dx / dy;
			if(m<=1)
			{
				Ey=L1_P1_y + m * (Rpx - L1_P1_x);
				error1=(Ey - Rpy) * (Ey - Rpy);
				Ey=L1_P1_y + m * (RpxN - L1_P1_x);
				error2=(Ey - RpyN) * (Ey - RpyN);
			}
			else
			{
				Ex=L1_P1_x + n * (Rpy - L1_P1_y);
				error1=(Ex - Rpx) * (Ex - Rpx);
				Ex=L1_P1_x + n * (RpyN - L1_P1_y);
				error2=(Ex - RpxN) * (Ex - RpxN);
			}
			if(error1<error2) error+=error1;
			else error+=error2;
		}

		if(cons[i].type==arcRules)
		{
			rad1=hypot(A1_Center_x - A1_Start_x , A1_Center_y - A1_Start_y);
			rad2=hypot(A1_Center_x - A1_End_x , A1_Center_y - A1_End_y);
			error+=(rad1 - rad2) * (rad1 - rad2);
			}

		if(cons[i].type==lineLength)
		{
			temp=hypot(L1_P2_x - L1_P1_x , L1_P2_y - L1_P1_y) - length;
			error += temp*temp;
		}

		if(cons[i].type==equalLegnth)
		{
			temp=hypot(L1_P2_x - L1_P1_x , L1_P2_y - L1_P1_y) - hypot(L2_P2_x - L2_P1_x , L2_P2_y - L2_P1_y);
			error += temp*temp;
		}

		if(cons[i].type==arcRadius)
		{
			rad1 = hypot(A1_Center_x - A1_Start_x , A1_Center_y - A1_Start_y);
			temp= rad1 - radius ;
			error += temp*temp;
		}

		if(cons[i].type==equalRadiusArcs)
		{
			rad1 = hypot(A1_Center_x - A1_Start_x , A1_Center_y - A1_Start_y);
			rad2 = hypot(A2_Center_x - A2_Start_x , A2_Center_y - A2_Start_y);
			temp = rad1-rad2;
			error += temp*temp;
		}

		if(cons[i].type==equalRadiusCircles)
		{
			temp = C1_rad - C2_rad;
			error += temp*temp;
		}

		if(cons[i].type==equalRadiusCircArc)
		{
			rad1 = hypot(A1_Center_x - A1_Start_x , A1_Center_y - A1_Start_y);
			temp = rad1-C1_rad;
			error += temp*temp;
		}

		if(cons[i].type==concentricArcs)
		{
			error += hypot(A1_Center_x - A2_Center_x , A1_Center_y - A2_Center_y);
		}

		if(cons[i].type==concentricCircles)
		{
			error += hypot(C1_Center_x - C2_Center_x , C1_Center_y - C2_Center_y);
		}

		if(cons[i].type==concentricCircArc)
		{
			error += hypot(A1_Center_x - C1_Center_x , A1_Center_y - C1_Center_y);
		}

		if(cons[i].type==circleRadius)
		{
			error += (C1_rad - radius)*(C1_rad - radius);
		}
		if(cons[i].type==internalAngle)
		{
			dx = L1_P2_x - L1_P1_x;
			dy = L1_P2_y - L1_P1_y;
			dx2 = L2_P2_x - L2_P1_x;
			dy2 = L2_P2_y - L2_P1_y;

			hyp1=hypot(dx,dy);
			hyp2=hypot(dx2,dy2);

			dx=dx/hyp1;
			dy=dy/hyp1;
			dx2=dx2/hyp2;
			dy2=dy2/hyp2;

			temp = dx*dx2+dy*dy2;
			temp2 = cos(angleP);
			error += (temp+temp2)*(temp+temp2);
		}

		if(cons[i].type==externalAngle)
		{
			dx = L1_P2_x - L1_P1_x;
			dy = L1_P2_y - L1_P1_y;
			dx2 = L2_P2_x - L2_P1_x;
			dy2 = L2_P2_y - L2_P1_y;

			hyp1=hypot(dx,dy);
			hyp2=hypot(dx2,dy2);

			dx=dx/hyp1;
			dy=dy/hyp1;
			dx2=dx2/hyp2;
			dy2=dy2/hyp2;

			temp = dx*dx2-dy*dy2;
			temp2 = cos(M_PI-angleP);
			error += (temp+temp2)*(temp+temp2);
		}

		if(cons[i].type==perpendicular)
		{
			dx = L1_P2_x - L1_P1_x;
			dy = L1_P2_y - L1_P1_y;
			dx2 = L2_P2_x - L2_P1_x;
			dy2 = L2_P2_y - L2_P1_y;

			hyp1=hypot(dx,dy);
			hyp2=hypot(dx2,dy2);

			dx=-dx/hyp1;
			dy=dy/hyp1;
			dx2=dx2/hyp2;
			dy2=dy2/hyp2;

			temp = dx*dx2-dy*dy2;
			error += (temp)*(temp);
		}

		if(cons[i].type==parallel)
		{
			dx = L1_P2_x - L1_P1_x;
			dy = L1_P2_y - L1_P1_y;
			dx2 = L2_P2_x - L2_P1_x;
			dy2 = L2_P2_y - L2_P1_y;

			hyp1=hypot(dx,dy);
			hyp2=hypot(dx2,dy2);

			dx=dx/hyp1;
			dy=dy/hyp1;
			dx2=dx2/hyp2;
			dy2=dy2/hyp2;

			temp = dx*dx2-dy*dy2;
			error += (temp)*(temp);
		}
	}
	return error;

}
