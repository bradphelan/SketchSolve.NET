/*
 * solve.h
 *
 *  Created on: May 4, 2009
 *      Author: Jonathan
 *      Copyright (c) 2009, Jonathan George
 *      This program is released under the BSD license. See the file COPYING for details.
 */
#include <iostream>

#ifndef SOLVE_H_
#define SOLVE_H_

#define pointOnPoint	 0
#define pointToLine 	 1
#define pointOnCurve     2
#define horizontal       3
#define vertical         4
#define angle			 5
#define radiusValue		 6
#define tangentToArc	 7
#define tangentToCircle	 8
#define arcRules 		 9
#define P2PDistance      10
#define P2PDistanceVert  11
#define P2PDistanceHorz  12
#define P2LDistance      13
#define P2LDistanceVert  14
#define P2LDistanceHorz  15
#define lineLength		 16
#define equalLegnth		 17
#define arcRadius		 18
#define equalRadiusArcs      19
#define equalRadiusCircles   20
#define equalRadiusCircArc   21
#define concentricArcs       22
#define concentricCircles    23
#define concentricCircArc    24
#define circleRadius         25


///////////////////////////////////////
/// BFGS Solver parameters
///////////////////////////////////////
#define pert              1e-10
#define XconvergenceRough 1e-4
#define XconvergenceFine  1e-10
#define smallF            1e-20

///////////////////////////////////////////////////////////////////////
/// constraint defines (these make writing constraint equations easier
///////////////////////////////////////////////////////////////////////
#define P1_x		   *cons[i].point1.x
#define P1_y	       *cons[i].point1.y
#define P2_x           *cons[i].point2.x
#define P2_y           *cons[i].point2.y
#define L1_P1_x        *cons[i].line1.p1.x
#define L1_P1_y        *cons[i].line1.p1.y
#define L1_P2_x        *cons[i].line1.p2.x
#define L1_P2_y        *cons[i].line1.p2.y
#define L2_P1_x        *cons[i].line2.p1.x
#define L2_P1_y        *cons[i].line2.p1.y
#define L2_P2_x        *cons[i].line2.p2.x
#define L2_P2_y        *cons[i].line2.p2.y
#define C1_Center_x    *cons[i].circle1.center.x
#define C1_Center_y    *cons[i].circle1.center.y
#define C1_rad         *cons[i].circle1.rad
#define C2_Center_x    *cons[i].circle2.center.x
#define C2_Center_y    *cons[i].circle2.center.y
#define C2_rad         *cons[i].circle2.rad
#define A1_Start_x     *cons[i].arc1.start.x
#define A1_Start_y     *cons[i].arc1.start.y
#define A1_End_x       *cons[i].arc1.end.x
#define A1_End_y       *cons[i].arc1.end.y
#define A1_Center_x    *cons[i].arc1.center.x
#define A1_Center_y    *cons[i].arc1.center.y
#define A2_Start_x     *cons[i].arc2.start.x
#define A2_Start_y     *cons[i].arc2.start.y
#define A2_End_x       *cons[i].arc2.end.x
#define A2_End_y       *cons[i].arc2.end.y
#define A2_Center_x    *cons[i].arc2.center.x
#define A2_Center_y    *cons[i].arc2.center.y
#define length		   *cons[i].parameter
#define distance	   *cons[i].parameter
#define radius		   *cons[i].parameter



struct point
{
	double * x;
	double * y;
};

struct line
{
	point p1;
	point p2;
};

struct arc
{
	point start;
	point end;
	point center;
};

struct circle
{
	point center;
	double *rad;
};

struct constraint
{
	int type;
	point point1;
	point point2;
	line line1;
	line line2;
	circle circle1;
	circle circle2;
	arc arc1;
	arc arc2;
	double *parameter; //radius, length, angle etc...
};





//Function Prototypes
void solve(double  x[],int xLength, constraint * cons, int consLength, int isFine);
double calc(constraint * cons, int consLength);

#endif /* SOLVE_H_ */