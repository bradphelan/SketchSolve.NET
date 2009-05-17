//============================================================================
// Name        : SolverPointers.cpp
// Author      : Jonathan George
// Version     :

//     Copyright (c) 2009, Jonathan George
//     This program is released under the BSD license. See the file COPYING for details.
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <math.h>
#include "solve.h"
using namespace std;




double parameters[100],constants[100];
double *P,*C;

point points[100];
line lines[100];

circle circles[100];
arc arcs[100];
constraint cons[100];

int main() {
	P= parameters;
	parameters[0]=0;//1x
	parameters[1]=0;//y
	parameters[2]=10;//x
	parameters[3]=2;//y
	parameters[4]=10;//xstart
	parameters[5]=.5;//y
	parameters[6]=10;//xend
	parameters[7]=5.5;//y
	parameters[8]=9;//xcenter
	parameters[9]=3;//y
	parameters[10]=0;

	point origin;
	double zero = 0;
	origin.x = &zero;
	origin.y = &zero;

	constants[0]=5;
	constants[1]=15;
	constants[2]=2;

	points[0].x = &parameters[0];
	points[0].y = &parameters[1];
	points[1].x = &parameters[2];
	points[1].y = &parameters[3];
	points[2].x = &parameters[4];
	points[2].y = &parameters[5];
	points[3].x = &parameters[6];
	points[3].y = &parameters[7];
	points[4].x = &parameters[8];
	points[4].y = &parameters[9];

	lines[0].p1 = points[0];
	lines[0].p2 = points[1];
	lines[1].p1 = points[3];
	lines[1].p2 = points[4];
	lines[2].p1 = points[4];
	lines[2].p2 = points[0];


	circles[0].center = points[2];
	circles[0].rad = &parameters[10];

	arcs[0].center = points[2];
	arcs[0].start = points[1];
	arcs[0].end = points[3];


	cons[0].type = pointOnPoint;
	cons[0].point1 = origin;
	cons[0].point2 = points[0];

	cons[1].type = horizontal;
	cons[1].line1 = lines[0];

	cons[2].type = horizontal;
	cons[2].line1 = lines[1];

	cons[3].type = vertical;
	cons[3].line1 = lines[2];

	cons[4].type = arcRules;
	cons[4].arc1 = arcs[0];

	cons[5].type = tangentToArc;
	cons[5].line1 = lines[0];
	cons[5].arc1 = arcs[0];

	cons[6].type = concentricCircArc;
	cons[6].arc1 = arcs[0];
	cons[6].circle1 = circles[0];

	cons[7].type = arcRadius;
	cons[7].arc1 = arcs[0];
	cons[7].parameter = &constants[0];

	cons[8].type = lineLength;
	cons[8].line1 = lines[0];
	cons[8].parameter = &constants[1];

	cons[9].type = circleRadius;
	cons[9].circle1 = circles[0];
	cons[9].parameter = &constants[2];

	cons[10].type = tangentToArc;
	cons[10].arc1 = arcs[0];
	cons[10].line1 = lines[1];



	//double x [5];
	//x[0]=45;
	for(int i=0;i<1;i++)
	{
	parameters[0]=1;//1x
	parameters[1]=1;//y
	parameters[2]=14;//x
	parameters[3]=-.5;//y
	parameters[4]=6;//xstart
	parameters[5]=12;//y
	parameters[6]=13.5;//xend
	parameters[7]=11;//y
	parameters[8]=-.5;//xcenter
	parameters[9]=9;//y
	parameters[10]=2.5;
	parameters[11]=2;

	int sol;

	sol=solve(parameters ,11,cons,11,rough);
	if(sol==succsess)
	{
		cout<<"A good Solution was found"<<endl;
	}
	else if(sol==noSolution)
	{
		cout<<"No valid Solutions were found from this start point"<<endl;
	}
	}
	//end
	return 0;
}


