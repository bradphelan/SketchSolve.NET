#include "stdafx.h"
#include "CppUnitTest.h"
#include "solve.h"
using namespace std;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
#define PARAM_COUNT 100

namespace SketchSolveTest
{        
    [TestClass]
	public class SketchSolveTest
    {
                [TestMethod]
        [DeploymentItem(@"xxx")]
        void TestMethod1()
        {

            double * parameters = new double[PARAM_COUNT];

            point points[PARAM_COUNT];
            line lines[PARAM_COUNT];

            constraint cons[PARAM_COUNT];

            parameters[0]=0;// First point X coordinate
            parameters[1]=0;// First point Y coordinate
            parameters[2]=10;// Second point X coordinate
            parameters[3]=2;// Second point Y coordinate

            //Define the first point
            points[0].x = &parameters[0]; //Pointer to the X coordinate;
            points[0].y = &parameters[1]; //Pointer to the Y coordinate;

            //Define the second point
            points[1].x = &parameters[2]; //Pointer to the X coordinate;
            points[1].y = &parameters[3]; //Pointer to the Y coordinate;

            //Define the line
            lines[0].p1 = points[0];
            lines[0].p2 = points[1];

            //Define the constraint
            cons[0].type = horizontal; //Make it a horizontal constraint
            cons[0].line1 = lines[0]; 

            int results = solve(&parameters, PARAM_COUNT , cons ,1, true );

        }

    };
}