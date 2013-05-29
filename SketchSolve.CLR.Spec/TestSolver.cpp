#include "stdafx.h"
#include "solve.h"

using namespace System;
using namespace System::Text;
using namespace System::Collections::Generic;
using namespace Microsoft::VisualStudio::TestTools::UnitTesting;

namespace SketchSolveCLRSpec
{
	[TestClass]
	public ref class UnitTest
	{
	private:
		TestContext^ testContextInstance;

	public: 
		/// <summary>
		///Gets or sets the test context which provides
		///information about and functionality for the current test run.
		///</summary>
		property Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^ TestContext
		{
			Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^ get()
			{
				return testContextInstance;
			}
			System::Void set(Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^ value)
			{
				testContextInstance = value;
			}
		};

		#pragma region Additional test attributes
		//
		//You can use the following additional attributes as you write your tests:
		//
		//Use ClassInitialize to run code before running the first test in the class
		//[ClassInitialize()]
		//static void MyClassInitialize(TestContext^ testContext) {};
		//
		//Use ClassCleanup to run code after all tests in a class have run
		//[ClassCleanup()]
		//static void MyClassCleanup() {};
		//
		//Use TestInitialize to run code before running each test
		//[TestInitialize()]
		//void MyTestInitialize() {};
		//
		//Use TestCleanup to run code after each test has run
		//[TestCleanup()]
		//void MyTestCleanup() {};
		//
		#pragma endregion 

		[TestMethod]
        void TestHorizontal()
		{
            double constants[100];
            double *parameters = new double[100];

            point points[100];
            line lines[100];

            constraint cons[100];

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

            int r = solve(&parameters, 4, cons, 1, true);

            Assert::AreEqual(*points[0].y, *points[1].y, 0.001);
            Assert::AreNotEqual(*points[0].x, *points[1].x, 0.001);


		};

		[TestMethod]
        void TestVertical()
		{
            double constants[100];
            double *parameters = new double[100];

            point points[100];
            line lines[100];

            constraint cons[100];

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
            cons[0].type = vertical; //Make it a horizontal constraint
            cons[0].line1 = lines[0]; 

            int r = solve(&parameters, 4, cons, 1, true);

            Assert::AreNotEqual(*points[0].y, *points[1].y, 0.001);
            Assert::AreEqual(*points[0].x, *points[1].x, 0.001);

		};
	};
}
