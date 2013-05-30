using System;
using NUnit.Framework;
using FluentAssertions;
using SketchSolve;
using System.Linq;

namespace SketchSolve.Spec
{
	[TestFixture()]
	public class Solver
	{
		[Test()]
		public void HorizontalConstraintShouldWork ()
		{

			var parameters = new Parameter[]{
				new Parameter(0),
				new Parameter(1),
				new Parameter(2),
				new Parameter(3)
			};

			var points = new point[]{
				new point(){x = parameters[0], y = parameters[1]},
				new point(){x = parameters[2], y = parameters[3]},
			};

			var lines = new line[]{
				new line(){p1 = points[0], p2 = points[1]}
			};

			var cons = new constraint[] {
				new constraint(){
					type =ConstraintEnum.horizontal,
					line1 = lines[0] 
				}
			};

			var r = SketchSolve.Solver.solve(parameters, cons, true);

			points [0].y.Value.Should ().BeInRange (points[1].y.Value-0.001, points[1].y.Value+0.001); 
			points [0].x.Value.Should ().NotBe (points[1].x.Value); 
		}
	}
}

