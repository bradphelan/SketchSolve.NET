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


			var line = new line () { p1 = new point(0,1), p2 = new point(2, 3) };
			var parameters = line.ToArray();

			var cons = new constraint[] {
				new constraint(){
					type =ConstraintEnum.horizontal,
					line1 = line 
				}
			};

			var r = SketchSolve.Solver.solve(parameters, cons, true);

			line.p1.y.Value.Should ().BeInRange (line.p2.y.Value-0.001, line.p2.y.Value+0.001); 
		}
	}
}

