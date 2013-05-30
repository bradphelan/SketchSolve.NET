using System;
using NUnit.Framework;
using FluentAssertions;
using SketchSolve;
using System.Linq;
using FluentAssertions.Numeric;

namespace SketchSolve.Spec
{
    public static class FluentExtensionsX {
        public static AndConstraint<NumericAssertions<double>> BeApproximately
            (this NumericAssertions<double> This, double val, double eps)
            {
                return This.BeInRange (val - eps/2, val + eps / 2);
            }
    }

    [TestFixture()]
        public class Solver
        {
            [Test()]
                public void HorizontalConstraintShouldWork ()
                {

                    var line = new line () { p1 = new point(0,1), p2 = new point(2, 3) };

                    SketchSolve.Solver.solve
                        ( true
                          , line.IsHorizontal()
                        );

                    line.p1.y.Value.Should ().BeApproximately (line.p2.y.Value, 0.001);
                }

            [Test()]
                public void VerticalConstraintShouldWork ()
                {

                    var line = new line () { p1 = new point(0,1), p2 = new point(2, 3) };

                    SketchSolve.Solver.solve
                        ( true
                          , line.IsVertical()
                        );

                    line.p1.x.Value.Should ().BeApproximately (line.p2.x.Value, 0.001);
                }


        }
}

