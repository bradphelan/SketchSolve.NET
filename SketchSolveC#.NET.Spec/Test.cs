using System;
using NUnit.Framework;
using FluentAssertions;
using SketchSolve;
using System.Linq;
using FluentAssertions.Numeric;
using NUnit.Framework.Constraints;

namespace SketchSolve.Spec
{
    public static class FluentExtensionsX
    {
        public static AndConstraint<NumericAssertions<double>> BeApproximately
            (this NumericAssertions<double> This, double val, double eps)
        {
            return This.BeInRange(val - eps / 2, val + eps / 2);
        }
    }

    [TestFixture()]
    public class Solver
    {
        [Test()]
        public void HorizontalConstraintShouldWork()
        {

            var line = new Line() { p1 = new Point(0, 1), p2 = new Point(2, 3) };

            SketchSolve.Solver.solve
                        (true
                          , line.IsHorizontal()
            );

            line.p1.y.Value.Should().BeApproximately(line.p2.y.Value, 0.001);
        }

        [Test()]
        public void VerticalConstraintShouldWork()
        {

            var line = new Line() { p1 = new Point(0, 1), p2 = new Point(2, 3) };

            SketchSolve.Solver.solve
                        (true
                          , line.IsVertical()
            );

            line.p1.x.Value.Should().BeApproximately(line.p2.x.Value, 0.001);
        }

        [Test()]
        public void PointOnPointConstraintShouldWork()
        {
            var line1 = new Line() { p1 = new Point(0, 1), p2 = new Point(2, 3) };
            var line2 = new Line() { p1 = new Point(10, 100), p2 = new Point(200, 300) };

            SketchSolve.Solver.solve
                (true
                 , line1.p1.IsColocated(line2.p2));

            line1.p1.x.Value.Should().BeApproximately(line2.p2.x.Value, 0.001);
            line1.p1.y.Value.Should().BeApproximately(line2.p2.y.Value, 0.001);

        }

        [Test()]
        public void InternalAngleConstraintShouldWork()
        {
            for (int i = 1; i < 10; i++)
            {
                var line1 = new Line() { p1 = new Point(0, 0, false), p2 = new Point(10, 0, false, true) };
                var line2 = new Line() { p1 = new Point(0, 0, false), p2 = new Point(10, -1, false) };

                Console.WriteLine(i);
                var a = Math.PI / 2 / 3;

                SketchSolve.Solver.solve
                    (true
                     , line1.HasInternalAngle(line2, new Parameter(a, false)));

                line1
                    .cosine(line2)
                    .Should()
                    .BeApproximately(Math.Cos(a), 0.001);

            }

        }


            }

        }

    }
}

