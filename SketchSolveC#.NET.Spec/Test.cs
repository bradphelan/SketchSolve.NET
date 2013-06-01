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

            var line = new Line(new Point(0, 1, false), new Point(2, 3, false, true) );

            var error = SketchSolve.Solver.solve
                        (true
                          , line.IsHorizontal()
            );

            error.Should ().BeApproximately (0, 0.0001);

            line.p1.y.Value.Should().BeApproximately(line.p2.y.Value, 0.001);
        }

        [Test()]
        public void VerticalConstraintShouldWork()
        {

            var line = new Line(new Point(0, 1, false), new Point(2, 3, true, false) );

            SketchSolve.Solver.solve
                        (true
                          , line.IsVertical()
            );

            Console.WriteLine (line);
            line.p1.x.Value.Should().BeApproximately(line.p2.x.Value, 0.001);
        }

        [Test()]
        public void PointOnPointConstraintShouldWork()
        {
            var line1 = new Line(new Point(0, 1),  new Point(2, 3, false) );
            var line2 = new Line(new Point(10, 100,false),new Point(200, 300, false) );

            SketchSolve.Solver.solve
                (true
                 , line1.p1.IsColocated(line2.p1));

            line1.p1.x.Value.Should().BeApproximately(line2.p1.x.Value, 0.001);
            line1.p1.y.Value.Should().BeApproximately(line2.p1.y.Value, 0.001);

        }

        [Test()]
        public void InternalAngleConstraintShouldWork()
        {
            for (int i = 1; i < 10; i++)
            {
                var line1 = new Line( new Point(0, 0, false), new Point(10, 0, false, true) );
                var line2 = new Line( new Point(0, 0, false), new Point(10, -1, false) );

                Console.WriteLine(i);
                var a = Math.PI / 2 / 3;

                SketchSolve.Solver.solve
                    (true
                     , line1.HasInternalAngle(line2, new Parameter(a, false)));

                line1
                    .Vector
                    .Cosine(line2.Vector)
                    .Should()
                    .BeApproximately(Math.Cos(a), 0.001);

            }

        }

        [Test()]
        public void ExternalAngleConstraintShouldWork()
        {
            for (int i = 1; i < 10; i++)
            {
                var line1 = new Line( new Point(0, 0, false),  new Point(10, 0, false, true) );
                var line2 = new Line(  new Point(0, 0, false),  new Point(10, -1, false) );

                Console.WriteLine(i);
                var a = Math.PI / 2 / 3;

                SketchSolve.Solver.solve
                    (true
                     , line1.HasExternalAngle(line2, new Parameter(a, false)));

                line1
                    .Vector
                    .Cosine(line2.Vector)
                    .Should()
                    .BeApproximately(Math.Cos(Math.PI - a), 0.001);

            }

        }

        [Test()]
        public void PerpendicularLineConstraintShouldWork()
        {
            for (int i = 1; i < 10; i++)
            {
                var line1 = new Line(new Point(0, 0, false), new Point(10, 0, false, true) );
                var line2 = new Line(new Point(0, 0, false), new Point(10, 10,true,  false) );

                var a = Math.PI / 2 / 3;

                SketchSolve.Solver.solve
                    (true
                     , line1.IsPerpendicularTo(line2));

                Console.WriteLine (line1);
                Console.WriteLine (line2);

                line1
                    .Vector
                    .Dot(line2.Vector)
                    .Should()
                    .BeApproximately(0, 0.001);

            }

        }

        [Test()]
        public void TangentToCircleConstraintShouldWork()
        {
            // Create a fully constrained circle at 0,0 with radius 1
            var circle = new Circle() { center = new Point(0, 0, false), rad = new Parameter(1, false) };

            var v = 1 / Math.Sin(Math.PI / 4);

            var line = new Line(new Point(0, -v, false, false), new Point(35, 0, true, false));

            var r = SketchSolve.Solver.solve
                ( true
                , line.IsTangentTo(circle));

            r.Should().BeApproximately(0,0.0001);

            line.p2.x.Value.Should().BeApproximately(v, 0.00001);

            
        }

        /// <summary>
        /// TODO
        /// Not sure how to fix this one. My gues is that we have a local maximum due
        /// to the initial conditions
        /// </summary>
        [Test()]
        public void TangentToCircleConstraintWithLineInitiallyThroughCircleCenterShouldWork()
        {
            // Create a fully constrained circle at 0,0 with radius 1
            var circle = new Circle() { center = new Point(0, 0, false), rad = new Parameter(1, false) };

            var v = 1 / Math.Sin(Math.PI / 4);

            var line = new Line(new Point(0, -v, false, false), new Point(0, 0, true, false));

            var r = SketchSolve.Solver.solve
                ( true
                , line.IsTangentTo(circle));

            r.Should().BeApproximately(0,0.0001);

            line.p2.x.Value.Should().BeApproximately(v, 0.00001);
        }

        /// <summary>
        /// TODO
        /// Not sure how to fix this one. My gues is that we have a local maximum due
        /// to the initial conditions
        /// </summary>
        [Test()]
        public void TangentToCircleConstraintWithLineInitiallyHorizontalShouldWork()
        {
            // Create a fully constrained circle at 0,0 with radius 1
            var circle = new Circle() { center = new Point(0, 0, false), rad = new Parameter(1, false) };

            var v = 1 / Math.Sin(Math.PI / 4);

            var line = new Line(new Point(-10, -v, false, false), new Point(10, -v, false, true));

            var r = SketchSolve.Solver.solve
                ( true
                , line.IsTangentTo(circle));

            r.Should().BeApproximately(0,0.0001);

            line.p2.x.Value.Should().BeApproximately(v, 0.00001);
        }

    }
}

