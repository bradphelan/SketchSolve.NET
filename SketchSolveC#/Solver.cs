using System;
using System.Linq;
using System.Collections.Generic;
using System.Collections;
using Accord.Math;
using Accord.Math.Optimization;
using Accord.Math.Differentiation;

namespace SketchSolve
{

    public enum Result
    {
        error = -1
,
        succsess = 0
,
        noSolution = 1
,
    }

    public static class Solver
    {
        ///////////////////////////////////////
        /// BFGS Solver parameters
        ///////////////////////////////////////
        const double pertMag = 1e-6;
        const double pertMin = 1e-10;
        const double XconvergenceRough = 1e-8;
        const double XconvergenceFine = 1e-10;
        const double smallF = 1e-20;
        const double validSolutionFine = 1e-12;
        const double validSoltuionRough = 1e-4;
        const double rough = 0;
        const double fine = 1;
        //Note that the total number of iterations allowed is MaxIterations *xLength;
        const double MaxIterations = 50 ;


        public static double solve (bool isFine, params Constraint[] cons)
        {
            return solve (isFine, (IEnumerable<Constraint>)cons);
        }

        public static double solve (bool isFine, IEnumerable<Constraint> cons)
        {
            var constraints = cons.ToArray ();

            // Get the parameters that need solving
            Parameter[] x = constraints.SelectMany (p=>p)
                .Distinct ()
                .Where(p=>p.free==true)
                .ToArray ();

            Func<double[], double> objective = args => {
                int i = 0;
                foreach (var arg in args) {
                    x [i].Value = arg;
                    i++;
                }
                var error = Constraint.calc (constraints);

                Console.WriteLine ( "o[" +  String.Join (", ", x.Select(y=>y.Value) ) + "]");
                Console.WriteLine(error);
                return error;
            };

            var gradient 
                = new FiniteDifferences (x.Length, objective); 

            var f = new NonlinearObjectiveFunction
                ( x.Length
                , objective
                 , args => { 
                    var g =  gradient.Compute(args);
                    Console.WriteLine ( "g[" +  String.Join (", ", g ) + "]");
                    return g;
                   }
                );

            // Now we can start stating the constraints 
            var nlConstraints = x.Select (p =>
                new NonlinearConstraint(f,
                    // 1st contraint: x should be greater than or equal to 0
                    function: (args) => p.Value, 
                    shouldBe: ConstraintType.GreaterThanOrEqualTo, 
                    value: 0
                )).ToList ();

            // Finally, we create the non-linear programming solver 
            var solver = new AugmentedLagrangianSolver(x.Length, new List<NonlinearConstraint>());

            // And attempt to solve the problem 
            double minValue = solver.Minimize(f);

            return minValue;

        }

        private static Random r = new Random();
    }
}

