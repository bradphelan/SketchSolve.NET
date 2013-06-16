
module SketchSolveFS.Spec
// NOTE: If warnings appear, you may need to retarget this project to .NET 4.0. Show the Solution
// Pad, right-click on the project node, choose 'Options --> Build --> General' and change the target
// framework to .NET 4.0 or .NET 4.5.

open NUnit.Framework
open FsUnit
open SketchSolveFS.Solver

[<TestFixture>]
type ``Tunable parameters`` ()=
   
   [<Test>] 
   member test. ``Give a fuck`` ()=
       let p0 = Free (ref 20.0)
       let p1 = Fixed 30.0
       
       let l = PointConstraint (fix 1.0,free 2.0)
       
       for v in FreeParameters [|l|] do
            printfn "%O" v
            
       
   
   


