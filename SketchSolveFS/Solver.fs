module SketchSolveFS.Solver


type Vector = {X:double; Y:double}
    with
    // Dot product
    static member (*) (v1 : Vector, v2 : Vector) = v1.X * v2.X + v1.Y * v2.Y
    // Scalar product
    static member (*) (s : double, v : Vector) = {X =v.X * s; Y=v.Y * s} 
    static member (*) (v : Vector, s : double) = {X =v.X * s; Y=v.Y * s} 
    // Negation
    static member inline (~-) (v : Vector) = -1.0 * v
    // Addition
    static member inline (+) (v1 : Vector, v2 : Vector) = {X = v1.X + v2.X; Y = v1.Y + v2.Y}
    // Subtraction
    static member inline (-) (v1 : Vector, v2 : Vector) = (-v1) + v2
    member m.LengthSquared with get() = m * m
    member m.Length with get() = sqrt m.LengthSquared
    member m.Unit 
        with get() =
            let len = m.Length
            { X = m.X/len; Y = m.Y/len }
    member m.Normal with get() = { X = -m.Y;  Y = m.X }
    
    // Project vector onto vector
    member m.ProjectOntoVector(v:Vector) = m * v.Unit * v.Unit 
    
    // Project  onto line
    member m.ProjectOntoLine(l:Line) = 
        (m - l.P0).ProjectOntoVector(l.Vector)  + l.P0
        
    // Project a point onto the circle radius
    member m.ProjectOntoCircle(c:Circle) =
        (m - c.Center).Unit * c.Radius + c.Center
        

and Line = {P0: Vector; P1: Vector }
    with 
    // Vector orientation from P0 to P1
    member m.Vector with get() = m.P1 - m.P0 
    // Shortest line between point and line
    member m.ShortestConnectingLine(v:Vector) = {P0= v.ProjectOntoLine(m); P1 = v}

and Circle = {Center: Vector; Radius: double} 
    with
    member m.ShortestConnectingLine(v:Vector) = {P0= v.ProjectOntoCircle(m); P1 = v}
        

type Parameter =
    | Fixed of double
    | Free of double ref
    with 
    member m.Val 
        with get() =
            match m with 
            | Fixed v -> v
            | Free v -> !v
    override m.ToString() = 
        match m with
        | Fixed v -> sprintf "%f" v
        | Free  v -> sprintf "$%f" !v

let free v = Free (ref v)     
let fix v = Fixed v
    
type IParameterizable =
    abstract member Parameters : unit -> seq<Parameter>
    
type VectorConstraint(X:Parameter, Y:Parameter) =
    override  m.ToString() = sprintf "%A %A" X Y
    member m.Vector with get() = Vector X.Val Y.Val
          
    interface IParameterizable with
        member this.Parameters() = seq {
                yield X
                yield Y
            }
            
type LineConstraint(A: VectorConstraint, B: VectorConstraint) =
    override  m.ToString() = sprintf "%A -> %A" A B
    member m.Line with get() = line A.Vector B.Vector 
    interface IParameterizable with
        member this.Parameters() = 
            seq { 
                yield! (A :> IParameterizable).Parameters() 
                yield! (B :> IParameterizable).Parameters() 
            }
            

type CircleConstraint(Center: VectorConstraint, Radius:Parameter) =
    override  m.ToString() = sprintf "%A -> %A" Center Radius 
    member m.Circle with get() = circle Center.Vector Radius.Val
    interface IParameterizable with
        member this.Parameters() = 
            seq { 
                yield! (Center :> IParameterizable).Parameters() 
                yield Radius 
            }
            
type Constraints  =
    | LineLineParallel of LineConstraint * LineConstraint
    | VectorVectorSame of VectorConstraint * VectorConstraint
    | CircleLineTangent of CircleConstraint * LineConstraint
    
let errorFn = function
    | LineLineParallel(l0, l1) -> fun ()->0.0
    | VectorVectorSame(p0,p1) -> fun()-> 0.0
    | CircleLineTangent(c0,l1) -> fun() -> 0.0
            
let FreeParameters (shapes:seq<IParameterizable>) = 
    shapes 
    |> Seq.collect (fun shape -> shape.Parameters())
    |> Seq.filter ( function Free _ -> true | _ -> false )
    
