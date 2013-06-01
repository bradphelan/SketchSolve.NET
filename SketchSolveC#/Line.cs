using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SketchSolve
{
    public class Vector : IEnumerable<Parameter>
    {
        public Parameter dx = new Parameter (0);
        public Parameter dy = new Parameter (0);

        public Vector (double dx, double dy, bool freex, bool freey)
        {
            this.dx = new Parameter (dx, freex);
            this.dy = new Parameter (dy, freey);
        }

        public Vector (double dx, double dy, bool free=true)
        {
            this.dx = new Parameter (dx, free);
            this.dy = new Parameter (dy, free);
        }

        public override string ToString()
        {
            return "-> " + dx.Value + ";" + dy.Value;
        }

        public double LengthSquared {
            get {
                return dx.Value * dx.Value + dy.Value * dy.Value;
            }
        }

        public double Length {
            get {
                return Math.Sqrt (LengthSquared);
            }
        }

        // The cosine of the angle between
        // the lines
        public double Cosine(Vector other){
            return this.Dot (other) /
                Length / 
                other.Length;
        }

        public double Dot ( Vector other){
            return dx.Value * other.dx.Value 
                + dy.Value * other.dy.Value; 
        }

        public Vector ProjectOnto(Vector other)
        {
            var unit = other.Unit;
            return this.Dot(unit) * unit;
        }

        #region basic operators
        public static Vector operator *(Vector a, double b) {
            return new Vector(a.dx.Value * b, a.dy.Value * b, false);
        }
        public static Vector operator *(double b , Vector a)
        {
            return a * b;
        }

        public static Vector operator /(Vector a, double b)
        {
            return a * (1 / b);
        }
        #endregion

        public Vector Unit
        {
            get
            {
                var l = Length;
                return new Vector(dx.Value / l, dy.Value / l);
            }
        }

        public Vector UnitNormal
        {
            get
            {
                var l = Length;
                return new Vector(-dy.Value / l, dx.Value / l);
            }
        }


#region IEnumerable implementation
        public IEnumerator<Parameter> GetEnumerator ()
        {
            yield return dx;
            yield return dy;
        }
#endregion
#region IEnumerable implementation
        IEnumerator IEnumerable.GetEnumerator ()
        {
            return this.GetEnumerator ();
        }
#endregion
    }

    public class Line : IEnumerable<Parameter>
    {
        public Point p1 = null;
        public Point p2 = null;
        public Vector v1 = null;

        public Line(Point p1, Point p2)
        {
            this.p1 = p1;
            this.p2 = p2;
        }

        public Line(Point p1, Vector v)
        {
            this.p1 = p1;
            this.v1 = v;
        }

        public override string ToString()
        {
            return "l " + p1 + " : " + p2;
        }

        public Vector Vector
        {
            get
            {
                if (v1==null)
                {
                    return new Vector(dx, dy, false, false);
                }
                else
                {
                    return v1;
                }

            }
        }

        private double dx {
            get { return p2.x.Value - p1.x.Value;}
        }

        private double dy {
            get { return p2.y.Value - p1.y.Value;}
        }



#region IEnumerable implementation
        public IEnumerator<Parameter> GetEnumerator ()
        {
            return new [] { p1, p2 }.SelectMany (p=>p).GetEnumerator ();
        }
#endregion
#region IEnumerable implementation
        IEnumerator IEnumerable.GetEnumerator ()
        {
            return this.GetEnumerator ();
        }
#endregion
    }
}
