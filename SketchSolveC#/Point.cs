using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SketchSolve
{
    public class Point : IEnumerable<Parameter>
    {
        public Parameter x = new Parameter (0);
        public Parameter y = new Parameter (0);

        public Point (double x, double y, bool freex, bool freey)
        {
            this.x = new Parameter (x, freex);
            this.y = new Parameter (y, freey);
        }
        public Point (double x, double y, bool free=true)
        {
            this.x = new Parameter (x, free);
            this.y = new Parameter (y, free);
        }

        public override string ToString()
        {
            return this.x.Value + ";" + this.y.Value;
        }

        public static Point operator +(Point a, Vector b) {
            return new Point(a.x.Value + b.dx.Value, a.y.Value + b.dy.Value);
        }

        public static Point operator -(Point a, Vector b) {
            return new Point(a.x.Value - b.dx.Value, a.y.Value - b.dy.Value);
        }

        public static Vector operator -(Point a, Point b) {
            return new Vector(a.x.Value - b.x.Value, a.y.Value - b.y.Value);
        }

#region IEnumerable implementation
        public IEnumerator<Parameter> GetEnumerator ()
        {
            yield return x;
            yield return y;
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
