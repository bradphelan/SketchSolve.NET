using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SketchSolve
{
    public class Line : IEnumerable<Parameter>
    {
        public Point p1 = new Point (0, 0);
        public Point p2 = new Point (0, 0);

        public double dx {
            get { return p2.x.Value - p1.x.Value;}
        }

        public double dy {
            get { return p2.y.Value - p1.y.Value;}
        }

        public double dot ( Line other){
            return dx * other.dx + dy * other.dy; 
        }

        public double lengthSquared {
            get {
                return dx * dx + dy * dy;
            }
        }

        public double length {
            get {
                return Math.Sqrt (lengthSquared);
            }
        }

        // The cosine of the angle between
        // the lines
        public double cosine(Line other){
            return this.dot (other) /
                length / 
                other.length;
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
