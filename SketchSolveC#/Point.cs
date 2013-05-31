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
