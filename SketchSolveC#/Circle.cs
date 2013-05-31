using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SketchSolve
{
    public class Circle : IEnumerable<Parameter>
    {
        public Point center = new Point (0, 0);
        public Parameter rad = new Parameter (0);
#region IEnumerable implementation

        public IEnumerator<Parameter> GetEnumerator ()
        {
            foreach (var p in center) {
                yield return p;
            }
            yield return rad;
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
