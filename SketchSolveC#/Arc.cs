using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SketchSolve
{
    public class Arc : IEnumerable<Parameter>
    {
        public Parameter startAngle = new Parameter (0);
        public Parameter endAngle = new Parameter (0);
        public Parameter rad = new Parameter (0);
        public Point center = new Point (0, 0);
#region IEnumerable implementation

        public IEnumerator<Parameter> GetEnumerator ()
        {
            yield return startAngle;
            yield return endAngle;
            yield return rad;
            foreach (var p in center) {
                yield return p;
            }
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
