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


        /// <summary>
        /// Returns a line normal to the circle and normal
        /// to the line 
        /// </summary>
        /// <returns></returns>
        /// <param name="line"></param>
        public Line CenterTo(Line line){
            var pCircCenter = center;
            var pLineP1 = line.p1;
            var vLine = line.Vector;

            var vLineStartToCircCenter = pCircCenter - pLineP1;

            var pProjection = pLineP1 + vLineStartToCircCenter.ProjectOnto(vLine);

            return new Line (center, pProjection);
        }

        public override string ToString ()
        {
            return "[c " + center + ", r " + rad.Value + "]";
        }

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
