using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SketchSolve
{
    public class Constraint : IEnumerable<Parameter>
    {
        public ConstraintEnum type;
        public Point point1;
        public Point point2;
        public Line line1;
        public Line line2;
        public Line SymLine;
        public Circle circle1;
        public Circle circle2;
        public Arc arc1;
        public Arc arc2;
        public Parameter parameter = null;
        //radius, length, angle etc...

#region IEnumerable implementation

        public IEnumerator<Parameter> GetEnumerator ()
        {
            List<IEnumerable<Parameter>> list = new List<IEnumerable<Parameter>> () {
                point1,
                point2,
                line1,
                line2,
                SymLine,
                circle1,
                circle2,
                arc1,
                arc2,
                new []{parameter}
            };
            return list
                .Where (p=>p!=null)
                .SelectMany (p=>p)
                .Where (p=>p!=null)
                .GetEnumerator ();
        }
#endregion

#region IEnumerable implementation

        IEnumerator IEnumerable.GetEnumerator ()
        {
            return this.GetEnumerator ();
        }
#endregion
    };
}
