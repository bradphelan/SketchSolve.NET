using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FluentAssertions;
using NUnit.Framework;


namespace SketchSolve.Spec
{
    [TestFixture()]
    class VectorTest
    {
        [Test()]
        public static void DotProductShouldWork()
        {
            var v0 = new Vector(0, 1);
            v0.Dot(v0.UnitNormal).Should().Be(0);
        }
    }
}
