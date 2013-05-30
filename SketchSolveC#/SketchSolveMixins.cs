using System;
using SketchSolve;

namespace SketchSolve
{
    public static class SketchSolveLineMixins
    {
        public static constraint IsHorizontal(this line This){
            return new constraint()
            {
                type = ConstraintEnum.horizontal,
                line1 = This
            };
        }
    }
}

