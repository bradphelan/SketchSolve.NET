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

        public static constraint IsVertical(this line This){
            return new constraint()
            {
                type = ConstraintEnum.vertical,
                line1 = This
            };
        }

        public static constraint IsColocated(this point This, point other){
            return new constraint()
            {
                type = ConstraintEnum.pointOnPoint,
                point1 = This,
                point2 = other
            };
        }

        public static constraint HasInternalAngle(this line This, line other, Parameter angle){
            return new constraint()
            {
                type = ConstraintEnum.internalAngle,
                line1 = This,
                line2 = other,
                parameter = angle
            };
        }

        public static constraint HasExternalAngle(this line This, line other, Parameter angle){
            return new constraint()
            {
                type = ConstraintEnum.externalAngle,
                line1 = This,
                line2 = other,
                parameter = angle
            };
        }
    }
}

