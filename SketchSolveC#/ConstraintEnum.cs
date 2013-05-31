using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SketchSolve
{
    public enum ConstraintEnum
    {
        // Done
        pointOnPoint
,
        pointOnLine
,
        // Done
        horizontal
,
        // Done
        vertical
,
        // Done 
        internalAngle
,
        // Done 
        externalAngle
,
        radiusValue
,
        tangentToArc
,
        tangentToCircle
,
        arcRules
,
        P2PDistance
,
        P2PDistanceVert
,
        P2PDistanceHorz
,
        P2LDistance
,
        P2LDistanceVert
,
        P2LDistanceHorz
,
        lineLength
,
        equalLegnth
,
        arcRadius
,
        equalRadiusArcs
,
        equalRadiusCircles
,
        equalRadiusCircArc
,
        concentricArcs
,
        concentricCircles
,
        concentricCircArc
,
        circleRadius
,
        parallel
,
        perpendicular
,
        colinear
,
        pointOnCircle
,
        pointOnArc
,
        pointOnLineMidpoint
,
        pointOnArcMidpoint
,
        pointOnCircleQuad
,
        symmetricPoints
,
        symmetricLines
,
        symmetricCircles
,
        symmetricArcs
,
    }
}
