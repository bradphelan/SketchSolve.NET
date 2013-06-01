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

        static double _hypot (double a, double b)
        {
            return Math.Sqrt (a*a+b*b);
        }


        public static double calc (Constraint [] cons)
        {
            int consLength = cons.Length;
            double error = 0;
            double temp, dx, dy, m, n, Ex, Ey, rad1, rad2, t, Xint, Yint, dx2, dy2, hyp1, hyp2, temp2;




            for (int i=0; i<consLength; i++) {


                // Crappy hack but it will get us going
                double P1_x = cons [i].point1 == null ? 0 : cons [i].point1.x.Value;
                double P1_y = cons [i].point1 == null ? 0 : cons [i].point1.y.Value;
                double P2_x = cons [i].point2 == null ? 0 : cons [i].point2.x.Value;
                double P2_y = cons [i].point2 == null ? 0 : cons [i].point2.y.Value;
                double L1_P1_x = cons [i].line1 == null ? 0 : cons [i].line1.p1.x.Value;
                double L1_P1_y = cons [i].line1 == null ? 0 : cons [i].line1.p1.y.Value;
                double L1_P2_x = cons [i].line1 == null ? 0 : cons [i].line1.p2.x.Value;
                double L1_P2_y = cons [i].line1 == null ? 0 : cons [i].line1.p2.y.Value;
                double L2_P1_x = cons [i].line2 == null ? 0 : cons [i].line2.p1.x.Value;
                double L2_P1_y = cons [i].line2 == null ? 0 : cons [i].line2.p1.y.Value;
                double L2_P2_x = cons [i].line2 == null ? 0 : cons [i].line2.p2.x.Value;
                double L2_P2_y = cons [i].line2 == null ? 0 : cons [i].line2.p2.y.Value;
                double C1_Center_x = cons [i].circle1 == null ? 0 : cons [i].circle1.center.x.Value;
                double C1_Center_y = cons [i].circle1 == null ? 0 : cons [i].circle1.center.y.Value;
                double C1_rad = cons [i].circle1 == null ? 0 : cons [i].circle1.rad.Value;
                double C2_Center_x = cons [i].circle2 == null ? 0 : cons [i].circle2.center.x.Value;
                double C2_Center_y = cons [i].circle2 == null ? 0 : cons [i].circle2.center.y.Value;
                double C2_rad = cons [i].circle2 == null ? 0 : cons [i].circle2.rad.Value;

                double A1_startA = cons [i].arc1 == null ? 0 : cons [i].arc1.startAngle.Value;
                double A1_endA = cons [i].arc1 == null ? 0 : cons [i].arc1.endAngle.Value;
                double A1_radius = cons [i].arc1 == null ? 0 : cons [i].arc1.rad.Value;
                double A1_Center_x = cons [i].arc1 == null ? 0 : cons [i].arc1.center.x.Value;
                double A1_Center_y = cons [i].arc1 == null ? 0 : cons [i].arc1.center.y.Value;
                double A2_startA = cons [i].arc2 == null ? 0 : cons [i].arc2.startAngle.Value;
                double A2_endA = cons [i].arc2 == null ? 0 : cons [i].arc2.endAngle.Value;
                double A2_radius = cons [i].arc2 == null ? 0 : cons [i].arc2.rad.Value;
                double A2_Center_x = cons [i].arc2 == null ? 0 : cons [i].arc2.center.x.Value;
                double A2_Center_y = cons [i].arc2 == null ? 0 : cons [i].arc2.center.y.Value;

                double A1_Start_x = (A1_Center_x + A1_radius * Math.Cos (A1_startA));
                double A1_Start_y = (A1_Center_y + A1_radius * Math.Sin (A1_startA));
                double A1_End_x = (A1_Center_x + A1_radius * Math.Cos (A1_endA));
                double A1_End_y = (A1_Center_y + A1_radius * Math.Sin (A1_endA));
                double A2_Start_x = (A1_Center_x + A2_radius * Math.Cos (A2_startA));
                double A2_Start_y = (A1_Center_y + A2_radius * Math.Sin (A2_startA));
                double A2_End_x = (A1_Center_x + A2_radius * Math.Cos (A2_endA));
                double A2_End_y = (A1_Center_y + A2_radius * Math.Sin (A2_endA));


                double length = cons [i].parameter == null ? 0 : cons [i].parameter.Value;
                double distance = length;
                double radius = length;
                double angleP = length;
                double quadIndex = length;

                double Sym_P1_x = cons [i].SymLine == null ? 0 : cons [i].SymLine.p1.x.Value;
                double Sym_P1_y = cons [i].SymLine == null ? 0 : cons [i].SymLine.p1.y.Value;

                double Sym_P2_x = cons [i].SymLine == null ? 0 : cons [i].SymLine.p2.x.Value;
                double Sym_P2_y = cons [i].SymLine == null ? 0 : cons [i].SymLine.p2.y.Value;



                if ((cons [i]).type == ConstraintEnum.pointOnPoint) {
                    //Hopefully avoid this constraint, make coincident points use the same parameters
                    var l2 = (cons [i].point1 - cons [i].point2).LengthSquared;
                    error += l2;
                }


                if (cons [i].type == ConstraintEnum.P2PDistance) {
                    error += (P1_x - P2_x) * (P1_x - P2_x) + (P1_y - P2_y) * (P1_y - P2_y) - distance * distance;

                }

                if (cons [i].type == ConstraintEnum.P2PDistanceVert) {
                    error += (P1_y - P2_y) * (P1_y - P2_y) - distance * distance;
                }

                if (cons [i].type == ConstraintEnum.P2PDistanceHorz) {
                    error += (P1_x - P2_x) * (P1_x - P2_x) - distance * distance;
                }


                if ((cons [i]).type == ConstraintEnum.pointOnLine) {
                    dx = L1_P2_x - L1_P1_x;
                    dy = L1_P2_y - L1_P1_y;

                    m = dy / dx; //Slope
                    n = dx / dy; //1/Slope

                    if (m <= 1 && m >= -1) {
                        //Calculate the expected y point given the x coordinate of the point
                        Ey = L1_P1_y + m * (P1_x - L1_P1_x);
                        error += (Ey - P1_y) * (Ey - P1_y);
                    } else {
                        //Calculate the expected x point given the y coordinate of the point
                        Ex = L1_P1_x + n * (P1_y - L1_P1_y);
                        error += (Ex - P1_x) * (Ex - P1_x);
                    }
                }

                if ((cons [i]).type == ConstraintEnum.P2LDistance) {
                    dx = L1_P2_x - L1_P1_x;
                    dy = L1_P2_y - L1_P1_y;

                    t = -(L1_P1_x * dx - P1_x * dx + L1_P1_y * dy - P1_y * dy) / (dx * dx + dy * dy);
                    Xint = L1_P1_x + dx * t;
                    Yint = L1_P1_y + dy * t;
                    temp = _hypot ((P1_x - Xint), (P1_y - Yint)) - distance;
                    error += temp * temp / 10;

                }
                if ((cons [i]).type == ConstraintEnum.P2LDistanceVert) {
                    dx = L1_P2_x - L1_P1_x;
                    dy = L1_P2_y - L1_P1_y;

                    t = (P1_x - L1_P1_x) / dx;
                    Yint = L1_P1_y + dy * t;
                    temp = Math.Abs ((P1_y - Yint)) - distance;
                    error += temp * temp;

                }
                if ((cons [i]).type == ConstraintEnum.P2LDistanceHorz) {
                    dx = L1_P2_x - L1_P1_x;
                    dy = L1_P2_y - L1_P1_y;

                    t = (P1_y - L1_P1_y) / dy;
                    Xint = L1_P1_x + dx * t;
                    temp = Math.Abs ((P1_x - Xint)) - distance;
                    error += temp * temp / 10;

                }


                if (cons [i].type == ConstraintEnum.vertical) {
                    double odx = L1_P2_x - L1_P1_x;
                    error += odx * odx;
                }

                if (cons [i].type == ConstraintEnum.horizontal) {
                    double ody = L1_P2_y - L1_P1_y;
                    error += ody * ody;
                }

                if (cons [i].type == ConstraintEnum.tangentToCircle) {

                    var l = cons [i].line1;
                    var c = cons [i].circle1;
                    temp = c.CenterTo (l).Vector.Length - c.rad.Value;
                    error += temp * temp;

                }

                if (cons [i].type == ConstraintEnum.tangentToArc) {
                    /*
               double dx,dy,Rpx,Rpy,RpxN,RpyN,hyp,error1,error2,rad;
               dx = L1_P2_x - L1_P1_x;
               dy = L1_P2_y - L1_P1_y;

               hyp=_hypot(dx,dy);

               double u = (A1_Center_x - L1_P1_x) * (L1_P2_x - L1_P1_x) + (A1_Center_y - L1_P1_y) * (L1_P2_y - L1_P1_y);
               u/=hyp*hyp;

               double x = L1_P1_x + u *(L1_P2_x - L1_P1_x);
               double y = L1_P1_y + u *(L1_P2_y - L1_P1_y);

               double dcsx = A1_Center_x - A1_Start_x;
               double dcsy = A1_Center_y - A1_Start_y;
               double dcex = A1_Center_x - A1_End_x;
               double dcey = A1_Center_y - A1_End_y;
               rad=(dcsx*dcsx + dcsy * dcsy);
            //  rad+=(dcex*dcex + dcey * dcey)/4;

            double dcx = A1_Center_x-x;
            double dcy = A1_Center_y-y;
            temp = (dcx * dcx + dcy * dcy) - rad;
            error += temp*temp*100;
            */

                    //#if defined(NEWARC)
                    dx = L1_P2_x - L1_P1_x;
                    dy = L1_P2_y - L1_P1_y;


                    double radsq = (A1_Center_x - A1_Start_x) * (A1_Center_x - A1_Start_x) + (A1_Center_y - A1_Start_y) * (A1_Center_y - A1_Start_y);
                    t = -(L1_P1_x * dx - A1_Center_x * dx + L1_P1_y * dy - A1_Center_y * dy) / (dx * dx + dy * dy);
                    Xint = L1_P1_x + dx * t;
                    Yint = L1_P1_y + dy * t;
                    temp = (A1_Center_x - Xint) * (A1_Center_x - Xint) + (A1_Center_y - Yint) * (A1_Center_y - Yint) - radsq;
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.arcRules) {
                    //rad1=_hypot(A1_Center_x - A1_Start_x , A1_Center_y - A1_Start_y);
                    //rad2=_hypot(A1_Center_x - A1_End_x , A1_Center_y - A1_End_y);
                    //error += (rad1-rad2)*(rad1-rad2);
                    //double dx,dy,Rpx,Rpy,RpxN,RpyN,hyp,error1,error2,rad;
                    //dx = A1_End_x - A1_Start_x;
                    //dy = A1_End_y - A1_Start_y;

                    //hyp=_hypot(dx,dy);

                    //double u = (A1_Center_x - A1_Start_x) * (A1_End_x - A1_Start_x) + (A1_Center_y - A1_Start_y) * (A1_End_y - A1_Start_y);
                    //u/=hyp*hyp;

                    //temp = Math.Sin(u - .5);
                    //error+=temp*temp*temp*temp*100000;
                    //error+=Math.Pow(-2*A1_Center_x*A1_End_y - 2*A1_Center_y*A1_End_y + A1_End_x*A1_End_y + Math.Pow(A1_End_y,2) + 2*A1_Center_x*A1_Start_x - 2*A1_Center_y*A1_Start_x - A1_End_x*A1_Start_x + 4*A1_End_y*A1_Start_x - 3*Math.Pow(A1_Start_x,2) +  2*A1_Center_y*A1_Start_y + A1_Start_x*A1_Start_y - Math.Pow(A1_Start_y,2),2)/(8*Math.Pow(A1_End_y,2) + 8*Math.Pow(A1_Start_x,2) - 8*A1_End_y*A1_Start_y -  8*A1_Start_x*A1_Start_y + 4*Math.Pow(A1_Start_y,2));
                    double a1endx2 = A1_End_x * A1_End_x;
                    double a1endy2 = A1_End_y * A1_End_y;
                    double a1startx2 = A1_Start_x * A1_Start_x;
                    double a1starty2 = A1_Start_y * A1_Start_y;
                    double num = -2 * A1_Center_x * A1_End_x + a1endx2 - 2 * A1_Center_y * A1_End_y + a1endy2 + 2 * A1_Center_x * A1_Start_x - a1startx2 + 2 * A1_Center_y * A1_Start_y - a1starty2;
                    error += num * num / (4.0 * a1endx2 + a1endy2 - 2 * A1_End_x * A1_Start_x + a1startx2 - 2 * A1_End_y * A1_Start_y + a1starty2);

                }

                if (cons [i].type == ConstraintEnum.lineLength) {
                    temp = Math.Sqrt (Math.Pow(L1_P2_x - L1_P1_x,2) + Math.Pow(L1_P2_y - L1_P1_y,2)) - length;
                    //temp=_hypot(L1_P2_x - L1_P1_x , L1_P2_y - L1_P1_y) - length;
                    error += temp * temp * 100;
                }

                if (cons [i].type == ConstraintEnum.equalLegnth) {
                    temp = _hypot (L1_P2_x - L1_P1_x, L1_P2_y - L1_P1_y) - _hypot (L2_P2_x - L2_P1_x, L2_P2_y - L2_P1_y);
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.arcRadius) {
                    rad1 = _hypot (A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
                    rad2 = _hypot (A1_Center_x - A1_End_x, A1_Center_y - A1_End_y);
                    temp = rad1 - radius;
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.equalRadiusArcs) {
                    rad1 = _hypot (A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
                    rad2 = _hypot (A2_Center_x - A2_Start_x, A2_Center_y - A2_Start_y);
                    temp = rad1 - rad2;
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.equalRadiusCircles) {
                    temp = C1_rad - C2_rad;
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.equalRadiusCircArc) {
                    rad1 = _hypot (A1_Center_x - A1_Start_x, A1_Center_y - A1_Start_y);
                    temp = rad1 - C1_rad;
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.concentricArcs) {
                    temp = _hypot (A1_Center_x - A2_Center_x, A1_Center_y - A2_Center_y);
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.concentricCircles) {
                    temp = _hypot (C1_Center_x - C2_Center_x, C1_Center_y - C2_Center_y);
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.concentricCircArc) {
                    temp = _hypot (A1_Center_x - C1_Center_x, A1_Center_y - C1_Center_y);
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.circleRadius) {
                    error += (C1_rad - radius) * (C1_rad - radius);
                }
                if (cons [i].type == ConstraintEnum.internalAngle) {
                    temp = cons [i].line1.Vector.Cosine (cons[i].line2.Vector);

                    temp2 = Math.Cos (angleP);
                    error += (temp - temp2) * (temp - temp2);
                }

                if (cons [i].type == ConstraintEnum.externalAngle) {
                    temp = cons [i].line1.Vector.Cosine (cons[i].line2.Vector);
                    temp2 = Math.Cos (Math.PI-angleP);
                    error += (temp - temp2) * (temp - temp2);
                }

                if (cons [i].type == ConstraintEnum.perpendicular) {
                    temp = cons[i].line1.Vector.Dot(cons[i].line2.Vector);
                    error += temp * temp;
                }

                if (cons [i].type == ConstraintEnum.parallel) {
                    dx = L1_P2_x - L1_P1_x;
                    dy = L1_P2_y - L1_P1_y;
                    dx2 = L2_P2_x - L2_P1_x;
                    dy2 = L2_P2_y - L2_P1_y;

                    hyp1 = _hypot (dx, dy);
                    hyp2 = _hypot (dx2, dy2);

                    dx = dx / hyp1;
                    dy = dy / hyp1;
                    dx2 = dx2 / hyp2;
                    dy2 = dy2 / hyp2;

                    temp = dy * dx2 - dx * dy2;
                    error += (temp) * (temp);
                }
                // Colinear constraint
                if (cons [i].type == ConstraintEnum.colinear) {
                    dx = L1_P2_x - L1_P1_x;
                    dy = L1_P2_y - L1_P1_y;

                    m = dy / dx;
                    n = dx / dy;
                    // Calculate the error between the expected intersection point
                    // and the true point of the second lines two end points on the
                    // first line
                    if (m <= 1 && m > -1) {
                        //Calculate the expected y point given the x coordinate of the point
                        Ey = L1_P1_y + m * (L2_P1_x - L1_P1_x);
                        error += (Ey - L2_P1_y) * (Ey - L2_P1_y);

                        Ey = L1_P1_y + m * (L2_P2_x - L1_P1_x);
                        error += (Ey - L2_P2_y) * (Ey - L2_P2_y);
                    } else {
                        //Calculate the expected x point given the y coordinate of the point
                        Ex = L1_P1_x + n * (L2_P1_y - L1_P1_y);
                        error += (Ex - L2_P1_x) * (Ex - L2_P1_x);

                        Ex = L1_P1_x + n * (L2_P2_y - L1_P1_y);
                        error += (Ex - L2_P2_x) * (Ex - L2_P2_x);
                    }
                }
                // Point on a circle
                if (cons [i].type == ConstraintEnum.pointOnCircle) {
                    //see what the current radius to the point is
                    rad1 = _hypot (C1_Center_x-P1_x, C1_Center_y - P1_y);
                    //Compare this radius to the radius of the circle, return the error squared
                    temp = rad1 - C1_rad;
                    error += temp * temp;
                    //cout<<"Point On circle error"<<temp*temp<<endl;
                }
                if (cons [i].type == ConstraintEnum.pointOnArc) {
                    //see what the current radius to the point is
                    rad1 = _hypot (A1_Center_x-P1_x, A1_Center_y - P1_y);
                    rad2 = _hypot (A1_Center_x-A1_Start_x, A1_Center_y - A1_Start_y);
                    //Compare this radius to the radius of the circle, return the error squared
                    temp = rad1 - rad2;
                    error += temp * temp;
                    //cout<<"Point On circle error"<<temp*temp<<endl;
                }
                if (cons [i].type == ConstraintEnum.pointOnLineMidpoint) {
                    Ex = (L1_P1_x + L1_P2_x) / 2;
                    Ey = (L1_P1_y + L1_P2_y) / 2;
                    temp = Ex - P1_x;
                    temp2 = Ey - P1_y;
                    error += temp * temp + temp2 * temp2;
                }
                if (cons [i].type == ConstraintEnum.pointOnArcMidpoint) {
                    rad1 = _hypot (A1_Center_x-A1_Start_x, A1_Center_y - A1_Start_y);
                    temp = Math.Atan2 (A1_Start_y-A1_Center_y, A1_Start_x - A1_Center_x);
                    temp2 = Math.Atan2 (A1_End_y-A1_Center_y, A1_End_x - A1_Center_x);
                    Ex = A1_Center_x + rad1 * Math.Cos ((temp2+temp)/2);
                    Ey = A1_Center_y + rad1 * Math.Sin ((temp2+temp)/2);
                    temp = (Ex - P1_x);
                    temp2 = (Ey - P1_y);
                    error += temp * temp + temp2 * temp2;
                }

                if (cons [i].type == ConstraintEnum.pointOnCircleQuad) {
                    Ex = C1_Center_x;
                    Ey = C1_Center_y;
                    switch ((int)quadIndex) {
                        case 0:
                        Ex += C1_rad;
                        break;
                        case 1:
                        Ey += C1_rad;
                        break;
                        case 2:
                        Ex -= C1_rad;
                        break;
                        case 3:
                        Ey -= C1_rad;
                        break;
                    }
                    temp = (Ex - P1_x);
                    temp2 = (Ey - P1_y);
                    error += temp * temp + temp2 * temp2;
                }
                if (cons [i].type == ConstraintEnum.symmetricPoints) {
                    dx = Sym_P2_x - Sym_P1_x;
                    dy = Sym_P2_y - Sym_P1_y;
                    t = -(dy * P1_x - dx * P1_y - dy * Sym_P1_x + dx * Sym_P1_y) / (dx * dx + dy * dy);
                    Ex = P1_x + dy * t * 2;
                    Ey = P1_y - dx * t * 2;
                    temp = (Ex - P2_x);
                    temp2 = (Ey - P2_y);
                    error += temp * temp + temp2 * temp2;
                }
                if (cons [i].type == ConstraintEnum.symmetricLines) {
                    dx = Sym_P2_x - Sym_P1_x;
                    dy = Sym_P2_y - Sym_P1_y;
                    t = -(dy * L1_P1_x - dx * L1_P1_y - dy * Sym_P1_x + dx * Sym_P1_y) / (dx * dx + dy * dy);
                    Ex = L1_P1_x + dy * t * 2;
                    Ey = L1_P1_y - dx * t * 2;
                    temp = (Ex - L2_P1_x);
                    temp2 = (Ey - L2_P1_y);
                    error += temp * temp + temp2 * temp2;
                    t = -(dy * L1_P2_x - dx * L1_P2_y - dy * Sym_P1_x + dx * Sym_P1_y) / (dx * dx + dy * dy);
                    Ex = L1_P2_x + dy * t * 2;
                    Ey = L1_P2_y - dx * t * 2;
                    temp = (Ex - L2_P2_x);
                    temp2 = (Ey - L2_P2_y);
                    error += temp * temp + temp2 * temp2;
                }
                if (cons [i].type == ConstraintEnum.symmetricCircles) {
                    dx = Sym_P2_x - Sym_P1_x;
                    dy = Sym_P2_y - Sym_P1_y;
                    t = -(dy * C1_Center_x - dx * C1_Center_y - dy * Sym_P1_x + dx * Sym_P1_y) / (dx * dx + dy * dy);
                    Ex = C1_Center_x + dy * t * 2;
                    Ey = C1_Center_y - dx * t * 2;
                    temp = (Ex - C2_Center_x);
                    temp2 = (Ey - C2_Center_y);
                    error += temp * temp + temp2 * temp2;
                    temp = (C1_rad - C2_rad);
                    error += temp * temp;
                }
                if (cons [i].type == ConstraintEnum.symmetricArcs) {
                    dx = Sym_P2_x - Sym_P1_x;
                    dy = Sym_P2_y - Sym_P1_y;
                    t = -(dy * A1_Start_x - dx * A1_Start_y - dy * Sym_P1_x + dx * Sym_P1_y) / (dx * dx + dy * dy);
                    Ex = A1_Start_x + dy * t * 2;
                    Ey = A1_Start_y - dx * t * 2;
                    temp = (Ex - A2_Start_x);
                    temp2 = (Ey - A2_Start_y);
                    error += temp * temp + temp2 * temp2;
                    t = -(dy * A1_End_x - dx * A1_End_y - dy * Sym_P1_x + dx * Sym_P1_y) / (dx * dx + dy * dy);
                    Ex = A1_End_x + dy * t * 2;
                    Ey = A1_End_y - dx * t * 2;
                    temp = (Ex - A2_End_x);
                    temp2 = (Ey - A2_End_y);
                    error += temp * temp + temp2 * temp2;
                    t = -(dy * A1_Center_x - dx * A1_Center_y - dy * Sym_P1_x + dx * Sym_P1_y) / (dx * dx + dy * dy);
                    Ex = A1_Center_x + dy * t * 2;
                    Ey = A1_Center_y - dx * t * 2;
                    temp = (Ex - A2_Center_x);
                    temp2 = (Ey - A2_Center_y);
                    error += temp * temp + temp2 * temp2;
                }
            }

            // Prevent symetry errors
            return error; 

        }

    };
}
