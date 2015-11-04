using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using TriangleRigging;

namespace TriangleRigging
{
    class TrianglePanel : List<List<Point>>, IDisposable
    {
        //FIELDS
        double Thickness;
        int Direction;
        bool disposed = false;

        //PROPERTIES**QUERY
        public double[] Angles { get; set; }
        public Point Centroid { get { return this[0][3]; } }
        public Point Circumcenter { get { return this[0][4]; } }
        public Point Incenter { get { return this[0][5]; } }
        public List<Point> ControlPoints { get { return this[0].GetRange(0, 3); } }
        public List<Point> PanelPoints { get { return this[1]; } }
        public List<List<Point>> ArcPoints { get { if (this.Count > 4) return this.GetRange(2, 3); else return null; } }
        public Vector[][] V { get; set; }

        //CONSTRUCTOR
        internal TrianglePanel(IEnumerable<Point> Points, double EdgeOffset, double CornerOffset, double Thickness, int Direction = 0, double MinRadius = 1/8)
        {
            this.Thickness = Thickness;
            this.Direction = Direction;
            // anchor points
            this.Add(Points.ToList());
            // side lengths of triangle opposite of point
            double[] s = { 
                             this[0][1].DistanceTo(this[0][2]), 
                             this[0][2].DistanceTo(this[0][0]), 
                             this[0][0].DistanceTo(this[0][1]) 
                         };
            // angle at point
            double[] angle = { 
                                 Math.Acos((s[1] * s[1] + s[2] * s[2] - s[0] * s[0]) / (2 * s[1] * s[2])), 
                                 Math.Acos((s[2] * s[2] + s[0] * s[0] - s[1] * s[1]) / (2 * s[2] * s[0])), 
                                 Math.Acos((s[0] * s[0] + s[1] * s[1] - s[2] * s[2]) / (2 * s[0] * s[1])) 
                             };
            Angles = angle;
            // centroid
            this[0].Add(Point.ByCoordinates(
                (this[0][0].X + this[0][1].X + this[0][2].X) / 3,
                (this[0][0].Y + this[0][1].Y + this[0][2].Y) / 3,
                (this[0][0].Z + this[0][1].Z + this[0][2].Z) / 3
                ));
            // circumcenter
            Circle c = Circle.ByThreePoints(this[0][0], this[0][1], this[0][2]);
            this[0].Add(c.CenterPoint);
            c.Dispose();
            // incenter
            this[0].Add(Point.ByCoordinates(
                (s[0] * this[0][0].X + s[1] * this[0][1].X + s[2] * this[0][2].X) / (s[0] + s[1] + s[2]),
                (s[0] * this[0][0].Y + s[1] * this[0][1].Y + s[2] * this[0][2].Y) / (s[0] + s[1] + s[2]),
                (s[0] * this[0][0].Z + s[1] * this[0][1].Z + s[2] * this[0][2].Z) / (s[0] + s[1] + s[2])
                ));
            // triangle vectors
            Vector[] v0 = {
                             Vector.ByTwoPoints(this[0][0], this[0][1]).Normalized(),
                             Vector.ByTwoPoints(this[0][0], this[0][2]).Normalized(),
                             null,
                             null
                         };
            v0[2] = v0[0].Add(v0[1]).Scale(EdgeOffset / Math.Sin(angle[0]));
            v0[3] = v0[0].Subtract(v0[1]).Normalized();
            Vector[] v1 = {
                             Vector.ByTwoPoints(this[0][1], this[0][2]).Normalized(),
                             Vector.ByTwoPoints(this[0][1], this[0][0]).Normalized(),
                             null,
                             null
                          };
            v1[2] = v1[0].Add(v1[1]).Scale(EdgeOffset / Math.Sin(angle[1]));
            v1[3] = v1[0].Subtract(v1[1]).Normalized();
            Vector[] v2 = {
                             Vector.ByTwoPoints(this[0][2], this[0][0]).Normalized(),
                             Vector.ByTwoPoints(this[0][2], this[0][1]).Normalized(),
                             null,
                             null
                          };
            v2[2] = v2[0].Add(v2[1]).Scale(EdgeOffset / Math.Sin(angle[2]));
            v2[3] = v2[0].Subtract(v2[1]).Normalized();
            Vector[][] v = { v0, v1, v2 };
            V = v;
            // panel points
            Point[] pts = {
                              this[0][0].Add(V[0][2]),
                              this[0][1].Add(V[1][2]),
                              this[0][2].Add(V[2][2])
                          };
            this.Add(pts.ToList());
            // corner offsets
            double[] r = {
                           (EdgeOffset + MinRadius) / Math.Sin(angle[0]) - MinRadius,
                           (EdgeOffset + MinRadius) / Math.Sin(angle[1]) - MinRadius,
                           (EdgeOffset + MinRadius) / Math.Sin(angle[2]) - MinRadius
                       };
            for (int i = 0; i < 3; i++) if (r[i] < CornerOffset) r[i] = CornerOffset;
            // arc at point 0
            Point[] arc0 = { null, this[0][0].Add(V[0][3].Normalized().Scale(r[0])), null };
            double r0 = (r[0] - this[0][0].DistanceTo(this[1][0])) * Math.Tan(angle[0] / 2);
            arc0[0] = arc0[1].Subtract(V[0][4].Scale(r0)).Add(V[0][1].Scale(r0));
            arc0[2] = arc0[1].Add(V[0][4].Scale(r0)).Add(V[0][0].Scale(r0));
            this.Add(arc0.ToList());
            // arc at point 1
            Point[] arc1 = { null, this[0][1].Add(V[1][3].Normalized().Scale(r[1])), null };
            double r1 = (r[1] - this[0][1].DistanceTo(this[1][1])) * Math.Tan(angle[1] / 2);
            arc1[0] = arc1[1].Subtract(V[1][4].Scale(r1)).Add(V[1][1].Scale(r1));
            arc1[2] = arc1[1].Add(V[1][4].Scale(r1)).Add(V[1][0].Scale(r1));
            this.Add(arc1.ToList());
            // arc at point 0
            Point[] arc2 = { null, this[0][2].Add(V[2][3].Normalized().Scale(r[2])), null };
            double r2 = (r[2] - this[0][2].DistanceTo(this[1][2])) * Math.Tan(angle[2] / 2);
            arc2[0] = arc2[1].Subtract(V[2][4].Scale(r2)).Add(V[2][1].Scale(r2));
            arc2[2] = arc2[1].Add(V[2][4].Scale(r2)).Add(V[2][0].Scale(r2));
            this.Add(arc2.ToList());
        }

        //METHOD**CREATE
        public static TrianglePanel ByPointsOffsetThickness(IEnumerable<Point> Points, double EdgeOffset, double CornerOffset, double Thickness)
        { return new TrianglePanel(Points, EdgeOffset, CornerOffset, Thickness); }
        public static TrianglePanel ByPointsOffsetThicknessDirectionRadius(IEnumerable<Point> Points, double EdgeOffset, double CornerOffset, double Thickness, int Direction, double MinRadius)
        { return new TrianglePanel(Points, EdgeOffset, CornerOffset, Thickness, Direction, MinRadius); }

        //METHODS**ACTION
        public PolyCurve GetPanelProfile()
        {
            Curve[] Curves = {
                          Arc.ByThreePoints(ArcPoints[0][0],ArcPoints[0][1], ArcPoints[0][2]),
                          Line.ByStartPointEndPoint(ArcPoints[0][2], ArcPoints[1][0]),
                          Arc.ByThreePoints(ArcPoints[1][0],ArcPoints[1][1], ArcPoints[1][2]),
                          Line.ByStartPointEndPoint(ArcPoints[1][2], ArcPoints[2][0]),
                          Arc.ByThreePoints(ArcPoints[2][0],ArcPoints[2][1], ArcPoints[2][2]),
                          Line.ByStartPointEndPoint(ArcPoints[2][2], ArcPoints[0][0])
                      };
            PolyCurve Profile = PolyCurve.ByJoinedCurves(Curves);
            Curves.ForEach(c => c.Dispose());
            return Profile;
        }
        public Surface GetPanelSurface()
        {
            Curve c = GetPanelProfile();
            Surface s = Surface.ByPatch(c);
            c.Dispose();
            return s;
        }
        public Solid GetPanelSolid()
        {
            Curve c = GetPanelProfile();
            Vector N = c.Normal;
            Point a = c.StartPoint;
            Point b = c.StartPoint;
            if (Direction == 0)
            {
                a = a.Add(N.Scale(-Thickness / 2));
                b = b.Add(N.Scale(Thickness / 2));
            }
            else b = b.Add(N.Scale(Direction * Thickness));
            Line l = Line.ByStartPointEndPoint(a,b);
            Solid s = c.SweepAsSolid(l);
            c.Dispose();
            N.Dispose();
            a.Dispose();
            b.Dispose();
            l.Dispose();
            return s;
        }
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
   
        //METHODS**PRIVATE
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                V.ForEach(v => v.Dispose());
            }
            disposed = true;
        }
    }



}
