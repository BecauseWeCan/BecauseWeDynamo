using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using Topology;

namespace Panelization
{
    public class TrianglePanel : IDisposable
    {
        //FIELDS
        double Thickness;
        int Direction;
        bool disposed = false;

        //PROPERTIES**QUERY
        public List<double> Angles { get; set; }
        public Point Centroid { get; private set; }
        public Point Circumcenter { get; private set; }
        public Point Incenter { get; private set; }
        public List<Point> ControlPoints { get; private set; }
        public List<Point> PanelPoints { get; private set; }
        public List<List<Point>> ArcPoints { get; private set; }
        public Vector[][] EdgeVectors { get; set; }
        public Vector Normal { get; set; }

        //CONSTRUCTOR
        internal TrianglePanel(IEnumerable<Point> Points, double EdgeOffset, double CornerOffset, double Thickness, int Direction = 0, double MinRadius = 1/8)
        {
            this.Thickness = Thickness;
            this.Direction = Direction;
            // anchor points
            ControlPoints = new List<Point>(Points.ToList());
            // side lengths of triangle opposite of point
            double[] s = { 
                             ControlPoints[1].DistanceTo(ControlPoints[2]), 
                             ControlPoints[2].DistanceTo(ControlPoints[0]), 
                             ControlPoints[0].DistanceTo(ControlPoints[1]) 
                         };
            // angle at point
            Angles = new List<double>{ 
                                 Math.Acos((s[1] * s[1] + s[2] * s[2] - s[0] * s[0]) / (2 * s[1] * s[2])), 
                                 Math.Acos((s[2] * s[2] + s[0] * s[0] - s[1] * s[1]) / (2 * s[2] * s[0])), 
                                 Math.Acos((s[0] * s[0] + s[1] * s[1] - s[2] * s[2]) / (2 * s[0] * s[1])) 
                             };
            // centroid
            Centroid = Point.ByCoordinates(
                (ControlPoints[0].X + ControlPoints[1].X + ControlPoints[2].X) / 3,
                (ControlPoints[0].Y + ControlPoints[1].Y + ControlPoints[2].Y) / 3,
                (ControlPoints[0].Z + ControlPoints[1].Z + ControlPoints[2].Z) / 3
                );
            // circumcenter
            Circle c = Circle.ByThreePoints(ControlPoints[0], ControlPoints[1], ControlPoints[2]);
            Normal = c.Normal;
            Circumcenter = c.CenterPoint;
            c.Dispose();
            // incenter
            Incenter = Point.ByCoordinates(
                (s[0] * ControlPoints[0].X + s[1] * ControlPoints[1].X + s[2] * ControlPoints[2].X) / (s[0] + s[1] + s[2]),
                (s[0] * ControlPoints[0].Y + s[1] * ControlPoints[1].Y + s[2] * ControlPoints[2].Y) / (s[0] + s[1] + s[2]),
                (s[0] * ControlPoints[0].Z + s[1] * ControlPoints[1].Z + s[2] * ControlPoints[2].Z) / (s[0] + s[1] + s[2])
                );
            // triangle vectors
            Vector[] v0 = {
                             Vector.ByTwoPoints(ControlPoints[0], ControlPoints[1]).Normalized(),
                             Vector.ByTwoPoints(ControlPoints[0], ControlPoints[2]).Normalized(),
                             null,
                             null
                         };
            v0[2] = v0[0].Add(v0[1]).Scale(EdgeOffset / Math.Sin(Angles[0]));
            v0[3] = v0[0].Subtract(v0[1]).Normalized();
            Vector[] v1 = {
                             Vector.ByTwoPoints(ControlPoints[1], ControlPoints[2]).Normalized(),
                             Vector.ByTwoPoints(ControlPoints[1], ControlPoints[0]).Normalized(),
                             null,
                             null
                          };
            v1[2] = v1[0].Add(v1[1]).Scale(EdgeOffset / Math.Sin(Angles[1]));
            v1[3] = v1[0].Subtract(v1[1]).Normalized();
            Vector[] v2 = {
                             Vector.ByTwoPoints(ControlPoints[2], ControlPoints[0]).Normalized(),
                             Vector.ByTwoPoints(ControlPoints[2], ControlPoints[1]).Normalized(),
                             null,
                             null
                          };
            v2[2] = v2[0].Add(v2[1]).Scale(EdgeOffset / Math.Sin(Angles[2]));
            v2[3] = v2[0].Subtract(v2[1]).Normalized();
            Vector[][] v = { v0, v1, v2 };
            EdgeVectors = v;
            // panel points
            PanelPoints = new List<Point> {
                              ControlPoints[0].Add(EdgeVectors[0][2]),
                              ControlPoints[1].Add(EdgeVectors[1][2]),
                              ControlPoints[2].Add(EdgeVectors[2][2])
                          };
            // corner offsets
            double[] r = {
                           (EdgeOffset + MinRadius) / Math.Sin(Angles[0]) - MinRadius,
                           (EdgeOffset + MinRadius) / Math.Sin(Angles[1]) - MinRadius,
                           (EdgeOffset + MinRadius) / Math.Sin(Angles[2]) - MinRadius
                       };
            for (int i = 0; i < 3; i++) if (r[i] < CornerOffset) r[i] = CornerOffset;

            Point a0,a1,a2;
            // arc at point 0
            double r0 = (r[0] - ControlPoints[0].DistanceTo(PanelPoints[0])) * Math.Tan(Angles[0] / 2);
            a1 = ControlPoints[0].Add(EdgeVectors[0][3].Normalized().Scale(r[0]));
            a0 = a1.Subtract(EdgeVectors[0][4].Scale(r0)).Add(EdgeVectors[0][1].Scale(r0));
            a2 = a1.Add(EdgeVectors[0][4].Scale(r0)).Add(EdgeVectors[0][0].Scale(r0));
            List<Point> arc0 =  new List<Point>{ a0, a1, a2 };
            // arc at point 1
            double r1 = (r[1] - ControlPoints[1].DistanceTo(PanelPoints[1])) * Math.Tan(Angles[1] / 2);
            a1 = ControlPoints[1].Add(EdgeVectors[1][3].Normalized().Scale(r[1]));
            a0 = a1.Subtract(EdgeVectors[1][4].Scale(r1)).Add(EdgeVectors[1][1].Scale(r1));
            a2 = a1.Add(EdgeVectors[1][4].Scale(r1)).Add(EdgeVectors[1][0].Scale(r1));
            List<Point> arc1 = new List<Point> { a0, a1, a2 };
            // arc at point 0
            double r2 = (r[2] - ControlPoints[2].DistanceTo(PanelPoints[2])) * Math.Tan(Angles[2] / 2);
            a1 = ControlPoints[2].Add(EdgeVectors[2][3].Normalized().Scale(r[2]));
            a0 = a1.Subtract(EdgeVectors[2][4].Scale(r2)).Add(EdgeVectors[2][1].Scale(r2));
            a2 = a1.Add(EdgeVectors[2][4].Scale(r2)).Add(EdgeVectors[2][0].Scale(r2));
            List<Point> arc2 = new List<Point> { a0, a1, a2 };
            a0.Dispose(); a1.Dispose(); a2.Dispose();
            ArcPoints = new List<List<Point>> { arc0, arc1, arc2};
            arc0 = arc1 = arc2 = null;
        }

        //METHOD**CREATE
        public static TrianglePanel ByMeshFace(Triangle Triangle, double EdgeOffset, double CornerOffset, double Thickness, int Direction, double MinRadius)
        { return new TrianglePanel(Triangle.GetVertexPoints(), EdgeOffset, CornerOffset, Thickness, Direction, MinRadius); }
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
            Point a = c.StartPoint;
            Point b = c.StartPoint;
            if (Direction == 0)
            {
                a = a.Add(Normal.Scale(-Thickness / 2));
                b = b.Add(Normal.Scale(Thickness / 2));
            }
            else b = b.Add(Normal.Scale(Direction * Thickness));
            Line l = Line.ByStartPointEndPoint(a,b);
            Solid s = c.SweepAsSolid(l);
            c.Dispose();
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
                Centroid.Dispose();
                Circumcenter.Dispose();
                Incenter.Dispose();
                ControlPoints.ForEach(p => p.Dispose());
                PanelPoints.ForEach(p => p.Dispose());
                ArcPoints.ForEach(p => p.ForEach(pt => pt.Dispose()));
                EdgeVectors.ForEach(v => v.Dispose());
                Normal.Dispose();
            }
            disposed = true;
        }
    }



}
