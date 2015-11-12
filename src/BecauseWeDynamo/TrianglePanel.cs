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
        internal double Thickness;
        internal int Direction;
        bool disposed = false;

        //PROPERTIES**QUERY
        public double EdgeOffset { get; set; }
        public Triangle Triangle { get; set; }
        public Point[][] ArcPoints { get; private set; }
        public Vector[][] EdgeVectors { get; private set; }

        //CONSTRUCTOR
        internal TrianglePanel(Triangle Triangle, double CornerOffset, double MinEdgeOffset, double Thickness, int Direction, double MinRadius)
        {
            this.Thickness = Thickness;
            this.Direction = Direction;
            this.Triangle = Triangle;
            double[] s = { Triangle.Edges[0].Length, Triangle.Edges[1].Length, Triangle.Edges[2].Length };
            List<Vector[]> eV = new List<Vector[]>(3);
            EdgeOffset = 0.5 * Thickness / Math.Tan(Triangle.MinEdgeAngle * Math.PI / 360);
            for (int i = 0; i < 3; i++)
            {
                double d = MinEdgeOffset / Math.Sin(Triangle.Angles[0] / 2) - 0.5 * Thickness / Math.Tan(Triangle.Angles[0]/2);
                if (EdgeOffset < d) EdgeOffset = d;
            }
            for (int i = 0; i < 3; i++)
            {
                List<Vector> V = new List<Vector> { Triangle.E[i].GetVector().Normalized(), Triangle.E[(i + 2) % 3].GetVector().Normalized().Reverse() };
                V.Add(V[0].Add(V[1]).Normalized());
                V.Add(V[0].Subtract(V[1]).Normalized());
                eV.Add(V.ToArray());
            }
            EdgeVectors = eV.ToArray();
            // corner offsets
            double[] r = {
                           (EdgeOffset + MinRadius) / Math.Sin(Triangle.Angles[0]/2) - MinRadius,
                           (EdgeOffset + MinRadius) / Math.Sin(Triangle.Angles[1]/2) - MinRadius,
                           (EdgeOffset + MinRadius) / Math.Sin(Triangle.Angles[2]/2) - MinRadius
                       };
            for (int i = 0; i < r.Length; i++) if (r[i] < CornerOffset) r[i] = CornerOffset;
            // arc at point 0
            List<Point[]> P = new List<Point[]>(3);
            for (int i = 0; i < 3; i++)
            {
                double rV = r[i] * Math.Tan(Triangle.Angles[i] / 2) - EdgeOffset/Math.Cos(Triangle.Angles[i]/2);
                Point a1 = Triangle.VertexPoints[i].Add(EdgeVectors[i][2].Scale(r[i]));
                Point a0 = a1.Add(EdgeVectors[i][1].Scale(rV)).Subtract(EdgeVectors[i][3].Scale(rV));
                Point a2 = a1.Add(EdgeVectors[i][0].Scale(rV)).Add(EdgeVectors[i][3].Scale(rV));
                Point[] arc0 = new Point[] { a0, a1, a2 };
                P.Add(arc0);
            }
            ArcPoints = P.ToArray();
        }

        //METHOD**CREATE
        public static TrianglePanel ByMeshFace(Triangle Triangle, double CornerOffset, double MinEdgeOffset, double Thickness, double MinRadius, int Direction = 0)
        { return new TrianglePanel(Triangle, CornerOffset, MinEdgeOffset, Thickness, Direction, MinRadius); }


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
                a = a.Add(Triangle.Normal.Scale(-Thickness / 2));
                b = b.Add(Triangle.Normal.Scale(Thickness / 2));
            }
            else b = b.Add(Triangle.Normal.Scale(Direction * Thickness));
            Line l = Line.ByStartPointEndPoint(a, b);
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
                ArcPoints.ForEach(p => p.ForEach(pt => pt.Dispose()));
                EdgeVectors.ForEach(v => v.Dispose());
            }
            disposed = true;
        }
    }



}
