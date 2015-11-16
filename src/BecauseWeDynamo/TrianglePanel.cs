using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using Topology;
using Text;

namespace Panelization
{
    public class TrianglePanel: IDisposable
    {
        internal double Thickness;
        internal int Direction;
        bool disposed = false;

        //PROPERTIES**QUERY
        public Triangle Triangle { get; private set; }
        public Point[][] ArcPoints { get; set; }
        public Vector[][] EdgeVectors { get; private set; }

        internal TrianglePanel(Triangle Triangle, double Thickness, double MinEdgeOffset, int Direction)
        {
            // initialize
            this.Thickness = Thickness;
            this.Direction = Direction;
            this.Triangle = Triangle;
            // edge vectors indexed by triangle vertex
            List<Vector[]> eV = new List<Vector[]>(3);
            for (int i = 0; i < 3; i++)
            {
                List<Vector> V = new List<Vector> { Triangle.E[i].GetVector().Normalized(), Triangle.E[(i + 2) % 3].GetVector().Normalized().Reverse() };
                V.Add(V[0].Add(V[1]).Normalized());
                V.Add(V[0].Subtract(V[1]).Normalized());
                eV.Add(V.ToArray());
            }
            EdgeVectors = eV.ToArray();
        }

        //METHODS**ACTION
        public PolyCurve GetPanelProfile()
        {
            if (ArcPoints.Equals(null)) return null;
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
            if (c.Equals(null)) { c.Dispose(); return null; }
            Surface s = Surface.ByPatch(c);
            c.Dispose();
            return s;
        }
        public Solid GetPanelSolid()
        {
            Curve c = GetPanelProfile();
            if (c.Equals(null)) { c.Dispose(); return null; }
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
        public PolyCurve[] GetEdgeLabels(double Scale)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < Triangle.E.Count; j++)
            {
                if (Triangle.E[j].Edge.E.Count == 1) continue;
                Point m = Triangle.E[j].Edge.MidPoint;
                Vector Y = Triangle.Normal.Cross(EdgeVectors[j][0]);
                Word w = Word.ByStringOriginVectors(Triangle.E[j].Edge.Name, m, EdgeVectors[j][0], Y);
                labels.AddRange(w.display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            return labels.ToArray();
        }
        public Circle[] GetHoles()
        {
            return null;
        }
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                if (!ArcPoints.Equals(null)) ArcPoints.ForEach(p => p.ForEach(pt => pt.Dispose()));
                if (!EdgeVectors.Equals(null)) EdgeVectors.ForEach(v => v.Dispose());
            }
            disposed = true;
        }
    }

    public class TrianglePanelC : TrianglePanel
    {
        //PROPERTIES
        public double EdgeOffset { get; set; }

        //CONSTRUCTOR
        internal TrianglePanelC(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerOffset, double MinRadius, int Direction) : base(Triangle, Thickness, MinEdgeOffset, Direction)
        {
            // edge offsets indexed by triangle halfedge angles
            EdgeOffset = 0.5 * Thickness / Math.Tan(Triangle.MinEdgeAngle * Math.PI / 360);
            for (int i = 0; i < 3; i++)
            {
                if (Triangle.E[i].Angle == 360) continue;
                double d = MinEdgeOffset / Math.Sin(Triangle.E[i].Angle * Math.PI / 360) - 0.5 * Thickness / Math.Tan(Triangle.E[i].Angle * Math.PI / 360);
                if (EdgeOffset < d) EdgeOffset = d;
            }
            // corner offsets based on triangle vertex
            double[] r = {
                           (EdgeOffset + MinRadius) / Math.Sin(Triangle.Angles[0]/2) - MinRadius,
                           (EdgeOffset + MinRadius) / Math.Sin(Triangle.Angles[1]/2) - MinRadius,
                           (EdgeOffset + MinRadius) / Math.Sin(Triangle.Angles[2]/2) - MinRadius
                       };
            for (int i = 0; i < r.Length; i++) if (r[i] < CornerOffset) r[i] = CornerOffset;
            // corner arcs based on triangle vertex
            List<Point[]> P = new List<Point[]>(3);
            for (int i = 0; i < 3; i++)
            {
                double rV = r[i] * Math.Tan(Triangle.Angles[i] / 2) - EdgeOffset / Math.Cos(Triangle.Angles[i] / 2);
                Point a1 = Triangle.VertexPoints[i].Add(EdgeVectors[i][2].Scale(r[i]));
                Point a0 = a1.Add(EdgeVectors[i][1].Scale(rV)).Subtract(EdgeVectors[i][3].Scale(rV));
                Point a2 = a1.Add(EdgeVectors[i][0].Scale(rV)).Add(EdgeVectors[i][3].Scale(rV));
                Point[] arc0 = new Point[] { a0, a1, a2 };
                P.Add(arc0);
            }
            ArcPoints = P.ToArray();
        }

        //METHOD**CREATE
        public static TrianglePanelC ByMeshFace(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerOffset, double MinRadius, int Direction = 0)
        { return new TrianglePanelC(Triangle, Thickness, MinEdgeOffset, CornerOffset, MinRadius, Direction); }
    }

    public class TrianglePanelE : TrianglePanel
    {
        //PROPERTIES**QUERY
        public double[] EdgeOffset { get; set; }

        //CONSTRUCTOR
        internal TrianglePanelE(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerRadius, int Direction): base(Triangle,Thickness,MinEdgeOffset,Direction)
        {
            // edge offsets indexed by triangle halfedge angles
            EdgeOffset = new double[] { 0, 0, 0 };
            for (int i = 0; i < 3; i++)
            {
                if (Triangle.E[i].Angle == 360) continue;
                EdgeOffset[i] = 0.5 * Thickness / Math.Tan(Triangle.E[i].Angle * Math.PI / 360);
                double OffsetAngle = MinEdgeOffset / Math.Sin(Triangle.E[i].Angle * Math.PI / 360) - 0.5 * Thickness / Math.Tan(Triangle.E[i].Angle * Math.PI / 360);
                if (EdgeOffset[i] < OffsetAngle) EdgeOffset[i] = OffsetAngle;
            }
            // corner arcs based on triangle vertex
            List<Point[]> P = new List<Point[]>(3);
            for (int i = 0; i < 3; i++)
            {
                double sinA = Math.Sin(Triangle.Angles[i]);
                double sinB = Math.Sin(Triangle.Angles[i] / 2);
                double tanB = Math.Tan(Triangle.Angles[i] / 2);
                int j = (i + 2) % 3;
                Point a0 = Triangle.VertexPoints[i].Add(EdgeVectors[i][1].Scale(EdgeOffset[i] / sinA + CornerRadius / tanB)).Add(EdgeVectors[i][0].Scale(EdgeOffset[j] / sinA));
                Point a1 = Triangle.VertexPoints[i].Add(EdgeVectors[i][1].Scale(EdgeOffset[i] / sinA)).Add(EdgeVectors[i][0].Scale(EdgeOffset[j] / sinA)).Add(EdgeVectors[i][2].Scale(CornerRadius / sinB - CornerRadius));
                Point a2 = Triangle.VertexPoints[i].Add(EdgeVectors[i][0].Scale(EdgeOffset[j] / sinA + CornerRadius / tanB)).Add(EdgeVectors[i][1].Scale(EdgeOffset[i] / sinA));
                Point[] arc = new Point[] { a0, a1, a2 };
                P.Add(arc);
            }
            ArcPoints = P.ToArray();
        }

        //METHOD**CREATE
        public static TrianglePanelE ByMeshFace(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerRadius, int Direction = 0)
        { return new TrianglePanelE(Triangle, Thickness, MinEdgeOffset, CornerRadius, Direction); }
    }
}
