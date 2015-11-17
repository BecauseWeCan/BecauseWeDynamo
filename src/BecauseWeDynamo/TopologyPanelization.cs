using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using Text;
using Topology;

namespace Topology.Panelization
{
    public class TrianglePanel : IDisposable
    {
        //**FIELDS
        internal double Thickness;
        internal int Direction;
        bool disposed = false;

        //**PROPERTIES**QUERY
        public Triangle Triangle { get; private set; }
        public Point[][] ArcPoints { get; set; }
        public Vector[][] VertexVectors { get; private set; }

        //**METHODS**CONSTRUCTOR
        internal TrianglePanel(Triangle Triangle, double Thickness, double MinEdgeOffset, int Direction)
        {
            // initialize
            this.Thickness = Thickness;
            this.Direction = Direction;
            this.Triangle = Triangle;
            // edge vectors indexed by triangle vertex
            VertexVectors = Triangle.VertexVectors;
        }

        //**METHODS**ACTION
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
                Vector Y = Triangle.Normal.Cross(VertexVectors[j][0]);
                Word w = Word.ByStringOriginVectors(Triangle.E[j].Edge.Name, m, VertexVectors[j][0], Y);
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
                if (!VertexVectors.Equals(null)) VertexVectors.ForEach(v => v.Dispose());
            }
            disposed = true;
        }
    }

    public class TrianglePanelC : TrianglePanel
    {
        //PROPERTIES
        public double EdgeOffset { get; set; }

        //CONSTRUCTOR
        internal TrianglePanelC(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerOffset, double MinRadius, int Direction)
            : base(Triangle, Thickness, MinEdgeOffset, Direction)
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
                Point a1 = Triangle.VertexPoints[i].Add(VertexVectors[i][2].Scale(r[i]));
                Point a0 = a1.Add(VertexVectors[i][1].Scale(rV)).Subtract(VertexVectors[i][3].Scale(rV));
                Point a2 = a1.Add(VertexVectors[i][0].Scale(rV)).Add(VertexVectors[i][3].Scale(rV));
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
        internal TrianglePanelE(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerRadius, int Direction)
            : base(Triangle, Thickness, MinEdgeOffset, Direction)
        {
            // edge offsets indexed by triangle halfedge angles
            EdgeOffset = new double[] { MinEdgeOffset, MinEdgeOffset, MinEdgeOffset };
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
                Point a0 = Triangle.VertexPoints[i].Add(VertexVectors[i][1].Scale(EdgeOffset[i] / sinA + CornerRadius / tanB)).Add(VertexVectors[i][0].Scale(EdgeOffset[j] / sinA));
                Point a1 = Triangle.VertexPoints[i].Add(VertexVectors[i][1].Scale(EdgeOffset[i] / sinA)).Add(VertexVectors[i][0].Scale(EdgeOffset[j] / sinA)).Add(VertexVectors[i][2].Scale(CornerRadius / sinB - CornerRadius));
                Point a2 = Triangle.VertexPoints[i].Add(VertexVectors[i][0].Scale(EdgeOffset[j] / sinA + CornerRadius / tanB)).Add(VertexVectors[i][1].Scale(EdgeOffset[i] / sinA));
                Point[] arc = new Point[] { a0, a1, a2 };
                P.Add(arc);
            }
            ArcPoints = P.ToArray();
        }

        //METHOD**CREATE
        public static TrianglePanelE ByMeshFace(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerRadius, int Direction = 0)
        { return new TrianglePanelE(Triangle, Thickness, MinEdgeOffset, CornerRadius, Direction); }
    }

    public class EdgeConnector : IDisposable
    {
        //FIELDS
        internal double Width;
        bool disposed = false;

        public Edge Edge { get; private set; }
        public HalfEdge[] HalfEdges { get; private set; }
        public double Inset { get; private set; }
        public List<Vector> Profile { get; private set; }
        public Vector[] Vectors { get; private set; }

        internal EdgeConnector(HalfEdge e1, HalfEdge e2, double Width, double PanelThickness, double PanelMinOffset)
        {
            //initialize
            Profile = new List<Vector>();
            Edge = e1.Edge;
            this.Width = Width;
            HalfEdges = new HalfEdge[] { e1, e2 };
            int i = 0;
            int i1 = Edge.E.IndexOf(e1);
            int i2 = Edge.E.IndexOf(e2);
            if (i1 + i2 != 2)
            {
                Vector N1 = e1.Face.Normal;
                Vector Z1 = e1.GetVector().Normalized();
                Vector N2 = e2.Face.Normal;
                Vector Z2 = e2.GetVector().Normalized();
                if (i1 + i2 == 3)
                {
                    i = 1;
                    if (Edge.E.IndexOf(e1) == 1) { N1 = N1.Reverse(); Z1 = Z1.Reverse(); }
                    if (Edge.E.IndexOf(e2) == 1) { N2 = N2.Reverse(); Z2 = Z2.Reverse(); }
                }
                double EdgeOffset = Math.Max(0.5 * PanelThickness / Math.Tan(Edge.Angle[i] * Math.PI / 360), PanelMinOffset / Math.Sin(Edge.Angle[i] * Math.PI / 360) - PanelThickness / Math.Tan(Edge.Angle[i] * Math.PI / 360));
                Inset = Math.Max(Width / 2 + EdgeOffset, (PanelThickness / 2 + Width) / Math.Tan(Edge.Angle[i] * Math.PI / 360));
                Vector eN = Edge.Normal[i];

                List<Vector> P1 = new List<Vector>();
                List<Vector> P2 = new List<Vector>();
                P1.Add(eN.Scale(-PanelThickness / 2 / Math.Sin(Edge.Angle[0] * Math.PI / 360)));
                P1.Add(P1[0].Add(N1.Cross(Z1).Normalized().Scale(Inset + 1.5 * Width)));
                P1.Add(P1[1].Add(N1.Scale(-Width / 2)));
                P1.Add(P1[2].Add(Z1.Cross(N1).Normalized().Scale(Width / 2)));
                if (Edge.Angle[0] < 180) P1.Add(P1[3].Add(eN.Scale(-Width / 2)));
                else
                {
                    P1.Add(P1[3].Add(N1.Scale(-Width / 2)));
                    P1.Add(P1[0].Add(eN.Scale(-Width / Math.Sin(Edge.Angle[0] * Math.PI / 360))));
                }
                P2.Add(P1[0].Add(N2.Cross(Z2).Scale(Inset + 1.5 * Width)));
                P2.Add(P2[0].Add(N2.Scale(-Width / 2)));
                P2.Add(P2[1].Add(Z2.Cross(N2).Scale(Width / 2)));
                if (Edge.Angle[0] < 180) P2.Add(P2[2].Add(eN.Scale(-Width / 2)));
                else P2.Add(P2[2].Add(N2.Scale(-Width / 2)));
                P2.Reverse();
                Profile.AddRange(P1);
                Profile.AddRange(P2);
                Vectors = new Vector[] { N1, Z1, N2, Z2, eN };
            }
        }
        internal EdgeConnector(HalfEdge e1, double Width, double PanelThickness, double Thickness)
        {
            Profile = new List<Vector>();
            Edge = e1.Edge;
        }


        //**METHOD**CREATE
        public static EdgeConnector ByHalfEdge(HalfEdge e1, HalfEdge e2, double Width, double PanelThickness, double PanelMinOffset) { return new EdgeConnector(e1, e2, Width, PanelThickness, PanelMinOffset); }

        //**METHOD**ACTIONS
        public PolyCurve GetConnectorProfile(Point Point)
        {
            if (!(Profile.Count > 8)) return null;
            List<Point> Points = new List<Point>(Profile.Count);
            for (int i = 0; i < Profile.Count; i++) Points.Add(Point.Add(Profile[i]));
            Curve[] Curves;
            if (Profile.Count == 9)
                Curves = new Curve[] {
                          Line.ByStartPointEndPoint(Points[0], Points[1]),
                          Line.ByStartPointEndPoint(Points[1], Points[2]),
                          Arc.ByCenterPointStartPointEndPoint(Points[3],Points[2],Points[4]),
                          Line.ByStartPointEndPoint(Points[4], Points[5]),
                          Arc.ByCenterPointStartPointEndPoint(Points[6],Points[5],Points[7]),
                          Line.ByStartPointEndPoint(Points[7], Points[8]),
                          Line.ByStartPointEndPoint(Points[8], Points[0])
                      };
            else
                Curves = new Curve[] {
                          Line.ByStartPointEndPoint(Points[0], Points[1]),
                          Line.ByStartPointEndPoint(Points[1], Points[2]),
                          Arc.ByCenterPointStartPointEndPoint(Points[3],Points[2],Points[4]),
                          Line.ByStartPointEndPoint(Points[4], Points[5]),
                          Line.ByStartPointEndPoint(Points[5], Points[6]),
                          Arc.ByCenterPointStartPointEndPoint(Points[7],Points[6],Points[8]),
                          Line.ByStartPointEndPoint(Points[8], Points[9]),
                          Line.ByStartPointEndPoint(Points[9], Points[0])
                      };
            PolyCurve Result = PolyCurve.ByJoinedCurves(Curves);
            Points.ForEach(p => p.Dispose());
            Curves.ForEach(c => c.Dispose());
            return Result;
        }
        public Surface GetConnectorSurface(Point Point)
        {
            Curve c = GetConnectorProfile(Point);
            if (c.Equals(null)) { c.Dispose(); return null; }
            Surface s = Surface.ByPatch(c);
            c.Dispose();
            return s;
        }
        public Solid GetConnectorSolid(Point Point)
        {
            Curve c = GetConnectorProfile(Point);
            if (c.Equals(null)) { c.Dispose(); return null; }
            Point a = c.StartPoint.Add(Vectors[1].Scale(Width / 2));
            Point b = c.StartPoint.Add(Vectors[3].Scale(Width / 2));
            Line l = Line.ByStartPointEndPoint(a, b);
            Solid s = c.SweepAsSolid(l);
            c.Dispose();
            a.Dispose();
            b.Dispose();
            l.Dispose();
            return s;
        }
        public PolyCurve[] GetEdgeLabels(Point Point, double Scale = 1)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            Vector X = Vectors[0].Cross(Vectors[1]).Normalized();
            Word w = Word.ByStringOriginVectors(Edge.Name, Point.Add(Profile[1].Subtract(Profile[0]).Scale(1 / 2)), X, Vectors[0]);
            labels.AddRange(w.display(Scale * Width / 32));
            X.Dispose(); w.Dispose();
            return labels.ToArray();
        }
        public Circle[] GetPockets(Point Point)
        {
            return null;
        }

        //**METHODS**DISPOSE
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
                if (!Profile.Equals(null)) Profile.ForEach(p => p.Dispose());
                if (!Vectors.Equals(null)) Vectors.ForEach(v => v.Dispose());
            }
            disposed = true;
        }
    }

    public class Panelization
    {
        Mesh m;

    }
}
