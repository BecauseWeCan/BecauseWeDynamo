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
        public List<Circle> Holes { get; set; }

        //**METHODS**CONSTRUCTOR
        internal TrianglePanel(Triangle Triangle, double Thickness, double MinEdgeOffset, int Direction)
        {
            // initialize
            this.Thickness = Thickness;
            this.Direction = Direction;
            this.Triangle = Triangle;
            Holes = new List<Circle>();
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
        public Solid GetPanelSolidHole()
        {
            Surface s = GetPanelSurface();
            Geometry[] G0 = s.Split(Holes[0]);
            Geometry[] G1 = G0[1].Split(Holes[1]);
            Geometry[] G2 = G1[1].Split(Holes[2]);
            Geometry[] G3 = G2[1].Split(Holes[3]);
            Geometry[] G4 = G3[1].Split(Holes[4]);
            Geometry[] G5 = G4[1].Split(Holes[5]);
            Solid S = (G5[1] as Surface).Thicken(Thickness);
            s.Dispose();
            G0.ForEach(g => g.Dispose());
            G1.ForEach(g => g.Dispose());
            G2.ForEach(g => g.Dispose());
            G3.ForEach(g => g.Dispose());
            G4.ForEach(g => g.Dispose());
            G5.ForEach(g => g.Dispose());
            return S;
        }
        public PolyCurve[] GetEdgeLabels(double Scale)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < Triangle.E.Count; j++)
            {
                if (Triangle.E[j].Edge.E.Count == 1) continue;
                Point m = Point.ByCoordinates(ArcPoints[j][2].X / 2 + ArcPoints[(j + 1) % 3][0].X / 2, ArcPoints[j][2].Y / 2 + ArcPoints[(j + 1) % 3][0].Y / 2, ArcPoints[j][2].Z / 2 + ArcPoints[(j + 1) % 3][0].Z / 2);
                Vector Y = VertexVectors[j][0].Cross(Triangle.Normal);
                Word w = Word.ByStringOriginVectors(Triangle.E[j].Edge.Name, m, VertexVectors[j][0], Y);
                labels.AddRange(w.display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            return labels.ToArray();
        }
        public PolyCurve[] GetLabels(double Scale, PanelSystem P)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < Triangle.E.Count; j++)
            {
                if (!P.E.ContainsKey(Triangle.E[j].Edge)) continue;
                Point m = Point.ByCoordinates(ArcPoints[j][2].X / 2 + ArcPoints[(j + 1) % 3][0].X / 2, ArcPoints[j][2].Y / 2 + ArcPoints[(j + 1) % 3][0].Y / 2, ArcPoints[j][2].Z / 2 + ArcPoints[(j + 1) % 3][0].Z / 2);
                Vector Y = VertexVectors[j][0].Cross(Triangle.Normal);
                m = m.Subtract(Y.Scale(P.E[Triangle.Edges[j]][0].Inset + 2 * P.E[Triangle.Edges[j]][0].Width));
                Word w = Word.ByStringOriginVectors(Triangle.E[j].Edge.Name, m, VertexVectors[j][0], Y);
                labels.AddRange(w.display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            Point c = Triangle.Center;
            Vector N = VertexVectors[0][0].Cross(Triangle.Normal);
            Word W = Word.ByStringOriginVectors(Triangle.Name, c, VertexVectors[0][0], N);
            labels.AddRange(W.display(2 * Scale));
            c.Dispose(); N.Dispose(); W.Dispose();
            return labels.ToArray();
        }
        public PolyCurve[] GetTriangleLabel(double Scale)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            Point c = Triangle.Center;
            Vector X = VertexVectors[0][0].Reverse();
            Vector Y = VertexVectors[0][0].Cross(Triangle.Normal);
            Word W = Word.ByStringOriginVectors(Triangle.Name, c, X, Y);
            labels.AddRange(W.display(Scale));
            c.Dispose(); X.Dispose(); Y.Dispose(); W.Dispose();
            return labels.ToArray();
        }
        public void AddHoles(Point Point, double Radius)
        {
            Holes.Add(Circle.ByCenterPointRadiusNormal(Point, Radius, Triangle.Normal));
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
                if (!Holes.Equals(null)) Holes.ForEach(h => h.Dispose());
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
        public double Angle { get; private set; }
        public List<Vector> Profile { get; private set; }
        public Vector[] Vectors { get; private set; }
        public List<Circle> Holes { get; set; }

        internal EdgeConnector(HalfEdge e1, HalfEdge e2, double Width, double PanelThickness, double PanelMinOffset)
        {
            //initialize
            Profile = new List<Vector>();
            Holes = new List<Circle>();
            HalfEdges = new HalfEdge[] { e1, e2 };
            Edge = e1.Edge;
            this.Width = Width;
            int i = 0;
            int i1 = Edge.E.IndexOf(e1);
            int i2 = Edge.E.IndexOf(e2);
            if (i1 + i2 != 2)
            {
                Vector X1 = e1.Face.Normal;
                Vector Z1 = e1.GetVector().Normalized();
                Vector Y1 = X1.Cross(Z1).Normalized();
                Vector X2 = e2.Face.Normal;
                Vector Z2 = e2.GetVector().Normalized();
                Vector Y2 = X2.Cross(Z2).Normalized();
                if (i1 + i2 == 3)
                {
                    i = 1;
                    if (Edge.E.IndexOf(e1) == 1) { X1 = X1.Reverse(); Z1 = Z1.Reverse(); }
                    if (Edge.E.IndexOf(e2) == 1) { X2 = X2.Reverse(); Z2 = Z2.Reverse(); }
                }
                Angle = Edge.Angle[0];
                Vector eN = Edge.Normal[i].Normalized();
                if (Edge.E.Count > 2)
                {
                    Angle = Math.Min(Edge.Angle[0], Edge.Angle[1]);
                    if (i == 0) eN = eN.Reverse();
                }
                double EdgeOffset = Math.Max(0.5 * PanelThickness / Math.Tan(Angle * Math.PI / 360), PanelMinOffset / Math.Sin(Angle * Math.PI / 360) - PanelThickness / Math.Tan(Angle * Math.PI / 360));
                Inset = Math.Max(Width / 2 + EdgeOffset, (PanelThickness / 2 + Width) / Math.Tan(Angle * Math.PI / 360));

                List<Vector> P1 = new List<Vector>();
                List<Vector> P2 = new List<Vector>();
                P1.Add(eN.Scale(-PanelThickness / 2 / Math.Sin(Edge.Angle[i] * Math.PI / 360)));
                P1.Add(X1.Scale(-PanelThickness / 2).Add(Y1.Normalized().Scale(Inset + 1.5 * Width)));
                P1.Add(P1[1].Add(X1.Scale(-Width / 2)));
                P1.Add(P1[2].Add(Y1.Normalized().Scale(-Width / 2)));
                if (Edge.Angle[i] < 180) P1.Add(P1[3].Add(eN.Scale(-Width / 2)));
                else
                {
                    P1.Add(P1[3].Add(X1.Scale(-Width / 2)));
                    P1.Add(P1[0].Add(eN.Scale(-Width / Math.Sin(Edge.Angle[i] * Math.PI / 360))));
                }
                P2.Add(X2.Scale(-PanelThickness / 2).Add(Y2.Scale(Inset + 1.5 * Width)));
                P2.Add(P2[0].Add(X2.Scale(-Width / 2)));
                P2.Add(P2[1].Add(Y2.Normalized().Scale(-Width / 2)));
                if (Edge.Angle[i] < 180) P2.Add(P2[2].Add(eN.Scale(-Width / 2)));
                else P2.Add(P2[2].Add(X2.Scale(-Width / 2)));
                P2.Reverse();
                Profile.AddRange(P1);
                Profile.AddRange(P2);
                Vectors = new Vector[] { X1, Y1, Z1, X2, Y2, Z2, eN };
            }
        }
        internal EdgeConnector(HalfEdge e1, double Width, double PanelThickness, double Thickness)
        {
            Profile = new List<Vector>();
            Edge = e1.Edge;
        }


        //**METHOD**CREATE
        public static EdgeConnector ByHalfEdge(HalfEdge e1, HalfEdge e2, double Width, double PanelThickness, double PanelMinOffset) { return new EdgeConnector(e1, e2, Width, PanelThickness, PanelMinOffset); }
        public static Object ByEdge(Edge e, double Width, double PanelThickness, double PanelMinOffset)
        {
            if (e.E.Count < 2) return null;
            if (e.E.Count == 2) return new EdgeConnector(e.E[0], e.E[1], Width, PanelThickness, PanelMinOffset);
            EdgeConnector[] result = {new EdgeConnector(e.E[0], e.E[1], Width, PanelThickness, PanelMinOffset),
                                         new EdgeConnector(e.E[1], e.E[2], Width, PanelThickness, PanelMinOffset) };
            return result;
        }


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
            Point a = c.StartPoint.Add(Vectors[2].Scale(Width / 2));
            Point b = c.StartPoint.Add(Vectors[5].Scale(Width / 2));
            Line l = Line.ByStartPointEndPoint(a, b);
            Solid s = c.SweepAsSolid(l);
            c.Dispose();
            a.Dispose();
            b.Dispose();
            l.Dispose();
            return s;
        }
        public Solid GetConnectorSolidHoles(Point Point)
        {
            Surface s = GetConnectorSurface(Point);
            Geometry[] G0 = s.Split(Holes[0]);
            Geometry[] G1 = G0[1].Split(Holes[1]);
            Geometry[] G2 = G1[1].Split(Holes[2]);
            Geometry[] G3 = G2[1].Split(Holes[3]);
            Solid S = (G3[1] as Surface).Thicken(Width, true);
            s.Dispose();
            G0.ForEach(g => g.Dispose());
            G1.ForEach(g => g.Dispose());
            G2.ForEach(g => g.Dispose());
            G3.ForEach(g => g.Dispose());
            return S;
        }
        public PolyCurve[] GetEdgeLabel(Point Point, double Scale = 1/32)
        {
            if (!(Profile.Count > 8)) return null;
            Point p4 = Point.Add(Profile[4]);
            Point p5 = Point.Add(Profile[5]);
            Point pt = p5;
            if (Profile.Count == 9) pt = Point.ByCoordinates(p4.X / 2 + p5.X / 2, p4.Y / 2 + p5.Y / 2, p4.Z / 2 + p5.Z / 2);
            Vector X = Vectors[3].Subtract(Vectors[0]);
            Vector Y = Vectors[6].Reverse();
            if (Angle > 180)
            {
                pt = Point.Add(Profile[0]);
                X = X.Reverse();
                Y = Y.Reverse();
            }
            Surface s = GetConnectorSurface(Point);
            CoordinateSystem cs = s.CoordinateSystemAtParameter(0.5, 0.5);
            Vector Z1 = cs.ZAxis.Normalized();
            Vector Z2 = X.Cross(Y).Normalized();
            if (!Z1.IsAlmostEqualTo(Z2)) X = X.Reverse();
            Word w = Word.ByStringOriginVectors(Edge.Name, pt, X, Y);
            PolyCurve[] label = w.display(Scale * Width).ToArray();
            p4.Dispose(); p5.Dispose(); pt.Dispose();
            X.Dispose(); Y.Dispose(); Z1.Dispose(); Z2.Dispose();
            s.Dispose(); cs.Dispose(); w.Dispose();
            return label;
        }
        public void AddPockets(Point Point, double Radius)
        {
            double a = Angle;
            if (Edge.Angle.Length > 1) a = Math.Min(Edge.Angle[0], Edge.Angle[1]);
            Point p1a = Point.Add(Profile[3]);
            Point p1b = p1a.Subtract(Vectors[1].Normalized().Scale(Width));
            Point p2a = Point.Add(Profile[Profile.Count - 3]);
            Point p2b = p2a.Subtract(Vectors[4].Normalized().Scale(Width));
            if (a == Angle || Edge.Angle[1] == a)
            {
                Holes.Add(Circle.ByCenterPointRadiusNormal(p1a, Radius, Vectors[5]));
                Holes.Add(Circle.ByCenterPointRadiusNormal(p1b, Radius, Vectors[5]));
            }
            if (a == Angle || Edge.Angle[0] == a)
            {
                Holes.Add(Circle.ByCenterPointRadiusNormal(p2a, Radius, Vectors[5]));
                Holes.Add(Circle.ByCenterPointRadiusNormal(p2b, Radius, Vectors[5]));
            }
            p1a.Dispose(); p1b.Dispose(); p2a.Dispose(); p2b.Dispose();
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
                if (!Holes.Equals(null)) Holes.ForEach(h => h.Dispose());
            }
            disposed = true;
        }
    }

    public class PanelSystem : IDisposable
    {
        //**FIELDS
        bool disposed = false;
        internal Dictionary<Triangle, TrianglePanel> T;
        internal Dictionary<Edge, EdgeConnector[]> E;
        internal TriangleMesh M;

        //**PROPERTIES**QUERY
        public List<TrianglePanel> Panels { get { return T.Values.ToList(); } }
        public List<EdgeConnector[]> Connectors { get { return E.Values.ToList(); } }

        //**CONSTRUCTOR
        internal PanelSystem(TriangleMesh Mesh) { M = Mesh; }
        internal PanelSystem(TriangleMesh Mesh, double Width, double Thickness, double MinEdgeOffset, double CornerRadius, double HoleRadius, double PocketRadius)
            : this(Mesh)
        { GetEdgeConnectors(Width, Thickness, MinEdgeOffset, PocketRadius); GetTrianglePanels(Thickness, MinEdgeOffset, CornerRadius, HoleRadius); }

        //**METHODS**CREATE
        public static PanelSystem ByMesh(TriangleMesh Mesh, double Width, double Thickness, double MinEdgeOffset, double CornerRadius, double HoleRadius, double PocketRadius)
        { return new PanelSystem(Mesh, Width, Thickness, MinEdgeOffset, CornerRadius, HoleRadius, PocketRadius); }

        public List<EdgeConnector[]> GetEdgeConnectors(double Width, double PanelThickness, double PanelMinOffset, double PocketRadius)
        {
            E = new Dictionary<Edge, EdgeConnector[]>(M.E2.Count + M.E3.Count);
            M.E2.ForEach(e => E.Add(e, new EdgeConnector[] { new EdgeConnector(e.E[0], e.E[1], Width, PanelThickness, PanelMinOffset) }));
            for (int i = 0; i < M.E3.Count; i++)
                E.Add(M.E3[i], new EdgeConnector[]{
                    new EdgeConnector(M.E3[i].E[0], M.E3[i].E[1], Width, PanelThickness, PanelMinOffset),
                    new EdgeConnector(M.E3[i].E[1], M.E3[i].E[2], Width, PanelThickness, PanelMinOffset)});
            Connectors.ForEach(C => C.ForEach(c => c.AddPockets(c.Edge.MidPoint, PocketRadius)));
            return Connectors;
        }
        public List<TrianglePanel> GetTrianglePanels(double Thickness, double MinEdgeOffset, double CornerRadius, double HoleRadius)
        {
            T = new Dictionary<Triangle, TrianglePanel>(M.Faces.Count);
            for (int i = 0; i < M.Faces.Count; i++)
            {
                TrianglePanel t = TrianglePanelE.ByMeshFace(M.Faces[i], Thickness, MinEdgeOffset, CornerRadius);
                T.Add(M.Faces[i], t);
                for (int j = 0; j < M.Faces[i].Edges.Length; j++)
                    if (E.ContainsKey(M.Faces[i].Edges[j]))
                    {
                        Point p = M.Faces[i].Edges[j].MidPoint;
                        Vector Y = M.Faces[i].Normal.Cross(M.Faces[i].E[j].GetVector()).Normalized();
                        p = p.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Inset));
                        t.AddHoles(p, HoleRadius);
                        p = p.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Width));
                        t.AddHoles(p, HoleRadius);
                        p.Dispose();
                        Y.Dispose();
                    }
            }
            return Panels;
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
                if (!Panels.Equals(null)) Panels.ForEach(t => t.Dispose());
                if (!E.Equals(null)) E.Values.ToArray().ForEach(e => e.Dispose());
            }
            disposed = true;
        }
    }
}
