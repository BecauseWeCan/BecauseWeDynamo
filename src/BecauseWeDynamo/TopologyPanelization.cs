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

namespace Panelization
{

    public class PanelBevel: IDisposable
    {
        //**FIELDS
        bool disposed = false;
        double Thickness;
        double EdgeThickness;
        double BevelAngle;


        //**PROPERTIES**
        public Topology.Face Face { get; private set; }
        /// <summary>
        /// gets point data for panel;
        /// i: index that refers to the vertex
        /// j: index of arc points at vertex with order inherited by triangle face (cc)
        /// </summary>
        public Point[][] ArcPoints { get; set; }
        public double[] EdgeOffset { get; set; }

        //**METHODS**CONSTRUCTOR
        internal PanelBevel(Topology.Face Face, double Thickness, double EdgeThickness, double MinEdgeOffset, double CornerRadius, double BevelAngle = 0)
        {
            // initialize
            this.Thickness = Thickness;
            this.Face = Face;
            // edge offsets indexed by triangle halfedge angles
            EdgeOffset = new double[Face.E.Count];
            for (int i = 0; i < Face.E.Count; i++)
            {
                if (Face.E[i].Angle == 360) continue;
                EdgeOffset[i] = 0.5 * Thickness / Math.Tan(Face.E[i].Angle * Math.PI / 360);
                double OffsetAngle = MinEdgeOffset / Math.Sin(Face.E[i].Angle * Math.PI / 360) - 0.5 * Thickness / Math.Tan(Face.E[i].Angle * Math.PI / 360);
                if (EdgeOffset[i] < OffsetAngle) EdgeOffset[i] = OffsetAngle;
            }
            // corner arcs based on triangle vertex
            List<Point[]> P = new List<Point[]>(3);
            for (int i = 0; i < 3; i++)
            {
                double sinA = Math.Sin(Face.Angles[i]);
                double sinB = Math.Sin(Face.Angles[i] / 2);
                double tanB = Math.Tan(Face.Angles[i] / 2);
                int j = (i + 2) % 3;
                Point a0 = Face.VertexPoints[i].Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA + CornerRadius / tanB)).Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA));
                Point a1 = Face.VertexPoints[i].Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA)).Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA)).Add(Face.VertexVectors[i][2].Scale(CornerRadius / sinB - CornerRadius));
                Point a2 = Face.VertexPoints[i].Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA + CornerRadius / tanB)).Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA));
                Point[] arc = new Point[] { a0, a1, a2 };
                P.Add(arc);
            }
            ArcPoints = P.ToArray();
        }

        //**METHOD**CREATE
        /// <summary>
        /// creates TrianglePanel with unique edge offsets determined by angle at edge, with given corner radius
        /// </summary>
        /// <param name="Triangle">reference triangular face on mesh</param>
        /// <param name="Thickness">thickness of panel in mesh units</param>
        /// <param name="MinEdgeOffset">minimum edge offset</param>
        /// <param name="CornerRadius">radius of corner fillet</param>
        /// <param name="Direction">direction of panel extrusion into a solid</param>
        /// <returns>TrianglePanel</returns>
        public static PanelBevel ByMeshFace(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerRadius, int Direction = 0)
        { return new PanelBevel(Triangle, Thickness, MinEdgeOffset, CornerRadius, Direction); }
    

        //**METHODS**ACTION
        /// <summary>
        /// returns profile curve of panel as single Polycurve
        /// </summary>
        /// <returns>Polycurve of joined Arcs and Lines</returns>
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
        /// <summary>
        /// returns flat panel on mesh face as single Surface
        /// </summary>
        /// <returns>Surface by Polycurve patch</returns>
        public Surface GetPanelSurface()
        {
            Curve c = GetPanelProfile();
            if (c.Equals(null)) { c.Dispose(); return null; }
            Surface s = Surface.ByPatch(c);
            c.Dispose();
            return s;
        }
        /// <summary>
        /// returns solid panel based on panel profile as single Solid
        /// </summary>
        /// <returns>Solid by Polycurve sweep</returns>
        public Solid GetPanelSolid()
        {
            Curve c = GetPanelProfile();
            if (c.Equals(null)) { c.Dispose(); return null; }
            double Bevel = (EdgeThickness-Thickness) * Math.Tan(BevelAngle * Math.PI /360);
            Vector In = Face.Normal.Cross(Face.E[0].GetVector().Normalized()).Normalized();
            Point[] pts = new Point[]{c.StartPoint.Add(Face.Normal.Scale(EdgeThickness)),
                              c.StartPoint,
                              c.StartPoint.Add(Face.Normal.Scale(EdgeThickness-Thickness)).Add(In.Scale(Bevel))};
            
            Curve l = PolyCurve.ByPoints(pts);
            Solid s = c.SweepAsSolid(l);
            c.Dispose(); l.Dispose(); pts.ForEach(p => p.Dispose());
            return s;
        }
        
        /// <summary>
        /// returns edge labels at edge on triangular mesh face backside as an Polycurve Array
        /// </summary>
        /// <param name="Scale">Scale = (letter height)/4; in mesh units</param>
        /// <returns>Polycurve Array</returns>
        public PolyCurve[] GetEdgeLabels(double Scale)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < Face.E.Count; j++)
            {
                if (Face.E[j].Edge.E.Count == 1) continue;
                Point m = Point.ByCoordinates(ArcPoints[j][2].X / 2 + ArcPoints[(j + 1) % 3][0].X / 2, ArcPoints[j][2].Y / 2 + ArcPoints[(j + 1) % 3][0].Y / 2, ArcPoints[j][2].Z / 2 + ArcPoints[(j + 1) % 3][0].Z / 2);
                Vector Y = Face.VertexVectors[j][0].Cross(Face.Normal);
                Word w = Word.ByStringOriginVectors(Face.E[j].Edge.Name, m, Face.VertexVectors[j][0], Y);
                labels.AddRange(w.display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            return labels.ToArray();
        }
        /// <summary>
        /// returns edge labels at edge connector and triangle label at triangle center 
        /// on triangular mesh face backside as an Polycurve Array
        /// </summary>
        /// <param name="Scale">Scale = (letter height)/4; in mesh units</param>
        /// <param name="PanelSystem">parent PanelSystem of Panel</param>
        /// <returns>Polycurve Array</returns>
        public PolyCurve[] GetLabels(double Scale, PanelSystem PanelSystem)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < Face.E.Count; j++)
            {
                if (!PanelSystem.E.ContainsKey(Face.E[j].Edge)) continue;
                Point m = Point.ByCoordinates(ArcPoints[j][2].X / 2 + ArcPoints[(j + 1) % 3][0].X / 2, ArcPoints[j][2].Y / 2 + ArcPoints[(j + 1) % 3][0].Y / 2, ArcPoints[j][2].Z / 2 + ArcPoints[(j + 1) % 3][0].Z / 2);
                Vector Y = Face.VertexVectors[j][0].Cross(Face.Normal);
                m = m.Subtract(Y.Scale(PanelSystem.E[Face.Edges[j]][0].Inset + 2 * PanelSystem.E[Face.Edges[j]][0].Width));
                Word w = Word.ByStringOriginVectors(Face.E[j].Edge.Name, m, Face.VertexVectors[j][0], Y);
                labels.AddRange(w.display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            Point c = Face.Center;
            Vector N = Face.VertexVectors[0][0].Cross(Face.Normal);
            Word W = Word.ByStringOriginVectors(Face.Name, c, Face.VertexVectors[0][0], N);
            labels.AddRange(W.display(2 * Scale));
            c.Dispose(); N.Dispose(); W.Dispose();
            return labels.ToArray();
        }
        /// <summary>
        /// returns triangle label at triangle center on triangular mesh face frontside as an Polycurve Array
        /// </summary>
        /// <param name="Scale"></param>
        /// <returns></returns>
        public PolyCurve[] GetTriangleLabel(double Scale)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            Point c = Face.Center;
            Vector X = Face.VertexVectors[0][0].Reverse();
            Vector Y = Face.VertexVectors[0][0].Cross(Face.Normal);
            Word W = Word.ByStringOriginVectors(Face.Name, c, X, Y);
            labels.AddRange(W.display(Scale));
            c.Dispose(); X.Dispose(); Y.Dispose(); W.Dispose();
            return labels.ToArray();
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
                if (!ArcPoints.Equals(null)) ArcPoints.ForEach(p => p.ForEach(pt => pt.Dispose()));
            }
            disposed = true;
        }
    }

    /// <summary>
    /// TrianglePanel: triangular panel base class with filleted edges based on existing triangular face with edge offsets, thickness, and extrusion direction
    /// </summary>
    public class TrianglePanel : IDisposable
    {
        //**FIELDS
        internal double Thickness;
        internal int Direction;
        bool disposed = false;

        //**PROPERTIES**QUERY
        /// <summary>
        /// gets reference triangular face on mesh
        /// </summary>
        public Triangle Triangle { get; private set; }
        /// <summary>
        /// gets point data for panel;
        /// i: index that refers to the vertex
        /// j: index of arc points at vertex with order inherited by triangle face (cc)
        /// </summary>
        public Point[][] ArcPoints { get; set; }
        /// <summary>
        /// gets hole information set by panelization of mesh
        /// </summary>
        public List<Circle> Holes { get; set; }

        //**METHODS**CONSTRUCTOR
        internal TrianglePanel(Triangle Triangle, double Thickness, double MinEdgeOffset, int Direction)
        {
            // initialize
            this.Thickness = Thickness;
            this.Direction = Direction;
            this.Triangle = Triangle;
            Holes = new List<Circle>();
        }

        //**METHODS**ACTION
        /// <summary>
        /// returns profile curve of panel as single Polycurve
        /// </summary>
        /// <returns>Polycurve of joined Arcs and Lines</returns>
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
        /// <summary>
        /// returns flat panel on mesh face as single Surface
        /// </summary>
        /// <returns>Surface by Polycurve patch</returns>
        public Surface GetPanelSurface()
        {
            Curve c = GetPanelProfile();
            if (c.Equals(null)) { c.Dispose(); return null; }
            Surface s = Surface.ByPatch(c);
            c.Dispose();
            return s;
        }
        /// <summary>
        /// returns solid panel based on panel profile as single Solid
        /// </summary>
        /// <returns>Solid by Polycurve sweep</returns>
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
        /// <summary>
        /// returns solid panel with holes based on flat panel as single Solid
        /// </summary>
        /// <returns>Solid by Surface Thickening</returns>
        public Solid GetPanelSolidHole()
        {
            Surface s = GetPanelSurface();
            Solid S = null;
            if (Holes.Count > 5)
            {
                Geometry[] G0 = s.Split(Holes[0]);
                Geometry[] G1 = G0[1].Split(Holes[1]);
                Geometry[] G2 = G1[1].Split(Holes[2]);
                Geometry[] G3 = G2[1].Split(Holes[3]);
                Geometry[] G4 = G3[1].Split(Holes[4]);
                Geometry[] G5 = G4[1].Split(Holes[5]);
                S = (G5[1] as Surface).Thicken(Thickness);
                G0.ForEach(g => g.Dispose());
                G1.ForEach(g => g.Dispose());
                G2.ForEach(g => g.Dispose());
                G3.ForEach(g => g.Dispose());
                G4.ForEach(g => g.Dispose());
                G5.ForEach(g => g.Dispose());
            }
            else if (Holes.Count > 3)
            {
                Geometry[] G0 = s.Split(Holes[0]);
                Geometry[] G1 = G0[1].Split(Holes[1]);
                Geometry[] G2 = G1[1].Split(Holes[2]);
                Geometry[] G3 = G2[1].Split(Holes[3]);
                S = (G3[1] as Surface).Thicken(Thickness);
                G0.ForEach(g => g.Dispose());
                G1.ForEach(g => g.Dispose());
                G2.ForEach(g => g.Dispose());
                G3.ForEach(g => g.Dispose());
            }
            else if (Holes.Count > 1)
            {
                Geometry[] G0 = s.Split(Holes[0]);
                Geometry[] G1 = G0[1].Split(Holes[1]);
                S = (G1[1] as Surface).Thicken(Thickness);
                G0.ForEach(g => g.Dispose());
                G1.ForEach(g => g.Dispose());
            }
            s.Dispose();
            return S;
        }
        /// <summary>
        /// returns edge labels at edge on triangular mesh face backside as an Polycurve Array
        /// </summary>
        /// <param name="Scale">Scale = (letter height)/4; in mesh units</param>
        /// <returns>Polycurve Array</returns>
        public PolyCurve[] GetEdgeLabels(double Scale)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < Triangle.E.Count; j++)
            {
                if (Triangle.E[j].Edge.E.Count == 1) continue;
                Point m = Point.ByCoordinates(ArcPoints[j][2].X / 2 + ArcPoints[(j + 1) % 3][0].X / 2, ArcPoints[j][2].Y / 2 + ArcPoints[(j + 1) % 3][0].Y / 2, ArcPoints[j][2].Z / 2 + ArcPoints[(j + 1) % 3][0].Z / 2);
                Vector Y = Triangle.VertexVectors[j][0].Cross(Triangle.Normal);
                Word w = Word.ByStringOriginVectors(Triangle.E[j].Edge.Name, m, Triangle.VertexVectors[j][0], Y);
                labels.AddRange(w.display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            return labels.ToArray();
        }
        /// <summary>
        /// returns edge labels at edge connector and triangle label at triangle center 
        /// on triangular mesh face backside as an Polycurve Array
        /// </summary>
        /// <param name="Scale">Scale = (letter height)/4; in mesh units</param>
        /// <param name="PanelSystem">parent PanelSystem of Panel</param>
        /// <returns>Polycurve Array</returns>
        public PolyCurve[] GetLabels(double Scale, PanelSystem PanelSystem)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < Triangle.E.Count; j++)
            {
                if (!PanelSystem.E.ContainsKey(Triangle.E[j].Edge)) continue;
                Point m = Point.ByCoordinates(ArcPoints[j][2].X / 2 + ArcPoints[(j + 1) % 3][0].X / 2, ArcPoints[j][2].Y / 2 + ArcPoints[(j + 1) % 3][0].Y / 2, ArcPoints[j][2].Z / 2 + ArcPoints[(j + 1) % 3][0].Z / 2);
                Vector Y = Triangle.VertexVectors[j][0].Cross(Triangle.Normal);
                m = m.Subtract(Y.Scale(PanelSystem.E[Triangle.Edges[j]][0].Inset + 2 * PanelSystem.E[Triangle.Edges[j]][0].Width));
                Word w = Word.ByStringOriginVectors(Triangle.E[j].Edge.Name, m, Triangle.VertexVectors[j][0], Y);
                labels.AddRange(w.display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            Point c = Triangle.Center;
            Vector N = Triangle.VertexVectors[0][0].Cross(Triangle.Normal);
            Word W = Word.ByStringOriginVectors(Triangle.Name, c, Triangle.VertexVectors[0][0], N);
            labels.AddRange(W.display(2 * Scale));
            c.Dispose(); N.Dispose(); W.Dispose();
            return labels.ToArray();
        }
        /// <summary>
        /// returns triangle label at triangle center on triangular mesh face frontside as an Polycurve Array
        /// </summary>
        /// <param name="Scale"></param>
        /// <returns></returns>
        public PolyCurve[] GetTriangleLabel(double Scale)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            Point c = Triangle.Center;
            Vector X = Triangle.VertexVectors[0][0].Reverse();
            Vector Y = Triangle.VertexVectors[0][0].Cross(Triangle.Normal);
            Word W = Word.ByStringOriginVectors(Triangle.Name, c, X, Y);
            labels.AddRange(W.display(Scale));
            c.Dispose(); X.Dispose(); Y.Dispose(); W.Dispose();
            return labels.ToArray();
        }
        /// <summary>
        /// add hole with center at closest point on panel to input point with given radius
        /// </summary>
        /// <param name="Point">Point is center of hole</param>
        /// <param name="Radius">Radius of hole</param>
        public void AddHole(Point Point, double Radius)
        {
            Surface s = GetPanelSurface();
            Point p = s.ClosestPointTo(Point);
            Holes.Add(Circle.ByCenterPointRadiusNormal(p, Radius, Triangle.Normal));
            s.Dispose();
            p.Dispose();
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
                if (!ArcPoints.Equals(null)) ArcPoints.ForEach(p => p.ForEach(pt => pt.Dispose()));
                if (!Holes.Equals(null)) Holes.ForEach(h => h.Dispose());
            }
            disposed = true;
        }
    }

    /// <summary>
    /// TrianglePanelEdgeBased: extends TrianglePanel so that edge offset is face-based and corner fillets are distance based
    /// </summary>
    public class TrianglePanelFaceBased : TrianglePanel
    {
        //**PROPERTIES
        /// <summary>
        /// gets offset of panel edge from triangle face edge
        /// </summary>
        public double EdgeOffset { get; set; }

        //**CONSTRUCTOR
        internal TrianglePanelFaceBased(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerOffset, double MinCornerRadius, int Direction)
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
                           (EdgeOffset + MinCornerRadius) / Math.Sin(Triangle.Angles[0]/2) - MinCornerRadius,
                           (EdgeOffset + MinCornerRadius) / Math.Sin(Triangle.Angles[1]/2) - MinCornerRadius,
                           (EdgeOffset + MinCornerRadius) / Math.Sin(Triangle.Angles[2]/2) - MinCornerRadius
                       };
            for (int i = 0; i < r.Length; i++) if (r[i] < CornerOffset) r[i] = CornerOffset;
            // corner arcs based on triangle vertex
            List<Point[]> P = new List<Point[]>(3);
            for (int i = 0; i < 3; i++)
            {
                double rV = r[i] * Math.Tan(Triangle.Angles[i] / 2) - EdgeOffset / Math.Cos(Triangle.Angles[i] / 2);
                Point a1 = Triangle.VertexPoints[i].Add(Triangle.VertexVectors[i][2].Scale(r[i]));
                Point a0 = a1.Add(Triangle.VertexVectors[i][1].Scale(rV)).Subtract(Triangle.VertexVectors[i][3].Scale(rV));
                Point a2 = a1.Add(Triangle.VertexVectors[i][0].Scale(rV)).Add(Triangle.VertexVectors[i][3].Scale(rV));
                Point[] arc0 = new Point[] { a0, a1, a2 };
                P.Add(arc0);
            }
            ArcPoints = P.ToArray();
        }

        //**METHOD**CREATE
        /// <summary>
        /// creates TrianglePanel with single edge offset determined by greatest edge offset determined by angle at edge, with corner fillets based on offset from corner
        /// </summary>
        /// <param name="Triangle">reference triangular face on mesh</param>
        /// <param name="Thickness">thickness of panel in mesh units</param>
        /// <param name="MinEdgeOffset">minimum edge offset</param>
        /// <param name="CornerOffset">offset of corner fillet from corner</param>
        /// <param name="MinCornerRadius">minimum radius of corner fillet</param>
        /// <param name="Direction">direction of panel extrusion into a solid</param>
        /// <returns>TrianglePanel</returns>
        public static TrianglePanelFaceBased ByMeshFace(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerOffset, double MinCornerRadius, int Direction = 0)
        { return new TrianglePanelFaceBased(Triangle, Thickness, MinEdgeOffset, CornerOffset, MinCornerRadius, Direction); }
    }

    /// <summary>
    /// TrianglePanelEdgeBased: extends TrianglePanel so that edge offset is face-based and corner fillets are distance based
    /// </summary>
    public class TrianglePanelEdgeBased : TrianglePanel
    {
        //**PROPERTIES**QUERY
        /// <summary>
        /// gets edge offsets in mesh face order starting at vertex 0
        /// </summary>
        public double[] EdgeOffset { get; set; }

        //**CONSTRUCTOR
        internal TrianglePanelEdgeBased(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerRadius, int Direction)
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
                Point a0 = Triangle.VertexPoints[i].Add(Triangle.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA + CornerRadius / tanB)).Add(Triangle.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA));
                Point a1 = Triangle.VertexPoints[i].Add(Triangle.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA)).Add(Triangle.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA)).Add(Triangle.VertexVectors[i][2].Scale(CornerRadius / sinB - CornerRadius));
                Point a2 = Triangle.VertexPoints[i].Add(Triangle.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA + CornerRadius / tanB)).Add(Triangle.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA));
                Point[] arc = new Point[] { a0, a1, a2 };
                P.Add(arc);
            }
            ArcPoints = P.ToArray();
        }

        //**METHOD**CREATE
        /// <summary>
        /// creates TrianglePanel with unique edge offsets determined by angle at edge, with given corner radius
        /// </summary>
        /// <param name="Triangle">reference triangular face on mesh</param>
        /// <param name="Thickness">thickness of panel in mesh units</param>
        /// <param name="MinEdgeOffset">minimum edge offset</param>
        /// <param name="CornerRadius">radius of corner fillet</param>
        /// <param name="Direction">direction of panel extrusion into a solid</param>
        /// <returns>TrianglePanel</returns>
        public static TrianglePanelEdgeBased ByMeshFace(Triangle Triangle, double Thickness, double MinEdgeOffset, double CornerRadius, int Direction = 0)
        { return new TrianglePanelEdgeBased(Triangle, Thickness, MinEdgeOffset, CornerRadius, Direction); }
    }

    /// <summary>
    /// EdgeConnector: geometry wrapper for connectors based on halfedges at edge that returns connectors at given point
    /// (connections are on backside of mesh except when there are an odd number of edges, ie there is branching in the mesh)
    /// </summary>
    public class EdgeConnector : IDisposable
    {
        //**FIELDS
        internal double Width;
        bool disposed = false;

        //**PROPERTIES
        /// <summary>
        /// gets reference edge in mesh
        /// </summary>
        public Topology.Edge Edge { get; private set; }
        /// <summary>
        /// gets halfedges of faces being connected by connector
        /// </summary>
        public HalfEdge[] HalfEdges { get; private set; }
        /// <summary>
        /// get distance from mesh edge to panel edge
        /// </summary>
        public double Inset { get; private set; }
        /// <summary>
        /// get angle between the faces being connected
        /// </summary>
        public double InsetAngle { get; private set; }
        /// <summary>
        /// get profile points as List of Vectors from given point on edge
        /// </summary>
        public List<Vector> Profile { get; private set; }
        /// <summary>
        /// get Vector Array {X1,Y1,Z1,X2,Y2,Z2,N} 
        /// (halfedge1 CS),(halfedge1 CS),edge normal
        /// </summary>
        public Vector[] Vectors { get; private set; }
        /// <summary>
        /// gets pocket information for dowel nuts
        /// </summary>
        public List<Circle> Pockets { get; set; }

        internal EdgeConnector(HalfEdge e1, HalfEdge e2, double Width, double PanelThickness, double PanelMinOffset)
        {
            // initialize properties
            Profile = new List<Vector>();
            Pockets = new List<Circle>();
            HalfEdges = new HalfEdge[] { e1, e2 };
            Edge = e1.Edge;
            this.Width = Width;
            // initial index for edge normal
            int i = 0;
            // determine HalfEdge index ( 0 - 1 - 2)
            int i1 = Edge.E.IndexOf(e1);
            int i2 = Edge.E.IndexOf(e2);
            // do not create connector in the case where HalfEdge index is (0,2)
            if (i1 + i2 != 2)
            {
                // generate orthonormal coordinate system based on right-hand rule
                // X: normal vector of face (front face normal)
                // Y: orthonormal vector on face plane (towards center)
                // Z: tangent vector of edge (halfedge direction)
                Vector X1 = e1.Face.Normal;
                Vector Z1 = e1.GetVector().Normalized();
                Vector Y1 = X1.Cross(Z1).Normalized();
                Vector X2 = e2.Face.Normal;
                Vector Z2 = e2.GetVector().Normalized();
                Vector Y2 = X2.Cross(Z2).Normalized();
                // case where there are three half edges to an edge
                // and where it is a front face to a back face connection
                // ie. halfedge1 and halfedge2
                if (i1 + i2 == 3)
                {
                    // edge normal is set to 1 which references normal between halfedge1 and halfedge2 
                    i = 1;
                    // change orthonormal coordinate system to front face based for halfedge with index 1
                    if (Edge.E.IndexOf(e1) == 1) { X1 = X1.Reverse(); Z1 = Z1.Reverse(); }
                    if (Edge.E.IndexOf(e2) == 1) { X2 = X2.Reverse(); Z2 = Z2.Reverse(); }
                }
                // set initial connector angle for inset calculation
                InsetAngle = Edge.Angle[0];
                // case halfedge count is two, edge normal array has one element that is average of face normals
                // case halfedge count is three, edge normal array has inverse average of face normals (0,1), (1,2), (2,0)
                // vector eN is edge normal vector pertinent to connector
                Vector eN = Edge.Normal[i].Normalized();
                // case where halfedge count is three
                if (Edge.E.Count > 2)
                {
                    // set inset calculation based on smallest angle pertinent to edgeconnectors at edge
                    InsetAngle = Math.Min(Edge.Angle[0], Edge.Angle[1]);
                    // reverse normal
                    if (i == 0) eN = eN.Reverse();
                }
                // calculate inset distance based on inset angle and panel offsets based on inset angle
                double EdgeOffset = Math.Max(0.5 * PanelThickness / Math.Tan(InsetAngle * Math.PI / 360), PanelMinOffset / Math.Sin(InsetAngle * Math.PI / 360) - PanelThickness / Math.Tan(InsetAngle * Math.PI / 360));
                Inset = Math.Max(Width / 2 + EdgeOffset, (PanelThickness / 2 + Width) / Math.Tan(InsetAngle * Math.PI / 360));
                // initialize Vector Lists to store affine geometry information
                List<Vector> P1 = new List<Vector>();
                List<Vector> P2 = new List<Vector>();
                // v0: point closest to edge
                P1.Add(eN.Scale(-PanelThickness / 2 / Math.Sin(Edge.Angle[i] * Math.PI / 360)));
                // v1: farthest point in first halfedge face direction (inset + 1.5 spacing)
                P1.Add(X1.Scale(-PanelThickness / 2).Add(Y1.Normalized().Scale(Inset + 1.5 * Width)));
                // v2: point to create square edge with face and beginning of fillet arc
                P1.Add(P1[1].Add(X1.Scale(-Width / 2)));
                // v3: center of fillet arc
                P1.Add(P1[2].Add(Y1.Normalized().Scale(-Width / 2)));
                // v4: end of fillet arc;
                // if angle is concave ie less than 180, connector is pie shaped and normal at arc end is edge normal
                if (Edge.Angle[i] < 180) P1.Add(P1[3].Add(eN.Scale(-Width / 2)));
                // if angle is convex ie greater than 180, connector is boomerang shaped and normal at arc end is halfedge face normal
                // v5: point where two sides of the connector meet; line (v0,v5) is line of symmetry
                else
                {
                    P1.Add(P1[3].Add(X1.Scale(-Width / 2)));
                    P1.Add(P1[0].Add(eN.Scale(-Width / Math.Sin(Edge.Angle[i] * Math.PI / 360))));
                }
                // w0: farthest point in second halfedge face direction (inset + 1.5 spacing)
                P2.Add(X2.Scale(-PanelThickness / 2).Add(Y2.Scale(Inset + 1.5 * Width)));
                // w1: point to create square edge with face and beginning of fillet arc
                P2.Add(P2[0].Add(X2.Scale(-Width / 2)));
                // w2: center of fillet arc
                P2.Add(P2[1].Add(Y2.Normalized().Scale(-Width / 2)));
                // w3: end of fillet arc
                // if angle is concave ie less than 180, connector is pie shaped and normal at arc end is edge normal
                if (Edge.Angle[i] < 180) P2.Add(P2[2].Add(eN.Scale(-Width / 2)));
                // if angle is convex ie greater than 180, connector is boomerang shaped and normal at arc end is halfedge face normal
                else P2.Add(P2[2].Add(X2.Scale(-Width / 2)));
                
                // reverse order of vectors generated by halfedge2 orthonormal coordinate system and combine vectors
                // to create a list of vectors in the direction from halfedge1 to halfedge2
                P2.Reverse();
                Profile.AddRange(P1);
                Profile.AddRange(P2);
                // store orthonormal coordinates in properties
                Vectors = new Vector[] { X1, Y1, Z1, X2, Y2, Z2, eN };
            }
        }
        internal EdgeConnector(HalfEdge e1, double Width, double PanelThickness, double Thickness)
        {
            Profile = new List<Vector>();
            Edge = e1.Edge;
        }


        //**METHOD**CREATE
        /// <summary>
        /// returns edge connector given halfedges, connector width, panel thickness
        /// </summary>
        /// <param name="HalfEdge1">halfedge of first face to connect</param>
        /// <param name="HalfEdge2">halfedge of second face to connect</param>
        /// <param name="Width">width, thickness, and hole spacing of edge connector</param>
        /// <param name="PanelThickness">thickness of panels that are being connected</param>
        /// <param name="PanelMinOffset">minimum distance of first hole to panel edge</param>
        /// <returns>EdgeConnector</returns>
        public static EdgeConnector ByHalfEdge(HalfEdge HalfEdge1, HalfEdge HalfEdge2, double Width, double PanelThickness, double PanelMinOffset) { return new EdgeConnector(HalfEdge1, HalfEdge2, Width, PanelThickness, PanelMinOffset); }
        /// <summary>
        /// returns edge connector given edges, connector width, panel thickness
        /// </summary>
        /// <param name="e">edge to place connector(s)</param>
        /// <param name="Width">width, thickness, and hole spacing of edge connector</param>
        /// <param name="PanelThickness">thickness of panels that are being connected</param>
        /// <param name="PanelMinOffset">minimum distance of first hole to panel edge</param>
        /// <returns>EdgeConnector Array</returns>
        public static EdgeConnector[] ByEdge(Topology.Edge e, double Width, double PanelThickness, double PanelMinOffset)
        {
            if (e.E.Count < 2) return null;
            if (e.E.Count == 2) return new EdgeConnector[] { new EdgeConnector(e.E[0], e.E[1], Width, PanelThickness, PanelMinOffset) };
            return new EdgeConnector[]{ new EdgeConnector(e.E[0], e.E[1], Width, PanelThickness, PanelMinOffset),
                                         new EdgeConnector(e.E[1], e.E[2], Width, PanelThickness, PanelMinOffset) };
        }


        //**METHOD**ACTIONS
        /// <summary>
        /// returns profile curve of panel as single Polycurve
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <returns>Polycurve of joined Arcs and Lines</returns>
        public PolyCurve GetConnectorProfile(Point Point)
        {
            // return null for faulty connectors
            if (!(Profile.Count > 8)) return null;
            // initialize and generate Point List using EdgeConnector geometry data at given point
            List<Point> Points = new List<Point>(Profile.Count);
            Profile.ForEach(v => Points.Add(Point.Add(v)));
            // initialize Curves for Profile
            Curve[] Curves;
            // generate curves for when angle is concave ie less than 180
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
            // generate curves for when angle is convex or greater than 180
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
            // create closed polycurve from curves
            PolyCurve Result = PolyCurve.ByJoinedCurves(Curves);
            // dispose unmanaged resources
            Points.ForEach(p => p.Dispose()); Curves.ForEach(c => c.Dispose());
            return Result;
        }
        /// <summary>
        /// returns flat panel on mesh face as single Surface
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <returns>Surface by Polycurve patch</returns>
        public Surface GetConnectorSurface(Point Point)
        {
            Curve c = GetConnectorProfile(Point);
            if (c.Equals(null)) { c.Dispose(); return null; }
            Surface s = Surface.ByPatch(c);
            c.Dispose();
            return s;
        }
        /// <summary>
        /// returns solid panel based on panel profile as single Solid
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <returns>Solid by Polycurve sweep</returns>
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
        /// <summary>
        /// returns solid panel with holes based on flat panel as single Solid
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <returns>Solid by Surface Thickening</returns>
        public Solid GetConnectorSolidHoles(Point Point)
        {
            Surface s = GetConnectorSurface(Point);
            Solid S = null;
            if (Pockets.Count > 3)
            {
                Geometry[] G0 = s.Split(Pockets[0]);
                Geometry[] G1 = G0[1].Split(Pockets[1]);
                Geometry[] G2 = G1[1].Split(Pockets[2]);
                Geometry[] G3 = G2[1].Split(Pockets[3]);
                S = (G3[1] as Surface).Thicken(Width, true);
                G0.ForEach(g => g.Dispose());
                G1.ForEach(g => g.Dispose());
                G2.ForEach(g => g.Dispose());
                G3.ForEach(g => g.Dispose());
            }
            else if (Pockets.Count > 1)
            {
                Geometry[] G0 = s.Split(Pockets[0]);
                Geometry[] G1 = G0[1].Split(Pockets[1]);
                S = (G1[1] as Surface).Thicken(Width, true);
                G0.ForEach(g => g.Dispose());
                G1.ForEach(g => g.Dispose());
            }
            s.Dispose();
            return S;
        }
        /// <summary>
        /// returns edge labels at edge on triangular mesh face backside as an Polycurve Array
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <param name="Scale">Scale = (letter height)/4; Scale is factor of width</param>
        /// <returns>Polycurve Array</returns>
        public PolyCurve[] GetEdgeLabel(Point Point, double Scale = 1/32)
        {
            if (!(Profile.Count > 8)) return null;
            Point p4 = Point.Add(Profile[4]);
            Point p5 = Point.Add(Profile[5]);
            Point pt = p5;
            if (Profile.Count == 9) pt = Point.ByCoordinates(p4.X / 2 + p5.X / 2, p4.Y / 2 + p5.Y / 2, p4.Z / 2 + p5.Z / 2);
            Vector X = Vectors[3].Subtract(Vectors[0]);
            Vector Y = Vectors[6].Reverse();
            if (InsetAngle > 180)
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
        /// <summary>
        /// returns coordinate system for transformation of edge
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <returns>CoordinateSystem</returns>
        public CoordinateSystem GetCS(Point Point)
        {
            if (!(Profile.Count > 8)) return null;
            Point p4 = Point.Add(Profile[4]);
            Point p5 = Point.Add(Profile[5]);
            Point pt = p5;
            if (Profile.Count == 9) pt = Point.ByCoordinates(p4.X / 2 + p5.X / 2, p4.Y / 2 + p5.Y / 2, p4.Z / 2 + p5.Z / 2);
            Vector X = Vectors[3].Subtract(Vectors[0]).Normalized();
            Vector Y = Vectors[6].Reverse().Normalized();
            if (InsetAngle > 180)
            {
                pt = Point.Add(Profile[0]);
                X = X.Reverse();
                Y = Y.Reverse();
            }
            CoordinateSystem CS = CoordinateSystem.ByOriginVectors(pt, X, Y);
            p4.Dispose(); p5.Dispose(); pt.Dispose();
            X.Dispose(); Y.Dispose();
            return CS;
        }
        /// <summary>
        /// add pockets based on edge condition and inset with given radius
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <param name="Radius">Radius of pocket</param>
        public void AddPockets(Point Point, double Radius)
        {
            int i = Edge.E.IndexOf(HalfEdges[0]) + Edge.E.IndexOf(HalfEdges[1]);
            int j = 0;
            if (i > 2) j = 1;
            double a = Edge.Angle[j];
            Point p1a = Point.Add(Profile[3]);
            Point p1b = p1a.Subtract(Vectors[1].Normalized().Scale(Width));
            Point p2a = Point.Add(Profile[Profile.Count - 3]);
            Point p2b = p2a.Subtract(Vectors[4].Normalized().Scale(Width));
            if (a == InsetAngle)
            {
                Pockets.Add(Circle.ByCenterPointRadiusNormal(p1a, Radius, Vectors[5]));
                Pockets.Add(Circle.ByCenterPointRadiusNormal(p1b, Radius, Vectors[5]));
                Pockets.Add(Circle.ByCenterPointRadiusNormal(p2a, Radius, Vectors[5]));
                Pockets.Add(Circle.ByCenterPointRadiusNormal(p2b, Radius, Vectors[5]));
            }
            else if (j == 1)
            {
                Pockets.Add(Circle.ByCenterPointRadiusNormal(p2a, Radius, Vectors[5]));
                Pockets.Add(Circle.ByCenterPointRadiusNormal(p2b, Radius, Vectors[5]));
            }
            else if (j == 0)
            {
                Pockets.Add(Circle.ByCenterPointRadiusNormal(p1a, Radius, Vectors[5]));
                Pockets.Add(Circle.ByCenterPointRadiusNormal(p1b, Radius, Vectors[5]));
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
                if (!Pockets.Equals(null)) Pockets.ForEach(h => h.Dispose());
            }
            disposed = true;
        }
    }

    /// <summary>
    /// PanelSystem: Mesh based system of EdgeConnectors and TrianglePanels
    /// </summary>
    public class PanelSystem : IDisposable
    {
        //**FIELDS
        bool disposed = false;
        internal Dictionary<Triangle, TrianglePanel> T;
        internal Dictionary<Topology.Edge, EdgeConnector[]> E;
        internal TriangleMesh M;

        //**PROPERTIES**QUERY
        /// <summary>
        /// gets TrianglePanel Array based on mesh
        /// </summary>
        public TrianglePanel[] Panels { get { return T.Values.ToArray(); } }
        /// <summary>
        /// gets two dimensional EdgeConnector Array indexed by mesh edge then by halfedge normals;
        /// each edge can have more than one EdgeConnector
        /// </summary>
        public EdgeConnector[][] Connectors { get { return E.Values.ToArray(); } }
        /// <summary>
        /// gets EdgeConnector List of mesh
        /// </summary>
        public List<EdgeConnector> ConnectorList
        {
            get
            {
                List<EdgeConnector> C = new List<EdgeConnector>();
                Connectors.ForEach(c => C.AddRange(c));
                return C;
            }
        }

        //**CONSTRUCTOR
        internal PanelSystem(TriangleMesh Mesh) { M = Mesh; }
        internal PanelSystem(TriangleMesh Mesh, double Width, double Thickness, double MinEdgeOffset, double CornerRadius, double HoleRadius, double PocketRadius)
            : this(Mesh)
        { GenerateEdgeConnectors(Width, Thickness, MinEdgeOffset, PocketRadius); GenerateTrianglePanels(Thickness, MinEdgeOffset, CornerRadius, HoleRadius); }

        //**METHODS**CREATE
        /// <summary>
        /// creates PanelSystem based on given Mesh and Panel parameters
        /// </summary>
        /// <param name="Mesh">BecauseWeDynamo Mesh</param>
        /// <param name="Width">Width, thickness and spacing of EdgeConnectors</param>
        /// <param name="Thickness">Thickness of TrianglePanels</param>
        /// <param name="MinEdgeOffset">Minimum edge offset from mesh face of TrianglePanels</param>
        /// <param name="CornerRadius">Radius of fillet for TrianglePanels</param>
        /// <param name="HoleRadius">Radius of holes for TrianglePanels</param>
        /// <param name="PocketRadius">Radius of pockets for EdgeConnectors</param>
        /// <returns>PanelSystem</returns>
        public static PanelSystem ByMesh(TriangleMesh Mesh, double Width, double Thickness, double MinEdgeOffset, double CornerRadius, double HoleRadius, double PocketRadius)
        { return new PanelSystem(Mesh, Width, Thickness, MinEdgeOffset, CornerRadius, HoleRadius, PocketRadius); }

        //**METHODS**ACTIONS
        /// <summary>
        /// generates new Edge/EdgeConnector Dictionary with given parameters
        /// </summary>
        /// <param name="Width">Width, thickness and spacing of EdgeConnectors</param>
        /// <param name="PanelThickness">Thickness of TrianglePanels</param>
        /// <param name="PanelMinOffset">Minimum edge offset from mesh face of TrianglePanels</param>
        /// <param name="PocketRadius">Radius of pockets for EdgeConnectors</param>
        /// <returns>EdgeConnector Array from new Dictionary values</returns>
        public EdgeConnector[][] GenerateEdgeConnectors(double Width, double PanelThickness, double PanelMinOffset, double PocketRadius)
        {
            // initialize mesh edge to edge connector dictionary
            E = new Dictionary<Topology.Edge, EdgeConnector[]>(M.E2.Count + M.E3.Count);
            // generate EdgeConnector for edges with two halfedges
            M.E2.ForEach(e => E.Add(e, new EdgeConnector[] { new EdgeConnector(e.E[0], e.E[1], Width, PanelThickness, PanelMinOffset) }));
            // generate EdgeConnector for edges with three halfedges
            M.E3.ForEach(e => E.Add(e, new EdgeConnector[]{ 
                new EdgeConnector(e.E[0], e.E[1], Width, PanelThickness, PanelMinOffset),
                new EdgeConnector(e.E[1], e.E[2], Width, PanelThickness, PanelMinOffset)}));
            // add pockets to EdgeConnectors
            Connectors.ForEach(C => C.ForEach(c => c.AddPockets(c.Edge.MidPoint, PocketRadius)));
            return Connectors;
        }
        /// <summary>
        /// generates new Triangle/TrianglePanel Dictionary with given parameters based on EdgeConnectors
        /// </summary>
        /// <param name="Thickness">Thickness of TrianglePanels</param>
        /// <param name="MinEdgeOffset">Minimum edge offset from mesh face of TrianglePanels</param>
        /// <param name="CornerRadius">Radius of fillet for TrianglePanels</param>
        /// <param name="HoleRadius">Radius of holes for TrianglePanels</param>
        /// <returns>TrianglePanel Array from new DIctionary values</returns>
        public TrianglePanel[] GenerateTrianglePanels(double Thickness, double MinEdgeOffset, double CornerRadius, double HoleRadius)
        {
            if (!(E.Count > 0)) return null;
            // initialize mesh face to panel dictionary
            T = new Dictionary<Triangle, TrianglePanel>(M.Faces.Count);
            // iterate through mesh faces to generate panels
            for (int i = 0; i < M.Faces.Count; i++)
            {
                // generate triangle panel with input values at given face
                TrianglePanel t = TrianglePanelEdgeBased.ByMeshFace(M.Faces[i], Thickness, MinEdgeOffset, CornerRadius);
                // add face/panel pair to dictionary
                T.Add(M.Faces[i], t);
                // iterate through edges and add holes based on edge connector calculations
                for (int j = 0; j < M.Faces[i].Edges.Length; j++)
                    if (E.ContainsKey(M.Faces[i].Edges[j]))
                    {
                        // midpoint of edge
                        Point p = M.Faces[i].Edges[j].MidPoint;
                        // orthogonal vector to edge that is on the plane of triangle face with direction towards center of face
                        Vector Y = M.Faces[i].Normal.Cross(M.Faces[i].E[j].GetVector()).Normalized();
                        // calculate first hole position based on edge using midpoint as origin
                        // distance between hole center and edge is equal to inset property of edge
                        p = p.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Inset));
                        // add hole to list of holes of TrianglePanel
                        t.AddHole(p, HoleRadius);
                        // calculate second hole position based on edge using midpoint as origin
                        // distance between hole center and edge is equal to the sum of the inset property of edge and the width property of edge
                        p = p.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Width));
                        // add hole to list of holes of TrianglePanel
                        t.AddHole(p, HoleRadius);
                        //dispose of unmanaged resources
                        p.Dispose(); Y.Dispose();
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
