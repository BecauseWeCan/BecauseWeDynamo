using System;
using System.Collections.Generic;
using System.Linq;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Runtime;
using Text;
using Topology;

namespace Topology.Panelization
{
    public class Panel : IDisposable
    {
        internal bool disposed = false;
        internal double[] Thickness;
        internal double BevelAngle;

        //**PROPERTIES**
        public Face Face { get; private set; }
        public Point[][] ArcPoints { get; set; }
        public double[] EdgeOffset { get; set; }

        //**CONSTRUCTOR
        internal Panel(Face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerRadius, double BevelAngle)
        {
            // initialize
            this.Face = Face;
            this.Thickness = new double[] { ThicknessFront, ThicknessBack };
            this.BevelAngle = BevelAngle;
            // edge offsets indexed by triangle halfedge angles
            EdgeOffset = new double[Face.E.Count];
            for (int i = 0; i < Face.E.Count; i++)
            {
                // check clearance in back
                EdgeOffset[i] = MinEdgeOffset;
                if (Face.E[i].Angle == 360 || Face.E[i].Angle == 180) continue;
                EdgeOffset[i] = 0;
                if (Face.E[i].Angle < 2 * (90 - BevelAngle))
                    EdgeOffset[i] = ThicknessBack / Math.Tan(Face.E[i].Angle * Math.PI / 360) - ThicknessBack * Math.Tan(BevelAngle * Math.PI / 180);
                // check clearance in front
                double OffsetAngle = MinEdgeOffset / Math.Sin(Face.E[i].Angle * Math.PI / 360) - ThicknessFront / Math.Tan(Face.E[i].Angle * Math.PI / 360);
                if (EdgeOffset[i] < OffsetAngle) EdgeOffset[i] = OffsetAngle;
            }
            // corner arcs based on triangle vertex
            List<Point[]> P = new List<Point[]>(Face.E.Count);
            for (int i = 0; i < Face.E.Count; i++)
            {
                double sinA = Math.Sin(Face.Angles[i] * Math.PI / 180);
                double sinB = Math.Sin(Face.Angles[i] * Math.PI / 360);
                double tanB = Math.Tan(Face.Angles[i] * Math.PI / 360);
                int j = (i + Face.E.Count - 1) % Face.E.Count;

                Point a0 = Face.VertexPoints[i].Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA + CornerRadius / tanB)).Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA));
                Point a1 = Face.VertexPoints[i].Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA)).Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA)).Add(Face.VertexVectors[i][2].Scale(CornerRadius / sinB - CornerRadius));
                Point a2 = Face.VertexPoints[i].Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA + CornerRadius / tanB)).Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA));
                Point[] arc = new Point[] { a0, a1, a2 };
                P.Add(arc);
            }
            ArcPoints = P.ToArray();
        }

        //**CREATE
        public static Panel ByFaceAndParameters(Face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerRadius, double BevelAngle = 0)
        { return new Panel(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, CornerRadius, BevelAngle); }

        public PolyCurve GetPanelProfile()
        {
            if (ArcPoints.Equals(null)) return null;
            Curve[] Curves = new Curve[2 * Face.E.Count];
            for (int i = 0; i < Face.E.Count; i++)
            {
                int j = (i + 1) % Face.E.Count;
                Curves[2 * i] = Arc.ByThreePoints(ArcPoints[i][0], ArcPoints[i][1], ArcPoints[i][2]);
                Curves[2 * i + 1] = Line.ByStartPointEndPoint(ArcPoints[i][2], ArcPoints[j][0]);
            }
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
        public Solid[] GetPanelSolid()
        {
            Curve[] C = new Curve[4];
            C[0] = GetPanelProfile();
            if (C[0].Equals(null)) { C[0].Dispose(); return null; }
            Solid[] S = new Solid[2];
            Point[] P1 = new Point[] { C[0].StartPoint, C[0].StartPoint.Add(Face.Normal.Scale(Thickness[0])) };
            C[1] = Line.ByStartPointEndPoint(P1[0], P1[1]);
            S[0] = C[0].SweepAsSolid(C[1]);
            List<Point> P2 = new List<Point>(Face.E.Count);
            for (int i = 0; i < Face.E.Count; i++)
            {
                Arc arc = Arc.ByThreePoints(ArcPoints[i][0], ArcPoints[i][1], ArcPoints[i][2]);
                P2.Add(arc.CenterPoint.Add(Face.VertexVectors[i][4].Scale((Thickness[1] * Math.Tan(BevelAngle * Math.PI / 180) - arc.Radius) / Math.Sin(Face.Angles[i] * Math.PI / 360))));
                arc.Dispose();
            }
            C[2] = Polygon.ByPoints(P2);
            P2.Add(P2[0].Add(Face.Normal.Scale(-Thickness[1])));
            C[3] = Line.ByStartPointEndPoint(P2[0], P2[P2.Count-1]);
            S[1] = C[2].SweepAsSolid(C[3]);
            P1.ForEach(p => p.Dispose()); P2.ForEach(p => p.Dispose()); C.ForEach(p => p.Dispose());
            return S;
        }
        public PolyCurve[] GetEdgeLabels(double Scale, double Offset=0, string LabelPrefix = "")
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int i = 0; i < Face.E.Count; i++)
            {
                if (Face.E[i].Edge.E.Count == 1) continue;
                int j = (i + 1) % Face.E.Count;
                Point m = Point.ByCoordinates(ArcPoints[i][2].X / 2 + ArcPoints[j][0].X / 2, ArcPoints[i][2].Y / 2 + ArcPoints[j][0].Y / 2, ArcPoints[i][2].Z / 2 + ArcPoints[j][0].Z / 2);
                Vector Y = Face.VertexVectors[i][0].Cross(Face.Normal);
                if (Offset > 0) m = m.Add(Y.Normalized().Scale(-Offset));
                Word w = Word.ByStringOriginVectors(LabelPrefix + Face.E[i].Edge.Name, m, Face.VertexVectors[i][0], Y);
                labels.AddRange(w.display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            return labels.ToArray();
        }
        public PolyCurve[] GetLabels(double Scale, double Offset = 0, string LabelPrefix = "")
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            labels.AddRange(GetEdgeLabels(Scale,Offset,LabelPrefix));
            Point c = Face.Center;
            Vector N = Face.VertexVectors[0][0].Cross(Face.Normal);
            Word W = Word.ByStringOriginVectors(Face.Name, c, Face.VertexVectors[0][0], N);
            labels.AddRange(W.display(2 * Scale));
            c.Dispose(); N.Dispose(); W.Dispose();
            return labels.ToArray();
        }
        public PolyCurve[] GetFaceLabel(double Scale)
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
        public PolyCurve[] GetFaceLabelBack(double Scale, string LabelPrefix = "")
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            Point c = Face.Center;
            Vector N = Face.VertexVectors[0][0].Cross(Face.Normal);
            Word W = Word.ByStringOriginVectors(LabelPrefix + Face.Name, c, Face.VertexVectors[0][0], N);
            labels.AddRange(W.display(Scale));
            c.Dispose(); N.Dispose(); W.Dispose();
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

    public class PanelHole : Panel
    {
        public List<Circle> Holes { get; set; }

        //**METHODS**CONSTRUCTOR
        internal PanelHole(Face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerRadius)
            : base(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, CornerRadius, 0)
        {
            // initialize
            Holes = new List<Circle>();
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
                S = ((G5[1] as Surface).Translate(Face.Normal, -Thickness[1]) as Surface).Thicken(Thickness[0] + Thickness[1], false);
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
                S = ((G3[1] as Surface).Translate(Face.Normal, -Thickness[1]) as Surface).Thicken(Thickness[0] + Thickness[1], false);
                G0.ForEach(g => g.Dispose());
                G1.ForEach(g => g.Dispose());
                G2.ForEach(g => g.Dispose());
                G3.ForEach(g => g.Dispose());
            }
            else if (Holes.Count > 1)
            {
                Geometry[] G0 = s.Split(Holes[0]);
                Geometry[] G1 = G0[1].Split(Holes[1]);
                S = ((G1[1] as Surface).Translate(Face.Normal, -Thickness[1]) as Surface).Thicken(Thickness[0] + Thickness[1], false);
                G0.ForEach(g => g.Dispose());
                G1.ForEach(g => g.Dispose());
            }
            s.Dispose();
            return S;
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
            Holes.Add(Circle.ByCenterPointRadiusNormal(p, Radius, Face.Normal));
            s.Dispose();
            p.Dispose();
        }
    }

    public class EdgeConnector : IDisposable
    {
        //**FIELDS
        internal int iAngle;
        internal double BevelAngle;
        internal Dictionary<char, double> Dim;
        bool disposed = false;
        internal double EdgeOffset;

        //**PROPERTIES
        /// <summary>
        /// gets reference edge in mesh
        /// </summary>
        public Edge Edge { get; private set; }
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

        internal EdgeConnector(HalfEdge e1, HalfEdge e2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        {
            // initialize properties
            this.BevelAngle = BevelAngle;
            Profile = new List<Vector>();
            HalfEdges = new HalfEdge[] { e1, e2 };
            Edge = e1.Edge;
            Dim = new Dictionary<char, double>();
            Dim.Add('H', Height); Dim.Add('D', Depth); Dim.Add('S', Spacing);
            // initial index for edge normal
            iAngle = 0;
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
                    iAngle = 1;
                    // change orthonormal coordinate system to front face based for halfedge with index 1
                    if (Edge.E.IndexOf(e1) == 1) { X1 = X1.Reverse(); Z1 = Z1.Reverse(); }
                    if (Edge.E.IndexOf(e2) == 1) { X2 = X2.Reverse(); Z2 = Z2.Reverse(); }
                }
                // set initial connector angle for inset calculation
                InsetAngle = Edge.Angle[0];
                // case halfedge count is two, edge normal array has one element that is average of face normals
                // case halfedge count is three, edge normal array has inverse average of face normals (0,1), (1,2), (2,0)
                // vector eN is edge normal vector pertinent to connector
                Vector eN = Edge.Normal[iAngle].Normalized();
                // case where halfedge count is three
                if (Edge.E.Count > 2)
                {
                    // reverse normal
                    if (iAngle == 0) eN = eN.Reverse();
                }
                // calculate inset distance based on inset angle and panel offsets based on inset angle
                // check clearance in back
                EdgeOffset = 0;
                if (InsetAngle == 360 || InsetAngle == 180) EdgeOffset = PanelMinOffset;
                if (InsetAngle < 2 * (90 - BevelAngle))
                    EdgeOffset = ThicknessBack / Math.Tan(InsetAngle * Math.PI / 360) - ThicknessBack * Math.Tan(BevelAngle * Math.PI / 180);
                // check clearance in front
                double OffsetAngle = PanelMinOffset / Math.Sin(InsetAngle * Math.PI / 360) - ThicknessFront / Math.Tan(InsetAngle * Math.PI / 360);
                if (EdgeOffset < OffsetAngle) EdgeOffset = OffsetAngle;
                Inset = Math.Max(Dim['H'] / 2 + EdgeOffset, (ThicknessBack + Dim['H']) / Math.Tan(InsetAngle * Math.PI / 360));
                // initialize Vector Lists to store affine geometry information
                List<Vector> P1 = new List<Vector>();
                List<Vector> P2 = new List<Vector>();
                // v0: point closest to edge
                P1.Add(eN.Scale(-ThicknessBack / Math.Sin(Edge.Angle[iAngle] * Math.PI / 360)));
                // v1: farthest point in first halfedge face direction (inset + 1.5 spacing)
                P1.Add(X1.Scale(-ThicknessBack).Add(Y1.Scale(Inset + Dim['S'] + Dim['H'] / 2)));
                // v2: point to create square edge with face and beginning of fillet arc
                P1.Add(P1[1].Add(X1.Scale(CornerRadius-Dim['H'])));
                // v3: center of fillet arc
                P1.Add(P1[2].Add(Y1.Scale(-CornerRadius)));
                // v4: end of fillet arc;
                P1.Add(P1[3].Add(X1.Scale(-CornerRadius)));
                // v5: point where two sides of the connector meet; line (v0,v5) is line of symmetry
                P1.Add(P1[0].Add(eN.Scale(-Dim['H'] / Math.Sin(Edge.Angle[iAngle] * Math.PI / 360))));

                // w0: farthest point in second halfedge face direction (inset + 1.5 spacing)
                P2.Add(X2.Scale(-ThicknessBack).Add(Y2.Scale(Inset + Dim['S'] + Dim['H'] / 2)));
                // w1: point to create square edge with face and beginning of fillet arc
                P2.Add(P2[0].Add(X2.Scale(CornerRadius-Dim['H'])));
                // w2: center of fillet arc
                P2.Add(P2[1].Add(Y2.Scale(-CornerRadius)));
                // w3: end of fillet arc
                P2.Add(P2[2].Add(X2.Scale(-CornerRadius)));

                // reverse order of vectors generated by halfedge2 orthonormal coordinate system and combine vectors
                // to create a list of vectors in the direction from halfedge1 to halfedge2
                P2.Reverse();
                Profile.AddRange(P1);
                Profile.AddRange(P2);
                // store orthonormal coordinates in properties
                Vectors = new Vector[] { X1, Y1, Z1, X2, Y2, Z2, eN };
            }
        }

        //**METHOD**CREATE
        public static EdgeConnector ByHalfEdgeHeight(HalfEdge HalfEdge1, HalfEdge HalfEdge2, double Height, double Depth, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        { return new EdgeConnector(HalfEdge1, HalfEdge2, Height, Depth, Height, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle); }
        public static EdgeConnector ByHalfEdge(HalfEdge HalfEdge1, HalfEdge HalfEdge2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        { return new EdgeConnector(HalfEdge1, HalfEdge2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle); }
        public static EdgeConnector ByEdge(Topology.Edge e, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        {
            if (e.E.Count < 2) return null;
            return new EdgeConnector(e.E[0], e.E[1], Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle);
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
            Curve[] Curves = new Curve[] {
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
        public Surface GetConnectorSurface(PolyCurve Profile)
        {
            if (Profile.Equals(null) || !Profile.IsClosed || !Profile.IsPlanar) return null; 
            Surface s = Surface.ByPatch(Profile);
            return s;
        }
        /// <summary>
        /// returns solid panel based on panel profile as single Solid
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <returns>Solid by Polycurve sweep</returns>
        public Solid GetConnectorSolid(PolyCurve Profile)
        {
            if (Profile.Equals(null) || !Profile.IsClosed || !Profile.IsPlanar) return null; 
            Point a = Profile.StartPoint.Add(Vectors[2].Scale(Dim['D'] / 2));
            Point b = Profile.StartPoint.Add(Vectors[5].Scale(Dim['D'] / 2));
            Line l = Line.ByStartPointEndPoint(a, b);
            Solid s = Profile.SweepAsSolid(l);
            a.Dispose();
            b.Dispose();
            l.Dispose();
            return s;
        }
        /// <summary>
        /// returns edge labels at edge on triangular mesh face backside as an Polycurve Array
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <param name="Scale">Scale = (letter height)/4; Scale is in project units</param>
        /// <returns>Polycurve Array</returns>
        public PolyCurve[] GetEdgeLabel(Point Point, double Scale, string LabelPrefix = "")
        {
            if (!(Profile.Count > 8)) return null;
            Point pt = Point.Add(Profile[5]);
            Vector X = Vectors[3].Subtract(Vectors[0]);
            Vector Y = Vectors[6].Reverse();
            if (Edge.Angle[iAngle] > 179.9999 && Edge.Angle[iAngle] < 180.0001) X = Profile[1].Subtract(Profile[0]);
            if (Edge.Angle[iAngle] > 180) { pt = Point.Add(Profile[0]); X = X.Reverse(); Y = Y.Reverse(); }
            Surface s = GetConnectorSurface(GetConnectorProfile(Point));
            CoordinateSystem cs = s.CoordinateSystemAtParameter(0.5, 0.5);
            Vector Z1 = cs.ZAxis.Normalized();
            Vector Z2 = X.Cross(Y).Normalized();
            if (!Z1.IsAlmostEqualTo(Z2)) X = X.Reverse();
            Word w = Word.ByStringOriginVectors(LabelPrefix+Edge.Name, pt, X, Y);
            List<PolyCurve> label = w.display(Scale);
            pt.Dispose();
            X.Dispose(); Y.Dispose(); Z1.Dispose(); Z2.Dispose();
            s.Dispose(); cs.Dispose(); w.Dispose();
            return label.ToArray();
        }
        public CoordinateSystem GetCS(Point Point) {
            Surface s = GetConnectorSurface(GetConnectorProfile(Point));
            CoordinateSystem cs = s.CoordinateSystemAtParameter(0.5, 0.5);
            s.Dispose();
            return cs; 
        }


        //**MODS
        public PolyCurve GetConnectorProfileBevel(Point Point)
        {
            // return default for angles smaller than 
            if (Edge.Angle[iAngle] <= 2 * (90 - BevelAngle)) return GetConnectorProfile(Point);
            // initialize and generate Point List using EdgeConnector geometry data at given point
            List<Point> Points = new List<Point>(Profile.Count + 7);
            double ThicknessBack = Profile[0].Length * Math.Sin(Edge.Angle[iAngle] * Math.PI / 360);
            double BevelInset = ThicknessBack * Math.Tan(BevelAngle * Math.PI / 180);
            //double DogBoneBevel = 2 * ToolRadius * Math.Sin(Math.PI / 4 - BevelAngle * Math.PI / 2);
            //double DogBoneChord = 2 * ToolRadius / Math.Sin(Math.PI / 4 + BevelAngle * Math.PI / 2);
            for (int i = 1; i < Profile.Count; i++) Points.Add(Point.Add(Profile[i]));
            //Points.Add(Point.Add(Vectors[4].Scale(EdgeOffset + BevelInset + 2 * ToolRadius)).Add(Vectors[3].Scale(-ThicknessBack)));
            //Points.Add(Point.Add(Vectors[4].Scale(EdgeOffset + BevelInset + ToolRadius)).Add(Vectors[3].Scale(-ThicknessBack - ToolRadius)));
            Points.Add(Point.Add(Vectors[4].Scale(EdgeOffset + BevelInset)).Add(Vectors[3].Scale(-ThicknessBack)));
            if (Edge.Angle[iAngle] == 2 * (180 - BevelAngle)) {}
            else if (Edge.Angle[iAngle] < 2 * (180 - BevelAngle))
            {
                Points.Add(Point.Add(Vectors[4].Scale(EdgeOffset + BevelInset / 2)).Add(Vectors[3].Scale(-ThicknessBack / 2)));
                Points.Add(Point.Add(Vectors[1].Scale(EdgeOffset + BevelInset / 2)).Add(Vectors[0].Scale(-ThicknessBack / 2)));
            }
            else if (Edge.Angle[iAngle] > 2 * (180 - BevelAngle))
            {
                //hyp = o/(sin(e/2)*tan(b) + cos(e/2))
                Points.Add(Point.Add(Vectors[6].Scale(EdgeOffset / (Math.Sin(Edge.Angle[iAngle] / 2) * Math.Tan(BevelAngle) + Math.Cos(Edge.Angle[iAngle] / 2)))));
            }
            Points.Add(Point.Add(Vectors[1].Scale(EdgeOffset + BevelInset)).Add(Vectors[0].Scale(-ThicknessBack)));
            //Points.Add(Point.Add(Vectors[1].Scale(EdgeOffset + BevelInset + ToolRadius)).Add(Vectors[0].Scale(-ThicknessBack - ToolRadius)));
            //Points.Add(Point.Add(Vectors[1].Scale(EdgeOffset + BevelInset + 2 * ToolRadius)).Add(Vectors[0].Scale(-ThicknessBack)));
            // initialize Curves for Profile
            /*List<Curve> Curves = new List<Curve> {
                Arc.ByThreePoints(Points[Points.Count-3],Points[Points.Count-2],Points[Points.Count-1]),
                Line.ByStartPointEndPoint(Points[Points.Count-1], Points[0]),
                Line.ByStartPointEndPoint(Points[0], Points[1]),
                Arc.ByCenterPointStartPointEndPoint(Points[2],Points[1],Points[3]),
                Line.ByStartPointEndPoint(Points[3], Points[4]),
                Line.ByStartPointEndPoint(Points[4], Points[5]),
                Arc.ByCenterPointStartPointEndPoint(Points[6],Points[5],Points[7]),
                Line.ByStartPointEndPoint(Points[7],Points[8]),
                Line.ByStartPointEndPoint(Points[8],Points[9]),
                Arc.ByThreePoints(Points[9],Points[10],Points[11])};
            for (int j = 11; j < Points.Count-3; j++) Curves.Add(Line.ByStartPointEndPoint(Points[j], Points[j + 1]));*/
            List<Curve> Curves = new List<Curve> {
                Line.ByStartPointEndPoint(Points[0], Points[1]),
                Arc.ByCenterPointStartPointEndPoint(Points[2],Points[1],Points[3]),
                Line.ByStartPointEndPoint(Points[3], Points[4]),
                Line.ByStartPointEndPoint(Points[4], Points[5]),
                Arc.ByCenterPointStartPointEndPoint(Points[6],Points[5],Points[7])};
            for (int j = 7; j < Points.Count; j++) Curves.Add(Line.ByStartPointEndPoint(Points[j], Points[(j + 1)%Points.Count]));
            // create closed polycurve from curves
            PolyCurve Result = PolyCurve.ByJoinedCurves(Curves);
            // dispose unmanaged resources
            Points.ForEach(p => p.Dispose()); Curves.ForEach(c => c.Dispose());
            return Result;
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

    public class PanelSystem : IDisposable
    {
        //**FIELDS
        bool disposed = false;
        internal Dictionary<Face, Panel> F;
        internal Dictionary<Edge, EdgeConnector> E;
        internal Mesh M;

        //**PROPERTIES**QUERY
        public Panel[] Panels { get { return F.Values.ToArray(); } }
        public EdgeConnector[] Connectors { get { return E.Values.ToArray(); } }

        //**CONSTRUCTOR
        internal PanelSystem(Mesh Mesh) { M = Mesh; }
        internal PanelSystem(Mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerRadius, double BevelAngle)
            : this(Mesh)
        { GenerateEdgeConnectors(Height, Depth, Spacing, ThicknessFront, ThicknessBack, MinEdgeOffset, CornerRadius, BevelAngle); GeneratePanels(ThicknessFront, ThicknessBack, MinEdgeOffset, CornerRadius, BevelAngle); }

        //**METHODS**CREATE
        public static PanelSystem ByMesh(Mesh Mesh) { return new PanelSystem(Mesh); }
        public static PanelSystem ByMeshParameters(Mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerRadius, double BevelAngle)
        { return new PanelSystem(Mesh, Height, Depth, Spacing, ThicknessFront, ThicknessBack, MinEdgeOffset, CornerRadius, BevelAngle); }

        //**METHODS**ACTIONS
        public void GenerateEdgeConnectors(double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        {
            // initialize mesh edge to edge connector dictionary
            E = new Dictionary<Topology.Edge, EdgeConnector>(M.E2.Count + M.E3.Count);
            // generate EdgeConnector
            M.E2.ForEach(e => E.Add(e, EdgeConnector.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle)));
            M.E3.ForEach(e => E.Add(e, EdgeConnector.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle)));
        }

        public bool GeneratePanels(double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerRadius, double BevelAngle)
        {
            if (!(E.Count > 0)) return false;
            // initialize mesh face to panel dictionary
            F = new Dictionary<Face, Panel>(M.Faces.Count);
            // iterate through mesh faces to generate panels
            for (int i = 0; i < M.Faces.Count; i++)
            {
                // generate triangle panel with input values at given face
                Panel t = Panel.ByFaceAndParameters(M.Faces[i], ThicknessFront, ThicknessBack, MinEdgeOffset, CornerRadius, BevelAngle);
                // add face/panel pair to dictionary
                F.Add(M.Faces[i], t);
                // iterate through edges and add holes based on edge connector calculations
                for (int j = 0; j < M.Faces[i].Edges.Length; j++)
                    if (E.ContainsKey(M.Faces[i].Edges[j]))
                    {
                        // midpoint of edge
                        Point p = M.Faces[i].Edges[j].MidPoint;
                        // orthogonal vector to edge that is on the plane of triangle face with direction towards center of face
                        Vector Y = M.Faces[i].Normal.Cross(M.Faces[i].E[j].GetVector()).Normalized();
                    }
            }
            return true;
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
                if (!Panels.Equals(null)) Panels.ForEach(p => p.Dispose());
                if (!E.Equals(null)) E.Values.ToArray().ForEach(e => e.Dispose());
            }
            disposed = true;
        }
    }
}
