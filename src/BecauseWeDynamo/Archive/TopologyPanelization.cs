using System;
using System.Collections.Generic;
using System.Linq;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Runtime;
using Fabrication.Text;
using Topology;
using Fabrication;

namespace Topology.Panelization
{
    /// <summary>
    /// panel object
    /// </summary>
    public class Panel : IDisposable
    {
        internal bool disposed = false;
        internal double[] Thickness;
        internal double BevelAngle;

        //**PROPERTIES**
        /// <summary>
        /// mesh face that defines base plane and panel geometry
        /// </summary>
        public Face Face { get; private set; }
        /// <summary>
        /// context coordinate system inherited from associated Mesh.Face
        /// </summary>
        public CoordinateSystem CS { get { return Face.CS; } }
        /// <summary>
        /// point array that defines panel geometry
        /// </summary>
        public Point[][] ArcPoints { get; set; }
        /// <summary>
        /// double array that stores edge offsets
        /// </summary>
        public double[] EdgeOffset { get; set; }

        //**CONSTRUCTOR
        internal Panel(Face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle)
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
                double sinA = Math.Sin(VectorMath.toRadians(Face.Angles[i]));
                double sinB = Math.Sin(VectorMath.toRadians(Face.Angles[i]/2));
                double tanB = Math.Abs(Math.Tan(VectorMath.toRadians(Face.Angles[i]/2)));
                double CornerRadius = MinCornerRadius * tanB;
                int j = (i + Face.E.Count - 1) % Face.E.Count;
                Point[] arc = new Point[] { Face.VertexPoints[i].Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA)).Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA)) };
                if (Face.Angles[i] < 180)
                {
                    Point a0 = arc[0].Add(Face.VertexVectors[i][1].Scale(CornerRadius / tanB));
                    Point a1 = arc[0].Add(Face.VertexVectors[i][2].Scale(CornerRadius / sinB - CornerRadius));
                    Point a2 = arc[0].Add(Face.VertexVectors[i][0].Scale(CornerRadius / tanB));
                    arc[0].Dispose();
                    arc = new Point[] { a0, a1, a2 };
                }
                P.Add(arc);
            }
            ArcPoints = P.ToArray();
        }

        //**CREATE
        /// <summary>
        /// creates panel system based on mesh face and parameters
        /// </summary>
        /// <param name="Face"></param>
        /// <param name="ThicknessFront"></param>
        /// <param name="ThicknessBack"></param>
        /// <param name="MinEdgeOffset"></param>
        /// <param name="CornerOffset"></param>
        /// <param name="MinFaceAngle"></param>
        /// <param name="BevelAngle"></param>
        /// <returns></returns>
        public static Panel ByFaceAndParameters(Face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double MinFaceAngle, double BevelAngle = 0)
        { return new Panel(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, MinFaceAngle, BevelAngle); }

        /// <summary>
        /// creates panel profile at mesh face as a polycurve
        /// </summary>
        /// <returns>panel profiles as polycurves</returns>
        public PolyCurve GetPanelProfile()
        {
            if (ArcPoints.Equals(null)) return null;
            List<Curve> Curves = new List<Curve>(2 * Face.E.Count);
            for (int i = 0; i < Face.E.Count; i++)
            {
                int j = (i + 1) % Face.E.Count;
                if (ArcPoints[i].Length > 1) Curves.Add(Arc.ByThreePoints(ArcPoints[i][0], ArcPoints[i][1], ArcPoints[i][2]));
                Curves.Add(Line.ByStartPointEndPoint(ArcPoints[i][ArcPoints[i].Length - 1], ArcPoints[j][0]));
            }
            PolyCurve C = PolyCurve.ByJoinedCurves(Curves);
            if (!C.IsPlanar)
            {
                Plane P = Plane.ByOriginNormal(Face.Center, Face.Normal);
                C = C.PullOntoPlane(P) as PolyCurve;
                P.Dispose();
            }
            Curves.ForEach(c => c.Dispose());
            return C;
        }
        /// <summary>
        /// creates panel face at mesh face as a surface
        /// </summary>
        /// <returns></returns>
        public Surface GetPanelSurface(PolyCurve c)
        {
            if (c.Equals(null)) { return null; }
            return Surface.ByPatch(c);
        }
        /// <summary>
        /// creates panel face extrusion
        /// </summary>
        /// <returns>panel face extrusion</returns>
        public Solid GetPanelSolid(Surface s) { return s.Thicken(Thickness[0],false); }
        /// <summary>
        /// creates edge labels as polylines on back of panel
        /// </summary>
        /// <param name="Scale">Scale</param>
        /// <param name="Offset">Label Offset from Edge</param>
        /// <param name="LabelPrefix">Label Prefix</param>
        /// <returns>Edge Labels on Panel Back</returns>
        public PolyCurve[] GetEdgeLabelsBack(double Scale = 1.0/12, double Offset = 0, string LabelPrefix = "")
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int i = 0; i < Face.E.Count; i++)
            {
                if (Face.E[i].Edge.E.Count == 1) continue;
                int j = (i + 1) % Face.E.Count;
                Point m = VectorMath.Center(ArcPoints[i][ArcPoints[i].Length - 1], ArcPoints[j][0]);
                Vector Y = Face.VertexVectors[i][0].Cross(Face.Normal);
                if (Offset > 0) m = m.Add(Y.Normalized().Scale(-Offset));
                Word w = Word.ByStringOriginVectors(LabelPrefix + Face.E[i].Edge.Name, m, Face.VertexVectors[i][0], Y);
                labels.AddRange(w.Display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            return labels.ToArray();
        }
        /// <summary>
        /// creates edge labels as polylines on front of panel
        /// </summary>
        /// <param name="Scale">Scale</param>
        /// <param name="Offset">Label Offset from Edge</param>
        /// <param name="LabelPrefix">Label Prefix</param>
        /// <returns>Edge Labels on Panel Front</returns>
        public PolyCurve[] GetEdgeLabelsFront(double Scale = 1.0/12, double Offset = 0, string LabelPrefix = "")
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int i = 0; i < Face.E.Count; i++)
            {
                if (Face.E[i].Edge.E.Count == 1) continue;
                int j = (i + 1) % Face.E.Count;
                Point m = VectorMath.Center(ArcPoints[i][ArcPoints[i].Length - 1], ArcPoints[j][0]);
                Vector Y = Face.VertexVectors[i][0].Cross(Face.Normal);
                if (Offset > 0) m = m.Add(Y.Normalized().Scale(-Offset));
                Word w = Word.ByStringOriginVectors(LabelPrefix + Face.E[i].Edge.Name, m, Face.VertexVectors[i][0].Reverse(), Y);
                labels.AddRange(w.Display(Scale));
                m.Dispose(); Y.Dispose(); w.Dispose();
            }
            return labels.ToArray();
        }
        /// <summary>
        /// creates edge and face labels as polylines on back of panel
        /// </summary>
        /// <param name="Scale">Scale</param>
        /// <param name="Offset">Label Offset from Edge</param>
        /// <param name="LabelPrefix">Label Prefix</param>
        /// <returns>Labels on Panel Back</returns>
        public PolyCurve[] GetLabelsBack(double Scale = 1.0/12, double Offset = 0, string LabelPrefix = "")
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            labels.AddRange(GetEdgeLabelsBack(Scale, Offset, LabelPrefix));
            labels.AddRange(GetFaceLabelBack(Scale, LabelPrefix));
            return labels.ToArray();
        }
        /// <summary>
        /// creates edge and face labels as polylines on front of panel
        /// </summary>
        /// <param name="Scale">Scale</param>
        /// <param name="Offset">Label Offset from Edge</param>
        /// <param name="LabelPrefix">Label Prefix</param>
        /// <returns>Labels on Panel Front</returns>
        public PolyCurve[] GetLabelsFront(double Scale = 1.0/12, double Offset = 0, string LabelPrefix = "")
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            labels.AddRange(GetEdgeLabelsFront(Scale, Offset, LabelPrefix));
            labels.AddRange(GetFaceLabelFront(Scale, LabelPrefix));
            return labels.ToArray();
        }
        /// <summary>
        /// creates face labels as polylines on front of panel
        /// </summary>
        /// <param name="Scale">Scale</param>
        /// <param name="LabelPrefix">Label Prefix</param>
        /// <returns>Face Labels on Panel Front</returns>
        public PolyCurve[] GetFaceLabelFront(double Scale = 1.0/12, string LabelPrefix = "")
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            Point c = Face.Center;
            Vector X = Face.VertexVectors[0][0].Reverse();
            Vector Y = Face.VertexVectors[0][0].Cross(Face.Normal);
            Word W = Word.ByStringOriginVectors(LabelPrefix + Face.Name, c, X, Y);
            labels.AddRange(W.Display(Scale));
            c.Dispose(); X.Dispose(); Y.Dispose(); W.Dispose();
            return labels.ToArray();
        }
        /// <summary>
        /// creates face labels as polylines on back of panel
        /// </summary>
        /// <param name="Scale">Scale</param>
        /// <param name="LabelPrefix">Label Prefix</param>
        /// <returns>Face Labels on Panel Back</returns>
        public PolyCurve[] GetFaceLabelBack(double Scale = 1.0/12, string LabelPrefix = "")
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            Point c = Face.Center;
            Vector N = Face.VertexVectors[0][0].Cross(Face.Normal);
            Word W = Word.ByStringOriginVectors(LabelPrefix + Face.Name, c, Face.VertexVectors[0][0], N);
            labels.AddRange(W.Display(Scale));
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

    /// <summary>
    /// panel object with holes that correspond to edge conector placement
    /// bevels are set to zero for hole addition
    /// </summary>
    public class PanelHole : Panel
    {
        /// <summary>
        /// double array list of hole geometry
        /// </summary>
        public List<double[]> Holes { get; set; }

        //**METHODS**CONSTRUCTOR
        internal PanelHole(Face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle)
            : base(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, MinCornerRadius, MinFaceAngle, BevelAngle)
        {
            // initialize
            Holes = new List<double[]>();
        }
        public static PanelHole ByFaceAndParameters(Face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle = 0)
        { return new PanelHole(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, MinCornerRadius, MinFaceAngle, BevelAngle); }

        //**MEDTHODS**ACTIONS
        /// <summary>
        /// get panel hole geometry data
        /// </summary>
        /// <returns></returns>
        public Circle[] GetPanelHoleProfile()
        {
            List<Circle> result = new List<Circle>();
            for (int i = 0; i < Holes.Count; i++)
            {
                Point center = Point.ByCoordinates(Holes[i][0], Holes[i][1], Holes[i][2]);
                result.Add(Circle.ByCenterPointRadiusNormal(center, Holes[i][3], Face.Normal));
                center.Dispose();
            }
            return result.ToArray();
        }
        /// <summary>
        /// get panel hole geometry data as polycurves
        /// </summary>
        /// <returns></returns>
        public PolyCurve[] GetPanelHoleProfileAsPolyCurve()
        {
            List<PolyCurve> result = new List<PolyCurve>();
            for (int i = 0; i < Holes.Count; i++)
            {
                Point center = Point.ByCoordinates(Holes[i][0], Holes[i][1], Holes[i][2]);
                Circle c = Circle.ByCenterPointRadiusNormal(center, Holes[i][3], Face.Normal);
                result.Add(PolyCurve.ByJoinedCurves(c.Patch().PerimeterCurves()));
                c.Dispose(); center.Dispose();
            }
            return result.ToArray();
        }
        /// <summary>
        /// returns solid panel with holes based on flat panel as single Solid
        /// </summary>
        /// <returns>Solid by Surface Thickening</returns>
        public Surface GetPanelSurfaceHole()
        {
            Surface s = GetPanelSurface(GetPanelProfile());
            List<Surface> S = new List<Surface>(){s};
            for (int i = 0; i < Holes.Count; i++)
            {
                Point center = Point.ByCoordinates(Holes[i][0], Holes[i][1], Holes[i][2]);
                Circle c = Circle.ByCenterPointRadiusNormal(center, Holes[i][3], Face.Normal);
                S.Add(S[S.Count-1].Trim(c,center)[0] as Surface);
                c.Dispose(); center.Dispose();
            }
            Surface result = S[S.Count - 1];
            for (int i = 0; i < S.Count - 1; i++) if (!result.Equals(S[i])) S[i].Dispose();
            S = null;
            s.Dispose();
            return result;
        }
        /// <summary>
        /// get panel solid hole
        /// </summary>
        /// <returns></returns>
        public Solid GetPanelSolidHole(Surface s)
        {
            Solid S = (s.Translate(Face.Normal, -Thickness[1]) as Surface).Thicken(Thickness[0] + Thickness[1], false);
            return S;
        }
        /// <summary>
        /// add hole with center at closest point on panel to input point with given radius
        /// </summary>
        /// <param name="Point">Point is center of hole</param>
        /// <param name="Radius">Radius of hole</param>
        public void AddHole(Point Point, double Radius)
        {
            Surface s = GetPanelSurface(GetPanelProfile());
            Point p = s.ClosestPointTo(Point);
            Holes.Add(new double[]{p.X,p.Y,p.Z,Radius});
            s.Dispose();
            p.Dispose();
        }

    }

    /// <summary>
    /// panel object with minimum edge offset determined by Z-axis elevation
    /// </summary>
    public class PanelGradient : Panel
    {
        internal PanelGradient(Face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle, double OffsetBaseZ)
            : base(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, MinCornerRadius, MinFaceAngle, BevelAngle)
        {
            for (int i = 0; i < Face.E.Count; i++)
            {
                double Offset = OffsetBaseZ * Face.E[i].Edge.MidPoint.Z;
                // check clearance in back
                if (EdgeOffset[i] < Offset) EdgeOffset[i] = Offset;
            }
            // corner arcs based on triangle vertex
            List<Point[]> P = new List<Point[]>(Face.E.Count);
            for (int i = 0; i < Face.E.Count; i++)
            {
                double sinA = Math.Sin(Face.Angles[i] * Math.PI / 180);
                double sinB = Math.Sin(Face.Angles[i] * Math.PI / 360);
                double tanB = Math.Abs(Math.Tan(Face.Angles[i] * Math.PI / 360));
                double CornerRadius = MinCornerRadius * tanB;
                int j = (i + Face.E.Count - 1) % Face.E.Count;
                Point[] arc = new Point[] { Face.VertexPoints[i].Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA)).Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA)) };
                if (Face.Angles[i] < 180)
                {
                    Point a0 = arc[0].Add(Face.VertexVectors[i][1].Scale(CornerRadius / tanB));
                    Point a1 = arc[0].Add(Face.VertexVectors[i][2].Scale(CornerRadius / sinB - CornerRadius));
                    Point a2 = arc[0].Add(Face.VertexVectors[i][0].Scale(CornerRadius / tanB));
                    arc[0].Dispose();
                    arc = new Point[] { a0, a1, a2 };
                }
                P.Add(arc);
            }
            ArcPoints = P.ToArray();
        }
        //**CREATE
        /// <summary>
        /// creates panel object with minimum edge offset determined by Z-axis elevation
        /// default bevel angle and z-axis contribution to offset is set to 0
        /// </summary>
        /// <param name="Face"></param>
        /// <param name="ThicknessFront"></param>
        /// <param name="ThicknessBack"></param>
        /// <param name="MinEdgeOffset"></param>
        /// <param name="MinCornerRadius"></param>
        /// <param name="MinFaceAngle"></param>
        /// <param name="BevelAngle"></param>
        /// <param name="OffsetBaseZ"></param>
        /// <returns></returns>
        public static PanelGradient ByFaceAndParameters(Face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle =0, double OffsetBaseZ = 0)
        { return new PanelGradient(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, MinCornerRadius, MinFaceAngle, BevelAngle, OffsetBaseZ); }

    }

    /// <summary>
    /// connectors at edges created for model
    /// </summary>
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
        public double Inset { get; set; }
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
        /// midpoint of edge
        /// </summary>
        public Point Midpoint { get { return Edge.MidPoint; } }

        internal EdgeConnector(HalfEdge e1, HalfEdge e2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        {
            // initialize properties
            this.BevelAngle = BevelAngle;
            Profile = new List<Vector>();
            HalfEdges = new HalfEdge[] { e1, e2 };
            Edge = e1.Edge;
            Dim = new Dictionary<char, double>();
            Dim.Add('H', Height); Dim.Add('D', Depth); Dim.Add('S', Spacing);
            double R = Math.Min(CornerRadius, Dim['H'] / 2);
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
                X1 = Z1.Cross(Y1).Normalized();
                Vector X2 = e2.Face.Normal;
                Vector Z2 = e2.GetVector().Normalized();
                Vector Y2 = X2.Cross(Z2).Normalized();
                X2 = Z2.Cross(Y2).Normalized();
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
                Vector eY = eN.Cross(Z1).Normalized();
                eN = Z1.Cross(eY).Normalized();
                eY.Dispose();
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
                P1.Add(P1[1].Add(X1.Scale(CornerRadius - Dim['H'])));
                // v3: center of fillet arc
                P1.Add(P1[2].Add(Y1.Scale(-CornerRadius)));
                // v4: end of fillet arc;
                P1.Add(P1[3].Add(X1.Scale(-CornerRadius)));
                // v5: point where two sides of the connector meet; line (v0,v5) is line of symmetry
                P1.Add(P1[0].Add(eN.Scale(-Dim['H'] / Math.Sin(Edge.Angle[iAngle] * Math.PI / 360))));

                // w0: farthest point in second halfedge face direction (inset + 1.5 spacing)
                P2.Add(X2.Scale(-ThicknessBack).Add(Y2.Scale(Inset + Dim['S'] + Dim['H'] / 2)));
                // w1: point to create square edge with face and beginning of fillet arc
                P2.Add(P2[0].Add(X2.Scale(CornerRadius - Dim['H'])));
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
        /// get connector profile at midpoint
        /// </summary>
        /// <returns></returns>
        public PolyCurve GetConnectorProfileAtMidpoint()
        {
            // return null for faulty connectors
            if (!(Profile.Count > 8)) return null;
            // initialize and generate Point List using EdgeConnector geometry data at given point
            List<Point> Points = new List<Point>(Profile.Count);
            Profile.ForEach(v => Points.Add(Edge.MidPoint.Add(v)));
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
        /// returns profile curve of panel as single Polycurve
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <returns>Polycurve of joined Arcs and Lines</returns>
        public PolyCurve GetConnectorProfileAtPoint(Point Point)
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
            Point pt = Point.Add(Profile[5]);
            Vector X = Vectors[3].Subtract(Vectors[0]);
            Vector Y = Vectors[6].Reverse();
            if (Edge.Angle[iAngle] > 179.999 && Edge.Angle[iAngle] < 180.001) X = Profile[1].Subtract(Profile[0]);
            if (Edge.Angle[iAngle] > 180) { pt = Point.Add(Profile[0]); X = X.Reverse(); Y = Y.Reverse(); }
            CoordinateSystem cs = GetCS(Point);
            Vector Z1 = cs.ZAxis.Normalized();
            Vector Z2 = X.Cross(Y).Normalized();
            if (!Z1.IsAlmostEqualTo(Z2)) X = X.Reverse();
            Word w = Word.ByStringOriginVectors(LabelPrefix + Edge.Name, pt, X.Normalized(), Y.Normalized());
            List<PolyCurve> label = w.Display(Scale);
            pt.Dispose();
            X.Dispose(); Y.Dispose(); Z1.Dispose(); Z2.Dispose();
            cs.Dispose(); w.Dispose();
            return label.ToArray();
        }
        public CoordinateSystem GetCS(Point Point)
        {
            Surface s = GetConnectorSurface(GetConnectorProfileAtPoint(Point));
            CoordinateSystem cs = s.CoordinateSystemAtParameter(0.5, 0.5);
            s.Dispose();
            return cs;
        }


        //**MODS
        /// <summary>
        /// get connector profile bevel at point
        /// </summary>
        /// <param name="Point"></param>
        /// <returns></returns>
        public PolyCurve GetConnectorProfileBevel(Point Point)
        {
            // return default for angles smaller than 
            if (Edge.Angle[iAngle] <= 2 * (90 - BevelAngle)) return GetConnectorProfileAtPoint(Point);
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
            if (Edge.Angle[iAngle] == 2 * (180 - BevelAngle)) { }
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
            for (int j = 7; j < Points.Count; j++) Curves.Add(Line.ByStartPointEndPoint(Points[j], Points[(j + 1) % Points.Count]));
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

    public class EdgeConnectorGradient : EdgeConnector
    {
        //**FIELDS

        internal EdgeConnectorGradient(HalfEdge e1, HalfEdge e2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double Z)
            : base(e1, e2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle)
        {
            int i1 = Edge.E.IndexOf(e1);
            int i2 = Edge.E.IndexOf(e2);
            // do not create connector in the case where HalfEdge index is (0,2)
            if (i1 + i2 != 2)
            {
                if (EdgeOffset < Edge.MidPoint.Z * Z) EdgeOffset = Edge.MidPoint.Z * Z;
                Inset = Math.Max(Dim['H'] / 2 + EdgeOffset, (ThicknessBack + Dim['H']) / Math.Tan(InsetAngle * Math.PI / 360));
                // initialize Vector Lists to store affine geometry information
                List<Vector> P1 = new List<Vector>();
                List<Vector> P2 = new List<Vector>();
                // v0: point closest to edge
                P1.Add(Vectors[6].Scale(-ThicknessBack / Math.Sin(Edge.Angle[iAngle] * Math.PI / 360)));
                // v1: farthest point in first halfedge face direction (inset + 1.5 spacing)
                P1.Add(Vectors[0].Scale(-ThicknessBack).Add(Vectors[1].Scale(Inset + Dim['S'] + Dim['H'] / 2)));
                // v2: point to create square edge with face and beginning of fillet arc
                P1.Add(P1[1].Add(Vectors[0].Scale(CornerRadius - Dim['H'])));
                // v3: center of fillet arc
                P1.Add(P1[2].Add(Vectors[1].Scale(-CornerRadius)));
                // v4: end of fillet arc;
                P1.Add(P1[3].Add(Vectors[0].Scale(-CornerRadius)));
                // v5: point where two sides of the connector meet; line (v0,v5) is line of symmetry
                P1.Add(P1[0].Add(Vectors[6].Scale(-Dim['H'] / Math.Sin(Edge.Angle[iAngle] * Math.PI / 360))));

                // w0: farthest point in second halfedge face direction (inset + 1.5 spacing)
                P2.Add(Vectors[3].Scale(-ThicknessBack).Add(Vectors[4].Scale(Inset + Dim['S'] + Dim['H'] / 2)));
                // w1: point to create square edge with face and beginning of fillet arc
                P2.Add(P2[0].Add(Vectors[3].Scale(CornerRadius - Dim['H'])));
                // w2: center of fillet arc
                P2.Add(P2[1].Add(Vectors[4].Scale(-CornerRadius)));
                // w3: end of fillet arc
                P2.Add(P2[2].Add(Vectors[3].Scale(-CornerRadius)));

                // reverse order of vectors generated by halfedge2 orthonormal coordinate system and combine vectors
                // to create a list of vectors in the direction from halfedge1 to halfedge2
                P2.Reverse();
                Profile.Clear();
                Profile.AddRange(P1);
                Profile.AddRange(P2);
            }

        }

        //**METHOD**CREATE
        /// <summary>
        /// 
        /// </summary>
        /// <param name="HalfEdge1"></param>
        /// <param name="HalfEdge2"></param>
        /// <param name="Height"></param>
        /// <param name="Depth"></param>
        /// <param name="Spacing"></param>
        /// <param name="ThicknessFront"></param>
        /// <param name="ThicknessBack"></param>
        /// <param name="PanelMinOffset"></param>
        /// <param name="CornerRadius"></param>
        /// <param name="OffsetBaseZ"></param>
        /// <returns></returns>
        public static EdgeConnectorGradient ByHalfEdge(HalfEdge HalfEdge1, HalfEdge HalfEdge2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double OffsetBaseZ)
        { return new EdgeConnectorGradient(HalfEdge1, HalfEdge2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, OffsetBaseZ); }
        public static EdgeConnectorGradient ByEdge(Topology.Edge e, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double OffsetBaseZ)
        {
            if (e.E.Count < 2) return null;
            return new EdgeConnectorGradient(e.E[0], e.E[1], Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, OffsetBaseZ);
        }
    }

    public class EdgeConnectorHole : EdgeConnector
    {
        //**FIELDS
        double Radius;

        internal EdgeConnectorHole(HalfEdge e1, HalfEdge e2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double PocketRadius)
            : base(e1, e2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle)
        {
            Radius = PocketRadius;
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
        public static EdgeConnectorHole ByHalfEdge(HalfEdge HalfEdge1, HalfEdge HalfEdge2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double PocketRadius) 
        { return new EdgeConnectorHole(HalfEdge1, HalfEdge2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, PocketRadius); }
        /// <summary>
        /// returns edge connector given edges, connector width, panel thickness
        /// </summary>
        /// <param name="e">edge to place connector(s)</param>
        /// <param name="Width">width, thickness, and hole spacing of edge connector</param>
        /// <param name="PanelThickness">thickness of panels that are being connected</param>
        /// <param name="PanelMinOffset">minimum distance of first hole to panel edge</param>
        /// <returns>EdgeConnector Array</returns>
        public static EdgeConnectorHole[] ByEdge(Topology.Edge e, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double PocketRadius)
        {
            if (e.E.Count < 2) return null;
            if (e.E.Count == 2) return new EdgeConnectorHole[] { new EdgeConnectorHole(e.E[0], e.E[1], Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, PocketRadius) };
            return new EdgeConnectorHole[]{ new EdgeConnectorHole(e.E[0], e.E[1], Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, PocketRadius),
                                         new EdgeConnectorHole(e.E[1], e.E[2], Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, PocketRadius) };
        }
        /// <summary>
        /// returns conector sufrace with holes at midpoint
        /// </summary>
        /// <returns></returns>
        public Surface GetConnectorSurfaceHoles()
        {
            Surface s = GetConnectorSurface(GetConnectorProfileAtMidpoint());
            Circle[] Pockets = GeneratePockets(Edge.MidPoint, Radius);
            List<Surface> S = new List<Surface>() { s };
            for (int i = 0; i < Pockets.Length; i++)  S.Add(S[S.Count - 1].Trim(Pockets[i], Pockets[i].CenterPoint)[0] as Surface);
            Surface result = S[S.Count - 1]; 
            for (int i = 0; i < S.Count - 1; i++) if (!result.Equals(S[i])) S[i].Dispose();
            Pockets.ForEach(c => c.Dispose());
            s.Dispose();
            S = null;
            return result;
        }
        /// <summary>
        /// returns solid panel with holes based on flat panel as single Solid
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <returns>Solid by Surface Thickening</returns>
        public Solid GetConnectorSolidHoles()
        {
            Surface s = GetConnectorSurfaceHoles();
            Solid S = s.Thicken(Dim['D'], true);
            s.Dispose();
            return S;
        }
        
        /// <summary>
        /// add pockets based on edge condition and inset with given radius
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <param name="Radius">Radius of pocket</param>
        public Circle[] GeneratePockets(Point Point, double Radius)
        {
            List<Circle> Pockets = new List<Circle>();
            int i = Edge.E.IndexOf(HalfEdges[0]) + Edge.E.IndexOf(HalfEdges[1]);
            int j = 0;
            if (i > 2) j = 1;
            double a = Edge.Angle[j];
            Point p1a = Point.Add(Profile[3]);
            Point p1b = p1a.Subtract(Vectors[1].Normalized().Scale(Dim['D']));
            Point p2a = Point.Add(Profile[Profile.Count - 3]);
            Point p2b = p2a.Subtract(Vectors[4].Normalized().Scale(Dim['D']));
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
            return Pockets.ToArray();
        }
    }

    /// <summary>
    /// panel system for mesh objects
    /// </summary>
    public class PanelSystem : IDisposable
    {
        //**FIELDS
        bool disposed = false;
        internal Dictionary<Face, Panel> F;
        internal Dictionary<Edge, EdgeConnector[]> E;
        internal Mesh M;

        //**PROPERTIES**QUERY
        /// <summary>
        /// panel array
        /// </summary>
        public Panel[] Panels { get { return F.Values.ToArray(); } }
        /// <summary>
        /// edgeconnector array
        /// </summary>
        public EdgeConnector[] Connectors
        {
            get
            {
                List<EdgeConnector> result = new List<EdgeConnector>();
                for (int i=0; i< E.Count; i++)
                {
                    result.AddRange(E.Values.ToArray()[i]);
                }
            return result.ToArray(); } }

        //**CONSTRUCTOR
        internal PanelSystem(Mesh Mesh) { M = Mesh; }
        internal PanelSystem(Mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double ConnectorFilletRadius, double BevelAngle)
            : this(Mesh)
        { GenerateEdgeConnectors(Height, Depth, Spacing, ThicknessFront, ThicknessBack, MinEdgeOffset, ConnectorFilletRadius, BevelAngle); GeneratePanels(ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, BevelAngle); }

        //**METHODS**CREATE
        /// <summary>
        /// returns default empty panel system with associated mesh
        /// </summary>
        /// <param name="Mesh">Mesh Object</param>
        /// <returns>Panel System</returns>
        public static PanelSystem ByMesh(Mesh Mesh) { return new PanelSystem(Mesh); }
        /// <summary>
        /// creates panel system with clearance bevels
        /// </summary>
        /// <param name="Mesh">Mesh Geometry</param>
        /// <param name="Height">Connector Height</param>
        /// <param name="Depth">Connector Depth</param>
        /// <param name="Spacing">Connector Hole Spacing</param>
        /// <param name="ThicknessFront">Panel Thickness Normal Direction</param>
        /// <param name="ThicknessBack">Panel Thickness Other Direction</param>
        /// <param name="MinEdgeOffset">Panel Minimum Edge Offset</param>
        /// <param name="CornerOffset">Panel Minimum Corner Offset</param>
        /// <param name="ConnectorFilletRadius">Connector Fillet Radius</param>
        /// <param name="BevelAngle">Panel Clearance Bevel Angle</param>
        /// <returns>Panel System with Clearance Bevels</returns>
        public static PanelSystem ByMeshParameters(Mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double ConnectorFilletRadius, double BevelAngle)
        {
            PanelSystem P = new PanelSystem(Mesh);
            P.GenerateEdgeConnectors(Height, Depth, Spacing, ThicknessFront, ThicknessBack, MinEdgeOffset, ConnectorFilletRadius, BevelAngle); 
            P.GeneratePanels(ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, BevelAngle);
            return P;
        }
        /// <summary>
        /// creates panel system with edge offset based on elevation
        /// </summary>
        /// <param name="Mesh">Mesh Geometry</param>
        /// <param name="Height">Connector Height</param>
        /// <param name="Depth">Connector Depth</param>
        /// <param name="Spacing">Connector Hole Spacing</param>
        /// <param name="ThicknessFront">Panel Thickness Normal Direction</param>
        /// <param name="ThicknessBack">Panel Thickness Other Direction</param>
        /// <param name="MinEdgeOffset">Panel Minimum Edge Offset</param>
        /// <param name="CornerOffset">Panel Minimum Corner Offset</param>
        /// <param name="ConnectorFilletRadius">Connector Fillet Radius</param>
        /// <param name="OffsetBaseZ">Panel Edge Elevation-Based Offset</param>
        /// <returns>Panel System with Elevation-Based Offsets</returns>
        public static PanelSystem ByMeshParametersGradient(Mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double ConnectorFilletRadius, double BevelAngle, double OffsetBaseZ)
        {
            PanelSystem P = new PanelSystem(Mesh);
            P.GenerateEdgeConnectorsGradient(Height, Depth, Spacing, ThicknessFront, ThicknessBack, MinEdgeOffset, ConnectorFilletRadius, BevelAngle, OffsetBaseZ);
            P.GeneratePanelsGradient(ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, BevelAngle, OffsetBaseZ);
            return P;
        }
        /// <summary>
        /// creates panel system with holes for bolt to dowel nut connection
        /// </summary>
        /// <param name="Mesh"></param>
        /// <param name="Height"></param>
        /// <param name="Depth"></param>
        /// <param name="Spacing"></param>
        /// <param name="ThicknessFront"></param>
        /// <param name="ThicknessBack"></param>
        /// <param name="MinEdgeOffset"></param>
        /// <param name="CornerOffset"></param>
        /// <param name="ConnectorFilletRadius"></param>
        /// <param name="BevelAngle"></param>
        /// <param name="PanelHoleRadius"></param>
        /// <param name="ConnectorPocketRadius"></param>
        /// <param name="NumConnectors"></param>
        /// <returns></returns>
        public static PanelSystem ByMeshParametersHole(Mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double ConnectorFilletRadius, double BevelAngle, double PanelHoleRadius, double ConnectorPocketRadius, int NumConnectors)
        {
            PanelSystem P = new PanelSystem(Mesh);
            P.GenerateEdgeConnectorsHole(Height, Depth, Spacing, ThicknessFront, ThicknessBack, MinEdgeOffset, ConnectorFilletRadius, BevelAngle, ConnectorPocketRadius);
            P.GeneratePanelsHole(ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, BevelAngle, PanelHoleRadius, NumConnectors);
            return P;
        }

        //**METHODS**ACTIONS
        /// <summary>
        /// generate edgeconnector geometry based on system inputs
        /// </summary>
        internal void GenerateEdgeConnectors(double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double ConnectorFilletRadius, double BevelAngle)
        {
            // initialize mesh edge to edge connector dictionary
            E = new Dictionary<Topology.Edge, EdgeConnector[]>(M.E2.Count + M.E3.Count);
            // generate EdgeConnector
            M.E2.ForEach(e => E.Add(e, new EdgeConnector[] { EdgeConnector.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle) }));
            M.E3.ForEach(e => E.Add(e, new EdgeConnector[] { EdgeConnector.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle) }));
        }
        /// <summary>
        /// generate edgeconnector geometry based on system inputs
        /// </summary>
        internal void GenerateEdgeConnectorsGradient(double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double ConnectorFilletRadius, double BevelAngle, double OffsetBaseZ)
        {
            // initialize mesh edge to edge connector dictionary
            E = new Dictionary<Topology.Edge, EdgeConnector[]>(M.E2.Count + M.E3.Count);
            // generate EdgeConnector
            M.E2.ForEach(e => E.Add(e, new EdgeConnector[] { EdgeConnectorGradient.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle, OffsetBaseZ) }));
            M.E3.ForEach(e => E.Add(e, new EdgeConnector[] { EdgeConnectorGradient.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle, OffsetBaseZ) }));
        }
        /// <summary>
        /// generates new Edge/EdgeConnector Dictionary with given parameters
        /// </summary>
        internal void GenerateEdgeConnectorsHole(double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double ConnectorFilletRadius, double BevelAngle, double ConnectorPocketRadius)
        {
            // initialize mesh edge to edge connector dictionary
            E = new Dictionary<Topology.Edge, EdgeConnector[]>(M.E2.Count + M.E3.Count);
            // generate EdgeConnector for edges with two halfedges
            M.E2.ForEach(e => E.Add(e, EdgeConnectorHole.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle, ConnectorPocketRadius)));
            // generate EdgeConnector for edges with three halfedges
            M.E3.ForEach(e => E.Add(e, EdgeConnectorHole.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle, ConnectorPocketRadius)));
        }
        /// <summary>
        /// generate panels based on system inputs
        /// </summary>
        /// <param name="ThicknessFront"></param>
        /// <param name="ThicknessBack"></param>
        /// <param name="MinEdgeOffset"></param>
        /// <param name="CornerOffset"></param>
        /// <param name="BevelAngle"></param>
        /// <returns></returns>
        internal bool GeneratePanels(double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double BevelAngle)
        {
            if (!(E.Count > 0)) return false;
            // initialize mesh face to panel dictionary
            F = new Dictionary<Face, Panel>(M.Faces.Count);
            // iterate through mesh faces to generate panels
            for (int i = 0; i < M.Faces.Count; i++)
            {
                // generate triangle panel with input values at given face
                Panel panel = Panel.ByFaceAndParameters(M.Faces[i], ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, M.MinFaceAngle, BevelAngle);
                // add face/panel pair to dictionary
                F.Add(M.Faces[i], panel);
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
        /// <summary>
        /// generate panels based on system inputs
        /// </summary>
        /// <param name="ThicknessFront"></param>
        /// <param name="ThicknessBack"></param>
        /// <param name="MinEdgeOffset"></param>
        /// <param name="CornerOffset"></param>
        /// <param name="OffsetBaseZ"></param>
        /// <returns></returns>
        internal bool GeneratePanelsGradient(double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double BevelAngle, double OffsetBaseZ)
        {
            if (!(E.Count > 0)) return false;
            // initialize mesh face to panel dictionary
            F = new Dictionary<Face, Panel>(M.Faces.Count);
            // iterate through mesh faces to generate panels
            for (int i = 0; i < M.Faces.Count; i++)
            {
                // generate triangle panel with input values at given face
                PanelGradient panel = PanelGradient.ByFaceAndParameters(M.Faces[i], ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, M.MinFaceAngle, BevelAngle, OffsetBaseZ);
                // add face/panel pair to dictionary
                F.Add(M.Faces[i], panel);
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
        /// <summary>
        /// generates new Triangle/TrianglePanel Dictionary with given parameters based on EdgeConnectors
        /// </summary>
        internal bool GeneratePanelsHole(double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double BevelAngle, double HoleRadius, int NumConnectors)
        {
            if (!(E.Count > 0)) return false;
            // initialize mesh face to panel dictionary
            F = new Dictionary<Face, Panel>(M.Faces.Count);
            // iterate through mesh faces to generate panels
            for (int i = 0; i < M.Faces.Count; i++)
            {
                // generate triangle panel with input values at given face
                PanelHole panel = PanelHole.ByFaceAndParameters(M.Faces[i], ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, M.MinFaceAngle, BevelAngle);
                // add face/panel pair to dictionary
                F.Add(M.Faces[i], panel);
                // iterate through edges and add holes based on edge connector calculations
                for (int j = 0; j < M.Faces[i].Edges.Length; j++)
                    if (E.ContainsKey(M.Faces[i].Edges[j]))
                    {
                        // orthogonal vector to edge that is on the plane of triangle face with direction towards center of face
                            Vector Y = M.Faces[i].Normal.Cross(M.Faces[i].E[j].GetVector()).Normalized();
                        if (NumConnectors.Equals(1))
                        {;
                            // calculate first hole position based on edge using midpoint as origin
                            // distance between hole center and edge is equal to inset property of edge
                            Point p = M.Faces[i].Edges[j].MidPoint.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Inset));
                            // add hole to list of holes of TrianglePanel
                            panel.AddHole(p, HoleRadius);
                            // calculate second hole position based on edge using midpoint as origin
                            // distance between hole center and edge is equal to the sum of the inset property of edge and the width property of edge
                            p = p.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Dim['S']));
                            // add hole to list of holes of TrianglePanel
                            panel.AddHole(p, HoleRadius);
                            //dispose of unmanaged resources
                            p.Dispose();
                        }
                        else
                        {
                            for (int k=0; k < NumConnectors; k++)
                            {
                                Point p = M.Faces[i].Edges[j].HalfEdges[0].V[0].Point;
                                Vector v = M.Faces[i].Edges[j].HalfEdges[0].GetVector();
                                p = p.Add(v.Scale(v.Length*(2*k+1)/(2*NumConnectors)));
                                p = p.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Inset));
                                // add hole to list of holes of TrianglePanel
                                panel.AddHole(p, HoleRadius);
                                // calculate second hole position based on edge using midpoint as origin
                                // distance between hole center and edge is equal to the sum of the inset property of edge and the width property of edge
                                p = p.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Dim['S']));
                                // add hole to list of holes of TrianglePanel
                                panel.AddHole(p, HoleRadius);
                                //dispose of unmanaged resources
                                p.Dispose(); v.Dispose();
                            }
                        }
                        Y.Dispose();
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
