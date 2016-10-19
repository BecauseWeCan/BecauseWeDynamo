using System;
using System.Collections.Generic;
using System.Linq;
using Autodesk.DesignScript.Geometry;
using Fabrication;
using Geometry;
using Topology;
using System.Diagnostics;


namespace Panelization
{
    /// <summary>
    /// panel object
    /// </summary>
    public class panel
    {
        internal bool disposed = false;
        internal double[] Thickness;
        internal double BevelAngle;

        //**PROPERTIES**
        /// <summary>
        /// mesh face that defines base plane and panel geometry
        /// </summary>
        public face Face { get; private set; }
        /// <summary>
        /// context coordinate system inherited from associated Mesh.Face
        /// </summary>
        public CoordinateSystem CS { get { return Face.CS; } }
        /// <summary>
        /// point array that defines panel geometry
        /// </summary>
        public point[][] ArcPoints { get; set; }
        /// <summary>
        /// double array that stores edge offsets
        /// </summary>
        public double[] EdgeOffset { get; set; }

        //**CONSTRUCTOR
        internal panel(face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle)
        {
            if (!(ThicknessFront > 0)) throw new ArgumentException("ThicknessFront must be greater than zero");
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
                if (math.Mod(Face.E[i].Angle, math.PI) == 0) continue;
                EdgeOffset[i] = 0;
                if (Face.E[i].Angle < 2 * math.toRadians(90 - BevelAngle))
                    EdgeOffset[i] = ThicknessBack / Math.Tan(Face.E[i].Angle / 2) - ThicknessBack * Math.Tan(math.toRadians(BevelAngle));
                // check clearance in front
                double OffsetAngle = MinEdgeOffset / Math.Sin(Face.E[i].Angle / 2) - ThicknessFront / Math.Tan(Face.E[i].Angle / 2);
                if (EdgeOffset[i] < OffsetAngle) EdgeOffset[i] = OffsetAngle;
            }
            // corner arcs based on triangle vertex
            List<point[]> P = new List<point[]>(Face.E.Count);
            for (int i = 0; i < Face.E.Count; i++)
            {
                double sinA = Math.Sin(Face.Angles[i]);
                double sinB = Math.Sin(Face.Angles[i] / 2);
                double tanB = Math.Abs(Math.Tan(Face.Angles[i] / 2));
                double CornerRadius = MinCornerRadius * tanB;
                int j = (i + Face.E.Count - 1) % Face.E.Count;
                point[] arc = new point[] { Face.Vertices[i].Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA)).Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA)) };
                if (Face.Angles[i] < math.PI)
                {
                    point a0 = arc[0].Add(Face.VertexVectors[i][1].Scale(CornerRadius / tanB));
                    point a1 = arc[0].Add(Face.VertexVectors[i][2].Scale(CornerRadius / sinB - CornerRadius));
                    point a2 = arc[0].Add(Face.VertexVectors[i][0].Scale(CornerRadius / tanB));
                    arc = new point[] { a0, a1, a2 };
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
        public static panel ByFaceAndParameters(face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double MinFaceAngle, double BevelAngle = 0)
        { return new panel(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, MinFaceAngle, BevelAngle); }

        /// <summary>
        /// gets panel profile at mesh face as Polycurve
        /// </summary>
        public PolyCurve GetPanelProfile()
        {
            if (ArcPoints.Equals(null)) return null;
            List<Curve> Curves = new List<Curve>(2 * Face.E.Count);
            for (int i = 0; i < Face.E.Count; i++)
            {
                int j = (i + 1) % Face.E.Count;
                Point A = ArcPoints[i][0].ToPoint();
                Point D = ArcPoints[j][0].ToPoint();
                if (ArcPoints[i].Length > 1)
                {
                    Point B = ArcPoints[i][1].ToPoint();
                    Point C = ArcPoints[i][2].ToPoint();
                    Curves.Add(Arc.ByThreePoints(A, B, C));
                    Curves.Add(Line.ByStartPointEndPoint(C,D));
                    B.Dispose(); C.Dispose();
                }
                else Curves.Add(Line.ByStartPointEndPoint(A,D));
                A.Dispose(); D.Dispose();
            }
            PolyCurve Profile = PolyCurve.ByJoinedCurves(Curves);
            Curves.ForEach(c => c.Dispose());
            return Profile;
        }
        /// <summary>
        /// gets panel surface at mesh face as Surface
        /// </summary>
        public Surface GetPanelSurface()
        {
            Curve c = GetPanelProfile();
            if (c.Equals(null)) { c.Dispose(); return null; }
            Point O = Face.Center.ToPoint();
            Vector N = Face.Normal.ToVector();
            Plane P = Plane.ByOriginNormal(O, N);
            Curve M = c.PullOntoPlane(P);
            Surface s = Surface.ByPatch(M);
            c.Dispose(); M.Dispose(); N.Dispose(); O.Dispose(); P.Dispose();
            return s;
        }
        /// <summary>
        /// gets panel shape as Solid
        /// </summary>
        public Solid GetPanelSolid() {
            Surface s = GetPanelSurface();
            Solid S=null;
            if (Thickness[0]>0) S = s.Thicken(Thickness[0],false);
            else if (Thickness[1] > 0) S = s.Thicken(-Thickness[1], false);
            s.Dispose();
            return S;
        }
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
                double x = ArcPoints[i][ArcPoints[i].Length - 1].X / 2 + ArcPoints[j][0].X / 2;
                double y = ArcPoints[i][ArcPoints[i].Length - 1].Y / 2 + ArcPoints[j][0].Y / 2;
                double z = ArcPoints[i][ArcPoints[i].Length - 1].Z / 2 + ArcPoints[j][0].Z / 2;
                Point m = Point.ByCoordinates(x,y,z);
                Vector X = Face.VertexVectors[i][0].ToVector();
                Vector Y = Face.VertexVectors[i][0].Cross(Face.Normal).ToVector();
                if (Offset > 0) m = m.Add(Y.Normalized().Scale(-Offset));
                word w = word.ByStringOriginVectors(LabelPrefix + Face.E[i].Edge.Name, m, X, Y);
                labels.AddRange(w.Display(Scale));
                m.Dispose(); X.Dispose(); Y.Dispose(); w.Dispose();
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
                double x = ArcPoints[i][ArcPoints[i].Length - 1].X / 2 + ArcPoints[j][0].X / 2;
                double y = ArcPoints[i][ArcPoints[i].Length - 1].Y / 2 + ArcPoints[j][0].Y / 2;
                double z = ArcPoints[i][ArcPoints[i].Length - 1].Z / 2 + ArcPoints[j][0].Z / 2;
                Point m = Point.ByCoordinates(x, y, z);
                Vector X = Face.VertexVectors[i][0].Reverse().ToVector();
                Vector Y = Face.VertexVectors[i][0].Cross(Face.Normal).ToVector();
                if (Offset > 0) m = m.Add(Y.Normalized().Scale(-Offset));
                word w = word.ByStringOriginVectors(LabelPrefix + Face.E[i].Edge.Name, m, X, Y);
                labels.AddRange(w.Display(Scale));
                m.Dispose(); X.Dispose(); Y.Dispose(); w.Dispose();
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
            Point c = Face.Center.ToPoint();
            Vector X = Face.VertexVectors[0][0].Reverse().ToVector();
            Vector Y = Face.VertexVectors[0][0].Cross(Face.Normal).ToVector();
            word W = word.ByStringOriginVectors(LabelPrefix + Face.Name, c, X, Y);
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
            Point c = Face.Center.ToPoint();
            Vector X = Face.VertexVectors[0][0].ToVector();
            Vector Y = Face.VertexVectors[0][0].Cross(Face.Normal).ToVector();
            word W = word.ByStringOriginVectors(LabelPrefix + Face.Name, c, X, Y);
            labels.AddRange(W.Display(Scale));
            c.Dispose(); X.Dispose(); Y.Dispose(); W.Dispose();
            return labels.ToArray();
        }
    }

    /// <summary>
    /// panel object with holes that correspond to edge conector placement
    /// bevels are set to zero for hole addition
    /// </summary>
    public class panelHole : panel
    {
        /// <summary>
        /// circle list that contains panel hole geometry
        /// </summary>
        public List<circle> Holes { get; set; }

        //**METHODS**CONSTRUCTOR
        internal panelHole(face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle)
            : base(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, MinCornerRadius, MinFaceAngle, BevelAngle)
        {
            // initialize
            Holes = new List<circle>();
        }
        /// <summary>
        /// creates panel object with holes for connectors, defualt bevel angle is 0
        /// </summary>
        /// <param name="Face"></param>
        /// <param name="ThicknessFront"></param>
        /// <param name="ThicknessBack"></param>
        /// <param name="MinEdgeOffset"></param>
        /// <param name="MinCornerRadius"></param>
        /// <param name="MinFaceAngle"></param>
        /// <param name="BevelAngle"></param>
        /// <returns></returns>
        new public static panelHole ByFaceAndParameters(face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle = 0)
        { return new panelHole(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, MinCornerRadius, MinFaceAngle, BevelAngle); }

        //**MEDTHODS**ACTIONS
        /// <summary>
        /// get panel hole geometry data
        /// </summary>
        /// <returns></returns>
        public Circle[] GetPanelHoleProfile()
        {
            List<Circle> result = new List<Circle>();
            Holes.ForEach(c => result.Add(c.ToCircle()));
            return result.ToArray();
        }
        /// <summary>
        /// get panel hole geometry data as polycurves
        /// </summary>
        /// <returns></returns>
        public PolyCurve[] GetPanelHoleProfileAsPolyCurve()
        {
            Circle[] C = GetPanelHoleProfile();
            List<PolyCurve> P = new List<PolyCurve>(C.Length);
            for (int i = 0; i < C.Length; i++) P.Add(PolyCurve.ByJoinedCurves(new Curve[]{C[i]}));
            C.ForEach(c => c.Dispose());
            return P.ToArray();
        }
        /// <summary>
        /// returns solid panel with holes based on flat panel as single Solid
        /// </summary>
        /// <returns>Solid by Surface Thickening</returns>
        public Surface GetPanelHoleSurface()
        {
            Surface s = GetPanelSurface();
            List<Surface> S = new List<Surface>() { s };
            Circle[] C = GetPanelHoleProfile();
            for (int i = 0; i < C.Length; i++) S.Add(S[S.Count - 1].Trim(C[i], C[i].CenterPoint)[0] as Surface);
            Surface result = S[S.Count - 1];
            for (int i = 0; i < S.Count - 1; i++) if (!result.Equals(S[i])) S[i].Dispose();
            S = null;
            s.Dispose();
            C.ForEach(c => c.Dispose());
            return result;
        }
        /// <summary>
        /// get panel solid hole from mesh to panel front
        /// </summary>
        /// <returns></returns>
        public Solid GetPanelHoleSolidFront()
        { 
            Surface s = GetPanelHoleSurface();
            Solid S = s.Thicken(Thickness[0], false);
            s.Dispose();
            return S;
        }
        /// <summary>
        /// get panel solid hole
        /// </summary>
        public Solid GetPanelHoleSolid()
        {
            Surface s1 = GetPanelHoleSurface();
            Vector v = Face.Normal.ToVector();
            Surface s2 = (Surface) s1.Translate(v, -Thickness[1]);
            Solid S = s2.Thicken(Thickness[0]+Thickness[1], false);
            s1.Dispose(); v.Dispose(); s2.Dispose();
            return S;
        }

        /// <summary>
        /// add hole with center at closest point on panel to input point with given radius
        /// </summary>
        /// <param name="Point">Point is center of hole</param>
        /// <param name="Radius">Radius of hole</param>
        public void AddHole(point Point, double Radius) { Holes.Add(new circle(Point, Radius, Face.Normal));}

    }

    /// <summary>
    /// panel object with minimum edge offset determined by Z-axis elevation
    /// </summary>
    public class panelGradient : panel
    {
        internal panelGradient(face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle, double OffsetBaseZ)
            : base(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, MinCornerRadius, MinFaceAngle, BevelAngle)
        {
            for (int i = 0; i < Face.E.Count; i++)
            {
                double Offset = OffsetBaseZ * Face.E[i].Edge.MidPoint.Z;
                // check clearance in back
                if (EdgeOffset[i] < Offset) EdgeOffset[i] = Offset;
            }
            // corner arcs based on triangle vertex
            List<point[]> P = new List<point[]>(Face.E.Count);
            for (int i = 0; i < Face.E.Count; i++)
            {
                double sinA = Math.Sin(Face.Angles[i]);
                double sinB = Math.Sin(Face.Angles[i]/2);
                double tanB = Math.Abs(Math.Tan(Face.Angles[i]/2));
                double CornerRadius = MinCornerRadius * tanB;
                int j = (i + Face.E.Count - 1) % Face.E.Count;
                point[] arc = new point[] { Face.Vertices[i].Add(Face.VertexVectors[i][1].Scale(EdgeOffset[i] / sinA)).Add(Face.VertexVectors[i][0].Scale(EdgeOffset[j] / sinA)) };
                if (Face.Angles[i] < math.PI)
                {
                    point a0 = arc[0].Add(Face.VertexVectors[i][1].Scale(CornerRadius / tanB));
                    point a1 = arc[0].Add(Face.VertexVectors[i][2].Scale(CornerRadius / sinB - CornerRadius));
                    point a2 = arc[0].Add(Face.VertexVectors[i][0].Scale(CornerRadius / tanB));
                    arc = new point[] { a0, a1, a2 };
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
        public static panelGradient ByFaceAndParameters(face Face, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double MinCornerRadius, double MinFaceAngle, double BevelAngle = 0, double OffsetBaseZ = 0)
        { return new panelGradient(Face, ThicknessFront, ThicknessBack, MinEdgeOffset, MinCornerRadius, MinFaceAngle, BevelAngle, OffsetBaseZ); }

    }

    /// <summary>
    /// connectors at edges created for model
    /// </summary>
    public class connector
    {
        //**FIELDS
        internal int iAngle;
        internal double BevelAngle;
        internal Dictionary<char, double> Dim;
        internal double EdgeOffset;

        //**PROPERTIES
        /// <summary>
        /// gets reference edge in mesh
        /// </summary>
        public edge Edge { get; private set; }
        /// <summary>
        /// gets halfedges of faces being connected by connector
        /// </summary>
        public halfedge[] HalfEdges { get; private set; }
        /// <summary>
        /// get distance from mesh edge to panel edge
        /// </summary>
        public double Inset { get; set; }
        /// <summary>
        /// get angle between the faces being connected (radians)
        /// </summary>
        public double InsetAngle { get; private set; }
        /// <summary>
        /// get profile points as List of Vectors from given point on edge
        /// </summary>
        public List<vector> Profile { get; private set; }
        /// <summary>
        /// get Vector Array {X1,Y1,Z1,X2,Y2,Z2,N} 
        /// (halfedge1 CS),(halfedge1 CS),edge normal
        /// </summary>
        public vector[] Vectors { get; private set; }
        /// <summary>
        /// midpoint of edge
        /// </summary>
        public Point Midpoint { get { return Edge.MidPoint.ToPoint(); } }

        internal connector(halfedge e1, halfedge e2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        {
            // initialize properties
            this.BevelAngle = BevelAngle;
            Profile = new List<vector>();
            HalfEdges = new halfedge[] { e1, e2 };
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
            if (i1+i2 != 2)
            {
                // generate orthonormal coordinate system based on right-hand rule
                // X: normal vector of face (front face normal)
                // Y: orthonormal vector on face plane (towards center)
                // Z: tangent vector of edge (halfedge direction)
                vector X1 = e1.Face.Normal;
                vector Z1 = e1.Normalized();
                vector Y1 = X1.NormalizedCross(Z1);
                X1 = Z1.NormalizedCross(Y1);
                vector X2 = e2.Face.Normal;
                vector Z2 = e2.Normalized();
                vector Y2 = X2.NormalizedCross(Z2);
                X2 = Z2.NormalizedCross(Y2);
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
                vector eN = Edge.Normal[iAngle].Normalized();
                vector eY = eN.NormalizedCross(Z1);
                eN = Z1.NormalizedCross(eY);
                // case where halfedge count is three
                if (Edge.E.Count > 2)
                {
                    // set inset calculation based on smallest angle pertinent to edgeconnectors at edge
                    InsetAngle = Math.Min(Edge.Angle[0], Edge.Angle[1]);
                    // reverse normal
                    if (iAngle == 0) eN = eN.Reverse();
                }
                // calculate inset distance based on inset angle and panel offsets based on inset angle
                // check clearance in back
                EdgeOffset = 0;
                if (InsetAngle == 2 * math.PI || InsetAngle == math.PI) EdgeOffset = PanelMinOffset;
                if (InsetAngle < 2 * math.toRadians(90 - BevelAngle))
                    EdgeOffset = ThicknessBack / Math.Tan(InsetAngle/2) - ThicknessBack * Math.Tan(BevelAngle * math.PI / 180);
                // check clearance in front
                double OffsetAngle = PanelMinOffset / Math.Sin(InsetAngle/2) - ThicknessFront / Math.Tan(InsetAngle/2);
                if (EdgeOffset < OffsetAngle) EdgeOffset = OffsetAngle;
                Inset = Math.Max(Dim['H'] / 2 + EdgeOffset, (ThicknessBack + Dim['H']) / Math.Tan(InsetAngle/2));
                // initialize Vector Lists to store affine geometry information
                List<vector> P1 = new List<vector>();
                List<vector> P2 = new List<vector>();
                // v0: point closest to edge
                P1.Add(eN.Scale(-ThicknessBack / Math.Sin(Edge.Angle[iAngle]/2)));
                // v1: farthest point in first halfedge face direction (inset + 1.5 spacing)
                P1.Add(X1.Scale(-ThicknessBack).Add(Y1.Scale(Inset + Dim['S'] + Dim['H'] / 2)));
                // v2: point to create square edge with face and beginning of fillet arc
                P1.Add(P1[1].Add(X1.Scale(CornerRadius - Dim['H'])));
                // v3: center of fillet arc
                P1.Add(P1[2].Add(Y1.Scale(-CornerRadius)));
                // v4: end of fillet arc;
                P1.Add(P1[3].Add(X1.Scale(-CornerRadius)));
                // v5: point where two sides of the connector meet; line (v0,v5) is line of symmetry
                P1.Add(P1[0].Add(eN.Scale(-Dim['H'] / Math.Sin(Edge.Angle[iAngle]/2))));

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
                Vectors = new vector[] { X1, Y1, Z1, X2, Y2, Z2, eN };
            }
        }

        //**METHOD**CREATE
        /// <summary>
        /// creates edgeConnector object with given halfedge objects and parameters with spacing being the same as height
        /// </summary>
        public static connector ByHalfEdgeHeight(halfedge HalfEdge1, halfedge HalfEdge2, double Height, double Depth, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        { return new connector(HalfEdge1, HalfEdge2, Height, Depth, Height, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle); }
        /// <summary>
        /// creates edgeConnector object with given halfedge objects and parameters
        /// </summary>
        public static connector ByHalfEdge(halfedge HalfEdge1, halfedge HalfEdge2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        { return new connector(HalfEdge1, HalfEdge2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle); }
        /// <summary>
        /// creates edgeConnector object with given edge object and parameters
        /// </summary>
        public static connector ByEdge(edge Edge, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle)
        {
            if (Edge.E.Count < 2) return null;
            halfedge he = Edge.E[0];
            return new connector(Edge.E[0], Edge.E[1], Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle);
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
            Profile.ForEach(v => Points.Add(Edge.MidPoint.ToPoint().Add(v.ToVector())));
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
            Profile.ForEach(v => Points.Add(Point.Add(v.ToVector())));
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
        public Surface GetConnectorSurface(PolyCurve Profile)
        {
            if (Profile.Equals(null)|| !Profile.IsClosed) return null;
            if (!Profile.IsPlanar)
            {
                Vector N = HalfEdges[0].ToVector();
                Point p = Profile.EndPoint;
                Plane P = Plane.ByOriginNormal(p,N);
                Curve C = Profile.PullOntoPlane(P);
                Surface s = Surface.ByPatch(C);
                N.Dispose(); p.Dispose(); P.Dispose(); C.Dispose();
                return s;

            }
            return Surface.ByPatch(Profile);
        }
        /// <summary>
        /// returns solid panel based on panel profile as single Solid
        /// </summary>
        public Solid GetConnectorSolid(PolyCurve Profile)
        {
            if (Profile.Equals(null) || !Profile.IsClosed) return null;
            Point a = Profile.StartPoint.Add(Vectors[2].Scale(Dim['D'] / 2).ToVector());
            Point b = Profile.StartPoint.Add(Vectors[5].Scale(Dim['D'] / 2).ToVector());
            Line l = Line.ByStartPointEndPoint(a, b);
            Solid s = Profile.SweepAsSolid(l);
            a.Dispose();
            b.Dispose();
            l.Dispose();
            return s;
        }
        /// <summary>
        /// gets edge labels at edge on triangular mesh face backside as an Polycurve array
        /// </summary>
        public PolyCurve[] GetEdgeLabel(Point Point, double Scale=1.0/12, string LabelPrefix = "")
        {
            point p = point.ByCoordinates(Point.X, Point.Y, Point.Z);
            point pt = p.Add(Profile[5]);
            vector X = Vectors[3].Subtract(Vectors[0]);
            vector Y = Vectors[6].Reverse();
            if (Edge.Angle[iAngle] > math.PI - 0.001 && Edge.Angle[iAngle] < math.PI + 0.001) X = Profile[1].Subtract(Profile[0]);
            if (Edge.Angle[iAngle] > math.PI) { pt = p.Add(Profile[0]); X = X.Reverse(); Y = Y.Reverse(); }
            Surface s = GetConnectorSurface(GetConnectorProfileAtPoint(Point));
            Vector Z1 = s.NormalAtParameter(0.5,0.5);
            Vector Z2 = X.Cross(Y).Normalized().ToVector();
            if (!Z1.IsAlmostEqualTo(Z2)) X = X.Reverse();
            word w = word.ByStringOriginVectors(LabelPrefix + Edge.Name, pt.ToPoint(), X.ToVector(), Y.ToVector());
            List<PolyCurve> label = w.Display(Scale);
            Z1.Dispose(); Z2.Dispose();  w.Dispose();
            return label.ToArray();
        }


        //**MODS
        /// <summary>
        /// get connector profile bevel at point
        /// </summary>
        /// <param name="Point"></param>
        /// <returns></returns>
        public PolyCurve GetConnectorProfileBevel(Point Point)
        {
            point p = point.ByCoordinates(Point.X, Point.Y, Point.Z);
            // return default for angles smaller than 
            if (Edge.Angle[iAngle] <= 2 * math.toRadians(90 - BevelAngle)) return GetConnectorProfileAtPoint(Point);
            // initialize and generate Point List using EdgeConnector geometry data at given point
            List<point> Points = new List<point>(Profile.Count + 7);
            double ThicknessBack = Profile[0].Length * Math.Sin(Edge.Angle[iAngle] /2);
            double BevelInset = ThicknessBack * Math.Tan(math.toRadians(BevelAngle));
            //double DogBoneBevel = 2 * ToolRadius * Math.Sin(math.PI / 4 - BevelAngle * math.PI / 2);
            //double DogBoneChord = 2 * ToolRadius / Math.Sin(math.PI / 4 + BevelAngle * math.PI / 2);
            for (int i = 1; i < Profile.Count; i++) Points.Add(p.Add(Profile[i]));
            //Points.Add(Point.Add(Vectors[4].Scale(EdgeOffset + BevelInset + 2 * ToolRadius)).Add(Vectors[3].Scale(-ThicknessBack)));
            //Points.Add(Point.Add(Vectors[4].Scale(EdgeOffset + BevelInset + ToolRadius)).Add(Vectors[3].Scale(-ThicknessBack - ToolRadius)));
            Points.Add(p.Add(Vectors[4].Scale(EdgeOffset + BevelInset)).Add(Vectors[3].Scale(-ThicknessBack)));
            if (Edge.Angle[iAngle] == 2 * math.toRadians(180-BevelAngle)) { }
            else if (Edge.Angle[iAngle] > 2 * math.toRadians(180-BevelAngle))
            {
                //hyp = o/(sin(e/2)*tan(b) + cos(e/2))
                //Points.Add(p.Add(Vectors[6].Scale(EdgeOffset / (Math.Sin(Edge.Angle[iAngle] / 2) * Math.Tan(BevelAngle) + Math.Cos(Edge.Angle[iAngle] / 2)))));
                Points.Add(p.Add(Vectors[6].Scale(EdgeOffset * Math.Sin(math.toRadians(90-BevelAngle)) / Math.Sin(Edge.Angle[iAngle] / 2 + math.toRadians(BevelAngle) - math.PI/2))));
            }
            else 
                //if (Edge.Angle[iAngle] < 2 * math.toRadians(180 - BevelAngle))
            {
                Points.Add(p.Add(Vectors[4].Scale(EdgeOffset + BevelInset / 2)).Add(Vectors[3].Scale(-ThicknessBack / 2)));
                Points.Add(p.Add(Vectors[1].Scale(EdgeOffset + BevelInset / 2)).Add(Vectors[0].Scale(-ThicknessBack / 2)));
            }
            Points.Add(p.Add(Vectors[1].Scale(EdgeOffset + BevelInset)).Add(Vectors[0].Scale(-ThicknessBack)));
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
            List<Point> Pts = new List<Point>(Points.Count);
            Points.ForEach(pt => Pts.Add(pt.ToPoint()));
            List<Curve> Curves = new List<Curve> {
                Line.ByStartPointEndPoint(Pts[0], Pts[1]),
                Arc.ByCenterPointStartPointEndPoint(Pts[2],Pts[1],Pts[3]),
                Line.ByStartPointEndPoint(Pts[3], Pts[4]),
                Line.ByStartPointEndPoint(Pts[4], Pts[5]),
                Arc.ByCenterPointStartPointEndPoint(Pts[6],Pts[5],Pts[7])};
            for (int j = 7; j < Points.Count; j++) Curves.Add(Line.ByStartPointEndPoint(Pts[j], Pts[(j + 1) % Pts.Count]));
            // create closed polycurve from curves
            PolyCurve Result = PolyCurve.ByJoinedCurves(Curves);
            // dispose unmanaged resources
            Pts.ForEach(pt => pt.Dispose()); Curves.ForEach(c => c.Dispose());
            return Result;
        }
        /// <summary>
        /// get bevel connector as points
        /// </summary>
        /// <param name="Point"></param>
        /// <returns></returns>
        public List<Point> GetConnectorProfileBevelPoint(Point Point)
        {
            point p = point.ByCoordinates(Point.X, Point.Y, Point.Z);
            if (Edge.Angle[iAngle] <= 2 * math.toRadians(90 - BevelAngle)) return null;
            List<point> Points = new List<point>(Profile.Count + 7);
            double ThicknessBack = Profile[0].Length * Math.Sin(Edge.Angle[iAngle] / 2);
            double BevelInset = ThicknessBack * Math.Tan(math.toRadians(BevelAngle));
            for (int i = 1; i < Profile.Count; i++) Points.Add(p.Add(Profile[i]));
            Points.Add(p.Add(Vectors[4].Scale(EdgeOffset + BevelInset)).Add(Vectors[3].Scale(-ThicknessBack)));
            if (Edge.Angle[iAngle] == 2 * math.toRadians(180 - BevelAngle)) { }
            else if (Edge.Angle[iAngle] < 2 * math.toRadians(180 - BevelAngle))
            {
                Points.Add(p.Add(Vectors[4].Scale(EdgeOffset + BevelInset / 2)).Add(Vectors[3].Scale(-ThicknessBack / 2)));
                Points.Add(p.Add(Vectors[1].Scale(EdgeOffset + BevelInset / 2)).Add(Vectors[0].Scale(-ThicknessBack / 2)));
            }
            else if (Edge.Angle[iAngle] > 2 * math.toRadians(180 - BevelAngle))
            {
                Points.Add(p.Add(Vectors[6].Scale(EdgeOffset / (Math.Sin(Edge.Angle[iAngle] / 2) * Math.Tan(BevelAngle) + Math.Cos(Edge.Angle[iAngle] / 2)))));
            }
            Points.Add(p.Add(Vectors[1].Scale(EdgeOffset + BevelInset)).Add(Vectors[0].Scale(-ThicknessBack)));
            List<Point> Pts = new List<Point>(Points.Count);
            Points.ForEach(pt => Pts.Add(pt.ToPoint()));
            Debug.WriteLine("Points:");
            for (int i = 0; i < Points.Count; i++)
            {
                Debug.WriteLine(Points[i]);
            }
            return Pts;
        }
    }

    /// <summary>
    /// edgeConnector object for panel system with edge offset determined by height in Z-axis
    /// </summary>
    public class connectorGradient : connector
    {
        internal connectorGradient(halfedge e1, halfedge e2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double Z)
            : base(e1, e2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle)
        {
            int i1 = Edge.E.IndexOf(e1);
            int i2 = Edge.E.IndexOf(e2);
            // do not create connector in the case where HalfEdge index is (0,2)
            if (i1 + i2 != 2)
            {
                if (EdgeOffset < Edge.MidPoint.Z * Z) EdgeOffset = Edge.MidPoint.Z * Z;
                Inset = Math.Max(Dim['H'] / 2 + EdgeOffset, (ThicknessBack + Dim['H']) / Math.Tan(InsetAngle / 2));
                // initialize Vector Lists to store affine geometry information
                List<vector> P1 = new List<vector>();
                List<vector> P2 = new List<vector>();
                // v0: point closest to edge
                P1.Add(Vectors[6].Scale(-ThicknessBack / Math.Sin(Edge.Angle[iAngle] / 2)));
                // v1: farthest point in first halfedge face direction (inset + 1.5 spacing)
                P1.Add(Vectors[0].Scale(-ThicknessBack).Add(Vectors[1].Scale(Inset + Dim['S'] + Dim['H'] / 2)));
                // v2: point to create square edge with face and beginning of fillet arc
                P1.Add(P1[1].Add(Vectors[0].Scale(CornerRadius - Dim['H'])));
                // v3: center of fillet arc
                P1.Add(P1[2].Add(Vectors[1].Scale(-CornerRadius)));
                // v4: end of fillet arc;
                P1.Add(P1[3].Add(Vectors[0].Scale(-CornerRadius)));
                // v5: point where two sides of the connector meet; line (v0,v5) is line of symmetry
                P1.Add(P1[0].Add(Vectors[6].Scale(-Dim['H'] / Math.Sin(Edge.Angle[iAngle] / 2))));

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
        /// creates edgeConnector object with given halfedge objects and parameters
        /// </summary>
        public static connectorGradient ByHalfEdge(halfedge HalfEdge1, halfedge HalfEdge2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double OffsetBaseZ)
        { return new connectorGradient(HalfEdge1, HalfEdge2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, OffsetBaseZ); }
        /// <summary>
        /// creates edgeConnector object with given edge object and parameters
        /// </summary>
        public static connectorGradient ByEdge(edge e, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double OffsetBaseZ)
        {
            if (e.E.Count < 2) return null;
            return new connectorGradient(e.E[0], e.E[1], Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, OffsetBaseZ);
        }
    }
    /// <summary>
    /// edgeConnector object for panel system using dowel nuts for assembly
    /// </summary>
    public class connectorHole : connector
    {
        //**FIELDS
        /// <summary>
        /// dowel nut radius
        /// </summary>
        double Radius;
        /// <summary>
        /// circle list of dowel nut geometry
        /// </summary>
        List<circle> Pockets;

        internal connectorHole(halfedge e1, halfedge e2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double PocketRadius)
            : base(e1, e2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle)
        {
            Radius = PocketRadius;
            Pockets = new List<circle>();
        }

        //**METHOD**CREATE
        /// <summary>
        /// creates edgeConnector object with given halfedge objects and parameters
        /// </summary>
        public static connectorHole ByHalfEdge(halfedge HalfEdge1, halfedge HalfEdge2, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double PocketRadius)
        { return new connectorHole(HalfEdge1, HalfEdge2, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, PocketRadius); }
        /// <summary>
        /// creates edgeConnector object with given edge object and parameters
        /// </summary>
        public static connectorHole[] ByEdge(edge e, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double CornerRadius, double BevelAngle, double PocketRadius)
        {
            if (e.E.Count < 2) return null;
            if (e.E.Count == 2) return new connectorHole[] { connectorHole.ByHalfEdge(e.E[0], e.E[1], Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, PocketRadius) };
            return new connectorHole[]{ new connectorHole(e.E[0], e.E[1], Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, CornerRadius, BevelAngle, PocketRadius)};
        }
        /// <summary>
        /// returns connector surface with holes at midpoint
        /// </summary>
        /// <returns></returns>
        public Surface GetConnectorHoleSurfaces()
        {
            Surface s = GetConnectorSurface(GetConnectorProfileAtMidpoint());
            List<Surface> S = new List<Surface>() { s };
            Circle[] C = GetPockets();
            for (int i = 0; i < C.Length; i++) S.Add(S[S.Count - 1].Trim(C[i], C[i].CenterPoint)[0] as Surface);
            Surface result = S[S.Count - 1];
            for (int i = 0; i < S.Count - 1; i++) if (!result.Equals(S[i])) S[i].Dispose();
            S = null;
            s.Dispose();
            C.ForEach(c => c.Dispose());
            return result;
        }
        /// <summary>
        /// gets solid panel with holes based on flat panel as single Solid
        /// </summary>
        public Solid GetConnectorHoleSolids()
        {
            Surface s = GetConnectorHoleSurfaces();
            Solid S = s.Thicken(Dim['D'], true);
            s.Dispose();
            return S;
        }
        /// <summary>
        /// adds pockets based on edge condition and inset with given radius
        /// </summary>
        /// <param name="Point">Point is origin of connector</param>
        /// <param name="Radius">Radius of pocket</param>
        public void AddPocket(point Point, double Radius)
        {
            int i = Edge.E.IndexOf(HalfEdges[0]) + Edge.E.IndexOf(HalfEdges[1]);
            int j = 0;
            if (i > 2) j = 1;
            double a = Edge.Angle[j];
            point p1a = Point.Add(Profile[3]);
            point p1b = p1a.Subtract(Vectors[1].Normalized().Scale(Dim['D']));
            point p2a = Point.Add(Profile[Profile.Count - 3]);
            point p2b = p2a.Subtract(Vectors[4].Normalized().Scale(Dim['D']));
            vector N = Vectors[5];
            if (a == InsetAngle)
            {
                Pockets.Add(circle.ByCenterRadiusNormal(p1a, Radius, N));
                Pockets.Add(circle.ByCenterRadiusNormal(p1b, Radius, N));
                Pockets.Add(circle.ByCenterRadiusNormal(p2a, Radius, N));
                Pockets.Add(circle.ByCenterRadiusNormal(p2b, Radius, N));
            }
            else if (j == 1)
            {
                Pockets.Add(circle.ByCenterRadiusNormal(p2a, Radius, N));
                Pockets.Add(circle.ByCenterRadiusNormal(p2b, Radius, N));
            }
            else if (j == 0)
            {
                Pockets.Add(circle.ByCenterRadiusNormal(p1a, Radius, N));
                Pockets.Add(circle.ByCenterRadiusNormal(p1b, Radius, N));
            }
        }
        /// <summary>
        /// gets pockets as Dynamo Circle object array
        /// </summary>
        public Circle[] GetPockets()
        {
            List<Circle> result = new List<Circle>(Pockets.Count);
            Pockets.ForEach(c => result.Add(c.ToCircle()));
            return result.ToArray();
        }
    }
    /// <summary>
    /// panel system
    /// </summary>
    public class panelSystem
    {
        //**FIELDS
        internal Dictionary<face, panel> F;
        internal Dictionary<edge, connector[]> E;
        internal mesh M;

        //**PROPERTIES**QUERY
        /// <summary>
        /// panel array
        /// </summary>
        public panel[] Panels { get { return F.Values.ToArray(); } }
        /// <summary>
        /// edgeconnector array
        /// </summary>
        public connector[] Connectors
        {
            get
            {
                List<connector> result = new List<connector>();
                for (int i = 0; i < E.Count; i++)
                {
                    result.AddRange(E.Values.ToArray()[i]);
                }
                return result.ToArray();
            }
        }

        //**CONSTRUCTOR
        internal panelSystem(mesh Mesh) { M = Mesh; }
        internal panelSystem(mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double ConnectorFilletRadius, double BevelAngle)
            : this(Mesh)
        { GenerateEdgeConnectors(Height, Depth, Spacing, ThicknessFront, ThicknessBack, MinEdgeOffset, ConnectorFilletRadius, BevelAngle); GeneratePanels(ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, BevelAngle); }

        //**METHODS**CREATE
        /// <summary>
        /// returns default empty panel system with associated mesh
        /// </summary>
        /// <param name="Mesh">Mesh Object</param>
        /// <returns>Panel System</returns>
        public static panelSystem ByMesh(mesh Mesh) { return new panelSystem(Mesh); }
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
        public static panelSystem ByMeshParameters(mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double ConnectorFilletRadius, double BevelAngle)
        {
            panelSystem P = new panelSystem(Mesh);
            P.GenerateEdgeConnectors(Height, Depth, Spacing, ThicknessFront, ThicknessBack, MinEdgeOffset, ConnectorFilletRadius, BevelAngle);
            P.GeneratePanels(ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, BevelAngle);
            return P;
        }
        /// <summary>
        /// creates panel system with edge offset based on elevation
        /// </summary>
        public static panelSystem ByMeshParametersGradient(mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double ConnectorFilletRadius, double BevelAngle, double OffsetBaseZ)
        {
            panelSystem P = new panelSystem(Mesh);
            P.GenerateEdgeConnectorsGradient(Height, Depth, Spacing, ThicknessFront, ThicknessBack, MinEdgeOffset, ConnectorFilletRadius, BevelAngle, OffsetBaseZ);
            P.GeneratePanelsGradient(ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, BevelAngle, OffsetBaseZ);
            return P;
        }
        /// <summary>
        /// creates panel system with holes for bolt to dowel nut connection
        /// </summary>
        public static panelSystem ByMeshParametersHole(mesh Mesh, double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double ConnectorFilletRadius, double BevelAngle, double PanelHoleRadius, double ConnectorPocketRadius, int NumConnectors)
        {
            panelSystem P = new panelSystem(Mesh);
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
            E = new Dictionary<edge, connector[]>(M.E2.Count + M.E3.Count);
            // generate EdgeConnector
            M.E2.ForEach(e => E.Add(e, new connector[] { connector.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle) }));
            M.E3.ForEach(e => E.Add(e, new connector[] { connector.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle) }));
        }
        /// <summary>
        /// generate edgeconnector geometry based on system inputs
        /// </summary>
        internal void GenerateEdgeConnectorsGradient(double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double ConnectorFilletRadius, double BevelAngle, double OffsetBaseZ)
        {
            // initialize mesh edge to edge connector dictionary
            E = new Dictionary<edge, connector[]>(M.E2.Count + M.E3.Count);
            // generate EdgeConnector
            M.E2.ForEach(e => E.Add(e, new connector[] { connectorGradient.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle, OffsetBaseZ) }));
            M.E3.ForEach(e => E.Add(e, new connector[] { connectorGradient.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle, OffsetBaseZ) }));
        }
        /// <summary>
        /// generates new Edge/EdgeConnector Dictionary with given parameters
        /// </summary>
        internal void GenerateEdgeConnectorsHole(double Height, double Depth, double Spacing, double ThicknessFront, double ThicknessBack, double PanelMinOffset, double ConnectorFilletRadius, double BevelAngle, double ConnectorPocketRadius)
        {
            // initialize mesh edge to edge connector dictionary
            E = new Dictionary<edge, connector[]>(M.E2.Count + M.E3.Count);
            // generate EdgeConnector for edges with two halfedges
            M.E2.ForEach(e => E.Add(e, connectorHole.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle, ConnectorPocketRadius)));
            // generate EdgeConnector for edges with three halfedges
            M.E3.ForEach(e => E.Add(e, connectorHole.ByEdge(e, Height, Depth, Spacing, ThicknessFront, ThicknessBack, PanelMinOffset, ConnectorFilletRadius, BevelAngle, ConnectorPocketRadius)));
            Connectors.ForEach(c => ((connectorHole) c).AddPocket(c.Edge.MidPoint, ConnectorPocketRadius));
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
            F = new Dictionary<face, panel>(M.Faces.Count);
            // iterate through mesh faces to generate panels
            for (int i = 0; i < M.Faces.Count; i++)
            {
                // generate triangle panel with input values at given face
                panel panel = panel.ByFaceAndParameters(M.Faces[i], ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, M.MinFaceAngle, BevelAngle);
                // add face/panel pair to dictionary
                F.Add(M.Faces[i], panel);
                // iterate through edges and add holes based on edge connector calculations
                for (int j = 0; j < M.Faces[i].Edges.Length; j++)
                    if (E.ContainsKey(M.Faces[i].Edges[j]))
                    {
                        // midpoint of edge
                        point p = M.Faces[i].Edges[j].MidPoint;
                        // orthogonal vector to edge that is on the plane of triangle face with direction towards center of face
                        vector Y = M.Faces[i].Normal.NormalizedCross(M.Faces[i].E[j]);
                    }
            }
            return true;
        }
        /// <summary>
        /// generate panels based on system inputs
        /// </summary>
        internal bool GeneratePanelsGradient(double ThicknessFront, double ThicknessBack, double MinEdgeOffset, double CornerOffset, double BevelAngle, double OffsetBaseZ)
        {
            if (!(E.Count > 0)) return false;
            // initialize mesh face to panel dictionary
            F = new Dictionary<face, panel>(M.Faces.Count);
            // iterate through mesh faces to generate panels
            for (int i = 0; i < M.Faces.Count; i++)
            {
                // generate triangle panel with input values at given face
                panelGradient panel = panelGradient.ByFaceAndParameters(M.Faces[i], ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, M.MinFaceAngle, BevelAngle, OffsetBaseZ);
                // add face/panel pair to dictionary
                F.Add(M.Faces[i], panel);
                // iterate through edges and add holes based on edge connector calculations
                for (int j = 0; j < M.Faces[i].Edges.Length; j++)
                    if (E.ContainsKey(M.Faces[i].Edges[j]))
                    {
                        // midpoint of edge
                        point p = M.Faces[i].Edges[j].MidPoint;
                        // orthogonal vector to edge that is on the plane of triangle face with direction towards center of face
                        vector Y = M.Faces[i].Normal.NormalizedCross(M.Faces[i].E[j]);
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
            F = new Dictionary<face, panel>(M.Faces.Count);
            // iterate through mesh faces to generate panels
            for (int i = 0; i < M.Faces.Count; i++)
            {
                // generate triangle panel with input values at given face
                panelHole panel = panelHole.ByFaceAndParameters(M.Faces[i], ThicknessFront, ThicknessBack, MinEdgeOffset, CornerOffset, M.MinFaceAngle, BevelAngle);
                // add face/panel pair to dictionary
                F.Add(M.Faces[i], panel);
                // iterate through edges and add holes based on edge connector calculations
                for (int j = 0; j < M.Faces[i].Edges.Length; j++)
                    if (E.ContainsKey(M.Faces[i].Edges[j]))
                    {
                        // orthogonal vector to edge that is on the plane of triangle face with direction towards center of face
                        vector Y = M.Faces[i].Normal.NormalizedCross(M.Faces[i].E[j]);
                        if (NumConnectors.Equals(1))
                        {
                            // calculate first hole position based on edge using midpoint as origin
                            // distance between hole center and edge is equal to inset property of edge
                            point p = M.Faces[i].Edges[j].MidPoint.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Inset));
                            // add hole to list of holes of TrianglePanel
                            panel.AddHole(p, HoleRadius);
                            // calculate second hole position based on edge using midpoint as origin
                            // distance between hole center and edge is equal to the sum of the inset property of edge and the width property of edge
                            vector V = Y.Scale(E[M.Faces[i].Edges[j]][0].Dim['S']);
                            p = p.Add(V);
                            // add hole to list of holes of TrianglePanel
                            panel.AddHole(p, HoleRadius);
                        }
                        else
                        {
                            for (int k = 0; k < NumConnectors; k++)
                            {
                                point p = M.Faces[i].Edges[j].HalfEdges[0].V[0];
                                vector v = M.Faces[i].Edges[j].HalfEdges[0];
                                vector v1 = v.Scale(v.Length * (2 * k + 1) / (2 * NumConnectors));
                                vector v2 = Y.Scale(E[M.Faces[i].Edges[j]][0].Inset);
                                p = p.Add(v1).Add(v2);
                                // add hole to list of holes of TrianglePanel
                                panel.AddHole(p, HoleRadius);
                                // calculate second hole position based on edge using midpoint as origin
                                // distance between hole center and edge is equal to the sum of the inset property of edge and the width property of edge
                                p = p.Add(Y.Scale(E[M.Faces[i].Edges[j]][0].Dim['S']));
                                // add hole to list of holes of TrianglePanel
                                panel.AddHole(p, HoleRadius);
                            }
                        }
                    }
            }
            return true;
        }
    }
}
