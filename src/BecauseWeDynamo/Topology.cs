using System;
using System.Collections.Generic;
using Autodesk.DesignScript.Runtime;
using Autodesk.DesignScript.Geometry;
using System.Linq;
using Geometry;

namespace Topology
{
    /// <summary>
    /// HalfEdge: ordered vertex array with reference to face, edge, normal, angle, and length 
    /// </summary>
    public class halfedge: vector
    {
        //**FIELD
        internal vertex[] V;

        //**PROPERTIES** //**QUERY**
        /// <summary>
        /// gets angle for halfedge at edge
        /// </summary>
        public double Angle { get; set; }
        /// <summary>
        /// get normal vector of halfedge face
        /// </summary>
        public vector Normal { get { return Face.Normal; } }
        /// <summary>
        /// gets Edge that contains this halfedge
        /// </summary>
        public edge Edge { get; private set; }
        /// <summary>
        /// get Face that contins this halfedge
        /// </summary>
        public face Face { get; private set; }

        //**CONSTRUCTOR**
        internal halfedge(vertex A, vertex B) : base(B.X - A.X, B.Y - A.Y, B.Z - A.Z)
        {
            Angle = 2*math.PI;
            V = new vertex[] { A, B };

        }
        internal halfedge(vertex A, vertex B, edge Edge, face Face)
            : this(A, B)
        { this.Edge = Edge; this.Face = Face; }
        internal halfedge(IEnumerable<vertex> Vertices) : this(Vertices.ElementAt(0), Vertices.ElementAt(1)) { }
        internal halfedge(IEnumerable<vertex> Vertices, edge Edge, face Face)
            : this(Vertices)
        { this.Edge = Edge; this.Face = Face; }

        //**METHODS** //**CREATE**
        /// <summary>
        /// creates HalfEdge instance
        /// </summary>
        /// <param name="Vertices">Vertices</param>
        /// <returns>HalfEdge</returns>
        public static halfedge ByVertices(IEnumerable<vertex> Vertices) { return new halfedge(Vertices); }
        /// <summary>
        /// adds reference edge if halfedge is part of edge
        /// and adds edge to vertices
        /// </summary>
        /// <param name="Edge">Mesh Edge</param>
        /// <returns>true if succeeded, false if failed</returns>
        public bool AddEdge(edge Edge)
        {
            if (Edge.E.Contains(this))
            {
                this.Edge = Edge;
                V[0].AddEdge(Edge);
                V[1].AddEdge(Edge);
                return true;
            }
            return false;
        }
        /// <summary>
        /// adds reference face if halfedge is part of face
        /// and adds face to vertices
        /// </summary>
        /// <param name="Face">Mesh Face</param>
        /// <returns>>true if succeeded, false if failed</returns>
        public bool AddFace(face Face)
        {
            if (Face.E.Contains(this))
            {
                this.Face = Face;
                V[0].AddFace(Face);
                V[1].AddFace(Face);
                return true;
            }
            return false;
        }
        /// <summary>
        /// flips direction of halfedge
        /// ie used when fliping face normals
        /// </summary>
        public void FlipDirection()
        {
            List<vertex> temp = new List<vertex>(V);
            V[0] = temp[1];
            V[1] = temp[0];
            temp = null;
            this.X = -X;
            this.Y = -Y;
            this.Z = -Z;
        }
    }

    /// <summary>
    /// Vertex: coordinate list with reference to face and edge sets
    /// </summary>
    public class vertex : point, IEquatable<vertex>
    {
        //**PROPERTIES** //**QUERY**
        /// <summary>
        /// Vertex Geometry as Point Object
        /// </summary>
        public Point Point { get { return Point.ByCoordinates(X, Y, Z); } }
        /// <summary>
        /// Edges connected to Vertex
        /// </summary>
        public HashSet<edge> Edges { get; private set; }
        /// <summary>
        /// Faces connected to Vertex
        /// </summary>
        public HashSet<face> Faces { get; private set; }

        //**CONSTRUCTOR**
        internal vertex(point Point) : base(Point.X, Point.Y, Point.Z)
        {
            Edges = new HashSet<edge>();
            Faces = new HashSet<face>();
        }
        internal vertex(point Point, IEnumerable<edge> Edges, IEnumerable<face> Faces)
            : base(Point.X, Point.Y, Point.Z)
        {
            this.Edges = new HashSet<edge>(Edges);
            this.Faces = new HashSet<face>(Faces);
        }
        internal vertex(double X, double Y, double Z): base(X,Y,Z)
        {
            Edges = new HashSet<edge>();
            Faces = new HashSet<face>();
        }
        internal vertex(double X, double Y, double Z, IEnumerable<edge> Edges, IEnumerable<face> Faces): base(X,Y,Z)
        {
            this.Edges = new HashSet<edge>(Edges);
            this.Faces = new HashSet<face>(Faces);
        }

        //**METHODS** //**CREAT**
        /// <summary>
        /// creates empty vertex at point
        /// </summary>
        /// <param name="Point">Point</param>
        /// <returns>Vertex</returns>
        public static vertex ByPoint(point Point) { return new vertex(Point); }

        //**METHODS** //**ACTION**
        /// <summary>
        /// return XYZ coordinates of Vertex
        /// </summary>
        /// <returns>XYZ Coordinates</returns>
        [MultiReturn(new[] { "X", "Y", "Z" })]
        public Dictionary<string, double> GetCoordinates() { return new Dictionary<string, double> { { "X", X }, { "Y", Y }, { "Z", Z } }; }
        /// <summary>
        /// adds edge to vertex if connected
        /// </summary>
        /// <param name="Edge">Edge to be added</param>
        public void AddEdge(edge Edge) { if (Edge.Vertices.Contains(this) && !Edges.Contains(Edge)) Edges.Add(Edge); }
        /// <summary>
        /// adds edges to vertex if connected
        /// </summary>
        /// <param name="Edges">Edges</param>
        public void AddEdges(IEnumerable<edge> Edges) { for (int i = 0; i < Edges.Count(); i++) AddEdge(Edges.ElementAt(i)); }
        /// <summary>
        /// adds face to vertex if connected
        /// </summary>
        /// <param name="Face">Face</param>
        public void AddFace(face Face) { if (Face.Vertices.Contains(this) && !Faces.Contains(Face)) Faces.Add(Face); }
        /// <summary>
        /// adds faces to vertex if connected
        /// </summary>
        /// <param name="Faces">Faces</param>
        public void AddFaces(IEnumerable<face> Faces) { for (int i = 0; i < Faces.Count(); i++) AddFace(Faces.ElementAt(i)); }
        /// <summary>
        /// returns distance to given vertex
        /// </summary>
        /// <param name="Vertex">Vertex</param>
        /// <returns>Distance</returns>
        public double DistanceTo(vertex Vertex)
        {
            double x = X - Vertex.X;
            double y = Y - Vertex.Y;
            double z = Z - Vertex.Z;
            return Math.Sqrt(x * x + y * y + z * z);
        }
        /// <summary>
        /// returns distance to given point
        /// </summary>
        /// <param name="Point">Point</param>
        /// <returns>Distance</returns>
        new public double DistanceTo(point Point)
        {
            double x = X - Point.X;
            double y = Y - Point.Y;
            double z = Z - Point.Z;
            return Math.Sqrt(x * x + y * y + z * z);
        }
        /// <summary>
        /// checks if vertex is located at point
        /// </summary>
        public bool IsAtPoint(point Point) { return (X == Point.X && Y == Point.Y && Z == Point.Z); }
        /// <summary>
        /// checks if vertex is located at Dynamo Point
        /// </summary>
        public bool IsAtPoint(Point Point) { return (X == Point.X && Y == Point.Y && Z == Point.Z); }
        /// <summary>
        /// checks if vertex is located at coordinates
        /// </summary>
        public bool IsAtPoint(double X, double Y, double Z) { return (X == this.X && Y == this.Y && Z == this.Z); }
        //**METHODS**IEQUATABLE
        /// <summary>
        /// gets object equality
        /// </summary>
        public override bool Equals(Object Object) { return this.Equals(Object as vertex); }
        /// <summary>
        /// gets vertex equality
        /// </summary>
        public bool Equals(vertex Vertex)
        {
            if (Object.ReferenceEquals(Vertex, null)) return false;
            if (Object.ReferenceEquals(this, Vertex)) return true;
            if (this.GetType() != Vertex.GetType()) return false;
            return (X == Vertex.X && Y == Vertex.Y && Z == Vertex.Z && Vertex.Edges == Edges && Vertex.Faces == Faces);
        }
        /// <summary>
        /// gets hashcode
        /// </summary>
        public override int GetHashCode() { return string.Format("{0}-{1}-{2}-{3}-{4}", X, Y, Z, Edges.GetHashCode(), Faces.GetHashCode()).GetHashCode(); }
    }

    /// <summary>
    /// Edge: HalfEdge List with reference to faces, vertices, length, name, angle, normal, midpoint, and length
    /// </summary>
    public class edge
    {
        //**FIELDS**
        internal List<halfedge> E;
        internal double[] N;

        //**PROPERTIES** //**QUERY**
        /// <summary>
        /// Edge Length
        /// </summary>
        public double Length { get { return E.ElementAt(0).Length; } }
        /// <summary>
        /// Edge Name for Label
        /// </summary>
        public string Name { get; set; }
        /// <summary>
        /// Angle between Faces at Edge (radians)
        /// </summary>
        public double[] Angle { get; set; }
        /// <summary>
        /// HalfEdge Array
        /// </summary>
        public halfedge[] HalfEdges { get { return E.ToArray(); } }
        /// <summary>
        /// Normal Vector
        /// </summary>
        public vector[] Normal
        {
            get
            {
                if (!(N.Length > 2)) return null;
                List<vector> V = new List<vector>();
                for (int i = 0; i < N.Length / 3; i++) V.Add(vector.ByCoordinates(N[i * 3], N[i * 3 + 1], N[i * 3 + 2]));
                return V.ToArray();
            }
        }
        /// <summary>
        /// Midpoint
        /// </summary>
        public point MidPoint
        {
            get
            {
                vertex[] V = Vertices; if (V.Length < 2 || V.Equals(null)) return null;
                return point.ByCoordinates(V[0].X / 2 + V[1].X / 2, V[0].Y / 2 + V[1].Y / 2, V[0].Z / 2 + V[1].Z / 2);
            }
        }
        /// <summary>
        /// Face List
        /// </summary>
        public List<face> Faces { get { List<face> F = new List<face>(E.Count); E.ToList().ForEach(e => F.Add(e.Face)); return F; } }
        /// <summary>
        /// Vertex Array
        /// </summary>
        public vertex[] Vertices { get { if (E.Count > 0) return E.ElementAt(0).V; return null; } }

        //**CONSTRUCTOR**
        internal edge() { E = new List<halfedge>(); Name = ""; Angle = new double[] { 2*Math.PI }; N = null; }
        internal edge(IEnumerable<halfedge> HalfEdges) : this() { E = new List<halfedge>(HalfEdges); Vertices.ForEach(v => v.AddEdge(this)); }
        internal edge(IEnumerable<halfedge> HalfEdges, string Name) : this(HalfEdges) { this.Name = Name; }

        //**METHODS**CREATE
        /// <summary>
        /// create and Edge object as an array of HalfEdges
        /// </summary>
        /// <param name="HalfEdges">HalfEdges</param>
        /// <returns>Edge</returns>
        public static edge ByHalfEdges(IEnumerable<halfedge> HalfEdges) { return new edge(HalfEdges); }

        //**METHODS** //**ACTION**
        /// <summary>
        /// get other vertex
        /// </summary>
        /// <param name="Vertex">Vertex</param>
        /// <returns>Other Vertex</returns>
        public vertex GetOtherVertex(vertex Vertex)
        {
            if (Vertices[0].Equals(Vertex)) return Vertices[1];
            if (Vertices[1].Equals(Vertex)) return Vertices[0];
            return null;
        }
        /// <summary>
        /// get normal vector to two faces that create an angle
        /// </summary>
        /// <param name="eA">HalfEdge A</param>
        /// <param name="eB">HalfEdge B</param>
        /// <returns></returns>
        internal double[] GetAngleNormal(halfedge eA, halfedge eB)
        {
            if (!E.Contains(eA) || !E.Contains(eB)) return null;
            vector aN = eA.Face.Normal; vector aX = eA; vector aY = aN.NormalizedCross(aX);
            vector bN = eB.Face.Normal; vector bX = eB; vector bY = bN.NormalizedCross(bX);
            vector eN = aN.Add(bN).Normalized(); vector eY = eN.Reverse();
            point M = MidPoint; point Ma = M.Add(aY); point Mb = M.Add(bY); point Me = M.Add(eY);
            arc arc = arc.ByThreePoints(Ma, Me, Mb);
            double[] result = new double[] { arc.SweepAngle, eN.X, eN.Y, eN.Z };
            return result;
        }
        /// <summary>
        /// creates line geometry based on edge
        /// </summary>
        /// <returns>Line</returns>
        public Line GetLine()
        {
            Point a = Vertices[0].Point;
            Point b = Vertices[1].Point;
            Line output = Line.ByStartPointEndPoint(a, b);
            a.Dispose(); b.Dispose();
            return output;
        }
        /// <summary>
        /// checks to see if edge has same geometric properties as given line
        /// </summary>
        /// <param name="Line">Line</param>
        /// <returns>Boolean</returns>
        public bool IsAtCurve(Curve Line)
        {
            return (Vertices[0].IsAtPoint(Line.EndPoint) && Vertices[1].IsAtPoint(Line.StartPoint)) || (Vertices[0].IsAtPoint(Line.StartPoint) && Vertices[1].IsAtPoint(Line.EndPoint));
        }
        /// <summary>
        /// rename edge
        /// </summary>
        /// <param name="Name"></param>
        public void Rename(string Name)
        {
            this.Name = Name;
        }
    }

    /// <summary>
    /// Spline: HalfEdge List with reference to vertices
    /// </summary>
    public class spline
    {
        /// <summary>
        /// Hlaf Edge List
        /// </summary>
        public List<edge> Edges { get; set; }
        /// <summary>
        /// Vertex List
        /// </summary>
        public List<vertex> Vertices { get; set; }

        internal spline() { Edges = new List<edge>(); Vertices = new List<vertex>(); }
        internal spline(edge[] Edges)
            : base()
        {
            this.Edges.AddRange(Edges);
            for (int i = 0; i < Edges.Length; i++)
            {
                if (!Vertices.Contains(Edges[i].Vertices[0])) Vertices.Add(Edges[i].Vertices[0]);
                if (!Vertices.Contains(Edges[i].Vertices[1])) Vertices.Add(Edges[i].Vertices[1]);

            }
        }

        //**METHODS**CREATE**
        /// <summary>
        /// creates spline from HalfEdgeList
        /// </summary>
        /// <param name="Edges"></param>
        /// <returns></returns>
        public static spline ByHalfEdges(edge[] Edges)
        {
            return new spline(Edges);
        }
    }

    /// <summary>
    /// Face: HalfEdge List with references to edges, vertices, center, name, angle, normal, midpoint, and length
    /// </summary>
    public class face
    {
        //**FIELDS**
        internal bool disposed = false;
        internal List<halfedge> E;

        //**PROPERTIES** //**QUERY**
        /// <summary>
        /// Name: face label
        /// </summary>
        public string Name { get; set; }
        /// <summary>
        /// CS: context coordinate system
        /// </summary>
        public CoordinateSystem CS
        {
            get
            {
                vector z = Normal.Normalized();
                vector x = E[0].Normalized();
                vector y = z.NormalizedCross(x);
                x = y.NormalizedCross(z);
                Point O = Point.ByCoordinates(Center.X,Center.Y, Center.Z);
                Vector X = Vector.ByCoordinates(x.X, x.Y, x.Z);
                Vector Y = Vector.ByCoordinates(y.X, y.Y, y.Z);
                Vector Z = Vector.ByCoordinates(z.X, z.Y, z.Z);
                CoordinateSystem CS = CoordinateSystem.ByOriginVectors(O,X,Y,Z);
                O.Dispose(); X.Dispose();  Y.Dispose(); Z.Dispose();
                return CS;
            }
        }
        /// <summary>
        /// Center: centerpoint of face
        /// </summary>
        public point Center { get; set; }
        /// <summary>
        /// Normal: face normal that defines outside
        /// </summary>
        public vector Normal { get; set; }
        /// <summary>
        /// HalfEdges: halfedge list
        /// </summary>
        public List<halfedge> HalfEdges { get { return E; } }
        /// <summary>
        /// Vertices: veterx array
        /// </summary>
        public vertex[] Vertices
        {
            get
            {
                vertex[] output = new vertex[E.Count];
                for (int i = 0; i < E.Count; i++) output[i] = E[i].V[0];
                return output;
            }
        }
        /// <summary>
        /// Edges: associated edge array
        /// </summary>
        public edge[] Edges
        {
            get
            {
                edge[] output = new edge[E.Count];
                for (int i = 0; i < E.Count; i++) output[i] = E[i].Edge;
                return output;
            }
        }
        /// <summary>
        /// VertexPoints: array of vertex point geometry 
        /// </summary>
        public Point[] VertexPoints
        {
            get
            {
                Point[] output = new Point[E.Count];
                for (int i = 0; i < E.Count; i++) output[i] = E[i].V[0].Point;
                return output;
            }
        }
        /// <summary>
        /// Vertex Vector Array:
        /// V1 is normalized vector from vertex in right-hand rule
        /// V2 is normalized vector from vertex in other direction
        /// N is face normal ie. V1 x V2 (cross product).
        /// returns array {V1, V2, V1+V2 (vertex bisector), V1-V2 (ON vector to biscetor and normal), NxV1 + V2xN(exterior bisector)}
        /// </summary>
        public vector[][] VertexVectors {get; internal set; }
        /// <summary>
        /// Angles: interior angles for face in radians
        /// </summary>
        public double[] Angles { get; set; }
        /// <summary>
        /// MinEdgeAngle: sharpest corner made with adjacent face
        /// </summary>
        public double MinEdgeAngle
        {
            get
            {
                double min = 2*math.PI;
                for (int i = 0; i < E.Count; i++) if (min > E[i].Angle) min = E[i].Angle;
                return min;
            }
        }

        //**CONSTRUCTOR**
        /// <summary>
        /// default face constructor
        /// </summary>
        public face() { Name = "";}
        /// <summary>
        /// face constructor from vertex array/list and normal
        /// </summary>
        public face(IEnumerable<vertex> Vertices, vector Normal)
            : this()
        {
            this.Normal = Normal;
            vertex[] V = Vertices.ToArray();
            VertexVectors = new vector[V.Length][];
            E = new List<halfedge>(Vertices.ToList().Count);
            for (int i = 0; i < V.Length; i++)
            {
                E.Add(new halfedge(V[i], V[(i + 1) % V.Length]));
                V[i].AddFace(this as face);
            }
            E.ForEach(he => he.AddFace(this));
            double[] xyz = { 0, 0, 0 };
            for (int i = 0; i < E.Count; i++)
            {
                xyz[0] += V[i].X / E.Count;
                xyz[1] += V[i].Y / E.Count;
                xyz[2] += V[i].Z / E.Count;
                int j = (i + E.Count - 1) % E.Count;
                VertexVectors[i] = vector.NormalizedVertexVectors(E[i], E[j].Reverse(), Normal);
            }
            Center = point.ByCoordinates(xyz[0], xyz[1], xyz[2]);
            Angles = new double[E.Count];
            for (int i = 0; i < E.Count; i++)
            {
                arc arc = arc.ByThreePoints(VertexVectors[i][0], VertexVectors[i][4], VertexVectors[i][1]);
                Angles[i] = arc.SweepAngle;
            }
        }

        //**METHODS**CREATE
        /// <summary>
        /// creates face from ordered vertices and normal
        /// </summary>
        /// <param name="Vertices"></param>
        /// <param name="Normal"></param>
        /// <returns></returns>
        public static face ByVertices(IEnumerable<vertex> Vertices, vector Normal) { return new face(Vertices, Normal); }

        //**METHODS** //**ACTION**
        /// <summary>
        /// reorders vertices based on given start vertex
        /// </summary>
        /// <param name="Start"></param>
        /// <returns></returns>
        public face ReOrderVertices(vertex Start)
        {
            if (E[0].V[0].Equals(Start)) return null;
            int index = 0;
            halfedge[] temp = new halfedge[E.Count];
            E.CopyTo(temp);
            for (int i = 1; i < E.Count; i++) if (E[i].Equals(Start)) { index = i; break; }
            for (int i = 0; i < E.Count; i++) E[i] = temp[(index + 1) % E.Count];
            temp = null;
            return this;
        }
        /// <summary>
        /// rename face
        /// </summary>
        /// <param name="Name">string</param>
        public void Rename(string Name)
        {
            this.Name = Name;
        }
    }
    /// <summary>
    /// triangle mesh face
    /// </summary>
    public class triangle : face
    {
        //**CONSTRUCTOR**
        /// <summary>
        /// default constructor
        /// </summary>
        public triangle() : base() { }
        /// <summary>
        /// constructor from vertex array/list and normal
        /// </summary>
        public triangle(IEnumerable<vertex> Vertices, vector Normal) : base(Vertices.Take(3), Normal) { }
        /// <summary>
        /// creates triangle mesh face from given vertex array/list and given normal
        /// </summary>
        public static triangle ByVerticesNormal(IEnumerable<vertex> Vertices, vector Normal) { return new triangle(Vertices, Normal); }
        /// <summary>
        /// creates triangle mesh face from given vertex array/list
        /// </summary>
        public static triangle ByVertices(IEnumerable<vertex> Vertices)
        {
            vector X = vector.ByTwoPoints(Vertices.ElementAt(0), Vertices.ElementAt(1));
            vector Y = vector.ByTwoPoints(Vertices.ElementAt(0), Vertices.ElementAt(2));
            vector N = X.Cross(Y);
            return new triangle(Vertices, N);
        }
        /// <summary>
        /// gets circumcenter
        /// </summary>
        public point GetCircumcenter()
        {
            arc c = arc.ByThreePoints(Vertices.ElementAt(0), Vertices.ElementAt(1), Vertices.ElementAt(2));
            return c.Center;
        }
        /// <summary>
        /// gets incenter
        /// </summary>
        public point GetIncenter()
        {
            double D = E[0].Length + E[1].Length + E[2].Length;
            double X = E[1].Length * E[0].V[0].X + E[2].Length * E[1].V[0].X + E[0].Length * E[2].V[0].X;
            double Y = E[1].Length * E[0].V[0].Y + E[2].Length * E[1].V[0].Y + E[0].Length * E[2].V[0].Y;
            double Z = E[1].Length * E[0].V[0].Z + E[2].Length * E[1].V[0].Z + E[0].Length * E[2].V[0].Z;
            return point.ByCoordinates(X / D, Y / D, Z / D);
        }
        /// <summary>
        /// gets other vertex
        /// </summary>
        public vertex GetOtherVertex(edge Edge)
        {
            List<vertex> V = new List<vertex>(Vertices);
            V.RemoveAll(v => Edge.Vertices.Contains(v));
            if (V.Count > 1) return null;
            return V[0];
        }
    }
    /// <summary>
    /// quad mesh face
    /// </summary>
    public class quad : face
    {
        //**PROPERTIES** //**QUERY**
        /// <summary>
        /// diagonal
        /// </summary>
        public int Diagonal { get; private set; }

        //**CONSTRUCTOR**
        internal quad(IEnumerable<vertex> Vertices, vector Normal)
            : base(Vertices.Take(4), Normal)
        {
            Diagonal = 0;
            if (E[0].V[0].DistanceTo(E[2].V[0]) > E[1].V[0].DistanceTo(E[3].V[0])) Diagonal = 1;
        }

        //**METHODS**CREATE
        /// <summary>
        /// creates quad object with given vertex list/array and given normal
        /// </summary>
        new public static quad ByVertices(IEnumerable<vertex> Vertices, vector Normal) { return new quad(Vertices, Normal); }

        //**METHODS**ACTION
        /// <summary>
        /// changes vertex ordering to flip diagonal
        /// </summary>
        public void FlipDiagonal()
        {
            Diagonal = (Diagonal + 1) % 2;
        }
        /// <summary>
        /// returns diagonal as Dynamo Line object
        /// </summary>
        public Line GetDiagonal()
        {
            Point a = E[Diagonal].V[0].Point;
            Point b = E[Diagonal + 2].V[0].Point;
            Line L = Line.ByStartPointEndPoint(a, b);
            a.Dispose(); b.Dispose();
            return L;
        }
    }


    /// <summary>
    /// Mesh Object with Topology and Geometry
    /// </summary>
    public class mesh
    {
        //**FIELDS**
        internal Dictionary<point, vertex> V;
        internal Dictionary<string, edge> E;
        internal List<edge> E1;
        internal List<edge> E2;
        internal List<edge> E3;
        internal List<edge> E0;
        internal int Df = 0;
        internal int De = 0;
        internal double MinFaceAngle = math.PI;

        //**PROPERTIES**QUERY
        /// <summary>
        /// Face List
        /// </summary>
        public List<face> Faces { get; set; }
        /// <summary>
        /// Edge List
        /// </summary>
        public List<edge> Edges { get; set; }
        /// <summary>
        /// Naked Edge List
        /// </summary>
        public List<edge> EdgesOuter { get { return E1; } }
        /// <summary>
        /// Vertex List
        /// </summary>
        public List<vertex> Vertices { get { return V.Values.ToList(); } }
        /// <summary>
        /// Spline List
        /// </summary>
        public List<spline> Splines { get; set; }
        /// <summary>
        /// Point List (Vertex Geometry)
        /// </summary>
        public point[] Points { get; set; }

        //**CONSTRUCTOR**
        internal mesh()
        {
            V = new Dictionary<point, vertex>();
            E = new Dictionary<string, edge>();
            Edges = new List<edge>();
            E1 = new List<edge>();
            E2 = new List<edge>();
            E3 = new List<edge>();
            E0 = new List<edge>();
            Splines = new List<spline>();
        }
        internal mesh(point[] Points)
            : this()
        {
            De = Points.Length;
            V = new Dictionary<point, vertex>(Points.Length);
            // create vertex lookup table from points
            for (int i = 0; i < Points.Length; i++)
            {
                if (V.ContainsKey(Points[i])) continue;
                V.Add(Points[i], new vertex(Points[i]));
            }
            this.Points = V.Keys.ToArray();
        }
        internal mesh(point[] Points, Surface[] Surfaces)
            : this(Points)
        {
            Df = Surfaces.Length.ToString().Length;
            De = (Surfaces.Length + Points.Length).ToString().Length;
            // initialize
            Faces = new List<face>(Surfaces.Length);
            //**CONSTRUCT MESH
            // create faces from surfaces
            int eCount = 1;
            for (int i = 0; i < Surfaces.Length; i++)
            {
                // find face vertices in lookup table
                // create face based on vertices
                if (Surfaces[i].Vertices.Length == 3)
                {
                    Vector N = Surfaces[i].NormalAtParameter(0.5, 0.5);
                    vector n = vector.ByCoordinates(N.X,N.Y,N.Z);
                    N.Dispose();
                    triangle t = triangle.ByVerticesNormal(FindFaceVertices(Surfaces[i]), n);
                    t.Name = "t" + (i + 1).ToString("D" + Df);
                    Faces.Add(t);
                    eCount = FindFaceEdges(t, eCount);
                }
                else
                {
                    Vector N = Surfaces[i].NormalAtParameter(0.5, 0.5);
                    vector n = vector.ByCoordinates(N.X, N.Y, N.Z);
                    N.Dispose();
                    face f = face.ByVertices(FindFaceVertices(Surfaces[i]), n);
                    f.Name = "f" + (i + 1).ToString("D" + Df);
                    Faces.Add(f);
                    eCount = FindFaceEdges(f, eCount);
                }
            }
            //**ITERATE THROUGH EDGES FOR ANGLE CALCS
            CalculateAngles();
        }

        //**METHODS**CREATE
        /// <summary>
        /// creates mesh object with given surfaces and derives vertices
        /// surfaces must be planar but can be concave as well as polygonal
        /// </summary>
        /// <param name="Surfaces">Surface Array</param>
        /// <returns>Mesh Object</returns>
        public static mesh BySurfaces(Surface[] Surfaces)
        {
            List<point> Points = new List<point>();
            for (int i = 0; i < Surfaces.Length; i++)
            {
                for (int j = 0; j < Surfaces[i].Vertices.Length; j++)
                {
                    Point P = Surfaces[i].Vertices[j].PointGeometry;
                    point p = point.ByCoordinates(P.X, P.Y, P.Z);
                    P.Dispose();
                    if (Points.Contains(p)) continue;
                    Points.Add(p);
                }
            }
            return new mesh(Points.ToArray(), Surfaces);
        }

        //**METHODS**ACTIONS
        /// <summary>
        /// returns topological entity (vertex) at point
        /// </summary>
        /// <param name="Point">Vertex Geometry</param>
        /// <returns>Vertex</returns>
        public vertex GetVertexAtPoint(Point Point)
        {
            point p = point.ByCoordinates(Point.X, Point.Y, Point.Z);
            vertex Result = null;
            for (int k = 0; k < V.Keys.Count; k++) if (p.Equals(Points[k])) { Result = V[Points[k]]; break; }
            return Result;
        }
        /// <summary>
        /// returns topological entity (edge) at line
        /// </summary>
        /// <param name="Line">Edge Geometry</param>
        /// <returns>Edge</returns>
        public edge GetEdgeAtLine(Curve Line)
        {
            edge Result = null;
            for (int k = 0; k < Edges.Count; k++) if (Edges[k].IsAtCurve(Line)) { Result = Edges[k]; break; }
            return Result;
        }
        /// <summary>
        /// calculates edge angles based on input geometry
        /// </summary>
        public void CalculateAngles()
        {
            //**ITERATE THROUGH EDGES FOR ANGLE CALCS
            for (int i = 0; i < Edges.Count; i++)
            {
                // case 0: is OuterEdge
                if (Edges[i].E.Count == 1) E1.Add(Edges[i]);
                // case 1: is toplogically consistent edge
                else if (Edges[i].E.Count == 2)
                {
                    double[] AngleNormal = Edges[i].GetAngleNormal(Edges[i].E.ElementAt(0), Edges[i].E.ElementAt(1));
                    Edges[i].Angle = new double[] { AngleNormal[0] };
                    Edges[i].N = new double[] { AngleNormal[1], AngleNormal[2], AngleNormal[3] };
                    Edges[i].E.ElementAt(0).Angle = AngleNormal[0];
                    Edges[i].E.ElementAt(1).Angle = AngleNormal[0];
                    E2.Add(Edges[i]);
                }
                // case 2: topologically inconsistent
                else if (Edges[i].E.Count == 3)
                {
                    vertex A = Edges[i].Vertices[0];
                    vertex B = Edges[i].Vertices[1];
                    List<halfedge> eA = new List<halfedge>(2);
                    List<halfedge> eB = new List<halfedge>(2);
                    for (int j = 0; j < 3; j++)
                    {
                        if (Edges[i].E.ElementAt(j).V[0].Equals(A)) eA.Add(Edges[i].E.ElementAt(j));
                        else if (Edges[i].E.ElementAt(j).V[0].Equals(B)) eB.Add(Edges[i].E.ElementAt(j));
                    }
                    if (eA.Count + eB.Count != 3) continue;
                    halfedge e0, e1, e2;
                    List<halfedge> e;
                    if (eA.Count == 1) { e0 = eA[0]; e = eB; } else { e0 = eB[0]; e = eA; }
                    double[] n0 = Edges[i].GetAngleNormal(e0, e[0]);
                    double[] n1 = Edges[i].GetAngleNormal(e0, e[1]);
                    double[] a1, a2;
                    if (n0[0] < n1[0]) { e1 = e[0]; e2 = e[1]; a1 = n0; a2 = n1; }
                    else { e1 = e[1]; e2 = e[0]; a1 = n1; a2 = n0; }
                    vector e1N = e1.Face.Normal, e2N = e2.Face.Normal;
                    vector a3 = e1N.Reverse().NormalizedAdd(e2N);
                    Edges[i].E = new List<halfedge> { e0, e1, e2 };
                    Edges[i].Angle = new double[] { a1[0], a2[0] - a1[0], 2*math.PI - a2[0] };
                    Edges[i].N = new double[] { -a1[1], -a1[2], -a1[3], -a3.X, -a3.Y, -a3.Z, a2[1], a2[2], a2[3] };
                    Edges[i].E.ElementAt(0).Angle = a1[0];
                    Edges[i].E.ElementAt(1).Angle = Math.Min(a1[0], a2[0] - a1[0]);
                    Edges[i].E.ElementAt(2).Angle = a2[0] - a1[0];
                    E3.Add(Edges[i]);
                }
                else E0.Add(Edges[i]);
            }
        }
        /// <summary>
        /// looks for Face Edge in Mesh Edge List,
        /// adds topological information if found
        /// otherwise creates Edge and adds to Mesh Edge List
        /// </summary>
        /// <param name="Face">Face</param>
        /// <param name="EdgeCount">ExistingEdgeCount</param>
        /// <returns>NewEdgeCount</returns>
        public int FindFaceEdges(face Face, int EdgeCount)
        {
            int eCount = EdgeCount;
            for (int j = 0; j < Face.E.Count; j++)
            {
                bool edgeFound = false;
                if (Edges.Count > 0)
                {
                    for (int k = 0; k < Edges.Count; k++)
                    {
                        if (Edges[k].Vertices.Contains(Face.E[j].V[0]) && Edges[k].Vertices.Contains(Face.E[j].V[1]))
                        {
                            Edges[k].E.Add(Face.E[j]);
                            Face.E[j].AddEdge(Edges[k]);
                            edgeFound = true;
                            break;
                        }
                    }
                }
                if (!edgeFound)
                {
                    edge e = new edge();
                    e.Name = eCount.ToString("D" + De);
                    e.E.Add(Face.E[j]);
                    Face.E[j].AddEdge(e);
                    Edges.Add(e);
                    eCount++;
                }
            }
            return eCount;
        }
        /// <summary>
        /// looks for Face Vertex in Mesh Vertex List,
        /// adds topological information if found
        /// otherwise creates Vertex and adds to Mesh Vertex List
        /// </summary>
        /// <param name="Surface"></param>
        /// <returns></returns>
        public List<vertex> FindFaceVertices(Surface Surface)
        {
            Autodesk.DesignScript.Geometry.Vertex[] vtx = Surface.Vertices;
            List<vertex> v = new List<vertex>(vtx.Length);
            for (int j = 0; j < vtx.Length; j++)
            {
                Point pt = vtx[j].PointGeometry;
                vertex search = GetVertexAtPoint(pt);
                point p = point.ByCoordinates(pt.X, pt.Y, pt.Z);
                pt.Dispose();
                if (!search.Equals(null)) v.Add(search);
                else
                {
                    V.Add(p,vertex.ByPoint(p));
                    v.Add(GetVertexAtPoint(pt));
                }
            }
            vtx.ForEach(x => x.Dispose());
            return v;
        }
        /// <summary>
        /// flips face normals
        /// </summary>
        /// <returns>Mesh</returns>
        public mesh FlipFaceNormals()
        {
            Faces.ForEach(f => f.E.ForEach(e => e.FlipDirection()));
            return this;
        }

        //**METHODS**IN PROGRESS
        internal mesh AddEdgeNames(Point Point, int SplineCount = 2, int EdgeCount = 3)
        {
            // find start point
            vertex v = GetVertexAtPoint(Point);
            if (v.Equals(null)) return null;
            List<edge> E0 = new List<edge>(Edges);
            spline S0 = new spline(); S0.Vertices.Add(v);
            HashSet<edge> vE = new HashSet<edge>(v.Edges); vE.IntersectWith(E1);
            halfedge e = null;
            for (int i = 0; i < vE.Count; i++) if (vE.ElementAt(i).E.ElementAt(0).V[0].Equals(v)) e = vE.ElementAt(i).E.ElementAt(0);
            if (e.Equals(null)) return null;
            int s = 0;
            int sN = 1;
            while (vE.Count > 0)
            {
                if (sN > 1) e = vE.ElementAt(0).E.ElementAt(0);
                e.Face.ReOrderVertices(S0.Vertices[sN - 1]);
                e.Edge.Name = "s" + s.ToString("D" + SplineCount) + "-" + sN.ToString("D" + EdgeCount);
                E.Add(e.Edge.Name, e.Edge);
                E0.Remove(e.Edge);
                S0.Edges.Add(e.Edge);
                vertex vS = e.Edge.GetOtherVertex(S0.Vertices[sN - 1]);
                if (!S0.Vertices.Contains(vS)) S0.Vertices.Add(vS);
                vE = new HashSet<edge>(S0.Vertices[sN].Edges);
                vE.IntersectWith(E1);
                vE.IntersectWith(E0);
                sN++;
            }
            Splines.Add(S0);
            s++;
            // splines
            while (E0.Count > 0 && this.E.Count < Edges.Count)
            {
                for (int i = 0; i < Splines[s].Edges.Count; i++)
                {
                    spline S = new spline();
                    sN = 1;
                    int eN = 1;
                    int tN = 1;
                    List<face> T0 = Splines[s].Edges[i].Faces.FindAll(f => f.Name.Equals(""));
                    if (T0.Count > 1) { }
                    S.Vertices.Add(T0[0].Vertices[2]);
                    for (int j = T0[0].E.Count - 1; j > 0; j--) if (T0[0].E[i].Edge.Name.Equals(""))
                        {
                            T0[0].E[i].Edge.Name = "e" + s.ToString("D" + SplineCount) + (s + 1).ToString("D" + SplineCount) + "-" + eN.ToString("D" + EdgeCount);
                            E.Add(T0[0].E[i].Edge.Name, T0[0].E[i].Edge);
                            E0.Remove(T0[0].E[i].Edge);
                            eN++;
                        }
                    T0[0].Name = "t" + s.ToString("D" + SplineCount) + (s + 1).ToString("D" + SplineCount) + "-" + tN.ToString("D" + EdgeCount);
                    tN++;
                    List<face> T1 = T0[0].E[1].Edge.Faces.FindAll(f => f.Name.Equals(""));
                    if (T0.Count > 1) { }
                    vertex v0 = (T1[0] as triangle).GetOtherVertex(T0[0].E[1].Edge);
                    T1[0].ReOrderVertices(v0);
                    T1[0].E[0].Edge.Name = "s" + s.ToString("D" + SplineCount) + "-" + sN.ToString("D" + EdgeCount);
                    E.Add(T1[0].E[0].Edge.Name, T1[0].E[0].Edge);
                    E0.Remove(T1[0].E[0].Edge);
                    S.Edges.Add(T1[0].E[0].Edge);
                    S.Vertices.Add(e.Edge.GetOtherVertex(S.Vertices[sN - 1]));
                    sN++;
                }

            }
            return this;
        }
        internal mesh GetSplineData(Curve[] Start, Curve[] End, Curve[] Lines)
        {
            List<edge> E = new List<edge>(Start.Length);
            List<Curve> Search = new List<Curve>(Lines);
            Search.RemoveAll(c => Start.Contains(c));
            Search.RemoveAll(c => End.Contains(c));
            for (int i = 0; i < Search.Count; i++)
            {
                edge e = GetEdgeAtLine(Search[i]);
                if (e.Equals(null))
                {
                    vertex v1 = GetVertexAtPoint(Lines[i].StartPoint);
                    vertex v2 = GetVertexAtPoint(Lines[i].EndPoint);
                    e = new edge(new HashSet<halfedge> { new halfedge(v1, v2), new halfedge(v2, v1) });
                    e.E[0].AddEdge(e);
                    e.E[1].AddEdge(e);
                    v1.Edges.Add(e);
                    v2.Edges.Add(e);
                }
                E.Add(e);

            }
            // case zero
            for (int i = 0; i < Start.Length; i++)
            {
                List<vertex> V = new List<vertex>();
                edge e0 = GetEdgeAtLine(Start[i]);
                V.Add(e0.E[0].V[0]);
                V.Add(e0.E[0].V[1]);
                HashSet<edge> H = V[V.Count - 1].Edges;
                H.IntersectWith(E);
                V.Add(H.ElementAt(0).GetOtherVertex(V[V.Count - 1]));
                E.Remove(H.ElementAt(0));



            }
            return this;
        }
    }
}
