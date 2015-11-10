using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using Text;

namespace Topology
{

    public class HalfEdge
    {
        internal List<Vertex> V;
        //**PROPERTIES** //**QUERY**
        public Edge Edge { get; private set; }
        public Face Face { get; private set; }

        //**CONSTRUCTOR**
        internal HalfEdge(Vertex A, Vertex B)
        {
            V = new List<Vertex>{A,B};
            Edge = null; Face = null; 
        }
        internal HalfEdge(Vertex A, Vertex B, Edge Edge, Face Face) : this(A,B)
        {this.Edge = Edge; this.Face = Face; }
        internal HalfEdge(IEnumerable<Vertex> Vertices)
        {V = Vertices.ToList(); this.Edge = null; Face = null; }
        internal HalfEdge(IEnumerable<Vertex> Vertices, Edge Edge, Face Face)
        {V = Vertices.ToList(); this.Edge = Edge; this.Face = Face; }

        //**METHODS** //**ACTION**
        public Vector GetVector()
        {
            Point A = V[0].Point;
            Point B = V[1].Point;
            Vector output = Vector.ByTwoPoints(A, B);
            A.Dispose(); B.Dispose();
            return output;
        }
        public HalfEdge FlipDirection()
        {
            List<Vertex> temp = new List<Vertex>(V);
            V[0] = temp[1];
            V[1] = temp[0];
            temp = null;
            return this;
        }
        public bool AddEdge(Edge Edge)
        {
            if (this.Edge == null)
            {
                this.Edge = Edge;
                return true;
            }
            return false;
        }
        public bool AddFace(Face Face)
        {
            if (this.Face == null)
            {
                this.Face = Face;
                return true;
            }
            return false;
        }
    }

    public class Vertex: IEquatable<Vertex>
    {
        //**PROPERTIES** //**QUERY**
        public List<double> Coordinates { get; private set; }
        public double X { get { return Coordinates[0]; } }
        public double Y { get { return Coordinates[1]; } }
        public double Z { get { return Coordinates[2]; } }
        public Point Point { get { return Point.ByCoordinates(X, Y, Z); } }
        public HashSet<Edge> Edges { get; private set; }
        public HashSet<Face> Faces { get; private set; }

        //**CONSTRUCTOR**
        internal Vertex(Point Point)
        {
            Coordinates = new List<double> { Point.X, Point.Y, Point.Z };
            Edges = new HashSet<Edge>();
            Faces = new HashSet<Face>();
        }
        internal Vertex(Point Point, IEnumerable<Edge> Edges, IEnumerable<Face> Faces)
        {
            Coordinates = new List<double> { Point.X, Point.Y, Point.Z };
            this.Edges = new HashSet<Edge>(Edges);
            this.Faces = new HashSet<Face>(Faces);
        }
        internal Vertex(double X, double Y, double Z)
        {
            Coordinates = new List<double> { X, Y, Z };
            Edges = new HashSet<Edge>();
            Faces = new HashSet<Face>();
        }
        internal Vertex(double X, double Y, double Z, IEnumerable<Edge> Edges, IEnumerable<Face> Faces)
        {
            Coordinates = new List<double> { X, Y, Z };
            this.Edges = new HashSet<Edge>(Edges);
            this.Faces = new HashSet<Face>(Faces);
        }

        //**METHODS** //**ACTION**
        [MultiReturn(new[] { "X", "Y", "Z" })]
        public Dictionary<string, double> GetCoordinates() { return new Dictionary<string, double> { { "X", X }, { "Y", Y }, { "Z", Z } }; }
        public void AddEdge(Edge Edge) { if (Edge.Vertices.Contains(this)) Edges.Add(Edge); }
        public void AddEdges(IEnumerable<Edge> Edges) { for (int i = 0; i < Edges.Count(); i++) AddEdge(Edges.ElementAt(i)); }
        public void AddFace(Face Face) { if (Face.Vertices.Contains(this)) Faces.Add(Face); }
        public void AddFaces(IEnumerable<Face> Faces) { for (int i = 0; i < Faces.Count(); i++) AddFace(Faces.ElementAt(i)); }

        public override bool Equals(Object Object) { return this.Equals(Object as Vertex); }
        public bool Equals(Vertex Vertex)
        {
            if (Object.ReferenceEquals(Vertex, null)) return false;
            if (Object.ReferenceEquals(this, Vertex)) return true;
            if (this.GetType() != Vertex.GetType()) return false;
            return (X == Vertex.X && Y == Vertex.Y && Z == Vertex.Z);
        }
        public bool IsAtPoint(Point Point) { return (X == Point.X && Y == Point.Y && Z == Point.Z); }
        public override int GetHashCode() { return string.Format("{0}-{1}-{2}", X, Y, Z).GetHashCode(); }
        public static bool operator ==(Vertex a, Vertex b)
        {
            if (Object.ReferenceEquals(a, null))
            {
                if (Object.ReferenceEquals(b, null)) return true;
                return false;
            }
            return a.Equals(b);
        }
        public static bool operator !=(Vertex a, Vertex b) { return !(a == b); }
    }


    public class Edge
    {
        //**FIELDS**
        internal HashSet<HalfEdge> E;

        //**PROPERTIES** //**QUERY**
        public string Name { get; set; }
        public double[] Angle { get; set; }
        public Point MidPoint { get { return Point.ByCoordinates(E.ToList()[0].V[0].X/2 + E.ToList()[0].V[1].X/2, E.ToList()[0].V[0].Y/2 + E.ToList()[0].V[1].Y/2, E.ToList()[0].V[0].Z/2 + E.ToList()[0].V[1].Z/2); } }
        public List<Face> Faces { get { List<Face> F = new List<Face>(E.Count); E.ToList().ForEach(e => F.Add(e.Face)); return F; } }
        public List<Vertex> Vertices { get { return E.ElementAt(0).V; } }
        //**CONSTRUCTOR**
        internal Edge() : base() { E = new HashSet<HalfEdge>(); Name = ""; Angle = null; }
        internal Edge(IEnumerable<HalfEdge> HalfEdges) { E = new HashSet<HalfEdge>(HalfEdges); Name = ""; Angle = null; }
        internal Edge(IEnumerable<HalfEdge> HalfEdges, string Name) { E = new HashSet<HalfEdge>(HalfEdges); this.Name = Name; Angle = null; }

        //**METHODS** //**ACTION**
        public Vertex GetOtherVertex(Vertex Vertex)
        {
            if (Vertices[0].Equals(Vertex)) return Vertices[1];
            if (Vertices[1].Equals(Vertex)) return Vertices[0];
            return null;
        }
        public Line GetLine()
        {
            Line output;
            using (Point a = Vertices[0].Point)
            {
                using (Point b = Vertices[1].Point)
                {
                    output = Line.ByStartPointEndPoint(a, b);
                }
            }
            return output;
        }
        public bool IsAtCurve(Curve Line)
        {
            Line ln = GetLine();
            if (ln.IsAlmostEqualTo(Line))
            {
                ln.Dispose();
                return true;
            }
            ln.Dispose();
            return false;
        }
    }

    public class Face: IDisposable
    {
        //**FIELDS**
        private bool disposed = false;
        internal List<HalfEdge> E;

        //**PROPERTIES** //**QUERY**
        public string Name { get; set; }
        public CoordinateSystem CS { get; set; }
        public Point Center { get { return CS.Origin; } }
        public Vector Normal { get { return CS.ZAxis; } }
        public Dictionary<string, Object> Parameters { get; set; }

        public Vertex[] Vertices
        {
            get
            {
                Vertex[] output = new Vertex[E.Count];
                for (int i = 0; i < E.Count; i++) output[i] = E[i].V[0];
                return output;
            }
        }
        public Edge[] Edges
        {
            get
            {
                Edge[] output = new Edge[E.Count];
                for (int i = 0; i < E.Count; i++) output[i] = E[i].Edge;
                return output;
            }
        }

        //**CONSTRUCTOR**
        internal Face() : base() { }
        internal Face(IEnumerable<Vertex> Vertices)
        {
            E = new List<HalfEdge>(Vertices.ToList().Count);
            for (int i = 0; i < Vertices.Count(); i++)
            {
                E.Add(new HalfEdge(Vertices.ElementAt(i), Vertices.ElementAt((i + 1) % E.Capacity)));
                Vertices.ElementAt(i).AddFace(this as Face);
            }
            E.ForEach(he => he.AddFace(this));
            double[] xyz = { 0, 0, 0 };
            for (int i = 0; i < E.Count; i++)
            {
                xyz[0] += Vertices.ElementAt(i).X / E.Count;
                xyz[1] += Vertices.ElementAt(i).Y / E.Count;
                xyz[2] += Vertices.ElementAt(i).Z / E.Count;
            }
            Point Center = Point.ByCoordinates(xyz[0], xyz[1], xyz[2]);
            SetCS(Center);
            Center.Dispose();
        }

        //**METHODS** //**ACTION**
        public Face ReOrderVertices(Vertex Start)
        {
            if (E[0].V[0].Equals(Start)) return null;
            int index = 0;
            HalfEdge[] temp = new HalfEdge[E.Count];
            E.CopyTo(temp);
            for (int i = 1; i < E.Count; i++) if (E[i].Equals(Start)) { index = i; break; }
            for (int i = 0; i < E.Count; i++) E[i] = temp[(index + 1) % E.Count];
            temp = null;
            SetCS(CS.Origin);
            return this;
        }
        internal void SetCS(Point Center)
        {
            Vector Y = E[E.Count - 1].GetVector();
            Y = Y.Reverse();
            Vector X = E[0].GetVector();
            Vector Z = X.Cross(Y);
            Y = Z.Cross(X);
            CS = CoordinateSystem.ByOriginVectors(Center, X, Y);
            X.Dispose(); Y.Dispose(); Z.Dispose();
        }
        public virtual void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                E.ForEach(e => e.Edge.E.Remove(e));
                E.Clear();
                if (Center != null) Center.Dispose();
                if (Normal != null) Normal.Dispose();
                if (CS != null) CS.Dispose();
                if (Parameters != null) for (int i = 0; i < Parameters.Count; i++)
                        if (Parameters.Values.ToArray()[i] is IDisposable) ((IDisposable)Parameters.Values.ToArray()[i]).Dispose();
            }
            disposed = true;
        }
    }

    public class Triangle : Face, IDisposable
    {
        //**FIELDS**
        private bool disposed = false;

        //**PROPERTIES** //**QUERY**
        public Point Circumcenter { get; private set; }

        //**CONSTRUCTOR**
        internal Triangle(IEnumerable<Vertex> Vertices)
            : base(Vertices)
        {
            Point[] pts = { Vertices.ElementAt(0).Point, Vertices.ElementAt(1).Point, Vertices.ElementAt(2).Point };
            Circle c = Circle.ByBestFitThroughPoints(pts);
            Circumcenter = c.CenterPoint;
            c.Dispose(); pts[0].Dispose(); pts[1].Dispose(); pts[2].Dispose(); pts = null;
        }

        //**METHODS** //**ACTION**
        public Point[] GetVertexPoints()
        {
            if (Vertices.Length < 3) { return null; }
            else
            {
                Point[] Points = { Vertices[0].Point, Vertices[1].Point, Vertices[2].Point };
                return Points;
            }
        }
        public override void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected new virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                if (Circumcenter != null) Circumcenter.Dispose();
                base.Dispose();
            }
            disposed = true;
        }
    }

    public class Quad : Face
    {
        //**FIELDS**
        private bool disposed = false;

        public override void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected new virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                base.Dispose();
            }
            disposed = true;
        }
    }

    public class Mesh : IDisposable
    {
        //**FIELDS**
        private bool disposed = false;
        private Dictionary<Point, Vertex> V;
        private Dictionary<string, Edge> E;

        //**PROPERTIES**QUERY
        public List<Face> Faces { get; set; }
        public List<Vertex> Vertices { get { return V.Values.ToList(); } }
        public List<Edge> Edges { get; set; }
        public Point[] Points { get; set; }

        internal Mesh(Surface[] Surfaces, Point[] Points)
        {
            // initialize
            Faces = new List<Face>(Surfaces.Length);
            V = new Dictionary<Point, Vertex>(Points.Length);
            E = new Dictionary<string, Edge>();
            Edges = new List<Edge>();

            // store input points
            this.Points = Points;
            // create vertex lookup table from points
            for (int i = 0; i < Points.Length; i++) V.Add(Points[i], new Vertex(Points[i]));
            // create faces from surfaces
            for (int i = 0; i < Surfaces.Length; i++)
            {
                // find face vertices in lookup table
                Autodesk.DesignScript.Geometry.Vertex[] vtx = Surfaces[i].Vertices;
                List<Vertex> v = new List<Vertex>(vtx.Length);
                for (int j = 0; j < vtx.Length; j++)
                {
                    Point pt = vtx[j].PointGeometry;
                    for (int k = 0; k < Points.Length; k++)
                    {
                        if (pt.IsAlmostEqualTo(Points[k]))
                        {
                            v.Add(V[Points[k]]);
                            break;
                        }
                    }
                    pt.Dispose();
                }
                vtx.ForEach(x => x.Dispose());
                // create face based on vertices
                Triangle t = new Triangle(v);
                Faces.Add(t);
                // create or find edges
                for (int j = 0; j < t.E.Count; j++)
                {
                    bool edgeFound = false;
                    if (Edges.Count > 0)
                    {
                        for (int k = 0; k < Edges.Count; k++)
                        {
                            if (Edges[k].Vertices.Contains(t.E[j].V[0]) && Edges[k].Vertices.Contains(t.E[j].V[1]))
                            {
                                Edges[k].E.Add(t.E[j]);
                                t.E[j].AddEdge(Edges[k]);
                                edgeFound = true;
                                break;
                            }
                        }
                    }
                    if (!edgeFound)
                    {
                        Edge e = new Edge();
                        e.E.Add(t.E[j]);
                        t.E[j].AddEdge(e);
                        Edges.Add(e);
                    }
                }
            }
        }

        public static Mesh BySurfacesPoints(Surface[] Surfaces, Point[] Points) { return new Mesh(Surfaces, Points); }

        public Mesh AddEdgeNames(PolyCurve[] Spline, int Digits = 3)
        {
            int D = (Spline.Length).ToString().Length;
            for (int i = 0; i < Spline.Length; i++)
            {
                for (int k = 0; k < Spline[i].Curves().Length; k++)
                {
                    for (int j = 0; j < Edges.Count; j++)
                    {
                        if (Edges[j].IsAtCurve(Spline[i].Curves()[k]))
                        {
                            Edges[j].Name = "s" + (i + 1).ToString("D" + D) + "-" + (k + 1).ToString("D" + Digits);
                            E.Add(Edges[j].Name, Edges[j]);
                            break;
                        }
                    }
                }
            }
            return this;
        }

        public void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                Faces.ForEach(f => f.Dispose());
            }
            disposed = true;
        }
    }
}
