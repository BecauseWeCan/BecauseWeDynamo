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

    public class HalfEdge<T>: List<T> where T: Vertex
    {
        //**PROPERTIES** //**QUERY**
        public Edge<T> Edge { get; private set; }
        public Face<T> Face { get; private set; }

        //**CONSTRUCTOR**
        internal HalfEdge(T A, T B): base(2)
        {
            this.Add(A); this.Add(B); 
            Edge = null; Face = null; 
        }
        internal HalfEdge(T A, T B, Edge<T> Edge, Face<T> Face) : this(A,B)
        {this.Edge = Edge; this.Face = Face; }
        internal HalfEdge(IEnumerable<T> Vertices): base(Vertices)
        {this.Edge = null; Face = null; }
        internal HalfEdge(IEnumerable<T> Vertices, Edge<T> Edge, Face<T> Face): base(Vertices)
        {this.Edge = Edge; this.Face = Face; }

        //**METHODS** //**ACTION**
        public Vector GetVector()
        {
            Point A = this[0].GetPoint();
            Point B = this[1].GetPoint();
            Vector output = Vector.ByTwoPoints(A, B);
            A.Dispose(); B.Dispose();
            return output;
        }
        public HalfEdge<T> FlipDirection()
        {
            HalfEdge<T> temp = new HalfEdge<T>(this);
            this[0] = temp[1];
            this[1] = temp[0];
            temp = null;
            return this;
        }
        public bool AddEdge(Edge<T> Edge)
        {
            if (this.Edge == null)
            {
                this.Edge = Edge;
                return true;
            }
            return false;
        }
        public bool AddFace(Face<T> Face)
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
        public double X { get; private set; }
        public double Y { get; private set; }
        public double Z { get; private set; }
        public HashSet<Edge<Vertex>> Edges { get; private set; }
        public HashSet<Face<Vertex>> Faces { get; private set; }

        //**CONSTRUCTOR**
        internal Vertex(Point Point)
        {
            X = Point.X; Y = Point.Y; Z = Point.Z;
            Edges = new HashSet<Edge<Vertex>>();
            Faces = new HashSet<Face<Vertex>>();
        }
        internal Vertex(Point Point, IEnumerable<Edge<Vertex>> Edges, IEnumerable<Face<Vertex>> Faces)
        {
            X = Point.X; Y = Point.Y; Z = Point.Z;
            this.Edges = new HashSet<Edge<Vertex>>(Edges);
            this.Faces = new HashSet<Face<Vertex>>(Faces);
        }
        internal Vertex(double X, double Y, double Z)
        {
            this.X = X; this.Y = Y; this.Z = Z;
            Edges = new HashSet<Edge<Vertex>>();
            Faces = new HashSet<Face<Vertex>>();
        }
        internal Vertex(double X, double Y, double Z, IEnumerable<Edge<Vertex>> Edges, IEnumerable<Face<Vertex>> Faces)
        {
            this.X = X; this.Y = Y; this.Z = Z;
            this.Edges = new HashSet<Edge<Vertex>>(Edges);
            this.Faces = new HashSet<Face<Vertex>>(Faces);
        }

        //**METHODS** //**ACTION**
        [MultiReturn(new[] { "X", "Y", "Z" })]
        public Dictionary<string, double> GetCoordinates() { return new Dictionary<string, double> { { "X", X }, { "Y", Y }, { "Z", Z } }; }
        public Point GetPoint() { return Point.ByCoordinates(X, Y, Z); }
        public void AddEdge(Edge<Vertex> Edge) { if (Edge.Vertices.Contains(this)) Edges.Add(Edge); }
        public void AddEdges(IEnumerable<Edge<Vertex>> Edges) { for (int i = 0; i < Edges.Count(); i++) AddEdge(Edges.ElementAt(i)); }
        public void AddFace(Face<Vertex> Face) { if (Face.Vertices.Contains(this)) Faces.Add(Face); }
        public void AddFaces(IEnumerable<Face<Vertex>> Faces) { for (int i = 0; i < Faces.Count(); i++) AddFace(Faces.ElementAt(i)); }

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


    public class Edge<T> : HashSet<HalfEdge<T>>, IDisposable where T : Vertex
    {
        //**FIELDS**
        private bool disposed = false;

        //**PROPERTIES** //**QUERY**
        public string Name { get; set; }
        public double[] Angle { get; set; }

        public Point MidPoint { get { return Point.ByCoordinates(this.ToList()[0][0].X + this.ToList()[0][1].X, this.ToList()[0][0].Y, this.ToList()[0][0].Z); } }
        public Face<T>[] Faces { get { return new Face<T>[2] { this.ToList()[0].Face, this.ToList()[1].Face }; } }
        public T[] Vertices { get { return new T[2] { this.ToList()[0][0], this.ToList()[1][0] }; } }

        //**CONSTRUCTOR**
        internal Edge() : base() { Name = ""; Angle = null; }
        internal Edge(IEnumerable<HalfEdge<T>> HalfEdges) : base(HalfEdges) { Name = ""; Angle = null; }
        internal Edge(IEnumerable<HalfEdge<T>> HalfEdges, string Name) : base(HalfEdges) { this.Name = Name; Angle = null; }

        //**METHODS** //**ACTION**
        public T GetOtherVertex(T Vertex)
        {
            if (this.ElementAt(0)[0].Equals(Vertex)) return this.ElementAt(0)[1];
            if (this.ElementAt(1)[0].Equals(Vertex)) return this.ElementAt(1)[1];
            return null;
        }
        public Line GetLine()
        {
            Line output;
            using (Point a = this.ToList()[0][0].GetPoint())
            {
                using (Point b = this.ToList()[0][1].GetPoint())
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
        public void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                this.Clear();
                if (MidPoint != null) MidPoint.Dispose();
            }
            disposed = true;
        }
    }

    public class Face<T> : List<HalfEdge<T>>, IDisposable where T: Vertex
    {
        //**FIELDS**
        private bool disposed = false;

        //**PROPERTIES** //**QUERY**
        public string Name { get; set; }
        public CoordinateSystem CS { get; set; }
        public Point Center { get { return CS.Origin; } }
        public Vector Normal { get { return CS.ZAxis; } }
        public Dictionary<string, Object> Parameters { get; set; }

        public T[] Vertices
        {
            get
            {
                T[] output = new T[Count];
                for (int i = 0; i < Count; i++) output[i] = this[i][0];
                return output;
            }
        }
        public Edge<T>[] Edges
        {
            get
            {
                Edge<T>[] output = new Edge<T>[Count];
                for (int i = 0; i < Count; i++) output[i] = this[i].Edge;
                return output;
            }
        }

        //**CONSTRUCTOR**
        internal Face() : base() { }
        internal Face(IEnumerable<T> Vertices)
            : base(Vertices.Count())
        {
            for (int i = 0; i < Vertices.Count(); i++)
            {
                this.Add(new HalfEdge<T>(Vertices.ElementAt(i), Vertices.ElementAt((i + 1) % Capacity)));
                Vertices.ElementAt(i).AddFace(this as Face<Vertex>);
            }
            this.ForEach(he => he.AddFace(this));
            double[] xyz = { 0, 0, 0 };
            for (int i = 0; i < Count; i++)
            {
                xyz[0] += Vertices.ElementAt(i).X / Count;
                xyz[1] += Vertices.ElementAt(i).Y / Count;
                xyz[2] += Vertices.ElementAt(i).Z / Count;
            }
            Point Center = Point.ByCoordinates(xyz[0], xyz[1], xyz[2]);
            Vector X = this[Count - 1].GetVector();
            X = X.Reverse();
            Vector Y = this[0].GetVector();
            Vector Z = X.Cross(Y);
            Y = Z.Cross(X);
            CS = CoordinateSystem.ByOriginVectors(Center, X, Y);
            X.Dispose(); Y.Dispose(); Z.Dispose(); Center.Dispose();
        }

        //**METHODS** //**ACTION**
        public Face<T> ReOrderVertices(T Start)
        {
            if (this[0][0].Equals(Start)) return null;
            int index = 0;
            HalfEdge<T>[] temp = new HalfEdge<T>[this.Count];
            this.CopyTo(temp);
            for (int i = 1; i < Count; i++) if (this[i].Equals(Start)) { index = i; break; }
            for (int i = 0; i < Count; i++) this[i] = temp[(index + 1) % Count];
            temp = null;
            return this;
        }
        public void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                this.ForEach(e => e.Edge.Remove(e));
                this.Clear();
                if (Center != null) Center.Dispose();
                if (Normal != null) Normal.Dispose();
                if (CS != null) CS.Dispose();
                if (Parameters != null) for (int i = 0; i < Parameters.Count; i++)
                        if (Parameters.Values.ToArray()[i] is IDisposable) ((IDisposable)Parameters.Values.ToArray()[i]).Dispose();
            }
            disposed = true;
        }
    }

    public class Triangle : Face<Vertex>, IDisposable
    {
        //**FIELDS**
        private bool disposed = false;

        //**PROPERTIES** //**QUERY**
        public Point Circumcenter { get; private set; }

        //**CONSTRUCTOR**
        internal Triangle(IEnumerable<Vertex> Vertices)
            : base(Vertices)
        {
            Point[] pts = { Vertices.ElementAt(0).GetPoint(), Vertices.ElementAt(1).GetPoint(), Vertices.ElementAt(2).GetPoint() };
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
                Point[] Points = { Vertices[0].GetPoint(), Vertices[1].GetPoint(), Vertices[2].GetPoint() };
                return Points;
            }
        }
        public new void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
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

    public class Quad : Face<Vertex>
    {
        //**FIELDS**
        private bool disposed = false;

        public new void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
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
        private Point circumcenter;
        private bool disposed = false;
        private Dictionary<Point, Vertex> V;
        private Dictionary<string, Edge<Vertex>> E;

        //**PROPERTIES**QUERY
        public List<Face<Vertex>> Faces { get; set; }
        public List<Vertex> Vertices { get { return V.Values.ToList(); } }
        public List<Edge<Vertex>> Edges { get; set; }
        public Point[] Points { get; set; }

        internal Mesh(Surface[] Surfaces, Point[] Points)
        {
            // initialize
            Faces = new List<Face<Vertex>>(Surfaces.Length);
            V = new Dictionary<Point, Vertex>(Points.Length);
            E = new Dictionary<string, Edge<Vertex>>();
            Edges = new List<Edge<Vertex>>();

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
                for (int j = 0; j < t.Count; j++)
                {
                    bool edgeFound = false;
                    if (Edges.Count > 0)
                    {
                        for (int k = 0; k < Edges.Count; k++)
                        {
                            if (Edges[k].Vertices.Contains(t[j][0]) && Edges[k].Vertices.Contains(t[j][1]))
                            {
                                Edges[k].Add(t[j]);
                                t[j].AddEdge(Edges[k]);
                                edgeFound = true;
                                break;
                            }
                        }
                    }
                    if (!edgeFound)
                    {
                        Edge<Vertex> e = new Edge<Vertex>();
                        e.Add(t[j]);
                        t[j].AddEdge(e);
                        Edges.Add(e);
                    }
                }

            }
        }

        public static Mesh BySurfacesPoints(Surface[] Surfaces, Point[] Points) { return new Mesh(Surfaces, Points); }

        public bool AddEdgeNames(PolyCurve[] Spline, int Digits = 3)
        {
            int D = (Spline.Length + 1).ToString().Length;
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
            return false;
        }

        public void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                Edges.ForEach(e => e = null);
                Faces.ForEach(f => f.Dispose());
                Vertices.ForEach(v => v = null);
            }
            disposed = true;
        }
    }
}
