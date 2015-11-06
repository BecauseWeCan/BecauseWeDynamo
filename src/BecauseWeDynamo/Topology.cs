using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using Text;
using Autodesk.Dynamo.MeshToolkit;

namespace Topology
{

    public class HalfEdge : List<Vertex>
    {
        //**FIELDS**
        private Edge edge;
        private Face face;
        private bool disposed = false;

        //**PROPERTIES** //**QUERY**
        public Edge Edge { get { return edge; } }
        public Face Face { get { return face; } }

        //**CONSTRUCTOR**
        internal HalfEdge(Vertex A, Vertex B) : base(2) { this.Add(A); this.Add(B); }
        internal HalfEdge(Vertex A, Vertex B, Edge Edge, Face Face) : this(A, B) { edge = Edge; face = Face; }
        internal HalfEdge(IEnumerable<Vertex> Vertices) : base(Vertices) { }
        internal HalfEdge(IEnumerable<Vertex> Vertices, Edge Edge, Face Face) : this(Vertices) { edge = Edge; face = Face; }

        //**METHODS** //**ACTION**
        public Vector GetVector()
        {
            Point A = this[0].GetPoint();
            Point B = this[1].GetPoint();
            Vector output = Vector.ByTwoPoints(A, B);
            A.Dispose(); B.Dispose();
            return output;
        }
        public HalfEdge FlipDirection()
        {
            List<Vertex> temp = new List<Vertex>(this);
            this[0] = temp[1];
            this[1] = temp[0];
            temp = null;
            return this;
        }
        public void Dispose() { edge = null; }
    }

    public class Vertex
    {
        //**FIELDS**


        //**PROPERTIES** //**QUERY**
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public HashSet<Edge> Edges { get; set; }
        public HashSet<Face> Faces { get; set; }

        //**CONSTRUCTOR**
        internal Vertex(Point Point)
        {
            X = Point.X; Y = Point.Y; Z = Point.Z;
            Edges = new HashSet<Edge>();
            Faces = new HashSet<Face>();
        }
        internal Vertex(Point Point, IEnumerable<Edge> Edges, IEnumerable<Face> Faces)
        {
            X = Point.X; Y = Point.Y; Z = Point.Z;
            this.Edges = new HashSet<Edge>(Edges);
            this.Faces = new HashSet<Face>(Faces);
        }
        internal Vertex(double X, double Y, double Z)
        {
            this.X = X; this.Y = Y; this.Z = Z;
            Edges = new HashSet<Edge>();
            Faces = new HashSet<Face>();
        }
        internal Vertex(double X, double Y, double Z, IEnumerable<Edge> Edges, IEnumerable<Face> Faces)
        {
            this.X = X; this.Y = Y; this.Z = Z;
            this.Edges = new HashSet<Edge>(Edges);
            this.Faces = new HashSet<Face>(Faces);
        }

        //**METHODS** //**ACTION**
        [MultiReturn(new[] { "X", "Y", "Z" })]
        public Dictionary<string, double> GetCoordinates() { return new Dictionary<string, double> { { "X", X }, { "Y", Y }, { "Z", Z } }; }
        public Point GetPoint() { return Point.ByCoordinates(X, Y, Z); }
        public void AddEdge(Edge Edge) { if (Edge.Vertices.Contains(this)) Edges.Add(Edge); }
        public void AddEdges(IEnumerable<Edge> Edges) { for (int i = 0; i < Edges.Count(); i++) AddEdge(Edges.ElementAt(i)); }
        public void AddFace(Face Face) { if (Face.Vertices.Contains(this)) Faces.Add(Face); }
        public void AddFaces(IEnumerable<Face> Faces) { for (int i = 0; i < Faces.Count(); i++) AddFace(Faces.ElementAt(i)); }

        public override bool Equals(System.Object Object)
        {
            if (Object == null) return false;
            Vertex vtx = Object as Vertex;
            if ((System.Object)vtx == null) return false;
            return (X == vtx.X && Y == vtx.Y && Z == vtx.Z);
        }
        public bool Equals(Vertex Vertex) { return (X == Vertex.X && Y == Vertex.Y && Z == Vertex.Z); }
        public bool IsAtPoint(Point Point) { return (X == Point.X && Y == Point.Y && Z == Point.Z); }
        public override int GetHashCode() { return string.Format("{0}-{1}-{2}", X, Y, Z).GetHashCode(); }
    }


    public class Edge : HashSet<HalfEdge>
    {
        //**FIELDS**
        private bool disposed = false;

        //**PROPERTIES** //**QUERY**
        public string Name { get; set; }
        public double[] Angle { get; set; }
        public Point MidPoint { get { return Point.ByCoordinates(this.ToList()[0][0].X + this.ToList()[0][1].X, this.ToList()[0][0].Y, this.ToList()[0][0].Z); } }
        public Face[] Faces { get { return new Face[2] { this.ToList()[0].Face, this.ToList()[1].Face }; } }
        public Vertex[] Vertices { get { return new Vertex[2] { this.ToList()[0][0], this.ToList()[1][0] }; } }

        //**CONSTRUCTOR**
        internal Edge() : base() { }
        internal Edge(IEnumerable<HalfEdge> HalfEdges) : base(HalfEdges) { }
        internal Edge(IEnumerable<HalfEdge> HalfEdges, string Name) : this(HalfEdges) { this.Name = Name; }
        internal Edge(IEnumerable<Vertex> Vertices) : base() { for (int i = 0; i < 2; i++) this.Add(new HalfEdge(Vertices.ElementAt(i), Vertices.ElementAt((i + 1) % 2))); }
        internal Edge(IEnumerable<Vertex> Vertices, string Name) : this(Vertices) { this.Name = Name; }

        //**METHODS** //**ACTION**
        public HalfEdge GetOtherHalfEdge(HalfEdge HalfEdge)
        {
            if (this.ElementAt(0).Equals(HalfEdge)) return this.ElementAt(1);
            if (this.ElementAt(1).Equals(HalfEdge)) return this.ElementAt(0);
            return null;
        }
        public Vertex GetOtherVertex(Vertex Vertex)
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
    }

    public class Face : List<HalfEdge>, IDisposable
    {
        //**FIELDS**
        private bool disposed = false;

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
                Vertex[] output = new Vertex[Count];
                for (int i = 0; i < Count; i++) output[i] = this[i][0];
                return output;
            }
        }
        public Edge[] Edges
        {
            get
            {
                Edge[] output = new Edge[Count];
                for (int i = 0; i < Count; i++) output[i] = this[i].Edge;
                return output;
            }
        }

        //**CONSTRUCTOR**
        internal Face() : base() { }
        internal Face(IEnumerable<Vertex> Vertices)
            : base(Vertices.Count())
        {
            for (int i = 0; i < Vertices.Count(); i++)
            {
                this.Add(new HalfEdge(Vertices.ElementAt(i), Vertices.ElementAt((i + 1) % Capacity)));
            }
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
        public Face ReOrderVertices(Vertex Start)
        {
            if (this[0][0].Equals(Start)) return this;
            int index = 0;
            HalfEdge[] temp = new HalfEdge[this.Count];
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
                for (int i = 0; i < Count; i++)
                {
                    if (this[i].Edge.Count < 2) this[i].Dispose();
                    if (this[i][0].Faces.Count < 2) this[i][0] = null;
                    this[i] = null;
                }
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

    public class Triangle : Face
    {
        //**FIELDS**
        private Point circumcenter;
        private bool disposed = false;

        //**PROPERTIES** //**QUERY**
        public Point Circumcenter { get { return circumcenter; } }

        //**CONSTRUCTOR**
        internal Triangle(IEnumerable<Vertex> Vertices)
            : base(Vertices)
        {
            Point[] pts = { Vertices.ElementAt(0).GetPoint(), Vertices.ElementAt(1).GetPoint(), Vertices.ElementAt(2).GetPoint() };
            Circle c = Circle.ByBestFitThroughPoints(pts);
            circumcenter = c.CenterPoint;
            c.Dispose(); pts[0].Dispose(); pts[1].Dispose(); pts[2].Dispose(); pts = null;
        }

        //**METHODS** //**ACTION**
        public void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected virtual void Dispose(bool disposing)
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
        private Point circumcenter;
        private bool disposed = false;

        public void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected virtual void Dispose(bool disposing)
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

        //**PROPERTIES**QUERY
        public List<Face> Faces { get; set; }
        public List<Vertex> Vertices { get { return V.Values.ToList(); } }
        public List<Edge> Edges { get; set; }
        public Point[] Points { get; set; }

        internal Mesh(Surface[] Surfaces, Point[] Points)
        {
            Faces = new List<Face>(Surfaces.Length);
            V = new Dictionary<Point, Vertex>(Points.Length);
            Edges = new List<Edge>();
            this.Points = Points;
            for (int i = 0; i < Points.Length; i++) V.Add(Points[i], new Vertex(Points[i]));
            for (int i = 0; i < Surfaces.Length; i++)
            {
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
                            break ;
                        }
                    }
                    pt.Dispose();
                }
                vtx.ForEach(x => x.Dispose());

                Triangle t = new Triangle(v);

                
            }
        }

    }
}
