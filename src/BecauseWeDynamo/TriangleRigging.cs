using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using Text;

namespace TriangleRigging
{
    public class TriangleEdgeComparer : IEqualityComparer<TriangleEdge>
    {
        public bool Equals(TriangleEdge e1, TriangleEdge e2)
        {
            return ((e1.A.AlmostEquals(e2.A) && e1.B.AlmostEquals(e2.B)) || (e1.A.AlmostEquals(e2.B) && e1.B.AlmostEquals(e2.A)));
        }
        public int GetHashCode(TriangleEdge e)
        {
            return e.GetHashCode();
        }
    }

    public class TriangleEdgeConnector
    {
        double radius;
        double spacing;
        public int[] numHoles { get; set; }
        public TriangleEdge edge { get; set; }

        internal TriangleEdgeConnector(TriangleEdge tEdge, double radiusHole, double spacingHole)
        {
            spacing = spacingHole;
            radius = radiusHole;
            edge = tEdge;
            numHoles = new int[3] { 0, 0, 0 };
            for (int i = 0; i < edge.triangles.Count; i++)
            {
                numHoles[i] = (int)(Math.Ceiling((edge.triangles[i].Geometry.ClosestPointTo(edge.midpoint).DistanceTo(edge.midpoint) + 3 * radius) / spacing)) + edge.triangles[i].numHoles - 1;
            }
            numHoles[2] = numHoles[0] + numHoles[1] + 1;
        }

        public static TriangleEdgeConnector ByEdge(TriangleEdge edge, double radiusHole, double spacingHole) { return new TriangleEdgeConnector(edge, radiusHole, spacingHole); }

        public List<Circle> PlaceHoles(double DisplayFactor = 1)
        {
            if (numHoles[2] < 2 || edge.isOuterEdge) return null;
            List<Circle> result = new List<Circle>(numHoles[2]);

            for (int i = 0; i < edge.triangles.Count; i++)
            {
                for (int k = 0; k < numHoles[i]; k++)
                {
                    result.Add(Circle.ByCenterPointRadiusNormal(
                        (Point)edge.midpoint.Translate(Vector.ByTwoPoints(edge.midpoint, edge.triangles[i].Geometry.ClosestPointTo(edge.midpoint)), spacing * (k + 1)),
                        DisplayFactor * radius,
                        edge.triangles[i].Normal
                        ));
                }
            }

            result.Add(Circle.ByCenterPointRadiusNormal(
                            (Point)edge.midpoint,
                            DisplayFactor * radius,
                            edge.Normal
                            ));
            return result;
        }
        public List<Circle> PlaceHolesByHoleCount(int numHoles, double DisplayFactor = 1)
        {
            if (this.numHoles[2] == numHoles) return PlaceHoles(DisplayFactor);
            return null;
        }
    }

    public class TriangleSpline
    {
        //**PROPERTIES
        List<TriangleEdge> edges;
        Dictionary<Point, TriangleVertex> vtxs;
        List<Point> pts;
        Point[] stpts;

        //**QUERY**PROPERTIES
        public List<TriangleEdge> Edges { get { return edges; } }
        public List<Point> Points { get { return pts; } }

        //**CONSTRUCTOR
        internal TriangleSpline(Curve[] Curves, Point[] StartPoints)
        {
            edges = new List<TriangleEdge>(Curves.Length);
            vtxs = new Dictionary<Point, TriangleVertex>();
            stpts = StartPoints;

            for (int i = 0; i < Curves.Length; i++)
            {
                Point stpt = Curves[i].StartPoint;
                Point endpt = Curves[i].EndPoint;
                bool hasStart = false;
                bool hasEnd = false;
                for (int j = 0; j < vtxs.Keys.ToList().Count; j++)
                {
                    if (vtxs.Keys.ToList()[j].IsAlmostEqualTo(Curves[i].StartPoint))
                    {
                        stpt = vtxs.Keys.ToList()[j];
                        hasStart = true;
                    }
                    if (vtxs.Keys.ToList()[j].IsAlmostEqualTo(Curves[i].EndPoint))
                    {
                        endpt = vtxs.Keys.ToList()[j];
                        hasEnd = true;
                    }
                    if (hasEnd && hasStart) break;
                }
                if (!hasStart) vtxs.Add(stpt, new TriangleVertex(stpt));
                if (!hasEnd) vtxs.Add(endpt, new TriangleVertex(endpt));
                TriangleEdge e = new TriangleEdge(vtxs[stpt], vtxs[endpt], "e" + i.ToString("D3"), new List<Triangle>(), true);
                vtxs[stpt].AddEdge(e);
                vtxs[endpt].AddEdge(e);
                edges.Add(e);
            }
            for (int i = 0; i < StartPoints.Length; i++)
            {
                bool inList = false;
                for (int j = 0; j < vtxs.Keys.ToList().Count; j++)
                {
                    if (vtxs.Keys.ToList()[j].IsAlmostEqualTo(StartPoints[i]))
                    {
                        inList = true;
                        break;
                    }
                }
                if (!inList) vtxs.Add(StartPoints[i], new TriangleVertex(StartPoints[i]));

            }
            pts = vtxs.Keys.ToList();
        }

        //**CREATE
        public static TriangleSpline ByCurvesAndPoints(Curve[] Curves, Point[] StartPoints) { return new TriangleSpline(Curves, StartPoints); }

        //**ACTIONS*METHODS
        public List<List<Point>> GetSplines()
        {
            List<List<Point>> result = new List<List<Point>>(stpts.Length);
            for (int i = 0; i < stpts.Length; i++)
            {
                List<Point> spline = new List<Point>();
                for (int n = 0; n < pts.Count; n++)
                {
                    if (pts[n].IsAlmostEqualTo(stpts[i]))
                    {
                        spline.Add(pts[n]);
                        break;
                    }
                }
                if (vtxs[spline[0]].edges.Count == 0)
                {
                    result.Add(spline);
                    continue;
                }

                bool added = true;
                while (added)
                {
                    if (vtxs[spline[spline.Count - 1]].edges.Count < 1) break;
                    added = false;
                    for (int k = 0; k < vtxs[spline[spline.Count - 1]].edges.Count; k++)
                    {
                        if (!spline.Contains(vtxs[spline[spline.Count - 1]].edges[k].A.point))
                        {
                            spline.Add(vtxs[spline[spline.Count - 1]].edges[k].A.point);
                            added = true;
                            break;
                        }
                        if (!spline.Contains(vtxs[spline[spline.Count - 1]].edges[k].B.point))
                        {
                            spline.Add(vtxs[spline[spline.Count - 1]].edges[k].B.point);
                            added = true;
                            break;
                        }
                    }
                }
                result.Add(spline);
            }
            return result;
        }
    }

    public class TriangleVertex : IDisposable
    {
        private List<Triangle> Triangles;
        private List<TriangleEdge> Edges;
        private Point id;
        bool disposed = false;

        //**QUERY**PROPERTIES
        /// <summary>
        /// point of vertex
        /// </summary>
        public Point point { get { return id; } }
        public List<Triangle> triangles { get { return Triangles; } }
        public List<TriangleEdge> edges { get { return Edges; } }
        public int SplineId { get; set; }

        //*CONSTRUCTOR
        internal TriangleVertex(Point pt)
        {
            id = pt;
            Triangles = new List<Triangle>();
            Edges = new List<TriangleEdge>();
            SplineId = -1;
        }

        //**CREATE
        /// <summary>
        /// make vertex from point
        /// </summary>
        /// <param name="pt">point</param>
        /// <returns></returns>
        public static TriangleVertex ByPoint(Point pt) { return new TriangleVertex(pt); }

        //**ACTIONS**METHODS
        /// <summary>
        /// adds Triangle information to vertex
        /// </summary>
        /// <param name="triangle">Triangle object to be added</param>
        public void AddTriangle(Triangle triangle)
        {
            if (!Triangles.Contains(triangle)) Triangles.Add(triangle);
        }
        public void AddEdge(TriangleEdge edge)
        {
            if (!Edges.Contains(edge)) Edges.Add(edge);
        }
        public Boolean AlmostEquals(TriangleVertex vtx)
        {
            return (vtx.point.IsAlmostEqualTo(id));
        }
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing) id.Dispose();
            disposed = true;
        }
    }

    public class TriangleEdge : IDisposable
    {

        IEqualityComparer<TriangleEdge> almostequal;
        private TriangleVertex a;
        private TriangleVertex b;
        private Point Midpoint;
        private List<Triangle> Triangles;
        private string id;
        private Boolean IsOuterEdge;
        private List<TriangleVertex> Vertices;
        bool disposed = false;

        //**QUERY**PROPERTIES
        public Vector Normal { get; set; }
        public TriangleVertex A { get { return a; } }
        public TriangleVertex B { get { return b; } }
        public List<Triangle> triangles { get { return Triangles; } }
        public Point midpoint { get { return Midpoint; } }
        public string name
        {
            get { return id; }
            set { id = value; }
        }
        public Boolean isOuterEdge
        {
            get { return IsOuterEdge; }
            set { IsOuterEdge = value; }
        }

        //*CONSTRUCTOR
        internal TriangleEdge(TriangleVertex pta, TriangleVertex ptb, string name, List<Triangle> triangles, Boolean isOuterEdge)
        {
            Vertices = new List<TriangleVertex>(2);
            a = pta;
            b = ptb;
            Midpoint = Point.ByCoordinates((pta.point.X + ptb.point.X) / 2, (pta.point.Y + ptb.point.Y) / 2, (pta.point.Z + ptb.point.Z) / 2);
            id = name;
            Triangles = triangles;
            IsOuterEdge = isOuterEdge;
            Vertices.Add(pta);
            Vertices.Add(ptb);
            Normal = null;
        }

        //**CREATE
        public static TriangleEdge ByPoints(TriangleVertex pta, TriangleVertex ptb, string name, List<Triangle> triangles) { return new TriangleEdge(pta, ptb, name, triangles, false); }
        public static TriangleEdge ByPoints(TriangleVertex pta, TriangleVertex ptb, string name, List<Triangle> triangles, Boolean boo) { return new TriangleEdge(pta, ptb, name, triangles, boo); }

        //**ACTIONS**METHODS
        /// <summary>
        /// adds triangle reference data to internal list
        /// </summary>
        /// <param name="t"></param>
        public void AddTriangle(Triangle t) { if (!Triangles.Contains(t)) Triangles.Add(t); }
        /// <summary>
        /// compare this edge to another
        /// </summary>
        /// <param name="e">edge</param>
        /// <returns></returns>
        public Boolean AlmostEquals(TriangleEdge e) { return ((e.A.AlmostEquals(a) && e.B.AlmostEquals(b)) || (e.A.AlmostEquals(b) && e.B.AlmostEquals(A))); }
        /// <summary>
        /// sees if input vertices are edge endpoints
        /// </summary>
        /// <param name="vtx1">first vertex</param>
        /// <param name="vtx2">second vertex</param>
        /// <returns></returns>
        public Boolean HasVertices(TriangleVertex vtx1, TriangleVertex vtx2) { return ((vtx1.AlmostEquals(a) && vtx2.AlmostEquals(b)) || (vtx1.AlmostEquals(b) && vtx2.AlmostEquals(a))); }
        public Line GetEdgeGeometry() { return Line.ByStartPointEndPoint(A.point, B.point); }
        public List<PolyCurve> GetEdgeLabels(double factor)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int i = 0; i < triangles.Count; i++)
            {
                for (int j = 0; j < triangles[i].edges.Count; j++)
                {
                    if (triangles[i].edges[j].name.Equals(this.name))
                    {
                        int a = triangles[i].vertices.IndexOf(triangles[i].edges[j].A);
                        int b = triangles[i].vertices.IndexOf(triangles[i].edges[j].B);
                        if (((a - b) % 3 + 3) % 3 == 1)
                        {
                            labels.AddRange(Word.ByStringOriginVectors(
                                    triangles[i].edges[j].name,
                                    triangles[i].Geometry.ClosestPointTo(triangles[i].edges[j].midpoint),
                                    Vector.ByTwoPoints(triangles[i].edges[j].midpoint, triangles[i].edges[j].B.point),
                                    Vector.ByTwoPoints(triangles[i].center, triangles[i].edges[j].midpoint)
                                    ).display(factor));
                        }
                        else if (((a - b) % 3 + 3) % 3 == 2)
                        {
                            labels.AddRange(Word.ByStringOriginVectors(
                                    triangles[i].edges[j].name, //string
                                    triangles[i].Geometry.ClosestPointTo(triangles[i].edges[j].midpoint), //cs point
                                    Vector.ByTwoPoints(triangles[i].edges[j].midpoint, triangles[i].edges[j].A.point), // X-axis
                                    Vector.ByTwoPoints(triangles[i].center, triangles[i].edges[j].midpoint) // Y-axis
                                    ).display(factor));
                        } // end switch
                    }
                } // end Edges
            }

            return labels;
        }
        public TriangleVertex GetOtherVertex(TriangleVertex vtx)
        {
            if (vtx.AlmostEquals(A)) return B;
            else if (vtx.AlmostEquals(B)) return A;
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
                midpoint.Dispose();
                Midpoint.Dispose();
                Normal.Dispose();
            }
            disposed = true;
        }

    }

    public class Triangle : IDisposable
    {
        //*PRIVATE**PROPERTIES
        public Dictionary<string, Object> Parameters;
        public Boolean HasName;

        private double Height;
        private List<Point> Points;
        private List<TriangleVertex> Vertices;
        private List<TriangleEdge> Edges;
        private Point Circumcenter;
        private Point Center;
        private CoordinateSystem CS;

        bool disposed = false;

        //**PROPERTIES - **QUERY**
        public CoordinateSystem contextCoordinateSystem { get { return CS; } }
        public List<TriangleEdge> edges { get { return Edges; } }
        public List<Point> points { get { return Points; } }
        public List<TriangleVertex> vertices { get { return Vertices; } }
        public Point circumcenter { get { return Circumcenter; } }
        public Point center { get { return Center; } }
        public Vector Normal { get; set; }
        public Geometry Geometry { get; set; }
        public string Name { get; set; }
        public int splineId { get; set; }
        public int numHoles { get; set; }
        public Curve[] PerimeterCurves
        {
            get
            {
                if (Geometry is Solid)
                {
                    Geometry[] g = Geometry.Explode();
                    Curve[] result = ((Surface)Geometry.Explode()[0]).PerimeterCurves();
                    g.ForEach(i => i.Dispose());
                    return result;
                }
                else if (Geometry is Surface) return ((Surface)Geometry).PerimeterCurves();
                else return null;
            }
        }
        public Boolean HasEdges
        {
            get
            {
                if (Edges.Count == 3) return true;
                else return false;
            }
        }

        //**CONSTRUCTOR
        internal Triangle(Geometry inputGeometry, List<Point> points)
        {
            Geometry = inputGeometry;
            Points = new List<Point>(3);
            Edges = new List<TriangleEdge>(3);
            Vertices = new List<TriangleVertex>(3);
            Name = "";
            HasName = false;
            Height = 0;
            numHoles = 1;

            Points.AddRange(points);
            Normal = Plane.ByBestFitThroughPoints(points).Normal;
            Circumcenter = Circle.ByThreePoints(points[0], points[1], points[2]).CenterPoint;
            Center = Point.ByCoordinates((points[0].X + points[1].X + points[2].X) / 3.0, (points[0].Y + points[1].Y + points[2].Y) / 3.0, (points[0].Z + points[1].Z + points[2].Z) / 3.0);
            CS = CoordinateSystem.ByOriginVectors(Center, Vector.ByTwoPoints(points[0], points[1]), Vector.ByTwoPoints(points[0], points[2]));
        }

        //**CREATE**
        /// <summary>
        /// 
        /// </summary>
        /// <param name="srf"></param>
        /// <param name="points"></param>
        /// <returns></returns>
        public static Triangle BySurfaceAndPoints(Surface srf, List<Point> points) { return new Triangle(srf, points); }

        //**ACTIONS**METHODS
        /// <summary>
        /// build internal geometry
        /// </summary>
        /// <param name="edge">base edge</param>
        public void BuildGeometricProperties(TriangleEdge edge)
        {
            int a = vertices.IndexOf(edge.A);
            int b = vertices.IndexOf(edge.B);
            int c = 3 - (a + b);
            if (((a - b) % 3 + 3) % 3 == 1) CS = CoordinateSystem.ByOriginVectors(center, Vector.ByTwoPoints(Points[b], Points[a]), Vector.ByTwoPoints(Points[b], Points[c]));
            else if (((a - b) % 3 + 3) % 3 == 2) CS = CoordinateSystem.ByOriginVectors(center, Vector.ByTwoPoints(Points[a], Points[b]), Vector.ByTwoPoints(Points[a], Points[c]));
            double w = points[a].DistanceTo(points[b]);
            double sa = points[c].DistanceTo(points[a]);
            double sb = points[c].DistanceTo(points[b]);
            double s = (w + sa + sb) / 2;
            Height = 2 * Math.Sqrt(s * (s - w) * (s - sa) * (s - sb)) / w;
        }
        /// <summary>
        /// adjusts CS
        /// </summary>
        public void BuildCS()
        {
            for (int i = 0; i < edges.Count; i++)
            {
                if (edges[i].name.Substring(1, 2).Equals(Name.Substring(3, 2)))
                {
                    CS = CS.Rotate(CS.Origin, CS.ZAxis, 180);
                }
            }
        }
        /// <summary>
        /// add edge reference to triangle
        /// </summary>
        /// <param name="edge">triangle edge</param>
        public void AddEdge(TriangleEdge edge) { if (!Edges.Contains(edge)) Edges.Add(edge); }
        /// <summary>
        /// add vertex reference to triangle
        /// </summary>
        /// <param name="vtx">triangle vertex</param>
        public void AddVertex(TriangleVertex vtx) { if (!Vertices.Contains(vtx)) Vertices.Add(vtx); }
        /// <summary>
        /// adds circles to Edges
        /// </summary>
        /// <param name="factor">text scale</param>
        /// <returns></returns>
        public List<PolyCurve> GetEdgeLabels(double factor)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < edges.Count; j++)
            {
                if (edges[j].isOuterEdge) continue;

                int a = vertices.IndexOf(edges[j].A);
                int b = vertices.IndexOf(edges[j].B);

                if (((a - b) % 3 + 3) % 3 == 1)
                {
                    labels.AddRange(Word.ByStringOriginVectors(
                                edges[j].name,
                                Geometry.ClosestPointTo(edges[j].midpoint),
                                Vector.ByTwoPoints(edges[j].midpoint, edges[j].B.point),
                                Vector.ByTwoPoints(center, edges[j].midpoint)
                                ).display(factor));
                }
                else if (((a - b) % 3 + 3) % 3 == 2)
                {
                    labels.AddRange(Word.ByStringOriginVectors(
                                edges[j].name, //string
                                Geometry.ClosestPointTo(edges[j].midpoint), //cs point
                                Vector.ByTwoPoints(edges[j].midpoint, edges[j].A.point), // X-axis
                                Vector.ByTwoPoints(center, edges[j].midpoint) // Y-axis
                                ).display(factor));
                } // end switch
            } // end Edges

            return labels;
        }// end action
        /// <summary>
        /// adds circles to Edges
        /// </summary>
        /// <param name="factor">text scale</param>
        /// <returns></returns>
        public List<Circle> PlaceHoles(int holes, double radius, double spacing)
        {
            List<Circle> circles = new List<Circle>();
            for (int j = 0; j < edges.Count; j++)
            {
                if (edges[j].isOuterEdge) continue;

                int a = vertices.IndexOf(edges[j].A);
                int b = vertices.IndexOf(edges[j].B);

                if (((a - b) % 3 + 3) % 3 == 1)
                {
                    double offset = spacing * Math.Ceiling((Geometry.ClosestPointTo(edges[j].midpoint).DistanceTo(edges[j].midpoint) + 3 * radius) / spacing);
                    for (int k = 0; k < holes; k++)
                    {
                        circles.Add(Circle.ByCenterPointRadiusNormal(
                            (Point)edges[j].midpoint.Translate(Vector.ByTwoPoints(edges[j].midpoint, Geometry.ClosestPointTo(edges[j].midpoint)), offset),
                            radius,
                            Normal
                            ));
                        offset += spacing;
                    }
                }
                else if (((a - b) % 3 + 3) % 3 == 2)
                {
                    double offset = spacing * Math.Ceiling((Geometry.ClosestPointTo(edges[j].midpoint).DistanceTo(edges[j].midpoint) + 3 * radius) / spacing);
                    for (int k = 0; k < holes; k++)
                    {
                        circles.Add(Circle.ByCenterPointRadiusNormal(
                            (Point)edges[j].midpoint.Translate(Vector.ByTwoPoints(edges[j].midpoint, Geometry.ClosestPointTo(edges[j].midpoint)), offset),
                            radius,
                            Normal
                            ));
                        offset += spacing;
                    }
                } // end switch
            } // end Edges
            return circles;
        }// end action
        /// <summary>
        /// get profiles and letters
        /// </summary>
        /// <param name="factor">letter size</param>
        /// <returns>curves</returns>
        public List<Curve> GetTriangleCurves(double factor)
        {
            List<Curve> curves = new List<Curve>();
            List<PolyCurve> letters = GetEdgeLabels(factor);
            curves.AddRange(letters);
            curves.AddRange(PerimeterCurves);
            return curves;
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
                circumcenter.Dispose();
                Circumcenter.Dispose();
                center.Dispose();
                Center.Dispose();
                contextCoordinateSystem.Dispose();
                CS.Dispose();
                Normal.Dispose();
                Geometry.Dispose();
                for (int i = 0; i < edges.Count; i++)
                {
                    if (edges[i].isOuterEdge) edges[i].Dispose();
                    else edges[i].isOuterEdge = true;
                }
                for (int i = 0; i < vertices.Count; i++)
                {
                    bool dispose = true;
                    for (int j = 0; j < vertices[i].triangles.Count; j++)
                    {
                        if (!vertices[i].triangles[j].Name.Equals(Name)) dispose = dispose && vertices[i].triangles[j].disposed;
                    }
                    if (dispose) vertices[i].Dispose();
                }
                PerimeterCurves.ForEach(c => c.Dispose());
            }
            disposed = true;
        }
    }

    public class TriangleRigging : IDisposable
    {
        //**GLOBAL**VARIABLES
        private Dictionary<Point, TriangleVertex> Vertices;
        bool disposed = false;
        //**PROPERTIES - **QUERY**
        public List<Point> VertexPoints { get; set; }
        public List<Triangle> Triangles { get; set; }
        public List<List<TriangleVertex>> Splines { get; set; }
        public List<TriangleEdge> Edges { get; set; }
        public List<TriangleVertex> vertices { get { return Vertices.Values.ToList(); } }

        //**CONSTRUCTORS
        internal TriangleRigging(Solid[][] Solids, Point[][] SolidsPoint, Point[][] OrderedSplines)
        {
            // initialize Vertices and Triangles
            Vertices = new Dictionary<Point, TriangleVertex>();
            Triangles = new List<Triangle>(Solids.Length);
            Edges = new List<TriangleEdge>();
            Splines = new List<List<TriangleVertex>>(OrderedSplines.Length);

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            for (int i = 0; i < Solids.Length; i++)
            {
                Triangles.Add(new Triangle(Solids[i][0], SolidsPoint[i].ToList()));
                for (int j = 0; j < SolidsPoint[i].Length; j++)
                {
                    if (!Vertices.ContainsKey(SolidsPoint[i][j]))
                    {
                        Vertices.Add(SolidsPoint[i][j], new TriangleVertex(SolidsPoint[i][j]));
                    }
                    Vertices[SolidsPoint[i][j]].AddTriangle(Triangles[i]);
                    Triangles[i].AddVertex(Vertices[SolidsPoint[i][j]]);
                }
            }

            VertexPoints = Vertices.Keys.ToList();

            BuildEdges(OrderedSplines);
        }
        internal TriangleRigging(Surface[] TriangleSurfaces, Point[][] OrderedSplines)
        {
            // initialize Vertices and Triangles
            Vertices = new Dictionary<Point, TriangleVertex>();
            Triangles = new List<Triangle>(TriangleSurfaces.Length);
            Edges = new List<TriangleEdge>();
            Splines = new List<List<TriangleVertex>>(OrderedSplines.Length);

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            for (int i = 0; i < TriangleSurfaces.Length; i++)
            {
                List<Point> pts = new List<Point>(3);
                TriangleSurfaces[i].Vertices.ForEach(vtx => pts.Add(vtx.PointGeometry));
                Triangles.Add(new Triangle(TriangleSurfaces[i], pts));
                for (int j = 0; j < pts.Count; j++)
                {
                    if (!Vertices.ContainsKey(pts[j]))
                    {
                        Vertices.Add(pts[j], new TriangleVertex(pts[j]));
                    }
                    Vertices[pts[j]].AddTriangle(Triangles[i]);
                    Triangles[i].AddVertex(Vertices[pts[j]]);
                }
            }

            VertexPoints = Vertices.Keys.ToList();

            BuildEdges(OrderedSplines);
        }
        internal TriangleRigging(Solid[][] Solids, Point[][] SolidsPoint, Point start)
        {
            // initialize Vertices and Triangles
            Vertices = new Dictionary<Point, TriangleVertex>();
            Triangles = new List<Triangle>(Solids.Length);
            Edges = new List<TriangleEdge>();
            Splines = new List<List<TriangleVertex>>();

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            for (int i = 0; i < Solids.Length; i++)
            {
                Triangle t = new Triangle((Surface)Solids[i][0].Explode()[0], SolidsPoint[i].ToList());
                Triangles.Add(t);
                for (int j = 0; j < SolidsPoint[i].Length; j++)
                {
                    if (!Vertices.ContainsKey(SolidsPoint[i][j]))
                    {
                        Vertices.Add(SolidsPoint[i][j], new TriangleVertex(SolidsPoint[i][j]));
                    }
                    Vertices[SolidsPoint[i][j]].AddTriangle(Triangles[i]);
                    Triangles[i].AddVertex(Vertices[SolidsPoint[i][j]]);
                    if (SolidsPoint[i][j].IsAlmostEqualTo(start)) Splines.Add(new List<TriangleVertex>(1) { Vertices[SolidsPoint[i][j]] });
                }
                for (int j = 0; j < SolidsPoint[i].Length; j++)
                {
                    TriangleEdge e = TriangleEdge.ByPoints(Vertices[SolidsPoint[i][j]], Vertices[SolidsPoint[i][(j + 1) % SolidsPoint[i].Length]], "e", new List<Triangle>(2) { t });
                    if (!Edges.Contains(e, new TriangleEdgeComparer()))
                    {
                        e.isOuterEdge = true;
                        Edges.Add(e);
                    }
                    else
                    {
                        TriangleEdge e0 = Edges.Find(f => f.AlmostEquals(e));
                        if (e0.isOuterEdge) e0.isOuterEdge = false;
                        e0.triangles.Add(t);
                        e.Dispose();
                    }
                }
            }

            VertexPoints = Vertices.Keys.ToList();

            //build unordered spline
            List<TriangleVertex> found = new List<TriangleVertex>() { Splines[0][0] };
            List<TriangleVertex> search = new List<TriangleVertex>(vertices);
            search.Remove(found[0]);
            int numSpline = 0;
            while (found.Count < vertices.Count)
            {
                List<TriangleVertex> spline = new List<TriangleVertex>();
                for (int i = 0; i < Splines[numSpline].Count; i++)
                    for (int j = 0; j < Splines[numSpline][i].edges.Count; j++)
                    {
                        bool AddVertex = false;
                        TriangleVertex vtx = Splines[numSpline][i].edges[j].GetOtherVertex(Splines[numSpline][i]);
                        if (found.Count < search.Count) { if (!found.Contains(vtx)) AddVertex = true; }
                        else { if (search.Contains(vtx)) AddVertex = true; }
                        if (AddVertex)
                        {
                            spline.Add(vtx);
                            found.Add(vtx);
                            search.Remove(vtx);
                        }
                    }

                Splines.Add(spline);
                numSpline++;
            }

            //order splines
            numSpline = Splines.Count.ToString().Length;
            List<List<TriangleVertex>> OrderedSplines = new List<List<TriangleVertex>>(Splines.Count);
            OrderedSplines.Add(Splines[0]);
            for (int i = 1; i < Splines.Count; i++)
            {

            }

        }


        //**PRIVATE**METHODS
        internal void BuildEdges(Point[][] OrderedSplines)
        {
            //build splines
            int numD = OrderedSplines.Length.ToString().Length;
            for (int i = 0; i < OrderedSplines.Length; i++)
            {
                List<TriangleVertex> Spline = new List<TriangleVertex>();
                int z = 0;
                int edgeCount = 1;
                int numDs = OrderedSplines[i].Length.ToString().Length;
                // iterate through spline vertices
                for (int j = 0; j < OrderedSplines[i].Length; j++)
                {
                    // search vertex points for match
                    for (int k = 0; k < VertexPoints.Count; k++)
                    {
                        if (VertexPoints[k].IsAlmostEqualTo(OrderedSplines[i][j]))
                        {
                            Spline.Add(Vertices[VertexPoints[k]]);
                            Vertices[VertexPoints[k]].SplineId = i;

                            if (j > 0)
                            {
                                // check Triangles
                                // find Triangles on edge
                                List<Triangle> listT = new List<Triangle>();
                                for (int t = 0; t < Vertices[VertexPoints[k]].triangles.Count; t++)
                                {
                                    if (Vertices[VertexPoints[z]].triangles.Contains(Vertices[VertexPoints[k]].triangles[t]))
                                        listT.Add(Vertices[VertexPoints[k]].triangles[t]);
                                }

                                if (listT.Count > 0)
                                {
                                    string edgeId = "s" + i.ToString("D" + numD) + "-" + edgeCount.ToString("D" + numDs);
                                    TriangleEdge e = new TriangleEdge(Vertices[VertexPoints[z]], Vertices[VertexPoints[k]], edgeId, listT, false);
                                    if (listT.Count == 1) e.isOuterEdge = true;
                                    Edges.Add(e);
                                    // add edge key to vertices
                                    for (int t = 0; t < listT.Count; t++)
                                    {
                                        listT[t].AddEdge(e); listT[t].BuildGeometricProperties(e);
                                        e.Normal = e.Normal.Add(listT[t].Normal);
                                    }
                                    Vertices[VertexPoints[z]].AddEdge(e);
                                    Vertices[VertexPoints[k]].AddEdge(e);
                                    edgeCount++;
                                }
                            }
                            z = k;
                            break;
                        }
                    } // end search
                } // end spline point iteration
                // check Triangles
                // find Triangles on edge
                if (OrderedSplines[i].Length > 2)
                {
                    List<Triangle> listTlast = new List<Triangle>();
                    for (int t = 0; t < Spline[0].triangles.Count; t++)
                    {
                        if (Vertices[VertexPoints[z]].triangles.Contains(Spline[0].triangles[t]))
                            listTlast.Add(Spline[0].triangles[t]);
                    }

                    if (listTlast.Count > 0)
                    {
                        string edgeId = "s" + i.ToString("D" + numD) + "-" + edgeCount.ToString("D" + numDs);
                        TriangleEdge e = new TriangleEdge(Vertices[VertexPoints[z]], Spline[0], edgeId, listTlast, false);
                        if (listTlast.Count == 1) e.isOuterEdge = true;
                        Edges.Add(e);
                        // add edge key to vertices
                        for (int t = 0; t < listTlast.Count; t++)
                        {
                            listTlast[t].AddEdge(e);
                            e.Normal = e.Normal.Add(listTlast[t].Normal);
                        }
                        Vertices[VertexPoints[z]].AddEdge(e);
                        Spline[0].AddEdge(e);
                        edgeCount++;
                    }
                }
                Splines.Add(Spline);
            } // end build spline


            //build Edges
            for (int s = 0; s < Splines.Count - 1; s++)
            {
                // for each spline
                int edgeCount = 0;
                int triangleCount = 1;
                int numDe = Splines[s].Count.ToString().Length;
                for (int v = 0; v < Splines[s].Count; v++)
                {
                    List<TriangleVertex> Search = new List<TriangleVertex>();
                    for (int t = 0; t < Splines[s][v].triangles.Count; t++)
                    {
                        if (!Splines[s][v].triangles[t].HasEdges)
                        {
                            for (int vt = 0; vt < Splines[s][v].triangles[t].vertices.Count; vt++)
                            {
                                if (!Search.Contains(Splines[s][v].triangles[t].vertices[vt]) && Splines[s][v].triangles[t].vertices[vt].SplineId == s + 1)
                                    Search.Add(Splines[s][v].triangles[t].vertices[vt]);
                            }
                        }
                    }
                    Search.Remove(Splines[s][v]);

                    for (int sv = 0; sv < Search.Count; sv++)
                    {
                        // check Triangles
                        // find Triangles on edge
                        List<Triangle> listT = new List<Triangle>();
                        for (int t = 0; t < Search[sv].triangles.Count; t++)
                        {
                            if (Splines[s][v].triangles.Contains(Search[sv].triangles[t]))
                                listT.Add(Search[sv].triangles[t]);
                        }

                        if (listT.Count > 0)
                        {
                            string edgeId = "e" + s.ToString("D" + numD) + (s + 1).ToString("D" + numD) + "-" + edgeCount.ToString("D" + numDe);
                            TriangleEdge e = new TriangleEdge(Splines[s][v], Search[sv], edgeId, listT, false);
                            if (listT.Count == 1) e.isOuterEdge = true;
                            Edges.Add(e);
                            // add edge key to vertices
                            for (int t = 0; t < listT.Count; t++)
                            {
                                listT[t].AddEdge(e);
                                e.Normal = e.Normal.Add(listT[t].Normal);
                                if (!listT[t].HasName)
                                {
                                    listT[t].Name = "t" + s.ToString("D" + numD) + (s + 1).ToString("D" + numD) + "-" + triangleCount.ToString("D" + numDe);
                                    listT[t].HasName = true;
                                    listT[t].splineId = s;
                                    listT[t].BuildCS();
                                    triangleCount++;
                                }
                            }
                            Splines[s][v].AddEdge(e);
                            Search[sv].AddEdge(e);
                            edgeCount++;
                        }
                    }


                } // end spline

            } // end build Edges 
        }

        //**STATIC**METHODS - **CREATE**
        public static TriangleRigging BySolidsPointsSplines(Solid[][] Solids, Point[][] SolidsPoints, Point[][] Splines) { return new TriangleRigging(Solids, SolidsPoints, Splines); }
        public static TriangleRigging BySolidsPointsPoints(Solid[][] Solids, Point[][] SolidsPoints, Point start) { return new TriangleRigging(Solids, SolidsPoints, start); }
        public static TriangleRigging ByTriangleSurfacesAndSplines(Surface[] Triangles, Point[][] Splines) { return new TriangleRigging(Triangles, Splines); }

        //**METHODS - **ACTION**
        public List<List<PolyCurve>> GetEdgeLabels(double factor)
        {
            List<List<PolyCurve>> labels = new List<List<PolyCurve>>(Triangles.Count);
            for (int i = 0; i < Triangles.Count; i++) { labels.Add(Triangles[i].GetEdgeLabels(factor)); }
            return labels;
        }
        /// <summary>
        /// gets index set of Triangles sorted by angles
        /// </summary>
        /// <returns></returns>
        public List<int> GetSortedTriangleIndexByName()
        {
            List<int> index = new List<int>(Triangles.Count);
            Dictionary<int, string> sorted = new Dictionary<int, string>();
            for (int i = 0; i < Triangles.Count; i++) { sorted.Add(i, Triangles[i].Name); }
            foreach (var item in sorted.OrderBy(i => i.Value)) { index.Add(item.Key); }

            return index;
        }

        public void SetHolesBySpline(int[] holes)
        {
            for (int t = 0; t < Triangles.Count; t++) Triangles[t].numHoles = holes[Triangles[t].splineId];
        }
        public List<TriangleEdgeConnector> PlaceEdgeConnectorHolesBySpline(int[] holes, double radius, double spacing)
        {
            SetHolesBySpline(holes);
            List<TriangleEdgeConnector> result = new List<TriangleEdgeConnector>();
            for (int i = 0; i < Edges.Count; i++)
            {
                if (Edges[i].isOuterEdge) continue;
                result.Add(new TriangleEdgeConnector(Edges[i], radius, spacing));
            }
            return result;
        }
        public List<TriangleEdgeConnector> GetEdgeConnectorByHoleCount(int HoleCount, List<TriangleEdgeConnector> input)
        {
            List<TriangleEdgeConnector> result = new List<TriangleEdgeConnector>();
            for (int i = 0; i < input.Count; i++) { if (input[i].numHoles[2] == HoleCount) result.Add(input[i]); }
            return result;
        }
        public void AddParamterToTriangles(string name, Object[] data)
        {
            if (data.Length == Triangles.Count) for (int i = 0; i < Triangles.Count; i++) { Triangles[i].Parameters.Add(name, data[i]); }
        }
        public List<int> GetTriangleIndexByParameterValue(string parameter, Object value)
        {
            if (!Triangles[0].Parameters.Keys.Contains(parameter)) return null;
            List<int> result = new List<int>();
            for (int i = 0; i < Triangles.Count; i++)
            {
                if (Triangles[i].Parameters[parameter].Equals(value)) result.Add(i);
            }
            return result;
        }

        public List<List<int>> GetTriangleRowIndexBySpline()
        {
            int numRow = Splines.Count - 1;
            int numD = Splines.Count.ToString().Length;
            List<int> Index = GetSortedTriangleIndexByName();
            List<int> Row = new List<int>(Index.Count);
            List<List<int>> result = new List<List<int>>(numRow);

            for (int i = 0; i < Index.Count; i++) Row.Add(Convert.ToInt32(Triangles[Index[i]].Name.Substring(1, numD)));
            for (int i = 0; i < numRow; i++)
            {
                List<int> list = new List<int>();
                for (int j = Row.IndexOf(i); j < Row.Count; j++)
                {
                    if (Row[j] > i) break;
                    list.Add(j);
                }
                result.Add(list);
            }
            return result;
        }

        public Mesh ConvertToMesh()
        {
            List<IndexGroup> indexGroup = new List<IndexGroup>(Triangles.Count);
            for (int t = 0; t < Triangles.Count; t++)
            {
                indexGroup.Add(IndexGroup.ByIndices(
                    (uint)VertexPoints.IndexOf(Triangles[t].points[0]),
                    (uint)VertexPoints.IndexOf(Triangles[t].points[1]),
                    (uint)VertexPoints.IndexOf(Triangles[t].points[2])
                    ));
            }
            return Mesh.ByPointsFaceIndices(VertexPoints, indexGroup);
        }
        public PolySurface ConvertToPolySurface()
        {
            List<Surface> TriangleSurfaces = new List<Surface>(Triangles.Count);
            for (int t = 0; t < Triangles.Count; t++) TriangleSurfaces.Add(Surface.ByPerimeterPoints(Triangles[t].points));
            PolySurface result = PolySurface.ByJoinedSurfaces(TriangleSurfaces);
            TriangleSurfaces.ForEach(s => s.Dispose());
            return result;
        }

        [MultiReturn(new[] { "VertexPoints", "Index" })]
        public static Dictionary<string, object> ConvertToMeshComponent(TriangleRigging triangleRigging)
        {
            List<int> indexGroup = new List<int>(triangleRigging.Triangles.Count);
            for (int t = 0; t < triangleRigging.Triangles.Count; t++) for (int p = 0; p < 3; p++) indexGroup.Add(triangleRigging.VertexPoints.IndexOf(triangleRigging.Triangles[t].points[p]));

            return new Dictionary<string, object>
            {
                { "VertexPoints",  triangleRigging.VertexPoints},
                { "Index", indexGroup }
            };
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
                Triangles.ForEach(t => t.Dispose());
                Splines.ForEach(s => s.ForEach(t => t.Dispose()));
                Edges.ForEach(e => e.Dispose());
                vertices.ForEach(v => v.Dispose());
                VertexPoints.ForEach(v => v.Dispose());
            }
            disposed = true;
        }



    } // end class


}