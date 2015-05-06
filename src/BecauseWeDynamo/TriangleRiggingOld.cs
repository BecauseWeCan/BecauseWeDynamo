using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;

/*
namespace TriangleRiggingOld
{
    public class TriangleVertex
    {
        private List<int> Triangles;
        private List<string> Edges;
        private Point id;
        private int Spline;
        private TriangleRigging Rig;

        public Point point {get { return id; } }
        public List<int> triangleIds { get { return Triangles; } }
        public List<string> edgeIds { get { return Edges; } }
        public int splineId { 
            get { return Spline; }
            set { Spline = value; }
        }

        internal TriangleVertex(Point pt, TriangleRigging rig)
        {
            id = pt;
            Triangles = new List<int>();
            Edges = new List<string>();
            Spline = -1;
            Rig = rig;
        }

        public static TriangleVertex ByPoint(Point pt, TriangleRigging rig)
        {
            return new TriangleVertex(pt, rig);
        }

        public void AddTriangle(int triangle)
        {
            if (!Triangles.Contains(triangle)) Triangles.Add(triangle);
        }

        public void AddEdge(string edge)
        {
            if (!Edges.Contains(edge)) Edges.Add(edge);
        }
    }

    public class TriangleEdge
    {
        private int a;
        private int b;
        private List<int> Triangles;
        private string id;
        private Boolean IsOuterEdge;

        public int A { get { return a; } }
        public int B { get { return b; } }
        public string name { 
            get { return id; }
            set { id = value; }
        }
        public Boolean isOuterEdge
        {
            get { return IsOuterEdge; }
            set { IsOuterEdge = value; }
        }
        public List<int> triangles { get { return Triangles; } }

        internal TriangleEdge(int pta, int ptb, string name, List<int> triangles, Boolean isOuterEdge)
        {
            a = pta;
            b = ptb;
            id = name;
            Triangles = triangles;
            IsOuterEdge = isOuterEdge;
        }

        public static TriangleEdge ByPoints(int pta, int ptb, string name, List<int> triangles)
        {
            return new TriangleEdge(pta, ptb, name, triangles, false);
        }

        public static TriangleEdge ByPoints(int pta, int ptb, string name, List<int> triangles, Boolean boo)
        {
            return new TriangleEdge(pta, ptb, name, triangles, boo);
        }

        public void AddTriangle(int t)
        {
            if (!Triangles.Contains(t)) Triangles.Add(t);
        }

    }

    public class Triangle
    {
        private List<Point> Vertices;
        private List<string> Edges;
        private Surface shape;
        public Boolean HasName;
        private string Name;

        public List<Point> vertices { get { return Vertices; }}
        public Surface surface { 
            get { return shape; }
            set { shape = value; }
        }
        public string name
        {
            get { return Name; }
            set { Name = value; }
        }

        internal Triangle(Surface srf, List<Point> points)
        {
            shape = srf;
            Vertices = new List<Point>();
            Edges = new List<string>();
            Name = "";
            HasName = false;
            Vertices.AddRange(points);
        }

        internal Triangle(List<Point> points)
        {
            Vertices = new List<Point>();
            Edges = new List<string>();
            Name = "";
            HasName = false;
            Vertices.AddRange(points);
        }

        public static Triangle BySurfaceAndPoints(Surface srf, List<Point> points)
        {
            return new Triangle(srf,points);
        }

        public static Triangle ByPoints(List<Point> points)
        {
            return new Triangle(points);
        }

        public void AddEdges(string edge)
        {
            if (!Edges.Contains(edge)) Edges.Add(edge);
        }

        public Boolean HasEdges()
        {
            if (Edges.Count == 3) return true;
            else return false;
        }
    }

    public class TriangleRigging
    {
        

        //**GLOBAL VARIABLES
        private Dictionary<Point, TriangleVertex> Vertices;
        private List<Point> VertexPoints;
        private List<Triangle> Triangles;
        private Dictionary<string, int> TriangleNames;
        private List<List<int>> Splines;
        private Dictionary<string, TriangleEdge> Edges;

        //**QUERY
        public List<Triangle> triangles { get { return Triangles; } }
        public List<string> triangleNames { get { return TriangleNames.Keys.ToList(); } }
        public List<TriangleEdge> edges { get { return Edges.Values.ToList(); } }
        public List<TriangleVertex> vertices{ get { return Vertices.Values.ToList(); } }
        public List<List<int>> splines { get { return Splines; } }


        //**CONSTRUCTORS
        internal TriangleRigging(Solid[][] Solids, Point[][] SolidsPoint, Point[][] OrderedSplines)
        {
            // initialize Vertices and Triangles
            Vertices = new Dictionary<Point, TriangleVertex>();
            Triangles = new List<Triangle>(Solids.Length);
            TriangleNames = new Dictionary<string, int>(Solids.Length);
            Edges = new Dictionary<string, TriangleEdge>();
            Splines = new List<List<int>>(OrderedSplines.Length);

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            for (int i = 0; i < Solids.Length; i++)
            {
                Triangles.Add(new Triangle((Surface) Solids[i][0].Explode()[0], SolidsPoint[i].ToList()));
                for (int j = 0; j < SolidsPoint[i].Length; j++)
                {
                    if (!Vertices.ContainsKey(SolidsPoint[i][j]))
                    {
                        Vertices.Add(SolidsPoint[i][j], new TriangleVertex(SolidsPoint[i][j], this));
                    }
                    Vertices[SolidsPoint[i][j]].AddTriangle(i);
                }
            }

            VertexPoints = Vertices.Keys.ToList();

            //build splines
            for (int i = 0; i < OrderedSplines.Length; i++)
            {
                List<int> Spline = new List<int>();
                int z = 0;
                int edgeCount = 1;
                // iterate through spline points
                for (int j = 0; j < OrderedSplines[i].Length; j++ )
                {
                    // search vertex points for match
                    for (int k = 0; k < VertexPoints.Count; k++)
                    {
                        if (VertexPoints[k].IsAlmostEqualTo(OrderedSplines[i][j]))
                        {
                            Spline.Add(k);
                            Vertices[VertexPoints[k]].splineId = i;

                            if (j > 0)
                            {
                                // check triangles
                                // find triangles on edge
                                List<int> triangles = new List<int>();
                                for (int t = 0; t < Vertices[VertexPoints[k]].triangleIds.Count; t++)
                                {
                                    if (Vertices[VertexPoints[z]].triangleIds.Contains(Vertices[VertexPoints[k]].triangleIds[t]))
                                        triangles.Add(Vertices[VertexPoints[k]].triangleIds[t]);
                                }

                                if (triangles.Count > 0)
                                {
                                    string edgeId = "s" + i + "-" + edgeCount;
                                    // add edge key to triangles
                                    for (int t = 0; t < triangles.Count; t++) { Triangles[triangles[t]].AddEdges(edgeId); }
                                    // create edges
                                    if (triangles.Count == 1) Edges.Add(edgeId, TriangleEdge.ByPoints(z, k, edgeId, triangles, true));
                                    else Edges.Add(edgeId, TriangleEdge.ByPoints(z, k, edgeId, triangles));
                                    // add edge key to vertices
                                    Vertices[VertexPoints[z]].AddEdge(edgeId);
                                    Vertices[VertexPoints[k]].AddEdge(edgeId);
                                    edgeCount++;
                                }
                            }
                                
                            z = k;
                            break;
                        }
                    } // end search
                } // end spline point iteration

                // find last edge if spline is circle
                // check triangles
                // find triangles on edge
                int v0 = Spline[0];
                z = Spline[Spline.Count - 1];
                List<int> ctriangles = new List<int>();
                for (int t = 0; t < Vertices[VertexPoints[v0]].triangleIds.Count; t++)
                {
                    if (Vertices[VertexPoints[z]].triangleIds.Contains(Vertices[VertexPoints[v0]].triangleIds[t]))
                        ctriangles.Add(Vertices[VertexPoints[v0]].triangleIds[t]);
                }

                if (ctriangles.Count > 0)
                {
                    string edgeId = "s" + i + "-" + edgeCount;
                    // add edge key to triangles
                    for (int t = 0; t < triangles.Count; t++) { Triangles[ctriangles[t]].AddEdges(edgeId); }
                    // create edges
                    if (triangles.Count == 1) Edges.Add(edgeId, TriangleEdge.ByPoints(z, v0, edgeId, ctriangles, true));
                    else Edges.Add(edgeId, TriangleEdge.ByPoints(z, v0, edgeId, ctriangles));
                    // add edge key to vertices
                    Vertices[VertexPoints[z]].AddEdge(edgeId);
                    Vertices[VertexPoints[v0]].AddEdge(edgeId);
                    edgeCount++;
                }


                Splines.Add(Spline);
            } // end build spline


            //build edges
            for (int i = 0; i < Splines.Count - 1; i++)
            {
                int j = i + 1;
                //spline point 
                int pti = 0;
                int ptj = 0;
                int num = 1;

                // build edges between splines
                // maximum number of edges is total number of points on splines - 1
                for (int e = 0; e < Splines[i].Count + Splines[j].Count - 1; e++)
                {
                    if (pti >= Splines[i].Count) break;
                    if (ptj >= Splines[j].Count) break;

                    int vi = Splines[i][pti];
                    int vj = Splines[j][ptj];

                    // check triangles
                    // find triangles on edge
                    List<int> triangles = new List<int>();
                    for (int t = 0; t < Vertices[VertexPoints[vj]].triangleIds.Count; t++)
                    {
                        if (Vertices[VertexPoints[vi]].triangleIds.Contains(Vertices[VertexPoints[vj]].triangleIds[t]))
                            triangles.Add(Vertices[VertexPoints[vj]].triangleIds[t]);
                    }

                    // create edges
                    string edgeId = "e" + i.ToString() + j.ToString() + "-" + e;
                    
                    // add edge key to triangles
                    for (int t = 0; t < triangles.Count; t++) 
                    { 
                        Triangles[triangles[t]].AddEdges(edgeId);
                        // tag triangle as complete if it has all edges, and it is last one, then add to both spline counters
                        if (Triangles[triangles[t]].HasEdges()) 
                        {
                            Triangles[triangles[t]].name = "t" + i.ToString() + j.ToString() + "-" + num;
                            TriangleNames.Add(Triangles[triangles[t]].name, triangles[t]);
                            Triangles[triangles[t]].HasName = true;
                            if (triangles.Count == 1)
                            {
                                pti++;
                                ptj++;
                            }
                            num++;
                        }
                        // otherwise figure out which spline nes point is on
                        else
                        {
                            for (int n = 0; n < Triangles[triangles[t]].vertices.Count; n++)
                            {
                                if (VertexPoints.IndexOf(Triangles[triangles[t]].vertices[n]) != vi 
                                    && VertexPoints.IndexOf(Triangles[triangles[t]].vertices[n]) != vj)
                                {
                                    if (Vertices[Triangles[triangles[t]].vertices[n]].splineId == i) { pti++; }
                                    else if (Vertices[Triangles[triangles[t]].vertices[n]].splineId == j) { ptj++; }
                                }
                            }
                        }
                    }

                    // create edges
                    if (triangles.Count == 1) Edges.Add(edgeId, TriangleEdge.ByPoints(vi, vj, edgeId, triangles, true));
                    else Edges.Add(edgeId, TriangleEdge.ByPoints(vi, vj, edgeId, triangles));
                    // add edge key to vertices
                    Vertices[VertexPoints[vi]].AddEdge(edgeId);
                    Vertices[VertexPoints[vj]].AddEdge(edgeId);

                } // end splines

            } // end build edges



        }

        //**CREATE
        public static TriangleRigging BySolidsPointsPoints(Solid[][] Solids, Point[][] SolidsPoints, Point[][] Splines)
        {
            return new TriangleRigging(Solids, SolidsPoints, Splines);
        }

        //**ACTION
        public TriangleVertex getVertex(Point pt)
        {
            if (VertexPoints.Contains(pt)) return Vertices[pt];
            return null;
        }
        public TriangleEdge getEdge(string name)
        {
            if (Edges.ContainsKey(name)) return Edges[name];
            return null;
        }
        public Triangle getTriangle(string name)
        {
            if (TriangleNames.ContainsKey(name)) return Triangles[TriangleNames[name]];
            return null;
        }

    }
}
*/