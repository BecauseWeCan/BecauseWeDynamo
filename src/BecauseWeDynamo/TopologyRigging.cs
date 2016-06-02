using System;
using System.Collections.Generic;
using System.Linq;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Runtime;
using Fabrication;
using Topology;
using Panelization;
using Geometry;

namespace Topology
{
    internal class rigging: mesh, IDisposable
    {
        //**PROPERTIES - **QUERY**
        new public List<List<vertex>> Splines { get; set; }
        public Dictionary<face, Dictionary<string, Object>> Parameters;
        private bool disposed = false;

        //**CONSTRUCTORS
        internal rigging(Solid[][] Solids, Point[][] SolidPoints, Point[][] OrderedSplines)
        {
            // initialize Vertices and Triangles
            V = new Dictionary<point, vertex>();
            Faces = new List<face>(Solids.Length);
            Edges = new List<edge>();
            Splines = new List<List<vertex>>(OrderedSplines.Length);
            Parameters = new Dictionary<face, Dictionary<string, object>>();

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            int eCount = 1;
            for (int i = 0; i < Solids.Length; i++)
            {
                List<vertex> TriangleV = new List<vertex>(3);
                for (int j = 0; j < SolidPoints[i].Length; j++)
                {
                    point v = point.ByPoint(SolidPoints[i][j]);
                    if (!V.ContainsKey(v))
                    {
                        V.Add(v, new vertex(v));
                    }
                    TriangleV.Add(V[v]);
                }
                triangle t = triangle.ByVertices(TriangleV);
                t.Name = "t" + (i + 1).ToString("D" + 4);
                Parameters.Add(t, new Dictionary<string,Object>());
                Parameters[t].Add("panel", Solids[i]);
                Faces.Add(t);
                eCount = FindFaceEdges(t, eCount);
            }
            Points = V.Keys.ToArray();
            BuildEdges(OrderedSplines);
        }
        internal rigging(Surface[] TriangleSurfaces, Point[][] OrderedSplines)
        {
            // initialize Vertices and Triangles
            V = new Dictionary<point, vertex>();
            Faces = new List<face>(TriangleSurfaces.Length);
            Edges = new List<edge>();
            Splines = new List<List<vertex>>(OrderedSplines.Length);

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            for (int i = 0; i < TriangleSurfaces.Length; i++)
            {
                Autodesk.DesignScript.Geometry.Vertex[] vtx = TriangleSurfaces[i].Vertices;
                List<vertex> v = new List<vertex>(vtx.Length);
                for (int j = 0; j < vtx.Length; j++)
                {
                    Point pt = vtx[j].PointGeometry;
                    vertex search = GetVertexAtPoint(pt);
                    if (!search.Equals(null)) v.Add(search);
                    else
                    {
                        point p = point.ByPoint(pt);
                        V.Add(p, new vertex(p));
                        v.Add(GetVertexAtPoint(pt));
                    }
                    pt.Dispose();
                }
                vtx.ForEach(x => x.Dispose());
            }

            Points = V.Keys.ToArray();
            BuildEdges(OrderedSplines);
        }
        

        //**PRIVATE**METHODS
        internal void BuildEdges(Point[][] OrderedSplines)
        {
            //build splines
            int numD = OrderedSplines.Length.ToString().Length;
            for (int i = 0; i < OrderedSplines.Length; i++)
            {
                List<vertex> Spline = new List<vertex>();
                int edgeCount = 1;
                int numDs = OrderedSplines[i].Length.ToString().Length;
                // iterate through spline Vertices
                for (int j = 0; j < OrderedSplines[i].Length; j++)
                {
                    if (!GetVertexAtPoint(OrderedSplines[i][j]).Equals(null))
                        {
                            Spline.Add(GetVertexAtPoint(OrderedSplines[i][j]));
                            if (j > 0)
                            {
                                Line L = Line.ByStartPointEndPoint(OrderedSplines[i][j-1],OrderedSplines[i][j]);
                                edge e = GetEdgeAtLine(L);
                                if (!e.Equals(null)) 
                                {
                                    e.Name = "s" + i.ToString("D" + numD) + "-" + edgeCount.ToString("D" + numDs);                                 
                                    edgeCount++;
                                }
                            }
                        }
                } // end search
                Splines.Add(Spline);
            } // end build spline


            //build Edges
            /*for (int s = 0; s < Splines.Count - 1; s++)
            {
                // for each spline
                int edgeCount = 0;
                int triangleCount = 1;
                int numDe = Splines[s].Count.ToString().Length;
                for (int v = 0; v < Splines[s].Count; v++)
                {
                    for (int sv = 0; sv < Search.Count; sv++)
                    {
                        // check Triangles
                        // find Triangles on edge
                        List<Triangle> listT = new List<Triangle>();
                        for (int t = 0; t < Search[sv].Triangles.Count; t++)
                        {
                            if (Splines[s][v].Triangles.Contains(Search[sv].Triangles[t]))
                                listT.Add(Search[sv].Triangles[t]);
                        }

                        if (listT.Count > 0)
                        {
                            string edgeId = "e" + s.ToString("D" + numD) + (s + 1).ToString("D" + numD) + "-" + edgeCount.ToString("D" + numDe);
                            Edge e = new Edge(Splines[s][v], Search[sv], edgeId, listT);
                            Edges.Add(e);
                            // add edge key to Vertices
                            for (int t = 0; t < listT.Count; t++)
                            {
                                listT[t].AddEdge(e);
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

            }*/ // end build Edges 
        }

        //**STATIC**METHODS - **CREATE**
        public static rigging BySolidsPointsSplines(Solid[][] Solids, Point[][] SolidsPoints, Point[][] Splines) { return new rigging(Solids, SolidsPoints, Splines); }

        //**METHODS - **ACTION**
        /// <summary>
        /// gets index set of Triangles sorted by angles
        /// </summary>
        /// <returns></returns>
        public List<int> GetSortedTriangleIndexByName()
        {
            List<int> index = new List<int>(Faces.Count);
            Dictionary<int, string> sorted = new Dictionary<int, string>();
            for (int i = 0; i < Faces.Count; i++) { sorted.Add(i, Faces[i].Name); }
            foreach (var item in sorted.OrderBy(i => i.Value)) { index.Add(item.Key); }

            return index;
        }
        public void AddParamterToTriangles(string name, Object[] data)
        {
            if (data.Length == Faces.Count) for (int i = 0; i < Faces.Count; i++) { Parameters[Faces[i]].Add(name, data[i]); }
        }
        public List<int> GetTriangleIndexByParameterValue(string parameter, Object value)
        {
            if (!Parameters[Faces[0]].Keys.Contains(parameter)) return null;
            List<int> result = new List<int>();
            for (int i = 0; i < Faces.Count; i++)
            {
                if (Parameters[Faces[i]][parameter].Equals(value)) result.Add(i);
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

            for (int i = 0; i < Index.Count; i++) Row.Add(Convert.ToInt32(Faces[Index[i]].Name.Substring(1, numD)));
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
        public void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                for (int j=0; j<Faces.Count;j++)
                {
                    face f = Faces[j];
                    for (int i = 0; i < Parameters[f].Values.Count; i++)
                    {
                        if (Parameters[f].Values.ToArray()[i] is IDisposable) ((IDisposable)Parameters.Values.ToArray()[i]).Dispose();
                    }
                }
            }
            disposed = true;
        }
    } // end class


}