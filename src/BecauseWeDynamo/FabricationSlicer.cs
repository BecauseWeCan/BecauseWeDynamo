using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using System;
using System.Collections.Generic;

namespace Fabrication
{
    /// <summary>
    /// slicing object
    /// </summary>
    public class Slicer : IDisposable
    {
        bool disposed = false;

        //**QUERY
        /// <summary>
        /// slice thickness
        /// </summary>
        public double Thickness { get; set; }
        /// <summary>
        /// spacing between slices
        /// </summary>
        public double Spacing { get; set; }
        /// <summary>
        /// solid to be sliced
        /// </summary>
        public Solid Solid { get; set; }
        /// <summary>
        /// List of CutPlane Arrays
        /// </summary>
        public List<Plane[]> CutPlanes { get; set; }

        //**CONSTRUCTORS
        internal Slicer(Solid solid, Plane plane, double thickness, double spacing)
        {
            Solid = solid;
            Thickness = thickness;
            Spacing = Math.Max(spacing, thickness);
            CutPlanes = new List<Plane[]>();
            CutPlanes.Add(GenerateCutPlanes(plane));
        }
        /// <summary>
        /// Internal constructor for the class that has all geometric inputs to create sliced model
        /// </summary>
        /// <param name="solid">Solid: geometry that is to be parsed</param>
        /// <param name="plane">Plane: primary cutplane</param>
        /// <param name="line1">Line1: (optional) defines secondary cutplane</param>
        /// <param name="line2">Line2: (optional) defines tertiary cutplane</param>
        /// <param name="thickness">Thickness: the thickness of the slices, or the thickness of the material to be used for the assembly</param>
        /// <param name="spacing">Spacing: the distance between each slice</param>
        internal Slicer(Solid solid, Plane plane, Line line1, Line line2, double thickness, double spacing)
            : this(solid, plane, thickness, spacing)
        {
            if (line1 != null)
            {
                Point pt1 = line1.StartPoint.Add(plane.Normal);
                Plane p1 = Plane.ByLineAndPoint(line1, pt1);
                CutPlanes.Add(GenerateCutPlanes(p1));
                pt1.Dispose(); p1.Dispose();
            }
            if (line2 != null)
            {
                Point pt2 = line1.StartPoint.Add(plane.Normal);
                Plane p2 = Plane.ByLineAndPoint(line2, pt2);
                CutPlanes.Add(GenerateCutPlanes(p2));
                pt2.Dispose(); p2.Dispose();
            }
        }
        /// <summary>
        /// creates slices using a curve as the spine
        /// </summary>
        /// <param name="solid">Solid: geometry that is to be parsed</param>
        /// <param name="curve">Curve: defines the normal used to create cut planes perpendicular to parameter "plane".
        /// If curve is too short, it will be extended using built-in extend function</param>
        /// <param name="thickness">Thickness: the thickness of the slices, or the thickness of the material to be used for the assembly</param>
        /// <param name="spacing">Spacing: the distance between each slice</param>
        /// <param name="origin">Origin</param>
        /// <returns>A newly-constructed Slicer object</returns>
        private Slicer(Solid solid, Curve curve, double thickness, double spacing, double origin)
        {
            Solid = solid;
            Thickness = thickness;
            Spacing = spacing;
            if (curve.IsPlanar)
            {
                Plane p1 = curve.PlaneAtParameter();
                CutPlanes.Add(GenerateCutPlanes(p1));
                p1.Dispose();
            }
            else
            {
                Point[] pts = curve.PointsAtEqualSegmentLength();
                Plane p1 = Plane.ByBestFitThroughPoints(pts);
                Curve crv = curve.PullOntoPlane(p1);
                CutPlanes.Add(GenerateCutPlanes(p1));
                CutPlanes.Add(GenerateCutPlanes(crv, origin));
                pts.ForEach(p => p.Dispose()); p1.Dispose(); crv.Dispose();
            }
        }

        //**CREATE
        /// <summary>
        /// 
        /// </summary>
        /// <param name="solid"></param>
        /// <param name="plane"></param>
        /// <param name="thickness"></param>
        /// <param name="spacing"></param>
        /// <returns></returns>
        public static Slicer ByPlane(Solid solid, Plane plane, double thickness, double spacing)
        {
            if (!solid.DoesIntersect(plane)) throw new ArgumentException("plane does not intersection solid");
            return new Slicer(solid, plane, thickness, spacing);
        }
        /// <summary>
        /// Construct an instance of slicer via a static method.
        /// This makes vertical stacks of given solid
        /// with each stack being of specified thickness
        /// and with given plane as the initial cutplane
        /// </summary>
        /// <param name="solid">Solid: geometry that is to be parsed</param>
        /// <param name="plane">Plane: primary cutplane</param>
        /// <param name="line1">Line1: (optional) defines secondary cutplane</param>
        /// <param name="line2">Line2: (optional) defines tertiary cutplane</param>
        /// <param name="thickness">Thickness: the thickness of the slices, or the thickness of the material to be used for the assembly</param>
        /// <param name="spacing">Spacing: the distance between each slice</param>
        /// <returns>A newly-constructed Slicer object</returns>
        public static Slicer ByPlaneAndLines(Solid solid, Plane plane, Line line1, Line line2 = null, double thickness = 0.25, double spacing = 1)
        {
            if (!solid.DoesIntersect(plane)) throw new ArgumentException("plane does not intersection solid");
            if (line1.Direction.IsParallel(plane.Normal)) throw new ArgumentException("line is perpendicular to plane");
            if (line2.Direction.IsParallel(plane.Normal)) throw new ArgumentException("line is perpendicular to plane");
            return new Slicer(solid, plane, line1, line2, thickness, spacing);
        }
        /// <summary>
        /// Construct an instance of slicer via a static method.
        /// This makes vertical stacks of given solid
        /// with each stack being of specified thickness
        /// and with given plane as the initial cutplane
        /// </summary>
        /// <param name="solid">Solid: geometry that is to be parsed</param>
        /// <param name="curve">Curve: defines the normal used to create cut planes perpendicular to parameter "plane".
        /// If curve is too short, it will be extended using built-in extend function</param>
        /// <param name="thickness">Thickness: the thickness of the slices, or the thickness of the material to be used for the assembly</param>
        /// <param name="spacing">Spacing: the distance between each slice</param>
        /// <param name="origin">Curve Parameter (0 - 1)</param>
        /// <returns>A newly-constructed Slicer object</returns>
        internal static Slicer ByCurve(Solid solid, Curve curve, double thickness = 0.25, double spacing = 1, double origin = 0)
        {
            return new Slicer(solid, curve, thickness, spacing, origin);
        }

        //**ACTIONS
        // *Note: needs to return results with intersection and assembly figured out;
        // *Insert code from 123Dmake;
        // *Results have placeholders
        /// <summary>
        /// Returns results as circles, planar surfaces, and solid slices.
        /// slices are within input solid so it can be skinned.
        /// </summary>
        /// <param name="slicer">slicer object</param>
        /// <returns>circles, planar surfaces, and solid slices listed by layer then by object</returns>
        [MultiReturn(new[] { "Profiles", "Surfaces", "Slices" })]
        public static Dictionary<string, object> GetResultsPlaneInscribed(Slicer slicer)
        {
            Plane[][] Shifted = slicer.ShiftCutPlanes(slicer.CutPlanes[0]);

            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();

            for (int i = 0; i < slicer.CutPlanes[0].Length; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] iUp = slicer.Solid.Intersect(Shifted[0][i]);
                Autodesk.DesignScript.Geometry.Geometry[] iDn = slicer.Solid.Intersect(Shifted[1][i]);

                List<Solid> cutsUp = new List<Solid>();
                List<Solid> cutsDn = new List<Solid>();
                for (int j = 0; j < iUp.Length; j++) { if (iUp[j] is Surface) cutsUp.Add(((Surface)iUp[j].Translate(slicer.CutPlanes[0][i].Normal, -slicer.Thickness / 2.0)).Thicken(slicer.Thickness, true)); else iUp[j].Dispose(); }
                for (int j = 0; j < iDn.Length; j++) { if (iDn[j] is Surface) cutsDn.Add(((Surface)iDn[j].Translate(slicer.CutPlanes[0][i].Normal, slicer.Thickness / 2.0)).Thicken(slicer.Thickness, true)); else iUp[j].Dispose(); }
                if (cutsUp.Count < 1 || cutsDn.Count < 1)
                {
                    iUp.ForEach(s => s.Dispose());
                    iDn.ForEach(s => s.Dispose());
                    cutsUp.ForEach(s => s.Dispose());
                    cutsDn.ForEach(s => s.Dispose());
                    continue;
                }
                else
                {
                    Solid[] Result = new Solid[] { Solid.ByUnion(cutsUp), Solid.ByUnion(cutsDn) };
                    iUp.ForEach(s => s.Dispose());
                    iDn.ForEach(s => s.Dispose());
                    cutsUp.ForEach(s => s.Dispose());
                    cutsDn.ForEach(s => s.Dispose());

                    Autodesk.DesignScript.Geometry.Geometry[] geometry = slicer.CutPlanes[0][i].IntersectAll(Result);
                    List<PolyCurve> profiles = new List<PolyCurve>();
                    List<Surface> surfaces = new List<Surface>();
                    List<Solid> slices = new List<Solid>();
                    for (int j = 0; j < geometry.Length; j++)
                    {
                        if (geometry[j] is Surface)
                        {
                            surfaces.Add((Surface)geometry[j]);
                            profiles.Add(PolyCurve.ByJoinedCurves(((Surface)geometry[j]).PerimeterCurves()));
                            slices.Add(((Surface)geometry[j]).Thicken(slicer.Thickness, true));
                        }
                        else geometry[j].Dispose();
                    }
                    Result.ForEach(r => r.Dispose());

                    Profiles.Add(profiles.ToArray());
                    Surfaces.Add(surfaces.ToArray());
                    Slices.Add(slices.ToArray());
                }
            }

            return new Dictionary<string, object>
            {
                { "Profiles", Profiles },
                { "Surfaces", Surfaces},
                { "Slices", Slices }
            };
        }
        /// <summary>
        /// Returns results as circles, planar surfaces, and solid slices.
        /// Slices contain solid so it can be sanded down.
        /// </summary>
        /// <param name="slicer">slicer object</param>
        /// <returns>circles, planar surfaces, and solid slices listed by layer then by object</returns>
        [MultiReturn(new[] { "Profiles", "Surfaces", "Slices" })]
        public static Dictionary<string, object> GetResultsPlaneCircumscribed(Slicer slicer)
        {
            Plane[][] Shifted = slicer.ShiftCutPlanes(slicer.CutPlanes[0]);

            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();

            for (int i = 0; i < slicer.CutPlanes[0].Length; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] iUp = slicer.Solid.Intersect(Shifted[0][i]);
                Autodesk.DesignScript.Geometry.Geometry[] iDn = slicer.Solid.Intersect(Shifted[1][i]);

                List<Solid> cuts = new List<Solid>();
                for (int j = 0; j < iUp.Length; j++) { if (iUp[j] is Surface) cuts.Add(((Surface)iUp[j].Translate(slicer.CutPlanes[0][i].Normal, -slicer.Thickness / 2.0)).Thicken(slicer.Thickness, true)); else iUp[j].Dispose(); }
                for (int j = 0; j < iDn.Length; j++) { if (iDn[j] is Surface) cuts.Add(((Surface)iDn[j].Translate(slicer.CutPlanes[0][i].Normal, slicer.Thickness / 2.0)).Thicken(slicer.Thickness, true)); else iUp[j].Dispose(); }
                Solid Result = Solid.ByUnion(cuts);
                iUp.ForEach(s => s.Dispose());
                iDn.ForEach(s => s.Dispose());
                cuts.ForEach(s => s.Dispose());

                Autodesk.DesignScript.Geometry.Geometry[] geometry = Result.Intersect(slicer.CutPlanes[0][i]);
                List<PolyCurve> profiles = new List<PolyCurve>();
                List<Surface> surfaces = new List<Surface>();
                List<Solid> slices = new List<Solid>();
                for (int j = 0; j < geometry.Length; j++)
                {
                    if (geometry[j] is Surface)
                    {
                        surfaces.Add((Surface)geometry[j]);
                        profiles.Add(PolyCurve.ByJoinedCurves(((Surface)geometry[j]).PerimeterCurves()));
                        slices.Add(((Surface)geometry[j]).Thicken(slicer.Thickness, true));
                    }
                    else geometry[j].Dispose();
                }
                Result.Dispose();

                Profiles.Add(profiles.ToArray());
                Surfaces.Add(surfaces.ToArray());
                Slices.Add(slices.ToArray());
            }

            return new Dictionary<string, object>
            {
                { "Profiles", Profiles },
                { "Surfaces", Surfaces},
                { "Slices", Slices }
            };
        }
        /// <summary>
        /// Returns results as circles, planar surfaces, and solid slices.
        /// </summary>
        /// <param name="slicer">slicer object</param>
        /// <returns>circles, planar surfaces, and solid slices listed by layer then by object</returns>
        [MultiReturn(new[] { "Profiles", "Surfaces", "Slices" })]
        public static Dictionary<string, object> GetResultsPlane(Slicer slicer)
        {
            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();

            for (int i = 0; i < slicer.CutPlanes[0].Length; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] geometry = slicer.Solid.Intersect(slicer.CutPlanes[0][i]);
                List<PolyCurve> profiles = new List<PolyCurve>();
                List<Surface> surfaces = new List<Surface>();
                List<Solid> slices = new List<Solid>();
                for (int j = 0; j < geometry.Length; j++)
                {
                    if (geometry[j] is Surface)
                    {
                        surfaces.Add((Surface)geometry[j]);
                        profiles.Add(PolyCurve.ByJoinedCurves(((Surface)geometry[j]).PerimeterCurves()));
                        slices.Add(((Surface)geometry[j]).Thicken(slicer.Thickness, true));
                    }
                    else geometry[j].Dispose();
                }
                Profiles.Add(profiles.ToArray());
                Surfaces.Add(surfaces.ToArray());
                Slices.Add(slices.ToArray());
            }

            return new Dictionary<string, object>
            {
                { "Profiles", Profiles },
                { "Surfaces", Surfaces},
                { "Slices", Slices }
            };
        }
        /// <summary>
        /// Returns surfaces by layer
        /// </summary>
        /// <param name="slicer">slicer object</param>
        /// <returns>list of polycurve arrays</returns>
        [MultiReturn(new[] { "Profiles", "Surfaces", "Slices" })]
        public static Dictionary<string, object> GetResultsLine1(Slicer slicer)
        {
            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();

            for (int i = 0; i < slicer.CutPlanes[1].Length; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] geometry = slicer.Solid.Intersect(slicer.CutPlanes[0][i]);
                List<PolyCurve> profiles = new List<PolyCurve>();
                List<Surface> surfaces = new List<Surface>();
                List<Solid> slices = new List<Solid>();
                for (int j = 0; j < geometry.Length; j++)
                {
                    if (geometry[j] is Surface)
                    {
                        surfaces.Add((Surface)geometry[j]);
                        profiles.Add(PolyCurve.ByJoinedCurves(((Surface)geometry[j]).PerimeterCurves()));
                        slices.Add(((Surface)geometry[j]).Thicken(slicer.Thickness, true));
                    }
                    else geometry[j].Dispose();
                }
                Profiles.Add(profiles.ToArray());
                Surfaces.Add(surfaces.ToArray());
                Slices.Add(slices.ToArray());
            }

            return new Dictionary<string, object>
            {
                { "Profiles", Profiles },
                { "Surfaces", Surfaces},
                { "Slices", Slices }
            };
        }
        /// <summary>
        /// Returns surfaces by layer
        /// </summary>
        /// <param name="slicer">slicer object</param>
        /// <returns>list of polycurve arrays</returns>
        [MultiReturn(new[] { "Profiles", "Surfaces", "Slices" })]
        public static Dictionary<string, object> GetResultsLine2(Slicer slicer)
        {
            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();

            for (int i = 0; i < slicer.CutPlanes[2].Length; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] geometry = slicer.Solid.Intersect(slicer.CutPlanes[0][i]);
                List<PolyCurve> profiles = new List<PolyCurve>();
                List<Surface> surfaces = new List<Surface>();
                List<Solid> slices = new List<Solid>();
                for (int j = 0; j < geometry.Length; j++)
                {
                    if (geometry[j] is Surface)
                    {
                        surfaces.Add((Surface)geometry[j]);
                        profiles.Add(PolyCurve.ByJoinedCurves(((Surface)geometry[j]).PerimeterCurves()));
                        slices.Add(((Surface)geometry[j]).Thicken(slicer.Thickness, true));
                    }
                    else geometry[j].Dispose();
                }
                Profiles.Add(profiles.ToArray());
                Surfaces.Add(surfaces.ToArray());
                Slices.Add(slices.ToArray());
            }

            return new Dictionary<string, object>
            {
                { "Profiles", Profiles },
                { "Surfaces", Surfaces},
                { "Slices", Slices }
            };
        }


        //**PRIVATE FUNCTIONS
        /// <summary>
        /// This function creates a list of cut planes;
        /// Initial cut plane is the plane specified in cosntructor
        /// and the rest of the planes are generated from the initial plane
        /// and the given spacing.
        /// </summary>
        /// <param name="plane">initial plane</param>
        /// <returns>List containing all cut planes</returns>
        internal Plane[] GenerateCutPlanes(Plane plane)
        {
            List<Plane> Planes = new List<Plane>();     // List for planes that will be returned
            List<Plane> upPlanes = new List<Plane>();   // List of planes in normal direction
            List<Plane> downPlanes = new List<Plane>(); // List of planes in reverse normal direction
            Plane cutplane = plane.Offset(Spacing);     // Plane object that is iterated to generate lists
            while (Solid.DoesIntersect(cutplane))
            {
                upPlanes.Add(cutplane);
                cutplane = cutplane.Offset(Spacing);
            }
            cutplane = plane.Offset(0.0 - Spacing);
            while (Solid.DoesIntersect(cutplane))
            {
                downPlanes.Add(cutplane);
                cutplane = cutplane.Offset(0.0 - Spacing);
            }
            downPlanes.Reverse();
            Planes.AddRange(downPlanes);
            Planes.Add(plane.Offset(0));
            Planes.AddRange(upPlanes);
            return Planes.ToArray();
        }
        /// <summary>
        /// This function creates a list of cut planes;
        /// Initial cut plane is normal plane to curve at given origin
        /// and the rest of the planes are generated from the initial plane
        /// and the given spacing.
        /// </summary>
        /// <param name="curve"></param>
        /// <param name="origin"></param>
        /// <returns>List containing all cut planes</returns>
        internal Plane[] GenerateCutPlanes(Curve curve, double origin)
        {
            List<Plane> Planes = new List<Plane>();
            if (origin > 0)
            {
                Point pt = curve.PointAtParameter(origin);
                Vector vec = curve.TangentAtParameter(origin);
                Planes.Add(Plane.ByOriginNormal(pt, vec));
                pt.Dispose();
                vec.Dispose();
            }
            else
            {
                Point pt = Solid.Intersect(curve)[0] as Point;
                Vector vec = curve.TangentAtParameter(curve.ParameterAtPoint(pt));
                Planes.Add(Plane.ByOriginNormal(pt, vec));
                pt.Dispose();
                vec.Dispose();
            }
            return Planes.ToArray();
        }
        /// <summary>
        /// Generates a list of arrays of closed result based on cutplanes 
        /// </summary>
        internal List<Surface[]> GenerateSurfaces(Plane[] planes)
        {
            List<Surface[]> Surfaces = new List<Surface[]>();
            for (int i = 0; i < planes.Length; i++)
            {
                List<Surface> surfaces = new List<Surface>();
                Autodesk.DesignScript.Geometry.Geometry[] intersection = Solid.Intersect(planes[i]);
                for (int j = 0; j < intersection.Length; j++)
                {
                    if (intersection[i] is Surface) surfaces.Add((Surface)intersection[i]);
                    else intersection[i].Dispose();
                }
                Surfaces.Add(surfaces.ToArray());
            }
            return Surfaces;
        }
        /// <summary>
        /// returns solid objects
        /// </summary>
        /// <param name="Surfaces"></param>
        /// <returns></returns>
        internal List<Solid[]> GenerateSlices(List<Surface[]> Surfaces)
        {
            List<Solid[]> Slices = new List<Solid[]>(Surfaces.Count);
            for (int i = 0; i < Surfaces.Count; i++) for (int j = 0; j < Surfaces[i].Length; j++)
                    Slices[i][j] = Surfaces[i][j].Thicken(Thickness, true);
            return Slices;
        }
        /// <summary>
        /// Takes a list of surface arrays, and returns perimeters
        /// </summary>
        internal List<PolyCurve[]> GenerateProfiles(List<Surface[]> Surfaces)
        {
            List<PolyCurve[]> Curves = new List<PolyCurve[]>();
            for (int i = 0; i < Surfaces.Count; i++) for (int j = 0; j < Surfaces[i].Length; j++)
                    Curves[i][j] = PolyCurve.ByJoinedCurves(Surfaces[i][j].PerimeterCurves());
            return Curves;
        }
        internal Plane[][] ShiftCutPlanes(Plane[] planes)
        {
            Plane[] shiftUp = new Plane[planes.Length];
            Plane[] shiftDown = new Plane[planes.Length];
            for (int i = 0; i < planes.Length; i++)
            {
                shiftUp[i] = planes[i].Offset(Thickness / 2.0);
                shiftDown[i] = planes[i].Offset(-Thickness / 2.0);
            }
            Plane[][] shifted = { shiftUp, shiftDown };
            return shifted;
        }
        //**METHODS**DISPOSE
        /// <summary>
        /// dispose function
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
        /// <summary>
        /// protected dispose function
        /// </summary>
        /// <param name="disposing"></param>
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                CutPlanes.ForEach(plns => plns.ForEach(pln => pln.Dispose()));
            }
            disposed = true;
        }

    }
}
