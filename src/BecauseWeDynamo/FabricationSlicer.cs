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
    public class Slicer: IDisposable
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
        /// primary cut planes
        /// </summary>
        public List<Plane> CutPlanesPrimary { get; set; }
        /// <summary>
        /// secondary cut planes
        /// </summary>
        public List<Plane> CutPlanesSecondary { get; set; }
        /// <summary>
        /// tertiary cut planes
        /// </summary>
        public List<Plane> CutPlanesTertiary { get; set; }
        /// <summary>
        /// inital geometry
        /// </summary>
        public List<Autodesk.DesignScript.Geometry.Geometry> InitialGeometry { get; set; }

        //**CONSTRUCTORS
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
        {
            Solid = solid;
            Thickness = thickness;
            Spacing = spacing;
            CutPlanesPrimary.AddRange(GenerateCutPlanes(plane));
            InitialGeometry = new List<Autodesk.DesignScript.Geometry.Geometry>();
            InitialGeometry.Add(plane);
            if (line1 != null)
            {
                CutPlanesSecondary.AddRange(GenerateCutPlanes(Plane.ByThreePoints(line1.StartPoint, line1.EndPoint, line1.StartPoint.Add(plane.Normal))));
                InitialGeometry.Add(line1);
            }
            if (line2 != null)
            {
                CutPlanesTertiary.AddRange(GenerateCutPlanes(Plane.ByThreePoints(line2.StartPoint, line2.EndPoint, line2.StartPoint.Add(plane.Normal))));
                InitialGeometry.Add(line2);
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
        internal Slicer(Solid solid, Curve curve, double thickness, double spacing, double origin)
        {
            Solid = solid;
            Thickness = thickness;
            Spacing = spacing;
            Plane plane = Plane.ByOriginNormal(curve.StartPoint, curve.Normal);
            Curve curvePlanar = curve;
            if(!curve.IsPlanar)
            {
                plane = Plane.ByBestFitThroughPoints(curve.ToNurbsCurve().ControlPoints());
                curvePlanar = curve.PullOntoPlane(plane);
            }
            CutPlanesPrimary.AddRange(GenerateCutPlanes(plane));
            CutPlanesSecondary.AddRange(GenerateCutPlanes(curvePlanar, origin));
            InitialGeometry = new List<Autodesk.DesignScript.Geometry.Geometry>(1) { curvePlanar };
        }

        //**CREATE
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
        public static Slicer ByPlaneAndLines(Solid solid, Plane plane, Line line1 = null, Line line2 = null, double thickness = 1, double spacing = 0)
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
        public static Slicer ByCurve(Solid solid, Curve curve, double thickness = 1, double spacing = 0, double origin = 0)
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
        public static Dictionary<string, object> GetPrimaryInscribedResults (Slicer slicer)
        {
            Plane[][] Shifted = slicer.ShiftCutPlanes(slicer.CutPlanesPrimary, slicer.Thickness);
            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();

            for (int i = 0; i < slicer.CutPlanesPrimary.Count; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] Up = slicer.Solid.Intersect(Shifted[0][i]);
                Autodesk.DesignScript.Geometry.Geometry[] Down = slicer.Solid.Intersect(Shifted[1][i]);
                List<Solid> solidsUp = new List<Solid>();
                List<Solid> solidsDown = new List<Solid>();

                for (int j = 0; j < Up.Length; j++)
                {
                    if ((Up[j] is Point) ^ (Up[j] is Curve)) { continue; }
                    else { solidsUp.Add(((Surface)Up[j].Translate(slicer.CutPlanesPrimary[i].Normal.Scale(-slicer.Thickness / 2.0))).Thicken(slicer.Thickness)); }
                }
                for (int j = 0; j < Down.Length; j++)
                {
                    if ((Down[j] is Point) ^ (Down[j] is Curve)) { continue; }
                    else { solidsDown.Add(((Surface)Down[j].Translate(slicer.CutPlanesPrimary[i].Normal.Scale(slicer.Thickness / 2.0))).Thicken(slicer.Thickness)); }
                }

                Solid solidUp = Solid.ByUnion(solidsUp);
                Solid solidDown = Solid.ByUnion(solidsDown);
                Autodesk.DesignScript.Geometry.Geometry[] results = slicer.CutPlanesPrimary[i].IntersectAll(solidUp.Intersect(solidDown));

                List<Surface> surfaces = new List<Surface>();
                List<Solid> slices = new List<Solid>();
                List<PolyCurve> profiles = new List<PolyCurve>();
                for (int j = 0; j < results.Length; j++)
                {
                    if ((results[j] is Point) ^ (results[j] is Curve)) { continue; }
                    else 
                    {
                        surfaces.Add((Surface) results[j]);
                        profiles.Add(PolyCurve.ByJoinedCurves(((Surface) results[j]).PerimeterCurves()));
                        slices.Add(((Surface) results[j]).Thicken(slicer.Thickness)); 
                    }
                }
                Surfaces.Add(surfaces.ToArray());
                Slices.Add(slices.ToArray());
                Profiles.Add(profiles.ToArray());
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
        public static Dictionary<string, object> GetPrimaryCircumscribedResults(Slicer slicer)
        {
            Plane[][] Shifted = slicer.ShiftCutPlanes(slicer.CutPlanesPrimary, slicer.Thickness);
            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();
            for (int i = 0; i < slicer.CutPlanesPrimary.Count; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] Up = slicer.Solid.Intersect(Shifted[0][i]);
                Autodesk.DesignScript.Geometry.Geometry[] Down = slicer.Solid.Intersect(Shifted[1][i]);
                List<Solid> solidsUp = new List<Solid>();
                List<Solid> solidsDown = new List<Solid>();

                for (int j = 0; j < Up.Length; j++)
                {
                    if ((Up[j] is Point) ^ (Up[j] is Curve)) { continue; }
                    else { solidsUp.Add(((Surface)Up[j].Translate(slicer.CutPlanesPrimary[i].Normal.Scale(-slicer.Thickness / 2.0))).Thicken(slicer.Thickness)); }
                }
                for (int j = 0; j < Down.Length; j++)
                {
                    if ((Down[j] is Point) ^ (Down[j] is Curve)) { continue; }
                    else { solidsDown.Add(((Surface)Down[j].Translate(slicer.CutPlanesPrimary[i].Normal.Scale(slicer.Thickness / 2.0))).Thicken(slicer.Thickness)); }
                }

                Solid solidUp = Solid.ByUnion(solidsUp);
                Solid solidDown = Solid.ByUnion(solidsDown);
                Solid result = solidUp.Union(solidDown);
                Autodesk.DesignScript.Geometry.Geometry[] results = result.Intersect(slicer.CutPlanesPrimary[i]);

                List<PolyCurve> profiles = new List<PolyCurve>();
                List<Surface> surfaces = new List<Surface>();
                List<Solid> slices = new List<Solid>();
                for (int j = 0; j < results.Length; j++)
                {
                    if ((results[j] is Point) ^ (results[j] is Curve)) { continue; }
                    else
                    {
                        surfaces.Add((Surface)results[j]);
                        profiles.Add(PolyCurve.ByJoinedCurves(((Surface)results[j]).PerimeterCurves()));
                        slices.Add(((Surface)results[j]).Thicken(slicer.Thickness));
                    }
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
        /// Returns results as circles, planar surfaces, and solid slices.
        /// </summary>
        /// <param name="slicer">slicer object</param>
        /// <returns>circles, planar surfaces, and solid slices listed by layer then by object</returns>
        [MultiReturn(new[] { "Profiles", "Surfaces", "Slices" })]
        public static Dictionary<string, object> GetPrimaryResults(Slicer slicer)
        {
            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();

            for (int i = 0; i < slicer.CutPlanesPrimary.Count; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] geometry = slicer.Solid.Intersect(slicer.CutPlanesPrimary[i]);
                List<PolyCurve> profiles = new List<PolyCurve>();
                List<Surface> surfaces = new List<Surface>();
                List<Solid> slices = new List<Solid>();
                for (int j = 0; j < geometry.Length; j++)
                {
                    if ((geometry[j] is Point) ^ (geometry[j] is Curve)) { continue; }
                    else
                    {
                        surfaces.Add((Surface)geometry[j]);
                        profiles.Add(PolyCurve.ByJoinedCurves(((Surface)geometry[j]).PerimeterCurves()));
                        slices.Add(((Surface)geometry[j]).Thicken(slicer.Thickness));
                    }
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
        [MultiReturn(new[] { "Primary", "Secondary", "Tertiary" })]
        public static Dictionary<string, object> GetSurfaces(Slicer slicer)
        {

            List<Solid[]> SliceResults = new List<Solid[]>(); // Container for Solid Results
            return new Dictionary<string, object>
            {
                { "Primary", slicer.GenerateSurfaces(slicer.CutPlanesPrimary) },
                { "Secondary", slicer.GenerateSurfaces(slicer.CutPlanesSecondary)},
                { "Tertiary", slicer.GenerateSurfaces(slicer.CutPlanesTertiary)}
            };
        }
        /// <summary>
        /// Returns solid slices by layer
        /// </summary>
        /// <param name="slicer">slicr object</param>
        /// <returns>list of solid arrays</returns>
        [MultiReturn(new[] { "Primary", "Secondary", "Tertiary" })]
        public static Dictionary<string, object> GetSlices(Slicer slicer)
        {
            return new Dictionary<string, object>
            {
                { "Primary", slicer.GenerateSlices(slicer.GenerateSurfaces(slicer.CutPlanesPrimary)) },
                { "Secondary", slicer.GenerateSlices(slicer.GenerateSurfaces(slicer.CutPlanesSecondary))},
                { "Tertiary", slicer.GenerateSlices(slicer.GenerateSurfaces(slicer.CutPlanesTertiary))}
            };
        }
        /// <summary>
        /// Returns profiles by layer
        /// </summary>
        /// <param name="slicer">slicer object</param>
        /// <returns>list of polycurve arrays</returns>
        [MultiReturn(new[] { "Primary", "Secondary", "Tertiary" })]
        public static Dictionary<string, object> GetProfiles(Slicer slicer)
        {
            List<Solid[]> SliceResults = new List<Solid[]>(); // Container for Solid Results
            return new Dictionary<string, object>
            {
                { "Primary", slicer.GenerateProfiles(slicer.CutPlanesPrimary) },
                { "Secondary", slicer.GenerateProfiles(slicer.CutPlanesSecondary)},
                { "Tertiary", slicer.GenerateProfiles(slicer.CutPlanesTertiary)}
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
        private List<Plane> GenerateCutPlanes(Plane plane)
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
            Planes.Add(plane);
            Planes.AddRange(upPlanes);
            return Planes;
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
        public List<Plane> GenerateCutPlanes(Curve curve, double origin)
        {
            List<Plane> Planes = new List<Plane>();
            if (origin > 0) Planes.Add(Plane.ByOriginNormal(curve.PointAtParameter(origin), curve.TangentAtParameter(origin)));
            else Planes.Add(Plane.ByOriginNormal((Point)Solid.Intersect(curve)[0], curve.TangentAtParameter(curve.ParameterAtPoint((Point)Solid.Intersect(curve)[0]))));
            return Planes;
        }
        /// <summary>
        /// Generates a list of arrays of closed result based on cutplanes 
        /// </summary>
        public List<Surface[]> GenerateSurfaces(List<Plane> planes)
        {
            List<Surface[]> Surfaces = new List<Surface[]>();
            for (int i = 0; i < planes.Count; i++ )
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
        public List<Solid[]> GenerateSlices(List<Surface[]> Surfaces)
        {
            List<Solid[]> Slices = new List<Solid[]>(Surfaces.Count);
            for (int i = 0; i < Surfaces.Count; i++) for (int j = 0; j < Surfaces[i].Length; j++)
                    Slices[i][j] = Surfaces[i][j].Thicken(Thickness);
            return Slices;
        }
        /// <summary>
        /// Takes a list of surface arrays, and returns perimeters
        /// </summary>
        private List<PolyCurve[]> GenerateProfiles(List<Plane> planes)
        {
            List<PolyCurve[]> Curves = new List<PolyCurve[]>();
            for (int i = 0; i < planes.Count; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] geometry = Solid.Intersect(planes[i]);
                List<PolyCurve> profiles = new List<PolyCurve>();
                for (int j = 0; j < geometry.Length; j++)
                {
                    if ((geometry[j] is Point) ^ (geometry[j] is Curve)) { continue; }
                    else { profiles.Add(PolyCurve.ByJoinedCurves(((Surface)geometry[j]).PerimeterCurves())); }
                }
                Curves.Add(profiles.ToArray());
            }
            return Curves;
        }
        private Plane[][] ShiftCutPlanes(List<Plane> planes, double thickness)
        {
            Plane[] shiftUp = new Plane[planes.Count];
            Plane[] shiftDown = new Plane[planes.Count];
            for (int i = 0; i < planes.Count; i++)
            {
                shiftUp[i] = planes[i].Offset(thickness / 2.0);
                shiftDown[i] = planes[i].Offset(- thickness / 2.0);
            }
            Plane[][] shifted = {shiftUp, shiftDown};
            return shifted;
        }
        private List<Surface[]> GenerateCircumscribedSurfaces(List<Plane> planes)
        {
            Plane[][] Shifted = ShiftCutPlanes(planes, Thickness);
            List<Surface[]> Surfaces = new List<Surface[]>();
            for (int i = 0; i < planes.Count; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] Up = Solid.Intersect(Shifted[0][i]);
                Autodesk.DesignScript.Geometry.Geometry[] Down = Solid.Intersect(Shifted[1][i]);
                List<Solid> solidsUp = new List<Solid>();
                List<Solid> solidsDown = new List<Solid>();
                
                for (int j = 0; j < Up.Length; j++)
                {
                    if ((Up[j] is Point) ^ (Up[j] is Curve)) { continue; }
                    else { solidsUp.Add(((Surface) Up[j].Translate(planes[i].Normal.Scale(-Thickness/2.0))).Thicken(Thickness)); }
                }
                for (int j = 0; j < Down.Length; j++)
                {
                    if ((Down[j] is Point) ^ (Down[j] is Curve)) { continue; }
                    else { solidsDown.Add(((Surface)Down[j].Translate(planes[i].Normal.Scale(Thickness/2.0))).Thicken(Thickness)); }
                }
                
                Solid solidUp = Solid.ByUnion(solidsUp);
                Solid solidDown = Solid.ByUnion(solidsDown);
                Solid result = solidUp.Union(solidDown);
                List<Surface> surfaces = new List<Surface>();
                Autodesk.DesignScript.Geometry.Geometry[] results = result.Intersect(planes[i]);
                for (int j = 0; j < results.Length; j++)
                {
                    if ((results[j] is Point) ^ (results[j] is Curve)) { continue; }
                    else { surfaces.Add((Surface) results[j]); }
                }
                Surfaces.Add(surfaces.ToArray());
            }
            return Surfaces;
        }
        /// <summary>
        /// returns inscribed surfaces at cutplanes
        /// </summary>
        /// <param name="planes"></param>
        /// <returns></returns>
        public List<Surface[]> GenerateInscribedSurfaces(List<Plane> planes)
        {
            Plane[][] Shifted = ShiftCutPlanes(planes, Thickness);
            List<Surface[]> Surfaces = new List<Surface[]>();
            for (int i = 0; i < planes.Count; i++)
            {
                Autodesk.DesignScript.Geometry.Geometry[] Up = Solid.Intersect(Shifted[0][i]);
                Autodesk.DesignScript.Geometry.Geometry[] Down = Solid.Intersect(Shifted[1][i]);
                List<Solid> solidsUp = new List<Solid>();
                List<Solid> solidsDown = new List<Solid>();

                for (int j = 0; j < Up.Length; j++)
                {
                    if ((Up[j] is Point) ^ (Up[j] is Curve)) { continue; }
                    else { solidsUp.Add(((Surface)Up[j].Translate(planes[i].Normal.Scale(-Thickness / 2.0))).Thicken(Thickness)); }
                }
                for (int j = 0; j < Down.Length; j++)
                {
                    if ((Down[j] is Point) ^ (Down[j] is Curve)) { continue; }
                    else { solidsDown.Add(((Surface)Down[j].Translate(planes[i].Normal.Scale(Thickness / 2.0))).Thicken(Thickness)); }
                }

                Solid solidUp = Solid.ByUnion(solidsUp);
                Solid solidDown = Solid.ByUnion(solidsDown);
                List<Surface> surfaces = new List<Surface>();
                Autodesk.DesignScript.Geometry.Geometry[] results = planes[i].IntersectAll(solidUp.Intersect(solidDown));
                for (int j = 0; j < results.Length; j++)
                {
                    if ((results[j] is Point) ^ (results[j] is Curve)) { continue; }
                    else { surfaces.Add((Surface)results[j]); }
                }
                Surfaces.Add(surfaces.ToArray());
            }
            return Surfaces;
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
                Solid.Dispose();
                CutPlanesPrimary.ForEach(p => p.Dispose());
                CutPlanesSecondary.ForEach(p => p.Dispose());
                CutPlanesTertiary.ForEach(p => p.Dispose());
                InitialGeometry.ForEach(p => p.Dispose());
            }
            disposed = true;
        }

    }
}
