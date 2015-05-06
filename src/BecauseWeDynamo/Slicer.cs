using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Runtime;
using System;
using System.Collections.Generic;

namespace Fabrication
{
    public class Slicer
    {
        // **GLOBAL VARIABLES
        private Solid _solid;
        private double _thickness;
        private double _spacing;
        private List<Plane> CutPlanesPrimary = new List<Plane>();
        private List<Plane> CutPlanesSecondary = new List<Plane>();
        private List<Plane> CutPlanesTertiary = new List<Plane>();


        //**QUERY
        /// <summary>
        /// List of initial cut planes from input plane
        /// </summary>
        public List<Plane> GetPlanesPrimary { get { return CutPlanesPrimary; } }
        /// <summary>
        /// List of secondary cut planes
        /// </summary>
        public List<Plane> GetPlanesSecondary { get { return CutPlanesSecondary; } }
        /// <summary>
        /// List of tertiary cut planes
        /// </summary>
        public List<Plane> GetPlanesTertiary { get { return CutPlanesTertiary; } }

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
        /// <param name="mode">mode: 0 - average, 1 - slices inside solid, 2 - slices outside solid</param>
        internal Slicer(Solid solid, Plane plane, Line line1, Line line2, double thickness, double spacing)
        {
            _solid = solid;
            _thickness = thickness;
            _spacing = spacing;

            CutPlanesPrimary.AddRange(GenerateCutPlanes(plane));
            if (line1 != null)
            {
                Line l = (Line) line1;
                CutPlanesSecondary.AddRange(GenerateCutPlanes(Plane.ByThreePoints(l.StartPoint, l.EndPoint, l.StartPoint.Add(plane.Normal))));
            }
            if (line2 != null)
            {
                Line l = (Line) line2;
                CutPlanesTertiary.AddRange(GenerateCutPlanes(Plane.ByThreePoints(l.StartPoint, l.EndPoint, l.StartPoint.Add(plane.Normal))));
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
        /// <returns>A newly-constructed Slicer object</returns>
        internal Slicer(Solid solid, Curve curve, double thickness, double spacing, double origin)
        {
            _solid = solid;
            _thickness = thickness;
            _spacing = spacing;

            Plane plane;
            Curve curvePlanar = curve;

            if(!curve.IsPlanar)
            {
                plane = Plane.ByBestFitThroughPoints(curve.ToNurbsCurve().ControlPoints());
                curvePlanar = curve.PullOntoPlane(plane);
            }
            else
            {
                plane = Plane.ByOriginNormal(curve.StartPoint, curve.Normal);
            }
            CutPlanesPrimary.AddRange(GenerateCutPlanes(plane));
            CutPlanesSecondary.AddRange(GenerateCutPlanes(curvePlanar, origin));

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
        public static Slicer ByPlane(Solid solid, Plane plane, double thickness, double spacing)
        {
            if (!solid.DoesIntersect(plane))
            {
                throw new ArgumentException("plane");
            }
            return new Slicer(solid, plane, null, null, thickness, spacing);
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
        public static Slicer ByLines(Solid solid, Plane plane, double thickness, double spacing, Line line1, Line line2)
        {
            if (!solid.DoesIntersect(plane))
            {
                throw new ArgumentException("plane");
            }
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
        /// <returns>A newly-constructed Slicer object</returns>
        public static Slicer ByCurve(Solid solid, Curve curve, double thickness, double spacing, double origin = 0)
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
            Plane[][] Shifted = slicer.ShiftCutPlanes(slicer.CutPlanesPrimary, slicer._thickness);
            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();

            for (int i = 0; i < slicer.CutPlanesPrimary.Count; i++)
            {
                Geometry[] Up = slicer._solid.Intersect(Shifted[0][i]);
                Geometry[] Down = slicer._solid.Intersect(Shifted[1][i]);
                List<Solid> solidsUp = new List<Solid>();
                List<Solid> solidsDown = new List<Solid>();

                for (int j = 0; j < Up.Length; j++)
                {
                    if ((Up[j] is Point) ^ (Up[j] is Curve)) { continue; }
                    else { solidsUp.Add(((Surface)Up[j].Translate(slicer.CutPlanesPrimary[i].Normal.Scale(-slicer._thickness / 2.0))).Thicken(slicer._thickness)); }
                }
                for (int j = 0; j < Down.Length; j++)
                {
                    if ((Down[j] is Point) ^ (Down[j] is Curve)) { continue; }
                    else { solidsDown.Add(((Surface)Down[j].Translate(slicer.CutPlanesPrimary[i].Normal.Scale(slicer._thickness / 2.0))).Thicken(slicer._thickness)); }
                }

                Solid solidUp = Solid.ByUnion(solidsUp);
                Solid solidDown = Solid.ByUnion(solidsDown);
                Geometry[] results = slicer.CutPlanesPrimary[i].IntersectAll(solidUp.Intersect(solidDown));

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
                        slices.Add(((Surface) results[j]).Thicken(slicer._thickness)); 
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
            Plane[][] Shifted = slicer.ShiftCutPlanes(slicer.CutPlanesPrimary, slicer._thickness);
            List<Surface[]> Surfaces = new List<Surface[]>();
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            List<Solid[]> Slices = new List<Solid[]>();
            for (int i = 0; i < slicer.CutPlanesPrimary.Count; i++)
            {
                Geometry[] Up = slicer._solid.Intersect(Shifted[0][i]);
                Geometry[] Down = slicer._solid.Intersect(Shifted[1][i]);
                List<Solid> solidsUp = new List<Solid>();
                List<Solid> solidsDown = new List<Solid>();

                for (int j = 0; j < Up.Length; j++)
                {
                    if ((Up[j] is Point) ^ (Up[j] is Curve)) { continue; }
                    else { solidsUp.Add(((Surface)Up[j].Translate(slicer.CutPlanesPrimary[i].Normal.Scale(-slicer._thickness / 2.0))).Thicken(slicer._thickness)); }
                }
                for (int j = 0; j < Down.Length; j++)
                {
                    if ((Down[j] is Point) ^ (Down[j] is Curve)) { continue; }
                    else { solidsDown.Add(((Surface)Down[j].Translate(slicer.CutPlanesPrimary[i].Normal.Scale(slicer._thickness / 2.0))).Thicken(slicer._thickness)); }
                }

                Solid solidUp = Solid.ByUnion(solidsUp);
                Solid solidDown = Solid.ByUnion(solidsDown);
                Solid result = solidUp.Union(solidDown);
                Geometry[] results = result.Intersect(slicer.CutPlanesPrimary[i]);

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
                        slices.Add(((Surface)results[j]).Thicken(slicer._thickness));
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
                Geometry[] geometry = slicer._solid.Intersect(slicer.CutPlanesPrimary[i]);
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
                        slices.Add(((Surface)geometry[j]).Thicken(slicer._thickness));
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
                { "Primary", slicer.GenerateSlices(slicer.CutPlanesPrimary) },
                { "Secondary", slicer.GenerateSlices(slicer.CutPlanesSecondary)},
                { "Tertiary", slicer.GenerateSlices(slicer.CutPlanesTertiary)}
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
            /// List for planes that will be returned
            List<Plane> Planes = new List<Plane>();

            /// List of planes in normal direction
            List<Plane> upPlanes = new List<Plane>();

            /// List of planes in reverse normal direction
            List<Plane> downPlanes = new List<Plane>();

            /// Plane object that is iterated to generate lists
            Plane cutplane = plane.Offset(_spacing);

            /// Conditional statement adding planes to list upPlanes
            while (_solid.DoesIntersect(cutplane))
            {
                upPlanes.Add(cutplane);
                cutplane = cutplane.Offset(_spacing);
            }

            /// Reset iterating plane
            cutplane = plane.Offset(0.0 - _spacing);

            /// Conditional statement adding planes to list downPlanes
            while (_solid.DoesIntersect(cutplane))
            {
                downPlanes.Add(cutplane);
                cutplane = cutplane.Offset(0.0 - _spacing);
            }

            /// Changes direction of downPlanes so indexing
            /// is consistent to the normal vector.
            downPlanes.Reverse();

            /// Lists are combined in relation to normal vector
            /// so list of planes starts at bottom and moves in normal direction 
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
        private List<Plane> GenerateCutPlanes(Curve curve, double origin)
        {
            // List for planes that will be returned
            List<Plane> Planes = new List<Plane>();

            //Initial Cut Plane
            if (origin > 0)
            {
                Planes.Add(Plane.ByOriginNormal(curve.PointAtParameter(origin), curve.TangentAtParameter(origin)));
            }
            else
            {
                Point pt = (Point)_solid.Intersect(curve)[0];
                Planes.Add(Plane.ByOriginNormal(pt, curve.TangentAtParameter(curve.ParameterAtPoint(pt))));
            }

            //***CRAWL THROUGH SOLID ALONG CURVE AT GIVEN SPACING
            
            return Planes;
        }
        /// <summary>
        /// Generates a list of arrays of closed result based on cutplanes 
        /// </summary>
        /// <param name="planes">A List containing cutplanes</param>
        /// <param name="solid">Solid being compared</param>
        /// <returns></returns>
        private List<Surface[]> GenerateSurfaces(List<Plane> planes)
        {
            List<Surface[]> Surfaces = new List<Surface[]>();
            for (int i = 0; i < planes.Count; i++ )
            {
                Geometry[] geometry = _solid.Intersect(planes[i]);
                List<Surface> surfaces = new List<Surface>();
                for (int j = 0; j < geometry.Length; j++) 
                {
                    if ((geometry[j] is Point) ^ (geometry[j] is Curve)) { continue; }
                    else { surfaces.Add((Surface) geometry[j]); }
                }
                Surfaces.Add(surfaces.ToArray());
            }
            return Surfaces;
        }
        /// <summary>
        /// This function creates a list of arrays of solids;
        /// Each layer is grouped by cutplane in an array,
        /// which is then stored in a list.
        /// </summary>
        /// <param name="planes">A List of arrays of closed nurbs result that outline slices</param>
        /// <returns>List containing all cut planes</returns>
        private List<Solid[]> GenerateSlices(List<Plane> planes)
        {
            List<Solid[]> Slices = new List<Solid[]>();
            for (int i = 0; i < planes.Count; i++)
            {
                Geometry[] geometry = _solid.Intersect(planes[i]);
                List<Solid> solids = new List<Solid>();
                for (int j = 0; j < geometry.Length; j++)
                {
                    if ((geometry[j] is Point) ^ (geometry[j] is Curve)) { continue; }
                    else { solids.Add(((Surface)geometry[j]).Thicken(_thickness)); }
                }
                Slices.Add(solids.ToArray());
            }
            return Slices;
        }
        /// <summary>
        /// Takes a list of surface arrays, and returns perimeters
        /// </summary>
        /// <param name="Surfaces"></param>
        /// <returns>A list of polycurve arrays with profiles </returns>
        private List<PolyCurve[]> GenerateProfiles(List<Plane> planes)
        {
            List<PolyCurve[]> Curves = new List<PolyCurve[]>();
            for (int i = 0; i < planes.Count; i++)
            {
                Geometry[] geometry = _solid.Intersect(planes[i]);
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
            Plane[][] Shifted = ShiftCutPlanes(planes, _thickness);
            List<Surface[]> Surfaces = new List<Surface[]>();
            for (int i = 0; i < planes.Count; i++)
            {
                Geometry[] Up = _solid.Intersect(Shifted[0][i]);
                Geometry[] Down = _solid.Intersect(Shifted[1][i]);
                List<Solid> solidsUp = new List<Solid>();
                List<Solid> solidsDown = new List<Solid>();
                
                for (int j = 0; j < Up.Length; j++)
                {
                    if ((Up[j] is Point) ^ (Up[j] is Curve)) { continue; }
                    else { solidsUp.Add(((Surface) Up[j].Translate(planes[i].Normal.Scale(-_thickness/2.0))).Thicken(_thickness)); }
                }
                for (int j = 0; j < Down.Length; j++)
                {
                    if ((Down[j] is Point) ^ (Down[j] is Curve)) { continue; }
                    else { solidsDown.Add(((Surface)Down[j].Translate(planes[i].Normal.Scale(_thickness/2.0))).Thicken(_thickness)); }
                }
                
                Solid solidUp = Solid.ByUnion(solidsUp);
                Solid solidDown = Solid.ByUnion(solidsDown);
                Solid result = solidUp.Union(solidDown);
                List<Surface> surfaces = new List<Surface>();
                Geometry[] results = result.Intersect(planes[i]);
                for (int j = 0; j < results.Length; j++)
                {
                    if ((results[j] is Point) ^ (results[j] is Curve)) { continue; }
                    else { surfaces.Add((Surface) results[j]); }
                }
                Surfaces.Add(surfaces.ToArray());
            }
            return Surfaces;
        }

        private List<Surface[]> GenerateInscribedSurfaces(List<Plane> planes)
        {
            Plane[][] Shifted = ShiftCutPlanes(planes, _thickness);
            List<Surface[]> Surfaces = new List<Surface[]>();
            for (int i = 0; i < planes.Count; i++)
            {
                Geometry[] Up = _solid.Intersect(Shifted[0][i]);
                Geometry[] Down = _solid.Intersect(Shifted[1][i]);
                List<Solid> solidsUp = new List<Solid>();
                List<Solid> solidsDown = new List<Solid>();

                for (int j = 0; j < Up.Length; j++)
                {
                    if ((Up[j] is Point) ^ (Up[j] is Curve)) { continue; }
                    else { solidsUp.Add(((Surface)Up[j].Translate(planes[i].Normal.Scale(-_thickness / 2.0))).Thicken(_thickness)); }
                }
                for (int j = 0; j < Down.Length; j++)
                {
                    if ((Down[j] is Point) ^ (Down[j] is Curve)) { continue; }
                    else { solidsDown.Add(((Surface)Down[j].Translate(planes[i].Normal.Scale(_thickness / 2.0))).Thicken(_thickness)); }
                }

                Solid solidUp = Solid.ByUnion(solidsUp);
                Solid solidDown = Solid.ByUnion(solidsDown);
                List<Surface> surfaces = new List<Surface>();
                Geometry[] results = planes[i].IntersectAll(solidUp.Intersect(solidDown));
                for (int j = 0; j < results.Length; j++)
                {
                    if ((results[j] is Point) ^ (results[j] is Curve)) { continue; }
                    else { surfaces.Add((Surface)results[j]); }
                }
                Surfaces.Add(surfaces.ToArray());
            }
            return Surfaces;
        }

    }
}
