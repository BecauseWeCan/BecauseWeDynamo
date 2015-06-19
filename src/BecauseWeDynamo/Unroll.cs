using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Fabrication
{
    class Unroll
    {
        //**TODOS
        //***Implement code copy from existing unfolding scripts from Inventor or 123Dmake

        //**GLOBAL VARIABLES
        private List<Surface> Surfaces;

        //**QUERY
        /// <summary>
        /// get Unrolled results
        /// </summary>
        public List<Surface> UnrolledSurface { get { return Surfaces; } }

        //**CONSTRUCTOR
        /// <summary>
        /// stores list of surfaces to unroll and unrollability
        /// </summary>
        /// <param name="surfaces">list of surfaces to unroll</param>
        /// <param name="canUnroll">true: singly curved, false: doubly curved</param>
        internal Unroll(List<Surface> surfaces, double thickness, bool canUnroll, double tolerance)
        {
            if (!canUnroll)
            {
                surfaces.GetRange(0, surfaces.Count);
                List<Surface> Surfaces = UnrollSurfaces(ApproximateWithUnrollableSurfaces(surfaces, tolerance), thickness);
            }
            else
            {
                List<Surface> Surfaces = UnrollSurfaces(surfaces, thickness);
            }
        }

        //**CREATE
        /// <summary>
        /// creates instance of unrolling
        /// </summary>
        /// <param name="surface">input surface</param>
        /// <returns>Unroll object</returns>
        public static Unroll BySurface(Surface surface, double thickness, double tolerance = 1)
        {
            bool canUnroll = true;
            for (float i = 0; i < 1; i += 0.1f)
            {
                for (float j = 0; j < 1; j += 0.1f)
                {
                    if (surface.GaussianCurvatureAtParameter(i, j) > 0)
                    {
                        canUnroll = false;
                        continue;
                    }
                }
            }
            Surface[] surfaces = {surface};
            return new Unroll(surfaces.ToList(), thickness, canUnroll, tolerance);
        }
        /// <summary>
        /// creates instance of unrolling
        /// </summary>
        /// <param name="solid">input solid</param>
        /// <returns>Unroll object</returns>
        public static Unroll BySolid(Solid solid, double thickness, double tolerance = 1)
        {
            bool canUnroll = true;
            List<Surface> surfaces = new List<Surface>();
            Geometry[] Parts = solid.Explode();
            for (int n = 0; n < Parts.Length; n++)
            {
                surfaces.Add((Surface) Parts[n]);
                for (float i = 0; i < 1; i += 0.1f)
                {
                    for (float j = 0; j < 1; j += 0.1f)
                    {
                        if (surfaces[n].GaussianCurvatureAtParameter(i, j) > 0)
                        {
                            canUnroll = false;
                            continue;
                        }
                    }
                }
            }
            return new Unroll(surfaces, thickness, canUnroll, tolerance);
        }
        /// <summary>
        /// creates instance of unrolling
        /// </summary>
        /// <param name="polysurface">input polysurface</param>
        /// <returns>Unroll object</returns>
        public static Unroll ByPolySurface(PolySurface polysurface, double thickness, double tolerance = 1)
        {
            bool canUnroll = true;
            List<Surface> surfaces = new List<Surface>();
            Geometry[] Parts = polysurface.Explode();
            for (int n = 0; n < Parts.Length; n++)
            {
                surfaces.Add((Surface)Parts[n]);
                for (float i = 0; i < 1; i += 0.1f)
                {
                    for (float j = 0; j < 1; j += 0.1f)
                    {
                        if (surfaces[n].GaussianCurvatureAtParameter(i, j) > 0)
                        {
                            canUnroll = false;
                            continue;
                        }
                    }
                }
            }
            return new Unroll(surfaces, thickness, canUnroll, tolerance);
        }

        //**ACTIONS

        //**INTERNAL FUNCTIONS
        /// <summary>
        /// maps surfaces to XYplane
        /// </summary>
        /// <param name="surfacesInput">list of closed planar surfaces</param>
        /// <returns>list of closed surface arrays on XYplane</returns>
        private List<Surface> MapToXYPlane(List<Surface> surfacesInput)
        {
            List<Surface> surfaces = new List<Surface>();
            for (int i = 0; i < surfacesInput.Count; i++)
            {
                surfaces.Add((Surface) surfacesInput[i].Transform(surfacesInput[i].CoordinateSystemAtParameter(0.5, 0.5), CoordinateSystem.Identity()));
            }
            return surfaces;
        }
        /// <summary>
        /// tiles vertically
        /// </summary>
        /// <param name="ListSurfaces">list of closed planar surface arrays</param>
        /// <param name="spacing">(double) spacing between result</param>
        /// <returns>list of closed planar circles</returns>
        private List<Surface> TileVertical(List<Surface> ListSurfaces, double spacing)
        {
            List<Surface> surfaces = new List<Surface>();

            Point LastMaxPoint = Point.Origin();
            Point LastMinPoint = Point.Origin();

            for (int i = 0; i < ListSurfaces.Count; i++)
            {
                Vector move = Vector.ByTwoPoints(ListSurfaces[i].BoundingBox.MinPoint, Point.ByCoordinates(spacing, LastMaxPoint.Y + spacing));
                surfaces.Add((Surface) ListSurfaces[i].Translate(move));
                LastMaxPoint = surfaces[i].BoundingBox.MaxPoint;
                LastMinPoint = surfaces[i].BoundingBox.MinPoint;
            }

            return surfaces;
        }
        //***NEEDS WORK***
        //***add code to approximate durface with singly curved surfaces within tolerance
        /// <summary>
        /// approximates unrollable surface with rollable ones within tolerance
        /// </summary>
        /// <param name="inputSurfaces">input surface list</param>
        /// <param name="tolerance">maximum difference between original and approximated surfaces</param>
        /// <returns></returns>
        private List<Surface> ApproximateWithUnrollableSurfaces(List<Surface> inputSurfaces, double tolerance)
        {
            List<Surface> surfaces = new List<Surface>();
            return surfaces;
        }
        //***NEEDS WORK***
        //***add code to unroll surface!
        /// <summary>
        /// unrolls surface with consideration to thickness
        /// </summary>
        /// <param name="inputSurfaces">input surface list</param>
        /// <param name="thickness">material thickness</param>
        /// <returns></returns>
        private List<Surface> UnrollSurfaces(List<Surface> inputSurfaces, double thickness)
        {
            List<Surface> surfaces = new List<Surface>();
            //***code to unrol surfaces
            return TileVertical(MapToXYPlane(surfaces), 2*thickness);
        }

    }
}
