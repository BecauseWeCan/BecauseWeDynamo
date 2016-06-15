using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using System;
using System.Collections.Generic;

namespace Fabrication
{
     class zNesting
    {
        //**TODO SUMMARY 
        /* ***Implement code to do actual nesting based on nesting code from existing freeware
         * ***Vertical Tile is set as default action for list of circles
         * ***Add option for rotation attempts to better optimize nesting
         * ***Add option of adding a sheet offset
         * ***Add option of preferred nesting direction
         * ***Potentially adding constructors for list of surfaces, or list of polycurve arrays
         */
 

        //**GLOBAL VARIABLES
        private List<PolyCurve> ListPolyCurve;
        private PolyCurve Sheet;
        private double Spacing;
        private double Width;
        private double Length;

        //**QUERY
        /// <summary>
        /// get internal list of circles
        /// </summary>
        public List<PolyCurve> PolyCurves { get { return ListPolyCurve; } }

        //**CONSTRUCTORS
        /// <summary>
        /// get nesting of circles within a set sheet
        /// </summary>
        /// <param name="circles">input result to be cut</param>
        /// <param name="spacing">spacing between input result</param>
        /// <param name="sheet">closed curve of sheet</param>
        internal zNesting(List<PolyCurve> polycurves, double spacing, PolyCurve sheet)
        {
            ListPolyCurve = TileVertical(MapToXYPlane(polycurves), spacing);
            Sheet = sheet;
            Spacing = spacing;
        }
        /// <summary>
        /// get nesting of circles within a set sheet
        /// </summary>
        /// <param name="surfaces">input planar surfaces to be cut</param>
        /// <param name="spacing">spacing between input result</param>
        /// <param name="sheet">closed curve of sheet</param>
        internal zNesting(List<Surface[]> surfaces, double spacing, PolyCurve sheet)
        {
            ListPolyCurve = TileVertical(MapToXYPlane(surfaces), spacing);
            Sheet = sheet;
            Spacing = spacing;
        }
        /// <summary>
        /// get nesting of circles within a set sheet
        /// </summary>
        /// <param name="surfaces">input planar surfaces to be cut</param>
        /// <param name="spacing">spacing between input result</param>
        /// <param name="width">width of sheet</param>
        /// <param name="length">length of sheet</param>
        internal zNesting(List<Surface[]> surfaces, double spacing, double width, double length)
        {
            Width = width;
            Length = length;
            Spacing = spacing;
            ListPolyCurve = TileVertical(MapToXYPlane(surfaces), spacing);
        }
        /// <summary>
        /// get nesting of circles within given dimensions
        /// </summary>
        /// <param name="circles">input planar surfaces to be cut</param>
        /// <param name="spacing">spacing between input result</param>
        /// <param name="width">width of sheet</param>
        /// <param name="length">length of sheet</param>
        internal zNesting(List<PolyCurve> polycurves, double spacing, double width, double length)
        {
            Width = width;
            Length = length;
            Spacing = spacing;
            ListPolyCurve = TileVertical(MapToXYPlane(polycurves), spacing);
        }

        //**CREATE
        /// <summary>
        /// create nesting layout of closed polycurve profiles within a closed polycurve sheet with given spacing
        /// </summary>
        /// <param name="circles">list of closed planar circles</param>
        /// <param name="spacing">(double) spacing between circles</param>
        /// <param name="sheet">closed planar polycurve of sheet outline</param>
        /// <returns>new instance of nesting</returns>
        public static zNesting ByPolyCurvesInPolyCurve(List<PolyCurve> polycurves, double spacing, PolyCurve sheet)
        {
            // error for sheet not closed or planar
            if (!sheet.IsClosed && !sheet.IsPlanar)
            {
                throw new ArgumentException("invalid sheet");
            }
            // error for empty list of input result to nest
            if (polycurves.Count <= 0)
            {
                throw new ArgumentException("no curves to array");
            }
            // error for input result not closed or planar
            for (int i = 0; i < polycurves.Count; i++)
            {
                if (!polycurves[i].IsClosed && !polycurves[i].IsPlanar)
                {
                    throw new ArgumentException("invalid curves");
                }
            }

            // NEEDS WORK:
            // error for input result not fitting inside sheet

            return new zNesting(polycurves, spacing, sheet);
        }
        /// <summary>
        /// create nesting layout of closed polycurve profiles within a closed polycurve sheet with given spacing
        /// </summary>
        /// <param name="Surfaces">list of closed planar surface arrays</param>
        /// <param name="spacing">(double) spacing between circles</param>
        /// <param name="sheet">closed planar polycurve of sheet outline</param>
        /// <returns>new instance of nesting</returns>
        public static zNesting BySurfacesInPolyCurve(List<Surface[]> Surfaces, double spacing, PolyCurve sheet)
        {
            // error for sheet not closed or planar
            if (!sheet.IsClosed && !sheet.IsPlanar)
            {
                throw new ArgumentException("invalid sheet");
            }
            // error for empty list of surfaces to nest
            if (Surfaces.Count <= 0)
            {
                throw new ArgumentException("no surfaces in array");
            }

            // NEEDS WORK:
            // error for input result not fitting inside sheet

            return new zNesting(Surfaces, spacing, sheet);
        }
        /// <summary>
        /// create nesting layout of closed polycurve profiles within a closed polycurve sheet with given spacing
        /// </summary>
        /// <param name="Surfaces">list of closed planar surface arrays</param>
        /// <param name="spacing">(double) spacing between circles</param>
        /// <param name="width">(double) width of sheet</param>
        /// <returns>new instance of nesting</returns>
        public static zNesting BySurfacesInWidth(List<Surface[]> Surfaces, double spacing, double width)
        {
            // error for width being too small
            if (!(width > 0))
            {
                throw new ArgumentException("invalid width");
            }
            // error for empty list of surfaces to nest
            if (Surfaces.Count <= 0)
            {
                throw new ArgumentException("no surfaces in array");
            }

            // NEEDS WORK:
            // error for input result not fitting inside sheet

            return new zNesting(Surfaces, spacing, width, 0);
        }
        /// <summary>
        /// create nesting layout of closed polycurve profiles within a closed polycurve sheet with given spacing
        /// </summary>
        /// <param name="Surfaces">list of closed planar surface arrays</param>
        /// <param name="spacing">(double) spacing between circles</param>
        /// <param name="width">(double) width of sheet</param>
        /// <returns>new instance of nesting</returns>
        public static zNesting BySurfacesInWidthAndLength(List<Surface[]> Surfaces, double spacing, double width, double length)
        {
            // error for width or length being too small
            if ((width <= 0) ^ (length <= 0))
            {
                throw new ArgumentException("invalid sheet dimension");
            }
            // error for empty list of surfaces to nest
            if (Surfaces.Count <= 0)
            {
                throw new ArgumentException("no surfaces in array");
            }

            // NEEDS WORK:
            // error for input result not fitting inside sheet

            return new zNesting(Surfaces, spacing, width, length);
        }

        //**ACTIONS
        /// <summary>
        /// get list of circles and list of lines from circles
        /// </summary>
        /// <param name="circles">list of closed planar polycurve profiles</param>
        /// <returns>list of circles and lines</returns>
        [MultiReturn(new[] { "Arcs", "Lines"})]
        public static Dictionary<string, object> ApproximateWithArcAndLineSegments(List<PolyCurve> polycurves)
        {
            // get list of circles and lines from circles
            List<Curve> Curves = new List<Curve>();
            // iterate over circles
            for (int i = 0; i < polycurves.Count; i++)
            {
                // iterate over result in circles
                for (int k = 0; k < polycurves[i].Curves().Length; k++)
                {
                    Curves.AddRange(polycurves[i].Curves()[k].ApproximateWithArcAndLineSegments());
                }
            }


            List<Line> Lines = new List<Line>();
            List<Arc> Arcs = new List<Arc>();
            // iterate over list of lines and circles and filter into list of circles and list of lines
            for (int i = 0; i < Curves.Count; i++)
            {
                double stIndex = Curves[i].StartParameter();
                double endIndex = Curves[i].EndParameter();
                Vector stTan = Curves[i].TangentAtParameter(stIndex);
                Vector endTan = Curves[i].TangentAtParameter(endIndex);
                Point stpt = Curves[i].StartPoint;
                Point endpt = Curves[i].EndPoint;

                if (stTan.IsParallel(endTan))
                {
                    Line line = (Line)Curves[i];
                    Lines.Add(line);
                }
                else
                {
                    Arc arc = Arc.ByStartPointEndPointStartTangent(stpt, endpt, stTan);
                    Arcs.Add(arc);
                }

            }

            return new Dictionary<string, object>
            {
                { "Arcs", Arcs },
                { "Lines", Lines},
            };
        }
        //**ACTIONS - NEEDS CODE
        /// <summary>
        /// takes a nesting object and returns circles on sheets
        /// </summary>
        /// <param name="nesting">nesting object</param>
        /// <returns>list of sheets and list of circles on sheet</returns>
        [MultiReturn(new[] { "Profiles", "Sheet" })]
        public static Dictionary<string, object> SheetLayout(zNesting nesting)
        {
            // list of polycurve arrays by sheet
            List<PolyCurve[]> Profiles = new List<PolyCurve[]>();
            // list of sheet
            List<PolyCurve> Sheets = new List<PolyCurve>();

            return new Dictionary<string, object>
            {
                { "Profiles", Profiles },
                { "Sheet", Sheets},
            };
        }

        //**INTERNAL FUNCTIONS
        /// <summary>
        /// maps surfaces to XYplane: iterate over list of surface arrays
        /// and map each surface onto XY plane while maintaining index structure.
        /// </summary>
        /// <param name="surfacesInput">list of closed planar surfaces</param>
        /// <returns>list of closed surface arrays on XYplane</returns>
        private List<Surface[]> MapToXYPlane(List<Surface[]> surfacesInput)
        {
            List<Surface[]> surfacesResult = new List<Surface[]>();
            for (int i = 0; i < surfacesInput.Count; i++ )
            {
                Surface[] surfaces = new Surface[surfacesInput[i].Length];
                for (int j = 0; j < surfacesInput[i].Length; j++)
                {
                    surfaces[j] = (Surface) surfacesInput[i][j].Transform(surfacesInput[i][j].CoordinateSystemAtParameter(0.5, 0.5), CoordinateSystem.Identity());
                }
                surfacesResult.Add(surfaces);
            }

            return surfacesResult;
        }
        /// <summary>
        /// maps surfaces to XYplane: iterate over list of surface arrays
        /// and map each surface onto XY plane while maintaining index structure.
        /// </summary>
        /// <param name="polycurvesInput">list of closed planar surfaces</param>
        /// <returns>list of closed surface arrays on XYplane</returns>
        private List<PolyCurve> MapToXYPlane(List<PolyCurve> polycurvesInput)
        {
            List<PolyCurve> polycurvesResult = new List<PolyCurve>();
            for (int i = 0; i < polycurvesInput.Count; i++)
            {
                polycurvesResult.Add( (PolyCurve) polycurvesInput[i].Transform( polycurvesInput[i].CoordinateSystemAtParameter(0.5), CoordinateSystem.Identity() ) );
            }
            return polycurvesResult;
        }
        /// <summary>
        /// tiles objects vertically
        /// </summary>
        /// <param name="ListSurfaces">list of closed planar surface arrays</param>
        /// <param name="spacing">(double) spacing between result</param>
        /// <returns>list of closed planar circles</returns>
        private List<PolyCurve> TileVertical(List<Surface[]> ListSurfaces, double spacing)
        {
            List<PolyCurve> Profiles = new List<PolyCurve>();

            Point LastMaxPoint = Point.Origin();
            Point LastMinPoint = Point.Origin();

            for (int i = 0; i < ListSurfaces.Count; i++)
            {
                for (int j = 0; j < ListSurfaces[i].Length; j++)
                {
                    Vector move = Vector.ByTwoPoints(ListSurfaces[i][j].BoundingBox.MinPoint, Point.ByCoordinates(spacing, LastMaxPoint.Y + spacing));
                    Surface surface = (Surface)ListSurfaces[i][j].Translate(move);
                        Profiles.Add(PolyCurve.ByJoinedCurves(surface.PerimeterCurves()));
                        LastMaxPoint = surface.BoundingBox.MaxPoint;
                        LastMinPoint = surface.BoundingBox.MinPoint;
                }
            }

            return Profiles;
        }
        /// <summary>
        /// tiles objects vertically
        /// </summary>
        /// <param name="ListPolyCurve">list of closed planar surface arrays</param>
        /// <param name="spacing">(double) spacing between result</param>
        /// <returns>list of closed planar circles</returns>
        private List<PolyCurve> TileVertical(List<PolyCurve> ListPolyCurve, double spacing)
        {
            List<PolyCurve> Profiles = new List<PolyCurve>();

            Point LastMaxPoint = Point.Origin();
            Point LastMinPoint = Point.Origin();

            for (int i = 0; i < ListPolyCurve.Count; i++)
            {
                    Vector move = Vector.ByTwoPoints(ListPolyCurve[i].BoundingBox.MinPoint, Point.ByCoordinates(spacing, LastMaxPoint.Y + spacing));
                    Surface surface = (Surface)ListPolyCurve[i].Translate(move);
                    Profiles.Add(PolyCurve.ByJoinedCurves(surface.PerimeterCurves()));
                    LastMaxPoint = surface.BoundingBox.MaxPoint;
                    LastMinPoint = surface.BoundingBox.MinPoint;
            }

            return Profiles;
        }
    }
}
