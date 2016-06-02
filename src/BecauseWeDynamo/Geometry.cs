using System;
using System.Collections.Generic;

namespace Geometry
{
    /// <summary>
    /// additional math functions
    /// </summary>
    public class math
    {
        /// <summary>
        /// PI as double
        /// </summary>
        public static double PI = 3.14159265358979323846;

        /// <summary>
        /// returns modulo operation on a double
        /// with given input number and modulus
        /// also the remainder after Euclidean division
        /// </summary>
        /// <param name="Number">Input Number</param>
        /// <param name="Modulus">Modulus</param>
        /// <returns>Remainder > 0</returns>
        public static double Mod(double Number, double Modulus)
        {
            return Number - Math.Abs(Modulus) * Math.Floor(Number / Math.Abs(Modulus));
        }

        /// <summary>
        /// takes angle in degrees and converts into radians
        /// </summary>
        /// <param name="Angle">Angle in Degrees</param>
        /// <returns>Angle in Radians</returns>
        public static double toRadians(double Angle)
        {
            return Angle * PI / 180;
        }

        /// <summary>
        /// takes angle in radians and converts into degrees
        /// </summary>
        /// <param name="Angle">Angle in Radians</param>
        /// <returns>Angle in Degrees</returns>
        public static double toDegrees(double Angle)
        {
            return Angle * 180 / PI;
        }

        /// <summary>
        /// scale SAT imports (m to ft)
        /// </summary>
        /// <param name="Geometry"></param>
        /// <returns></returns>
        public static Autodesk.DesignScript.Geometry.Geometry ScaleSAT(Autodesk.DesignScript.Geometry.Geometry Geometry)
        {
            return Geometry.Scale(25 / 2.54 / 3);
        }

        /// <summary>
        /// find center point of point cluster
        /// </summary>
        /// <param name="Points"></param>
        /// <returns></returns>
        public static point Center(point[] Points)
        {
            double[] XYZ = new double[] { 0, 0, 0 };
            for (int i = 0; i < Points.Length; i++)
            {
                XYZ[0] += Points[i].X / Points.Length;
                XYZ[1] += Points[i].Y / Points.Length;
                XYZ[2] += Points[i].Z / Points.Length;
            }
            return point.ByCoordinates(XYZ[0], XYZ[1], XYZ[2]);
        }

        /// <summary>
        /// find weighted center of point cluster based on given points and weights
        /// </summary>
        /// <param name="Points"></param>
        /// <param name="Weights"></param>
        /// <returns></returns>
        public static point CenterWeighted(point[] Points, double[] Weights)
        {
            double[] XYZ = new double[] { 0, 0, 0 };
            double W = 0;
            for (int i = 0; i < Points.Length; i++)
            {
                XYZ[0] += Weights[i] * Points[i].X / Points.Length;
                XYZ[1] += Weights[i] * Points[i].Y / Points.Length;
                XYZ[2] += Weights[i] * Points[i].Z / Points.Length;
                W += Weights[i];
            }
            return point.ByCoordinates(XYZ[0] / W, XYZ[1] / W, XYZ[2] / W);
        }

        /// <summary>
        /// find center vector of vector cluster
        /// </summary>
        /// <param name="Vectors"></param>
        /// <returns></returns>
        public static vector Center(vector[] Vectors)
        {
            double[] XYZ = new double[] { 0, 0, 0 };
            for (int i = 0; i < Vectors.Length; i++)
            {
                XYZ[0] += Vectors[i].X / Vectors.Length;
                XYZ[1] += Vectors[i].Y / Vectors.Length;
                XYZ[2] += Vectors[i].Z / Vectors.Length;
            }
            return vector.ByCoordinates(XYZ[0], XYZ[1], XYZ[2]);
        }

        /// <summary>
        /// find midpoint of two points
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static point Center(point A, point B)
        {
            return point.ByCoordinates(A.X / 2 + B.X / 2, A.Y / 2 + B.Y / 2, A.Z / 2 + B.Z / 2);
        }
    }
    /// <summary>
    /// coordinates
    /// </summary>
    public class coordinates
    {
        /// <summary>
        /// x-coordinate
        /// </summary>
        public double X { get; set; }
        /// <summary>
        /// y-coordinate
        /// </summary>
        public double Y { get; set; }
        /// <summary>
        /// z-coordinate
        /// </summary>
        public double Z { get; set; }

        internal coordinates()
        { X = 0; Y = 0; Z = 0; }
        internal coordinates(double X, double Y, double Z)
        { this.X = X; this.Y = Y; this.Z = Z; }
    }
    /// <summary>
    /// vector object from coordinates
    /// </summary>
    public class vector: coordinates
    {
        /// <summary>
        /// vector length
        /// </summary>
        public double Length { get { return Math.Sqrt(X * X + Y * Y + Z * Z); } }

        internal vector() : base() { }
        internal vector(double X, double Y, double Z) : base(X, Y, Z) { }

        /// <summary>
        /// vector from coordinate object
        /// </summary>
        public static vector ByCoordinates(coordinates Coordinates) { return new vector(Coordinates.X, Coordinates.Y, Coordinates.Z); }
        /// <summary>
        /// vector from coordinates
        /// </summary>
        public static vector ByCoordinates(double X, double Y, double Z = 0) { return new vector(X, Y, Z); }
        /// <summary>
        /// vector from two points (tail=base, head=arrow)
        /// </summary>
        /// <returns></returns>
        public static vector ByTwoPoints(point Tail, point Head) { return new vector(Head.X - Tail.X, Head.Y - Tail.Y, Head.Z - Tail.Z); }
        /// <summary>
        /// global x-axis (1,0,0)
        /// </summary>
        public static vector Xaxis() { return new vector(1, 0, 0); }
        /// <summary>
        /// global y-axis (0,1,0)
        /// </summary>
        public static vector Yaxis() { return new vector(0, 1, 0); }
        /// <summary>
        /// global z-axis (0,0,1)
        /// </summary>
        public static vector Zaxis() { return new vector(0, 0, 1); }
        /// <summary>
        /// orthogonal x-vector to given normal (z-vector) on xy-plane by righthand rule
        /// </summary>
        /// <param name="Normal"></param>
        /// <returns></returns>
        public static vector NormalizedZeroVector(vector Normal) { return vector.ByCoordinates(-Math.Sign(Normal.Z) * Normal.Y, Math.Sign(Normal.Z) * Normal.X, 0).Normalized(); }
        /// <summary>
        /// creates vector basis based on normal vector with normal as Z-axis
        /// where the X-Axis is the right-handed orthogonal vector to normal
        /// that lies in the global XY plane with the normal being the Z-axis;
        /// the Y-axis is generated by right-hand rule from Z cross X
        /// </summary>
        public static vector[] NormalizedBasisXY(vector Normal)
        {
            vector[] result = new vector[3];
            result[2] = Normal.Normalized();
            result[0] = NormalizedZeroVector(Normal);
            result[1] = result[2].NormalizedCross(result[0]);
            return result;
        }
        /// <summary>
        /// creates orthonormal vector basis
        /// based on given X-Axis and Y-Axis Vectors
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <returns></returns>
        public static vector[] NormalizedBasis(vector X, vector Y) { return new vector[] { X.Normalized(), X.NormalizedCross(Y).NormalizedCross(X.Normalized()), X.NormalizedCross(Y) }; }        
        /// <summary>
        /// creates vector basis based on normal vector as defined in dxf OCS
        /// if x-coordinate and y-coordinate of normal vector is within 1/64 of 0
        /// then X-Axis is Vector(0,1,0) cross normal vector
        /// otherwise, X-Axis is (0,0,1) cross normal vector
        /// the Y-axis is generated by right-hand rule from Z cross X
        /// </summary>
        public static vector[] NormalizedBasisDXF(vector Normal)
        {
            vector[] XYZ = new vector[3];
            XYZ[2] = Normal.Normalized();
            if (Math.Abs(XYZ[2].X) < 1.0 / 64 && Math.Abs(XYZ[2].Y) < 1.0 / 64) XYZ[0] = vector.Yaxis().NormalizedCross(XYZ[2]);
            else XYZ[0] = vector.Zaxis().NormalizedCross(XYZ[2]);
            XYZ[1] = XYZ[2].NormalizedCross(XYZ[0]);
            return XYZ;
        }
        /// <summary>
        /// given two vectors defining shape vertex,
        /// V1 is vector from vertex in right-hand rule
        /// V2 is vector from vertex in other direction
        /// where face normal is V1 x V2 (cross product).
        /// returns array {V1, V2, V1+V2 (bisector), V1-V2 (ON vector to biscetor and normal), NxV1 + V2xN(exterior bisector)}
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="V2"></param>
        /// <param name="N"></param>
        /// <returns></returns>
        public static vector[] NormalizedVertexVectors(vector V1, vector V2, vector N)
        {
            return new vector[] { V1.Normalized(), V2.Normalized(), V1.Normalized().NormalizedAdd(V2.Normalized()), V1.Normalized().NormalizedSubtract(V2.Normalized()), N.NormalizedCross(V1.Normalized()).NormalizedAdd(V2.Normalized().NormalizedCross(N)) };
        }
        
        /// <summary>
        /// angle between vectors in radians (always less than PI)
        /// </summary>
        /// <param name="Vector"></param>
        /// <returns></returns>
        public double AngleBetween(vector Vector)
        {
            vector N1 = this.Normalized();
            vector N2 = Vector.Normalized();
            return Math.Acos(N1.Dot(N2));
        }
        /// <summary>
        /// dot product of two vectors (order invariant)
        /// </summary>
        public double Dot(vector Vector) { return X * Vector.X + Y * Vector.Y + Z * Vector.Z; }
        /// <summary>
        /// vector addition
        /// </summary>
        public vector Add(vector Vector) { return new vector(X + Vector.X, Y + Vector.Y, Z + Vector.Z); }
        /// <summary>
        /// vector subtraction
        /// </summary>
        public vector Subtract(vector Vector) { return new vector(X - Vector.X, Y - Vector.Y, Z - Vector.Z); }
        /// <summary>
        /// uniform scalar multiplication
        /// </summary>
        public vector Scale(double Factor) { return new vector(Factor * X, Factor * Y, Factor * Z); }
        /// <summary>
        /// scalar multiplication by axis
        /// </summary>
        public vector Scale(double xFactor, double yFactor, double zFactor) { return new vector(xFactor * X, yFactor * Y, zFactor * Z); }
        /// <summary>
        /// vector cross product
        /// </summary>
        public vector Cross(vector Vector) { return vector.ByCoordinates(Y * Vector.Z - Vector.Y * Z, Z * Vector.X - Vector.Z * X, X * Vector.Y - Vector.X * Y); }
        /// <summary>
        /// reverse vector direction
        /// </summary>
        public vector Reverse() { return vector.ByCoordinates( -X, -Y, -Z); }
        /// <summary>
        /// normalize vector
        /// </summary>
        public vector Normalized() { return new vector(X / Length, Y / Length, Z / Length); }
        /// <summary>
        /// normalized cross product
        /// </summary>
        public vector NormalizedCross(vector Vector) {  return this.Normalized().Cross(Vector.Normalized()).Normalized(); }
        /// <summary>
        /// takes an orthogonal basis or object coordinate system (OCS) 
        /// and returns orthonormal basis as object coordinate
        /// </summary>
        /// <param name="CS">Object Coordinate System</param>
        /// <returns>Normalized OCS</returns>
        public static vector[] NormalizeCS(vector[] CS) { return new vector[] { CS[0].Normalized(), CS[1].Normalized(), CS[2].Normalized() }; }
        /// <summary>
        /// returns normalized vector addition
        /// </summary>
        public vector NormalizedAdd(vector Vector) { return this.Add(Vector).Normalized(); }
        /// <summary>
        /// returns normalized vector subtraction
        /// </summary>
        public vector NormalizedSubtract(vector V2) { return this.Subtract(V2).Normalized(); }
        /// <summary>
        /// calculates coefficients of vector from universal basis
        /// to vector basis with given vectors as basis
        /// </summary>
        /// <param name="X">X-axis of Basis</param>
        /// <param name="Y">Y-axis of Basis</param>
        /// <param name="Z">Z-axis of Basis</param>
        /// <returns>Coefficients in Given Basis</returns>
        public double[] ChangeBasis(vector X, vector Y, vector Z) { return new double[] { this.Dot(X) / X.Dot(X), this.Dot(Y) / Y.Dot(Y), this.Dot(Z) / Z.Dot(Z) }; }
        /// <summary>
        /// calculates coefficients of vector from universal basis
        /// to vector basis with given vectors as basis
        /// </summary>
        /// <param name="XYZ">Orthogonal Basis as Vector Array [X,Y,Z]</param>
        /// <returns>Coefficients for Decomposition (a, b, c) where V = aX + bY + cZ</returns>
        public double[] ChangeBasis(vector[] XYZ) { return this.ChangeBasis( XYZ[0], XYZ[1], XYZ[2]); }
        /// <summary>
        /// calculates new coordinates of point from universal basis
        /// to vector basis derived from given normal using NormalizedBasisDXF
        /// </summary>
        public double[] ChangeBasisDXF(vector Normal)
        {
            return this.ChangeBasis(NormalizedBasisDXF(Normal));
        }
        /// <summary>
        /// checks if two vectors are parallel
        /// using methods in cross product
        /// </summary>
        public bool IsParallel(vector V2)
        {
            double x = Y * V2.Z - V2.Y * Z;
            double y = Z * V2.X - V2.Z * X;
            double z = X * V2.Y - V2.X * Y;
            return (x + y + z).Equals(0);
        }
        /// <summary>
        /// checks if two vectors are perpendicular using dot product
        /// (dot product of perpendicular vectors is 0)
        /// </summary>
        /// <param name="V2">Vector 2</param>
        /// <returns></returns>
        public bool IsPerpendicular(vector V2)
        {
            return this.Dot(V2).Equals(0);
        }
        /// <summary>
        /// returns Dynamo Vector equivalent
        /// </summary>
        /// <returns></returns>
        public Autodesk.DesignScript.Geometry.Vector ToVector()
        {
            return Autodesk.DesignScript.Geometry.Vector.ByCoordinates(X, Y, Z);
        }
        /// <summary>
        /// gets object equality
        /// </summary>
        public override bool Equals(Object Object) { return this.Equals(Object as vector); }
        /// <summary>
        /// gets vector equality
        /// </summary>
        public bool Equals(vector Vector)
        {
            if (Object.ReferenceEquals(Vector, null)) return false;
            if (Object.ReferenceEquals(this, Vector)) return true;
            if (this.GetType() != Vector.GetType()) return false;
            return (X == Vector.X && Y == Vector.Y && Z == Vector.Z);
        }
        /// <summary>
        /// gets hashcode
        /// </summary>
        public override int GetHashCode() { return string.Format("{0}-{1}-{2}-{3}", X, Y, Z, "v").GetHashCode(); }

    }

    /// <summary>
    /// point object from coordinates
    /// </summary>
    public class point: coordinates
    {
        /// <summary>
        /// vector equivalent
        /// </summary>
        public vector AsVector { get { return vector.ByCoordinates(X, Y, Z); } }

        internal point() : base() { }
        internal point(double X, double Y, double Z) : base(X, Y, Z) { }
        /// <summary>
        /// creates point from given coordinate object
        /// </summary>
        public static point ByCoordinates(coordinates C) { return new point(C.X, C.Y, C.Z); }
        /// <summary>
        /// creates point from given coordinates
        /// </summary>
        public static point ByCoordinates(double X, double Y, double Z=0) { return new point(X, Y, Z); }
        /// <summary>
        /// creates point from Dynamo Point object
        /// </summary>
        public static point ByPoint(Autodesk.DesignScript.Geometry.Point Point) { return point.ByCoordinates(Point.X,Point.Y,Point.Z); }
        /// <summary>
        /// distance to point
        /// </summary>
        public double DistanceTo(point point)
        {
            double x = X - point.X;
            double y = Y - point.Y;
            double z = Z - point.Z;
            return Math.Sqrt(x * x + y * y + z * z);
        }
        /// <summary>
        /// new coordinates for point in given basis
        /// </summary>
        public double[] ChangeBasis(vector X, vector Y, vector Z) { return new double[] { this.AsVector.Dot(X) / X.Dot(X), this.AsVector.Dot(Y) / Y.Dot(Y), this.AsVector.Dot(Z) / Z.Dot(Z) }; }
        /// <summary>
        /// new coordinates for point in given basis
        /// </summary>
        public double[] ChangeBasis(vector[] XYZ) { return this.ChangeBasis(XYZ[0], XYZ[1], XYZ[2]); }
        /// <summary>
        /// scalar mulitplication
        /// </summary>
        public point Multiply(double Factor) { return point.ByCoordinates(Factor * X, Factor * Y, Factor * Z);}
        /// <summary>
        /// point addition
        /// </summary>
        /// <param name="Point"></param>
        /// <returns></returns>
        public point Add(point Point) { return point.ByCoordinates(X + Point.X, Y + Point.Y, Z + Point.Z); }
        /// <summary>
        /// vector addition for point
        /// </summary>
        public point Add(vector Vector) { return point.ByCoordinates(X + Vector.X, Y + Vector.Y, Z + Vector.Z); }
        /// <summary>
        /// vector subtraction for point
        /// </summary>
        public point Subtract(vector Vector) { return point.ByCoordinates(X - Vector.X, Y - Vector.Y, Z - Vector.Z); }
        /// <summary>
        /// returns equivalent Dynamo Point object
        /// </summary>
        public Autodesk.DesignScript.Geometry.Point ToPoint()
        {
            return Autodesk.DesignScript.Geometry.Point.ByCoordinates(X, Y, Z);
        }
        /// <summary>
        /// gets object equality
        /// </summary>
        /// <param name="Object"></param>
        /// <returns></returns>
        public override bool Equals(Object Object) { return this.Equals(Object as point); }
        /// <summary>
        /// gets point equality
        /// </summary>
        /// <param name="Point"></param>
        /// <returns></returns>
        public bool Equals(point Point)
        {
            if (Object.ReferenceEquals(Point, null)) return false;
            if (Object.ReferenceEquals(this, Point)) return true;
            if (this.GetType() != Point.GetType()) return false;
            return (X == Point.X && Y == Point.Y && Z == Point.Z);
        }
        /// <summary>
        /// gets hashcode
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode() { return string.Format("{0}-{1}-{2}", X, Y, Z).GetHashCode(); }
    }

    /// <summary>
    /// coordinate system
    /// </summary>
    public class coordinatesystem
    {
        /// <summary>
        /// x-axis
        /// </summary>
        public vector X { get; set; }
        /// <summary>
        /// y-axis
        /// </summary>
        public vector Y { get; set; }
        /// <summary>
        /// z-axis
        /// </summary>
        public vector Z { get; set; }
        /// <summary>
        /// origin
        /// </summary>
        public point Origin { get; set; }

        internal coordinatesystem()
        {
            Origin = point.ByCoordinates(0, 0, 0);
            X = vector.Xaxis();
            Y = vector.Yaxis();
            Z = vector.Zaxis();
        }
        internal coordinatesystem(point origin): base()
        {
            Origin = origin;
        }
        internal coordinatesystem(point origin, vector X, vector Y, vector Z)
        {
            Origin = origin;
            this.X = X;
            this.Y = Y;
            this.Z = Z;
        }
        /// <summary>
        /// creates universal coordinate system with origin at given point
        /// </summary>
        public static coordinatesystem ByOrigin(point origin) { return new coordinatesystem(origin); }
        /// <summary>
        /// creates universal coordinate system with origin at given point
        /// </summary>
        public static coordinatesystem ByOrigin(double x = 0, double y = 0, double z = 0) { return new coordinatesystem(point.ByCoordinates(x,y,z)); }
        /// <summary>
        /// creates coordinate system base on given point and vectors
        /// </summary>
        public static coordinatesystem ByOriginVectors(point origin, vector xAxis, vector yAxis, vector zAxis) { return new coordinatesystem(origin, xAxis, yAxis, zAxis); }
    }

    /// <summary>
    /// arc object
    /// </summary>
    public class arc
    {
        /// <summary>
        /// arc radius
        /// </summary>
        public double Radius { get; set; }
        /// <summary>
        /// arc center
        /// </summary>
        public point Center { get; set; }
        /// <summary>
        /// arc sweep angle in radians
        /// </summary>
        public double SweepAngle { get; set; }
        /// <summary>
        /// arc normal by righthand rule
        /// </summary>
        public vector Normal { get; set; }
        /// <summary>
        /// context coordinate system with arc start point as x-axis
        /// </summary>
        public coordinatesystem CS { get; set; }

        internal arc(coordinates A, coordinates B, coordinates C)
        {
            vector vA = vector.ByCoordinates(A);
            vector vB = vector.ByCoordinates(B);
            vector vC = vector.ByCoordinates(C);
            point pA = point.ByCoordinates(A);
            point pB = point.ByCoordinates(B);
            point pC = point.ByCoordinates(C);
            double denom = 2 * Math.Pow(vA.Subtract(vB).Cross(vB.Subtract(vC)).Length, 2);
            Radius = vA.Subtract(vB).Length * vB.Subtract(vC).Length * vC.Subtract(vA).Length / denom;
            double nA = vB.Subtract(vC).Dot(vB.Subtract(vC)) * vA.Subtract(vB).Dot(vA.Subtract(vC)) / denom;
            double nB = vA.Subtract(vC).Dot(vA.Subtract(vC)) * vB.Subtract(vA).Dot(vB.Subtract(vC)) / denom;
            double nC = vA.Subtract(vB).Dot(vA.Subtract(vB)) * vC.Subtract(vA).Dot(vC.Subtract(vB)) / denom;
            Center = pA.Multiply(nA).Add(pB.Multiply(nB)).Add(pC.Multiply(nC));

            vector cA = vector.ByTwoPoints(Center, pA);
            vector cB = vector.ByTwoPoints(Center, pB);
            vector cC = vector.ByTwoPoints(Center, pC);

            vector nAB = cA.NormalizedCross(cB);
            vector nBC = cB.NormalizedCross(cC);
            vector nCA = cC.NormalizedCross(cA);

            if (nAB.Dot(nBC) > 0)
            {
                Normal = nAB;
                SweepAngle = cA.AngleBetween(cB) + cB.AngleBetween(cC);
                CS = coordinatesystem.ByOriginVectors(Center, cA.Normalized(), Normal.NormalizedCross(cA), Normal);
            }
            else if (nBC.Dot(nCA) > 0)
            {
                Normal = nBC;
                SweepAngle = 2*Math.PI - cA.AngleBetween(cC);
                CS = coordinatesystem.ByOriginVectors(Center, cA.Normalized(), Normal.NormalizedCross(cA), Normal);
            }
            else if (nCA.Dot(nAB) > 0)
            {
                Normal = nAB;
                SweepAngle = 2*Math.PI - cA.AngleBetween(cC);
                CS = coordinatesystem.ByOriginVectors(Center, cA.Normalized(), Normal.NormalizedCross(cA), Normal);
            }
        }
        /// <summary>
        /// create arc from three points;
        /// first point determines x-axis with normal derived by righthand rule
        /// </summary>
        public static arc ByThreePoints(coordinates A, coordinates B, coordinates C) { return new arc(A, B, C); }
    }

    /// <summary>
    /// circle object
    /// </summary>
    public class circle
    {
        /// <summary>
        /// circle radius
        /// </summary>
        public double Radius { get; set; }
        /// <summary>
        /// circle center
        /// </summary>
        public point Center { get; set; }
        /// <summary>
        /// circle normal set by input
        /// </summary>
        public vector Normal { get; set; }

        internal circle(point Center, double Radius, vector Normal)
        {
            this.Radius = Radius;
            this.Center = Center;
            this.Normal = Normal;
        }
        /// <summary>
        /// creates circle with center at given point with given radius on plane orthogonal to normal
        /// </summary>
        public static circle ByCenterRadiusNormal(point Center, double Radius, vector Normal) { return new circle(Center, Radius, Normal); }
        /// <summary>
        /// gets Dynamo Circle object equivalent
        /// </summary>
        public Autodesk.DesignScript.Geometry.Circle ToCircle()
        {
            return Autodesk.DesignScript.Geometry.Circle.ByCenterPointRadiusNormal(Center.ToPoint(), Radius, Normal.ToVector());
        }
    }
}
