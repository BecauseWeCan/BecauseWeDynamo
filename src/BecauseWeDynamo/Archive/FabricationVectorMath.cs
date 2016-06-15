using System;
using Autodesk.DesignScript.Geometry;

namespace Fabrication
{
    public static class VectorMath
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
        /// returns angle between vectors in radians
        /// angle is equal to or less than 180
        /// </summary>
        /// <param name="V1">Vector 1</param>
        /// <param name="V2">Vector 2</param>
        /// <returns>Angle in Radians</returns>
        public static double AngleBetween(Vector V1, Vector V2)
        {
            Vector N1 = V1.Normalized();
            Vector N2 = V2.Normalized();
            double a = Math.Acos(VectorMath.Dot(N1, N2));
            N1.Dispose(); N2.Dispose();
            return a;
        }

        /// <summary>
        /// returns dot product between two vectors
        /// </summary>
        /// <param name="V1">Vector 1</param>
        /// <param name="V2">Vector 2</param>
        /// <returns>Dot Product</returns>
        public static double Dot(Vector V1, Vector V2)
        {
            return V1.X * V2.X + V1.Y * V2.Y + V1.Z * V2.Z;
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
        public static Point Center(Point[] Points)
        {
            double[] XYZ = new double[] { 0, 0, 0 };
            for (int i = 0; i < Points.Length; i++)
            {
                XYZ[0] += Points[i].X / Points.Length;
                XYZ[1] += Points[i].Y / Points.Length;
                XYZ[2] += Points[i].Z / Points.Length;
            }
            return Point.ByCoordinates(XYZ[0], XYZ[1], XYZ[2]);
        }

        /// <summary>
        /// find weighted center of point cluster based on given points and weights
        /// </summary>
        /// <param name="Points"></param>
        /// <param name="Weights"></param>
        /// <returns></returns>
        public static Point CenterWeighted(Point[] Points, double[] Weights)
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
            return Point.ByCoordinates(XYZ[0]/W, XYZ[1]/W, XYZ[2]/W);
        }

        /// <summary>
        /// find center vector of vector cluster
        /// </summary>
        /// <param name="Vectors"></param>
        /// <returns></returns>
        public static Vector Center(Vector[] Vectors)
        {
            double[] XYZ = new double[] { 0, 0, 0 };
            for (int i = 0; i < Vectors.Length; i++)
            {
                XYZ[0] += Vectors[i].X / Vectors.Length;
                XYZ[1] += Vectors[i].Y / Vectors.Length;
                XYZ[2] += Vectors[i].Z / Vectors.Length;
            }
            return Vector.ByCoordinates(XYZ[0], XYZ[1], XYZ[2]);
        }

        /// <summary>
        /// find midpoint of two points
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static Point Center(Point A, Point B)
        {
            return Point.ByCoordinates(A.X / 2 + B.X / 2, A.Y / 2 + B.Y / 2, A.Z / 2 + B.Z / 2);
        }

        /// <summary>
        /// takes an orthogonal basis or object coordinate system (OCS) 
        /// and returns orthonormal basis as object coordinate
        /// </summary>
        /// <param name="CS">Object Coordinate System</param>
        /// <returns>Normalized OCS</returns>
        public static Vector[] NormalizeCS(Vector[] CS)
        {
            return new Vector[] { CS[0].Normalized(), CS[1].Normalized(), CS[2].Normalized() };
        }
        
        /// <summary>
        /// returns normalized vector addition
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="V2"></param>
        /// <returns></returns>
        public static Vector NormalizedAdd(Vector V1, Vector V2)
        {
            return V1.Add(V2).Normalized();
        }

        /// <summary>
        /// returns normalized vector subtraction
        /// </summary>
        /// <param name="V1"></param>
        /// <param name="V2"></param>
        /// <returns></returns>
        public static Vector NormalizedSubtract(Vector V1, Vector V2)
        {
            return V1.Subtract(V2).Normalized();
        }

        /// <summary>
        /// return normalized cross vector of two vectors
        /// </summary>
        /// <param name="V1">Vector 1</param>
        /// <param name="V2">Vector 2</param>
        /// <returns></returns>
        public static Vector NormalizedCross(Vector V1, Vector V2)
        {
            double v = Math.Sqrt(Math.Pow(V1.Y * V2.Z - V2.Y * V1.Z, 2) + Math.Pow(V1.Z * V2.X - V2.Z * V1.X, 2) + Math.Pow(V1.X * V2.Y - V2.X * V1.Y, 2));
            Vector V = Vector.ByCoordinates((V1.Y * V2.Z - V2.Y * V1.Z)/v, (V1.Z * V2.X - V2.Z * V1.X)/v, (V1.X * V2.Y - V2.X * V1.Y)/v);
            return V;
        }

        /// <summary>
        /// creates right-handed orthogonal vector to normal
        /// that lies in the global XY plane
        /// </summary>
        /// <param name="Normal"></param>
        /// <returns></returns>
        public static Vector NormalizedZeroVector(Vector Normal)
        {
            return Vector.ByCoordinates(-Math.Sign(Normal.Z) * Normal.Y, Math.Sign(Normal.Z) * Normal.X, 0).Normalized();
        }

        /// <summary>
        /// creates vector basis based on normal vector with normal as Z-axis
        /// where the X-Axis is the right-handed orthogonal vector to normal
        /// that lies in the global XY plane with the normal being the Z-axis;
        /// the Y-axis is generated by right-hand rule from Z cross X
        /// </summary>
        /// <param name="Normal"></param>
        /// <returns></returns>
        public static Vector[] NormalizedBasisXY(Vector Normal)
        {
            Vector[] result = new Vector[3];
            result[2] = Normal.Normalized();
            result[0] = NormalizedZeroVector(Normal);
            result[1] = NormalizedCross(result[2],result[0]);
            return result;
        }

        /// <summary>
        /// creates orthonormal vector basis
        /// based on given X-Axis and Y-Axis Vectors
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <returns></returns>
        public static Vector[] NormalizedBasis(Vector X, Vector Y)
        {
            Vector[] result = new Vector[3];
            result[0] = X.Normalized();
            result[2] = NormalizedCross(X, Y);
            result[1] = NormalizedCross(result[2], result[0]);
            return result;
        }

        /// <summary>
        /// creates vector basis based on normal vector as defined in dxf OCS
        /// if x-coordinate and y-coordinate of normal vector is within 1/64 of 0
        /// then X-Axis is Vector(0,1,0) cross normal vector
        /// otherwise, X-Axis is (0,0,1) cross normal vector
        /// the Y-axis is generated by right-hand rule from Z cross X
        /// </summary>
        /// <param name="Normal"></param>
        /// <returns></returns>
        public static Vector[] NormalizedBasisDXF(Vector Normal)
        {
            Vector[] XYZ = new Vector[3];
            XYZ[2] = Normal.Normalized();
            if (Math.Abs(XYZ[2].X) < 1.0/64 && Math.Abs(XYZ[2].Y) < 1.0/ 64) XYZ[0] = NormalizedCross(Vector.YAxis(), XYZ[2]);
            else XYZ[0] = NormalizedCross(Vector.ZAxis(), XYZ[2]);
            XYZ[1] = NormalizedCross(XYZ[2], XYZ[0]);
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
        public static Vector[] NormalizedVertexVectors(Vector V1, Vector V2, Vector N)
        {
            Vector V1N = V1.Normalized();
            Vector V2N = V2.Normalized();
            Vector[] V = new Vector[]{ V1N, V2N, NormalizedAdd(V1N,V2N), NormalizedSubtract(V1N,V2N), NormalizedAdd(NormalizedCross(N, V1N), NormalizedCross(V2N, N)) };
            N.Dispose();
            return V;
        }

        /// <summary>
        /// calculates coefficients of vector from universal basis
        /// to vector basis with given vectors as basis
        /// </summary>
        /// <param name="V">Vector in Universal Basis</param>
        /// <param name="X">X-axis of Basis</param>
        /// <param name="Y">Y-axis of Basis</param>
        /// <param name="Z">Z-axis of Basis</param>
        /// <returns>Coefficients in Given Basis</returns>
        public static double[] ChangeBasis(Vector V, Vector X, Vector Y, Vector Z)
        {
            if (X.Length * Y.Length * Z.Length == 0) return new double[] { 0, 0, 0 };
            return new double[] { Dot(V, X) / Dot(X, X), Dot(V, Y) / Dot(Y, Y), Dot(V, Z) / Dot(Z, Z) };
        }

        /// <summary>
        /// calculates coefficients of vector from universal basis
        /// to vector basis with given vectors as basis
        /// </summary>
        /// <param name="V">Vector to Decompose</param>
        /// <param name="XYZ">Orthogonal Basis as Vector Array [X,Y,Z]</param>
        /// <returns>Coefficients for Decomposition (a, b, c) where V = aX + bY + cZ</returns>
        public static double[] ChangeBasis(Vector V, Vector[] XYZ)
        {
            return ChangeBasis(V, XYZ[0], XYZ[1], XYZ[2]);
        }

        /// <summary>
        /// calculates coefficients of vector from universal basis
        /// to vector basis provided by given coordinate system
        /// </summary>
        /// <param name="V">Vector in Universal Basis</param>
        /// <param name="CS">Coordinate System</param>
        /// <returns>Coefficients in Coordinate System Basis</returns>
        public static double[] ChangeBasis(Vector V, CoordinateSystem CS)
        {
            return ChangeBasis(V, CS.XAxis, CS.YAxis, CS.ZAxis);
        }

        /// <summary>
        /// calculates new coordinates of point from universal basis
        /// to vector basis with given vectors as basis
        /// </summary>
        /// <param name="Point">Point in Universal Basis</param>
        /// <param name="XYZ">Orthogonal Basis as Vector Array [X,Y,Z]</param>
        /// <returns>Point Coordinates in Orthogonal Basis (a, b, c) where Point = aX + bY + cZ</returns>
        public static double[] ChangeBasis(Point Point, Vector[] XYZ)
        {
            Vector V = Vector.ByTwoPoints(Point.Origin(), Point);
            double[] result = ChangeBasis(V, XYZ);
            V.Dispose();
            return result;
        }

        /// <summary>
        /// calculates new coordinates of point from universal basis
        /// to vector basis provided by given coordinate system
        /// </summary>
        /// <param name="Point">Point in Universal Basis</param>
        /// <param name="CS">Coordinate System</param>
        /// <returns>Point Coordinates in Coordinate System Basis</returns>
        public static double[] ChangeBasis(Point Point, CoordinateSystem CS)
        {
            Vector V = Vector.ByTwoPoints(Point.Origin(), Point);
            double[] result = ChangeBasis(V, CS);
            V.Dispose();
            return result;
        }

        /// <summary>
        /// calculates new coordinates of point from universal basis
        /// to vector basis derived from given normal using NormalizedBasisDXF
        /// </summary>
        /// <param name="Point">Point in Universal Basis</param>
        /// <param name="Normal">Normal Vector</param>
        /// <returns>Point Coordinates in Coordinate System Basis</returns>
        public static double[] ChangeBasisDXF(Point Point, Vector Normal)
        {
            return ChangeBasis(Vector.ByCoordinates(Point.X,Point.Y,Point.Z), NormalizedBasisDXF(Normal));
        }

        /// <summary>
        /// calculates new coordinates of point from universal basis
        /// to vector basis derived from given normal using NormalizedBasisDXF
        /// </summary>
        /// <param name="V">Vector in Universal Basis</param>
        /// <param name="Normal">Normal Vector</param>
        /// <returns>Coefficients in Coordinate System Basis</returns>
        public static double[] ChangeBasisDXF(Vector Vector, Vector Normal)
        {
            return ChangeBasis(Vector, NormalizedBasisDXF(Normal));
        }

        /// <summary>
        /// checks if two vectors are parallel
        /// using methods in cross product
        /// </summary>
        /// <param name="V1">Vector 1</param>
        /// <param name="V2">Vector 2</param>
        /// <returns></returns>
        public static bool IsParallel(Vector V1, Vector V2)
        {
            double x = V1.Y * V2.Z - V2.Y * V1.Z;
            double y = V1.Z * V2.X - V2.Z * V1.X;
            double z = V1.X * V2.Y - V2.X * V1.Y;
            return (x + y + z).Equals(0);
        }

        /// <summary>
        /// checks if two vectors are perpendicular using dot product
        /// (dot product of perpendicular vectors is 0)
        /// </summary>
        /// <param name="V1">Vector 1</param>
        /// <param name="V2">Vector 2</param>
        /// <returns></returns>
        public static bool IsPerpendicular(Vector V1, Vector V2)
        {
            return Dot(V1, V2).Equals(0);
        }
    }

}
