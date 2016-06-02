using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Geometry;

namespace BecauseWeDynamoTest
{
    [TestClass]
    public class vectorTest
    {
        [TestMethod]
        public void NormalizedCross()
        {
            vector X = vector.ByCoordinates(2, 2, 0);
            vector Y = vector.ByCoordinates(-2, 2, 0);
            vector[] V = vector.NormalizedBasis(X, Y);
            vector XN = X.Normalized();
            vector YN = Y.Normalized();
            vector ZN = XN.Cross(YN).Normalized();
            vector[] expected = new vector[] { XN, YN, ZN };
            // assert
            Assert.AreEqual(expected[0], V[0]);
            //Assert.AreEqual(expected[1], V[1]);
            //Assert.AreEqual(expected[2], V[2]);
        }

        [TestMethod]
        public void ArcSweepAngle270()
        {
            point A = point.ByCoordinates(1, 0, 0);
            point B = point.ByCoordinates(0, 1, 0);
            point C = point.ByCoordinates(0, -1, 0);
            arc a = arc.ByThreePoints(A,B,C);

            Assert.AreEqual(a.SweepAngle, math.toRadians(270), 0.00001);
        }
        [TestMethod]
        public void ArcNormal270()
        {
            point A = point.ByCoordinates(1, 0, 0);
            point B = point.ByCoordinates(0, 1, 0);
            point C = point.ByCoordinates(0, -1, 0);
            arc a = arc.ByThreePoints(A, B, C);

            Assert.AreEqual(a.Normal, vector.ByCoordinates(0, 0, 1));
        }
        [TestMethod]
        public void ArcSweepAngle180()
        {
            point A = point.ByCoordinates(1, 0, 0);
            point B = point.ByCoordinates(0, 1, 0);
            point C = point.ByCoordinates(-1, 0, 0);
            arc a = arc.ByThreePoints(A, B, C);

            Assert.AreEqual(a.SweepAngle, Math.PI, 0.000000000000001);
        }
        [TestMethod]
        public void ArcNormal180()
        {
            point A = point.ByCoordinates(1, 0, 0);
            point B = point.ByCoordinates(0, 1, 0);
            point C = point.ByCoordinates(0, -1, 0);
            arc a = arc.ByThreePoints(A, B, C);

            Assert.AreEqual(a.Normal, vector.ByCoordinates(0, 0, 1));
        }
        [TestMethod]
        public void ArcSweepAngle271()
        {
            point A = point.ByCoordinates(1, 0, 0);
            point B = point.ByCoordinates(-1, 0, 0);
            point C = point.ByCoordinates(0, -1, 0);
            arc a = arc.ByThreePoints(A, B, C);

            Assert.AreEqual(a.SweepAngle, math.toRadians(270), 0.000000000000001);
        }
        [TestMethod]
        public void ArcNormal271()
        {
            point A = point.ByCoordinates(1, 0, 0);
            point B = point.ByCoordinates(-1, 0, 0);
            point C = point.ByCoordinates(0, -1, 0);
            arc a = arc.ByThreePoints(A, B, C);

            Assert.AreEqual(a.Normal, vector.ByCoordinates(0, 0, 1));
        }

    }
}
