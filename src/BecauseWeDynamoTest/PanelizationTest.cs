using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Geometry;
using Topology;
using Panelization;

namespace BecauseWeDynamoTest
{
    [TestClass]
    public class PanelTest
    {
        [TestMethod]
        [ExpectedException(typeof(ArgumentException), "ThicknessFront must be greater than zero")]
        public void CreatePanel()
        {
            face Face = new face();
            Panel p = Panel.ByFaceAndParameters(Face, 0, 0, 0, 0, 0);
        }
    }
}
