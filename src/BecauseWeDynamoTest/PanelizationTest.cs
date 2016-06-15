using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Autodesk.DesignScript.Geometry;
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
            panel p = panel.ByFaceAndParameters(Face, 0, 0, 0, 0, 0);
        }
        [TestMethod]
        [ExpectedException(typeof(ArgumentException), "ThicknessFront must be greater than zero")]
        public void CreatePanelHole()
        {
            face Face = new face();
            panelHole p = panelHole.ByFaceAndParameters(Face, 0, 0, 0, 0, 0);
        }
    }
}
