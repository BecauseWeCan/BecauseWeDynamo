using NUnit.Framework;
using RevitTestServices;
using RTF.Framework;

namespace Test
{
    [TestFixture]
    public class SystemTestExample : RevitSystemTestBase
    {
        [Test, TestModel(@".\Models\Test.rvt")]
        public void Location()
        {
            OpenAndRunDynamoDefinition(@".\test.dyn");

            // Your test logic goes here.
        }
    }
}