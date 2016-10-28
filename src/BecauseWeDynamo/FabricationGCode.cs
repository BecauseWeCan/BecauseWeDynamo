using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Geometry;
using Autodesk.DesignScript.Geometry;

namespace Fabrication
{
    /// <summary>
    /// list of autodesk curve objects with added gcode properties
    /// </summary>
    public class listDraw : List<Curve>
    {
        //**FIELD

        //**PROPERTIES** //**QUERY**
        /// <summary>
        /// 
        /// </summary>
        public double feedSpeed {get; set;}
        /// <summary>
        /// 
        /// </summary>
        public double toolDiameter {get; set;}
        /// <summary>
        /// 
        /// </summary>
        public double moveSpeed {get; set;}
        /// <summary>
        /// 
        /// </summary>
        public double jogSpeed {get; set;}
        /// <summary>
        /// 
        /// </summary>
        public double plungeDepth {get; set;}

        //**CONSTRUCTOR**
        internal listDraw(double feedSpeed, double toolDiameter, double moveSpeed, double jogSpeed, double plungeDepth) : base()
        {
            this.feedSpeed = feedSpeed;
            this.toolDiameter = toolDiameter;
            this.moveSpeed = moveSpeed;
            this.jogSpeed = jogSpeed;
            this.plungeDepth = plungeDepth;
        }
        internal listDraw(int capacity, double feedSpeed, double toolDiameter, double moveSpeed, double jogSpeed, double plungeDepth) : base(capacity)
        {
            this.feedSpeed = feedSpeed;
            this.toolDiameter = toolDiameter;
            this.moveSpeed = moveSpeed;
            this.jogSpeed = jogSpeed;
            this.plungeDepth = plungeDepth;
        }
        internal listDraw(IEnumerable<Curve> collection, double feedSpeed, double toolDiameter, double moveSpeed, double jogSpeed, double plungeDepth) : base(collection)
        {
            this.feedSpeed = feedSpeed;
            this.toolDiameter = toolDiameter;
            this.moveSpeed = moveSpeed;
            this.jogSpeed = jogSpeed;
            this.plungeDepth = plungeDepth;
        }

        //**STATIC METHOD**
        /// <summary>
        /// create a listDraw object;
        /// extends list class with added tooling properties, double feedSpeed, double toolDiameter, double moveSpeed, double jogSpeed, double plungeDepth
        /// </summary>
        /// <param name="feedSpeed">spped of plunge</param>
        /// <param name="toolDiameter">diameter of tool (in)</param>
        /// <param name="moveSpeed">speed of movement (ips)</param>
        /// <param name="jogSpeed">speed of jog (ips)</param>
        /// <param name="zHeight">z-height</param>
        /// <returns></returns>
        public static listDraw BySettings(double feedSpeed, double toolDiameter, double moveSpeed, double jogSpeed, double zHeight)
        {
            return new listDraw(feedSpeed, toolDiameter, moveSpeed, jogSpeed, zHeight);
        }
    }

    /// <summary>
    /// module that exports g-code as a text file
    /// </summary>
    public class gCode
    {
        //**FIELD

        //**PROPERTIES** //**QUERY**

        //**CONSTRUCTOR**
        internal gCode(string prefixFilename, string postfixFilename, listDraw[] orderedDrawList)
        {

        }

        //**METHODS** //**CREATE**
        /// <summary>
        /// 
        /// </summary>
        /// <param name="prefixFile"></param>
        /// <param name="postfixFile"></param>
        /// <param name="orderedDrawList"></param>
        /// <returns></returns>
        public static gCode ByDefault(string prefixFile, string postfixFile, listDraw[] orderedDrawList)
        { return new gCode(prefixFile, postfixFile, orderedDrawList); }

        //**METHODS** //**ACTIONS
        /// <summary>
        /// 
        /// </summary>
        /// <param name="filename"></param>
        /// <returns></returns>
        public static string[] readFile(string filename)
        {
            List<string> lines = new List<string>();
            try
            {
                //Pass the file path and file name to the StreamReader constructor
                StreamReader sr = new StreamReader(filename);

                //Read the first line of text
                string line = sr.ReadLine();

                //Continue to read until you reach end of file
                while (line != null)
                {
                    //write the lie to console window
                    lines.Add(line);
                    //Read the next line
                    line = sr.ReadLine();
                }

                //close the file
                sr.Close();
                Console.ReadLine();
            }
            catch (Exception e)
            {
                Console.WriteLine("Exception: " + e.Message);
            }
            finally
            {
                Console.WriteLine("Executing finally block.");
            }
            return lines.ToArray();
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="lines"></param>
        /// <returns></returns>
        public static bool writeFile(string filename, string[] lines)
        {
            try
            {

                //Pass the filepath and filename to the StreamWriter Constructor
                StreamWriter sw = new StreamWriter(filename);
                lines.ForEach(n => sw.WriteLine(n + "/n"));
                //Close the file
                sw.Close();
                return true;
            }
            catch (Exception e)
            {
                Console.WriteLine("Exception: " + e.Message);
            }
            finally
            {
                Console.WriteLine("Executing finally block.");
            }
            return false;
        }
    }
}
