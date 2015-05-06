using System;

namespace DXFLibrary
{
	/// <summary>
	/// Representation of ENTITIES section in dxf.
	/// </summary>
	class Entities : DXFLibrary.Section
	{
		public Entities():base("ENTITIES")
		{
			// 
			// TODO: Add constructor logic here
			//
		}
		public void AddEntity(Entity e)
		{
			this.AddElement(e);
		}
	}
}
