namespace DXFLibrary
{
	/// <summary>
	/// C# representation of data in a dxf file:
	/// a dxf code paired with a unit of data.
	/// </summary>
	struct Data
	{
		public int code;
		public object data;
		public Data(int code, object data)
		{
			this.code = code;
			this.data = data;
		}
	}
}