<h1>BecauseWeDynamo</h1>
A series of digital fabrication Nodes for Autodesk Dynamo / Autodesk Revit.

It includes custom nodes for:
<ul>
<li>Building and transversing Mesh Object with topological and geometric properties (Half-Edge based)</li>
<li>Labeling faces and edges with curve based font appropriate for milling / etching</li>
<li>Isolating and arraying curve profiles for subtractive fabrication</li>
<li>Exporting DXF of curve-based objects as 3D arcs and lines on designated layer and color</li>
<li>Sectioning larger solids into layers of a controllable thickness using native Dynamo geometry</li>
</ul>

This project includes references to:
<ul><li>DynamoServices</li>
<li>Protogeometry</li>
</ul>

<h2>Installing</h2>
<p>A pre-compiled Package is available within the Dynamo Package Manager for easy installtion. It may be slightly out of date. The latest version will always be posted here to GitHub.</p>
<p>If you build BecauseWeDynamo yourself, once it's done navigate to BecauseWeDynamo/bin/ and import BecauseWeDynamo.dll into Dynamo either using library import or addding the .dll file to the root folder of Dynamo. If only interested in the necessary .dll's to import through Zero-Touch Import into Dynamo, check the Zero-Touch Node folder. These .dll's will be updated periodically and will be based on the latest build. You may need to right-click on these .dll's and give them the rights to run, as Windows will restrict them by default.</p>

<h2>License</h2>
BecauseWeDynamo: Copyright 2016 Because We Can
</br>Dynamo: Copyright 2016 Autodesk
</br>Revit: Copyright 2016 Autodesk
<p>Licensed under The MIT License (MIT); you may not use this file except in compliance with the License. You may obtain a copy of the License <a href="https://github.com/BecauseWeCan/BecauseWeDynamo/blob/master/LICENSE.md">here</a>. Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
