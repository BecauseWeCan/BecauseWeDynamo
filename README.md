<h1>BecauseWeDynamo</h1>
A series of digital fabrication Nodes for Autodesk Dynamo.

It includes such things as:
<ul>
<li>Nodes for building and walking triangular mesh topology to label edges for assembly.</li>
<li>Nodes for isolating, arraying, and quasi-nesting parts for subtractive fabrication.</li>
<li>A DXF (with mathmatical curves) exporting Node.</li>
<li>A single-line-based font approprate for milling / etching for part labeling.</li>
<li>Auto-sectioning tools that slice larger solids into layers of a controllable thickness and will produce clean outline vectors of each slice.</li>
</ul>

This project includes references to:
<ul><li>DynamoCore</li>
<li>Protogeometry</li></ul>

<h2>Installing</h2>
<p>A pre-compiled Package is available within the Dynamo Package Manager for easy installtion. It may be slightly out of date. The latest version will always be posted here to GitHub.</p>
<p>If you build BecauseWeDynamo yourself, once it's done navigate to BecauseWeDynamo/bin/ and import BecauseWeDynamo.dll into Dynamo either using library import or addding the .dll file to the root folder of Dynamo. If only interested in the necessary .dll's to import through Zero-Touch Import into Dynamo, check the Zero-Touch Node folder. These .dll's will be updated periodically and will be based on the latest build. You may need to right-click on these .dll's and give them the rights to run, as Windows will restrict them by default.</p>

<h2>Known Issues</h2>
Error messages include:
- External component has thrown an exception
- Unable to get FACES range: ACIS_EXCEPTION_ACCESS_VIOLATION -- access violation
- Unable to copy the face! ACIS_EXCEPTION_ACCESS_VIOLATION -- access violation
- Unable to make sheet double-sided: CHANGE_LOST_ENT -- attempt to change record marked as deleted
- Unable to make sheet double-sided: ILLEGAL_ENTITY_MODIFICATION -- entity modification outside API BEGIN/END block

What helps to prevent errors:
- Close and reopen Revit when opening different scripts
- Interacting with GUI before running script
- Making sure that there is only one active process for Revit

<h2>License</h2>
BecauseWeDynamo: Copyright 2015 Because We Can
</br>Dynamo: Copyright 2015 Autodesk
</br>Revit: Copyright 2015 Autodesk
<p>Licensed under The MIT License (MIT); you may not use this file except in compliance with the License. You may obtain a copy of the License <a href="https://github.com/BecauseWeCan/BecauseWeDynamo/blob/master/LICENSE.md">here</a>. Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
