<h1>BecauseWeDynamo</h1>
Library for building topology from Revit adaptive components based on three points and sets of ordering points to label and array parts for fabrication. This project includes references to:
<ul><li>DynamoCore</li>
<li>Protogeometry</li></ul>

<h2>Installing</h2>
<p>Once BecauseWeDynamo is finished building, navigate to BecauseWeDynamo/bin/ and import BecauseWeDynamo.dll into Dynamo either using library import or addding the .dll file to the root folder of Dynamo. If only interested in the necessary .dll's to import through Zero-Touch Import into Dynamo, check the Zero-Touch Node folder. These .dll's will be updated periodically and will be based on the latest build.</p>

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
