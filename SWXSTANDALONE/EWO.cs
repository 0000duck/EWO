using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;
using SolidWorks.Interop.sldworks;

namespace SWXSTANDALONE
{
    internal class EWO
    {
        public Matrix<Double> TessalationPoints = null;
        public Matrix<Double> JoiningNormals = null;
        public Matrix<Double> BaseNormals = null;
        private ModelDoc2 swModel;
        private ModelDocExtension swModelDocExt;
        private SldWorks swApp;
        private SelectionMgr swSelectionMgr;
        private int numberOfTessalationPoints;
        private Double[] curveTessalationPoints;



        public static Vector<Double> Cross(Vector<Double> left, Vector<Double> right)
        {
            if ((left.Count != 3 || right.Count != 3))
            {
                string message = "Vectors must have a length of 3.";
                throw new Exception(message);
            }
            Vector<Double> result = new DenseVector(3);
            result[0] = left[1] * right[2] - left[2] * right[1];
            result[1] = -left[0] * right[2] + left[2] * right[0];
            result[2] = left[0] * right[1] - left[1] * right[0];

            return result;
        }

        public EWO(CosmeticWeldBeadFeatureData weldBeadData, ModelDoc2 Model, ModelDocExtension DocExt, SldWorks App)
        {
            swModel = Model;
            swModelDocExt = DocExt;
            swApp = App;
            swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
            Object[] weldmentEdgesObj = (Object[])weldBeadData.GetReferenceEdges();
            for (int edgeIndex = 0; edgeIndex < weldmentEdgesObj.Length; edgeIndex++)
            {
                Edge weldmentEdge = (Edge)weldmentEdgesObj[edgeIndex];
                Curve weldmentCurve = (Curve)weldmentEdge.GetCurve();
                Double curveStartParameter, curveEndParameter;
                bool circularCurve, periodicCurve;
                bool successfulOperation = weldmentCurve.GetEndParams(out curveStartParameter, out curveEndParameter, out circularCurve, out periodicCurve);
                Double[] curveStartPoint = (Double[])weldmentCurve.Evaluate2(curveStartParameter, 0);
                Double[] curveEndPoint = (Double[])weldmentCurve.Evaluate2(curveEndParameter, 0);
                curveTessalationPoints = (Double[])weldmentCurve.GetTessPts(0.0001, 0.0001, curveStartPoint, curveEndPoint);
                
                
                numberOfTessalationPoints = curveTessalationPoints.Length / 3;
                int columnIndex;
                if (TessalationPoints == null)
                {
                    columnIndex = 0;
                    TessalationPoints = Matrix<Double>.Build.Random(3, curveTessalationPoints.Length / 3);
                }
                else
                {
                    columnIndex = TessalationPoints.ColumnCount;
                    Matrix<Double> auxMatrix = Matrix<Double>.Build.Random(3, curveTessalationPoints.Length / 3 + columnIndex);
                    for (int auxIndex = 0; auxIndex < columnIndex; auxIndex++)
                    {
                        auxMatrix[0, auxIndex] = TessalationPoints[0, auxIndex];
                        auxMatrix[1, auxIndex] = TessalationPoints[1, auxIndex];
                        auxMatrix[2, auxIndex] = TessalationPoints[2, auxIndex];
                    }
                    TessalationPoints = auxMatrix.Clone();
                }

                for (int tessIndex = 0; tessIndex < curveTessalationPoints.Length; tessIndex += 3, columnIndex++)
                {
                    //TessalationPoints[0, columnIndex] = curveTessalationPoints[curveTessalationPoints.Length - tessIndex - 3];
                    //TessalationPoints[1, columnIndex] = curveTessalationPoints[curveTessalationPoints.Length - tessIndex - 2];
                    //TessalationPoints[2, columnIndex] = curveTessalationPoints[curveTessalationPoints.Length - tessIndex - 1];
                    TessalationPoints[0, columnIndex] = curveTessalationPoints[tessIndex];
                    TessalationPoints[1, columnIndex] = curveTessalationPoints[tessIndex + 1];
                    TessalationPoints[2, columnIndex] = curveTessalationPoints[tessIndex + 2];
                }

                Double[] edgeSelectionPoint = null, edgeTangentPoint1 = null;
                if (curveTessalationPoints.Length > 6)
                {
                    int curveMid = ((curveTessalationPoints.Length / 3) / 2) * 3;
                    edgeSelectionPoint = new Double[] { curveTessalationPoints[curveMid], curveTessalationPoints[curveMid + 1], curveTessalationPoints[curveMid + 2] };
                    edgeTangentPoint1 = new Double[] { curveTessalationPoints[curveMid - 3], curveTessalationPoints[curveMid - 2], curveTessalationPoints[curveMid - 1] };
                }
                else
                {
                    edgeSelectionPoint = new Double[] { (curveTessalationPoints[0] + curveTessalationPoints[3]) / 2, (curveTessalationPoints[1] + curveTessalationPoints[4]) / 2, (curveTessalationPoints[2] + curveTessalationPoints[5]) / 2 };
                    edgeTangentPoint1 = new Double[] { curveTessalationPoints[3], curveTessalationPoints[4], curveTessalationPoints[5] };
                }

                Vector<double> edgeTangentVector = Vector<double>.Build.DenseOfArray(new double[] { edgeSelectionPoint[0] - edgeTangentPoint1[0], edgeSelectionPoint[1] - edgeTangentPoint1[1], edgeSelectionPoint[2] - edgeTangentPoint1[2] });
                edgeTangentVector = edgeTangentVector.Normalize(1);
                Double[] edgeTangentPoint2 = edgeSelectionPoint;
                successfulOperation = swModelDocExt.SelectByID2("", "EDGE", edgeSelectionPoint[0], edgeSelectionPoint[1], edgeSelectionPoint[2], false, 0, null, 0);
                weldmentEdge = (Edge)swSelectionMgr.GetSelectedObject6(1, 0);
                Face2 edgeFace1, edgeFace2, joiningFace = null, baseFace = null;
                weldmentEdge.IGetTwoAdjacentFaces2(out edgeFace1, out edgeFace2);
                int faceCookies = swApp.RegisterTrackingDefinition("WeldExporter");
                edgeFace1.SetTrackingID(faceCookies, 1);
                edgeFace2.SetTrackingID(faceCookies, 2);
                object trackCookiesObject = null;
                baseFace = getBaseFace(edgeFace1, edgeFace2, edgeTangentVector, edgeSelectionPoint, faceCookies);
                float[] baseFaceTessalationNorms = (float[])baseFace.GetTessNorms();
                float[] edgeFace1TessalationNorms = (float[])edgeFace1.GetTessNorms();
                float[] edgeFace2TessalationNorms = (float[])edgeFace2.GetTessNorms();
                int tessNo = baseFaceTessalationNorms.Length;
                bool term1 = Math.Round(baseFaceTessalationNorms[0], 4) == Math.Round(-edgeFace1TessalationNorms[0], 4);
                bool term2 = Math.Round(baseFaceTessalationNorms[1], 4) == Math.Round(-edgeFace1TessalationNorms[1], 4);
                bool term3 = Math.Round(baseFaceTessalationNorms[2], 4) == Math.Round(-edgeFace1TessalationNorms[2], 4);
                bool term4 = Math.Round(baseFaceTessalationNorms[tessNo - 3], 4) == Math.Round(-edgeFace1TessalationNorms[0], 4);
                bool term5 = Math.Round(baseFaceTessalationNorms[tessNo - 2], 4) == Math.Round(-edgeFace1TessalationNorms[1], 4);
                bool term6 = Math.Round(baseFaceTessalationNorms[tessNo - 1], 4) == Math.Round(-edgeFace1TessalationNorms[2], 4);
                bool term7 = Math.Round(baseFaceTessalationNorms[0], 4) == Math.Round(-edgeFace2TessalationNorms[0], 4);
                bool term8 = Math.Round(baseFaceTessalationNorms[1], 4) == Math.Round(-edgeFace2TessalationNorms[1], 4);
                bool term9 = Math.Round(baseFaceTessalationNorms[2], 4) == Math.Round(-edgeFace2TessalationNorms[2], 4);
                bool term10 = Math.Round(baseFaceTessalationNorms[tessNo - 3], 4) == Math.Round(-edgeFace2TessalationNorms[0], 4);
                bool term11 = Math.Round(baseFaceTessalationNorms[tessNo - 2], 4) == Math.Round(-edgeFace2TessalationNorms[1], 4);
                bool term12 = Math.Round(baseFaceTessalationNorms[tessNo - 1], 4) == Math.Round(-edgeFace2TessalationNorms[2], 4);
                if ((term1 && term2 && term3) || (term4 && term5 && term6))
                {
                    joiningFace = edgeFace2;
                } else if ((term7 && term8 && term9) || (term10 && term11 && term12))
                {
                    joiningFace = edgeFace1;
                }

                if (JoiningNormals == null)
                {
                    JoiningNormals = getDiscretizationOfNormals(joiningFace);
                }
                else
                {
                    JoiningNormals = appendMatrix(JoiningNormals, getDiscretizationOfNormals(joiningFace));
                }

                if (BaseNormals == null)
                {
                    BaseNormals = getDiscretizationOfNormals(baseFace);
                }
                else
                {
                    BaseNormals = appendMatrix(BaseNormals, getDiscretizationOfNormals(baseFace));
                }

            }
        }

        public EWO(Face2[] fromFaces, Face2[] toFaces, CosmeticWeldBeadFeatureData weldBeadData, ModelDoc2 Model, ModelDocExtension DocExt, SldWorks App)
        {
            swModel = Model;
            swModelDocExt = DocExt;
            swApp = App;
            swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
            Object[] weldmentEdgesObj = (Object[])weldBeadData.GetReferenceEdges();
            for (int edgeIndex = 0; edgeIndex < weldmentEdgesObj.Length; edgeIndex++)
            {               
                Edge weldmentEdge = (Edge)weldmentEdgesObj[edgeIndex];
                Curve weldmentCurve = (Curve)weldmentEdge.GetCurve();
                Double curveStartParameter, curveEndParameter;
                bool circularCurve, periodicCurve;
                bool successfulOperation = weldmentCurve.GetEndParams(out curveStartParameter, out curveEndParameter, out circularCurve, out periodicCurve);
                Double[] curveStartPoint = (Double[])weldmentCurve.Evaluate2(curveStartParameter, 0);
                Double[] curveEndPoint = (Double[])weldmentCurve.Evaluate2(curveEndParameter, 0);
                curveTessalationPoints = (Double[])weldmentCurve.GetTessPts(0.0001, 0.0001, curveStartPoint, curveEndPoint);
                int curveMid = ((curveTessalationPoints.Length / 3) / 2) * 3;
                Double[] edgeSelectionPoint = { curveTessalationPoints[curveMid], curveTessalationPoints[curveMid + 1], curveTessalationPoints[curveMid + 2] };
                Double[] edgeTangentPoint1 = { curveTessalationPoints[curveMid - 3], curveTessalationPoints[curveMid - 2], curveTessalationPoints[curveMid - 1] };
                numberOfTessalationPoints = curveTessalationPoints.Length / 3;
                int columnIndex;
                if (TessalationPoints == null)
                {
                    columnIndex = 0;
                    TessalationPoints = Matrix<Double>.Build.Random(3, curveTessalationPoints.Length / 3);
                }
                else
                {
                    columnIndex = TessalationPoints.ColumnCount;
                    Matrix<Double> auxMatrix = Matrix<Double>.Build.Random(3, curveTessalationPoints.Length / 3 + columnIndex);
                    for (int auxIndex = 0; auxIndex < columnIndex; auxIndex++)
                    {
                        auxMatrix[0, auxIndex] = TessalationPoints[0, auxIndex];
                        auxMatrix[1, auxIndex] = TessalationPoints[1, auxIndex];
                        auxMatrix[2, auxIndex] = TessalationPoints[2, auxIndex];
                    }
                    TessalationPoints = auxMatrix.Clone();
                }
                
                for (int tessIndex = 0; tessIndex < curveTessalationPoints.Length; tessIndex+=3, columnIndex++)
                {
                    //TessalationPoints[0, columnIndex] = curveTessalationPoints[curveTessalationPoints.Length - tessIndex - 3];
                    //TessalationPoints[1, columnIndex] = curveTessalationPoints[curveTessalationPoints.Length - tessIndex - 2];
                    //TessalationPoints[2, columnIndex] = curveTessalationPoints[curveTessalationPoints.Length - tessIndex - 1];
                    TessalationPoints[0, columnIndex] = curveTessalationPoints[tessIndex ];
                    TessalationPoints[1, columnIndex] = curveTessalationPoints[tessIndex + 1];
                    TessalationPoints[2, columnIndex] = curveTessalationPoints[tessIndex + 2];
                }

                Vector<double> edgeTangentVector = Vector<double>.Build.DenseOfArray(new double[] { edgeSelectionPoint[0] - edgeTangentPoint1[0], edgeSelectionPoint[1] - edgeTangentPoint1[1], edgeSelectionPoint[2] - edgeTangentPoint1[2] });
                edgeTangentVector = edgeTangentVector.Normalize(1);
                Double[] edgeTangentPoint2 = edgeSelectionPoint;
                successfulOperation = swModelDocExt.SelectByID2("", "EDGE", edgeSelectionPoint[0], edgeSelectionPoint[1], edgeSelectionPoint[2], false, 0, null, 0);
                weldmentEdge = (Edge)swSelectionMgr.GetSelectedObject6(1, 0);
                Face2 edgeFace1, edgeFace2, joiningFace = null, baseFace = null;
                weldmentEdge.IGetTwoAdjacentFaces2(out edgeFace1, out edgeFace2);
                int faceCookies = swApp.RegisterTrackingDefinition("WeldExporter");
                edgeFace1.SetTrackingID(faceCookies, 1);
                edgeFace2.SetTrackingID(faceCookies, 2);
                object trackCookiesObject = null;
                for (int faceIndex = 0; faceIndex < toFaces.Length; faceIndex++)
                {
                    toFaces[faceIndex].GetTrackingIDs(faceCookies, out trackCookiesObject);
                    if (trackCookiesObject != null)
                    {
                        int[] trackCookies = (int[])trackCookiesObject;
                        if (trackCookies[0] == 1)
                        {
                            joiningFace = edgeFace1;
                            baseFace = getBaseFace(joiningFace, edgeTangentVector, edgeSelectionPoint, faceCookies);
                        }
                        else
                        {
                            joiningFace = edgeFace2;
                            baseFace = getBaseFace(joiningFace, edgeTangentVector, edgeSelectionPoint, faceCookies);
                        }
                        edgeFace1.RemoveTrackingID(faceCookies);
                        edgeFace2.RemoveTrackingID(faceCookies);

                        break;
                    }

                }
                if (joiningFace == null && baseFace == null)
                {
                    for (int faceIndex = 0; faceIndex < fromFaces.Length; faceIndex++)
                    {
                        fromFaces[faceIndex].GetTrackingIDs(faceCookies, out trackCookiesObject);
                        if (trackCookiesObject != null)
                        {
                            int[] trackCookies = (int[])trackCookiesObject;
                            if (trackCookies[0] == 1)
                            {
                                joiningFace = edgeFace1;
                                baseFace = getBaseFace(joiningFace, edgeTangentVector, edgeSelectionPoint, faceCookies);
                            }
                            else
                            {
                                joiningFace = edgeFace2;
                                baseFace = getBaseFace(joiningFace, edgeTangentVector, edgeSelectionPoint, faceCookies);
                            }
                            edgeFace1.RemoveTrackingID(faceCookies);
                            edgeFace2.RemoveTrackingID(faceCookies);

                            break;
                        }

                    }
                }
                if (JoiningNormals == null)
                {
                    JoiningNormals = getDiscretizationOfNormals(joiningFace);
                }
                else
                {
                    JoiningNormals = appendMatrix(JoiningNormals, getDiscretizationOfNormals(joiningFace));                    
                }

                if (BaseNormals == null)
                {
                    BaseNormals = getDiscretizationOfNormals(baseFace);
                }
                else
                {
                    BaseNormals = appendMatrix(BaseNormals, getDiscretizationOfNormals(baseFace));
                }
            }
        }

        private Matrix<Double> appendMatrix(Matrix<Double> originalMatrix, Matrix<Double> appendedMatrix)
        {
            Matrix<Double> newMatrix = Matrix<Double>.Build.Dense(3, originalMatrix.ColumnCount + appendedMatrix.ColumnCount);
            for (int columnIndex = 0; columnIndex < originalMatrix.ColumnCount + appendedMatrix.ColumnCount; columnIndex++)
            {
                if (columnIndex < originalMatrix.ColumnCount)
                {
                    newMatrix[0, columnIndex] = originalMatrix[0, columnIndex];
                    newMatrix[1, columnIndex] = originalMatrix[1, columnIndex];
                    newMatrix[2, columnIndex] = originalMatrix[2, columnIndex];
                }
                else
                {
                    newMatrix[0, columnIndex] = appendedMatrix[0, columnIndex - originalMatrix.ColumnCount];
                    newMatrix[1, columnIndex] = appendedMatrix[1, columnIndex - originalMatrix.ColumnCount];
                    newMatrix[2, columnIndex] = appendedMatrix[2, columnIndex - originalMatrix.ColumnCount];
                }
            }

            return newMatrix;
            
        }

        private Matrix<Double> getDiscretizationOfNormals(Face2 face)
        {

            Matrix<Double> faceDiscreteNormals = Matrix<Double>.Build.Dense(3, numberOfTessalationPoints);
            MathTransform faceTranform = getTransform(face);
            float[] faceTessalationNorms = (float[])face.GetTessNorms();
            int midIndex, endIndex = faceTessalationNorms.Length;
            Vector<Double> faceFirstNormal = computeTransform(faceTranform, Vector<double>.Build.DenseOfArray(new double[] { faceTessalationNorms[0], faceTessalationNorms[1], faceTessalationNorms[2] }));
            midIndex = ((endIndex / 3) / 2) * 3;
            Vector<Double> faceMidNormal = computeTransform(faceTranform, Vector<double>.Build.DenseOfArray(new double[] { faceTessalationNorms[midIndex], faceTessalationNorms[midIndex + 1], faceTessalationNorms[midIndex + 2] }));
            Vector<Double> faceLastNormal = computeTransform(faceTranform, Vector<double>.Build.DenseOfArray(new double[] { faceTessalationNorms[endIndex - 3], faceTessalationNorms[endIndex - 2], faceTessalationNorms[endIndex - 1] }));
            Vector<Double> DifNormal = (faceFirstNormal - faceMidNormal) / (numberOfTessalationPoints / 2);
            //int sign_x = Math.Sign(Math.Round(curveTessalationPoints[(curveTessalationPoints.Length / 3 / 2) * 3] - curveTessalationPoints[0], 6));
            //int sign_y = Math.Sign(Math.Round(curveTessalationPoints[(curveTessalationPoints.Length / 3 / 2) * 3 + 1] - curveTessalationPoints[1], 6));
            //int sign_z = Math.Sign(Math.Round(curveTessalationPoints[(curveTessalationPoints.Length / 3 / 2) * 3 + 2] - curveTessalationPoints[2], 6));
            //int normal_x = Math.Sign(faceMidNormal[0] - faceFirstNormal[0]);
            //int normal_y = Math.Sign(faceMidNormal[1] - faceFirstNormal[1]);
            //int normal_z = Math.Sign(faceMidNormal[2] - faceFirstNormal[2]);
            for (int i = 0; i < 3; i++)
            {
                if (faceFirstNormal[i] > faceMidNormal[i])
                    DifNormal[i] = -Math.Abs(DifNormal[i]);    
                else
                    DifNormal[i] = Math.Abs(DifNormal[i]);

            }
            faceDiscreteNormals[0, 0] = faceFirstNormal[0];
            faceDiscreteNormals[1, 0] = faceFirstNormal[1];
            faceDiscreteNormals[2, 0] = faceFirstNormal[2];

            for (int indexTess = 1; indexTess < numberOfTessalationPoints / 2; indexTess++)
            {
                faceDiscreteNormals[0, indexTess] = (faceDiscreteNormals[0, indexTess - 1] + DifNormal[0]);
                faceDiscreteNormals[1, indexTess] = (faceDiscreteNormals[1, indexTess - 1] + DifNormal[1]);
                faceDiscreteNormals[2, indexTess] = (faceDiscreteNormals[2, indexTess - 1] + DifNormal[2]);

            }

            int signX = Math.Sign(Math.Round(curveTessalationPoints[3] - curveTessalationPoints[0], 6));
            int signY = Math.Sign(Math.Round(curveTessalationPoints[4] - curveTessalationPoints[1], 6));
            int signZ = Math.Sign(Math.Round(curveTessalationPoints[5] - curveTessalationPoints[2], 6));
            int signNormalX = Math.Sign(faceDiscreteNormals[0, 0] - faceDiscreteNormals[0, 1]);
            int signNormalY = Math.Sign(faceDiscreteNormals[1, 0] - faceDiscreteNormals[1, 1]);
            int signNormalZ = Math.Sign(faceDiscreteNormals[2, 0] - faceDiscreteNormals[2, 1]);
            int signNormal = 0;
            int signPosition = 0;
            if (signNormalX != 0 || signNormalY != 0 || signNormalZ != 0)
            {
                if (signNormalX == 0)
                    signNormal = signNormalY * signNormalZ;
                else if (signNormalY == 0)
                    signNormal = signNormalX * signNormalZ;
                else if (signNormalZ == 0)
                    signNormal = signNormalX * signNormalY;
                else
                    signNormal = signNormalX * signNormalY * signNormalZ;
            }
            if (signX != 0 || signY != 0 || signZ != 0)
            {
                if (signX == 0)
                    signPosition = signY * signZ;
                else if (signNormalY == 0)
                    signPosition = signX * signZ;
                else if (signNormalZ == 0)
                    signPosition = signX * signY;
                else
                    signPosition = signX * signY * signZ;
            }
            int signTotal = signNormal * signPosition;


            DifNormal = (faceMidNormal - faceLastNormal) / ((numberOfTessalationPoints) / 2);
            for (int i = 0; i < 3; i++)
            {
                if (faceMidNormal[i] > faceLastNormal[i])
                    DifNormal[i] = -Math.Abs(DifNormal[i]);
                else
                    DifNormal[i] = Math.Abs(DifNormal[i]);

            }
            faceDiscreteNormals[0, numberOfTessalationPoints / 2] = faceMidNormal[0];
            faceDiscreteNormals[1, numberOfTessalationPoints / 2] = faceMidNormal[1];
            faceDiscreteNormals[2, numberOfTessalationPoints / 2] = faceMidNormal[2];

            for (int indexTess = numberOfTessalationPoints / 2 + 1; indexTess < numberOfTessalationPoints; indexTess++)
            {
                faceDiscreteNormals[0, indexTess] = (faceDiscreteNormals[0, indexTess - 1] + DifNormal[0]);
                faceDiscreteNormals[1, indexTess] = (faceDiscreteNormals[1, indexTess - 1] + DifNormal[1]);
                faceDiscreteNormals[2, indexTess] = (faceDiscreteNormals[2, indexTess - 1] + DifNormal[2]);

            }

            if (signTotal < 0)
            {
                Matrix<Double> auxMatrix = faceDiscreteNormals.Clone();
                for (int i = 0; i < faceDiscreteNormals.ColumnCount; i++)
                {
                    faceDiscreteNormals[0, i] = auxMatrix[0, auxMatrix.ColumnCount - 1 - i];
                    faceDiscreteNormals[1, i] = auxMatrix[1, auxMatrix.ColumnCount - 1 - i];
                    faceDiscreteNormals[2, i] = auxMatrix[2, auxMatrix.ColumnCount - 1 - i];
                }
            }

            return faceDiscreteNormals;

        }

        private static Vector<Double> computeTransform (MathTransform normalTransform, Vector<Double> normalVector)
        {
            Double[] joiningTransformArray = (Double[])normalTransform.ArrayData;
            Double[] joiningNormalDouble = { joiningTransformArray[0] * normalVector[0] + joiningTransformArray[3] * normalVector[1] + joiningTransformArray[6] * normalVector[2],
                                       joiningTransformArray[1] * normalVector[0] + joiningTransformArray[4] * normalVector[1] + joiningTransformArray[7] * normalVector[2],
                                       joiningTransformArray[2] * normalVector[0] + joiningTransformArray[5] * normalVector[1] + joiningTransformArray[8] * normalVector[2]};
            Vector<Double> transformedNormal = Vector<double>.Build.DenseOfArray(new double[] { joiningNormalDouble[0], joiningNormalDouble[1], joiningNormalDouble[2] });
            return transformedNormal;
        }

        private MathTransform getTransform (Face2 face)
        {
            Feature joiningFeature = (Feature)face.GetFeature();
            string featureName = joiningFeature.GetNameForSelection(out featureName);
            int pFrom = featureName.IndexOf("@") + "@".Length;
            int pTo = featureName.LastIndexOf("@");
            String componentName = featureName.Substring(pFrom, pTo - pFrom);
            bool feature_successfuly_selected = swModelDocExt.SelectByID2(componentName, "COMPONENT", 0, 0, 0, false, 0, null, 0);
            Component2 comp = (Component2)swSelectionMgr.GetSelectedObject6(1, 0);
            comp.DeSelect();
            return comp.Transform;
            
        }

        private Face2 getBaseFace (Face2 joiningFace, Vector<Double> edgeTangentVector, Double[] edgeSelectionPoint, int faceCookies)
        {
            object trackCookiesObject = null;
            Face2 baseFace = null;

            MathTransform joiningTransform = getTransform(joiningFace);

            float[] joiningTessalationNorms = (float[])joiningFace.GetTessNorms();
            int index = ((joiningTessalationNorms.Length / 3) / 2) * 3;
            Vector<Double> joiningNormal = computeTransform(joiningTransform, Vector<double>.Build.DenseOfArray(new double[] { joiningTessalationNorms[index], joiningTessalationNorms[index + 1], joiningTessalationNorms[index + 2] }));

            
            
            Vector<Double> parallel = (joiningNormal.DotProduct(edgeTangentVector) / edgeTangentVector.DotProduct(edgeTangentVector)) * edgeTangentVector;
            Vector<Double> perpendicular = joiningNormal - parallel;
            Vector<Double> w = Cross(edgeTangentVector, perpendicular);

            for (Double indexPi = 0; indexPi < Math.PI * 2; indexPi += Math.PI / 180)
            {
                Double x1 = Math.Cos(indexPi) / perpendicular.Norm(1);
                Double x2 = Math.Sin(indexPi) / w.Norm(1);
                Vector<Double> apb = perpendicular.Norm(1) * (x1 * perpendicular + x2 * w);
                Vector<Double> result = (apb + parallel) * 0.001;
                Double[] faceSelectionPoint = { result[0] + edgeSelectionPoint[0], result[1] + edgeSelectionPoint[1], result[2] + edgeSelectionPoint[2] };
                bool successfulOperation = swModelDocExt.SelectByID2("", "FACE", faceSelectionPoint[0], faceSelectionPoint[1], faceSelectionPoint[2], false, 0, null, 0);
                if (successfulOperation)
                {
                    Face2 face = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                    face.GetTrackingIDs(faceCookies, out trackCookiesObject);
                    if (trackCookiesObject == null)
                    {
                        baseFace = face;
                        break;
                    }
                }
            }

            return baseFace;
        }

        private Face2 getBaseFace (Face2 edgeFace1, Face2 edgeFace2, Vector<Double> edgeTangentVector, Double[] edgeSelectionPoint, int faceCookies)
        {
            object trackCookiesObject = null;
            Face2 baseFace = null;

            MathTransform joiningTransform = getTransform(edgeFace1);

            float[] joiningTessalationNorms = (float[])edgeFace1.GetTessNorms();
            int index = ((joiningTessalationNorms.Length / 3) / 2) * 3;
            Vector<Double> joiningNormal = computeTransform(joiningTransform, Vector<double>.Build.DenseOfArray(new double[] { joiningTessalationNorms[index], joiningTessalationNorms[index + 1], joiningTessalationNorms[index + 2] }));



            Vector<Double> parallel = (joiningNormal.DotProduct(edgeTangentVector) / edgeTangentVector.DotProduct(edgeTangentVector)) * edgeTangentVector;
            Vector<Double> perpendicular = joiningNormal - parallel;
            Vector<Double> w = Cross(edgeTangentVector, perpendicular);

            for (Double indexPi = 0; indexPi < Math.PI * 2; indexPi += Math.PI / 180)
            {
                Double x1 = Math.Cos(indexPi) / perpendicular.Norm(1);
                Double x2 = Math.Sin(indexPi) / w.Norm(1);
                Vector<Double> apb = perpendicular.Norm(1) * (x1 * perpendicular + x2 * w);
                Vector<Double> result = (apb + parallel) * 0.001;
                Double[] faceSelectionPoint = { result[0] + edgeSelectionPoint[0], result[1] + edgeSelectionPoint[1], result[2] + edgeSelectionPoint[2] };
                bool successfulOperation = swModelDocExt.SelectByID2("", "FACE", faceSelectionPoint[0], faceSelectionPoint[1], faceSelectionPoint[2], false, 0, null, 0);
                if (successfulOperation)
                {
                    Face2 face = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                    face.GetTrackingIDs(faceCookies, out trackCookiesObject);
                    if (trackCookiesObject == null)
                    {
                        baseFace = face;
                        break;
                    }
                }
            }

            if (baseFace == null)
            {
                joiningTransform = getTransform(edgeFace2);

                joiningTessalationNorms = (float[])edgeFace2.GetTessNorms();
                index = ((joiningTessalationNorms.Length / 3) / 2) * 3;
                joiningNormal = computeTransform(joiningTransform, Vector<double>.Build.DenseOfArray(new double[] { joiningTessalationNorms[index], joiningTessalationNorms[index + 1], joiningTessalationNorms[index + 2] }));



                parallel = (joiningNormal.DotProduct(edgeTangentVector) / edgeTangentVector.DotProduct(edgeTangentVector)) * edgeTangentVector;
                perpendicular = joiningNormal - parallel;
                w = Cross(edgeTangentVector, perpendicular);

                for (Double indexPi = 0; indexPi < Math.PI * 2; indexPi += Math.PI / 180)
                {
                    Double x1 = Math.Cos(indexPi) / perpendicular.Norm(1);
                    Double x2 = Math.Sin(indexPi) / w.Norm(1);
                    Vector<Double> apb = perpendicular.Norm(1) * (x1 * perpendicular + x2 * w);
                    Vector<Double> result = (apb + parallel) * 0.001;
                    Double[] faceSelectionPoint = { result[0] + edgeSelectionPoint[0], result[1] + edgeSelectionPoint[1], result[2] + edgeSelectionPoint[2] };
                    bool successfulOperation = swModelDocExt.SelectByID2("", "FACE", faceSelectionPoint[0], faceSelectionPoint[1], faceSelectionPoint[2], false, 0, null, 0);
                    if (successfulOperation)
                    {
                        Face2 face = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                        face.GetTrackingIDs(faceCookies, out trackCookiesObject);
                        if (trackCookiesObject == null)
                        {
                            baseFace = face;
                            break;
                        }
                    }
                }
            }

            return baseFace;
        }
    }
}
