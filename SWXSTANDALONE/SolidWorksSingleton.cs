using SolidWorks.Interop.sldworks;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Threading;
using SolidWorks.Interop.swconst;
using System.IO;
using System.Numerics;
using System.Xml;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;






namespace SWXSTANDALONE
{
    class SolidWorksSingleton
    {
        public static SldWorks swApp;

        static Double computeAngle(Double x1, Double y1, Double z1, Double x2, Double y2, Double z2)
        {
            Double xx = x1 * x2 + y1 * y2 + z1 * z2;
            Double yy = Math.Sqrt(x1 * x1 + y1 * y1 + z1 * z1);
            Double yyy = Math.Sqrt(x2 * x2 + y2 * y2 + z2 * z2);
            Double angle = Math.Acos(xx / yy * yyy);
            Double angle_degrees = (180 / Math.PI) * angle;
            return angle_degrees;
        }

        static Matrix<Double> AbsoluteOrientation(Matrix<Double> B, Matrix<Double> A)
        {
            int col_count = A.ColumnCount;
            Double lc_x = 0, lc_y = 0, lc_z = 0, rc_x = 0, rc_y = 0, rc_z = 0;
            for (int index = 0; index < col_count; index++)
            {
                lc_x = lc_x + A[0, index];
                lc_y = lc_y + A[1, index];
                lc_z = lc_z + A[2, index];
                rc_x = rc_x + B[0, index];
                rc_y = rc_y + B[1, index];
                rc_z = rc_z + B[2, index];
            }
            lc_x = lc_x / col_count;
            lc_y = lc_y / col_count;
            lc_z = lc_z / col_count;
            rc_x = rc_x / col_count;
            rc_y = rc_y / col_count;
            rc_z = rc_z / col_count;

            Matrix<double> lc = DenseMatrix.OfArray(new double[,] { { lc_x }, { lc_y }, { lc_z } });
            Matrix<double> rc = DenseMatrix.OfArray(new double[,] { { rc_x }, { rc_y }, { rc_z } });

            for (int index = 0; index < col_count; index++)
            {
                A[0, index] = A[0, index] - lc_x;
                A[1, index] = A[1, index] - lc_y;
                A[2, index] = A[2, index] - lc_z;
                B[0, index] = B[0, index] - rc_x;
                B[1, index] = B[1, index] - rc_y;
                B[2, index] = B[2, index] - rc_z;
            }



            Matrix<double> M = A.Multiply(B.Transpose());
            Double Sxx = M[0, 0];
            Double Syx = M[1, 0];
            Double Szx = M[2, 0];
            Double Sxy = M[0, 1];
            Double Syy = M[1, 1];
            Double Szy = M[2, 1];
            Double Sxz = M[0, 2];
            Double Syz = M[1, 2];
            Double Szz = M[2, 2];

            Matrix<Double> N = DenseMatrix.OfArray(new double[,] { { Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx }, { Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz }, { Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy }, { Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz } });


            Evd<double> eigen = N.Evd();
            MathNet.Numerics.LinearAlgebra.Vector<Complex> D = eigen.EigenValues;
            Matrix<Double> V = eigen.EigenVectors;
            

            int max_index = 3;
            Double max_eigen = D[0].Magnitude;
            for (int index = 1; index < D.Count; index++)
            {
                if (D[index].Magnitude > max_eigen)
                {
                    max_eigen = D[index].Magnitude;
                    max_index = index;
                }
            }

            Double max_val = Math.Abs(V[0, max_index]);
            int max_index2 = 0;
            MathNet.Numerics.LinearAlgebra.Vector<Double> q = DenseVector.OfArray(new double[] { 0, 0, 0, 0 });
            q[0] = V[0, max_index];
            for (int index = 1; index < V.RowCount; index++)
            {
                q[index] = V[index, max_index];
                if (Math.Abs(V[index, max_index]) > max_val)
                {
                    max_val = Math.Abs(V[index, max_index]);
                    max_index2 = index;
                }
            }
            Double sgn = 0;
            if (V[max_index2, max_index] > 0)
                sgn = 1;
            else
                sgn = -1;

            q = q * sgn;
            Double nrm = 0;
            nrm = q.L2Norm();
            q = q / nrm;
            Matrix<Double> Z = DenseMatrix.OfArray(new double[,] { { q[0], -q[3], q[2] }, { q[3], q[0], -q[1] }, { -q[2], q[1], q[0] } });
            Matrix<Double> qT = DenseMatrix.OfArray(new double[,] { { q[1] }, { q[2] }, { q[3] } });
            Matrix<Double> qR = DenseMatrix.OfArray(new double[,] { { q[1], q[2], q[3] } });
            Matrix<Double> w = qT.Multiply(qR);
            Matrix<Double> R = w + Z.Multiply(Z);

            //R = R.Transpose();
            return R;
        }

        static Vector3 getNormal(Face2 frFace, ModelDocExtension swModelDocExt, ModelDoc2 swModel)
        {
            bool feature_successfuly_selected = false;
            SelectionMgr swSelectionMgr = default(SelectionMgr);
            String name;

            Double[] face2_normal = (Double[])frFace.Normal;

            Feature ft = (Feature)frFace.GetFeature();
            name = ft.GetNameForSelection(out name);
            int pFrom = name.IndexOf("@") + "@".Length;
            int pTo = name.LastIndexOf("@");
            String result = name.Substring(pFrom, pTo - pFrom);

            feature_successfuly_selected = swModelDocExt.SelectByID2(result, "COMPONENT", 0, 0, 0, false, 0, null, 0);
            swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
            Component2 comp = (Component2)swSelectionMgr.GetSelectedObject6(1, 0);

            MathTransform trsf = comp.Transform2;
            Double[] arr = (Double[])trsf.ArrayData;

            int edges_count = frFace.GetEdgeCount();
            Object[] frFaceEdges = (Object[])frFace.GetEdges();
            Edge edge1 = (Edge)frFaceEdges[0];
            Edge edge2 = (Edge)frFaceEdges[2];
            Vertex vertex1 = (Vertex)edge2.GetStartVertex();
            Vertex vertex2 = (Vertex)edge2.GetEndVertex();
            Vertex vertex3 = (Vertex)edge1.GetEndVertex();

            Double[] point1 = (Double[])vertex1.GetPoint();
            Double[] point2 = (Double[])vertex2.GetPoint();
            Double[] point3 = (Double[])vertex3.GetPoint();


            Double x1 = arr[0] * point1[0] + Math.Abs(arr[1] * point1[1]) + arr[2] * point1[2] + arr[9] * 1;
            Double y1 = arr[3] * point1[0] + arr[4] * point1[1] + Math.Abs(arr[5] * point1[2]) + arr[10] * 1;
            Double z1 = arr[6] * point1[0] + arr[7] * point1[1] + arr[8] * point1[2] + arr[11] * 1;

            Double x2 = arr[0] * point2[0] + arr[1] * point2[1] + arr[2] * point2[2] + arr[9] * 1;
            Double y2 = arr[3] * point2[0] + arr[4] * point2[1] + arr[5] * point2[2] + arr[10] * 1;
            Double z2 = arr[6] * point2[0] + arr[7] * point2[1] + arr[8] * point2[2] + arr[11] * 1;

            Double x3 = arr[0] * point3[0] + arr[1] * point3[1] + arr[2] * point3[2] + arr[9] * 1;
            Double y3 = arr[3] * point3[0] + arr[4] * point3[1] + arr[5] * point3[2] + arr[10] * 1;
            Double z3 = arr[6] * point3[0] + arr[7] * point3[1] + arr[8] * point3[2] + arr[11] * 1;

            Double x4 = arr[0] * face2_normal[0] + arr[1] * face2_normal[1] + arr[2] * face2_normal[2];
            Double y4 = arr[3] * face2_normal[0] + arr[4] * face2_normal[1] + arr[5] * face2_normal[2];
            Double z4 = arr[6] * face2_normal[0] + arr[7] * face2_normal[1] + arr[8] * face2_normal[2];

            Double[] r1 = { x1 - x2, y1 - y2, z1 - z2 };
            Double[] r2 = { x1 - x3, y1 - y3, z1 - z3 };
            Vector3 v1 = new Vector3(Convert.ToSingle(r1[0]), Convert.ToSingle(r1[1]), Convert.ToSingle(r1[2]));
            Vector3 v2 = new Vector3(Convert.ToSingle(r2[0]), Convert.ToSingle(r2[1]), Convert.ToSingle(r2[2]));
            Vector3 normalvv = Vector3.Cross(v1, v2);
            return normalvv;

        }

        static Vector3 getNormal2(Face2 frFace, ModelDocExtension swModelDocExt, ModelDoc2 swModel, ref MathTransform componentTransform)
        {
            bool feature_successfuly_selected = false;
            SelectionMgr swSelectionMgr = default(SelectionMgr);
            String name;

            Double[] face2_normal = (Double[])frFace.Normal;

            Feature ft = (Feature)frFace.GetFeature();
            name = ft.GetNameForSelection(out name);
            int pFrom = name.IndexOf("@") + "@".Length;
            int pTo = name.LastIndexOf("@");
            String result = name.Substring(pFrom, pTo - pFrom);

            feature_successfuly_selected = swModelDocExt.SelectByID2(result, "COMPONENT", 0, 0, 0, false, 0, null, 0);
            swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
            Component2 comp = (Component2)swSelectionMgr.GetSelectedObject6(1, 0);

            componentTransform = comp.Transform2;
            Double[] arr = (Double[])componentTransform.ArrayData;

            Double x4 = arr[0] * face2_normal[0] + arr[1] * face2_normal[1] + arr[2] * face2_normal[2];
            Double y4 = arr[3] * face2_normal[0] + arr[4] * face2_normal[1] + arr[5] * face2_normal[2];
            Double z4 = arr[6] * face2_normal[0] + arr[7] * face2_normal[1] + arr[8] * face2_normal[2];
            comp.DeSelect();
            Vector3 normalvv = new Vector3(Convert.ToSingle(x4), Convert.ToSingle(y4), Convert.ToSingle(z4));
            return normalvv;

        }

        static Vector3 computeTransform(MathTransform trsf, Double normal_x, Double normal_y, Double normal_z)
        {
            Double[] arr = (Double[])trsf.ArrayData;
            Double x4 = arr[0] * normal_x + arr[1] * normal_y + arr[2] * normal_z;
            Double y4 = arr[3] * normal_x + arr[4] * normal_y + arr[5] * normal_z;
            Double z4 = arr[6] * normal_x + arr[7] * normal_y + arr[8] * normal_z;
            Vector3 normalvv = new Vector3(Convert.ToSingle(x4), Convert.ToSingle(y4), Convert.ToSingle(z4));
            return normalvv;
        }

        static MathNet.Numerics.LinearAlgebra.Vector<Double> computeTransform2(MathTransform trsf, Double normal_x, Double normal_y, Double normal_z)
        {
            Double[] arr = (Double[])trsf.ArrayData;            
            return DenseVector.OfArray(new double[] { arr[0] * normal_x + arr[3] * normal_y + arr[6] * normal_z,
                                                      arr[1] * normal_x + arr[4] * normal_y + arr[7] * normal_z,
                                                      arr[2] * normal_x + arr[5] * normal_y + arr[8] * normal_z});
        }

        static Double Distance(MathNet.Numerics.LinearAlgebra.Vector<Double> A, MathNet.Numerics.LinearAlgebra.Vector<Double> B)
        {
            return Math.Sqrt((B[0] - A[0])*(B[0] - A[0]) + (B[1] - A[1])*(B[1] - A[1]) + (B[2] - A[2])* (B[2] - A[2]));
        }


        internal static void Holla(int file_index)
        {


            
            ModelDoc2 swModel;
            ModelDocExtension swModelDocExt;
            SelectionMgr swSelectionMgr = default(SelectionMgr);
            WeldmentBeadFeatureData swWeldBead;
            Feature swFeat;
            CosmeticWeldBeadFeatureData swCosmeticWeldBeadFeatureData = default(CosmeticWeldBeadFeatureData);
            SelectData swSelData;
            Feature swFeature;
            MathTransform swXform;
            Edge[] edges;
            FeatureManager swFeatureManager;
            CosmeticWeldBeadFolder swWeldFolder = default(CosmeticWeldBeadFolder);
            IFeatureFolder swFeatFolder = default(IFeatureFolder);            
            XmlDocument XMLFile = new XmlDocument();
            XmlDocument XMLLog = new XmlDocument();

            //XMLFile.Load("C:\\Users\\Matei\\OneDrive\\Desktop\\ElementaryOperations\\weldments.xml");
            XmlElement RootNode = XMLFile.CreateElement("Weldments");
            XmlElement LogRootNode = XMLLog.CreateElement("PoseData");
            XMLFile.AppendChild(RootNode);
            XMLLog.AppendChild(LogRootNode);
            Double stickout_D = 0.018;
            float stickout = Convert.ToSingle(stickout_D);

            Double[] weldment_start_coordinates;
            Double[] weldment_end_coordinates;

            Object[] weldment_edges_obj;

            
            bool feature_successfuly_selected = false;



            if (swApp == null)
            {
                swApp = Activator.CreateInstance(Type.GetTypeFromProgID("SldWorks.Application")) as SldWorks;
            }

            swModel = (ModelDoc2)swApp.ActiveDoc;
            String file_address = swModel.GetPathName();


            swModelDocExt = swModel.Extension;


            if (((swModel != null)))
            {
                
                RootNode.SetAttribute("MeasurementUnit", "mm");
                RootNode.SetAttribute("SWFile", file_address);
                file_address = file_address.Replace(".SLDASM", ".STEP");
                RootNode.SetAttribute("STEPFile", file_address);



                int weld_index = 1;
                while (true)
                {
                    bool added_approach = false;


                    feature_successfuly_selected = swModelDocExt.SelectByID2("Weld Bead" + weld_index.ToString(), "COSMETIC_WELD", 0, 0, 0, false, 0, null, 0);
                    if (!feature_successfuly_selected)
                        break;

                    XmlElement WeldmentElement = XMLFile.CreateElement("Weldment");
                    RootNode.AppendChild(WeldmentElement);
                    WeldmentElement.SetAttribute("No", Convert.ToString(weld_index));

                    swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
                    swFeature = (Feature)swSelectionMgr.GetSelectedObject6(1, 0);
                    swXform = swModelDocExt.GetCoordinateSystemTransformByName(swFeature.Name);

                    swCosmeticWeldBeadFeatureData = (CosmeticWeldBeadFeatureData)swFeature.GetDefinition();
                    swCosmeticWeldBeadFeatureData.AccessSelections(swModel, null);
                    int d;
                    Object[] FromFaceObj = (Object[])swCosmeticWeldBeadFeatureData.GetEntitiesWeldFrom(out d);
                    Object[] ToFacesObj = (Object[])swCosmeticWeldBeadFeatureData.GetEntitiesWeldTo(out d);
                    Face2[] FromFace = null;
                    Face2[] ToFaces = null;
                    if (FromFaceObj != null && ToFacesObj != null)
                    {
                        FromFace = new Face2[FromFaceObj.Length];
                        ToFaces = new Face2[ToFacesObj.Length];
                        for (int i = 0; i < FromFaceObj.Length; i++)
                        {
                            FromFace[i] = FromFaceObj[i] as Face2;
                        }
                        for (int i = 0; i < ToFacesObj.Length; i++)
                        {
                            ToFaces[i] = ToFacesObj[i] as Face2;
                        }
                    }


                    //Face2[] FromFace = FromFaceObj as Face2[];
                    //Face2[] ToFaces = ToFacesObj as Face2[];



                    weldment_edges_obj = (Object[])swCosmeticWeldBeadFeatureData.GetReferenceEdges();



                    Double[] pointOnEdgeCoordinates = { 0, 0, 0 };
                    int ewoindex = 1;
                    for (int edge_index = 0; edge_index < weldment_edges_obj.Length; edge_index++)
                    {
                        XmlElement EdgeElement = XMLLog.CreateElement("Edge");
                        LogRootNode.AppendChild(EdgeElement);
                        EdgeElement.SetAttribute("No", Convert.ToString(edge_index));

                        Edge weldment_edge = (Edge)weldment_edges_obj[edge_index];


                        bool circle_edge = false;

                        if (weldment_edge.GetStartVertex() != null)
                        {

                            Vertex weldment_edge_start_vertex = (Vertex)weldment_edge.GetStartVertex();
                            Vertex weldment_edge_end_vertex = (Vertex)weldment_edge.GetEndVertex();


                            weldment_start_coordinates = (Double[])weldment_edge_start_vertex.GetPoint();
                            weldment_end_coordinates = (Double[])weldment_edge_end_vertex.GetPoint();

                        }
                        else
                        {
                            object pointOnEdgeObj = weldment_edge.GetClosestPointOn(0, 0, 0);
                            pointOnEdgeCoordinates = pointOnEdgeObj as Double[];
                            weldment_start_coordinates = pointOnEdgeCoordinates;
                            weldment_end_coordinates = pointOnEdgeCoordinates;
                            circle_edge = true;

                        }

                        Curve weldment_curve = (Curve)weldment_edge.GetCurve();
                        object weldment_curve_tessalation_points_Obj = weldment_curve.GetTessPts(0.0001, 0.0001, weldment_start_coordinates, weldment_end_coordinates);
                        Double[] weldment_curve_tessalation_points = weldment_curve_tessalation_points_Obj as Double[];

                        int point_index = 0;

                        Matrix<double> TessalationPoints = Matrix<double>.Build.Random(3, weldment_curve_tessalation_points.Length / 3);
                        int column_index = 0;

                        for (point_index = 0; point_index < weldment_curve_tessalation_points.Length; point_index += 3)
                        {
                            TessalationPoints[0, column_index] = weldment_curve_tessalation_points[point_index];
                            TessalationPoints[1, column_index] = weldment_curve_tessalation_points[point_index + 1];
                            TessalationPoints[2, column_index] = weldment_curve_tessalation_points[point_index + 2];

                            column_index += 1;
                        }

                        MathNet.Numerics.LinearAlgebra.Vector<Double> SelectionPoint = MathNet.Numerics.LinearAlgebra.Vector<Double>.Build.Random(3);
                        if (circle_edge)
                        {
                            SelectionPoint[0] = weldment_start_coordinates[0];
                            SelectionPoint[1] = weldment_start_coordinates[1];
                            SelectionPoint[2] = weldment_start_coordinates[2];
                        }
                        else if (TessalationPoints.ColumnCount > 2)
                        {
                            SelectionPoint[0] = TessalationPoints[0, TessalationPoints.ColumnCount / 2];
                            SelectionPoint[1] = TessalationPoints[1, TessalationPoints.ColumnCount / 2];
                            SelectionPoint[2] = TessalationPoints[2, TessalationPoints.ColumnCount / 2];
                        }
                        else
                        {
                            SelectionPoint[0] = (TessalationPoints[0, 0] + TessalationPoints[0, 1]) / 2;
                            SelectionPoint[1] = (TessalationPoints[1, 0] + TessalationPoints[1, 1]) / 2;
                            SelectionPoint[2] = (TessalationPoints[2, 0] + TessalationPoints[2, 1]) / 2;

                        }

                        feature_successfuly_selected = swModelDocExt.SelectByID2("", "EDGE", SelectionPoint[0], SelectionPoint[1], SelectionPoint[2], false, 0, null, 0);
                        swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
                        weldment_edge = (Edge)swSelectionMgr.GetSelectedObject6(1, 0);
                        weldment_edge.IGetTwoAdjacentFaces2(out Face2 face1, out Face2 face2);

                        Debug.Print("Face1 area " + face1.GetArea().ToString());
                        Debug.Print("Face2 area " + face2.GetArea().ToString());

                        Face2 JoiningFace = null, BaseFace = null;
                        bool curvedJoiningFace = false;
                        MathTransform FaceTransform = null, BaseTransform = null;
                        MathNet.Numerics.LinearAlgebra.Vector<Double> BaseNormal;
                        MathNet.Numerics.LinearAlgebra.Vector<Double> JoiningNormal;
                        bool found_base = false;
                        if (FromFace != null && ToFaces != null)
                        {
                            
                            for (int index = 0; index < ToFaces.Length; index++)
                            {
                                if (ToFaces[index].GetArea() == face1.GetArea() || ToFaces[index].GetArea() == face2.GetArea())
                                {
                                    JoiningFace = ToFaces[index];
                                    BaseFace = FromFace[0];
                                    found_base = true;
                                    break;
                                }
                            }
                            if (!found_base)
                            {
                                JoiningFace = FromFace[0];
                                Feature FaceFeature1 = (Feature)JoiningFace.GetFeature();
                                String FeatureName;
                                FeatureName = FaceFeature1.GetNameForSelection(out FeatureName);
                                int pFrom = FeatureName.IndexOf("@") + "@".Length;
                                int pTo = FeatureName.LastIndexOf("@");
                                String result = FeatureName.Substring(pFrom, pTo - pFrom);
                                feature_successfuly_selected = swModelDocExt.SelectByID2(result, "COMPONENT", 0, 0, 0, false, 0, null, 0);
                                swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
                                Component2 comp = (Component2)swSelectionMgr.GetSelectedObject6(1, 0);
                                FaceTransform = comp.Transform;

                                comp.DeSelect();
                                Double[] FaceNormalDouble;
                                object strips;
                                if (JoiningFace.GetArea() == face1.GetArea())
                                {
                                    FaceNormalDouble = (Double[])face2.Normal;
                                    strips = (object)face2.GetTessNorms();
                                }
                                else
                                {
                                    FaceNormalDouble = (Double[])face1.Normal;
                                    strips = (object)face1.GetTessNorms();
                                }
                                MathNet.Numerics.LinearAlgebra.Vector<Double> FaceNormal = computeTransform2(FaceTransform, FaceNormalDouble[0], FaceNormalDouble[1], FaceNormalDouble[2]);
                                if (FaceNormal.Average() != 0)
                                {
                                    MathNet.Numerics.LinearAlgebra.Vector<Double> RayPoint = DenseVector.OfArray(new double[] { (weldment_end_coordinates[0] + weldment_end_coordinates[0]) / 2, (weldment_end_coordinates[1] + weldment_end_coordinates[1]) / 2, (weldment_end_coordinates[2] + weldment_end_coordinates[2]) });
                                    swModelDocExt.SelectByRay(SelectionPoint[0], SelectionPoint[1], SelectionPoint[2], FaceNormal[0], FaceNormal[1], FaceNormal[2], 0.1, 2, true, 0, 0);
                                }
                                else
                                {
                                    
                                    float[] stripsa = strips as float[];
                                    Double first_normal_x, first_normal_y, first_normal_z, last_normal_x, last_normal_y, last_normal_z;

                                    if (stripsa.Length % 3 == 0)
                                    {
                                        first_normal_x = stripsa[3];
                                        first_normal_y = stripsa[4];
                                        first_normal_z = stripsa[5];
                                        last_normal_x = stripsa[stripsa.Length - 3];
                                        last_normal_y = stripsa[stripsa.Length - 2];
                                        last_normal_z = stripsa[stripsa.Length - 1];

                                    }
                                    else
                                    {
                                        first_normal_x = stripsa[2];
                                        first_normal_y = stripsa[3];
                                        first_normal_z = stripsa[4];
                                        last_normal_x = stripsa[stripsa.Length - 3];
                                        last_normal_y = stripsa[stripsa.Length - 2];
                                        last_normal_z = stripsa[stripsa.Length - 1];
                                    }
                                    MathNet.Numerics.LinearAlgebra.Vector<Double> firstNormal, lastNormal;
                                    
                                    firstNormal = computeTransform2(FaceTransform, first_normal_x, first_normal_y, first_normal_z);
                                    lastNormal = computeTransform2(FaceTransform, last_normal_x, last_normal_y, last_normal_z);
                                    MathNet.Numerics.LinearAlgebra.Vector<Double> avgNormal = (firstNormal + lastNormal) / 2;
                                    
                                    swModelDocExt.SelectByRay(SelectionPoint[0], SelectionPoint[1], SelectionPoint[2], avgNormal[0], avgNormal[1], avgNormal[2], 0.1, 2, true, 0, 0);
                                    
                                }
                                BaseFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                swModel.ClearSelection2(true);
                                swModel.ClearSelection2(false);
                                //FaceFeature1 = (Feature)BaseFace.GetFeature();
                                
                                //FeatureName = FaceFeature1.GetNameForSelection(out FeatureName);
                                //pFrom = FeatureName.IndexOf("@") + "@".Length;
                                //pTo = FeatureName.LastIndexOf("@");
                                //result = FeatureName.Substring(pFrom, pTo - pFrom);
                                //feature_successfuly_selected = swModelDocExt.SelectByID2(result, "COMPONENT", 0, 0, 0, false, 0, null, 0);
                                //swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
                                //comp = (Component2)swSelectionMgr.GetSelectedObject6(1, 0);
                                //BaseTransform = comp.Transform;

                                //comp.DeSelect();


                            }
                        }

                        if (JoiningFace == null && BaseFace == null)
                        {

                            Feature FaceFeature1 = (Feature)face1.GetFeature();
                            String FeatureName;
                            FeatureName = FaceFeature1.GetNameForSelection(out FeatureName);
                            int pFrom = FeatureName.IndexOf("@") + "@".Length;
                            int pTo = FeatureName.LastIndexOf("@");
                            String result = FeatureName.Substring(pFrom, pTo - pFrom);
                            feature_successfuly_selected = swModelDocExt.SelectByID2(result, "COMPONENT", 0, 0, 0, false, 0, null, 0);
                            swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
                            Component2 comp = (Component2)swSelectionMgr.GetSelectedObject6(1, 0);
                            FaceTransform = comp.Transform;
                            comp.DeSelect();
                            Double[] FaceNormalDouble1 = (Double[])face1.Normal;
                            Double[] FaceNormalDouble2 = (Double[])face2.Normal;
                            MathNet.Numerics.LinearAlgebra.Vector<Double> FaceNormal1 = computeTransform2(FaceTransform, FaceNormalDouble1[0], FaceNormalDouble1[1], FaceNormalDouble1[2]);
                            MathNet.Numerics.LinearAlgebra.Vector<Double> FaceNormal2 = computeTransform2(FaceTransform, FaceNormalDouble2[0], FaceNormalDouble2[1], FaceNormalDouble2[2]);


                            MathNet.Numerics.LinearAlgebra.Vector<Double> RayPoint = DenseVector.OfArray(new double[] { (weldment_end_coordinates[0] + weldment_start_coordinates[0]) / 2, (weldment_end_coordinates[1] + weldment_start_coordinates[1]) / 2, (weldment_end_coordinates[2] + weldment_start_coordinates[2]) });
                            Face2 selectedFace;

                            BaseNormal = MathNet.Numerics.LinearAlgebra.Vector<Double>.Build.Random(3);
                            JoiningNormal = MathNet.Numerics.LinearAlgebra.Vector<Double>.Build.Random(3);
                            Double[] JoiningNormalDouble = { 0, 0, 0 };

                            if (FaceNormal1.Average() == 0 || (FaceNormal2.Average() == 0))
                            {
                                curvedJoiningFace = true;
                                if (FaceNormal1.Average() == 0)
                                {
                                    JoiningFace = face1;
                                    swModelDocExt.SelectByRay(SelectionPoint[0], SelectionPoint[1], SelectionPoint[2], FaceNormal2[0], FaceNormal2[1], FaceNormal2[2], 0.1, 2, true, 0, 0);
                                    selectedFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                    Double[] selectedFaceNormalDouble = (Double[])selectedFace.Normal;
                                    MathNet.Numerics.LinearAlgebra.Vector<Double> selectedFaceNormal = computeTransform2(FaceTransform, selectedFaceNormalDouble[0], selectedFaceNormalDouble[1], selectedFaceNormalDouble[2]);
                                    if (selectedFaceNormal[0] == FaceNormal2[0] & selectedFaceNormal[1] == FaceNormal2[1] & selectedFaceNormal[2] == FaceNormal2[2])
                                    {
                                        swModel.ClearSelection2(true);
                                        swModel.ClearSelection2(false);
                                        swModelDocExt.SelectByRay(RayPoint[0], RayPoint[1], RayPoint[2], -FaceNormal2[0], -FaceNormal2[1], -FaceNormal2[2], 0.1, 2, true, 0, 0);
                                        selectedFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                    }

                                    BaseFace = selectedFace;
                                }
                                else
                                {
                                    JoiningFace = face2;
                                    swModelDocExt.SelectByRay(RayPoint[0], RayPoint[1], RayPoint[2], FaceNormal1[0], FaceNormal1[1], FaceNormal1[2], 0.1, 2, true, 0, 0);

                                    selectedFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                    Double[] selectedFaceNormalDouble = (Double[])selectedFace.Normal;
                                    MathNet.Numerics.LinearAlgebra.Vector<Double> selectedFaceNormal = computeTransform2(FaceTransform, selectedFaceNormalDouble[0], selectedFaceNormalDouble[1], selectedFaceNormalDouble[2]);
                                    if (selectedFaceNormal[0] == FaceNormal1[0] & selectedFaceNormal[1] == FaceNormal1[1] & selectedFaceNormal[2] == FaceNormal1[2])
                                    {
                                        swModel.ClearSelection2(true);
                                        swModel.ClearSelection2(false);
                                        swModelDocExt.SelectByRay(RayPoint[0], RayPoint[1], RayPoint[2], -FaceNormal1[0], -FaceNormal1[1], -FaceNormal1[2], 0.1, 2, true, 0, 0);
                                        selectedFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                    }
                                    BaseFace = selectedFace;
                                }
                            }
                            else
                            {
                                MathNet.Numerics.LinearAlgebra.Vector<Double> bisector = (FaceNormal1 + FaceNormal2) / 2;
                                MathNet.Numerics.LinearAlgebra.Vector<Double> point1 = 0.01 * bisector;
                                point1 = point1 + RayPoint;
                                MathNet.Numerics.LinearAlgebra.Vector<Double> point2 = point1 - 0.04 * FaceNormal1;
                                MathNet.Numerics.LinearAlgebra.Vector<Double> point3 = point1 - 0.04 * FaceNormal2;

                                RayPoint = RayPoint - point1;
                                Face2 jFace = null, bFace = null;
                                MathNet.Numerics.LinearAlgebra.Vector<Double> jNormal = MathNet.Numerics.LinearAlgebra.Vector<Double>.Build.Random(3);
                                MathNet.Numerics.LinearAlgebra.Vector<Double> bNormal = MathNet.Numerics.LinearAlgebra.Vector<Double>.Build.Random(3);
                                swModelDocExt.SelectByRay(point2[0], point2[1], point2[2], -FaceNormal1[0], -FaceNormal1[1], -FaceNormal1[2], 0.001, 2, true, 0, 0);
                                selectedFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                if (selectedFace != null)
                                {
                                    FaceFeature1 = (Feature)selectedFace.GetFeature();
                                    FeatureName = FaceFeature1.GetNameForSelection(out FeatureName);
                                    if (FeatureName.Contains(result))
                                    {
                                        jFace = selectedFace;
                                        jNormal = FaceNormal1;
                                        bNormal = -FaceNormal2;
                                        swModel.ClearSelection2(true);
                                        swModel.ClearSelection2(false);
                                        swModelDocExt.SelectByRay(point2[0], point2[1], point2[2], -bNormal[0], -bNormal[1], -bNormal[2], 0.001, 2, true, 0, 0);
                                        bFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                    }
                                }
                                swModel.ClearSelection2(true);
                                swModel.ClearSelection2(false);
                                swModelDocExt.SelectByRay(point2[0], point2[1], point2[2], -FaceNormal2[0], -FaceNormal2[1], -FaceNormal2[2], 0.001, 2, true, 0, 0);
                                selectedFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                if (selectedFace != null)
                                {
                                    FaceFeature1 = (Feature)selectedFace.GetFeature();
                                    FeatureName = FaceFeature1.GetNameForSelection(out FeatureName);
                                    if (FeatureName.Contains(result))
                                    {
                                        jFace = selectedFace;
                                        jNormal = FaceNormal2;
                                        bNormal = -FaceNormal1;
                                        swModel.ClearSelection2(true);
                                        swModel.ClearSelection2(false);
                                        swModelDocExt.SelectByRay(point2[0], point2[1], point2[2], -bNormal[0], -bNormal[1], -bNormal[2], 0.001, 2, true, 0, 0);
                                        bFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                    }
                                }
                                swModel.ClearSelection2(true);
                                swModel.ClearSelection2(false);
                                swModelDocExt.SelectByRay(point3[0], point3[1], point3[2], -FaceNormal1[0], -FaceNormal1[1], -FaceNormal1[2], 0.001, 2, true, 0, 0);
                                selectedFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                if (selectedFace != null)
                                {
                                    FaceFeature1 = (Feature)selectedFace.GetFeature();
                                    FeatureName = FaceFeature1.GetNameForSelection(out FeatureName);
                                    if (FeatureName.Contains(result))
                                    {
                                        jFace = selectedFace;
                                        jNormal = FaceNormal1;
                                        bNormal = -FaceNormal2;
                                        swModel.ClearSelection2(true);
                                        swModel.ClearSelection2(false);
                                        swModelDocExt.SelectByRay(point3[0], point3[1], point3[2], -bNormal[0], -bNormal[1], -bNormal[2], 0.001, 2, true, 0, 0);
                                        bFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                    }
                                }
                                swModel.ClearSelection2(true);
                                swModel.ClearSelection2(false);
                                swModelDocExt.SelectByRay(point3[0], point3[1], point3[2], -FaceNormal2[0], -FaceNormal2[1], -FaceNormal2[2], 0.001, 2, true, 0, 0);
                                selectedFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                if (selectedFace != null)
                                {
                                    FaceFeature1 = (Feature)selectedFace.GetFeature();
                                    FeatureName = FaceFeature1.GetNameForSelection(out FeatureName);
                                    if (FeatureName.Contains(result))
                                    {
                                        jFace = selectedFace;
                                        jNormal = FaceNormal2;
                                        bNormal = -FaceNormal1;
                                        swModel.ClearSelection2(true);
                                        swModel.ClearSelection2(false);
                                        swModelDocExt.SelectByRay(point3[0], point3[1], point3[2], -bNormal[0], -bNormal[1], -bNormal[2], 0.001, 2, true, 0, 0);
                                        bFace = (Face2)swSelectionMgr.GetSelectedObject6(1, 0);
                                    }
                                }
                                JoiningFace = jFace;
                                BaseNormal = bNormal;
                                JoiningNormal = jNormal;




                            }

                        }
                        else
                        {
                            Feature JoiningFaceFeature = (Feature)JoiningFace.GetFeature();
                            String FeatureName;
                            FeatureName = JoiningFaceFeature.GetNameForSelection(out FeatureName);
                            int pFrom = FeatureName.IndexOf("@") + "@".Length;
                            int pTo = FeatureName.LastIndexOf("@");
                            String result = FeatureName.Substring(pFrom, pTo - pFrom);
                            feature_successfuly_selected = swModelDocExt.SelectByID2(result, "COMPONENT", 0, 0, 0, false, 0, null, 0);
                            swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
                            Component2 comp = (Component2)swSelectionMgr.GetSelectedObject6(1, 0);
                            FaceTransform = comp.Transform;
                            comp.DeSelect();
                            Double[] JoiningNormalDouble = (Double[])JoiningFace.Normal;
                            JoiningNormal = computeTransform2(FaceTransform, JoiningNormalDouble[0], JoiningNormalDouble[1], JoiningNormalDouble[2]);

                            Feature BaseFaceFeature = (Feature)BaseFace.GetFeature();
                            
                            FeatureName = BaseFaceFeature.GetNameForSelection(out FeatureName);
                            pFrom = FeatureName.IndexOf("@") + "@".Length;
                            pTo = FeatureName.LastIndexOf("@");
                            result = FeatureName.Substring(pFrom, pTo - pFrom);
                            feature_successfuly_selected = swModelDocExt.SelectByID2(result, "COMPONENT", 0, 0, 0, false, 0, null, 0);
                            swSelectionMgr = (SelectionMgr)swModel.SelectionManager;
                            comp = (Component2)swSelectionMgr.GetSelectedObject6(1, 0);
                            BaseTransform = comp.Transform;
                            comp.DeSelect();
                            Double[] BaseNormalDouble = (Double[])JoiningFace.Normal;
                            BaseNormal = computeTransform2(BaseTransform, BaseNormalDouble[0], BaseNormalDouble[1], BaseNormalDouble[2]);

                        }



                        MathTransform JoiningTransform = FaceTransform;


                        if (curvedJoiningFace)
                        {

                            Double first_normal_x, first_normal_y, first_normal_z, last_normal_x, last_normal_y, last_normal_z, mid_normal_x, mid_normal_y, mid_normal_z;
                            object strips = (object)JoiningFace.GetTessTriStripNorms();
                            float[] stripsa = strips as float[];

                            object norms = (object)JoiningFace.GetTessNorms();
                            float[] normss = norms as float[];

                            MathNet.Numerics.LinearAlgebra.Vector < Double > FirstNormal = computeTransform2(JoiningTransform, Convert.ToSingle(normss[0]), Convert.ToSingle(normss[1]), Convert.ToSingle(normss[2]));
                            MathNet.Numerics.LinearAlgebra.Vector < Double > LastNormal = computeTransform2(JoiningTransform, Convert.ToSingle(normss[normss.Length - 3]), Convert.ToSingle(normss[normss.Length - 2]), Convert.ToSingle(normss[normss.Length - 1]));
                            MathNet.Numerics.LinearAlgebra.Vector<Double> MidNormal = computeTransform2(JoiningTransform, Convert.ToSingle(normss[(normss.Length/6)*3]), Convert.ToSingle(normss[(normss.Length/6)*3+1]), Convert.ToSingle(normss[(normss.Length/6)*3+2]));

                            MathNet.Numerics.LinearAlgebra.Vector<Double> FirstWeldmentNormal = (FirstNormal + BaseNormal) / 2;
                            MathNet.Numerics.LinearAlgebra.Vector<Double> LastWeldmentNormal = (LastNormal + BaseNormal) / 2;

                            MathNet.Numerics.LinearAlgebra.Vector<Double> DifNormal = (FirstNormal - MidNormal) / ((TessalationPoints.ColumnCount-1) / 2);

                            for (int i = 0; i<3; i++)
                            {
                                if (FirstNormal[i] > MidNormal[i])                                
                                    DifNormal[i] = -Math.Abs(DifNormal[i]);
                                else
                                    DifNormal[i] = Math.Abs(DifNormal[i]);
                                
                            }

                            Matrix<Double> WeldmentNormals = Matrix<Double>.Build.Dense(3, TessalationPoints.ColumnCount - 1);
                            
                            WeldmentNormals[0, 0] = FirstNormal[0];
                            WeldmentNormals[1, 0] = FirstNormal[1];
                            WeldmentNormals[2, 0] = FirstNormal[2];

                            for (int index_tess = 1; index_tess < TessalationPoints.ColumnCount / 2; index_tess++)
                            {
                                WeldmentNormals[0, index_tess] = (WeldmentNormals[0, index_tess - 1] + DifNormal[0]);
                                WeldmentNormals[1, index_tess] = (WeldmentNormals[1, index_tess - 1] + DifNormal[1]);
                                WeldmentNormals[2, index_tess] = (WeldmentNormals[2, index_tess - 1] + DifNormal[2]);

                            }

                            DifNormal = (MidNormal - LastNormal) / ((TessalationPoints.ColumnCount) / 2);

                            for (int i = 0; i < 3; i++)
                            {
                                if (MidNormal[i] > LastNormal[i])
                                    DifNormal[i] = -Math.Abs(DifNormal[i]);
                                else
                                    DifNormal[i] = Math.Abs(DifNormal[i]);

                            }

                            WeldmentNormals[0, TessalationPoints.ColumnCount / 2] = MidNormal[0];
                            WeldmentNormals[1, TessalationPoints.ColumnCount / 2] = MidNormal[1];
                            WeldmentNormals[2, TessalationPoints.ColumnCount / 2] = MidNormal[2];

                            for (int index_tess = TessalationPoints.ColumnCount / 2 + 1; index_tess < TessalationPoints.ColumnCount - 1; index_tess++)
                            {
                                WeldmentNormals[0, index_tess] = (WeldmentNormals[0, index_tess - 1] + DifNormal[0]);
                                WeldmentNormals[1, index_tess] = (WeldmentNormals[1, index_tess - 1] + DifNormal[1]);
                                WeldmentNormals[2, index_tess] = (WeldmentNormals[2, index_tess - 1] + DifNormal[2]);

                            }

                            WeldmentNormals = WeldmentNormals / 100;

                            //for (int index = 0; index < WeldmentNormals.ColumnCount; index++)
                            //{
                            //    WeldmentNormals[0, index] = WeldmentNormals[0, index] * stickout;
                            //    WeldmentNormals[1, index] = WeldmentNormals[1, index] * stickout;
                            //    WeldmentNormals[2, index] = WeldmentNormals[2, index] * stickout;

                            //}

                            int norm_index;
                            MathNet.Numerics.LinearAlgebra.Vector<Double> firstNormal, lastNormal, midNormal;
                            if (stripsa.Length % 3 == 0)
                            {
                                first_normal_x = stripsa[3];
                                first_normal_y = stripsa[4];
                                first_normal_z = stripsa[5];
                                norm_index = stripsa.Length / 3;
                                if (norm_index % 2 != 0)
                                    norm_index++;
                                norm_index = (norm_index / 2) * 3;
                                mid_normal_x = stripsa[norm_index];
                                mid_normal_y = stripsa[norm_index + 1];
                                mid_normal_z = stripsa[norm_index + 2];

                                last_normal_x = stripsa[stripsa.Length - 3];
                                last_normal_y = stripsa[stripsa.Length - 2];
                                last_normal_z = stripsa[stripsa.Length - 1];

                            }
                            else
                            {
                                first_normal_x = stripsa[2];
                                first_normal_y = stripsa[3];
                                first_normal_z = stripsa[4];
                                norm_index = stripsa.Length / 3;
                                if (norm_index % 2 != 0)
                                    norm_index++;
                                norm_index = (norm_index / 2) * 3;
                                mid_normal_x = stripsa[norm_index-1];
                                mid_normal_y = stripsa[norm_index];
                                mid_normal_z = stripsa[norm_index + 1];
                                last_normal_x = stripsa[stripsa.Length - 3];
                                last_normal_y = stripsa[stripsa.Length - 2];
                                last_normal_z = stripsa[stripsa.Length - 1];
                            }
                            firstNormal = computeTransform2(JoiningTransform, first_normal_x, first_normal_y, first_normal_z);
                            lastNormal = computeTransform2(JoiningTransform, last_normal_x, last_normal_y, last_normal_z);
                            midNormal = computeTransform2(JoiningTransform, mid_normal_x, mid_normal_y, mid_normal_z);


                            Double weldment_angle = computeAngle(BaseNormal[0], BaseNormal[1], BaseNormal[2], first_normal_x, first_normal_y, first_normal_z);

                            Double angle_degre = computeAngle(firstNormal[0], firstNormal[1], firstNormal[2], lastNormal[0], lastNormal[1], lastNormal[2]);
                            Double[] tessalation_normals_x;
                            Double[] tessalation_normals_y;
                            Double[] tessalation_normals_z;
                            Double angle_check;
                            if (angle_degre != 0)
                            {
                                int points_no = weldment_curve_tessalation_points.Length / 3;

                                Double x_difference = (lastNormal[0] - firstNormal[0]) / points_no;
                                Double y_difference = (lastNormal[1] - firstNormal[1]) / points_no;
                                Double z_difference = (lastNormal[2] - firstNormal[2]) / points_no;

                                tessalation_normals_x = new Double[points_no - 1];
                                tessalation_normals_y = new Double[points_no - 1];
                                tessalation_normals_z = new Double[points_no - 1];

                                tessalation_normals_x[0] = firstNormal[0] + x_difference / 2;
                                tessalation_normals_y[0] = firstNormal[1] + y_difference / 2;
                                tessalation_normals_z[0] = firstNormal[2] + z_difference / 2;

                                for (int index = 1; index < points_no - 1; index++)
                                {
                                    tessalation_normals_x[index] = tessalation_normals_x[index - 1] + x_difference;
                                    tessalation_normals_y[index] = tessalation_normals_y[index - 1] + y_difference;
                                    tessalation_normals_z[index] = tessalation_normals_z[index - 1] + z_difference;


                                }
                                Vector3[] vectorNormals = new Vector3[tessalation_normals_x.Length];
                                
                                for (int index = 0; index < tessalation_normals_x.Length; index++)
                                {
                                    vectorNormals[index] = new Vector3(Convert.ToSingle(tessalation_normals_x[index]),
                                        Convert.ToSingle(tessalation_normals_y[index]),
                                        Convert.ToSingle(tessalation_normals_z[index]));

                                    vectorNormals[index] = vectorNormals[index] * stickout;
                                }

                                int index_normals = 0;
                                for (int index = 0; index < TessalationPoints.ColumnCount - 1; index += 1)

                                {

                                    XmlElement EWOLogElement = XMLLog.CreateElement("EWO");
                                    LogRootNode.AppendChild(EWOLogElement);
                                    EWOLogElement.SetAttribute("No", Convert.ToString(index));
                                    Double midX = (TessalationPoints[0, index] + TessalationPoints[0, index + 1]) / 2;
                                    Double midY = (TessalationPoints[1, index] + TessalationPoints[1, index + 1]) / 2;
                                    Double midZ = (TessalationPoints[2, index] + TessalationPoints[2, index + 1]) / 2;


                                    Matrix<Double> A = DenseMatrix.OfArray(new double[,] {
                                        { TessalationPoints[0, index], midX, TessalationPoints[0, index + 1], TessalationPoints[0, index] + WeldmentNormals[0, index],  midX + WeldmentNormals[0, index], TessalationPoints[0, index + 1] + WeldmentNormals[0, index]},
                                        { TessalationPoints[1, index], midY, TessalationPoints[1, index + 1], TessalationPoints[1, index] + WeldmentNormals[1, index],  midY + WeldmentNormals[1, index], TessalationPoints[1, index + 1] + WeldmentNormals[1, index]},
                                        { TessalationPoints[2, index], midZ, TessalationPoints[2, index + 1], TessalationPoints[2, index] + WeldmentNormals[2, index],  midZ + WeldmentNormals[2, index], TessalationPoints[2, index + 1] + WeldmentNormals[2, index]}
                                        });

                                    String A_string = "[";

                                    for(int row_index = 0; row_index < A.RowCount; row_index++)
                                    {
                                        
                                        for (int col_index = 0; col_index < A.ColumnCount; col_index++)
                                        {
                                            A_string = A_string + Convert.ToString(A[row_index, col_index]);
                                            if (col_index < A.ColumnCount - 1)
                                                A_string = A_string + ", ";
                                            

                                        }
                                        if (row_index < A.RowCount - 1)
                                            A_string = A_string + "; ";
                                        else
                                            A_string = A_string + "]";
                                    }

                                    EWOLogElement.SetAttribute("A", A_string);

                                    Matrix<Double> A_aux = DenseMatrix.OfArray(new double[,] {
                                        { TessalationPoints[0, index], midX, TessalationPoints[0, index + 1], TessalationPoints[0, index] + WeldmentNormals[0, index],  midX + WeldmentNormals[0, index], TessalationPoints[0, index + 1] + WeldmentNormals[0, index]},
                                        { TessalationPoints[1, index], midY, TessalationPoints[1, index + 1], TessalationPoints[1, index] + WeldmentNormals[1, index],  midY + WeldmentNormals[1, index], TessalationPoints[1, index + 1] + WeldmentNormals[1, index]},
                                        { TessalationPoints[2, index], midZ, TessalationPoints[2, index + 1], TessalationPoints[2, index] + WeldmentNormals[2, index],  midZ + WeldmentNormals[2, index], TessalationPoints[2, index + 1] + WeldmentNormals[2, index]}
                                        });

                                    Matrix<Double> B = DenseMatrix.OfArray(new double[,]
                                        {
                                            { 0, 0, 0, 1, 1, 1},
                                            { 0, Distance(A.Column(0), A.Column(1)), Distance(A.Column(0), A.Column(2)), 0, Distance(A.Column(0), A.Column(1)), Distance(A.Column(0), A.Column(2)) },
                                            { 0, 0, 0, 0, 0, 0}
                                        });

                                    String B_string = "[";

                                    for (int row_index = 0; row_index < B.RowCount; row_index++)
                                    {

                                        for (int col_index = 0; col_index < B.ColumnCount; col_index++)
                                        {
                                            B_string = B_string + Convert.ToString(B[row_index, col_index]);
                                            if (col_index < B.ColumnCount - 1)
                                                B_string = B_string + ", ";
                                            

                                        }
                                        if (row_index < B.RowCount - 1)
                                            B_string = B_string + "; ";
                                        else
                                            B_string = B_string + "]";
                                    }

                                    EWOLogElement.SetAttribute("B", B_string);

                                    Matrix<Double> B_aux = DenseMatrix.OfArray(new double[,]
                                        {
                                            { 0, 0, 0, 1, 1, 1},
                                            { 0, Distance(A.Column(0), A.Column(1)), Distance(A.Column(0), A.Column(2)), 0, Distance(A.Column(0), A.Column(1)), Distance(A.Column(0), A.Column(2)) },
                                            { 0, 0, 0, 0, 0, 0}
                                        });

                                    Matrix<Double> R = AbsoluteOrientation(A_aux, B_aux);


                                    String R_string = "[";

                                    for (int row_index = 0; row_index < R.RowCount; row_index++)
                                    {

                                        for (int col_index = 0; col_index < R.ColumnCount; col_index++)
                                        {
                                            R_string = R_string + Convert.ToString(R[row_index, col_index]);
                                            if (col_index < R.ColumnCount - 1)
                                                R_string = R_string + ", ";
                                            

                                        }
                                        if (row_index < R.RowCount - 1)
                                            R_string = R_string + "; ";
                                        else
                                            R_string = R_string + "]";
                                    }

                                    EWOLogElement.SetAttribute("R", R_string);



                                    String StartPose = "[[" + Convert.ToString(R[0, 0]) + ", " + Convert.ToString(R[0, 1]) + ", " + Convert.ToString(R[0, 2]) + ", " + Convert.ToString(TessalationPoints[0, index]) + "], [" +
                                        Convert.ToString(R[1, 0]) + ", " + Convert.ToString(R[1, 1]) + ", " + Convert.ToString(R[1, 2]) + ", " + Convert.ToString(TessalationPoints[1, index]) + "], [" +
                                        Convert.ToString(R[2, 0]) + ", " + Convert.ToString(R[2, 1]) + ", " + Convert.ToString(R[2, 2]) + ", " + Convert.ToString(TessalationPoints[2, index]) + "], " +
                                        "[0, 0, 0, 1]]";

                                    String EndPose = "[[" + Convert.ToString(R[0, 0]) + ", " + Convert.ToString(R[0, 1]) + ", " + Convert.ToString(R[0, 2]) + ", " + Convert.ToString(A[0, 4]) + "], [" +
                                        Convert.ToString(R[1, 0]) + ", " + Convert.ToString(R[1, 1]) + ", " + Convert.ToString(R[1, 2]) + ", " + Convert.ToString(A[1, 4]) + "], [" +
                                        Convert.ToString(R[2, 0]) + ", " + Convert.ToString(R[2, 1]) + ", " + Convert.ToString(R[2, 2]) + ", " + Convert.ToString(A[2, 4]) + "], " +
                                        "[0, 0, 0, 1]]";

                                    //String SeamAngle = Convert.ToString(computeAngle(baseFaceNormalDouble[0], baseFaceNormalDouble[1], baseFaceNormalDouble[2], joiningFaceNormalDouble[0], joiningFaceNormalDouble[1], joiningFaceNormalDouble[1]));

                                    XmlElement EWOElement = XMLFile.CreateElement("EWO");
                                    WeldmentElement.AppendChild(EWOElement);
                                    EWOElement.SetAttribute("No", Convert.ToString(ewoindex));
                                    EWOElement.SetAttribute("StartPose", StartPose);
                                    EWOElement.SetAttribute("EndPose", EndPose);
                                    EWOElement.SetAttribute("SeamAngle", Convert.ToString(weldment_angle));
                                    EWOElement.SetAttribute("Stickout", Convert.ToString(stickout));
                                    EWOElement.SetAttribute("BasePlateThickness", "0");
                                    EWOElement.SetAttribute("JoiningPlateThickness", "0");

                                    ewoindex += 1;
                                    index_normals += 1;

                                }



                                //Double mid_point_x = (weldment_curve_tessalation_points[5] - weldment_curve_tessalation_points[2]) / 2;
                                //Double mid_point_y = (weldment_curve_tessalation_points[4] - weldment_curve_tessalation_points[1]) / 2;
                                //Double mid_point_z = (weldment_curve_tessalation_points[3] - weldment_curve_tessalation_points[0]) / 2;

                                //angle_check = computeAngle(mid_point_x, mid_point_y, mid_point_z, tessalation_normals_x[tessalation_normals_x.Length - 1], tessalation_normals_y[tessalation_normals_x.Length - 1], tessalation_normals_z[tessalation_normals_x.Length - 1]);
                                //mid_point_x = (weldment_curve_tessalation_points[8] - weldment_curve_tessalation_points[5]) / 2;
                                //mid_point_y = (weldment_curve_tessalation_points[7] - weldment_curve_tessalation_points[4]) / 2;
                                //mid_point_z = (weldment_curve_tessalation_points[6] - weldment_curve_tessalation_points[3]) / 2;
                                //angle_check = computeAngle(mid_point_x, mid_point_y, mid_point_z, tessalation_normals_x[tessalation_normals_x.Length - 2], tessalation_normals_y[tessalation_normals_x.Length - 2], tessalation_normals_z[tessalation_normals_x.Length - 2]);

                            }

                        }
                        else
                        {

                            String SeamAngle = Convert.ToString(computeAngle(BaseNormal[0], BaseNormal[1], BaseNormal[2], JoiningNormal[0], JoiningNormal[1], JoiningNormal[2]));

                            MathNet.Numerics.LinearAlgebra.Vector<Double> WeldmentNormal = (BaseNormal + JoiningNormal) / 2;


                            MathNet.Numerics.LinearAlgebra.Vector<Double> ApproachVector = WeldmentNormal * 0.3;

                            WeldmentNormal = WeldmentNormal * stickout;

                            Matrix <Double> A = DenseMatrix.OfArray(new double[,] {
                                { weldment_start_coordinates[0], weldment_end_coordinates[0], (weldment_end_coordinates[0] + weldment_start_coordinates[0]) / 2, weldment_start_coordinates[0] + JoiningNormal[0] * 0.01,  weldment_end_coordinates[0] + JoiningNormal[0] * 0.01, (weldment_end_coordinates[0] + weldment_start_coordinates[0]) / 2 + JoiningNormal[0] * 0.01},
                                { weldment_start_coordinates[1], weldment_end_coordinates[1], (weldment_end_coordinates[1] + weldment_start_coordinates[1]) / 2, weldment_start_coordinates[1] + JoiningNormal[1] * 0.01,  weldment_end_coordinates[1] + JoiningNormal[1] * 0.01, (weldment_end_coordinates[1] + weldment_start_coordinates[1]) / 2 + JoiningNormal[1] * 0.01 },
                                { weldment_start_coordinates[2], weldment_end_coordinates[2], (weldment_end_coordinates[2] + weldment_start_coordinates[2]) / 2, weldment_start_coordinates[2] + JoiningNormal[2] * 0.01,  weldment_end_coordinates[2] + JoiningNormal[2] * 0.01, (weldment_end_coordinates[2] + weldment_start_coordinates[2]) / 2 + JoiningNormal[2] * 0.01 }
                            });

                            Matrix<Double> A_aux = DenseMatrix.OfArray(new double[,] {
                                { weldment_start_coordinates[0], weldment_end_coordinates[0], (weldment_end_coordinates[0] + weldment_start_coordinates[0]) / 2, weldment_start_coordinates[0] + JoiningNormal[0] * 0.01,  weldment_end_coordinates[0] + JoiningNormal[0] * 0.01, (weldment_end_coordinates[0] + weldment_start_coordinates[0]) / 2 + JoiningNormal[0] * 0.01},
                                { weldment_start_coordinates[1], weldment_end_coordinates[1], (weldment_end_coordinates[1] + weldment_start_coordinates[1]) / 2, weldment_start_coordinates[1] + JoiningNormal[1] * 0.01,  weldment_end_coordinates[1] + JoiningNormal[1] * 0.01, (weldment_end_coordinates[1] + weldment_start_coordinates[1]) / 2 + JoiningNormal[1] * 0.01 },
                                { weldment_start_coordinates[2], weldment_end_coordinates[2], (weldment_end_coordinates[2] + weldment_start_coordinates[2]) / 2, weldment_start_coordinates[2] + JoiningNormal[2] * 0.01,  weldment_end_coordinates[2] + JoiningNormal[2] * 0.01, (weldment_end_coordinates[2] + weldment_start_coordinates[2]) / 2 + JoiningNormal[2] * 0.01 }
                            });

                            Matrix<Double> B = DenseMatrix.OfArray(new double[,]
                            {
                                { 0, 0, 0, 1, 1, 1},
                                { 0, Distance(A.Column(0), A.Column(1)), Distance(A.Column(0), A.Column(2)), 0, Distance(A.Column(0), A.Column(1)), Distance(A.Column(0), A.Column(2)) },
                                { 0, 0, 0, 0, 0, 0}
                            });

                            Matrix<Double> B_aux = DenseMatrix.OfArray(new double[,]
                            {
                                { 0, 0, 0, 1, 1, 1},
                                { 0, Distance(A.Column(0), A.Column(1)), Distance(A.Column(0), A.Column(2)), 0, Distance(A.Column(0), A.Column(1)), Distance(A.Column(0), A.Column(2)) },
                                { 0, 0, 0, 0, 0, 0}
                            });

                            

                            Matrix<Double> R = AbsoluteOrientation(A_aux, B_aux);

                            

                            String StartPose = "[[" + Convert.ToString(R[0, 0]) + ", " + Convert.ToString(R[0, 1]) + ", " + Convert.ToString(R[0, 2]) + ", " + Convert.ToString(A[0, 0]) + "], [" +
                                Convert.ToString(R[1, 0]) + ", " + Convert.ToString(R[1, 1]) + ", " + Convert.ToString(R[1, 2]) + ", " + Convert.ToString(A[1, 0]) + "], [" +
                                Convert.ToString(R[2, 0]) + ", " + Convert.ToString(R[2, 1]) + ", " + Convert.ToString(R[2, 2]) + ", " + Convert.ToString(A[2, 0]) + "], " +
                                "[0, 0, 0, 1]]";

                            String EndPose = "[[" + Convert.ToString(R[0, 0]) + ", " + Convert.ToString(R[0, 1]) + ", " + Convert.ToString(R[0, 2]) + ", " + Convert.ToString(A[0, 1]) + "], [" +
                                Convert.ToString(R[1, 0]) + ", " + Convert.ToString(R[1, 1]) + ", " + Convert.ToString(R[1, 2]) + ", " + Convert.ToString(A[1, 1]) + "], [" +
                                Convert.ToString(R[2, 0]) + ", " + Convert.ToString(R[2, 1]) + ", " + Convert.ToString(R[2, 2]) + ", " + Convert.ToString(A[2, 1]) + "], " +
                                "[0, 0, 0, 1]]";


                            if (!added_approach)
                            {

                                

                                String ApproachPose = "[[" + Convert.ToString(R[0, 0]) + ", " + Convert.ToString(R[0, 1]) + ", " + Convert.ToString(R[0, 2]) + ", " + Convert.ToString((A[0, 5] + ApproachVector[0])) + "], [" +
                                Convert.ToString(R[1, 0]) + ", " + Convert.ToString(R[1, 1]) + ", " + Convert.ToString(R[1, 2]) + ", " + Convert.ToString(A[1, 5] + ApproachVector[1]) + "], [" +
                                Convert.ToString(R[2, 0]) + ", " + Convert.ToString(R[2, 1]) + ", " + Convert.ToString(R[2, 2]) + ", " + Convert.ToString(A[2, 5] + ApproachVector[2]) + "], " +
                                "[0, 0, 0, 1]]";

                                String DeparturePose = "[[" + Convert.ToString(R[0, 0]) + ", " + Convert.ToString(R[0, 1]) + ", " + Convert.ToString(R[0, 2]) + ", " + Convert.ToString((A[0, 5] + ApproachVector[0])) + "], [" +
                                Convert.ToString(R[1, 0]) + ", " + Convert.ToString(R[1, 1]) + ", " + Convert.ToString(R[1, 2]) + ", " + Convert.ToString(A[1, 5] + ApproachVector[1]) + "], [" +
                                Convert.ToString(R[2, 0]) + ", " + Convert.ToString(R[2, 1]) + ", " + Convert.ToString(R[2, 2]) + ", " + Convert.ToString(A[2, 5] + ApproachVector[2]) + "], " +
                                "[0, 0, 0, 1]]";


                                WeldmentElement.SetAttribute("ApproachPose", ApproachPose);
                                
                                WeldmentElement.SetAttribute("DeparturePose", DeparturePose);
                                
                                added_approach = true;

                            }

                            XmlElement EWOElement = XMLFile.CreateElement("EWO");
                            WeldmentElement.AppendChild(EWOElement);
                            EWOElement.SetAttribute("No", Convert.ToString(ewoindex));
                            EWOElement.SetAttribute("StartPose", StartPose);
                            EWOElement.SetAttribute("EndPose", EndPose);                            
                            EWOElement.SetAttribute("SeamAngle", SeamAngle);
                            EWOElement.SetAttribute("Stickout", Convert.ToString(stickout));
                            EWOElement.SetAttribute("BasePlateThickness", "10");
                            EWOElement.SetAttribute("JoiningPlateThickness", "10");

                            ewoindex += 1;
                        }

                    }

                    weld_index++;
                }
                String directory = System.AppDomain.CurrentDomain.BaseDirectory;
                String settings_path = Path.Combine(directory, "EWOsettings.xml");
                bool path_exists = File.Exists(settings_path);
                String file_path;
                if (!path_exists)
                {
                    System.Windows.Forms.SaveFileDialog saveFileDialog1;
                    saveFileDialog1 = new System.Windows.Forms.SaveFileDialog();
                    var path = System.Environment.GetFolderPath(System.Environment.SpecialFolder.MyDocuments);
                    saveFileDialog1.InitialDirectory = path;
                    saveFileDialog1.Title = "Select where to save the EWO file";
                    saveFileDialog1.DefaultExt = "xml";
                    saveFileDialog1.Filter = "xml files (*.xml)|*.xml|All files (*.*)|*.*";
                    saveFileDialog1.FileName = "EWO";
                    XmlDocument XMLSettings = new XmlDocument();
                    XmlElement SettingsRootNode = XMLSettings.CreateElement("EWO_File");
                    XMLSettings.AppendChild(SettingsRootNode);
                    if (saveFileDialog1.ShowDialog() == System.Windows.Forms.DialogResult.OK)
                    {
                        file_path = saveFileDialog1.FileName;
                        XMLFile.Save(file_path);
                        SettingsRootNode.SetAttribute("Address", file_path);
                        XMLSettings.Save(settings_path);
                    }
                }
                else
                {

                    XmlDocument reader = new XmlDocument();
                    reader.Load(settings_path);

                    XmlNodeList element = reader.GetElementsByTagName("EWO_File");
                    String address = element.Item(0).Attributes["Address"].Value;
                    XMLFile.Save(address);
                    

                }
                //XMLFile.Save("C:\\Users\\simat\\OneDrive\\Desktop\\ElementaryOperations\\weldments.xml");
                //XMLLog.Save("C:\\Users\\simat\\OneDrive\\Desktop\\ElementaryOperations\\weldments_log.xml");
                //var path = System.Environment.GetFolderPath(System.Environment.SpecialFolder.MyDocuments);
                //var subFolderPath = Path.Combine(path, "EWO");
                //var filePath = Path.Combine(path, "EWO", "weldments.xml");
                //System.IO.Directory.CreateDirectory(subFolderPath);
                //XMLFile.Save(filePath);
                //Debug.Print("HOPA");

            }

        }

    }

}
