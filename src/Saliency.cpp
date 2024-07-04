#include "Saliency.h"
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCurvatures.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkAutoInit.h>

#include <vtkSmartPointer.h>
#include <vtkCurvatures.h>
#include <vtkPolyDataReader.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkColorSeries.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkScalarBarActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>


std::vector<double>* Env::calculate_mean_curvature(cinolib::DrawableTrimesh<> m) {
    std::vector<std::vector<uint>> faces = m.vector_polys();
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();

    vtkNew<vtkPolyData> geometry;
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> polys;

    for (size_t i = 0; i < vertexs.size(); ++i) {
        double point[3] = { vertexs[i][0],vertexs[i][1] ,vertexs[i][2] };  
        points->InsertPoint(i, point);
    }
    for (size_t i = 0; i < faces.size(); i++){
        vtkIdType face[3] = { faces[i][0],faces[i][1] ,faces[i][2] };
        polys->InsertNextCell(3, face);
    }
    geometry->SetPoints(points);
    geometry->SetPolys(polys);

    vtkSmartPointer<vtkCurvatures> curvaturesFilter= vtkSmartPointer<vtkCurvatures>::New();
    curvaturesFilter->SetInputData(geometry);
    curvaturesFilter->SetCurvatureTypeToMean(); // 设置计算平均曲率
    curvaturesFilter->Update();

    vtkSmartPointer<vtkFloatArray> meanCurvatureArray = vtkFloatArray::SafeDownCast(curvaturesFilter->GetOutput()->GetPointData()->GetArray("Mean_Curvature"));


    std::vector<double>* meanCurvatures = new std::vector<double>;

    for (uint i = 0; i < vertexs.size(); i++) {
        double meanCurvature = meanCurvatureArray->GetValue(i);
        meanCurvatures->push_back(meanCurvature);
    }
    return meanCurvatures;
}


std::vector<double>* Env::get_saliency_score(cinolib::DrawableTrimesh<> m, unv::data_type type) {

    std::vector<std::vector<uint>> faces = m.vector_polys();
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    std::vector<cinolib::vec3d> normals = m.vector_vert_normals();
    int vertexCnt = vertexs.size();
    int FaceCount = faces.size();

    double xMin = INT_MAX, yMin = INT_MAX, zMin = INT_MAX;
    double xMax = INT_MIN, yMax = INT_MIN, zMax = INT_MIN;
    for (uint i = 0; i < vertexs.size(); i++) {
        xMin = std::min(xMin, vertexs[i][0]);
        xMax = std::max(xMax, vertexs[i][0]);
        yMin = std::min(yMin, vertexs[i][1]); 
        yMax = std::max(yMax, vertexs[i][1]);
        zMin = std::min(zMin, vertexs[i][2]); 
        zMax = std::max(zMax, vertexs[i][2]);
    }


    // Calculate each vertecies' shape operator
    MFUN::mat3* shapeOperators = NULL;
    float* vertexArea = NULL;
    shapeOperators = (MFUN::mat3*)malloc(vertexCnt * sizeof(MFUN::mat3));
    vertexArea = (float*)malloc(vertexCnt * sizeof(float));
    for (int i = 0; i < vertexCnt; i++) {
        vertexArea[i] = 0.0f;
        for (int j = 0; j < 9; j++)
            shapeOperators[i].m[j] = 0.0f;
    }
    for (int k = 0; k < FaceCount; k++) {
        // Calculate the face's area

        float faceArea = m.poly_area(k);

        for (int idx = 0; idx < 3; idx++) {
            int i = faces[k][idx];
            int j = faces[k][(idx + 1) % 3];
            // Get vertex i and j's normal vectors.
            MFUN::vec3 Ni = MFUN::vec3(normals[i][0], normals[i][1], normals[i][2]);
            MFUN::vec3 Nj = MFUN::vec3(normals[j][0], normals[j][1], normals[j][2]);
            // Get vertex i and j's location.
            MFUN::vec3 Vi = MFUN::vec3(vertexs[i][0], vertexs[i][1], vertexs[i][2]);
            MFUN::vec3 Vj = MFUN::vec3(vertexs[j][0], vertexs[j][1], vertexs[j][2]);

            // For vertex i, update the relative part of its shape operator
            MFUN::vec3 Tij = (MFUN::identity_mat3() - MFUN::wedge(Ni, Ni)) * (Vi - Vj);
            Tij = MFUN::normalise(Tij);
            float kappa_ij = 2 * dot(Ni, Vj - Vi);
            kappa_ij /= MFUN::get_squared_dist(Vi, Vj);
            // Maintain vi's shape operator
            shapeOperators[i] = shapeOperators[i] + (MFUN::wedge(Tij, Tij) * (kappa_ij * faceArea));
            vertexArea[i] += faceArea;

            // For vertex j, update the relative part of its shape operator
            MFUN::vec3 Tji = (MFUN::identity_mat3() - MFUN::wedge(Nj, Nj)) * (Vj - Vi);
            Tji = MFUN::normalise(Tji);
            float kappa_ji = 2 * dot(Nj, Vi - Vj);
            kappa_ji /= get_squared_dist(Vi, Vj);
            // Maintain vj's shape operator
            shapeOperators[j] = shapeOperators[j] + (MFUN::wedge(Tji, Tji) * (kappa_ji * faceArea));

            vertexArea[j] += faceArea;
        }
    }

    for (int i = 0; i < vertexCnt; i++) {
        shapeOperators[i] = shapeOperators[i] * (1.0f / vertexArea[i]);// * 10000000.0f;
        //print(shapeOperators[i]);
    }
    free(vertexArea);

    // Diagonalize the shape operator, and get the mean curvature
    float *meanCurvature_T = (float*)malloc(vertexCnt * sizeof(float));
    for (int k = 0; k < vertexCnt; k++) {
        MFUN::vec3 E1 = MFUN::vec3(1.0f, 0.0f, 0.0f);
        MFUN::vec3 Nk = MFUN::vec3(normals[k][0], normals[k][1], normals[k][2]);
        bool isMinus = MFUN::get_squared_dist(E1, Nk) > MFUN::get_squared_dist(E1 * (-1.0f), Nk);
        MFUN::vec3 Wk;
        // Diagnoalization by the Householder transform
        if (!isMinus)
            Wk = E1 + Nk;
        else
            Wk = E1 - Nk;
        Wk = normalise(Wk);
        MFUN::mat3 Qk = MFUN::identity_mat3() - (MFUN::wedge(Wk, Wk) * 2.0f);
        MFUN::mat3 Mk = MFUN::transpose(Qk) * shapeOperators[k] * Qk;
        // Calculate the mean curvature by M_k's trace;
        meanCurvature_T[k] = (float)(Mk.m[4] + Mk.m[8]);
    }
    free(shapeOperators);






    

    std::vector<double>* meanCurvature = new std::vector<double>;
    for (int k = 0; k < vertexCnt; k++)
        meanCurvature->push_back(meanCurvature_T[k]);
    free(meanCurvature_T);

    // Calculate the incident matrix ( as linked list )
    int* first = NULL;
    int* next = NULL;
    int* incidentVertex = NULL;


    // Calculate the mesh saliency by BFS
    float diagonalLength = sqrt((xMax - xMin) * (xMax - xMin) + (yMax - yMin) * (yMax - yMin) + (zMax - zMin) * (zMax - zMin));
    float sigma = 0.003 * diagonalLength;
    float* saliency[7];
    float maxSaliency[7];
    for (int i = 2; i <= 6; i++) {
        saliency[i] = NULL;;
        saliency[i] = (float*)malloc(vertexCnt * sizeof(float));
        maxSaliency[i] = INT_MIN;
    }

    // Labeled the vertecies whether covered or not.
    bool* used = NULL;
    used = (bool*)malloc(vertexCnt * sizeof(bool));
    for (int k = 0; k < vertexCnt; k++) {
        printf("\r 1/2 [%.2f%%]\t>", float(k * 100.0 / (vertexCnt - 1)));
        for (int j = 1; j <= k * 40 / vertexCnt; j++)
            std::cout << "▉";

        // Initialize the saliency and its local counter.
        for (int i = 2; i <= 6; i++)
            saliency[i][k] = 0.0f;
        // Initialize the saliency's Gaussian filter.
        float gaussianSigma1[7], gaussianSigma2[7], sumSigma1[7], sumSigma2[7];
        for (int i = 2; i <= 6; i++)
            gaussianSigma1[i] = gaussianSigma2[i] = 0.0f,
            sumSigma1[i] = sumSigma2[i] = 0.0f;
        // Get the current vertex's information.
        cinolib::vec3d vVec = vertexs[k];
        // Initialize the queue to find neighbourhood.
  //      memset(used, 0, sizeof(used));
        for (int i = 0; i < vertexCnt; i++)
            used[i] = false;
        std::queue<int> Q;
        Q.push(k);
        used[k] = true;
        // Frsit BFS
        while (!Q.empty()) {
            // Get the front element in the queue.
            int idx = Q.front(); Q.pop();
            cinolib::vec3d idxVec = vertexs[idx];
            // Put the next level vertecies into the queue.
            std::vector<uint> adj = m.vert_verts_link(idx);
            for (uint e = 0; e < adj.size(); e++) {
                if (used[adj[e]]) continue;
                cinolib::vec3d nextVec = vertexs[adj[e]];
                if ((vVec - nextVec).norm() <= 36 * sigma * sigma) {
                    Q.push(adj[e]);
                    used[adj[e]] = true;
                }
            }

            // Update Gaussian filter
            double dist = (vVec-idxVec).norm();
            for (int i = 2; i <= 6; i++) {
                float sigmaHere = i * i * sigma * sigma;
                if (dist <= sigmaHere) {
                    float factor = exp(-dist / (2 * sigmaHere));
                    gaussianSigma1[i] += meanCurvature->at(idx) * factor;
                    sumSigma1[i] += factor;
                }
                if (dist <= 2 * sigmaHere) {
                    float factor = exp(-dist / (8 * sigma * sigma));
                    gaussianSigma2[i] += meanCurvature->at(idx) * factor;
                    sumSigma2[i] += factor;
                }
            }
        }
        for (int i = 2; i <= 6; i++) {
            saliency[i][k] = fabs(gaussianSigma1[i] / sumSigma1[i]
                - gaussianSigma2[i] / sumSigma2[i]);
            maxSaliency[i] = std::max(maxSaliency[i], saliency[i][k]);
        }
    }

    // Second BFS and get the non-linear normailization of suppressian's saliency.
    std::vector<double>* smoothSaliency = new std::vector<double>;
    printf("\n");

    for (int k = 0; k < vertexCnt; k++) {
        printf("\r 2/2 [%.2f%%]\t>", float(k * 100.0 / (vertexCnt - 1)));
        for (int j = 1; j <= k * 40 / vertexCnt; j++)
            std::cout << "▉";
        smoothSaliency->push_back(0);
        float localMaxSaliency[7];//, localCntSaliency[7];
        for (int i = 2; i <= 6; i++)
            localMaxSaliency[i] = INT_MIN;
        // Get the current vertex's information.
        cinolib::vec3d vVec = vertexs[k];
        // Initialize the queue to find neighbourhood.
        for (int i = 0; i < vertexCnt; i++)
            used[i] = false;
        std::queue<int> Q;
        Q.push(k);
        used[k] = true;
        while (!Q.empty()) {
            // Get the front element in the queue.
            int idx = Q.front(); Q.pop();
            //aiVec = &(mesh.mVertices[idx]);
            //vec3 idxVec = vec3(aiVec.x, aiVec.y, aiVec.z);
            // Put the next level vertecies into the queue.
            std::vector<uint> adj = m.vert_verts_link(idx);
            for (uint e = 0; e < adj.size(); e++) {
                if (used[adj[e]]) continue;
                cinolib::vec3d nextVec = vertexs[adj[e]];
                if ((vVec - nextVec).norm() <= 36 * sigma * sigma) {
                    Q.push(adj[e]);
                    used[adj[e]] = true;
                }
            }
            // Update Gaussian filter
            for (int i = 2; i <= 6; i++)
                localMaxSaliency[i] = std::max(localMaxSaliency[i], saliency[i][idx]);
        }
        // Calculate the weighted saliency
        float saliencySum = 0.0f;
        for (int i = 2; i <= 6; i++) {
            float factor = (maxSaliency[i] - localMaxSaliency[i]) * (maxSaliency[i] - localMaxSaliency[i]);
            smoothSaliency->at(k) += (GLfloat)saliency[i][k] * factor;
            saliencySum += factor;
        }
        smoothSaliency->at(k) /= saliencySum;
    }
    printf("\n");
    if (type == unv::POINT_DATA) {
        return unv::normalize(smoothSaliency);
    }
    else {
        std::vector<double>* smoothSaliency_face = new std::vector<double>;
        smoothSaliency_face = unv::point2face(*smoothSaliency, m);
        delete smoothSaliency;
        return unv::normalize(smoothSaliency_face);
    }
    return smoothSaliency;
    // Clean up resources
    
}
