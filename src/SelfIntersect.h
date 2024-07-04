#pragma once
#include <fstream>
#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/ARAP.h>
#include <cinolib/gl/glcanvas.h>
#include <omp.h>
#include<vector>
#include <set>
#include <cinolib/gl/surface_mesh_controls.h>



#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
//#include <CGAL/IO/OFF_reader.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <vtkActor.h>
#include <vtkCellData.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkSTLReader.h>
#include <vtkLookupTable.h>
#include <vtkRenderWindow.h>
#include <vtkPolyDataMapper.h>
#include <vtkUnsignedCharArray.h>
#include <vtkColorTransferFunction.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

namespace Eva {
    bool SelfIntersection(const std::vector<cinolib::vec3d> vertexs, const std::vector<std::vector<uint>> faces);
    std::set<uint>* IntersectionList(const std::vector<cinolib::vec3d> vertexs, const std::vector<std::vector<uint>> faces);
}