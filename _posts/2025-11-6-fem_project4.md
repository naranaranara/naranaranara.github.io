---
layout: post
title:  "Paper Rebuild 시리즈 #4: haning node 원인: 국부정련"
toc: true          # 우측 목차
toc_sticky: true   # 스크롤해도 고정
date: 2025-11-6
tags: [markdown, blog]
use_math: true
mathjax: true
render_with_liquid: false
---
## 국부 정련 끈 그래프끼리 비교
<img width="400" height="357" alt="nohanging_(20,40,40)vs(24,48,48)" src="https://github.com/user-attachments/assets/256fd86e-6f0e-496a-ba8c-9fc7b8ec4eb1" />  

1\. haning node는 국부 정련을 했을 때, 큰 셀의 모서리/변 중간에만 걸려 있는(맞은편 셀 꼭짓점과 안 맞는) ‘매달린’ 절점이다.
   
2\. 그래프 해석: 그래프가 끊기는 현상은 발생하지 않았지만 값이 이론값과 많이 달라졌고 x=1로 가면서의 값이 두 경우 모두 비슷한 값이여야 하지만 아예 다른 경향을 보인다. $\rightarrow$ 뭔가 잘못됨.

## ubuntu를 멀티 부팅으로 다운
<img width="325" height="391" alt="coarsettilde" src="https://github.com/user-attachments/assets/10461b1f-0af6-4fa8-a984-2e1096423b9a" />  

1 \. 멀티 부팅으로 ubuntu 다운 받으니까 계산 속도가 빨라졌다. 그래서 30,60,60을 넣었는데 매쉬가 딱 맞게 들어가서 그런지 haning node가 발생하지 않았다. $\rightarrow$ Ny=60이면 hy=20/60=1/3 $\rightarrow$ y<3이 딱 9칸(정수)라 경계가 깔끔하게 떨어져서 못 읽는 격자가 없어졌다.  

<img width="400" height="357" alt="haning_(20,40,40)vs(30,60,60)" src="https://github.com/user-attachments/assets/cea697e0-2a3d-4b99-8d2f-30aa24af071c" />  
2\. 한 그래프에서 한번에 비교해봤다. 근데 수렴되는 값이 서로 달랐다. x=1일 때 (패치 바로 뒷면)에서는 단열 조건이므로 비슷한 값으로 수렴해야한다. &\rightarrow$ T+C가 있어도 미분하면 상수항이 없어지는 문제가 있다. 따라서 메시/정련에 따라서 몇 개가 빠지거나  더 걸리거나 하면 곡선이 평행이동하는 현상이 생긴다. 따라서 id를 사용해서 상수 오프셋을 없애주면 전체 face의 모든 dof가 항상 동일하게 고정된다. 매쉬 수가 달라져도 같은 값으로 수렴할 수 있다.

<img width="400 " height="357" alt="bc_edit_(20,40,40)vs(30,60,60)" src="https://github.com/user-attachments/assets/7e8462d1-a411-4006-a4f7-5ea56cdac23a" />  
3\. 최종 결과 그래프


  - 경계를 id로 저장

```
template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements)
{
  const double x_min = 0.0,  x_max = 1.0;   // L=1
  const double y_min = 0.0,  y_max = 20.0;  // b=20
  const double z_min = 0.0,  z_max = 20.0;

  Point<dim,double> min(x_min,y_min,z_min), max(x_max,y_max,z_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);

  constexpr types::boundary_id PATCH= 2;  
  constexpr types::boundary_id YMAX = 3;  
  constexpr types::boundary_id ZMAX = 4;   
  const double tol = 1e-12;

  for (auto cell : triangulation.active_cell_iterators()){
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
      if (!cell->face(f)->at_boundary()) continue;
      const auto c = cell->face(f)->center();

      if (std::abs(c[1]-y_max) < tol){
        cell->face(f)->set_boundary_id(YMAX);
        continue;
      }
      if (std::abs(c[2]-z_max) < tol){
        cell->face(f)->set_boundary_id(ZMAX);
        continue;
      }

      if (std::abs(c[0]-x_min) < tol){
        if (c[1] >= 0.0 && c[1] <= 2.0 &&
            c[2] >= 0.0 && c[2] <= 2.0){
          cell->face(f)->set_boundary_id(PATCH); // 지정 열유속 q0
        }
        continue;
      }
    }
  }
  unsigned faces_ymax=0, faces_zmax=0, faces_patch=0;
  for (auto &cell : triangulation.active_cell_iterators()){
    for (unsigned f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
      if (!cell->face(f)->at_boundary()) continue;
      const auto bid = cell->face(f)->boundary_id();
      if (bid==YMAX)  ++faces_ymax;
      if (bid==ZMAX)  ++faces_zmax;
      if (bid==PATCH) ++faces_patch;
    }
  }
  std::cout << "[mesh] faces YMAX=" << faces_ymax
        << ", ZMAX=" << faces_zmax
        << ", PATCH=" << faces_patch << std::endl;
}
```
  - dirichlet b.c.

```
template <int dim>
void FEM<dim>::define_boundary_conds()
{
  boundary_values_of_D.clear();
  boundary_values_of_V.clear();

  Functions::ZeroFunction<dim> zero; 

  constexpr types::boundary_id YMAX = 3;
  constexpr types::boundary_id ZMAX = 4;

  VectorTools::interpolate_boundary_values(dof_handler, YMAX, zero, boundary_values_of_D);
  VectorTools::interpolate_boundary_values(dof_handler, ZMAX, zero, boundary_values_of_D);

  boundary_values_of_V = boundary_values_of_D;

  std::cout << "[bc] Dirichlet DOFs = " << boundary_values_of_D.size() << std::endl;

  std::set<types::boundary_id> ids = {YMAX, ZMAX};
  const ComponentMask cmask(fe.n_components(), true);

  std::vector<bool> on_id(dof_handler.n_dofs());
  DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), on_id, ids);

  const unsigned n_fixed = std::count(on_id.begin(), on_id.end(), true);
  std::cout << "[bc] Dirichlet DOFs (check) = " << n_fixed << std::endl;
}
```
  - 확인 코드: dirichlet dofs(check)와 dirichlet dofs의 수가 값게 출력되어야 한다.

```
std::cout << "[mesh] faces YMAX=" << faces_ymax
          << ", ZMAX=" << faces_zmax
          << ", PATCH=" << faces_patch << std::endl;

std::cout << "[bc] Dirichlet DOFs = " << boundary_values_of_D.size() << std::endl;
std::cout << "[bc] Dirichlet DOFs (check) = " << n_fixed << std::endl;
```
3\. 결과

  - (20,40,40)
  
    [mesh] faces YMAX=800, ZMAX=800, PATCH=16
    [bc] Dirichlet DOFs = 1701  
    [bc] Dirichlet DOFs (check) = 1701  
    Number of active elems:       41772  
    Number of degrees of freedom: 46489  
    sum(F) = -4  
    Ttilde_L(eps->0) = 7.54319  
    **오차 2.9568558%**
  
- (30,60,60)
  
  [mesh] faces YMAX=1800, ZMAX=1800, PATCH=36  
  [bc] Dirichlet DOFs = 3751 
  [bc] Dirichlet DOFs (check) = 3751  
  Number of active elems:       140109  
  Number of degrees of freedom: 150523 
  sum(F) = -4  
  Ttilde_L(eps->0) = 7.34
  **오차: 0.2830205%**
4\. 고찰
  
  생각보다 wsl과 멀티부팅 우분투의 계산 속도 차이는 확연했다. wsl로 할 떄는 격자수가 조금만 커져도 reconnecting이 뜨면서 연결이 불안정했는데 멀티부팅으로 우분투를 여니까 격자수가 커져도 버텼다. 만약에 계속 wsl로 했으면 이 프로젝트는 못 끝냈을 것 같다. 그리고 값이 나왔어도 paraview를 통해서 그래프를 그려보고 설정한 경계조건이 제대로 들어갔는지 나온 값이 유의미한 결과인지 확인해보는 과정도 중요하다. haning node 현상을 처음 봤는데 격자수가 많으면 많을수록 좋긴 하지만 정수값으로 딱 맞아 떨어지는지 확인하고 격자수를 늘리는게 중요한 것 같다. 상수 오프셋 되는 현상은 면을 고정해서 상수 자유도를 잡아주면 매쉬를 바꿔도 값이 비슷하게 떨어진다. 만약에 두 케이스의 값이 너무 다르면 의심해봐야 한다.
  
# 최종 코드
- fem.h

```
#ifndef FEM4_H
#define FEM4_H
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <vector>
#include <set>
#include <algorithm>

using namespace dealii;

const unsigned int order    = 1;  
const unsigned int quadRule = 2;  
template <int dim>
class FEM
{
public:
  FEM (double Alpha);
  ~FEM();

  void generate_mesh(std::vector<unsigned int> numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
  void assemble_neumann_rhs_steady(double q0);
  void solve_steady();
  void apply_initial_conditions();
  void solve_trans();
  void output_steady_results();
  void output_trans_results(unsigned int index);
  double l2norm();

  Triangulation<dim> triangulation;
  FESystem<dim>      fe;
  DoFHandler<dim>    dof_handler;
  QGauss<dim>        quadrature_formula;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> M, K, system_matrix;
  Vector<double>       D_steady, D_trans, V_trans, F, RHS;

  Table<2,double>      nodeLocation; 
  std::map<unsigned int,double> boundary_values_of_D;
  std::map<unsigned int,double> boundary_values_of_V;

  std::vector<double> l2norm_results;
  double              alpha;

  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    nodal_data_component_interpretation;
};

template <int dim>
FEM<dim>::FEM (double Alpha)
: fe (FE_Q<dim>(order), 1),
  dof_handler (triangulation),
  quadrature_formula(quadRule)
{
  alpha = Alpha;
  nodal_solution_names.push_back("D");
  nodal_data_component_interpretation.push_back(
    DataComponentInterpretation::component_is_scalar);
}

template <int dim>
FEM<dim>::~FEM () { dof_handler.clear (); }

template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements)
{
  const double x_min = 0.0,  x_max = 1.0;  
  const double y_min = 0.0,  y_max = 20.0; 
  const double z_min = 0.0,  z_max = 20.0;

  Point<dim,double> min(x_min,y_min,z_min), max(x_max,y_max,z_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);

  constexpr types::boundary_id PATCH= 2; 
  constexpr types::boundary_id YMAX = 3;  
  constexpr types::boundary_id ZMAX = 4;   
  const double tol = 1e-12;

  for (auto cell : triangulation.active_cell_iterators()){
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
      if (!cell->face(f)->at_boundary()) continue;
      const auto c = cell->face(f)->center();

      if (std::abs(c[1]-y_max) < tol){
        cell->face(f)->set_boundary_id(YMAX);
        continue;
      }
      if (std::abs(c[2]-z_max) < tol){
        cell->face(f)->set_boundary_id(ZMAX);
        continue;
      }

      if (std::abs(c[0]-x_min) < tol){
        if (c[1] >= 0.0 && c[1] <= 2.0 &&
            c[2] >= 0.0 && c[2] <= 2.0){
          cell->face(f)->set_boundary_id(PATCH); // 지정 열유속 q0
        }

        continue;
      }

    }
  }
  unsigned faces_ymax=0, faces_zmax=0, faces_patch=0;
  for (auto &cell : triangulation.active_cell_iterators()){
    for (unsigned f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
      if (!cell->face(f)->at_boundary()) continue;
      const auto bid = cell->face(f)->boundary_id();
      if (bid==YMAX)  ++faces_ymax;
      if (bid==ZMAX)  ++faces_zmax;
      if (bid==PATCH) ++faces_patch;
    }
  }
  std::cout << "[mesh] faces YMAX=" << faces_ymax
        << ", ZMAX=" << faces_zmax
        << ", PATCH=" << faces_patch << std::endl;
}

template <int dim>
void FEM<dim>::define_boundary_conds()
{
  boundary_values_of_D.clear();
  boundary_values_of_V.clear();

  Functions::ZeroFunction<dim> zero;

  constexpr types::boundary_id YMAX = 3;
  constexpr types::boundary_id ZMAX = 4;

  VectorTools::interpolate_boundary_values(dof_handler, YMAX, zero, boundary_values_of_D);
  VectorTools::interpolate_boundary_values(dof_handler, ZMAX, zero, boundary_values_of_D);

  boundary_values_of_V = boundary_values_of_D;

  std::cout << "[bc] Dirichlet DOFs = " << boundary_values_of_D.size() << std::endl;

  std::set<types::boundary_id> ids = {YMAX, ZMAX};
  const ComponentMask cmask(fe.n_components(), true);

  std::vector<bool> on_id(dof_handler.n_dofs());
  DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), on_id, ids);

  const unsigned n_fixed = std::count(on_id.begin(), on_id.end(), true);
  std::cout << "[bc] Dirichlet DOFs (check) = " << n_fixed << std::endl;
}

template <int dim>
void FEM<dim>::assemble_neumann_rhs_steady(const double q0)
{
  QGauss<dim-1> qf_face(quadRule);
  FEFaceValues<dim> fe_face(fe, qf_face, update_values | update_JxW_values);

  const unsigned int dofs_per_elem = fe.dofs_per_cell;
  std::vector<types::global_dof_index> local(dofs_per_elem);
  Vector<double> Flocal(dofs_per_elem);

  for (auto cell : dof_handler.active_cell_iterators()){
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
      if (!cell->face(f)->at_boundary()) continue;
      if (cell->face(f)->boundary_id() != 2) continue; 

      fe_face.reinit(cell, f);
      Flocal = 0.0;

      for (unsigned int q=0; q<qf_face.size(); ++q){
        const double JxW = fe_face.JxW(q);
        for (unsigned int A=0; A<dofs_per_elem; ++A){
          const double NA = fe_face.shape_value(A, q);
          Flocal[A] -= q0 * NA * JxW;
        }
      }

      cell->get_dof_indices(local);
      for (unsigned int A=0; A<dofs_per_elem; ++A)
        F[local[A]] += Flocal[A];
    }
  }
}

template <int dim>
void FEM<dim>::setup_system()
{
  dof_handler.distribute_dofs (fe);

  // DoF -> coordinates
  MappingQ1<dim,dim> mapping;
  std::vector<Point<dim,double>> dof_coords(dof_handler.n_dofs());
  nodeLocation.reinit(dof_handler.n_dofs(),dim);
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++)
    for(unsigned int j=0; j<dim; j++)
      nodeLocation[i][j] = dof_coords[i][j];

  define_boundary_conds();

  sparsity_pattern.reinit (dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  K.reinit (sparsity_pattern);
  M.reinit (sparsity_pattern);
  system_matrix.reinit (sparsity_pattern);
  D_steady.reinit(dof_handler.n_dofs());
  D_trans.reinit(dof_handler.n_dofs());
  V_trans.reinit(dof_handler.n_dofs());
  RHS.reinit(dof_handler.n_dofs());
  F.reinit(dof_handler.n_dofs());

  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
}

template <int dim>
void FEM<dim>::assemble_system()
{
  M=0; K=0; F=0;

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients | update_JxW_values);

  const unsigned int dofs_per_elem = fe.dofs_per_cell;
  const unsigned int num_quad_pts  = quadrature_formula.size();

  FullMatrix<double> Mlocal(dofs_per_elem, dofs_per_elem);
  FullMatrix<double> Klocal(dofs_per_elem, dofs_per_elem);
  std::vector<unsigned int> local_dof_indices(dofs_per_elem);

  const double rho = 3.8151e6; 

  for (auto elem = dof_handler.begin_active(); elem != dof_handler.end(); ++elem)
  {
    elem->get_dof_indices(local_dof_indices);
    fe_values.reinit(elem);

    Mlocal = 0.0;
    Klocal = 0.0;

    for (unsigned int q=0; q<num_quad_pts; ++q){
      const double JxW = fe_values.JxW(q);
      for (unsigned int A=0; A<dofs_per_elem; ++A){
        const double NA = fe_values.shape_value(A,q);
        for (unsigned int B=0; B<dofs_per_elem; ++B){
          const double NB = fe_values.shape_value(B,q);
          Mlocal[A][B] += rho * NA * NB * JxW;
        }
      }
    }

    FullMatrix<double> kappa(dim,dim); kappa = 0.0;
    for (unsigned int d=0; d<dim; ++d) kappa[d][d] = 385.0;

    for (unsigned int q=0; q<num_quad_pts; ++q){
      const double JxW = fe_values.JxW(q);
      for (unsigned int A=0; A<dofs_per_elem; ++A){
        const Tensor<1,dim> gradA = fe_values.shape_grad(A,q);
        for (unsigned int B=0; B<dofs_per_elem; ++B){
          const Tensor<1,dim> gradB = fe_values.shape_grad(B,q);
          double accum = 0.0;
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              accum += gradA[i]*kappa[i][j]*gradB[j];
          Klocal[A][B] += accum * JxW;
        }
      }
    }

    for (unsigned int i=0; i<dofs_per_elem; ++i){
      const unsigned int I = local_dof_indices[i];
      for (unsigned int j=0; j<dofs_per_elem; ++j){
        const unsigned int J = local_dof_indices[j];
        K.add(I, J, Klocal[i][j]);
        M.add(I, J, Mlocal[i][j]);
      }
    }
  }
}

template <int dim>
void FEM<dim>::solve_steady()
{
  F = 0.0;
  const double q0 = 1.0;
  assemble_neumann_rhs_steady(q0);

  MatrixTools::apply_boundary_values (boundary_values_of_D, K, D_steady, F, false);

  double sumF = 0.0;
  for (unsigned int i=0; i<F.size(); ++i) sumF += F[i];
  std::cout << "sum(F) = " << sumF << std::endl; 

  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D_steady, F);

  output_steady_results();
}

template <int dim>
void FEM<dim>::apply_initial_conditions()
{
  const unsigned int N = dof_handler.n_dofs();
  for(unsigned int i=0; i<N; i++){
    const double x = nodeLocation[i][0];
    if(x < 0.5) D_trans[i] = 300.0;
    else        D_trans[i] = 300.0 + 20.0*(x - 0.5);
  }

  system_matrix.copy_from(M);
  K.vmult(RHS,D_trans);
  RHS *= -1.;
  RHS.add(1.,F);
  MatrixTools::apply_boundary_values (boundary_values_of_V, system_matrix, V_trans, RHS, false);

  SparseDirectUMFPACK  A;
  A.initialize(system_matrix);
  A.vmult (V_trans, RHS);

  output_trans_results(0);
  l2norm_results.push_back(l2norm());
}

template <int dim>
void FEM<dim>::solve_trans()
{
  apply_initial_conditions();

  const double delta_t = 1.0;
  const unsigned int N = dof_handler.n_dofs();
  Vector<double> D_tilde(N);

  for(unsigned int t_step=1; t_step<3001; ++t_step){
    D_tilde.equ(1.0, D_trans);
    D_tilde.add(delta_t*(1.0 - alpha), V_trans);

    system_matrix.copy_from(M);
    system_matrix.add(alpha*delta_t, K);
    K.vmult(RHS, D_tilde);
    RHS *= -1.0;
    MatrixTools::apply_boundary_values (boundary_values_of_V, system_matrix, V_trans, RHS, false);

    SparseDirectUMFPACK  A;
    A.initialize(system_matrix);
    A.vmult (V_trans, RHS);

    D_trans = D_tilde;
    D_trans.add(alpha*delta_t, V_trans);

    if(t_step%100 == 0){
      output_trans_results(t_step);
      l2norm_results.push_back(l2norm());
    }
  }
}

template <int dim>
void FEM<dim>::output_steady_results ()
{
  std::ofstream output1 ("solution.vtk");
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (D_steady, nodal_solution_names,
                            DataOut<dim>::type_dof_data,
                            nodal_data_component_interpretation);
  data_out.build_patches (); data_out.write_vtk (output1); output1.close();
}

template <int dim>
void FEM<dim>::output_trans_results (unsigned int index)
{
  char filename[100];
  std::snprintf(filename, 100, "solution_%u.vtk", index);
  std::ofstream output1 (filename);
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (D_trans, nodal_solution_names,
                            DataOut<dim>::type_dof_data,
                            nodal_data_component_interpretation);
  data_out.build_patches (); data_out.write_vtk (output1); output1.close();
}

template <int dim>
double FEM<dim>::l2norm()
{
  double l2 = 0.0;
  FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

  const unsigned int dofs_per_elem = fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  const unsigned int num_quad_pts = quadrature_formula.size();
  double u_steady = 0.0, u_trans = 0.0;

  for (auto elem = dof_handler.begin_active(); elem != dof_handler.end(); ++elem){
    elem->get_dof_indices (local_dof_indices);
    fe_values.reinit(elem);

    for(unsigned int q=0; q<num_quad_pts; ++q){
      u_steady = 0.0; u_trans = 0.0;
      for(unsigned int A=0; A<dofs_per_elem; ++A){
        const double NA = fe_values.shape_value(A,q);
        const unsigned int I = local_dof_indices[A];
        u_steady += D_steady[I] * NA;
        u_trans  += D_trans [I] * NA;
      }
      const double diff = (u_steady - u_trans);
      l2 += diff*diff * fe_values.JxW(q);
    }
  }
  return std::sqrt(l2);
}

#endif 
```
- main.cc

```
#define BOOST_ALLOW_DEPRECATED_HEADERS
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "FEM4.h"
#include <deal.II/base/point.h>
#include <deal.II/fe/mapping_q1.h>         
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/grid_tools.h>
using namespace dealii;

#define USE_LOCAL_REFINE 1

static inline double probe_T_at(const DoFHandler<3>& dof,
                                const Vector<double>& D,
                                const Point<3>& p)
{
  MappingQ1<3> mapping;
  return VectorTools::point_value(mapping, dof, D, p);
}

static inline double boundary_limit_T_center(const DoFHandler<3>& dof,
                                             const Vector<double>& D,
                                             double eps1, double eps2)
{
  const double T1 = probe_T_at(dof, D, Point<3>(eps1, 1.0, 1.0));
  const double T2 = probe_T_at(dof, D, Point<3>(eps2, 1.0, 1.0));
  return T1 + (T1 - T2) * (eps1 / (eps2 - eps1)); // linear extrapolation to eps->0
}

static inline double local_hx_at_patch_center(const Triangulation<3>& tria)
{
  const double y0 = 1.0, z0 = 1.0;
  double best = 1e300;
  Triangulation<3>::active_cell_iterator best_cell;

  for (auto cell : tria.active_cell_iterators()){
    for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f){
      if (!cell->face(f)->at_boundary()) continue;
      if (cell->face(f)->boundary_id()!=2) continue; 
      const auto c = cell->face(f)->center();
      const double d = std::hypot(c[1]-y0, c[2]-z0);
      if (d < best){ best = d; best_cell = cell; }
    }
  }

  if (best == 1e300) return GridTools::minimal_cell_diameter(tria);

  double xmin = 1e300, xmax = -1e300;
  for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_cell; ++v){
    const auto &pt = best_cell->vertex(v);
    xmin = std::min(xmin, pt[0]);
    xmax = std::max(xmax, pt[0]);
  }
  return (xmax - xmin);
}

int main(){
  try{
    deallog.depth_console(0);
    std::vector<unsigned int> num_of_elems(3);
    num_of_elems[0] = 20;  // Nx
    num_of_elems[1] = 40;  // Ny
    num_of_elems[2] = 40;  // Nz

    FEM<3> problem(/*alpha*/1.0);
    problem.generate_mesh(num_of_elems);

#if USE_LOCAL_REFINE
    for (unsigned int r=0; r<2; ++r){
      for (auto cell : problem.triangulation.active_cell_iterators()){
        const auto c = cell->center();
        if (c[0] < 0.2 && c[1] < 3.0 && c[2] < 3.0)
          cell->set_refine_flag();
      }
      problem.triangulation.execute_coarsening_and_refinement();
    }
#endif

    problem.setup_system();
    problem.assemble_system();
    problem.solve_steady(); // prints sum(F)

    const double k = 385.0, q0 = 1.0, L = 1.0;
    auto Ttilde_L = [&](double T){ return -(k/(q0*L))*T; };

    const double hx   = local_hx_at_patch_center(problem.triangulation);
    const double eps1 = 0.5*hx;
    const double eps2 = 1.5*hx;

    const double T_eps0      = boundary_limit_T_center(problem.dof_handler, problem.D_steady, eps1, eps2);
    const double Ttilde_eps0 = Ttilde_L(T_eps0);
    std::cout << "Ttilde_L(eps->0) = " << Ttilde_eps0 << std::endl;

    return 0;
  }
  catch (std::exception &exc){
    std::cerr << "\nException: " << exc.what() << std::endl; return 1;
  }
  catch (...){
    std::cerr << "\nUnknown exception!" << std::endl; return 1;
  }
}
```
  
