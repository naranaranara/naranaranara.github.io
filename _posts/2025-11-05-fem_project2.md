---
layout: post
title:  "Paper Rebuild 시리즈 #2:문제 상황에 맞는 코드로 변경"
toc: true          # 우측 목차
toc_sticky: true   # 스크롤해도 고정
date: 2025-11-5
tags: [markdown, blog]
use_math: true
mathjax: true
---
# Paper Rebuild 시리즈 #2:문제 상황에 맞는 코드로 변경

## 매쉬 & 패치 라벨링
- 도메인 $[0,1]x[0,20]x[0,20]$ 생성 후, $x=0$ (표면)에 boundardy_id=2 부여
<br>
```
GridGenerator::subdivided_hyper_rectangle(triangulation, numberOfElements,
                                          Point<dim>(0,0,0), Point<dim>(1,20,20));
for (auto cell : triangulation.active_cell_iterators()){
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
    if (!cell->face(f)->at_boundary()) continue;
    const auto c = cell->face(f)->center();
    if (std::abs(c[0]) < 1e-12){             // x=1 face
      if (c[1] >= 0.0 && c[1] <= 2.0 &&
          c[2] >= 0.0 && c[2] <= 2.0)
        cell->face(f)->set_boundary_id(2);
    }
  }
}
```

## Dirichlet B.C. 경계 설정  

- $ y=20 or z=20 $ 면에 $T=0$ 적용  

```
boundary_values_of_D.clear(); boundary_values_of_V.clear();
for (unsigned int i=0; i<dof_handler.n_dofs(); ++i){
  const double y = nodeLocation[i][1], z = nodeLocation[i][2];
  if (std::abs(y-20.0) < 1e-12 || std::abs(z-20.0) < 1e-12){
    boundary_values_of_D[i] = 0.0; 
    boundary_values_of_V[i] = 0.0; 
  }
}
```
## Neuman B.C. 면적분 
1\.$ \int _{\Gamma_N}(-q_0)N_A \, d\Gamma$  
<br>
2\.검증 방법: sum(F) 값이 -4 근처
<br>
```
QGauss<dim-1> qf_face(quadRule);
FEFaceValues<dim> fe_face(fe, qf_face, update_values | update_JxW_values);
for (auto cell : dof_handler.active_cell_iterators()){
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
    if (!cell->face(f)->at_boundary()) continue;
    if (cell->face(f)->boundary_id() != 2) continue; 
    fe_face.reinit(cell,f);
    Vector<double> Flocal(fe.dofs_per_cell); Flocal = 0.0;
    for (unsigned int q=0; q<qf_face.size(); ++q){
      const double JxW = fe_face.JxW(q);
      for (unsigned int A=0; A<fe.dofs_per_cell; ++A){
        const double NA = fe_face.shape_value(A,q);
        Flocal[A] -= q0 * NA * JxW;
      }
    }
    std::vector<types::global_dof_index> local(fe.dofs_per_cell);
    cell->get_dof_indices(local);
    for (unsigned int A=0; A<fe.dofs_per_cell; ++A)
      F[local[A]] += Flocal[A];
  }
}
```
## 측정값은 0보다 살짝 안쪽에서 평가
1\.패치부분에서 기울기가 확 바뀌기 때문에 경계에 정확히 찍으면 매시에 따라서 값이 튐  
- main.cc 파일에서 수정

```
const double k=385.0, q0=1.0, L=1.0, a=2.0, b=20.0;
auto Ttilde_L = [&](double T){ return -(k/(q0*L))*T; };
auto Ttilde_a = [&](double T){ return -(k/(q0*a))*T; };
auto Ttilde_b = [&](double T){ return -(k/(q0*b))*T; };

const double eps = 1e-8;
Point<3> p_surf_in(eps, 1.0, 1.0); 
double T_in = VectorTools::point_value(problem.dof_handler, problem.D_steady, p_surf_in);
std::cout << "[surf in] Ttilde_L=" << Ttilde_L(T_in) << "\n";

Point<3> p_surf(1.0, 1.0, 1.0);
double T_s = VectorTools::point_value(problem.dof_handler, problem.D_steady, p_surf);
std::cout << "[surf]    Ttilde_L=" << Ttilde_L(T_s) << "\n";
```
## 국부 정련 
1\.x=0 부분에서 패치 엣지의 영향이 x 방향으로 전파되는 것을 막기 위해서 패치 둘레에 얇은 띠를 정밀하게 만드는 과정  
<br>
2\.해상도(공간 정확도) 높이기  
```
for (unsigned int r=0; r<2; ++r){ 
  for (auto cell : problem.triangulation.active_cell_iterators()){
    const auto c = cell->center();
    if (c[0] > 0.9){                                       
      const double dy = std::min(c[1], std::abs(c[1]-2.0)); 
      const double dz = std::min(c[2], std::abs(c[2]-2.0)); 
      if (std::min(dy,dz) < 1.2) cell->set_refine_flag(); 
    }
  }
  problem.triangulation.execute_coarsening_and_refinement();
}
```
## 실행 방법 (빌드 -> 실행)
- ./ 뒤에는 파일명 넣으면 됨
```
cmake -S .. -B build -DDEAL_II_DIR=/usr/lib/cmake/deal.II -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./main4
```
## 실행 결과 및 고찰
1\. sum(F)= -4 확인
2\. [surface patch center] Ttilde_L ≈ 7.529  
$\rightarrow$ 논문 값보다 +2.8%  
3\.지금은 mesh를 x,y,z 방향으로 각각 20,40,40으로 나눠놨는데 이게 너무 작아서 그런건 아닐까 의심중
