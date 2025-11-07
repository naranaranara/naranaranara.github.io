---
layout: post
title:  "Paper Rebuild 시리즈 #3: 더 정교하게 하기"
toc: true          # 우측 목차
toc_sticky: true   # 스크롤해도 고정
date: 2025-11-6
tags: [markdown, blog]
use_math: true
mathjax: true
render_with_liquid: false
---
# Paper Rebuild 시리즈 #3:값 더 정교하게 -> 실패

## 평가점을 안쪽으로 
1 \. 현재 평가점이 경계와 너무 가까워서 값이 이상하다고 생각해서 eps=$0.25h_x$ 로 정의

- 결과: **6.91189**에서 **6.90041** 로 값이 더 줄어들었다.
$\rightarrow$ 표면에서 한 칸 안쪽을 찍었으므로 자연스러운 반응
<br>

2 \. 정확도 개선 시도
- 2차함수, 점 3개로 계산하여 보려고 함
```
const unsigned int order    = 2;
const unsigned int quadRule = 3; 
```
$\rightarrow$ 계산이 매우 느려지면서 wsl로 열었던 우분투가 응답없음 상태 (일단 보류)  

## 제대로 열유속이 적용되었는지 확인
1 \. 가열 패치에서 열유속($q_0$)가 정확한 면적에, 정확한 크기로 적용됐는지 확인하기 위해서 적용 유속을 면적분 한 값을 확인  
$ \rightarrow \int (-k \nabla Tn) ,\ dA = q_0(2a)2=4 $

- #include 아래, mian() 함수 위쪽에 추가

```
double integrate_solution_flux(const dealii::DoFHandler<3>& dof,
                               const dealii::FiniteElement<3>& fe,
                               const dealii::Vector<double>& D_steady,
                               double k /*=385.0*/)
{
  using namespace dealii;
  const unsigned int face_q = 2;
  QGauss<2> qf_face(face_q);
  FEFaceValues<3> fe_face(fe, qf_face,
                          update_gradients | update_normal_vectors | update_JxW_values);

  double flux = 0.0;
  std::vector<Tensor<1,3>> grads(qf_face.size());

  for (auto cell : dof.active_cell_iterators())
    for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f){
      if (!cell->face(f)->at_boundary()) continue;
      if (cell->face(f)->boundary_id()!=2) continue; 
      fe_face.reinit(cell, f);

      fe_face.get_function_gradients(D_steady, grads); // 면 위 ∇T
      for (unsigned int q=0; q<qf_face.size(); ++q){
        const Tensor<1,3> n = fe_face.normal_vector(q);
        const double qn = -k * grads[q] * n; 
        flux += qn * fe_face.JxW(q); 
      }
    }
  return flux; 
}
```
- solver_steady() 바로 뒤에 추가
  
```
double q_from_sol = integrate_solution_flux(problem.dof_handler, problem.fe,
                                            problem.D_steady, /*k=*/385.0);
std::cout << "[check] solution_flux (∫-k∇T·n dA) = " << q_from_sol
          << "  (target ≈ 4.0)\n";
```
2 \. 결과: 항상 sum(F)=4로 정확히 나와서 이 문제는 아닌 걸로 판단

## 외삽법 써보기
1 \. 개념  
  경계 정확히 위를 찍으면 패치 윗면은 기울기가 급해서 경계 위에 잡으면 값이 부정확할 수 있다. 따라서 2점 외삽법을 사용하기로 한다. 2점 외삽법은 경계에서 아주 작게 떨어진 두점을 잡고 직선으로 0까지 만드는 방법이다
  $$
  T(0) \approx 2T(\Delta) -T(2\Delta)
  $$  
  <br>
2 \. 방법  
  값이 바껴도 적용될 수 있도록 비율을 사용한다. (0.25h,205h) 두 점을 쓴다.

- main.cc #include 바로 뒤

```
double probe_T_at(const dealii::DoFHandler<3>& dof,
                  const dealii::Vector<double>& D,
                  const dealii::Point<3>& p)
{
  return dealii::VectorTools::point_value(dof, D, p);
}

double boundary_limit_T_center(const dealii::DoFHandler<3>& dof,
                               const dealii::Vector<double>& D,
                               double eps1, double eps2)
{
  using dealii::Point;
  const double T1 = probe_T_at(dof, D, Point<3>(eps1, 1.0, 1.0)); // (x=eps1, y=1, z=1)
  const double T2 = probe_T_at(dof, D, Point<3>(eps2, 1.0, 1.0)); // (x=eps2, y=1, z=1)
  return T1 + (T1 - T2) * (eps1 / (eps2 - eps1));
}
```

- 호출 위치는 main에서 solve_steady()다음

```
const double k = 385.0, q0 = 1.0, L = 1.0;
auto Ttilde_L = [&](double T){ return -(k/(q0*L)) * T; };

const double h_min = dealii::GridTools::minimal_cell_diameter(problem.triangulation);
const double eps1  = 0.25 * h_min;
const double eps2  = 0.50 * h_min;

{
  dealii::Point<3> p1(eps1, 1.0, 1.0);
  double T = probe_T_at(problem.dof_handler, problem.D_steady, p1);
  std::cout << "[QoI] center @ eps=0.25h  Ttilde_L = " << Ttilde_L(T) << "\n";
}

{
  double T_eps0 = boundary_limit_T_center(problem.dof_handler, problem.D_steady, eps1, eps2);
  std::cout << "[QoI] center (eps->0)  Ttilde_L = " << Ttilde_L(T_eps0) << "\n";
}
```
3\.결과 : 7.55711으로 목표값 대비 +3.148% $\rightarrow$ 3점 외삽법 쓰면 더 정확해지려나?
## 3점 외삽법
1\. 결과: 7.9로 오차가 더욱 커지는 것을 확인했다.  
<br>
2\. 고찰: 3점 외삽은 2차 함수(곡성)으로 예측하는 것이다. 만약에 데이터가 부드러운 곡선이라면 3점 외삽법이 더 정확할 수 있다. 띠리사 해석 결과가 완벽하게 매끄럽지 않고 노이즈(울퉁불퉁함)이 섞여있을 수 있다는 것을 유추할 수 있었다. $\rightarrow$ 그렇다면 해상도를 더 높여서 mesh를 25,50,50으로 해볼까?

## 매쉬 더 잘게 자르기  
1\. 해상도를 높이는 방법을 선택했다. 더 정확한 값이면 이론값과 비슷해질거라고 생각했다. 그리고 25,50,50으로 나눈 이유는 이때까지만 해도 그냥 x와 y,z가 2배만 되면 다 괜찮다고 생각했다.
```
    std::vector<unsigned int> num_of_elems(3);
    num_of_elems[0] = 25;  // Nx
    num_of_elems[1] = 50;  // Ny
    num_of_elems[2] = 50;  // Nz
```
2\. 결과가 7.4정도로 나와서 정확한 값인줄 알고 paraview에 플롯팅 해봤다.  
<img width="383" height="205" alt="first_cliff" src="https://github.com/user-attachments/assets/15379cdd-88b1-4941-a29c-5d9ddf344825" />


중간에 값이 없어서 절벽마냥 뚝 끊긴 그래프가 생겼다. $\rightarrow$ hanging node 문제 발생
