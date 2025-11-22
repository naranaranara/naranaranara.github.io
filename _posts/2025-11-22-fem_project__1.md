---
layout: post
title:  "FSI 시리즈 #1"
toc: true          # 우측 목차
toc_sticky: true   # 스크롤해도 고정
date: 2025-11-22
tags: [markdown, blog]
use_math: true
mathjax: true
---
## FSI 
유체가 구조에 미치는 영향 = 구조 표면에 작용하는 압력 또는 전단력(하중)  
즉, 유체가 구조에 영향을 미치면 바뀐 구조가 다시 그 주변 유체에 영향을 주고 그게 다시 구조에 영향을 주는 순환구조랄까,,  

이번에는 로켓 핀에 작용하는 응력을 확인해보기 위해서 직사각형 패널을 모델링하고 fsi의 효과를 테스트 해보기!  

## 구조  

\- 지배방정식 (quasi-static)  

$$
-\nabla \cdot \sigma(\mathbf{u}) = \mathbf{0} \quad \text{in } \Omega 
$$  

-헷갈리는 부분-  

Quasi-static이라고 했는데 지배방정식은 왜 static인지?  
왜냐면 순간순간은 정적 평형 상태이지만, 노이만 경계조건에서 압력에 따른 외력이 $sin(t)$로 시간에 따라 변하기 때문이다. 즉, 정적인 식을 시간에 따라 값을 바꿔가면서 여러 번 풀기 때문에 준정적이다.  

\- 구성방정식 (훅의 법칙)  

$$
\sigma_{ij} = C_{ijkl} \varepsilon_{kl}
$$  

```c++

template <int dim>
double FEM<dim>::C(unsigned int i,unsigned int j,
                   unsigned int k,unsigned int l)
{
  double E  = 7.0e4;
  double nu = 0.33;

  double lambda=(E*nu)/((1.+nu)*(1.-2.*nu));
  double mu    =E/(2.*(1.+nu));

  return lambda*(i==j)*(k==l)
       + mu*((i==k)*(j==l) + (i==l)*(j==k));
}
```

\- 경계조건  

왼쪽($x=0$) : 로켓 몸통에 붙은 곳 $\rightarrow$ 디리클레 경계조건  

$$
u=0
$$  

오른쪽($x=L_x$): 유체가 미는 면 $\rightarrow$ 노이만 경계조건(압력)  

$$
\boldsymbol{\sigma} \mathbf{n} = \mathbf{t}(t) = -p(t) \mathbf{n}
$$  

위/아래($y=0,y=L_y$): 하중 없음 $\rightarrow$ 노이만 경계조건 (걍 내비둠)  

$$
\boldsymbol{\sigma} \mathbf{n} = \mathbf{0}
$$  

구조 문제 전체:  

$$
\begin{cases} 
-\nabla \cdot \boldsymbol{\sigma}(\mathbf{u}) = \mathbf{0} & \text{in } \Omega \\
\boldsymbol{\sigma} = \mathbb{C} : \boldsymbol{\varepsilon}(\mathbf{u}) & \\
\mathbf{u} = \mathbf{0} & \text{on } \Gamma_{\text{left}} \\
\boldsymbol{\sigma} \mathbf{n} = -p(t)\mathbf{n} & \text{on } \Gamma_{\text{right}} \\
\boldsymbol{\sigma} \mathbf{n} = \mathbf{0} & \text{on } \Gamma_{\text{top,bottom}}
\end{cases}
$$  

코드:  

**\<디리클레 경계조건\>**  

```c++
template <int dim>
void FEM<dim>::define_boundary_conds()
{
  boundary_values.clear();
  const unsigned int totalDOFs = dof_handler.n_dofs();
  const double tol = 1e-10;

  for (unsigned int g = 0; g < totalDOFs; ++g)
  {
    const double x = dofLocation[g][0];

    if (use_bc1 && std::abs(x - 0.0) < tol)
    {
      boundary_values[g] = 0.0;
    }
  }
}
```

**\<노이만 경계조건\>**  

```c++
if (use_pressure_right)
    {
      for (unsigned int f=0; f < faces_per_elem; ++f)
      {
        if (elem->face(f)->at_boundary())
        {
          if (std::abs(elem->face(f)->center()[0] - Lx) < 1e-10)
          {
            fe_face_values.reinit (elem, f);

            for (unsigned int q=0; q<num_face_quad_pts; ++q)
            {
              const Tensor<1,dim> normal   = fe_face_values.normal_vector(q);
              Tensor<1,dim>       traction = -p_t * normal;

              for (unsigned int A=0; A<nodes_per_elem; ++A)
              {
                for(unsigned int i=0; i<dim; ++i)
                {
                  const unsigned int a  = dim*A + i;
                  const double Ni = fe_face_values.shape_value(a,q);

                  Flocal[a] += Ni * traction[i] * fe_face_values.JxW(q);
                }
              }
            }
          }
        }
      }
    }
```

## 유체  

\- 지배방정식  

$$
\tau_f \dot{p}(t) + p(t) = -K_v v(t) + p_{\text{ext}}(t)
$$  

$\tau_f \dot{p}(t)$ : 유체는 점성이 있으므로 압력이 갑자기 못 변하고 $\tau_f$ 스케일로 부드럽게 따라간다.  
$p(t)$ : 초기 압력  
$-K_v v(t)$ : 유체 압력이 그 움직임을 저지하는 방향(-)로 바뀜  
$p_{\text{ext}}(t)$: 외부에서 주어지는 압력  

\- 구성방정식 ($p$가 $v$와 $p_{ext}$를 어떤 형태로 따른다고 보는가?)  

$$
\dot{p} = \frac{-p - K_v v + p_{\text{ext}}}{\tau_f}
$$  

\- 초기조건  

$$
p(0)=p_0
$$  

```c++
double time  = 0.0;
double p     = 0.0;
```

코드:  

```c++
template <int dim>
void FEM<dim>::run_fsi()
{
  setup_system(); 

  const double dt    = 0.01;   
  const double T_end = 5.0;   
  const double Kv    = 2.0;  
  const double tau_f = 0.05; 
  const double P0   = 1.0e3; 
  const double fHz  = 1.0; 
  const double two_pi_f = 2.0 * 3.141592653589793 * fHz; // 2πf

  double time  = 0.0;
  double p     = 0.0;     
  double q_old = 0.0;        
  unsigned int step = 0;

  while (time <= T_end)
  {
    assemble_system(p);
    solve();

    double q = get_tip_displacement_y();  
    double v = (q - q_old) / dt;         

    const double p_ext = P0 * std::sin(2.0 * M_PI * fHz * time);

    const double alpha = dt / tau_f;
    double p_new = ( p - alpha * Kv * v + alpha * p_ext )
                   / ( 1.0 + alpha );

    p     = p_new;
    q_old = q;

    time += dt;
    ++step;
  }
}
```

**\<정리\>**  

1\. 유체 $\rightarrow$ 구조  
\- 유체쪽 1-DOF ODE 가 계산한 압력 $p(t)$를 구조한테 넘겨줌  
\- 구조 FEM 코드에서 이 $p(t)$를 오른쪽 노이만 하중으로 사용  

```c++
assemble_system(p);
solve();        
```

2\. 구조 $\rightarrow$ 유체  

\- 구조 해석 결과에서 tip 변위 $q(t)$를 뽑고  
\- 타임 스텝 차이 속도 term을 만들어서 유체 ODE에 넣어줌  

```c++
double q = get_tip_displacement_y();
double v = (q - q_old) / dt;   
```

<img width="704" height="384" alt="Gemini_Generated_Image_7c4z657c4z657c4z" src="https://github.com/user-attachments/assets/ccc347f8-f26c-4bcd-8763-3e07adf741f0" />{: .centered }   
 

## 결과  

![fsi_x](https://github.com/user-attachments/assets/6a003fd0-5b60-4369-b74b-6fba56473153){: .centered }   
  

![fsi_y](https://github.com/user-attachments/assets/65e66196-95f7-4687-86b7-32e8d7e3330e){: .centered }  

## 검증  

검증은 매트랩 코드를 통해서 했다.  

```matlab
clear; clc; close all;

target_axis = 'y';  

filename = 'fsi_result.csv'; 
col_idx_x = 6; 
col_idx_y = 7;  

Lx = 10.0; Ly = 50.0;
E  = 7.0e4; nu = 0.33;
dt = 0.01; T_end = 5.0;
Kv = 2.0; tau_f = 0.05;
P0 = 1.0e3; fHz = 1.0;

E_eff = E / (1 - nu^2);    
K_eq  = (E_eff * Ly) / Lx;   

t = 0:dt:T_end;
N = length(t);
p = zeros(1, N);      
u_x_model = zeros(1, N); % X방향 변위 (Model)
v_tip = zeros(1, N);  

for i = 1:N-1
    force = p(i) * Ly; 
    u_x_model(i) = - force / K_eq; % 압축이므로 음수
    
    if i > 1, v_tip(i) = (u_x_model(i) - u_x_model(i-1)) / dt; else, v_tip(i) = 0; end
    p_ext = P0 * sin(2 * pi * fHz * t(i));
    alpha = dt / tau_f;
    p_new = (p(i) - alpha * Kv * v_tip(i) + alpha * p_ext) / (1 + alpha);
    p(i+1) = p_new;
end
u_x_model(N) = -(p(N) * Ly) / K_eq;

if strcmpi(target_axis, 'x')
    fprintf('>> 모드: X방향 변위 검증 (압축 거동)\n');

    theory_data = u_x_model;

    target_col = col_idx_x;

    y_label_str = 'Tip Displacement u_x (mm)';
    title_str   = 'Validation: X-Direction (Compression)';
    
elseif strcmpi(target_axis, 'y')

    fprintf('>> 모드: Y방향 변위 검증 (푸아송 효과)\n');

    nu_eff = nu / (1 - nu);
    % Y 변위 = |X변위| * nu_eff * (Ly/Lx) * 0.5 (상단 노드 기준)
    u_y_model = -u_x_model * nu_eff * (Ly/Lx) * 0.5;
    
    theory_data = u_y_model;

    target_col = col_idx_y;

    y_label_str = 'Top Corner Displacement u_y (mm)';
    title_str   = 'Validation: Y-Direction (Poisson Expansion)';
else
    error('target_axis 변수는 ''x'' 또는 ''y''만 가능합니다.');
end
```

1\. x 방향 변위 ($u_0$)  

\- 기본 식 (훅의 법칙)  

$$
F = Kx \quad \rightarrow \quad x = \frac{F}{K}
$$  

\- 하중 조건  

$$
-F = P L_y
$$  

\- 스프링 상수($K$)  

$$
K = \frac{EA}{L}
$$  

- 조건: 평면변형률이여서 원래는 압축되면 y,z 축으로 변형이 생기는데 지금은 z축으로의 변형이 없다. 따라서 재료가 더 딱딱하게 느껴진다. 그래서 $E$ 대신에 유효탄성계수 $E_{eff}$를 사용한다.

$$
E_{eff} = \frac{E}{1 - \nu^2}
$$  

- 증명:

$$
\epsilon_x = \frac{1}{E}[\sigma_x - \nu(\sigma_y + \sigma_z)]
\epsilon_z = \frac{1}{E}[\sigma_z - \nu(\sigma_x + \sigma_y)] = 0 \quad \rightarrow \quad \sigma_z = \nu(\sigma_x + \sigma_y)
$$  

**가정**  
매트랩 1D 모델은 x방향으로만 스프링을 누르는 상황으로 가정한다. y방향으로 누르는 힘이 없어서 자유롭게 팽창한다.  
FEM 해석에서는 왼쪽 벽이 막혀있다고 가정했는데 매트랩에서는 뚫려있다고 가정할 수 있는 이유는 생베낭의 원리로 $x=L$ 부근에서는 y방향 응력이 없다는 가정을 적용가능하다.  

위의 가정에 따라서 $\sigma_y=0$이고 $\sigma_z=\nu \sigma_x$ 이므로,  

$$
\epsilon_x = \frac{1}{E}[\sigma_x - \nu(0 + \sigma_z)]
\sigma_z = \nu\sigma_x
E\epsilon_x = \sigma_x - \nu^2\sigma_x = (1 - \nu^2)\sigma_x \quad \rightarrow \quad \sigma_x = \left(\frac{E}{1 - \nu^2}\right) \epsilon_x  
$$

2\. y방향 변위 ($u_1$) $\rightarrow$ 푸아송의 원리  

이 경우에도 위의 경우와 같이 z축 방향으로 변형이 못 일어나니까 y축으로 더 많이 튀어나오니까 유효 포아송비를 써줘야 한다.  

$$
\nu_{eff} = \frac{\nu}{1 - \nu}
$$  

- 증명:
  $$
  \nu_{eff} = -\frac{\epsilon_y}{\epsilon_x}
\epsilon_x = \frac{\sigma_x}{E}(1 - \nu^2)
\epsilon_y = \frac{1}{E}[\sigma_y - \nu(\sigma_x + \sigma_z)] = -\frac{\nu\sigma_x}{E}(1 + \nu)
\therefore \nu_{eff} = \frac{-\frac{\nu\sigma_x}{E}(1 + \nu)}{\frac{\sigma_x}{E}(1 - \nu^2)} = \frac{\nu}{1 - \nu}
  $$

매트랩에 넣는 식은 최종 변위를 계산해주는 식을 넣어주고, 재료는 위아래 대칭으로 뚱뚱해지니까 절반을 나눠준다.

$$
\Delta L_y = \epsilon_y L_y = \left(-\nu_{eff} \frac{\Delta L_x}{L_x}\right) L_y = -\Delta L_x \cdot \nu_{eff} \cdot \frac{L_y}{L_x}
\Delta L_y = -\Delta L_x \cdot \nu_{eff} \cdot \frac{L_y}{L_x} \times \frac{1}{2}
$$

## 그래프 비교  

매트랩에서 가정할 때는 $x=L$인 지점에서는 왼쪽 모서리의 구속 효과가 작용하지 않는 것으로 하였는데 원래 모델은 ($L_x=10, L_y=50$)으로 x 방향의 길이가 y 방향보다 짧다. 따라서 $x=L$인 지점에서도 구속 효과가 작용한다. 따라서 결과 값을 비교해보면 오차가 심한 것을 확인할 수 있다. 따라서 코드 검증을 위해서 x 방향 길이를 늘려서 계산해보았을 때는 오차가 거의 없는 것을 확인할 수 있었다.  

**\<Case 1\>**  

<img width="450" height="285" alt="image" src="https://github.com/user-attachments/assets/4f51ad53-3264-48fb-9930-59150cfd0fae" />{: .centered}  

<img width="450" height="285" alt="image (1)" src="https://github.com/user-attachments/assets/69a77672-00f0-421b-96b4-609972ea340b" />{: .centered}   
  
**\<Case 2\>**  

<img width="450" height="285" alt="image (2)" src="https://github.com/user-attachments/assets/f512b51f-0273-4268-b1b4-e0225cdd91e5" />{: .centered}    

<img width="450" height="285" alt="image (3)" src="https://github.com/user-attachments/assets/84f00be2-8ac0-4996-a2c1-e454d8f6e62a" />{: .centered} 





















