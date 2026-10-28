---
layout: single
title: "IMLS 연구 3 — 2D FEM, 자유 DOF와 접촉 Newton 해법"
date: 2026-10-28 12:00:00 +0900
last_modified_at: 2026-10-28
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/03-fem2d-solver/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 검증 기준일: 2026-10-28

## 목적과 두 가지 실험

Stage 4는 2D 변형 표면과 스프링에 연결된 질의점을 결합한다. 이후 `block_sliding.py`는 고정 평면 위의 변형 블록을 다룬다. 구동 방식과 자유도가 다르다.

| 항목 | Stage 4 | 블록 슬라이딩 |
|---|---|---|
| 접촉 | 질의점과 변형 표면 | 블록 아래 경계 적분점과 고정 평면 |
| 구동 | 질의점 스프링 | 상단 x 변위 지정, 외부 수직하중 |
| 스프링 | 있음 | 없음 |
| FEM | 소변형 평면응력 선형 삼각형 | 같은 종류의 FEM |

## FEM의 B, D, K와 Hessian은 어떻게 연결되나?

삼각형 절점당 $(u_x,u_y)$ 두 자유도가 있다. 선형 형상함수로 요소 내부 변위를 보간하고, 그 기울기를 조립한 $B$로 변형률을 계산한다.

$$
\epsilon=Bu_e,\qquad\sigma=D\epsilon,\qquad
D=\frac{E}{1-\nu^2}\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix}.
$$

변형률 순서는 $(\epsilon_{xx},\epsilon_{yy},\gamma_{xy})$이며 $\gamma_{xy}=2\epsilon_{xy}$이다. 요소 면적 $A$, 두께 $t$에서

$$
K_e=tA B^TDB,\qquad E_{\rm elastic}=\tfrac12u^TKu,
\quad\nabla E_{\rm elastic}=Ku,\quad\nabla^2E_{\rm elastic}=K.
$$

따라서 이 **선형 탄성 에너지의** Hessian이 강성행렬 $K$다. 전체 에너지에는 접촉 Hessian도 추가된다. 대변형 비선형 재료까지 상수 $K$ 하나로 표현되는 것은 아니다.

## 자유 DOF와 BC는 무엇인가?

변위 또는 위치가 지정된 성분을 제외하고 풀어야 하는 성분이 자유 DOF다. Stage 4 기본 메시에는 FEM 절점 6개가 있고 아래 절점 3개의 x·y가 고정된다. 위 절점의 6성분에 질의점 2성분을 더해 미지수는 8개다. FEM만 8자유도라는 뜻은 아니다.

블록 슬라이딩은 상단 절점의 x좌표만 기준 위치+지정 이동량으로 고정한다. 수직방향은 하중과 접촉 평형으로 결정한다. 검증용 압축 문제는 다른 BC를 사용하므로 결과를 섞지 않는다.

## base_terms와 전체 에너지는?

`stage4_solver.py`의 `base_terms`는 탄성+스프링을 반환하고 접촉을 별도로 더한다. 블록의 `base_terms`는 탄성-외력이 한 일이다. 일반적인 정적 문제는

$$
\Pi=E_{\rm elastic}-F_{\rm ext}^Tu+E_{\rm contact}
$$

로 쓰고 Stage 4에는 스프링 에너지가 추가된다. `base_terms`는 현재 문제에서 접촉을 제외한 부분을 모으는 인터페이스다.

## 경계 적분점은 절점과 다른가?

`surface.py`의 `boundary_quadrature`는 각 선분에 두 Gauss 점을 만든다. 기준 구간 $\xi\in[-1,1]$에서 두 점은 $\pm1/\sqrt3$이고 보간 비율은 $t=(1+\xi)/2$다.

$$
q=(1-t)x_a+tx_b,\qquad E_c\approx\sum_k w_k b(s(q_k)),\qquad w_k=L_{\rm ref}/2.
$$

$N_a=1-t,N_b=t$가 형상함수이고 행렬 $S$에 넣으면 전체 적분점 위치를 `S @ nodes`로 만든다. 적분점에 독립 자유도를 추가하는 것이 아니다. 각 점의 에너지 미분이 형상함수를 통해 절점으로 전달된다. 현재 접촉 적분은 **기준 선분 길이**를 쓰며 현재 변형된 길이를 사용하는 적분과 다르다.

## solver, FEM, run 파일은 무엇을 맡나?

`mesh.py`는 연결 관계·경계·절점, `fem.py`는 $B,D,K$와 탄성 에너지, `contact.py`는 접촉 에너지의 미분을 맡는다. solver는 이를 더해 평형을 찾고 run 파일은 예제를 구성해 실행한다. `if __name__ == '__main__':`는 파일을 직접 실행했을 때만 예제 코드를 호출한다. import만 했을 때 자동 실행되는 것을 막는다.

Newton 반복은 자유 성분에서 $H_{ff}\Delta z_f=-g_f$를 푼다. 음의 고유값 등으로 방향이 불안정할 때 최소 고유값을 확인해 대각 shift를 더한다. 이는 탐색용 행렬의 SPD 보정이며 원래 재료 모델을 바꾸는 것은 아니다.

Armijo 조건은 $\Pi(z+\alpha\Delta z)\le\Pi(z)+10^{-4}\alpha g^T\Delta z$다. $\alpha=1$부터 실패하면 절반으로 줄인다. 양의 요소 면적, 커널 지지 영역, 양의 gap도 검사한다. 현재 RIMLS는 후보 끝 상태의 feasibility 검사이고 연속 경로 전체의 CCD 보장은 없다.

## 실행과 검증

```bash
python fem2d/run_stage4.py
python -m unittest discover -s tests -p test_stage4_reference.py -v
python fem2d/block_sliding.py
python tests/test_fem_compression.py
python tests/test_fem_contact_compression.py
```

Stage 4는 float64 Torch 자동미분을 기본으로 사용한다. 2026-10-28 독립 압축 재검증에서는 FEM 단독의 상단 수직변위 $-0.0005$, 전체 반력 $0.06$을 확인했고, 결합 검증에서는 gap $0.158606842714$와 같은 수직 반력을 확인했다. 상세 BC와 오차는 다음 기준해 검증 편에 정리했다. 검증 범위는 소변형 평면·균일 압축이다. 다음에는 슬라이딩 에너지와 수평력을 구분해 읽는다.

이번 Stage 4 예제에서는 Full 4회, GN·filtered GN 5회 업데이트로 모두 수렴했다. 질의점은 약 `(0.099757, 0.168206)`, gap은 `0.1758048`, 최대 절점 이동은 약 `0.008502`였다. `test_stage4_reference.py`의 탄성 미분, 표면 변형과 힘 균형, 접촉 결합, 평형·구속 검사 4개도 통과했다.
