---
layout: single
title: "IMLS 연구 1 — 정오각형 샘플에서 IMLS·RIMLS 표면까지"
date: 2026-10-28 12:00:00 +0900
last_modified_at: 2026-10-28
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/01-geometry-imls/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 검증 기준일: 2026-10-28

## 목적과 핵심 질문

표면을 선분·삼각형으로 직접 나누면 요소 경계에서 법선이나 접촉 조합이 바뀔 수 있다. 위치와 법선 샘플로 암시적 함수를 만들고, 그 함수를 접촉 간격으로 사용하는 준비부터 시작했다. 문답은 학습 중 다룬 질문을 주제별로 재구성한 설명이며 대화의 직접 인용은 아니다.

## 정오각형 좌표는 왜 cos와 sin인가?

`fem2d/geometry.py`의 핵심은 다음과 같다.

```python
angles = np.pi / 2 + np.arange(5) * 2 * np.pi / 5
vertices = radius * np.column_stack((np.cos(angles), np.sin(angles)))
edges = np.roll(vertices, -1, axis=0) - vertices
edges_unit = edges / np.linalg.norm(edges, axis=1, keepdims=True)
edge_normals = np.column_stack((edges_unit[:, 1], -edges_unit[:, 0]))
```

반지름 $r$인 원의 점은 $(r\cos\theta,r\sin\theta)$이다. $\pi/2$에서 시작하면 첫 꼭짓점이 위에 놓이고, $2\pi/5$씩 증가시키면 다섯 꼭짓점이 반시계 방향으로 생긴다. `np.arange(5)`는 0부터 4까지이므로 시작점을 중복하지 않는다.

`np.roll(vertices, -1, axis=0)`은 각 행에 다음 꼭짓점을 대응시키고 마지막에는 첫 꼭짓점을 붙인다. 따라서 `edges[i]`는 $v_{i+1}-v_i$이다. 반대로 빼면 변 방향이 뒤집힌다. 변 길이로 나눈 벡터가 단위 접선 $t$이고, 반시계 다각형의 외향 법선은 오른쪽 회전 $n=(t_y,-t_x)$이다. 열린 경계에서는 변 순서와 선택한 접촉 면에 따라 회전 부호가 달라지므로 방향부터 확인해야 한다.

## 꼭짓점 법선과 변 위 샘플은 어떻게 만드나?

모서리에는 인접한 두 변 법선이 있다. 이 구현은 두 법선을 더하고 정규화한 평균 방향을 꼭짓점 법선으로 정한다. 이는 샘플 법선의 선택이며 원래 다각형의 모서리가 유일한 법선을 갖게 된다는 뜻은 아니다.

변 위 위치는 $p_{i,j}=v_i+t(v_{i+1}-v_i)$, $t=j/n_{\rm per\_edge}$이다. `range(n_per_edge)`에서 $t=1$은 제외한다. 그 끝점은 다음 변의 $t=0$으로 들어가므로 중복을 피한다. 변마다 5점을 만들면 전체는 25점이다.

| 이름 | 역할 |
|---|---|
| 기하 꼭짓점 | 원래 다각형의 모양을 정의한다. |
| 표면 샘플 | IMLS에 위치와 법선을 제공한다. 꼭짓점과 변 내부점을 포함할 수 있다. |
| FEM 절점 | 변위 자유도를 갖고 요소를 연결한다. |
| 접촉 적분점 | 경계 에너지를 적분한다. 위치는 FEM 절점에서 보간한다. |

샘플을 늘린다고 FEM 자유도가 자동으로 늘어나지는 않는다.

## 수평선과 커널의 의미는?

수평선에서는 $p_i=(x_i,0)$, $n_i=(0,1)$이다. 질의점 $q$와 샘플 $p_i$의 거리로 공간 가중치를 정한다.

$$
t_i=\max\left(1-\frac{\|q-p_i\|^2}{R^2},0\right),\qquad \phi_i=t_i^4.
$$

$R$은 영향을 주는 이웃의 반경이다. $R$ 밖의 점은 가중치가 0이다. 여기의 $t_i$는 무차원 커널 값이며 선분 보간 비율과 같은 변수가 아니다. 반경 안에 유효한 샘플이 없어 분모가 거의 0이면 함수를 평가할 수 없다.

샘플의 접평면까지의 부호 있는 값은 $s_i=n_i^T(q-p_i)$이다. IMLS는 이를 공간 가중 평균한다.

$$
f(q)=\frac{\sum_i\phi_i(q)s_i(q)}{\sum_i\phi_i(q)}.
$$

$f=0$이 암시적 표면이고 법선이 향하는 쪽을 양의 간격으로 사용한다. 일반 곡면에서 $f$를 정확한 최단 유클리드 거리라고 가정해서는 안 된다. 모든 샘플이 같은 수평선에 있으면 $s_i=y$이므로 분모가 유효한 영역에서 $f(q)=y$가 정확히 성립한다. 가중치가 달라져도 같은 값의 평균이기 때문이다.

## RIMLS는 모서리를 어떻게 다루나?

현재 `fem2d/rimls.py`는 기본 IMLS 값 $f_0$와 기울기 $g_0$로 한 번 재가중한다.

$$
w_i=\phi_i\exp\left[-\frac{(f_0-s_i)^2}{(\sigma_rR)^2}\right]
\exp\left[-\frac{\|g_0-n_i\|^2}{\sigma_n^2}\right],\qquad
f_r=\frac{\sum_i w_i s_i}{\sum_i w_i}.
$$

공간 가중치는 가까운 샘플을, 잔차 가중치는 현재 추정과 잘 맞는 접평면을, 법선 가중치는 현재 기울기와 비슷한 법선을 강조한다. 서로 다른 면의 법선을 무조건 섞는 효과를 줄일 수 있지만, 어떤 모서리든 원래 각도로 정확히 복원한다는 보장은 없다. $R,\sigma_r,\sigma_n$과 샘플 법선에 영향을 받는다. 현재 코드는 반복 robust fitting 전체가 아니라 위 한 번의 재가중 형태다.

## 구현과 확인

NumPy 구현은 값을 확인하는 데 사용하고 Torch 구현은 이후 접촉 미분 그래프를 연결한다. 코드는 `fem2d/geometry.py`, `imls.py`, `rimls.py`에 나뉜다. 아래 명령은 연구 코드 루트에서 실행한다. 블로그 저장소에 연구 코드 본체를 복사한 것은 아니다.

```bash
python fem2d/sliding_plane.py
```

2026-10-28 재실행에서 샘플 9·17·33개, 높이 0.12의 평면 스캔은 IMLS와 RIMLS 모두 통과했다. 최대 간격 오차는 약 $2.78\times10^{-17}$, 최대 $|F_x|$는 약 $3.56\times10^{-18}$이었다. 이 검사는 평면 재현과 계산 일관성을 확인한다. 다음에는 이 간격 함수를 barrier 에너지에 넣어 힘과 Hessian을 구한다.

## 참고와 범위

Du, Li, Coros, Thomaszewski의 [Robust and Artefact-Free Deformable Contact with Smooth Surface Representations (2024)](https://doi.org/10.1111/cgf.15187), §4.1, 식 (4)–(6)을 확인했다. 논문의 IMLS 접촉 아이디어와 현재 공부용 RIMLS 재가중 구현은 구분한다. 논문의 self-contact와 CCD는 이 단계에 구현하지 않았다.
