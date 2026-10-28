---
layout: single
title: "IMLS 연구 2 — 접촉 barrier, Full·GN Hessian과 자동미분"
date: 2026-10-28 12:00:00 +0900
last_modified_at: 2026-10-28
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/02-barrier-hessian/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 검증 기준일: 2026-10-28

## 목적: 같은 에너지에 서로 다른 Hessian을 쓰기

Stage 3는 고정 표면 앞의 질의점 하나를 목표점에 연결하는 스프링 문제다. 접촉 에너지와 힘은 유지하면서 Full, GN, filtered GN이 어떤 Hessian을 사용하는지 비교한다. 현재 방법은 scalar gap을 합성한 단순 근사이며 GIPC 전체 구현이 아니다.

## barrier 에너지와 힘은 무엇인가?

간격을 $s>0$, 활성 거리 임계값을 $\hat d$라고 두면 `fem2d/contact.py`의 에너지는

$$
b(s)=\begin{cases}
-(s-\hat d)^2\log(s/\hat d),&0<s<\hat d,\\
0,&s\ge\hat d.
\end{cases}
$$

$s\le0$은 허용하지 않아 에너지를 무한대로 처리한다. $\hat d$는 커널 반경 $R$과 다르다. $R$은 표면 함수의 지지 영역, $\hat d$는 접촉 에너지의 활성 영역이다.

$$
\nabla E_c=b'(s)\nabla s,\qquad F_c=-\nabla E_c.
$$

힘은 에너지가 감소하는 방향이다. 수평 평면의 $s=y$에서는 $F_y=-b'(y)>0$이므로 위쪽 접촉력이 생긴다.

## Full과 GN은 어떤 항이 다른가?

chain rule을 한 번 더 적용하면

$$
H_{\rm full}=b''(s)\nabla s\nabla s^T+b'(s)\nabla^2s,
\qquad H_{\rm GN}=b''(s)\nabla s\nabla s^T.
$$

GN은 gap의 곡률 $\nabla^2s$가 들어가는 두 번째 항을 생략한다. 첫 항은 $b''\ge0$이면 양의 준정부호지만, 접촉 Hessian만으로 모든 자유도에서 양의 정부호인 것은 아니다. 스프링·탄성 강성과 경계조건이 함께 작용한다. Full은 정확한 2차 미분이라는 의미이며 매 반복에서 더 빠른 수렴을 보장하지 않는다.

현재 filtered GN은 $\rho=0.1$, $s_{\rm thr}=\rho\hat d$로 두고

$$
H_{\rm filtered}=b''\big(\max(s,s_{\rm thr})\big)\nabla s\nabla s^T
$$

를 쓴다. 작은 gap에서 계수가 커지는 것을 제한한다. 에너지와 gradient는 원래 $b(s),b'(s)$를 유지한다. 따라서 새로운 평형식을 만드는 것이 아니라 탐색 방향에 쓰는 곡률 근사를 바꾼다. 필터가 활성화되지 않으면 GN과 같다.

## 유한차분의 h와 자동미분은 어떻게 다른가?

`fn` 인자는 미분할 **함수**다. 예를 들어 `gap_fn(state)`가 scalar gap을 반환한다. `h`는 유한차분에서 입력을 움직이는 작은 양이다.

$$
\frac{\partial f}{\partial x_j}\approx\frac{f(x+he_j)-f(x-he_j)}{2h}.
$$

$h$가 크면 절단 오차가, 지나치게 작으면 뺄셈 소거와 부동소수점 오차가 문제가 된다. 현재 gradient 기본값은 $10^{-5}$, Hessian 기본값은 $10^{-4}$다. 모든 문제에 최적인 상수는 아니다.

자동미분은 알맞은 $h$를 찾아주는 방식이 아니다. 실행한 계산 그래프의 연산별 도함수를 chain rule로 연결한다. 이산화 오차나 부동소수점 오차까지 없어지는 것은 아니다. 미분 경로 안에서 Tensor를 NumPy나 Python `float`로 변환하면 연결이 끊길 수 있으므로 Torch RIMLS에서는 가중치와 보간도 Tensor로 유지한다. 제곱 잔차를 직접 사용해 `abs`·norm 표현의 0 근처 미분 문제도 피했다.

## Stage 3의 전체 문제와 실행

$$
\Pi(q)=\frac{k}{2}\|q-q_{\rm target}\|^2+b(s(q)).
$$

```bash
python fem2d/run_stage3.py
python tests/test_plane_contact_reference.py
```

첫 명령은 평면과 오각형 모서리의 세 방법을 실행한다. 두 번째는 평면 $s=y$와 barrier 식 자체를 Torch로 미분한 독립 기준을 사용한다. $k=3$, 목표 $y=-0.05$, $\hat d=0.4$에서 이분법 기준해는

$$
y^*=0.173808700231,\qquad F_y=0.671426100693.
$$

2026-10-28 재실행에서 세 방법은 초기 높이 0.02·0.12 모두 이 기준에 일치했다. 평면에서는 $\nabla^2s=0$이므로 Full과 GN이 수치 오차를 제외하면 같은 것이 정상이다. 이를 곡면의 근사 정확성 증거로 확대할 수는 없다.

## GIPC와 현재 코드의 차이

Huang, Chitalu, Lin, Komura의 [GIPC: Fast and Stable Gauss-Newton Optimization of IPC Barrier Energy (2024)](https://doi.org/10.1145/3643028)는 접촉 기하의 재정식화, 해석적 eigensystem, filtering·stiffening, nearly parallel edge-edge 처리와 GPU 구현 등을 다룬다. 논문 첫 부분과 기여 설명을 확인했다. 현재 scalar-gap GN과 계수 제한은 그 아이디어에서 출발한 공부용 실험이다. 논문의 전체 방법이나 보고된 가속을 재현했다고 말할 수 없다. 다음 단계는 FEM 강성과 접촉 Hessian을 한 에너지 안에서 결합하는 것이다.

이번 `run_stage3.py` 재실행에서 평면은 세 방법 모두 4회 업데이트로 수렴했다. 오각형 모서리도 모두 5회로 수렴했고, 위치는 약 `(0.435376, 1.230267)`, gap은 `0.1989185`, 에너지는 `0.0673452`였다. 이 한 예제의 결과가 일반 곡면 전체의 수렴 보장은 아니다.
