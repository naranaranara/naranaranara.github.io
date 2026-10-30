---
layout: single
title: "IMLS 연구 13 — 보폭과 모드 분석에서 Gap SR1까지"
date: 2026-10-29 23:00:00 +0900
last_modified_at: 2026-10-30
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact, optimization]
permalink: /imls/13-step-modes-and-sr1/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 2026-10-29 연구 및 학습, 2026-10-30 정리

## 어떤 곡률 정보를 저렴하게 되찾을 것인가

GN은 간격의 이차미분을 생략해 반복당 비용을 줄인다. 그러나 구 중심이나 원기둥 축 근처를 목표로 둔 단일 질의점 실험에서는 반복이 크게 늘었다. 바로 새 Hessian을 제안하기보다, 먼저 이동 길이만 잘 고르면 해결되는지, 방향별 강성 오차가 남는지 분석했다.

접촉 에너지가 $E_c=b(s(q))$일 때 전체 Hessian을 다음처럼 나눈다.

$$
H_{\mathrm{Full}}=\underbrace{H_{\mathrm{base}}+b''(s)\nabla s\nabla s^\mathsf T}_{G=H_{\mathrm{GN}}}+\underbrace{b'(s)\nabla^2s}_{C}.
$$

목표는 알고 있는 $G$를 버리는 것이 아니라, 계산비용이 큰 $C$의 유용한 효과를 저렴하게 얻는 것이다. 현재 scalar-gap GN은 GIPC 논문 전체 알고리즘과 같지 않다.

## 질문 1 보폭만 키우면 되는가

GN 방향 $p=-G^{-1}g$를 유지하고 기본 선 탐색의 시작 보폭 1보다 큰 값도 허용했다. 이어 같은 방향에서 에너지를 줄이는 보폭을 적극적으로 찾는 수치적 선상 최소화와, 이전 이동의 강성 비율을 이용하는 적응형 보폭을 비교했다.

$$
d_k=q_{k+1}-q_k,\qquad r_k=\frac{d_k^\mathsf T(g_{k+1}-g_k)}{d_k^\mathsf TG_kd_k}.
$$

작은 이동에서는 gradient 차이가 실제 Hessian 작용을 반영한다. $r_k$가 작으면 해당 방향에서 GN이 강성을 크게 예측했다는 단서다. 다만 secant는 구간 평균 정보이고, 한 방향만 측정한다. 모든 방향의 곡률 비율이라고 해석하면 안 된다.

저장된 `gn_step_results/matched_comparison.csv`에서 구 중심 근처는 확대 GN 50회, 적응형 74회, 선상 최소화 258회였고 원기둥 축 근처는 각각 52·44·186회였다. 선상 최소화는 수치 탐색 절차의 결과이며 수학적으로 정확한 전역 최적 보폭이라는 뜻은 아니다. 비싼 보폭 탐색이 항상 반복 수를 줄이지 않았다. 이 초기 실험 시간은 뒤의 통제된 7회 benchmark와 섞지 않는다.

## 질문 2 방향마다 오차가 다른가

같은 위치에서 일반화 고유값 문제를 진단했다.

$$
H_{\mathrm{Full}}v_i=\lambda_iGv_i,\qquad v_i^\mathsf TGv_j=\delta_{ij}.
$$

$G=LL^\mathsf T$인 Cholesky 분해 후 $L^{-1}H_{\mathrm{Full}}L^{-\mathsf T}$를 대칭 고유값 분해한다. $\lambda_i$는 일반 좌표의 강성 자체가 아니라 그 모드에서 GN 대비 실제 강성의 비율이다. 작은 $\lambda_i$는 GN이 해당 모드의 이동을 지나치게 억제할 수 있다는 뜻이다.

구 중심 근처에서는 약 0.0313·0.0531·0.9978, 원기둥 축 근처에서는 약 0.0499·0.9989·1.0039인 국소 상태가 관찰됐다. 모든 모드에 같은 보폭 배율을 적용해서 이런 서로 다른 비율을 동시에 맞추기는 어렵다. 이 숫자는 특정 상태의 진단이며 경로 전체에 고정된 물성값이 아니다. 초기 비볼록 상태의 음의 모드와 안정화도 별도로 다뤘다.

## 질문 3 모드 몇 개만 정확히 보정하면 되는가

정확한 Full Hessian과 고유모드를 알고 있다는 Oracle 조건에서 모드 1개·2개·전체 보정을 비교했다. 반복 수 감소는 확인했지만 Full 미분과 고유값 분해 비용을 내므로 실험 자체가 저비용 알고리즘은 아니다. 중요한 곡률 효과가 저차원일 수 있는지 확인하는 진단이다.

초기 모드 점수는 강성 오차와 gradient 활성도를 결합했다. 이후 방향 변화 크기와 에너지 개선량을 각각 고려한 점수를 시험했다.

$$
S_i^{\mathrm{direction}}=\left|\frac1{\lambda_i}-1\right||v_i^\mathsf Tg|,
\qquad S_i^E=\frac{(\lambda_i-1)^2}{\lambda_i}(v_i^\mathsf Tg)^2.
$$

첫 점수는 $G$ 노름에서 GN 방향이 얼마나 달라지는지, 두 번째는 양의 정부호 이차모델에서의 개선량과 관계가 있다. 비선형 수렴을 보장하는 공식은 아니다. 음수·0 근처 고유값에는 역수를 그대로 적용하지 않고 안정화와 대체 처리가 필요하다. 몇몇 문제에서는 점수를 바꿔도 선택 모드와 경로가 같았다. 다음 과제는 점수 자체보다 정확한 모드를 계산하지 않고 곡률 정보를 얻는 것이 됐다.

## Gap SR1은 무엇을 학습하는가

| 방법 | secant에 넣는 변화량 | 사용하는 행렬 |
|---|---|---|
| 일반 SR1 | $\Delta g$ | 전체 Hessian 근사 |
| Structured SR1 | $\Delta g-G_{\mathrm{new}}d$ | $G+B_C$ |
| Gap SR1 | $\Delta(\nabla s)$ | $G+b'(s)B_s$ |

Gap SR1은 현재 이미 계산한 $b'(s)$를 그대로 사용하고 $B_s\approx\nabla^2s$만 추정한다. Structured 방식은 현재 GN을 anchor로 쓰므로 유한 이동 중 GN 자체의 변화가 잔차에 섞일 수 있다.

$$
u=y-Bd,\qquad B_{\mathrm{new}}=B+\frac{uu^\mathsf T}{u^\mathsf Td}.
$$

분모가 너무 작으면 갱신을 건너뛴다. SR1은 음의 곡률도 보존할 수 있지만 최종 방향의 하강성을 위해 합산 행렬의 안정화가 필요하다. $B_s=0$으로 시작하므로 첫 방향은 GN이다. 현재는 밀집행렬을 사용하는 단일-gap 연구용이며 다중 접촉·제한 메모리·대규모 병렬 성능까지 입증한 것은 아니다.

## 코드 읽는 순서

`fem3d/solvers/gn_step.py`, `fem3d/diagnostics/modal_curvature.py`, `fem3d/solvers/modal_correction.py`, `fem3d/solvers/sr1.py` 순서로 보면 분석이 새 근사법으로 이어지는 이유를 볼 수 있다. solver는 방향과 반복 규칙을 구현하고 comparison 실험은 같은 문제로 여러 방법을 돌려 결과를 집계한다.

[다음: 14개 케이스와 곡률 근사 검증]({{ '/imls/14-sr1-validation/' | relative_url }})

