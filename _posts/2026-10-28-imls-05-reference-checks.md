---
layout: single
title: "IMLS 연구 5 — 독립 기준해로 2D FEM과 접촉 검증하기"
date: 2026-10-28 12:00:00 +0900
last_modified_at: 2026-10-28
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/05-reference-checks/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 검증 기준일: 2026-10-28

## 목적: 코드끼리의 일치보다 독립적인 답 찾기

같은 조립 코드를 두 번 호출해 값이 같다는 것만으로는 공통 오류를 찾기 어렵다. 탄성 단독, 접촉 단독, 결합 문제를 나누어 독립적인 해석식 또는 scalar 평형식과 비교했다. 아래는 2026-10-28 현재 코드를 재실행한 결과다.

## 1. 2D 평면응력 압축

`tests/test_fem_compression.py`는 폭 $W=2$, 높이 $H=0.5$, 두께 $t=1$, $E=30$, $\nu=0.3$, 압력 $p=0.03$인 직사각형을 사용한다. 바닥의 $u_y$만 고정하고 바닥 왼쪽 한 점의 $u_x$만 고정해 Poisson 팽창을 허용한다. 바닥 전체 x를 고정하는 문제와 다르다.

상단 균일 traction은 선분 양 끝에 $p\,t\,L/2$씩 consistent load로 분배한다. 정확한 단축응력 변위장은 좌하단을 기준으로

$$
u_x=\nu\frac{p}{E}x,\qquad u_y=-\frac{p}{E}y.
$$

$$
u_y^{\rm top}=-0.0005,\quad u_x^{\rm right}=0.0006,\quad
R_y=pWt=0.06,\quad E_e=\frac{p^2}{2E}WHt=1.5\times10^{-5}.
$$

```bash
python tests/test_fem_compression.py
```

기본 조건의 최대 변위 오차는 $2.168\times10^{-19}$, 반력과 에너지도 기준에 일치했다. 테스트는 `nx=2,3,7`과 $(\nu,t)=(0.3,1),(0,0.7),(0.45,2)$의 9조건을 포함하며 모두 통과했다.

## 2. 평면과 스프링의 접촉

`tests/test_plane_contact_reference.py`는 정확한 평면 $s=y$를 쓰고 barrier 식 자체를 Torch로 미분한다. 기준은 접촉 코드의 `barrier_first_second`를 호출하지 않는다. 평형 조건

$$
k(y-y_{\rm target})+b'(y)=0
$$

을 구간 $(0,\hat d)$의 이분법으로 푼다. $k=3$, $y_{\rm target}=-0.05$, $\hat d=0.4$에서

$$
y^*=0.173808700231,\qquad -b'(y^*)=0.671426100693.
$$

```bash
python tests/test_plane_contact_reference.py
```

3개 테스트가 통과했다. 접촉장 검사는 샘플 수 9·17·33, 여러 x·y, 세 Hessian 방법을 포함한다. Stage 3 평형은 세 방법과 두 초기값에서 기준에 일치했다. filtered GN의 Hessian은 필터 계수에 대응하는 기준과 비교하며 Full Hessian과 같다고 요구하지 않는다.

## 3. FEM과 경계 적분 접촉의 결합

`tests/test_fem_contact_compression.py`는 2D FEM·RIMLS·Gauss 적분·기존 solver를 결합한다. 슬라이딩 BC를 바꾸어 수평 한 점만 고정하고 바닥 수직변위는 자유롭게 둔다. 상단은 압력 0.03의 consistent load다. 단위 두께에서 독립적인 접촉 평형은

$$
-b'(s^*)=p,\qquad s^*=0.158606842714\quad(\hat d=0.2).
$$

이 scalar 식을 이분법으로 풀고, 강체적인 상하 이동과 탄성 압축을 합쳐 기준 배치를 만든다. 바닥 gap과 블록 높이 변화는 별개의 양이다. 높이 변화는 $-pH/E=-0.0005$, 위쪽 접촉력 합계는 $pW=0.06$이다.

기준 배치 자체의 잔차·gap·힘·전체 에너지를 먼저 확인하고, 초기 상태에서 solver가 그 배치로 수렴하는지 검사한다.

```bash
python tests/test_fem_contact_compression.py
```

2개 테스트가 통과했다. solver 검사는 `(nx, initial_gap)=(3,0.12),(5,0.12),(3,0.01)`과 세 방법의 9조합이다. 최대 좌표 오차는 약 $8.07\times10^{-12}$였다. 초기 gap 0.01은 필터 임계값 $0.1\hat d=0.02$보다 작아 filtered GN이 실제로 한 번 이상 활성화되는지도 확인한다. 필터가 꺼진 상황에서 세 방법이 같은 결과를 내는 것만 검사한 것은 아니다.

## 검증 범위와 다음 단계

이 검사는 평면·균일 변형에 대한 요소 강성, 하중, 반력, 접촉 미분, 경계 적분과 비선형 평형의 일관성을 뒷받침한다. 선형 요소가 정확히 표현할 수 있는 아핀 변위장이므로 작은 오차가 나오는 이유도 분명하다.

일반 곡면의 거리 근사, 법선 변화, 대변형, 연속 경로 충돌 회피, self-contact·friction, 성능 우위를 증명하는 결과는 아니다. 다음에는 접촉을 붙이기 전에 3D 메시와 FEM을 독립적으로 검증한다.
