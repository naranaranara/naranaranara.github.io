---
layout: single
title: "IMLS 연구 11 — 여러 곡면과 변형 표면의 접촉 실험"
date: 2026-10-28 16:50:00 +0900
last_modified_at: 2026-10-28
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/11-curved-contact-cases/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 2026-10-28 추가 학습 기록

## 목적: 구면에서 본 현상을 일반화해도 되는가?

[구면 진단]({{ '/imls/10-sphere-curvature/' | relative_url }}) 뒤에 원통·타원체·오목 구면·안장면을 추가했다. 같은 `sphere_contact.solve_point`로 고정 표면 문제를 풀고, 별도로 변형되는 FEM 표면에 질의점을 결합했다. 모든 곡면이 GN에 불리하다는 전제를 두지 않았다.

## 해석적 기준장은 모두 정확한 거리인가?

`fem3d/geometry.py`는 샘플과 법선을, `curved_contact.py`는 각 장과 초기화를 구성한다.

| 표면 | 해석적 장과 샘플의 의미 |
|---|---|
| 구면 | 바깥쪽 정확한 Euclidean SDF, 반지름 1 |
| 원통 | 무한 원통의 방사 SDF; RIMLS 샘플의 z범위는 [-2,2]로 유한 |
| 타원체 | 반축 1.4·1·0.7, 바깥쪽 정확한 최단거리 |
| 오목 구면 | 안쪽이 허용 영역, $s=1-\|q\|$, 법선은 안쪽 |
| 안장면 | 정규화한 그래프 정의 함수로 정확한 Euclidean SDF는 아님 |

타원체는 최근접점의 Lagrange multiplier root를 detached 이분법으로 먼저 찾고 4회의 미분 가능한 Newton 갱신을 연결한다. 내부 trial의 음수 반환은 비허용 배치를 거절하기 위한 것이며 내부 SDF 구현으로 사용하지 않는다. 바깥쪽에서는 $\|\nabla s\|=1$, $(\nabla^2s)\nabla s=0$도 검사한다.

안장면은 $h=0.35(x^2-y^2)$에 대해

$$
s(q)=\frac{z-h(x,y)}{\sqrt{1+0.49(x^2+y^2)}}
$$

를 쓴다. 부호와 영등위면을 정의하는 함수이지만, 이를 정확한 거리 기준이라고 부르면 안 된다. 해석적 장과 RIMLS 장은 별개의 에너지이며 각 장 안에서 세 Hessian 방법을 비교한다. 같은 모델 gap 초기화도 서로 다른 좌표를 가질 수 있다.

## 고정 곡면 72조합의 결과

12조건×2개 장×3방법, 총 72조합 중 64개가 수렴했고 8개는 300회 예산에 걸렸다. 구면 중심 근처와 원통 축 근처의 GN·filtered GN이며 실패 행을 제외하지 않았다. 아래는 `curved_contact_results/summary.csv`의 **RIMLS 업데이트 횟수**다.

| 고정 표면 / 조건 | Full | GN | filtered GN |
|---|---:|---:|---:|
| 구면 / 접선 목표 | 7 | 20 | 20 |
| 구면 / 중심 근처 | 18 | 300 (예산 종료) | 300 (예산 종료) |
| 원통 / 축 방향 | 6 | 7 | 7 |
| 원통 / 둘레 방향 | 7 | 19 | 19 |
| 원통 / 축 근처 | 18 | 300 (예산 종료) | 300 (예산 종료) |
| 타원체 / x 방향 | 5 | 14 | 14 |
| 타원체 / y 방향 | 7 | 19 | 19 |
| 오목 구면 / 방사 방향 | 6 | 11 | 11 |
| 오목 구면 / 접선 목표 | 16 | 24 | 20 |
| 안장면 / x 방향 | 6 | 13 | 13 |
| 안장면 / y 방향 | 5 | 9 | 9 |
| 안장면 / 혼합 방향 | 5 | 10 | 10 |


오목 구면의 접선 목표에서는 backtracking이 Full 1회, GN 22회, filtered GN 6회였다. 보정 부호와 선 탐색이 달라질 수 있음을 보여 주지만 곡률 부호 하나만으로 전체 난이도를 결정하지는 않는다. 수렴한 고정 예제의 GN과 Full 해 차이는 대체로 $10^{-9}$ 수준이었다. 일반 비볼록 문제에서 Full 해와의 거리가 독립적인 참값 오차는 아니다.

```bash
python -B fem3d/curved_contact_sweep.py --max-iter 300 --outdir results/curved-contact
```

이 명령은 숫자 CSV·JSON을 생성하며 이미지는 저장하지 않는다. 이번 기록에서는 보존된 결과를 읽었고 sweep 전체를 다시 실행하지 않았다.

## 표면 자체가 변형되면 무엇을 추가해야 하나?

`deformable_curved_contact.py`는 3×3×3 절점 사면체 패치의 윗면을 휘게 만들고 바닥을 고정한다. 자유 FEM 54성분과 질의점 3성분을 합쳐 **57 DOF**다. 질의점 하나를 스프링 $k=3$으로 목표 `(0.4, 0.1, -0.1)`에 연결하며 초기 질의점은 `(-0.4, 0, 0.18)`이다. $R=0.9$, $\hat d=0.2$를 사용한다.

상단 8삼각형에서 24개 위치 샘플을 보간하고, 변형된 면의 면적 가중 법선→절점 법선→샘플 법선을 매번 계산한다. 위치뿐 아니라 법선까지 Torch 그래프로 미분해 질의점과 FEM 사이의 교차 Hessian 블록이 생긴다. 고정 표면처럼 법선을 상수로 두는 실험과 다르다.

| 재료 E | Full 업데이트 | GN / filtered GN | 세 방법의 결과 |
|---|---:|---:|---|
| 10 | 6 | 13 / 13 | 모두 수렴, 해 차이 약 $10^{-9}$ |
| 30 | 6 | 10 / 10 | 모두 수렴, 해 차이 약 $10^{-9}$ |

`deformable_curved_results/summary.csv`의 6행을 확인했다. 테스트는 교차 블록이 실제로 0이 아닌지, 임의 방향 에너지 2차 유한차분과 일치하는지, 세 방법의 배치가 같은 허용오차 안에 드는지 검사한다.

```bash
python -B fem3d/deformable_curved_contact.py --outdir results/deformable-curved
python -B tests/test_curved_contact.py
```

이번 자동검증에서 `test_curved_contact.py`의 5개 테스트는 통과했다. 접촉 Full Hessian의 직접 AD·유한차분, 곡률 부호, 타원체 거리 성질, 배치·scalar RIMLS 미분 일치와 변형 패치 결합을 포함한다.

## 한계와 다음 단계

변형 예제도 스프링 질의점 하나와 변형 표면의 접촉이다. 두 변형체의 다중 접촉, 대회전 재료, self-contact·마찰·CCD를 검증한 것이 아니다. 현재 양의 사면체 체적과 양의 gap은 끝 상태에서 검사한다. 고정 곡면과 변형 패치에서 반복 수 차이를 보았으므로, 다음에는 그 차이가 실제 비용 증가인지 총 실행시간으로 판단해야 한다.
