---
layout: single
title: "IMLS 연구 7 — 3D 표면 적분과 FEM·접촉 결합 검증"
date: 2026-10-28 16:50:00 +0900
last_modified_at: 2026-10-28
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/07-fem3d-contact/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 2026-10-28 추가 학습 기록

## 목적: 3D FEM에 실제 접촉을 결합하기

[3D FEM 기초 편]({{ '/imls/06-fem3d/' | relative_url }}) 다음으로, 삼각형 표면 적분과 RIMLS 접촉을 연결했다. 먼저 정확한 평면 압축을 검증한 뒤 정적 슬라이딩과 공식 IPC 비교로 확장했다. 이 단계의 압축·지정변위 문제는 동적 시뮬레이션이 아니다.

## S와 B는 둘 다 형상함수인데 무엇이 다른가?

`fem3d/surface.py`는 바닥에 놓인 외부 삼각형 face를 찾고 각 면에 세 적분점을 만든다. barycentric 좌표는 $(2/3,1/6,1/6)$의 세 순열이다. 면적 $A$인 각 면에서 가중치는 $A/3$이다.

$$
q_k=\sum_i N_i(\xi_k)x_i,\qquad Q=SX,\qquad
E_c\approx\sum_k w_k b(s(q_k)),\quad w_k=A_{\rm ref}/3.
$$

S는 적분점에서 평가한 **형상함수 값**을 전체 절점 열에 배치한 보간행렬이다. 각 행의 합은 1이며 `samples = S @ nodes`로 위치를 얻는다. B는 **형상함수 미분**으로 절점 변위에서 변형률을 만든다. S가 접촉 위치를, B가 재료 변형률을 계산한다. 현재 접촉 적분은 고정된 기준 면적을 사용한다.

`St`, `pt`, `nt`의 t는 Tensor를 구분하는 이름이다. 시간 변수가 아니다. `contact.py`의 평면 RIMLS 반경 $R=0.9$는 이 실험의 선택값이며 논문에서 자동으로 정해 주는 값이 아니다. 자동미분도 이 R이나 유한차분 h를 찾아주는 알고리즘이 아니다.

## 3D 압축의 하중과 자유 DOF는?

폭 2·깊이 1·높이 0.5, 선형 등방성 $E=30,\nu=0.3$, 평면 $z=0$, 상단 압력 $p=0.03$을 사용한다. 전체 아래쪽 하중은 $0.03\times2=0.06$이다. 스프링 질의점 없이 블록과 평면이 접촉한다.

절점 i의 성분 번호는 `3*i`, `3*i+1`, `3*i+2`가 각각 x·y·z이다. 바닥 z를 고정하면 접촉 지지가 검증되지 않으므로 모든 수직 DOF를 자유롭게 둔다. 원점 ux·uy와 바닥 +x 모서리점 uy의 수평 3성분만 anchor로 지정한다. 상단 삼각형의 각 절점에는 아래 방향 힘 $pA/3$을 준다. 면적 A에 곱하는 압력 방향과 절점 분배 계수는 구분해야 한다.

## 독립 기준해와 실제 검증

평면에서는 $s=z$이고 단위 면적당 위쪽 접촉력이 $-b'(s)$다. 따라서

$$
-b'(s^*)=p,\qquad s^*=0.158606842714\quad(\hat d=0.2).
$$

이 식을 독립 이분법으로 풀고 강체 수직 이동+균일 탄성 압축의 기준 배치를 만든다. 높이 변화는 $-pH/E=-0.0005$, 접촉력 합은 $F_z=0.06$이다.

```bash
python -B tests/test_fem3d_contact_compression.py
python -B fem3d/compression.py
```

이번 전체 검증에서 해당 2개 테스트는 통과했다. 기준 배치의 잔차·에너지·적분점 Jacobian과 접촉 Hessian을 확인하고, 세 방법×초기 gap 0.12·0.01의 6조합을 풀었다. 최대 좌표 오차는 약 $8.07\times10^{-12}$였다. 작은 초기 gap에서 filtered GN의 활성화도 확인했다.

## 정적 슬라이딩은 압축과 구속이 같나?

`fem3d/block_sliding.py`는 먼저 이동량 0에서 preload하고 상단 ux를 0~0.4로 지정한다. 상단 한 점 uy를 anchor로 두고 모든 uz는 자유롭다. 저장된 `block_sliding_3d_results/summary.csv`의 세 방법×5위치는 모두 수렴했다. 접선력은 대략 $10^{-19}$ 이하이고 $F_z\approx0.06$이다.

```bash
python -B fem3d/block_sliding.py
```

상단 x 구속은 Poisson 팽창 일부를 제한한다. 따라서 자유 압축 기준해의 gap과 같아야 하는 문제가 아니다. 저장된 슬라이딩의 최소 **적분점 gap**은 약 `0.158601518`이다. 비교 드라이버가 보고하는 최소 절점 높이와도 구분한다.

## 공식 3D IPC와 비교하기

`ipc/ipc_worker.py`는 외부 삼각형 face와 중복 없는 edge로 collision mesh를 구성한다. 내부 FEM 절점은 충돌 정점에서 빼고 결과 gradient·Hessian을 전체 DOF에 다시 배치한다. 27절점 메시의 내부 절점 1개에 직접 접촉 미분이 0인지, 전체 Hessian이 81×81인지 검사했다. 내부 절점은 탄성을 통해서는 영향을 받는다.

공식 ipctk 1.3.0은 지속적인 별도 프로세스·Toolkit 1스레드로 호출한다. 블록–바닥 접촉만 허용하며 마찰·자기접촉은 제외한다. IPC는 local PSD clamp와 공식 CCD를 사용하고 RIMLS는 끝 상태 feasibility 검사다.

공유 바닥은 9×7=63정점, 96삼각형, 간격 0.5다. RIMLS도 같은 바닥 샘플을 사용한다. 미변형 기준 상태에서 전체 법선력을 한 번 맞춘 scale `0.05931098177337765`를 유지한다. 평형 뒤에 매번 재보정하지 않는다.

| 저장된 5위치 결과 | 최대 접선력 x의 절댓값 | 수직력 |
|---|---:|---:|
| IPC | 0.00543046 | 약 0.06 |
| RIMLS 세 방법 | $7.99\times10^{-20}$ 이하 | 약 0.06 |

`block_sliding_3d_ipc_results/summary.csv`의 20행은 모두 수렴했다. 공식 gradient의 유한차분, PSD, 작용·반작용, 바닥 관통 CCD 제한, 내부 절점 매핑 검사는 이번 재검증에서도 통과했다. 현재 전체 검증 상태도 함께 구분한다.

> 최종 검증: 현재 코드에서 전체 51개 테스트가 통과했다. 이 글의 긴 실행 수치는 보존된 CSV·JSON을 재확인한 값이며, 게시 작업에서 모든 실험을 다시 실행했다는 뜻은 아니다.

정식 실행 진입점은 `python -B ipc/ipc3d_comparison.py`다. 이전 공부용 test_ 접두사 대신 현재 정식 파일 이름을 사용한다.

## 범위와 다음 질문

이것은 공식 접촉 라이브러리를 사용한 자체 소변형 FEM simulator다. 원래 IPC 전체 solver의 논문 장면을 재현한 것은 아니다. 에너지 정의·적분·CCD·native 프로세스 통신과 Torch 계산이 달라 절대 에너지와 실행시간으로 방법의 우열을 정하지 않는다. 다음 질문은 지정변위 대신 동일한 힘으로 구동했을 때 어떤 운동이 나오는가이다.
