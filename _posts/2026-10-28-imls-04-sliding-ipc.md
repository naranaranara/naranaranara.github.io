---
layout: single
title: "IMLS 연구 4 — 평면 슬라이딩과 공식 IPC 비교 읽기"
date: 2026-10-28 12:00:00 +0900
last_modified_at: 2026-10-28
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/04-sliding-ipc/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 검증 기준일: 2026-10-28

## 목적: 에너지와 힘을 구분해 비교하기

평면 위에서 접선으로 움직일 때 물리적인 마찰을 넣지 않았다면 접선 접촉력은 0이어야 한다. 고정 높이의 질의점 스캔과, 하중 아래 정적 평형을 매번 푸는 변형 블록 실험을 구분해 살펴본다.

## 고정 평면 스캔에서 에너지가 0인가?

`fem2d/sliding_plane.py`는 높이를 0.12로 고정하고 $x\in[-1,1]$을 스캔한다. 스프링 평형이나 시간 적분을 푸는 실험이 아니다. 평면이면 $s=y$이므로 $E_c=b(y)$는 x에 무관한 **0이 아닌 상수**다.

```bash
python fem2d/sliding_plane.py
```

2026-10-28 재실행에서 $\hat d=0.4$일 때 에너지는 `0.0943914679`, $F_y$는 `1.3275581038`이었다. 샘플 9·17·33개에 대한 최대 $|F_x|$는 약 $3.56\times10^{-18}$ 이하로 부동소수점 오차 수준이다. $F_x\approx0$은 x에 따른 미분이 0이라는 뜻이지 에너지 자체가 0이라는 뜻은 아니다.

## 변형 블록은 어떻게 움직이나?

`block_sliding.py`는 상단 x좌표를 지정하고 아래로 전체 하중 0.3을 가한다. 아래 경계의 두 점 Gauss 적분으로 접촉 에너지를 더한다. 스프링은 없다. 각 지정 위치마다 Newton으로 정적 평형을 찾고 다음 위치의 초기값에 앞선 해를 사용한다.

```bash
python fem2d/block_sliding.py --compare-ipc --show
```

`--show`는 화면에만 표시하며 그림 파일을 저장하지 않는다. 수치 파일이 필요할 때만 `--outdir`을 지정한다. 기본 비교는 이동 0~0.4의 **5개 위치**다. 별도 IPC 드라이버의 기본값과 혼동하지 않도록 위 진입점을 기준으로 삼는다.

평면 gap은 적분점의 y좌표이고 절점 위치에 선형이므로 gap Hessian은 0이다. Full과 GN은 겹칠 수 있고 필터 임계값보다 멀면 filtered GN도 같다. 같은 곡선이 나온다고 세 코드를 실행하지 않았다는 뜻은 아니다.

## 공식 IPC는 어떻게 분리했나?

실제 비교 구현은 `ipc/test_ipc_comparison.py`, `ipc/test_ipc_worker.py`, 의존성은 `ipc/requirements-ipc.txt`에 분리했다. 공식 `ipctk==1.3.0`을 지속적인 별도 subprocess에서 사용해 이 환경의 Torch·native 라이브러리 충돌을 피한다. worker는 공식 collision constraints, barrier 미분, PSD projection과 `compute_collision_free_stepsize`를 호출한다. 블록과 고정 바닥 사이 접촉만 허용하며 self-contact와 friction은 제외한다.

IPC는 공식 CCD를 사용하고 RIMLS는 후보 상태의 양의 gap·면적·지지 영역을 확인한다. RIMLS도 IPC와 같은 CCD 보장을 갖는다고 쓰면 안 된다.

초기 상태에서 RIMLS와 IPC의 전체 법선력을 맞추어 IPC scale을 한 번 정한 뒤 고정했다. 기본 저장 결과의 scale은 `0.8218807474310905`다. 이 보정은 초기 법선력 한 조건을 맞추며 두 모델의 에너지 정의와 적분을 동일하게 만드는 것은 아니다. IPC는 traditional discrete formulation을 선택했다. 다른 IPC formulation 전반에 대한 결론으로 확대할 수 없다.

## 기본 5개 위치의 실제 결과

기존 `block_sliding_ipc_results/summary.csv`·`metadata.json`을 읽고 같은 명령을 그림 없이 재실행했다. 네 방법 모두 5개 위치에서 수렴했다.

| 모델 | 최대 $\lvert F_x\rvert$ | 전체 $F_y$ | 관찰 |
|---|---:|---:|---|
| IPC | 약 0.06606 | 약 0.3 | 위치에 따라 수평력이 바뀜 |
| RIMLS Full | 약 $2.48\times10^{-18}$ | 약 0.3 | 수치 오차 수준의 수평력 |
| RIMLS GN / filtered GN | 약 $4.14\times10^{-18}$ | 약 0.3 | 평면에서 거의 동일한 결과 |

IPC 블록의 기울어짐은 앞선 결과 그림에서 관찰되었지만 절점별 힘과 모멘트를 분해해 원인을 확정한 상태는 아니다. 전체 수평력만으로 특정 절점이 회전의 원인이라고 단정하지 않는다.

## steps, 바닥 샘플 수, FEM nx는 어떻게 다른가?

41개 위치는 40개 이동 구간이다. 바닥 점 13개는 12개 선분을 만든다. `nx`는 블록의 가로 FEM 절점 수이며 현재 메시에는 가로 2행이 있다. 셋은 서로 다른 해상도다.

```bash
python ipc/test_ipc_comparison.py --increments 41 --max-displacement 0.4 --plane-samples 13
```

기존 `floor13_steps41`, `floor25_steps41`, `floor49_steps41`의 CSV에서 IPC만 집계한 값은 다음과 같다. 이번에 CSV를 재집계한 결과이며 41개 위치 실험 전체를 새로 실행한 것은 아니다.

| 바닥 점 / 선분 | 위치 수 | 접촉 에너지 max−min | 샘플링된 최대 $\lvert F_x\rvert$ |
|---|---:|---:|---:|
| 13 / 12 | 41 | 0.00437456070716 | 0.0829900864934 |
| 25 / 24 | 41 | 0.00156514416796 | 0.0526077330563 |
| 49 / 48 | 41 | 0.0000447674708353 | 0.00111667492402 |

이 데이터에서는 바닥 세분화에 따라 에너지 변동과 샘플링된 수평력 최대값이 함께 감소했다. 하지만 에너지 진폭 감소만으로 그 미분인 힘의 감소를 일반적으로 보장할 수는 없다. 작은 진폭이어도 더 짧은 길이에서 바뀌면 기울기가 클 수 있다. 또한 표의 힘은 전체 평형 경로를 따라 얻은 에너지 곡선의 수치 미분과 동일한 양이라고 가정하지 않는다.

위치를 촘촘하게 평가하면 곡선의 샘플링이 개선된다. 물리 모델 자체에 smoothing을 추가한 것은 아니다. 41개 위치에서 측정한 최대값도 연속 위치 전체의 엄밀한 최대값은 아니다.

## 한계와 다음 단계

절대 접촉 에너지는 정의·적분 방식이 달라 우열 비교의 기준이 아니다. C++ IPC와 Python/Torch 계산, 프로세스 통신 비용이 섞인 시간도 공정한 속도 benchmark가 아니다. 현재 결과는 특정 2D 평면 문제의 관찰이다. 연구 아이디어의 노벨티나 GIPC 대비 속도 우월성은 아직 입증하지 않았다. 다음에는 독립 기준해로 조립 정확성을 확인하고 일반 곡면과 3D 접촉은 별도로 검증해야 한다.

배경은 Du et al.의 [smooth surface contact 논문 (2024)](https://doi.org/10.1111/cgf.15187) §3–4에서 확인했다. 공식 도구의 출처는 [IPC Toolkit 저장소](https://github.com/ipc-sim/ipc-toolkit)다.
