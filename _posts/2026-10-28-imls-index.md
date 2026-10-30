---
layout: single
title: "IMLS·RIMLS 접촉 연구 — 학습 기록 인덱스"
date: 2026-10-28 12:00:00 +0900
last_modified_at: 2026-10-30
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/index/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

## 연구 목적과 현재 단계

매끄러운 IMLS/RIMLS 표면으로 표면 이산화의 인공 접촉력을 줄이고, 접촉 Hessian의 유용한 곡률 효과를 저렴하게 반영할 수 있는지 공부한다. GIPC에서 영감을 얻었지만 현재 scalar-gap GN·filtered GN은 GIPC 전체 방법이 아니다.

2026-10-29에는 보폭·고유모드 진단에서 Gap SR1으로 이어지는 분석, 14개 케이스 반복 benchmark와 오프라인 곡률 검증을 진행했다. 또 두 변형체의 다중 접촉, 희소 FEM, 독립 `src` 코드와 공식 BVH·CCD 연결을 구현했다. 마지막 질문은 CPU 기반 전체 계산의 비용과 GPU 확장 방향이었다. GPU 가속은 아직 구현·성능 검증하지 않았다.

## 학습 순서와 글

각 문답은 실제 질문 주제를 재구성한 설명이며 대화의 직접 인용은 아니다. 1–12편은 10월 28일까지의 기록이다. 과거 단회 시간·구조 설명을 최신 반복 측정이나 현재 경로와 혼동하지 않도록, 최신 상태는 13–16편과 아래 표를 기준으로 읽는다.

### 기초 1–6편

- [IMLS 연구 1 — 정오각형 샘플에서 IMLS·RIMLS 표면까지]({{ '/imls/01-geometry-imls/' | relative_url }})
- [IMLS 연구 2 — 접촉 barrier, Full·GN Hessian과 자동미분]({{ '/imls/02-barrier-hessian/' | relative_url }})
- [IMLS 연구 3 — 2D FEM, 자유 DOF와 접촉 Newton 해법]({{ '/imls/03-fem2d-solver/' | relative_url }})
- [IMLS 연구 4 — 평면 슬라이딩과 공식 IPC 비교 읽기]({{ '/imls/04-sliding-ipc/' | relative_url }})
- [IMLS 연구 5 — 독립 기준해로 2D FEM과 접촉 검증하기]({{ '/imls/05-reference-checks/' | relative_url }})
- [IMLS 연구 6 — 3D 사면체 메시와 선형 FEM 압축 검증]({{ '/imls/06-fem3d/' | relative_url }})

### 2026-10-28 3D 접촉과 곡률 7–12편

- [IMLS 연구 7 — 3D 표면 적분과 FEM·접촉 결합 검증]({{ '/imls/07-fem3d-contact/' | relative_url }})
- [IMLS 연구 8 — 지정변위에서 힘 구동 동역학으로, 힘의 부호 검증]({{ '/imls/08-force-driven-dynamics/' | relative_url }})
- [IMLS 연구 9 — 시간 간격을 줄여 무엇을 검증했나]({{ '/imls/09-time-convergence/' | relative_url }})
- [IMLS 연구 10 — 구면 곡률과 GN의 느린 접선 수렴]({{ '/imls/10-sphere-curvature/' | relative_url }})
- [IMLS 연구 11 — 여러 곡면과 변형 표면의 접촉 실험]({{ '/imls/11-curved-contact-cases/' | relative_url }})
- [IMLS 연구 12 — 반복 횟수보다 총 시간, 코드 정리와 재검증]({{ '/imls/12-runtime-and-refactor/' | relative_url }})

### 2026-10-29 SR1 검증과 두 변형체 13–16편

- [IMLS 연구 13 — 보폭과 모드 분석에서 Gap SR1까지]({{ '/imls/13-step-modes-and-sr1/' | relative_url }})
- [IMLS 연구 14 — Gap SR1의 시간과 곡률 근사 검증]({{ '/imls/14-sr1-validation/' | relative_url }})
- [IMLS 연구 15 — 두 변형체 접촉과 희소 FEM 구현]({{ '/imls/15-two-body-sparse-fem/' | relative_url }})
- [IMLS 연구 16 — 독립 코드와 BVH 및 CCD 연결]({{ '/imls/16-src-bvh-ccd/' | relative_url }})

## 현재 구현과 검증 범위

| 분야 | 구현·관측 | 현재 한계 |
|---|---|---|
| 2D·3D FEM | 선형 탄성·압축 기준해·희소 조립 | 소변형 solid 중심, cloth 미구현 |
| 단일-gap 최적화 | 보폭·BFGS·모드 Oracle·Structured/Gap SR1 | 밀집 SR1이며 다중 접촉 SR1은 미구현 |
| 반복 benchmark | 14케이스×5방법×7회=490회 모두 수렴 | CPU 단일스레드·작은 모델, 모든 경우의 가속 아님 |
| 곡률 근사 | 정확한 gap Hessian은 오프라인에서만 비교 | 변형 패치의 새 방향에서 큰 오차도 관측 |
| 두 변형체 | 양방향 표면 적분·작용반작용·정적 압축 | 선형 탄성·무마찰, 대변형/동역학 검증 아님 |
| 자동 탐색 | BVH FV/EE 후보→적분점별 RIMLS 패치 | 복잡한 접힘·자기접촉·패치 교체 연속성 미검증 |
| CCD | 공식 IPC Toolkit의 삼각형 경로 검사 연결 | RIMLS 영점 표면 비관통 보장과 다름 |
| 새 개발 코드 | `src/` 아래 기능 폴더, `fem3d`와 독립 | 현재 CPU, 지역 접촉 AD는 밀집 |

Gap SR1은 어려운 고정 곡면에서 크게 빨라졌지만 변형 패치 E=10에서는 시간 차이가 거의 없었다. 수렴·미분 일치는 해당 이산 모델의 수치 검증이며 범용 물리 정확성이나 신규성을 확정하는 근거는 아니다.

## 실행 구조와 재현

연구 코드 루트 기준 상대 경로다. 블로그에는 논문 PDF·그림, 화면 캡처·영상, 로컬 사용자 경로와 원본 실험 결과를 복사하지 않는다.

- `fem2d/`: 기존 2D 기준 구현.
- `fem3d/core`, `problems`, `solvers`, `experiments`, `diagnostics`: 독립 3D 기준 구현과 연구 실험. 새 `src`를 부르는 wrapper가 아니다.
- `src/geometry`, `mechanics`, `collision`, `contact`, `optimization`, `simulation`: 새 독립 개발 코드. 중간 `rimls_contact/` 폴더는 없다.
- `ipc/`: 기존 IPC 비교·동역학 실험.
- `tests/`: 수치·연결·독립성 검증.

```bash
python -B -m unittest discover -s tests -v
python -B fem3d/run.py --list
python -B fem3d/run.py sr1_validate --repeats 7
python -B fem3d/run.py two_body --methods gn --steps 20 --show
python -B run.py two-body --contact auto --methods gn --steps 20 --show
```

`--show`는 화면 표시이며 그림·영상 저장이 아니다. 새 두 블록 CLI는 정적 하중 단계의 애니메이션이다. 기존 `ipc`의 힘 구동 backward Euler 실험과 구분한다. 현재 환경의 PyTorch·SciPy·ipctk 계산 경로는 CPU다.

## 검증 기록과 숫자의 출처

2026-10-30 게시 점검에서 현재 연구 코드 전체 unittest를 독립 재실행해 **116개 모두 통과**했다(57.939초). 기록 작업은 연구 코드를 수정하지 않았다. 테스트 개수는 메서드 수이며 임의 형상의 모든 조건을 검증했다는 뜻은 아니다.

14개 케이스 시간·곡률·프로파일 표는 보존된 `sr1_validation_results/`의 CSV·README에서 확인했다. 490회 benchmark 자체를 게시 작업에서 다시 실행한 것은 아니다. 앞선 보폭·모드 실험도 해당 결과 파일을 읽었다. 10월 28일 글의 BLAS 스레드 조건이 다른 단회 시간과 최신 중앙값을 섞지 않는다.

플랫폼 선택 probe의 임시 산출물은 남기지 않았고 기존 연구 결과와 frozen baseline은 보존했다. SOFA·libuipc는 설치·성능 검증 완료 대상으로 쓰지 않는다.

## 읽은 논문

- Y. Du, Y. Li, S. Coros, B. Thomaszewski. [Robust and Artefact-Free Deformable Contact with Smooth Surface Representations](https://doi.org/10.1111/cgf.15187). Computer Graphics Forum 43(8), 2024.
- K. Huang, F. M. Chitalu, H. Lin, T. Komura. [GIPC: Fast and Stable Gauss-Newton Optimization of IPC Barrier Energy](https://doi.org/10.1145/3643028). ACM Transactions on Graphics 43(2), Article 23, 2024.

논문의 줄·매듭·자기접촉 장면은 현재 코드의 검증 범위가 아니다. 논문의 가속률을 우리 방법의 성능 근거로 전용하지 않는다.

## 업데이트와 다음 연구

- 2026-10-28: 1–12편. 표면 표현·2D/3D FEM·접촉·동역학·곡률과 초기 시간 관측.
- 2026-10-29 작업, 2026-10-30 정리: 13–16편. 보폭과 모드 진단, SR1·정식 반복 측정, 두 변형체·희소화, 독립 코드·자동 탐색·CCD.

다음은 탐색·RIMLS 미분·조립·선형 풀이를 큰 모델에서 나눠 측정하고, CPU↔GPU 전송까지 포함한 실제 전체 시간 감소를 검증하는 것이다. 범용 FEM 툴 개발이 아니라 연구 접촉 방법을 비교할 수 있는 작은 in-house 코드가 목적이다. 반복 예약은 설정하지 않았다.
