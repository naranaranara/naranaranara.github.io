---
layout: single
title: "IMLS·RIMLS 접촉 연구 — 학습 기록 인덱스"
date: 2026-10-28 12:00:00 +0900
last_modified_at: 2026-10-28
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/index/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[처음 공부한 기하·표면 표현부터]({{ '/imls/01-geometry-imls/' | relative_url }}) · [오늘 추가한 3D 접촉부터]({{ '/imls/07-fem3d-contact/' | relative_url }})

## 연구 목적과 오늘의 핵심

매끄러운 IMLS/RIMLS 표면으로 표면 이산화의 인공 접촉력을 줄이고, GIPC에서 영감을 얻은 Hessian 근사로 계산을 효율화할 수 있는지 공부한다. 현재는 3D FEM·접촉의 정적/동적 검증, 구면 곡률 진단, 여러 곡면과 변형 표면 실험까지 진행했다.

오늘의 핵심은 세 가지다. 지정변위와 힘 구동 동역학을 구분했고, 평면 RIMLS의 시간 1차 수렴을 확인했으며, GN의 반복이 많아도 전체 실행시간은 더 짧을 수 있음을 수치로 확인했다. 목표는 논문 그림과 비슷한 장면 만들기가 아니라 수식·미분·평형·적분을 검증하는 것이다.

## 학습 순서와 글

### 기초: 1–6편

- [IMLS 연구 1 — 정오각형 샘플에서 IMLS·RIMLS 표면까지]({{ '/imls/01-geometry-imls/' | relative_url }})
- [IMLS 연구 2 — 접촉 barrier, Full·GN Hessian과 자동미분]({{ '/imls/02-barrier-hessian/' | relative_url }})
- [IMLS 연구 3 — 2D FEM, 자유 DOF와 접촉 Newton 해법]({{ '/imls/03-fem2d-solver/' | relative_url }})
- [IMLS 연구 4 — 평면 슬라이딩과 공식 IPC 비교 읽기]({{ '/imls/04-sliding-ipc/' | relative_url }})
- [IMLS 연구 5 — 독립 기준해로 2D FEM과 접촉 검증하기]({{ '/imls/05-reference-checks/' | relative_url }})
- [IMLS 연구 6 — 3D 사면체 메시와 선형 FEM 압축 검증]({{ '/imls/06-fem3d/' | relative_url }})

### 2026-10-28 추가: 7–12편

- [IMLS 연구 7 — 3D 표면 적분과 FEM·접촉 결합 검증]({{ '/imls/07-fem3d-contact/' | relative_url }})
- [IMLS 연구 8 — 지정변위에서 힘 구동 동역학으로, 힘의 부호 검증]({{ '/imls/08-force-driven-dynamics/' | relative_url }})
- [IMLS 연구 9 — 시간 간격을 줄여 무엇을 검증했나]({{ '/imls/09-time-convergence/' | relative_url }})
- [IMLS 연구 10 — 구면 곡률과 GN의 느린 접선 수렴]({{ '/imls/10-sphere-curvature/' | relative_url }})
- [IMLS 연구 11 — 여러 곡면과 변형 표면의 접촉 실험]({{ '/imls/11-curved-contact-cases/' | relative_url }})
- [IMLS 연구 12 — 반복 횟수보다 총 시간, 코드 정리와 재검증]({{ '/imls/12-runtime-and-refactor/' | relative_url }})


각 글의 문답은 실제 학습 주제를 재구성한 설명이며 대화의 직접 인용은 아니다. 초기 단계의 글은 남기되 현재 경로와 진행 상태는 갱신했다.

## 현재 구현과 검증 범위

| 단계 | 구현·확인 | 한계 |
|---|---|---|
| 2D 기초 | IMLS/RIMLS, 스프링 질의점, 평면응력 FEM·경계 적분 | 소변형·평면 기준해 중심 |
| 3D 정적 | 사면체 FEM+기준 면적 접촉, 압축·지정변위 슬라이딩 | RIMLS는 끝 상태 검사, CCD 없음 |
| 공식 IPC | 3D 표면 매핑·gradient·PSD·작용반작용·CCD | ipctk 기반 자체 simulator, 원 논문 전체 구현 아님 |
| 3D 동역학 | lumped mass, 힘 구동 backward Euler, 충격량·영 외력 대조군 | 소변형 선형 탄성, 논문 장면 파라미터 불일치 |
| 시간 검증 | RIMLS 평면 문제 위치 오차의 1차 감소 | IPC의 엄밀한 시간 수렴은 세 결과만으로 미확정 |
| 고정 곡면 | 72조합 중 64 수렴, 8은 300회 예산 종료 | 단일 질의점+스프링, 모델별 거리장 차이 |
| 변형 곡면 | 57 DOF 패치+질의점, 위치·법선까지 미분 | 두 변형체의 다중 접촉 검증은 아님 |
| 비용 평가 | 반복 수와 단회 총 시간을 별도 비교 | 정식 반복 benchmark·새 방법의 가속은 미입증 |

현재 GN/filtered GN은 scalar gap 기반 근사이며 GIPC 전체 방법이 아니다. 일반적인 노벨티·속도 우월성·논문 예제 완전 재현을 주장하지 않는다. 많은 수렴 사례에서는 GN이 이미 Full보다 시간상 유리했다.

## 실행 구조와 재현

연구 코드 루트 기준 상대 경로다. 블로그는 연구 소스 본체·논문 PDF·이미지·임시 검증 산출물을 복사하지 않는다.

- `fem2d/`: 기존 2D 구현.
- `fem3d/`: 정식 `mesh.py`, `fem.py`, `surface.py`, `contact.py`, `compression.py`, `block_sliding.py`, `sphere_contact.py`와 곡면·최적화 모듈.
- `ipc/`: 정식 `ipc_comparison.py`, `ipc3d_comparison.py`, `ipc3d_dynamics.py`, `ipc_worker.py`.
- `tests/`: 자동 검증. `ipc/test_time_convergence.py`와 `fem3d/test_sphere_curvature_diagnostic.py`는 아직 공부용 진입점.

```bash
python -m pip install -r requirements-autodiff.txt
python -m pip install matplotlib
python -m pip install -r ipc/requirements-ipc.txt
python -B -m unittest discover -s tests -v
python -B fem3d/compression.py
python -B fem3d/block_sliding.py
python -B ipc/ipc3d_comparison.py
python -B ipc/ipc3d_dynamics.py --steps 50 --dt 0.1 --force-x 0.002
python -B fem3d/sphere_contact.py --max-iter 300
```

확인 환경은 NumPy·Torch 1.9 CPU·matplotlib, 공식 ipctk 1.3.0이다. 임의의 Python/OS 전체를 검증한 것은 아니다. `-B`는 로컬 bytecode 캐시를 쓰지 않도록 한 옵션이다. 선택적으로 `--show`를 지원하는 스크립트는 화면 표시만 하며 그림 파일을 저장하지 않는다. 숫자 출력 경로는 각 글의 `--outdir` 예를 따른다.

## 이번에 새로 확인한 것과 기존 결과의 구분

게시 작업에서 수정된 현재 코드의 전체 자동검증을 재실행해 **51개 모두 통과**했다. 3D 압축 최대 좌표 오차 약 $8.07\times10^{-12}$를 재확인했다. 저장된 동적 배치를 이용한 공식 IPC 힘 부호 유한차분과 누적 접촉력 변위 계산도 다시 확인했다. 이 기록 작업에서 연구 코드는 수정하지 않았다.

50스텝 동역학, 3개 시간 간격, 72조합 곡면 sweep, 구면 확장 예산·곡률 진단, 변형 패치 결과는 보존된 CSV·JSON을 직접 읽어 재집계했다. 긴 실험 전체를 이번에 다시 수행했다고 표현하지 않는다. 테스트 51은 메서드 수이며 각 메서드는 여러 조건을 포함한다. 중간 import 누락은 연구 채팅에서 복구 후 전체 재검증을 마쳤다.

## 읽은 논문

- Y. Du, Y. Li, S. Coros, B. Thomaszewski. [Robust and Artefact-Free Deformable Contact with Smooth Surface Representations](https://doi.org/10.1111/cgf.15187). Computer Graphics Forum 43(8), 2024. 매끄러운 표면의 접촉, 자기접촉·CCD 연구. 현재 RIMLS의 적분·재료·충돌 검사와 차이를 구분한다.
- K. Huang, F. M. Chitalu, H. Lin, T. Komura. [GIPC: Fast and Stable Gauss-Newton Optimization of IPC Barrier Energy](https://doi.org/10.1145/3643028). ACM Transactions on Graphics 43(2), Article 23, 2024. 해석적 접촉 Hessian 근사와 최적화 연구. 논문의 가속률을 현재 코드의 성능 근거로 전용하지 않는다.

## 업데이트와 다음 연구

- 2026-10-28 최초: 1–6편과 인덱스 게시. 기하·2D 접촉부터 3D FEM 압축까지 정리.
- 2026-10-28 추가: 7–12편 게시. 3D 접촉·동역학·시간 수렴·곡률·여러 곡면·실행시간·코드 이관을 정리. 2·4·6편의 연결과 현재 경로·진행 상태 갱신.

다음 연구는 GN의 총 시간이나 안정성이 나빠지는 조건을 먼저 찾고, 그때 유용한 곡률을 저렴하게 반영할 수 있는지 검토하는 것이다. 반복 수만 줄이는 것으로 가속을 주장하지 않는다. 새 진행사항에 맞춰 해당 글과 검증일·조건·한계를 이어 갱신한다. 반복 예약은 설정하지 않았다.
