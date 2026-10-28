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

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 검증 기준일: 2026-10-28

## 연구 목적

IMLS/RIMLS의 매끄러운 표면 표현으로 불연속 표면 이산화의 인공 접촉력을 줄일 수 있을까? GIPC에서 영감을 얻은 Hessian 근사로 변형체 접촉 계산을 가속할 수 있을까? 이 질문을 위해 2D 접촉의 작은 문제부터 3D 선형 FEM 기초까지 공부하고 구현했다.

## 학습 순서와 글

- [IMLS 연구 1 — 정오각형 샘플에서 IMLS·RIMLS 표면까지]({{ '/imls/01-geometry-imls/' | relative_url }})
- [IMLS 연구 2 — 접촉 barrier, Full·GN Hessian과 자동미분]({{ '/imls/02-barrier-hessian/' | relative_url }})
- [IMLS 연구 3 — 2D FEM, 자유 DOF와 접촉 Newton 해법]({{ '/imls/03-fem2d-solver/' | relative_url }})
- [IMLS 연구 4 — 평면 슬라이딩과 공식 IPC 비교 읽기]({{ '/imls/04-sliding-ipc/' | relative_url }})
- [IMLS 연구 5 — 독립 기준해로 2D FEM과 접촉 검증하기]({{ '/imls/05-reference-checks/' | relative_url }})
- [IMLS 연구 6 — 3D 사면체 메시와 선형 FEM 압축 검증]({{ '/imls/06-fem3d/' | relative_url }})


문답은 실제 학습 주제를 재구성한 설명으로 대화의 직접 인용은 아니다. 구현·검증한 사실과 아직 확인하지 못한 주장을 구분한다.

## 현재까지의 진행

| 단계 | 구현·확인한 내용 | 범위 |
|---|---|---|
| 표면 | 2D IMLS, 한 번의 robust 재가중, 샘플과 법선 | 일반 형상의 정확한 최단거리는 아님 |
| Stage 3 | 고정 표면+스프링 질의점, Full/GN/filtered GN | scalar gap 근사 비교 |
| Stage 4 | 2D 소변형 평면응력 FEM+질의점, Torch 자동미분 | 기본 문제 합계 8 DOF |
| 슬라이딩 | 고정 평면 스캔, 변형 블록, 공식 ipctk 비교 | 마찰·self-contact 제외, RIMLS CCD 없음 |
| 독립 검증 | 2D 압축, 평면 접촉, FEM+접촉 결합 | 평면·균일 변형 |
| 3D 기초 | 사면체 메시·선형 탄성·압축 | 3D contact 구현 전 |

현재 GN/filtered GN은 GIPC 전체 방법의 구현이 아니다. 논문 노벨티와 속도 우월성은 입증되지 않았다. 논문의 3D 예제를 재현 완료한 것도 아니다.

## 실행과 재현의 입구

경로는 연구 코드 루트 기준이다. 현재 2D 코드는 `fem2d/`, 3D는 `fem3d/`, 검증은 `tests/`, 공식 IPC 비교는 `ipc/`로 나누었다. 블로그에는 연구 소스 본체·논문 PDF·그림·임시 백업을 복사하지 않았다. 아래 명령은 연구 코드를 가진 환경에서 실행한다.

```bash
python -m pip install -r requirements-autodiff.txt
python -m pip install matplotlib
python -m pip install -r ipc/requirements-ipc.txt
python tests/test_fem_compression.py
python tests/test_plane_contact_reference.py
python tests/test_fem_contact_compression.py
python tests/test_fem3d_compression.py
python fem2d/sliding_plane.py
python fem2d/block_sliding.py --compare-ipc --show
```

검증 환경은 NumPy·Torch·matplotlib, 공식 ipctk 1.3.0을 사용한다. 임의의 Python 버전이나 모든 OS에서의 재현을 보장한 것은 아니다. 위 4개 검증 스크립트는 이번에 모두 재실행했고 총 9개 테스트 메서드가 통과했다. 각 메서드는 여러 조건을 포함한다. 평면 스캔과 기본 5위치 IPC 비교도 재실행했다. 41위치 바닥 세분화 결과는 기존 CSV의 재집계로 구분했다.

## 읽은 논문

- Y. Du, Y. Li, S. Coros, B. Thomaszewski. [Robust and Artefact-Free Deformable Contact with Smooth Surface Representations](https://doi.org/10.1111/cgf.15187). Computer Graphics Forum 43(8), 2024. 접선력 문제, IMLS 접촉·self-contact·CCD를 다룬다. 현재 구현은 후자의 self-contact와 CCD를 포함하지 않는다.
- K. Huang, F. M. Chitalu, H. Lin, T. Komura. [GIPC: Fast and Stable Gauss-Newton Optimization of IPC Barrier Energy](https://doi.org/10.1145/3643028). ACM Transactions on Graphics 43(2), Article 23, 2024. 접촉 Hessian 근사와 최적화 연구다. 논문의 가속률을 현재 코드의 속도 근거로 사용하지 않는다.

## 업데이트와 다음 단계

2026-10-28: 최초 기록. 2D 표면부터 3D 압축 검증까지 6편으로 정리하고 현재 코드·수치 결과와 대조했다.

다음 학습은 3D 메시·FEM의 이해를 다진 뒤의 3D 접촉이다. 이 인덱스를 시리즈 입구로 유지한다. 새 진행사항이 오면 해당 글을 보완하거나 새 회차를 추가하고, 검증일·조건·결과·미해결 사항을 함께 갱신한다. 자동 게시 예약이나 반복 실행은 설정하지 않았다.
