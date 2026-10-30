---
layout: single
title: "IMLS 연구 16 — 독립 코드와 BVH 및 CCD 연결"
date: 2026-10-29 23:00:00 +0900
last_modified_at: 2026-10-30
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact, optimization]
permalink: /imls/16-src-bvh-ccd/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 2026-10-29 연구 및 학습, 2026-10-30 정리

## 새 코드와 기존 검증 코드를 분리했다

새 개발 코드를 `src/` 바로 아래에 역할별로 배치했다. 중간에 두었던 `rimls_contact/` 하위 폴더는 제거했다. 기존 `fem3d`를 새 코드의 wrapper로 남기는 방식도 철회하고 자체 구현을 복원했다. `src`는 `fem3d`를 import하지 않고 `fem3d` 역시 `src`를 호출하지 않는 독립 기준 코드다. 같은 기능을 유지하는 비용은 있지만 이후 수치 비교의 기준이 된다.

```text
src/
├── geometry/       표면 샘플·법선·RIMLS·패치 선택
├── mechanics/      소변형 solid FEM
├── collision/      IPC Toolkit BVH·CCD 연결
├── contact/        barrier와 접촉 미분
├── optimization/   방향·안정화·선 탐색·단일-gap SR1
└── simulation/     문제 설정·하중 단계·표시
```

`mechanics/cloth.py`는 미구현 인터페이스다. 관성 포텐셜 보조 함수도 있지만 새 두 물체 CLI는 정적 하중 단계이지 완전한 동역학 드라이버가 아니다.

## 플랫폼은 무엇을 선택했나

현재 선택은 Python·PyTorch·SciPy 연구 코어와 공식 IPC Toolkit 1.3.0 충돌 백엔드다. IPC Toolkit은 FEM 엔진 전체가 아니라 충돌 후보·CCD·IPC 관련 계산을 제공한다. 우리 코드에서 FEM과 RIMLS 접촉 에너지는 자체 구현한다.

SOFA는 프레임워크와 사용자 force field 확장 경로, libuipc는 GPU 시뮬레이션과 FEM·천·줄 예제를 제공하는 후보로 검토했다. 그러나 이번에 설치·실행 성능을 검증한 것은 아니다. RIMLS Hessian과 반복 루프를 교체하는 경로를 확인하지 않고 GPU 예제 실행만으로 연구 방법이 이식됐다고 할 수는 없다.

공식 참고: [IPC Toolkit](https://github.com/ipc-sim/ipc-toolkit), [ipctk 1.3.0](https://pypi.org/project/ipctk/1.3.0/), [SOFA Python force field 예](https://github.com/sofa-framework/SofaPython3/blob/master/examples/emptyForceField.py), [libuipc](https://github.com/spiriMirror/libuipc).

기존 Torch 1.9 CPU와 ipctk의 native 라이브러리 충돌을 피하려고 ipctk는 Torch를 import하지 않는 지속 실행 child process에 둔다. 부모가 좌표와 요청을 보내고 숫자 결과를 받는다. `collision/_worker.py`가 공식 라이브러리를 호출하고 `collision/ipc_toolkit.py`가 연결을 담당한다. `import math`는 Python 기본 수학 모듈이고 IPC import는 `ipctk`다.

## 탐색과 접촉 응답은 다른 계산이다

처음 두 블록은 접촉면을 직접 지정했다. 양쪽이 변형돼도 지정한 면 안에서 힘은 계산할 수 있지만 어디가 닿는지 찾는 일반 문제를 해결한 것은 아니었다. 이를 전체 경계 삼각형 후보 탐색과 연결했다.

```text
변형된 전체 경계면
 → BVH의 face–vertex / edge–edge 후보
 → 원래 FEM 절점·삼각형 번호 복원
 → incidence로 삼각형 쌍 구성·중복 제거
 → 적분점별 상대 표면 패치 선택
 → RIMLS gap·barrier·접촉 미분
```

BVH는 계층적 bounding box로 멀리 있는 부분을 제외한다. 후보 번호를 compact 표면 메시에서 전체 FEM 번호로 되돌려야 하며 모서리에 연결된 면들로 확장할 때 중복을 제거한다. 후보 수를 접촉 에너지 개수로 쓰면 같은 접촉을 여러 번 더할 수 있다.

기준 형상에서 연결되고 법선 차이가 30도 이내인 면들을 smooth patch로 묶는다. 블록 모서리에서 서로 다른 면 법선을 섞지 않기 위한 선택이다. 현재 위치의 대향 법선·상대 위치로 필터한 뒤 가우스 적분점마다 상대 물체의 가까운 패치 한 개를 고른다. 복잡한 접힘에서 여러 패치가 필요한 상황은 아직 검증하지 않았다.

## 0.3과 0.2와 R은 어떻게 다른가

| 값 | 현재 기본값 | 역할 |
|---|---:|---|
| 탐색 거리 $\hat d+\mathrm{skin}$ | 0.3 | 가까운 패치를 후보로 고려 |
| Barrier 활성 gap $\hat d$ | 0.2 | RIMLS gap이 이보다 작으면 접촉 에너지 발생 |
| RIMLS 커널 반경 $R$ | 1.5 | 선택된 패치에서 gap에 기여하는 샘플 범위 |

IPC bounding box에는 검색 거리의 절반인 0.15를 inflation으로 전달하고 이후 실제 점–삼각형 거리 0.3 조건으로 다시 거른다. 후보라고 바로 힘이 생기지 않는다. 탐색용 삼각형 거리와 RIMLS gap도 서로 다른 양이다.

`geometry/contact_patches.py`의 `AutomaticRIMLSContact.refresh()`에서 선택을 갱신한다. 자동미분 밖에서 Newton 선형화 직전에 선택하고 방향 계산과 선 탐색 동안 고정한다. skin/4의 보수적 절점 이동 제한과 패치 교체 에너지 jump 기록을 추가했지만 이 규칙만으로 전역적인 에너지 매끄러움을 보장하지 않는다.

## CCD는 이동 경로를 검사한다

탐색은 어떤 표면을 쓸지 고르고 CCD는 제안한 $q\to q+p$의 선형 절점 경로에서 삼각형이 통과하는지 검사한다. 공식 `compute_collision_free_stepsize()`가 반환한 보폭 이하에서 Armijo 선 탐색을 시작한다.

$$
q_{\mathrm{new}}=q+\alpha p,\qquad \alpha\leq\alpha_{\mathrm{CCD}}.
$$

도착점만 확인하는 것과 다르다. crossing probe에서 충돌 없는 한계는 약 0.4999466으로 해석적 교차 0.5보다 작았다. 하지만 삼각형 CCD가 RIMLS 영점 표면의 비관통을 증명하는 것은 아니다. 두 표현 사이 차이는 연구 과제다.

## 무엇으로 검증했나

평면 압축에서 자동·수동 접촉 선택의 Full·GN·filtered GN 최종 좌표 차이는 최대 약 $3.1\times10^{-10}$였다. 이 예제의 패치 교체 에너지 jump는 0이었다. 옆으로 분리했다 되돌리는 선택 해제·재생성, 비연속 표면 번호 매핑, FV/EE 중복 제거, 고정 패치의 유한차분 미분, 방법별 패치 상태 분리도 테스트했다.

한 압축 상태에서 BVH 후보 488개가 288개 face pair와 최종 48개 접촉 샘플로 정리됐다. 개수는 해당 메시 상태에 한정된다. 2026-10-30 게시 점검에서 전체 unittest를 독립 재실행해 116개 모두 통과했다. 연구 코드는 수정하지 않았다. 복잡한 곡면 전환·자기접촉·대변형까지 검증한 숫자는 아니다.

```bash
# 새 독립 코드
python -B run.py two-body --contact auto --methods gn --steps 20 --show
# 수동 선택 대조군
python -B run.py two-body --contact manual --methods gn --steps 20 --show
# 기존 독립 기준 코드
python -B fem3d/run.py two_body --methods gn --steps 20 --show
```

선택 실험의 임시 산출물은 남기지 않았다. 기존 연구 결과·frozen baseline은 유지했고 이 기록에는 논문 그림이나 화면 캡처·자동생성 영상을 게시하지 않는다.

## CPU라 느리다면 다음은 GPU인가

현재 NumPy·SciPy 조립과 선형 풀이, Torch 텐서·자동미분, 선택한 IPC BVH·CCD 경로는 모두 CPU다. SciPy는 물리를 정하는 것이 아니라 희소행렬과 수치 선형대수를 제공한다. GPU 장치가 있어도 CPU Torch에서 자동으로 GPU를 사용하지 않는다.

IPC Toolkit 전체에 GPU 가능성이 없다고 단정할 수는 없지만 현재 연결 경로는 CPU다. RIMLS 텐서만 CUDA로 옮기고 매번 NumPy로 되돌리면 전송이 병목이 될 수 있다. 지금은 GPU 속도 개선 결과가 없다.

다음 목표는 기능을 더 쌓는 것보다 탐색·gap 미분·조립·풀이 시간을 분해하고 GPU에 계산과 데이터를 유지할 경로를 검증하는 것이다. 성공 기준은 GPU에서 실행됐다는 사실이 아니라 같은 정확도의 전체 해석시간 감소다. 기존 `fem3d`는 CPU 수치 기준으로 남긴다.

