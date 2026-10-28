---
layout: single
title: "IMLS 연구 12 — 반복 횟수보다 총 시간, 코드 정리와 재검증"
date: 2026-10-28 16:50:00 +0900
last_modified_at: 2026-10-28
categories: [IMLS]
tags: [IMLS, RIMLS, FEM, contact]
permalink: /imls/12-runtime-and-refactor/
toc: true
toc_sticky: true
use_math: true
mathjax: true
locale: ko-KR
---

[시리즈 인덱스]({{ '/imls/index/' | relative_url }}) · 2026-10-28 추가 학습 기록

## 오늘의 교정: 반복이 많으면 느리다는 결론은 성립하지 않는다

[구면 곡률]({{ '/imls/10-sphere-curvature/' | relative_url }})과 [여러 곡면 비교]({{ '/imls/11-curved-contact-cases/' | relative_url }})에서 GN의 업데이트 횟수가 늘어나는 이유를 분석했다. 그러나 사용자가 지적한 대로 **GN은 한 반복이 싸서 총 시간으로 더 빠를 수 있다.** 실제 저장된 CPU 측정에서도 수렴한 예제 다수는 그랬다.

## 같은 허용오차에 수렴한 실행시간 비교

아래는 `curved_contact_results/summary.csv`와 `deformable_curved_results/summary.csv`의 단회 solve-loop 시간이다. warm-up·실행 순서·반복 측정이 통제된 정식 benchmark가 아니다. 비율은 해당 실행에서 Full 시간이 GN 시간의 몇 배였는지 보여 준다.

| 저장된 수렴 예제 | Full 업데이트 / 시간 | GN 업데이트 / 시간 | 관측한 Full 시간 ÷ GN 시간 |
|---|---:|---:|---:|
| 구면 접선 목표 | 7 / 59.5 ms | 20 / 35.1 ms | 1.70 |
| 원통 둘레 방향 | 7 / 90.7 ms | 19 / 46.8 ms | 1.94 |
| 타원체 x 방향 | 5 / 48.2 ms | 14 / 25.7 ms | 1.88 |
| 오목 구면 접선 목표 | 16 / 122.0 ms | 24 / 47.5 ms | 2.57 |
| 안장면 x 방향 | 6 / 63.5 ms | 13 / 27.4 ms | 2.32 |
| 변형 패치 E=10 | 6 / 0.905 s | 13 / 0.156 s | 5.80 |
| 변형 패치 E=30 | 6 / 0.888 s | 10 / 0.116 s | 7.63 |


예를 들어 57 DOF 변형 패치 E=30은 Full 6회보다 GN 10회가 반복은 많지만 시간은 약 0.888초 대 0.116초였다. 반복 수만으로 계산이 악화됐다고 쓰면 이 결과를 반대로 읽게 된다.

반대로 `sphere_contact_results`의 깊은 목표는 Full 14회·0.119초, GN 91회·0.149초로 Full이 조금 빨랐다. 중심 근처의 GN 300회 예산 종료와 Full 수렴 시간을 비교해 같은 정확도의 속도비라고 해석할 수는 없다. 충분한 예산의 수렴 결과, 잔차뿐 아니라 좌표 정확도도 함께 맞춰야 한다.

## 총 시간은 어떤 항으로 구성되나?

대략적인 비용은 반복마다 gap·gradient·Hessian 평가, 선형계 풀이, 선 탐색 중 추가 평가의 합이다. Full은 gap Hessian을 구하지만 GN은 이를 생략한다. 작은 3 DOF 문제와 57 DOF 변형 표면의 비용 구성이 다르며, 곡률 보정을 추가하면 반복 수는 줄어도 비싼 미분 때문에 전체는 느려질 수 있다.

다음 성능 비교는 같은 초기 상태·모델·허용오차에서 warm-up 후 실행 순서를 바꾸어 여러 번 측정하고 중앙값과 분산을 보고해야 한다. 현재 수치는 연구 방향을 정하는 관찰이지 일반적인 가속률 약속이 아니다. IPC native worker+통신 시간과 Torch 모델 시간을 직접 비교하는 문제도 별도로 남아 있다.

## 모델을 바꾸지 않는 코드 정리는 무엇을 했나?

- `geometry.py`: 평면·구면·원통·타원체·안장면 샘플과 법선, 곡면 패치 좌표 변환을 모았다. 체적 연결 관계는 `mesh.py`에 남긴다. `itertools.permutations`는 Python 표준 라이브러리의 축 순열 생성이며 Torch 연산이 아니다.
- `optimization.py`: Newton 방향과 고유값 shift 계산을 공통화했다. 호출별 고유값 바닥값을 유지하고 shift가 0이면 불필요한 단위행렬 할당을 피한다. 새 결합 예제는 공통 정적 반복기를 사용한다.
- `contact.py`: 여러 질의점의 RIMLS를 Tensor batch로 평가한다. scalar와 같은 공간·잔차·법선 가중식이며 값·Jacobian·Hessian 일치를 검사했다. 비교 코드가 공급한 바닥 샘플도 받아 중복 생성을 줄였다.
- `fem.py`: 질량 계산에 필요한 사면체 체적만 벡터화해 얻는다. 체적만 필요할 때 형상함수 역행렬과 B까지 구하지 않는다.

학습 중 24질의·63평면 샘플의 Jacobian 계산에서 5회 측정한 지역 비용은 scalar 약 0.2096초, batch 약 0.02054초였다는 실행 기록이 있다. 약 10.2배라는 숫자는 이 미세 구간의 관찰이며 전체 solver 가속률이 아니다. 해당 임시 benchmark JSON은 정리 과정에서 삭제되어 현재 파일을 다시 집계한 수치로 제시하지 않는다.

## 파일 이관과 __main__의 의미

3D 정식 구현은 이제 `mesh.py`, `fem.py`, `surface.py`, `contact.py`, `compression.py`, `block_sliding.py`, `sphere_contact.py`다. 예전의 `test_mesh.py`·`test_fem.py`를 import하는 설명은 더 이상 현재 구조가 아니다. IPC도 `ipc_comparison.py`, `ipc3d_comparison.py`, `ipc3d_dynamics.py`, `ipc_worker.py`가 정식 이름이다.

`ipc/test_time_convergence.py`와 `fem3d/test_sphere_curvature_diagnostic.py`는 아직 공부용이다. 빈 `ipc/time_convergence.py`를 실행하면 시간 수렴 실험이 수행되는 것은 아니다. `if __name__ == '__main__':` 아래 예제만 직접 실행 시 호출되므로, 함수 import 자체를 결과 검증으로 오해하지 않는다.

## 정리 후 전체 검증과 발견한 누락

게시 점검 중 불필요 import를 정리한 뒤 `ipc3d_comparison.py`에서 필요한 `torch`까지 빠진 것을 발견했다. 이 중간 상태는 51개 중 48개 통과·3개 `NameError`였다. 연구 채팅에서 import를 복구한 뒤 전체를 재실행해 **51개 전부 통과**했다. 기존 연구 결과와 현재 실행 상태를 섞지 않기 위해 최종 검증 결과를 기준으로 기록한다.

```bash
python -B -m unittest discover -s tests -v
```

검증 범위는 기존 2D·3D FEM, 3D 접촉 압축과 슬라이딩, IPC gradient·PSD·작용반작용·CCD·내부 절점 매핑, 질량·힘 구동·영 수평하중 대조군, 구면 미분, 여러 곡면·변형 표면 결합을 포함한다. 51은 테스트 메서드 수이고 각 메서드에는 여러 조건이 들어 있다. 긴 72조합 sweep이나 50·100·200스텝 전체를 이 게시 작업에서 새로 실행한 것은 아니다.

## 다음 연구의 판단 기준

현재 증거는 GN이 많은 조건에서 이미 비용 효율적이라는 것이다. 따라서 연구 질문은 모든 곡률 항을 항상 복원하는 것이 아니라, GN의 **총 시간 또는 안정성**이 실제로 나빠지는 조건에서만 유용한 곡률을 싸게 반영할 수 있는가로 좁혀진다. 반복 수·총 시간·허용오차·수렴 실패를 함께 비교해야 한다. 이는 후속 가설이며 새 방법의 성능이나 논문 노벨티는 아직 입증하지 않았다. scalar gap의 filtered GN을 실제 GIPC 전체 방법으로 부르지도 않는다.
