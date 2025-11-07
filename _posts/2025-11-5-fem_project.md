---
layout: post
title:  "Paper Rebuild 시리즈 #1: 논문 상황 정의"
toc: true          # 우측 목차
toc_sticky: true   # 스크롤해도 고정
date: 2025-11-5
tags: [markdown, blog]
use_math: true
mathjax: true
---
# Paper Rebuild 시리즈 #1: 논문 상황 정의
## 동기 및 논문 선정 이유
FEM강의를 끝내고 직접 짠 코드가 진짜 되는지 논문에 적용과 검증을 해보고 싶었다. 
3차원 열해석 과정 그리고 결과 값이 명시되어 있는 논문을 찾았다. 이 논문은 FEM 문제는 아니지만 결과 값이 분명하게 명시되어 있어 선택하였다. 
## 논문 케이스
- $x \in [0,1] \quad y,z \in [0,20] $
- 균일 열유속 $q_0 = 1$ (패치)
- 등온(Dirichlet B.C.) : $T=0$ on $y=20, z=20$
- 나머지 면은 단열 (Neuman B.C.)
  
<img width="524" height="328" alt="Screenshot from 2025-11-05 19-23-04" src="https://github.com/user-attachments/assets/6e5a0bb3-6912-42b5-a6ef-6f4dc88fc281" />{: .centered }
<img width="509" height="179" alt="Screenshot from 2025-11-05 19-23-39" src="https://github.com/user-attachments/assets/bc5cb465-afa3-4207-8213-35d471022eca" />{: .centered }

## 목표 문제 정의 
- 도메인: $x \in [0,1], y,z \in [0,20] \(b=20)$  
- Neuman B.C. (가열 패치) : $ x=1 \quad \& \quad y,z \in [0,2] (a=2)$에서 유속 $q_0=+1$로 들어감
- Dirichlet B.C : $T=0 on y=20 or z=20$ 두 경우 중 하나라도 해당되면 0 // 그 외 바깥면은 단열  
- 비교 지표: 가열면 중심 ($ x=0, y=a/2 , z= a/2 $)의
$$
\tilde{T}_L =\frac{k}{q_0L}(-T), k=385, q_0=1,L=1
$$
<img width="995" height="139" alt="Screenshot from 2025-11-05 19-36-51" src="https://github.com/user-attachments/assets/b0c6380b-e1d4-443b-b73a-0e3a451ae697" />{: .centered }

### 참고문헌  

Beck, J. V., Transient three-dimensional heat conduction problems with partial heating, International Journal of Heat and Mass Transfer, 54 (2011), 2479–2489. DOI: 10.1016/j.ijheatmasstransfer.2011.02.014.


