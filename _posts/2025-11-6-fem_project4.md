---
layout: post
title:  "Paper Rebuild 시리즈 #4: haning node 원인: 국부정련"
toc: true          # 우측 목차
toc_sticky: true   # 스크롤해도 고정
date: 2025-11-6
tags: [markdown, blog]
use_math: true
mathjax: true
---
## 국부 정련 끈 그래프끼리 비교
<img width="801" height="714" alt="nohanging_(20,40,40)vs(24,48,48)" src="https://github.com/user-attachments/assets/256fd86e-6f0e-496a-ba8c-9fc7b8ec4eb1" />
1. haning node는 국부 정련을 했을 때, 큰 셀의 모서리/변 중간에만 걸려 있는(맞은편 셀 꼭짓점과 안 맞는) ‘매달린’ 절점이다. 원래는25,50,50으로 매쉬를 나눴는데
2. 그래프 해석: 그래프가 끊기는 현상은 발생하지 않았지만 값이 이론값과 많이 달라졌고 x=1로 가면서의 값이 두 경우 모두 비슷한 값이여야 하지만 아예 다른 경향을 보인다. $\rightarrow$ 뭔가 잘못됨.

## ubuntu를 멀티 부팅으로 다운
<img width="651" height="782" alt="coarsettilde" src="https://github.com/user-attachments/assets/10461b1f-0af6-4fa8-a984-2e1096423b9a" />  

1. 멀티 부팅으로 ubuntu 다운 받으니까 계산 속도가 빨라졌다. 그래서 30,60,60을 넣었는데 매쉬가 딱 맞게 들어가서 그런지 haning node가 발생하지 않았다.

<img width="801" height="714" alt="bc_edit_(20,40,40)vs(30,60,60)" src="https://github.com/user-attachments/assets/69063397-3443-4c09-8ef7-a024cc393927" />  

2. 30,60,60이 괜찮았던 이유:
 - 정련 경계의 위치가 격자 분할 수에 따라 달라졌다.
   - Ny=60이면 hy=20/60=1/3 $\rightarrow$ y<3이 딱 9칸(정수)라 경계가 깔끔
