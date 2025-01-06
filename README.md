# How to Draw Pole figure

## 1. Basic Preparation
Install [`texture3`](https://github.com/youngung/texture3) from github using clone or download
```
$ pip install .
```
Execute [```setup.py```](setup.py) in folder `texture3` 

## 2. Run command
```
$ ipython --pylab 
```
Enter python in cmd
```
$ from TX import upf
```
Load ```upf``` from module ```TX```
```
$ mypf = upf.polefigure(fnsx='fnsx file name', filename='TEX_PH1.OUT')
```


Save data in `mypf`,
Use `fnsx`file used in evpsc simulation and `TEX_PH1.OUT` created after the simulation as Arguments
```
$ mypf.pf_new(poles=[[1,0,0],[1,1,0]],mn=0.5,mx=3.5,ismooth=10)
```
```pf_new```  Draw Polefigure 

| Arguments                |    Description                            |
|:-------------------------|:------------------------------------------|
| poles                    | Miller index of Polefigure                |
| mn                       | Minimun level of contour                  |
| mx                       | Maximun level of contour                  |
| ismooth                  | Contours become smooth                    |















# Pole figure 그리기

## 1. 기본 준비
 git hub로부터 texture3를 clone 혹은 download를 통해 설치
```
$ pip install .
```
texture3 폴더내의 ```setup.py```실행
## 2. 명령어 실행
```
$ ipython --pylab 
```
cmd 창에서 파이썬실행
```
$ from TX import upf
```
module ```TX```에서 ```upf```불러오기
```
$ mypf = upf.polefigure(fnsx='fnsx file name', filename='TEX_PH1.OUT')
```
```mypf```에 데이터 저장, evpsc 시뮬레이션에 사용한 ```fnsx```파일, 시뮬레이션 이후  생성된 ```TEX_PH1.OUT```를 Arguments로 사용
```
$ mypf.pf_new(poles=[[1,0,0],[1,1,0]],mn=0.5,mx=3.5,ismooth=10)
```
polefigure 생성
```pf_new``` 
| Arguments                |    Description                            |
|:-------------------------|:------------------------------------------|
|poles                     | Polefigure의 Miller index를 지정           |
|mn                        | 등고선의 최솟값을 지정                      |
|mx                        | 등고선의 최댓값을 지정                      |
|ismooth                   | 등고선이 각지지 않고 부드럽게 나타냄         |
