# How to Draw Pole figure

## 1. Basic Preparation
Install [`texture3`](https://github.com/youngung/texture3) from github using clone or download
```shell
$ pip install .
```
Execute [```setup.py```](setup.py) in folder `texture3`

## 2. Run command
```shell
$ ipython --pylab
```
Enter python in cmd
```python
from TX import upf
```
Load ```upf``` from module ```TX```
```python
mypf = upf.polefigure(fnsx='fnsx file name', filename='TEX_PH1.OUT')
```


Save data in `mypf`,
Use [`fnsx`](README/F_voce.sx) file used in evpsc simulation and [`TEX_PH1.OUT`](README/TEX_PH1_5000.OUT) created after the simulation as Arguments
```
$ mypf.pf_new(poles=[[1,0,0],[1,1,0]],mn=0.5,mx=3.5,ismooth=10)
```
```pf_new``` do Draw Polefigure

| Arguments      |    Description                                                    | default value   |
|:---------------|:------------------------------------------------------------------|:----------------|
| ifig           | figure index(ifig and axs cannot be used at the same time.)       | None            |
| axs            | Arranging an axis or subplot in a graph or plot                   | None            |
| proj           | proj can be either 'pf' or 'ipf'                                  | 'pf'            |
| poles          | Miller index of Polefigure                                        |[[1,0,0],[1,1,0]]|
| ix, iy         | x and y tick labels appended to each pole figure                  | ix='1',iy='2'   | 
| dph            | Grid of tilting angle                                             | 10              |
| dth            | Grid of in-plane rotation angle                                   | 10              |
| rot            | in-plane rotatation (radian)                                      |                 |
| n_rim          | The number of 'central' rims to be *averaged*                     | 2               |
| cdim           | crystal dimension                                                 | None            |
| ires           | If True, indicate the grid                                        | True            |
| mn             | Minimun level of contour                                          | None            |
| mx             | Maximun level of contour                                          | None            |
| lev_norm_log   | If True, use logarithmic scales. If False, linear scale.          | True            |
| nlev           | Level of iso contour bins                                         | 7               |
| cmap           | Color map used to color-code the contour levels.                  |                 |
| iline_khi80    | Whether or not to draw a line of chi=80                           |                 |
| mode           | ex. Contour modes, dot modes                                      | 'line'          |
| ilev           | Miller index of Polefigure                                        | 1               |
| levels         | level options have 0 or 1                                         | None            |
| transform      | transformation matrix applied to the entire polycrystal aggregate | 'magma'         |
| ideco_lev      | switch to turn on or off the levels                               | True            |
| ismooth        | Contours become smooth                                            | 1               |

![example image](README/Figure_1.png)
*Use arument poles only*
![example image](README/Figure_2.png)
*Use additional arguments mn, mx, ismooth*








# Pole figure 그리기

## 1. 기본 준비
 git hub로부터 texture3를 clone 혹은 download를 통해 설치
```shell
$ pip install .
```
texture3 폴더내의 ```setup.py```실행
## 2. 명령어 실행
```shell
$ ipython --pylab
```
cmd 창에서 파이썬실행
```python
from TX import upf
```
module ```TX```에서 ```upf```불러오기
```python
mypf = upf.polefigure(fnsx='fnsx file name', filename='TEX_PH1.OUT')
```
```mypf```에 데이터 저장, evpsc 시뮬레이션에 사용한 ```fnsx```파일, 시뮬레이션 이후  생성된 ```TEX_PH1.OUT```를 Arguments로 사용
```python
mypf.pf_new(poles=[[1,0,0],[1,1,0]],mn=0.5,mx=3.5,ismooth=10)
```
polefigure 생성
```pf_new```
| Arguments                |    Description                            |
|:-------------------------|:------------------------------------------|
|poles                     | Polefigure의 Miller index를 지정           |
|mn                        | 등고선의 최솟값을 지정                      |
|mx                        | 등고선의 최댓값을 지정                      |
|ismooth                   | 등고선이 각지지 않고 부드럽게 나타냄         |

![example image](README/Figure_1.png)
*poles arument만 사용*

![example image](README/Figure_2.png)
*mn, mx, ismooth arguments를 추가로 사용*
