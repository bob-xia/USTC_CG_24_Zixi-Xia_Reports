# Overview of My Course Study of Computer Graphics
###### Zixi Xia, teacher: Ligang Liu

##### Tips: To see the .gif animated pictures, please click [https://github.com/bob-xia/USTC_CG_24_Zixi-Xia_Reports/blob/main/brief_intro.md](https://github.com/bob-xia/USTC_CG_24_Zixi-Xia_Reports/blob/main/brief_intro.md).

## Selected Homeworks

### Hw2: Image Warping

- Two types of interpolation function: IDW (Inverse distance-weighted) and RBF (Radial basis functions);
- Method to fill the black seams: reversed transform, average number.

| Method            | Control Points         | IDW  | RBF-1 | RBF-2 |
| ----------------- | ---------------------- | ---- | ----- | ----- |
| No filling        | ![](figures/1.1.1.png) | ![](figures/1.1.2.png) | ![](figures/1.1.3.png)|![](figures/1.1.4.png)|
| Reverse Transform | ![](figures/1.2.1.png) | ![](figures/1.2.2.png) | ![](figures/1.2.3.png)|![](figures/1.2.4.png)|

### Hw3: Poisson Image Editing [Perez et al., Siggraph 2003]

- Seamless image cloning by solving Poisson Equation;
- Two major types of algorithm: Seamless and Mixed.

|            Seamless            |             Mixed              |
| :----------------------------: | :----------------------------: |
| <img src="figures/2.1.png"  /> | <img src="figures/2.1.png"  /> |

### Hw4: Tutte Parameterization

- Minimal Surface: different weight-uniform weights, cotangent weights, Floater's shape-preserving weights;
- Parameterization: minimal surface, where the boundary are set to a plane.

|                        |                  Uniform weights                   |                 Cotangent weights                  |              Shape-preserving weights              |
| :--------------------: | :------------------------------------------------: | :------------------------------------------------: | :------------------------------------------------: |
| map boundary to circle | <img src="figures/3.1.1.png" style="zoom: 20%;" /> | <img src="figures/3.1.2.png" style="zoom: 20%;" /> | <img src="figures/3.1.3.png" style="zoom: 20%;" /> |
| map boundary to square | <img src="figures/3.2.1.png" style="zoom: 20%;" /> | <img src="figures/3.2.2.png" style="zoom: 20%;" /> | <img src="figures/3.2.3.png" style="zoom: 20%;" /> |

|        |                  Uniform weights                  |                 Cotangent weights                 |             Shape-preserving weights              |
| :----: | :-----------------------------------------------: | :-----------------------------------------------: | :-----------------------------------------------: |
| Circle | <img src="figures/3.3.1.gif" style="zoom:20%;" /> | <img src="figures/3.3.2.gif" style="zoom:20%;" /> | <img src="figures/3.3.3.gif" style="zoom:20%;" /> |
| Square | <img src="figures/3.4.1.gif" style="zoom:20%;" /> | <img src="figures/3.4.2.gif" style="zoom:20%;" /> | <img src="figures/3.4.3.gif" style="zoom:20%;" /> |

### Hw5: ARAP Parameterization [Liu et al., Siggraph 2008]

- The relationship between texture coordination and its current position came up to two types: ARAP (As-rigid-as-possible) and ASAP (As-similar-as-possible).
- Solve nonlinear equations by Local/Global Approach.

|                 |                       ARAP                        |                   ARAP-filling                    |                       ASAP                        |                   ASAP-filling                    |
| :-------------: | :-----------------------------------------------: | :-----------------------------------------------: | :-----------------------------------------------: | :-----------------------------------------------: |
|  `Balls.usda`   | <img src="figures/4.1.1.png" style="zoom:33%;" /> | <img src="figures/4.1.2.png" style="zoom:33%;" /> | <img src="figures/4.1.3.png" style="zoom:33%;" /> | <img src="figures/4.1.4.png" style="zoom:33%;" /> |
| `Cow_dABF.usda` | <img src="figures/4.2.1.png" style="zoom:33%;" /> | <img src="figures/4.2.2.png" style="zoom:33%;" /> | <img src="figures/4.2.3.png" style="zoom:33%;" /> | <img src="figures/4.2.4.png" style="zoom:33%;" /> |

|                     |          ARAP          |                    ARAP-filling                    |          ASAP          |      ASAP-filling      |
| :-----------------: | :--------------------: | :------------------------------------------------: | :--------------------: | :--------------------: |
| `Gargoyle_ABF.usda` | ![](figures/4.3.1.png) |               ![](figures/4.3.2.png)               | ![](figures/4.3.3.png) | ![](figures/4.3.4.png) |
|  `Isis_dABF.usda`   | ![](figures/4.4.1.png) | <img src="figures/4.4.2.png" style="zoom: 50%;" /> | ![](figures/4.4.3.png) | ![](figures/4.4.4.png) |

### Hw6: Shader Programming

- Realization of Blinn-Phong lighting model, as well as normal map, and gamma correlation. 

	| Normal by Interpolation | Normal Map             |
	| ----------------------- | ---------------------- |
	| ![](figures/5.1.1.png)  | ![](figures/5.1.2.png) |

- Realization of Shadow Mapping, and soften the edge of shadow.

	| <img src="figures/5.2.1.gif" style="zoom: 33%;" /> | ![](figures/5.2.2.png) |
	| :--: | ---- |
	| ![](figures/5.2.3.png)  | ![](figures/5.2.4.png) |


### Hw8: Mass Spring Model

- Mass-Spring Model with realization of Euler equation (explicit, implicit, semi-implicit);
- Fast simulation method of Mass-Spring model [Liu et al., Siggraph 2008].

| Different Mesh Density: 10x10 | 20x20                  | 40x40                  |
| ----------------------------- | ---------------------- | ---------------------- |
| ![](figures/6.1.1.gif)        | ![](figures/6.1.2.gif) | ![](figures/6.1.3.gif) |

| Fixed points: A vertex of the square | Diagonal vertices of the square | Four vertices of the square |
| ------------------------------------ | ------------------------------- | --------------------------- |
| ![](figures/6.2.1.gif)               | ![](figures/6.2.2.gif)          | ![](figures/6.2.3.gif)      |

| External force $a_w=10r=10||X_i||$ | External force $a_w=\cfrac{-(X_i-X_0)}{||X_i-X_0||^3}$，where $X_0=(0,0,-2)$ |
| ---------------------------------- | ------------------------------------------------------------ |
| ![](figures/6.3.1.gif)             | ![](figures/6.3.2.gif)                                       |

| Penalty Force method for collision | Optimized method       | Example: collision between two balls |
| ---------------------------------- | ---------------------- | ------------------------------------ |
| ![](figures/6.4.1.gif)             | ![](figures/6.4.2.gif) | ![](figures/6.4.3.gif)               |

### Hw9: SPH fluid

- Solving Navier-Stokes equation by Operator Splitting, as well as getting physical quantities from neighboring particles, called WCSPH.
- Realization of incompressible fluid, which is called IISPH.

|                      WCSPH                       |                      IISPH                       |
| :----------------------------------------------: | :----------------------------------------------: |
| <img src="figures/7.1.gif" style="zoom: 33%;" /> | <img src="figures/7.2.gif" style="zoom: 50%;" /> |

## Final Project

### Simulation of visual effect with Einstein's Special Theory of Relativity, Author: Yutian Zhu, Zichao Liu, Zixi Xia

#### Results and Report

##### See [https://github.com/bob-xia/USTC_CG_24_Zixi-Xia_Reports/tree/main/project](https://github.com/bob-xia/USTC_CG_24_Zixi-Xia_Reports/tree/main/project) for seven video tapes and report (Chinese version and English version translated by ChatGPT).

