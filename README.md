<p align="center">
<b>Usage</b>
</p>

1 - Prepare the workdir with `./run.sh`, or set up a box mesh yourself (explained in run.sh)

2 - Set up the field by editing  `scripts/StreamlineAux.py`, or write yourself a plugin to consume an external `.csv`

3 - Generate N streamlines by running `python3 scripts/Streamline.py N`

<p align="center">
<b>Example 1</b>
</p>

Given the following almost linear field

<img src="https://render.githubusercontent.com/render/math?math=V_x= 0.01 sin(\pi x)  cos(\pi y) + 0.1">

<img src="https://render.githubusercontent.com/render/math?math=V_y=-0.01  cos(\pi  x) sin(\pi  y) + 0.1">

<img src="https://render.githubusercontent.com/render/math?math=V_z=0.1">

<p align="center">
<img src="gallery/linearSource/toy.png" width="350" title="hover text">
</p>

Results visualized in [Meshlab](https://www.meshlab.net/) are

<div><img src="gallery/linearSource/1.png" height="150"><img src="gallery/linearSource/2.png" height="150"><img src="gallery/linearSource/3.png" height="150"><img src="gallery/linearSource/4.png" height="150">
</div>

<p align="center">
<b>Example 2</b>
</p>

Given the following almost 'whirly' field

<img src="https://render.githubusercontent.com/render/math?math=V_x=  sin(\pi x)  cos(\pi y) + 0.1">

<img src="https://render.githubusercontent.com/render/math?math=V_y=- cos(\pi  x) sin(\pi  y) + 0.1">

<img src="https://render.githubusercontent.com/render/math?math=V_z=0.1">

<p align="center">
<img src="gallery/oscillating/toy.png" width="350" title="hover text">
</p>

Results visualized in [Meshlab](https://www.meshlab.net/) are

<div><img src="gallery/oscillating/1.png" height="170"><img src="gallery/oscillating/2.png" height="170"><img src="gallery/oscillating/3.png" height="170">
</div>

