### four-points-controls

Allows you to orient an object in 3D using positions of four points of its local XY plane in the camera plane. Get it from [NPM](https://www.npmjs.com/package/four-points-controls) or [CDN](https://cdn.jsdelivr.net/npm/four-points-controls/build/fpc.js) and then just:

```js
const controls = FourPointsControls( camera, renderer.domElement );
controls.add( mesh ); // your 3D object to orient
scene.add( controls );
```

#### Demo

https://makc.github.io/four-points-controls/

#### API

```js
// camera: an active instance of PerspectiveCamera
// element: optional HTMLElement to add pointer listeners to; if
// not specified, no handles are created and you are responsible
// for updating controls.points values
const controls = FourPointsControls( camera, element );

// four points in the camera plane, set by controls UI if element
// was specified, or by you (both x and y should be normalized to
// -1...1 range
controls.points = [...];

// currently active corner, 0...3 (see exact3 comment below)
const index = controls.currentPoint;

// a method of pose calculation (see below for supported values)
controls.method = ...;
```

#### Methods

* `FourPointsControls.exact3` \
calculates exact fit of three points and minimizes fourth point reprojection error (this is similar to opencv SOLVEPNP_P3P method); three points used are those around the active corner

* `FourPointsControls.exact3.average` \
*experimental*; calculates approximate fit of four points by averaging the results of `FourPointsControls.exact3` for every point

* `FourPointsControls.exact4` \
*default*; calculates exact fit of four points, but sacrifices the transformation orthonormality (this works best with flat or near-flat objects)

* `FourPointsControls.exact4.diamond` \
*experimental*; tweaks the results of `FourPointsControls.exact4` to equalize the lengths of X and Y axes

* `FourPointsControls.opencv` \
calculates approximate fit of four points using opencv SOLVEPNP_IPPE_SQUARE method (you will need to preload [~9 MB opencv.js build](https://docs.opencv.org/4.6.0/opencv.js) to use this option)
