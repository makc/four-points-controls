import {
    BufferGeometry,
    Float32BufferAttribute,
    Group,
    Points,
    PointsMaterial,
    Raycaster,
    TextureLoader,
    Vector2,
    Vector3
} from 'three';

import exact3 from './methods/exact3';
import exact4 from './methods/exact4';
import opencv from './methods/opencv';

export function FourPointsControls(perspectiveCamera, domElement) {
    perspectiveCamera.clearViewOffset();

    const object = new Group();
    object.clone = function() {
        return FourPointsControls(perspectiveCamera, domElement);
    };

    object.points = [
        new Vector2( 0.0,  0.0), new Vector2( 0.6,  0.1),
        new Vector2( 0.8, -0.5), new Vector2(-0.2, -0.5)
    ];

    const handles = [];
    const Q = object.points.map(function() { return new Vector3() });
    const X = new Vector3(), Y = new Vector3(), Z = new Vector3(), T = new Vector3();

    object.method = exact4;
    object.currentPoint = 1;

    let needsReprojection = false;

    object.updateMatrix = function() {
        
        const f = 1 / Math.tan(0.5 * perspectiveCamera.fov * Math.PI / 180);

        Q.forEach(function(point, index) {
            const sourceIndex = (index + object.currentPoint) % 4;
            point.x = object.points[sourceIndex].x * perspectiveCamera.aspect;
            point.y = object.points[sourceIndex].y;
            point.z = -f;
        });

        const P = this.method(Q);

        if(!P) {
            console.warn('Failed to find a pose for', JSON.stringify(object.points), 'current point:', object.currentPoint);
            return;
        }

        for(let i = 0; i < object.currentPoint; i++) {
            P.unshift(P.pop());
            Q.unshift(Q.pop());
        }

        T.copy(P[0]).add(P[2]).multiplyScalar(0.5);
        X.copy(P[1]).add(P[2]).multiplyScalar(0.5).sub(T);
        Y.copy(P[0]).add(P[1]).multiplyScalar(0.5).sub(T);

        // re-scale back to 1
        const s = 2 / (X.length() + Y.length());
        T.multiplyScalar(s);
        X.multiplyScalar(s);
        Y.multiplyScalar(s);

        Z.copy(X).cross(Y);


        this.matrix.makeBasis(X, Y, Z);
        this.matrix.setPosition(T);

        perspectiveCamera.updateMatrixWorld();
        this.matrix.premultiply(perspectiveCamera.matrix);

        this.matrixWorldNeedsUpdate = true;

        // update transformation params for the user's reference 
        this.matrix.decompose(this.position, this.quaternion, this.scale);


        if(handles.length) {
            handles.forEach(function(point, index) {
                if(needsReprojection) {
                    const zoom = -f / P[index].z;
                    Q[index].x = zoom * P[index].x; object.points[index].x = Q[index].x / perspectiveCamera.aspect;
                    Q[index].y = zoom * P[index].y; object.points[index].y = Q[index].y;
                }
                point.position.copy(Q[index]).multiplyScalar(
                    (1.000001 * perspectiveCamera.near) / f
                );
                point.updateMatrix();
                point.matrix.premultiply(perspectiveCamera.matrix);
                point.matrixWorld.copy(point.matrix);
            });

            needsReprojection = false;
        }
    };

    if(domElement) {
        const loader = new TextureLoader();
        const makeMaterial = function(url) {
            return new PointsMaterial({ size: 16, sizeAttenuation: false, transparent: true, depthTest: false, map: loader.load(url) });
        };

        const material = makeMaterial("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32' width='32' height='32'%3E%3Ccircle cx='16' cy='16' r='13.5' stroke='%23F80' stroke-width='3' fill='%23F80' fill-opacity='.8' /%3E%3C/svg%3E%0A");
        const materialActive = makeMaterial("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32' width='32' height='32'%3E%3Ccircle cx='16' cy='16' r='15' fill='%23FC0' /%3E%3C/svg%3E%0A");

        const geometry = new BufferGeometry();
        geometry.setAttribute('position', new Float32BufferAttribute([0, 0, 0], 3));

        object.points.forEach(function(v, index) {
            const point = new Points(geometry, material);
            point.frustumCulled = false;
            //point.position.set([-1, 1, 1, -1][index], [1, 1, -1, -1][index], 0);
            point.updateMatrixWorld = function() {};
            object.add(point);
            handles.push(point);
        });

        const event2xy = function(event) {
            const rect = domElement.getBoundingClientRect();
            return {
                x: ((event.clientX - rect.left) / rect.width) * 2 - 1,
                y: -((event.clientY - rect.top) / rect.height) * 2 + 1
            }
        };

        let pointerIsDown = false;
        domElement.addEventListener('pointerdown', function(event) {
            pointerIsDown = true;
        });
        domElement.addEventListener('pointerup', function() {
            pointerIsDown = false;
            needsReprojection = true;
        });

        const raycaster = new Raycaster();
        raycaster.params.Points.threshold =
            // try to come up with something reasonable here...
            20 * perspectiveCamera.near * Math.tan(0.5 * perspectiveCamera.fov * Math.PI / 180) / domElement.offsetHeight;

        domElement.addEventListener('pointermove', function(event) {
            const xy = event2xy(event);
            if(pointerIsDown) {
                // drag
                const index = handles.findIndex(function(point) {
                    return point.material === materialActive;
                });
                if(index > -1) {
                    object.points[object.currentPoint = index].copy(xy);
                }
            } else {
                // hover
                handles.forEach(function(point) {
                    point.material = material;
                });
                domElement.style.cursor = '';

                raycaster.setFromCamera(event2xy(event), perspectiveCamera);
                const results = raycaster.intersectObjects(handles);
                if(results[0]) {
                    results[0].object.material = materialActive;
                    domElement.style.cursor = 'move';
                }
            }
        });

        needsReprojection = true;
    }

    return object;
}

FourPointsControls.exact3 = exact3;
FourPointsControls.exact4 = exact4;
FourPointsControls.opencv = opencv;