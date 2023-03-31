import { Matrix4, Group, Vector2, Vector3, TextureLoader, BufferGeometry, Float32BufferAttribute, Points, Raycaster, PointsMaterial } from 'three';

/**
	 * A class for complex numbers.
	 * @author makc
	 */
	class Complex {
		/** Real part. */
		x;

		/** Imaginary part. */
		y;

		/** Principal argument. */
		get arg () {
			return Math.atan2 (this.y, this.x);
		}

		set arg (angle) {
			const _r = this.r; x = _r * Math.cos (angle); y = _r * Math.sin (angle);
		}

		/** Modulus. */
		get r () {
			return Math.sqrt (this.x ** 2 + this.y ** 2);
		}

		set r (length) {
			const scale = length / this.r; this.x *= scale; this.y *= scale;
		}

		/** Returns the conjugate. */
		conj () {
			return Complex.C (this.x, -this.y);
		}

		/** Returns the sum of this number and the argument. */
		add (b) {
			return Complex.C (this.x + b.x, this.y + b.y);
		}

		/** Returns the difference of this number and the argument. */
		sub (b) {
			return Complex.C (this.x - b.x, this.y - b.y);
		}

		/** Returns the product of this number and the argument. */
		mul (b) {
			return Complex.C (this.x * b.x - this.y * b.y, this.y * b.x + this.x * b.y);
		}

		/** Returns the ratio of this number and the argument. */
		div (b) {
			const d = b.x * b.x + b.y * b.y; return Complex.C ((this.x * b.x + this.y * b.y) / d, (this.y * b.x - this.x * b.y) / d);
		}

		/** Returns integer power of this number. */
		ipow (p) {
			const c = Complex.C (1, 0);
			if (p > -1) {
				for (let i = 0; i < p; i++) {
					const _x = this.x * c.x - this.y * c.y;
					const _y = this.y * c.x + this.x * c.y;
					c.x = _x;
					c.y = _y;
				}
				return c;
			}
			return c.div (this).ipow (-p);
		}

		/** Returns exponential function value for this number. */
		exp () {
			const ex = Math.exp (this.x); return Complex.C (ex * Math.cos (this.y), ex * Math.sin (this.y));
		}

		/** Returns the logarithm of this number to base e. */
		ln (k = 0) {
			return Complex.C (Math.log (this.r), this.arg + 2 * Math.PI * k);
		}

		/** Returns complex power of this number. */
		pow (p, k = 0) {
			return p.mul (this.ln (k)).exp ();
		}

		/** Returns the logarithm of this number to complex base. */
		log (b, k = 0, kb = 0) {
			return this.ln (k).div (b.ln (kb));
		}

		/**
		 * Creates permanent instance.
		 * @param	x Real part.
		 * @param	y Imaginary part.
		 */
		constructor (x = 0, y = 0) {
			this.x = x;
			this.y = y;
		}

		static #pool = [];
		static #index = 0;

		/**
		 * Returns temporary instance.
		 * You could use it for a constant until ReleaseTemporaries is called.
		 * @param	x Real part.
		 * @param	y Imaginary part.
		 */
		static C (x = 0, y = 0) {
			if (Complex.#index === Complex.#pool.length) Complex.#pool [Complex.#index] = new Complex;
			const c = Complex.#pool [Complex.#index++]; c.x = x; c.y = y; return c;
		}

		/**
		 * Returns all temporary instances back to the pool.
		 */
		static ReleaseTemporaries () {
			Complex.#index = 0;
		}

		/**
		 * Makes this instsance a copy of given number.
		 * @param	from Input instance to copy.
		 * @return This instance, modified.
		 */
		copy (from) {
			this.x = from.x; this.y = from.y; return this;
		}

		/**
		 * Saves calculation result (or makes a copy of this number).
		 * @param	result Output instance; will be created if not provided.
		 * @param	releaseTemporaries Returns all temporary instances back to the pool.
		 * @return Output instance.
		 */
		save (result = null, releaseTemporaries = false) {
			result = result || new Complex;
			result.copy (this);
			if (releaseTemporaries) Complex.ReleaseTemporaries ();
			return result;
		}

		/**
		 * Returns string representation of this number.
		 * @param	p Required precision.
		 */
		toString (p = 4) {
			return this.x.toPrecision (p) + " " + ((this.y > 0) ? "+" : "") + this.y.toPrecision (p) + "i";
		}
	}

const A$1 = new Complex;
const B$1 = new Complex;
const C$1 = new Complex;
const D$1 = new Complex;

const p = new Complex;
const q = new Complex;
const r = new Complex;
const s = new Complex;

const P = new Complex;
const Q = new Complex;
const R = new Complex;
const S = new Complex;

/** z⁴ + Az³ + Bz² + Cz + D = 0 */
function f ( z ) {
    return ( z.ipow (4) )
      .add ( z.ipow (3).mul (A$1) )
      .add ( z.ipow (2).mul (B$1) )
      .add ( z.mul (C$1) )
      .add ( D$1 );
}

function safeDiv ( a, b ) {
    const c = a.div ( b );
    if( isNaN (c.x) ) {
        c.x = Math.random () - Math.random ();
        c.y = Math.random () - Math.random ();
    }
    return c;
}

/**
 * Solves 4th power polynomials using Durand-Kerner method.
 * @see https://en.wikipedia.org/wiki/Durand-Kerner_method#Explanation
 */
function solve4 (a, b, c, d, releaseTemporaries = true) {
    // accept duck-typed arguments
    A$1.copy (a); B$1.copy (b); C$1.copy (c); D$1.copy (d);

    const seed = Complex.C (0.4, 0.9);

    seed.ipow (0).save (p);
    seed.ipow (1).save (q);
    seed.ipow (2).save (r);
    seed.ipow (3).save (s);

    while (true) {
        // P = p-f(p)/((p-q)(p-r)(p-s))
        p.sub (
            safeDiv (f (p),
                ( p.sub (q) ).mul( p.sub (r) ).mul( p.sub (s) )
            )
        ).save (P);

        // Q = q-f(q)/((q-p)(q-r)(q-s))
        q.sub (
            safeDiv (f (q),
                ( q.sub (p) ).mul( q.sub (r) ).mul( q.sub (s) )
            )
        ).save (Q);

        // R = r-f(r)/((r-p)(r-q)(r-s))
        r.sub (
            safeDiv (f (r),
                ( r.sub (p) ).mul( r.sub (q) ).mul( r.sub (s) )
            )
        ).save (R);

        // S = s-f(s)/((s-p)(s-q)(s-r))
        s.sub (
            safeDiv (f (s),
                ( s.sub (p) ).mul( s.sub (q) ).mul( s.sub (r) )
            )
        ).save (S);

        if (p.sub (P).r + q.sub (Q).r + r.sub (R).r + s.sub (S).r < 1e-8) {
            // converged
            break;
        } else {
            // on to next iteration
            P.save (p); Q.save (q); R.save (r); S.save (s, releaseTemporaries);
        }
    }

    if (releaseTemporaries) {
        // we no longer need temporary variables - return them all to the pool
        Complex.ReleaseTemporaries ();
    }

    return [P, Q, R, S];
}

const metric = new WeakMap();

function Pose(k1, k3, q) {
    const p = [];
    p[0] = q[0].clone();
    p[1] = q[1].clone().multiplyScalar (k1);
    p[3] = q[3].clone().multiplyScalar (k3);

    // p2 = p0 + (p1 - p0) + (p3 - p0) = p1 + p3 - p0
    p[2] = p[1].clone().add(p[3]).sub(p[0]);

    // reproject q2 from p2
    var zoom = q[2].z / p[2].z;
    var q2_x = p[2].x * zoom;
    var q2_y = p[2].y * zoom;

    metric.set(p,
        1 / ( 1e-6 +
            // reprojection error
            (q2_x - q[2].x) ** 2 + (q2_y - q[2].y) ** 2
        )
    );

    return p;
}

Pose.Best = function (a, b) {
    if(!a || (b && (metric.get(a) < metric.get(b)))) return +1;
    if(!b || (a && (metric.get(a) > metric.get(b)))) return -1;
    return 0;
};

const A = { x: 0, y: 0 }, B = { x: 0, y: 0 }, C = { x: 0, y: 0 }, D = { x: 0, y: 0 };

function exact3(q) {
    // qi are points in a camera plane
    // pi are quad vertices in 3D
    // we seek a set of coefs k0...k3 such that pi = ki * qi
    // we fix k0 to be 1 because scale is arbitrary
    // we find k3 as a function of k1 from (p1-p0) dot (p3-p0) = 0 [1]
    // and then from (p1-p0) dot (p1-p0) = (p3-p0) dot (p3-p0) [2]
    // which gives this 4th degree polynomial:
    // 0 = -e*c^2*k1^4
    // +2*c*(a*c + e*d)*k1^3
    // +(a*(f*a -6*c*d) -e*d^2)*k1^2
    // +2*(a*(2*d^2 -f*b) +b*c*d)*k1
    // +b*(f*b -2*d^2)
    const a = q[0].dot (q[1]);
    const b = q[0].dot (q[0]);
    const c = q[3].dot (q[1]);
    const d = q[3].dot (q[0]); // k3 = (a k1 - b)/(c k1 - d) [1]
    const e = q[1].dot (q[1]);
    const f = q[3].dot (q[3]); // k1^2 * e - 2a * k1 = k3^2 * f - 2d * k3 [2]

    // simplify into k1^4 + Ak1^3 + Bk1^2 + Ck1 + D
    const g = -e * c * c;
    A.x = 2 * c * (a * c + e * d) / g;
    B.x = (a * (f * a - 6 * c * d) - e * d * d) / g;
    C.x = 2 * (a * (2 * d * d - f * b) + b * c * d) / g;
    D.x = b * (f * b - 2 * d * d) / g;

    const roots = solve4 (A, B, C, D);

    // threshold real roots to ensure at least one k1 value
    for (let i = 0; i < 4; i++) {
        roots[i].y = Math.abs (roots[i].y);
    }

    const threshold = Math.max (1e-3,
        Math.min (roots[0].y, roots[1].y, roots[2].y, roots[3].y)
    );

    return roots.map (function (r) {
        if (r.y < threshold) {
            const k1 = r.x;
            const k3 = (a * k1 - b) / (c * k1 - d);
            if ((k1 > 0) && (k3 > 0)) {
                return new Pose(k1, k3, q);
            }
        }
    }).sort (Pose.Best) [0];
}

exact3.average = function (q) {
    // this produces bogus pose by averaging exact3() results for each qi
    // however, exact3() results change once pi-s get reprojected, and so
    // average() result does, too - therefore it is repeated here a few
    // times, hoping that it will converge to something pretty
    let iterations = 3, qz = q[0].z, p;
    do {
        let r = [], w = 0;

        for (let i = 0; i < 4; i++) {
            let result = exact3 (q.map (function(_, j) { return q[(i + j) % 4] }));
            if (result) {
                // weigh results with their metrics
                let m = metric.get (result); w += m; r.push (result);

                for (let j = 0; j < 4; j++) {
                    result[j].multiplyScalar (m);
                }

                for (let j = 0; j < i; j++) {
                    result.unshift (result.pop ());
                }
            }
        }

        if (w > 0) {
            p = q.map (function(qi) {
                return qi.clone().multiplyScalar(0);
            });

            for (let i = 0; i < 4; i++) {
                for (let j = 0; j < r.length; j++) {
                    p[i].add (r[j][i].multiplyScalar (1 / w));
                }
            }

            // try to keep it square - the method here is also bogus, the only
            // sound condition it satisfies is that square input is not changed
            const T = p[0].clone ().add (p[1]).add (p[2]).add (p[3]).multiplyScalar (0.25);
            const TP0 = p[0].clone ().sub (T);
            const TP1 = p[1].clone ().sub (T);
            const Z = TP1.clone ().cross (TP0).normalize ();
            TP0.applyAxisAngle (Z, -Math.PI / 4);
            TP1.applyAxisAngle (Z,  Math.PI / 4);
            const Y = TP0.add (TP1).multiplyScalar (0.5 * Math.SQRT1_2); // Y candidate
            const X = TP1.copy (Y).applyAxisAngle (Z, -Math.PI / 2); // X candidate
            p[0].copy (T).add (Y).sub (X);
            p[1].copy (T).add (Y).add (X);
            p[2].copy (T).sub (Y).add (X);
            p[3].copy (T).sub (Y).sub (X);

            // reproject p and repeat
            q = p.map (function(pi) {
                return pi.clone ().multiplyScalar (qz / pi.z);
            });
        } else {
            // there's nothing to iterate
            break;
        }
    } while (iterations --> 0);

    return p;
};

const m = new Matrix4();

function exact4(q) {
    // qi are points in a camera plane
    // pi are quad vertices in 3D
    // we seek a set of coefs k0...k3 such that pi = ki * qi
    // we fix k0 to be 1 because scale is arbitrary
    // pi must be in the same plane, so p3 + p1 = p2 + p0
    // or -q3 * k3 +q2 * k2 -q1 * k1 +q0 * 1 = 0
    const p = q.map (function (qi) { return qi.clone(); });
    p[3].multiplyScalar (-1);
    p[1].multiplyScalar (-1);

    m.makeBasis (p[3], p[2], p[1]);
    m.setPosition (p[0]);
    m.invert ();

    for (let i = 3; i > 0; i--) {
        p[i].copy (q[i]).multiplyScalar (m.elements[15 - i]);
    }

    return p;
}

let cmat, dist, ipts, mpts, rvec, tvec, rmat;

window.cv?.then(function() {
    cmat = cv.Mat.zeros(3, 3, cv.CV_64FC1); // camera matrix
    dist = cv.Mat.zeros(4, 1, cv.CV_64FC1); // distortion coef-s
    ipts = cv.Mat.zeros(4, 2, cv.CV_64FC1); // image points
    mpts = cv.matFromArray(4, 3, cv.CV_64FC1, [ // model points
        -1,  1,  0,
         1,  1,  0,
         1, -1,  0,
        -1, -1,  0
    ]);
    rvec = new cv.Mat({ width: 1, height: 3 }, cv.CV_64FC1);
    tvec = new cv.Mat({ width: 1, height: 3 }, cv.CV_64FC1);
    rmat = new cv.Mat({ width: 3, height: 3 }, cv.CV_64FC1);
});

function opencv(q) {
    if (cmat) {
        cmat.data64F[0] = cmat.data64F[4] = -q[0].z; // fx, fy
        cmat.data64F[8] = 1;

        for (let i = 0; i < 4; i++) {
            ipts.data64F[2 * i] = q[i].x; ipts.data64F[2 * i + 1] = -q[i].y;
        }

        // initial values for *vec
        tvec.data64F[0] = 0; tvec.data64F[1] = 0; tvec.data64F[2] = -q[0].z * 2;
        rvec.data64F[0] = 0; rvec.data64F[1] = 0; rvec.data64F[2] = 0;

        try {
            if (cv.solvePnP (
                mpts,
                ipts,
                cmat,
                dist,
                rvec,
                tvec,
                true, // useExtrinsicGuess
                cv.SOLVEPNP_IPPE_SQUARE
            )) {
                // success
                cv.Rodrigues (rvec, rmat);

                const X = q[0].clone ().set (rmat.data64F[0], -rmat.data64F[3], -rmat.data64F[6]);
                const Y = q[0].clone ().set (rmat.data64F[1], -rmat.data64F[4], -rmat.data64F[7]);
                const T = q[0].clone ().set (tvec.data64F[0], -tvec.data64F[1], -tvec.data64F[2]);

                return [
                    T.clone ().add (Y).sub (X),
                    T.clone ().add (Y).add (X),
                    T.clone ().sub (Y).add (X),
                    T.clone ().sub (Y).sub (X)
                ];
            }
        } catch (oopsie) {
            // would you look at that - opencv throws errors sometimes
        }
    }
}

function FourPointsControls(perspectiveCamera, domElement) {
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

export { FourPointsControls };
