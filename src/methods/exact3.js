import { solve4 } from 'cdk4';

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

export default function exact3(q) {
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
}