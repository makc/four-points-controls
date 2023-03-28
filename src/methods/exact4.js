import { Matrix4 } from 'three';

const m = new Matrix4();

export default function exact4(q) {
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