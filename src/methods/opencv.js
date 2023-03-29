
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

export default function opencv(q) {
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