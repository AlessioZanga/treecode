pub type Real = f32;

pub fn rsqr(x: Real) -> Real {
    x * x
}

pub fn rqbe(x: Real) -> Real {
    x * x * x
}

pub fn rlog2(x: Real) -> Real {
    x.log2()
}

pub fn rexp2(x: Real) -> Real {
    x.exp2()
}

pub fn rdex(x: Real) -> Real {
    10.0f32.powf(x)
}

pub fn fcbrt(x: f32) -> f32 {
    x.cbrt()
}

pub fn xrandom(xl: f64, xh: f64) -> f64 {
    extern "C" {
        fn random() -> libc::c_long;
    }
    xl + (xh - xl) * (unsafe { random() } as f64) / 2147483647.0
}

pub fn grandom(mean: f64, sdev: f64) -> f64 {
    let mut v1;
    let mut v2;
    let mut s;
    loop {
        v1 = xrandom(-1.0, 1.0);
        v2 = xrandom(-1.0, 1.0);
        s = v1 * v1 + v2 * v2;
        if s < 1.0 {
            break;
        }
    }
    mean + sdev * v1 * (-2.0 * s.ln() / s).sqrt()
}

pub fn pickshell(vec: &mut [Real], ndim: usize, rad: Real) {
    loop {
        let mut rsq: Real = 0.0;
        for v in vec.iter_mut().take(ndim) {
            *v = xrandom(-1.0, 1.0) as Real;
            rsq += *v * *v;
        }
        if rsq <= 1.0 {
            let rscale = rad / rsq.sqrt();
            for v in vec.iter_mut().take(ndim) {
                *v *= rscale;
            }
            return;
        }
    }
}

pub fn pickball(vec: &mut [Real], ndim: usize, rad: Real) {
    loop {
        let mut rsq: Real = 0.0;
        for v in vec.iter_mut().take(ndim) {
            *v = xrandom(-1.0, 1.0) as Real;
            rsq += *v * *v;
        }
        if rsq <= 1.0 {
            for v in vec.iter_mut().take(ndim) {
                *v *= rad;
            }
            return;
        }
    }
}

pub fn pickbox(vec: &mut [Real], ndim: usize, size: Real) {
    for v in vec.iter_mut().take(ndim) {
        *v = xrandom(-(size as f64), size as f64) as Real;
    }
}
