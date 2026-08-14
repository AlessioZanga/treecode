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
