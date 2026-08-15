pub fn rsqr(x: f32) -> f32 {
    x * x
}

pub fn rqbe(x: f32) -> f32 {
    x * x * x
}

pub fn rlog2(x: f32) -> f32 {
    x.log2()
}

pub fn rexp2(x: f32) -> f32 {
    x.exp2()
}

pub fn rdex(x: f32) -> f32 {
    10.0f32.powf(x)
}

pub fn fcbrt(x: f32) -> f32 {
    x.cbrt()
}
