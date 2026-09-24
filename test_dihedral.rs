fn dihedral(
    p1: (f64, f64, f64),
    p2: (f64, f64, f64),
    p3: (f64, f64, f64),
    p4: (f64, f64, f64),
) -> f64 {
    let b1 = (p2.0 - p1.0, p2.1 - p1.1, p2.2 - p1.2);
    let b2 = (p3.0 - p2.0, p3.1 - p2.1, p3.2 - p2.2);
    let b3 = (p4.0 - p3.0, p4.1 - p3.1, p4.2 - p3.2);

    let cross = |a: (f64, f64, f64), b: (f64, f64, f64)| -> (f64, f64, f64) {
        (
            a.1 * b.2 - a.2 * b.1,
            a.2 * b.0 - a.0 * b.2,
            a.0 * b.1 - a.1 * b.0,
        )
    };

    let dot = |a: (f64, f64, f64), b: (f64, f64, f64)| -> f64 {
        a.0 * b.0 + a.1 * b.1 + a.2 * b.2
    };

    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);

    // m1 = b2 x b3
    let m1 = cross(b2, b3);
    // n1 = b1 x b2
    // formula: atan2( |b2| * b1 . (b2 x b3) , (b1 x b2) . (b2 x b3) )
    let b2_mag = (b2.0 * b2.0 + b2.1 * b2.1 + b2.2 * b2.2).sqrt();
    let x = dot(n1, n2);
    let y = dot(b1, n2) * b2_mag;
    
    y.atan2(x)
}

fn main() {
    let a = dihedral((0.,0.,0.), (1.,0.,0.), (1.,1.,0.), (1.,1.,1.));
    println!("Angle in rad: {}, in deg: {}", a, a.to_degrees());
}
