
use nalgebra::SVector;

fn lower_triangle<const N: usize>(m : &nalgebra::SMatrix<f64, N, N>, b_vec : &nalgebra::SVector<f64, N>) ->
(nalgebra::SMatrix<f64, N, N>,nalgebra::SVector<f64, N>) {
    
    let mut mat = nalgebra::SMatrix::<f64, N, N>::zeros();
    mat.copy_from(m);
    let mut b = nalgebra::SVector::<f64, N>::zeros();
    b.copy_from(b_vec);

    let epsilon = 1e-8;
    let mut pivot = (0,0);
    
    for _ in 0..N {
        if mat[pivot] != 0.0 {
            for i in pivot.0+1..N {
                let factor = mat[(i, pivot.1)] / mat[pivot];

                let mut row = mat.row(i).into_owned();
                row -= mat.row(pivot.0)*factor;
                // let row = row.map(|x| if x.abs() < epsilon { 0.0 } else { x });
                mat.row_mut(i).copy_from(&row);

                b[i] -= factor * b[pivot.0];
            }
            println!("{} {}", mat, b);

            pivot.0 += 1;
            pivot.1 += 1;
        } else {
            let mut found = false;
            for i in 1..N-pivot.0 {
                if mat[(pivot.0+i, pivot.1+i)].abs() > epsilon {
                    mat.swap_rows(pivot.0, pivot.0+i);
                    b.swap_rows(pivot.0, pivot.0+i);
                    found = true;
                }
            }
            if !found {
                panic!("singular matrix")
            }
        }
    }
    return (mat, b);
}

fn gaussian<const N: usize>(m: &nalgebra::SMatrix<f64, N, N>, b_vec : &nalgebra::SVector<f64, N>) -> nalgebra::SVector<f64, N>{

    let mut mat = nalgebra::SMatrix::<f64, N, N>::zeros();
    mat.copy_from(m);
    let mut b = nalgebra::SVector::<f64, N>::zeros();
    b.copy_from(b_vec);

    let (mat, b) = lower_triangle::<N>(&mat, &b);

    let mut sol = SVector::<f64, N>::zeros();

    for i in (0..N).rev() {
        let mut sum = 0.0;
        for j in i+1..N {
            sum += mat[(i, j)] * sol[j];
        }
        sol[i] = (b[i] - sum) / mat[(i, i)];
    }
    
    println!("backsub {} {} {}",mat,b, sol);
    return sol;
}


pub fn gaussian_main() {

    let mat = nalgebra::SMatrix::<f64, 3, 3>::new(
        2.0, 6.0, 3.0,
        4.0, 6.0, 6.0,
        7.0, 8.0, 9.0
    );
    let b = nalgebra::SVector::<f64, 3>::new(1.0, 2.0, 3.0);
    let sol = gaussian::<3>(&mat, &b);
    println!("{} {} {}", mat, b, mat*sol);
}


