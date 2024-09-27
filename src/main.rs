use std::fmt::{Display, Formatter, Write};

fn main() {
    let a = Matrix::from_rows(vec![
        vec![1.0, 1.0, 0.0],
        vec![-2.0, 2.0, 1.0],
        vec![0.0, 1.0, -1.0],
    ]);
    // let a = Matrix::from_rows(vec![
    //     vec![0.0, 1.0, 0.0, 0.0],
    //     vec![2.0, 3.0, 2.0, 0.0],
    //     vec![0.0, 1.0, 1.0, 0.0],
    //     vec![0.0, 0.0, 0.0, 1.0],
    // ]);
    println!("{}", a);
    println!("{}", a.determinant());
}

struct Matrix {
    v: Vec<Vec<f64>>,
    m: usize,
    n: usize,
    aug: Option<usize>,
}

#[allow(dead_code)]
impl Matrix {
    fn from_rows(v: Vec<Vec<f64>>) -> Self {
        let m = v.len();
        assert!(m != 0);
        let n = v[0].len();
        assert!(v.iter().all(|row| row.len() == n));
        Self { v, m, n, aug: None }
    }

    fn identity(m: usize, n: usize) -> Matrix {
        let mut v = vec![vec![0.0; n]; m];
        for i in 0..(std::cmp::min(m, n)) {
            v[i][i] = 1.0;
        }
        Matrix { v, m, n, aug: None }
    }

    fn transpose(&self) -> Matrix {
        let m = self.n;
        let n = self.m;
        let mut v = vec![vec![0.0; n]; m];
        for (i, row) in self.v.iter().enumerate() {
            for (j, a) in row.iter().enumerate() {
                v[j][i] = *a;
            }
        }
        Matrix { v, m, n, aug: None }
    }

    fn minor(&self, i: usize, j: usize) -> Matrix {
        assert!(self.m > 1);
        assert!(self.n > 1);
        let m = self.m - 1;
        let n = self.n - 1;
        let mut v = self.v.clone();
        v.remove(i);
        v.iter_mut().for_each(|row| _ = row.remove(j));
        Matrix { v, m, n, aug: None }
    }

    fn augment(&self, a: Matrix) -> Matrix {
        assert!(a.m == self.m);
        let m = a.m;
        let n = self.n + a.n;
        let mut v = Vec::with_capacity(m);
        for i in 0..m {
            let mut row = Vec::with_capacity(n);
            row.extend_from_slice(&self.v[i]);
            row.extend_from_slice(&a.v[i]);
            v.push(row);
        }
        let aug = Some(self.n);
        Matrix { v, m, n, aug }
    }

    fn determinant(&self) -> f64 {
        assert!(self.n > 0);
        assert!(self.m == self.n);
        if self.n == 1 {
            return self.v[0][0];
        }
        self.v[0]
            .iter()
            .enumerate()
            .map(|(j, a)| {
                let cofactor = {
                    let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
                    sign * self.minor(0, j).determinant()
                };
                a * cofactor
            })
            .sum()
    }

    fn swap_rows(&mut self, row_i: usize, row_j: usize) {
        assert!(row_i < self.m);
        assert!(row_j < self.m);
        self.v.swap(row_i, row_j);
    }

    fn scale_row(&mut self, i: usize, scalar: f64) {
        assert!(i < self.m);
        assert!(scalar != 0.0);
        for x in self.v[i].iter_mut() {
            *x *= scalar;
        }
    }

    fn add_scaled_row(&mut self, target_i: usize, source_i: usize, scalar: f64) {
        assert!(target_i < self.m);
        assert!(source_i < self.m);
        assert!(scalar != 0.0);
        assert!(target_i != source_i);
        for i in 0..self.m {
            self.v[target_i][i] += scalar * self.v[source_i][i];
        }
    }

    fn row_reduce(&mut self) {
        let mut h = 0;
        let mut k = 0;
        while h <= self.m && k <= self.n {
            if let Some(i_max) = self
                .v
                .iter()
                .map(|row| row[k])
                .enumerate()
                .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
                .map(|(i, _)| i)
            {
                self.swap_rows(h, i_max);
                for i in h + 1..self.m {
                    let _f = self.v[i][k] / self.v[h][k];
                    unimplemented!()
                    // https://en.wikipedia.org/wiki/Gaussian_elimination#Example_of_the_algorithm
                }
                h += 1;
            }
            k += 1;
        }
    }
}

impl Display for Matrix {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        if self.m == 0 {
            return write!(f, "[]",);
        }
        if self.m == 1 {
            write!(f, "[")?;
            for (j, x) in self.v[0].iter().enumerate() {
                if j != 0 {
                    write!(f, " ")?;
                }
                write!(f, "{}", x)?;
            }
            write!(f, "]")?;
            return Ok(());
        }
        let max_widths = self
            .v
            .iter()
            .try_fold(vec![0; self.n], |mut max_widths, row| {
                row.iter().enumerate().try_for_each(|(i, x)| {
                    let mut s = String::new();
                    if let Some(p) = f.precision() {
                        write!(&mut s, "{:.*}", p, x)?;
                    } else {
                        write!(&mut s, "{}", x)?;
                    }
                    let width = s.chars().count();
                    max_widths[i] = std::cmp::max(max_widths[i], width);
                    Ok(())
                })?;
                Ok(max_widths)
            })?;

        for (i, row) in self.v.iter().enumerate() {
            let (row_start, row_end) = match i {
                0 => ('┌', '┐'),
                j if j == self.m - 1 => ('└', '┘'),
                _ => ('│', '│'),
            };
            write!(f, "{}", row_start)?;
            for (j, x) in row.iter().enumerate() {
                let column_delimiter = if let Some(_) = self.aug.filter(|&x| x == j) {
                    " │ "
                } else {
                    " "
                };
                if j != 0 {
                    write!(f, "{}", column_delimiter)?;
                }
                if let Some(p) = f.precision() {
                    write!(f, "{:>width$.*}", p, x, width = max_widths[j],)?;
                } else {
                    write!(f, "{:>width$}", x, width = max_widths[j])?;
                }
            }
            writeln!(f, "{}", row_end)?;
        }

        Ok(())
    }
}
