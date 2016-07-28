//! # bluenoisers
//!
//! Implementation of blue noise in Rust using Fast Poisson Disk Sampling.
//! This implementation can be used for images (2-dimensional), for volumes (3-dimensional) and in higher dimensions.
//!
//! For background information see here: https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
//!
//! Documentation is here: https://connyonny.github.io/bluenoisers/bluenoisers

extern crate rand;
use rand::Rng;

use std::cmp::min;

#[derive(Debug)]
struct BackgroundGrid {
    data: Vec<usize>,
    dimensions: Vec<f64>,
    min_dst_sqr: f64,
    cell_size: f64,
    cell_count: Vec<usize>,
    cell_multiplicators: Vec<usize>,
}

impl BackgroundGrid {
    pub fn new(dimensions: Vec<f64>, min_distance: f64) -> BackgroundGrid {
        assert!(min_distance > 0.0);
        let dimension = dimensions.len();
        let cell_size = min_distance / (dimension as f64).sqrt();
        let cell_count: Vec<usize> = dimensions.iter()
                                               .map(|x| (x / (cell_size as f64)).ceil() as usize)
                                               .collect();
        let data_size = cell_count.iter().fold(1_usize, |accu, x| accu * x);
        let mut cell_multiplicators = Vec::new();
        let mut multi_accu = 1_usize;
        for i in 0..dimension {
            cell_multiplicators.push(multi_accu);
            multi_accu *= cell_count[i];
        }
        BackgroundGrid {
            data: vec![0; data_size],
            dimensions: dimensions,
            min_dst_sqr: min_distance * min_distance,
            cell_size: cell_size,
            cell_count: cell_count,
            cell_multiplicators: cell_multiplicators,
        }
    }
    pub fn dst_sqr(x: &Vec<f64>, y: &Vec<f64>) -> f64 {
        debug_assert_eq!(x.len(), y.len());
        x.iter().zip(y.iter()).fold(0_f64, |accu, (xx, yx)| {
            let diff = xx - yx;
            accu + diff * diff
        })
    }
    fn calc_idx(&self, cell_id: &Vec<usize>) -> usize {
        self.cell_multiplicators
            .iter()
            .zip(cell_id.iter())
            .skip(1)
            .fold(cell_id[0], |accu, (multi, cell)| accu + multi * cell)
    }
    pub fn insert(&mut self,
                  sample_position: Vec<f64>,
                  samples: &mut Vec<Vec<f64>>)
                  -> Result<usize, ()> {
        if sample_position.iter()
                          .zip(self.dimensions.iter())
                          .any(|(samp_x, dim_x)| *samp_x < 0_f64 || samp_x >= dim_x) {
            return Err(());
        }
        let dimension = self.dimensions.len();
        debug_assert_eq!(sample_position.len(), dimension);
        let cell_id: Vec<usize> = sample_position.iter()
                                                 .map(|x| (*x / self.cell_size) as usize)
                                                 .collect();
        let samp_idx = self.calc_idx(&cell_id);
        debug_assert!(cell_id.iter().zip(self.cell_count.iter()).all(|(cid, cc)| cid < cc));
        let cell_offs = (self.min_dst_sqr / (self.cell_size as f64)).ceil() as usize;
        let min_cell: Vec<usize> = cell_id.iter().map(|x| x.saturating_sub(cell_offs)).collect();
        let max_cell: Vec<usize> = cell_id.iter()
                                          .zip(self.cell_count.iter())
                                          .map(|(x, size_x)| min(x + cell_offs, size_x - 1))
                                          .collect();
        debug_assert!(min_cell.iter()
                              .zip(max_cell.iter())
                              .zip(cell_id.iter())
                              .all(|((cmin, cmax), c)| cmin <= c && c <= cmax));
        let mut indices = min_cell.clone();
        let mut checked_own_idx = false;
        loop {
            debug_assert!(min_cell.iter()
                                  .zip(max_cell.iter())
                                  .zip(indices.iter())
                                  .all(|((cmin, cmax), c)| cmin <= c && c <= cmax));
            let idx = self.calc_idx(&indices);
            if idx == samp_idx {
                checked_own_idx = true;
            }
            match self.data[idx] {
                0 => (),
                other_id => {
                    let other_sample = &samples[other_id - 1];
                    if BackgroundGrid::dst_sqr(&sample_position, other_sample) < self.min_dst_sqr {
                        return Err(());
                    }
                }
            }
            // loop exit check
            if indices == max_cell {
                break;
            }
            // iterate indices
            for i in 0..dimension {
                if indices[i] == max_cell[i] {
                    indices[i] = min_cell[i];
                } else {
                    indices[i] += 1;
                    break;
                }
            }
        }
        // no collission found
        debug_assert!(checked_own_idx,
                      format!("Didn't check own idx.\n\tMin cells: {:?}\n\tMax cells: \
                               {:?}\n\tself cells: {:?}",
                              min_cell,
                              max_cell,
                              cell_id));
        samples.push(sample_position);
        debug_assert_eq!(self.data[samp_idx], 0);
        self.data[samp_idx] = samples.len();
        Ok(samples.len())
    }
}

fn polar_to_cartesian(radius: f64, angles: Vec<f64>) -> Vec<f64> {
    let mut cart = vec![0_f64; angles.len()+1];
    let sines: Vec<f64> = angles.iter().map(|x| x.sin()).collect();
    for (i, xi) in cart.iter_mut().enumerate() {
        // formula from https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates
        *xi = sines.iter().take(i).fold(radius, |accu, sine| accu * sine) *
              if let Some(ang) = angles.get(i) {
            ang.cos()
        } else {
            1_f64
        };
    }
    cart
}

/// Generates blue noise samples
///
/// # Arguments
///
/// * `dimensions` - How broad in each dimension the space to be filled is. Its length determines the dimensionality, i.e. length `n` generates an `n`-dimensional problem space.
/// * `min_distance` - How far away from each other should samples at least be, euclidean distance
/// * `k_abort` - How often should the generator try to generate a valid new neighbor of an existing point before giving that existing point up as starting point. A value of 30 is recommended.
///
/// The samples returned are in order of generation.
/// Each sample is at most `2 * min_distance` away from a previous sample (except the first sample, of course).
pub fn blue_noise(dimensions: Vec<f64>, min_distance: f64, k_abort: usize) -> Vec<Vec<f64>> {
    let dimension = dimensions.len();
    let mut rng = rand::thread_rng();
    let initial_sample = dimensions.iter().map(|x| rng.gen_range::<f64>(0_f64, *x)).collect();
    let mut samples: Vec<Vec<f64>> = Vec::new();
    let mut bggrid = BackgroundGrid::new(dimensions, min_distance);
    let initial_sample_id = bggrid.insert(initial_sample, &mut samples).unwrap();
    debug_assert_eq!(initial_sample_id, 1);
    let mut active: Vec<usize> = vec![initial_sample_id];
    while !active.is_empty() {
        let mut next_active: Vec<usize> = Vec::new();
        for current_id in active.iter() {
            let current_samp = samples[current_id - 1].clone();
            let mut found_something = false;
            for _ in 0..k_abort {
                let radius = rng.gen_range::<f64>(min_distance, 2_f64 * min_distance);
                let angles = (0..dimension - 1)
                                 .map(|_| {
                                     rng.gen_range::<f64>(0_f64, 2_f64 * std::f64::consts::PI)
                                 })
                                 .collect();
                let samp_offs = polar_to_cartesian(radius, angles);
                debug_assert_eq!(samp_offs.len(), dimension);
                // if polar_to_cartesian would return an iterator, this might be more efficient
                let samp = samp_offs.into_iter()
                                    .zip(current_samp.iter())
                                    .map(|(offs, x)| x + offs)
                                    .collect();
                match bggrid.insert(samp, &mut samples) {
                    Ok(new_samp_id) => {
                        next_active.push(new_samp_id);
                        found_something = true;
                    }
                    Err(_) => {
                        // wait for the next iteration
                    }
                }
            }
            if found_something {
                next_active.push(*current_id)
            }
        }
        active = next_active;
    }
    samples
}

#[test]
fn grid_corners() {
    let mut grid = BackgroundGrid::new(vec![35_f64, 9_f64], 4.0);
    let mut samples = Vec::new();
    assert_eq!(grid.cell_count.len(), 2);
    assert_eq!(grid.insert(vec![0., 9.], &mut samples), Err(()));
    assert_eq!(samples.len(), 0);
    assert_eq!(grid.insert(vec![0., 0.], &mut samples), Ok(1));
    assert_eq!(samples.len(), 1);
    assert_eq!(grid.insert(vec![34., 0.], &mut samples), Ok(2));
    assert_eq!(samples.len(), 2);
    assert_eq!(grid.insert(vec![0., 8.], &mut samples), Ok(3));
    assert_eq!(samples.len(), 3);
    assert_eq!(grid.insert(vec![34., 8.], &mut samples), Ok(4));
    assert_eq!(samples.len(), 4);
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate rand;
    use rand::Rng;
    use std::f64;
    #[test]
    fn sanity_3d() {
        sanity_nd(3, 15., 25.);
    }
    #[test]
    #[ignore]
    fn sanity_6d() {
        sanity_nd(6, 6., 6.5);
    }
    fn sanity_nd(dimension: usize, minr: f64, maxr: f64) {
        let mut rng = rand::thread_rng();
        let radius = 3.;
        let mut dimensions: Vec<f64> = Vec::new();
        for _ in 0..dimension {
            dimensions.push(rng.gen_range::<f64>(minr, maxr));
        }
        assert_eq!(dimensions.len(), dimension);
        let samples = blue_noise(dimensions, radius, 30);
        println!("there are {} samples.", samples.len());
        for s1 in samples.iter() {
            let mut mindst = f64::INFINITY;
            for s2 in samples.iter() {
                if s1 == s2 {
                    continue;
                }
                let dst = super::BackgroundGrid::dst_sqr(s1, s2).sqrt();
                if dst < mindst {
                    mindst = dst;
                }
            }
            assert!(mindst >= radius); // distance constraint violated
            assert!(mindst < 2_f64 * radius); // not nicely spread in the room
        }
    }
    fn get_image(radius: f64, size: usize) -> Vec<Vec<bool>> {
        let samples = blue_noise(vec![size as f64, size as f64], radius, 30);
        let mut image = vec![vec![false; size]; size];
        for s in samples {
            image[s[1] as usize][s[0] as usize] = true;
        }
        image
    }
    #[test]
    fn sanity_2d() {
        let size: isize = 128;
        let radius: isize = 8;
        let image = get_image(radius as f64, size as usize);
        for y in 0..size {
            for x in 0..size {
                if image[y as usize][x as usize] {
                    for dy in 0..radius {
                        for dx in 0..radius {
                            if dx == 0 && dy == 0 {
                                continue;
                            }
                            image.get((y + dy) as usize).map(|line| {
                                match line.get((x + dx) as usize) {
                                    Some(&true) => {
                                        // the -1 is to accomodate for rounding errors
                                        assert!(dx * dx + dy * dy >= (radius - 1) * (radius - 1));
                                    }
                                    _ => {}
                                }
                            });
                        }
                    }
                }
            }
        }
    }
}
