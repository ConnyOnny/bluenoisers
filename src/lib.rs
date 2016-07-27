extern crate rand;
use rand::Rng;

use std::cmp::min;

#[derive(Clone,Copy,PartialEq,Eq,Debug)]
pub struct Sample2D {
    x: usize,
    y: usize,
}

impl Sample2D {
    pub fn dstsqr(&self, other: &Sample2D) -> usize {
        let dx = (self.x as isize - other.x as isize).abs() as usize;
        let dy = (self.y as isize - other.y as isize).abs() as usize;
        dx*dx+dy*dy
    }
}

struct BackgroundGrid2D {
    data: Vec<Option<usize>>, // stores indices
    cells_x: usize, // for linearization of the 2D array
    cells_y: usize,
    width: usize,
    height: usize,
    min_distance_sqr: usize,
    cell_size: usize,
}

impl BackgroundGrid2D {
    pub fn new(width: usize, height: usize, min_distance: f64) -> BackgroundGrid2D {
        assert!(min_distance > 0.0_f64);
        let max_cell_size = (min_distance / 2.0_f64.sqrt()).floor();
        let cells_x : usize = (width as f64 / max_cell_size).ceil() as usize;
        let cells_y : usize = (height as f64 / max_cell_size).ceil() as usize;
        BackgroundGrid2D {
            data: vec![None; cells_x*cells_y],
            cells_x: cells_x,
            cells_y: cells_y,
            width: width,
            height: height,
            min_distance_sqr: (min_distance * min_distance).ceil() as usize,
            cell_size: max_cell_size as usize,
        }
    }
    pub fn insert(&mut self, sample_x: usize, sample_y: usize, samples: &mut Vec<Sample2D>) -> Result<usize,()> {
        let new_sample = Sample2D { x: sample_x, y: sample_y };
        if sample_x >= self.width || sample_y >= self.height {
            return Err(());
        }
        let cell_id_x = sample_x / self.cell_size;
        let cell_id_y = sample_y / self.cell_size;
        assert!(cell_id_x < self.cells_x);
        assert!(cell_id_y < self.cells_y);
        let min_cell_id_x = cell_id_x.saturating_sub(2);
        let max_cell_id_x = min((cell_id_x + 2),self.cells_x-1);
        let min_cell_id_y = cell_id_y.saturating_sub(2);
        let max_cell_id_y = min((cell_id_y + 2),self.cells_y-1);
        debug_assert!(min_cell_id_x <= cell_id_x);
        debug_assert!(max_cell_id_x >= cell_id_x);
        debug_assert!(min_cell_id_y <= cell_id_y);
        debug_assert!(max_cell_id_y >= cell_id_y);
        for y in min_cell_id_y..max_cell_id_y+1 {
            for x in min_cell_id_x..max_cell_id_x+1 {
                match self.data[y*self.cells_x+x] {
                    None => (),
                    Some(other_id) => {
                        let other_sample = &samples[other_id];
                        let dstsqr = new_sample.dstsqr(other_sample);
                        if dstsqr < self.min_distance_sqr {
                            return Err(());
                        }
                    }
                }
            }
        }
        // no collission found
        samples.push(new_sample);
        debug_assert_eq!(self.data[cell_id_y*self.cells_x+cell_id_x], None);
        self.data[cell_id_y*self.cells_x+cell_id_x] = Some(samples.len()-1);
        debug_assert_eq!(samples[samples.len()-1], new_sample);
        Ok(samples.len()-1)
    }
}

pub fn blue_noise (dimensions: [usize; 2], min_distance: f64, k_abort: usize) -> Vec<Sample2D> {
    let mut rng = rand::thread_rng();
    let initial_sample = Sample2D {
        x: rng.gen::<usize>() % dimensions[0],
        y: rng.gen::<usize>() % dimensions[1],
    };
    let mut samples : Vec<Sample2D> = Vec::new();
    let mut bggrid = BackgroundGrid2D::new(dimensions[0], dimensions[1], min_distance);
    let initial_sample_id = bggrid.insert(initial_sample.x, initial_sample.y, &mut samples).unwrap();
    debug_assert_eq!(initial_sample_id, 0);
    let mut active : Vec<usize> = vec![initial_sample_id];
    while !active.is_empty() {
        let mut next_active = Vec::new();
        for current_id in active.iter() {
            let current_sample = samples[*current_id].clone();
            let mut found_something = false;
            for _ in 0..k_abort {
                let distance = rng.gen_range::<f64>(min_distance, 2_f64 * min_distance);
                let angle = rng.gen_range::<f64>(0_f64, 2_f64 * std::f64::consts::PI);
                let xoffs = angle.cos() * distance;
                let yoffs = angle.sin() * distance;
                let test_samp_x = current_sample.x as f64 + xoffs;
                let test_samp_y = current_sample.y as f64 + yoffs;
                if test_samp_x < 0. || test_samp_y < 0. {
                    continue;
                }
                match bggrid.insert(test_samp_x.round() as usize, test_samp_y.round() as usize, &mut samples) {
                    Ok(new_id) => {
                        next_active.push(new_id);
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
    let mut grid = BackgroundGrid2D::new(35,9,4.0);
    let mut samples = Vec::new();
    println!("cell size: {}", grid.cell_size);
    println!("grid dimensions: {} x {}", grid.cells_x, grid.cells_y);
    assert_eq!(grid.insert(0,9,&mut samples), Err(()));
    assert_eq!(grid.insert(0,0,&mut samples), Ok(0));
    assert_eq!(grid.insert(34,0,&mut samples), Ok(1));
    assert_eq!(grid.insert(0,8,&mut samples), Ok(2));
    assert_eq!(grid.insert(34,8,&mut samples), Ok(3));
}

#[cfg(test)]
mod tests {
    use super::*;
    fn get_image(radius: f64, size: usize) -> Vec<Vec<bool>> {
        let samples = blue_noise([size,size], radius, 30);
        let mut image = vec![vec![false; size]; size];
        for s in samples {
        	image[s.y][s.x] = true;
        }
        image
    }
    fn output_ppm(image: &Vec<Vec<bool>>) {
        let size = image.len();
        assert_eq!(image[0].len(),size);
        println!("P1");
        println!("{} {}", size, size);
        for y in 0..size {
            for x in 0..size {
                print!("{} ", if image[y][x] { "1" } else { "0" });
            }
            println!("");
        }
    }
    #[test]
    fn output_sanity() {
        let size : isize = 128;
        let radius : isize = 8;
        let image = get_image(radius as f64, size as usize);
        for y in 0..size {
            for x in 0..size {
                if image[y as usize][x as usize] {
                    for dy in -radius..radius {
                        for dx in -radius..radius {
                            if dx == 0 && dy == 0 {
                                continue;
                            }
                            if -dy > y as isize || -dx > x as isize {
                                continue;
                            }
                            image.get((y+dy) as usize).map(|line| match line.get((x+dx) as usize) {
                                        Some(&true) => {
                                            assert!(dx*dx+dy*dy >= radius*radius);
                                        },
                                        _ => {}
                            });
                        }
                    }
                }
            }
        }
    }
    #[test]
    #[ignore]
    fn ppm() {
        let image = get_image(8.0, 128);
        output_ppm(&image);
    }
}
