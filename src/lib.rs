extern crate rand;
use rand::Rng;

use std::cmp::min;

pub trait NoiseContainer2D {
    fn width(&self) -> usize;
    fn height(&self) -> usize;
    fn put_sample(&mut self, x: usize, y: usize);
}

impl NoiseContainer2D for Vec<Vec<bool>> {
    fn width(&self) -> usize {
        self[0].len()
    }
    fn height(&self) -> usize {
        self.len()
    }
    fn put_sample(&mut self, x: usize, y: usize) {
        self[y][x] = true;
    }
}

#[derive(Clone,Copy,PartialEq,Eq,Debug)]
struct Sample2D {
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
    min_distance_sqr: usize,
    cell_size: usize,
}

impl BackgroundGrid2D {
    pub fn new(width: usize, height: usize, min_distance: f64) -> BackgroundGrid2D {
        assert!(min_distance > 0.0_f64);
        let max_cell_size = min_distance / 2.0_f64.sqrt();
        let cells_x : usize = (width as f64 / max_cell_size).ceil() as usize;
        let cells_y : usize = (height as f64 / max_cell_size).ceil() as usize;
        BackgroundGrid2D {
            data: vec![None; cells_x*cells_y],
            cells_x: cells_x,
            cells_y: cells_y,
            min_distance_sqr: (min_distance * min_distance).ceil() as usize,
            cell_size: max_cell_size as usize,
        }
    }
    pub fn insert(&mut self, sample_x: usize, sample_y: usize, samples: &mut Vec<Sample2D>) -> Result<usize,()> {
        let new_sample = Sample2D { x: sample_x, y: sample_y };
        let cell_id_x = sample_x / self.cell_size;
        let cell_id_y = sample_y / self.cell_size;
        if cell_id_x >= self.cells_x || cell_id_y >= self.cells_y {
            return Err(());
        }
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

pub fn blue_noise (image: &mut NoiseContainer2D, min_distance: f64, k_abort: usize) {
    let mut rng = rand::thread_rng();
    let initial_sample = Sample2D {
        x: rng.gen::<usize>() % image.width(),
        y: rng.gen::<usize>() % image.height(),
    };
    let mut samples : Vec<Sample2D> = Vec::new();
    let mut bggrid = BackgroundGrid2D::new(image.width(), image.height(), min_distance);
    let initial_sample_id = bggrid.insert(initial_sample.x, initial_sample.y, &mut samples).unwrap();
    debug_assert_eq!(initial_sample_id, 0);
    let mut active : Vec<usize> = vec![initial_sample_id];
    while !active.is_empty() {
        let mut next_active = Vec::new();
        for current_id in active.iter() {
            let current_sample = samples[*current_id].clone();
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
                    }
                    Err(_) => {
                        // wait for the next iteration
                    }
                }
            }
        }
        active = next_active;
    }
    for samp in samples {
        image.put_sample(samp.x, samp.y);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn get_image(radius: i32, size: usize) -> Vec<Vec<bool>> {
        let mut image = vec![vec![false; size]; size];
        blue_noise(&mut image, radius as f64, 30);
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
        let size : usize = 128;
        let radius : i32 = 8;
        let image = get_image(radius, size);
        for y in 0..size {
            for x in 0..size {
                if image[y][x] {
                    for dy in -radius..radius {
                        for dx in -radius..radius {
                            if dx == 0 && dy == 0 {
                                continue;
                            }
                            image.get(y).map(|line| match line.get(x) {
                                        Some(&true) => {
                                            if !(dx*dx+dy*dy >= radius*radius) {
                                                println!("dx: {}\ndy: {}",dx,dy);
                                            }
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
        let size : usize = 128;
        let radius : i32 = 8;
        let image = get_image(radius, size);
        output_ppm(&image);
    }
}
